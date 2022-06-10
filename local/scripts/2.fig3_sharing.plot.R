##### GMK April 2019 #####
#### packages ####
rm(list=ls())
dev.off()

library(ggplot2)
library(dplyr)
library(tidyr)
library(qqman)
library(cowplot)
library(data.table)
library(stringr)
library(devtools)
library(ggman)
library(patchwork)
options(scipen = 999)

##### functions ##### 
# modify ggmanHighlight so that it runs even if no matches with list (some chr no matches)
check.input.ggmanHighlight <- function(){
  ## ggmanPlot input
  if(! any(class(ggmanPlot) == "ggman")){
    stop("The ggmanPlot input is not a ggman object")
  }
  ## highlight input
  if(! is.vector(highlight)){
    stop(paste0("The highlight input ",highlight," is not a vector object"))
  }
  
  ## colour input
  ## lineColor input
  ## Thanks to Sacha Epskamp for isColor function. 
  ## Reference:http://stackoverflow.com/questions/13289009/check-if-character-string-is-a-valid-color-representation/13290832
  isColor <- function(x)
  {
    res <- try(col2rgb(x),silent=TRUE)
    return(!"try-error"%in%class(res))
  }
  if(! isColor(colour)){
    stop(paste0("\'",colour,"\'"," is not a valid color"))
  }
  
}

ggmanHighlight <- function(ggmanPlot, highlight,colour = "#56B4E9", ...){
  ##input checks
  environment(check.input.ggmanHighlight) <- environment()
  check.input.ggmanHighlight()
  
  dfm <- ggmanPlot[[1]]
  dfm <- dfm[dfm$snp %in% highlight,]
  ggmanPlot +
    scale_colour_grey(start = 0.5,end = 0.6) +
    geom_point(data = dfm,colour= colour,...)
}

zscore <- function(x) { 
  x<- as.numeric(as.character(x))
  (x-mean(x))/sd(x) }

# expands dataset +-x windows
bp.list.expanded <- function(dat, x){
  BP.wg <- dat[["BP.wg"]]
  # x= number of 1kb windows, convert to bp
  y <- x*1000
  seq(from=as.numeric(as.character(BP.wg))-y, to=as.numeric(as.character(BP.wg))+y, by=1000)}

## shift one row (either dir)
rowShift <- function(x, shiftLen = 1L) {
  r <- (1L + shiftLen):(length(x) + shiftLen)
  r[r<1] <- NA
  return(x[r])
}



#################################### 0. data prep #####
setwd("git/2021-altitude-heliconius/")
ref.scaff.era <- read.table(" local/data/ref/Heliconius_erato_demophoon_v1_-_scaffolds.fa.fai", row.names = NULL)
ref.scaff.mel <- read.table(" local/data/ref/Hmel2.5.scaffolds.fa.fai", row.names = NULL)

era.all.pbs <- read.csv("local/data/pbs.out/era.all.pbs.recomb.csv")
mel.all.pbs <- read.csv("local/data/pbs.out/mel.all.pbs.recomb.csv")

head(era.all.pbs)
era.all.pbs$CHR <- as.numeric(as.character(substr(era.all.pbs$scaff,7,8))); unique(era.all.pbs$CHR)
mel.all.pbs$CHR <- as.numeric(as.character(substr(mel.all.pbs$scaff,6,7))); unique(mel.all.pbs$CHR)

# create df with start ends of shared outlier regions
names(era.all.pbs)
# this is not ideal as some across sites get split by those regions that dont overlap, but wont see it on the plot
era.shared.regions <-  summarise(group_by(era.all.pbs, sharing.type, shdr.para.east.id, shdr.para.west.id, CHR),
          start.BP.wg=min(BP.wg),
          end.BP.wg=max(BP.wg)); head(era.shared.regions)

mel.shared.regions <-  summarise(group_by(mel.all.pbs, sharing.type, shdr.para.east.id, shdr.para.west.id, CHR),
                                 start.BP.wg=min(BP.wg),
                                 end.BP.wg=max(BP.wg)); head(mel.shared.regions)

# add whether region forms part of a hdr

era.all.pbs$is.hdr<- if_else(!(is.na(era.all.pbs$hdr.CoE))|!(is.na(era.all.pbs$hdr.EcE))|!(is.na(era.all.pbs$hdr.CoW))|!(is.na(era.all.pbs$hdr.EcW)), "hdr", "background" )
mel.all.pbs$is.hdr<- if_else(!(is.na(mel.all.pbs$hdr.CoE))|!(is.na(mel.all.pbs$hdr.EcE))|!(is.na(mel.all.pbs$hdr.CoW))|!(is.na(mel.all.pbs$hdr.EcW)), "hdr", "background" )

## recomb rates
mel.rec.summ.S <-  read.csv("local/data/recomb/mel.rec.summ.scaled.csv")
era.rec.summ.S <-  read.csv("local/data/recomb/era.rec.summ.scaled.csv")
str(mel.rec.summ.S ); str(era.rec.summ.S )

############### prep recomb splines ###############
# Fit a smooth spline across rec rates per chromosome
chr.splines.era <- list()
chr.splines.pred.era <- list()
for( i in 1:length(unique(era.rec.summ.S$chromosome))) {
  chr.splines.era[[i]] <- smooth.spline(subset(era.rec.summ.S, chromosome==i)$BP.wg, subset(era.rec.summ.S, chromosome==i)$meanS, spar = 0.5)
  chr.splines.pred.era[[i]] <- predict(chr.splines.era[[i]], subset(era.rec.summ.S, chromosome==i)[,"BP.wg"])$y 
}
era.rec.summ.S$meanS.spline <- unlist(rbindlist(lapply(chr.splines.pred.era, as.data.table)))


chr.splines.mel <- list()
chr.splines.pred.mel <- list()
for( i in 1:length(unique(mel.rec.summ.S$chr))) {
  chr.splines.mel[[i]] <- smooth.spline(subset(mel.rec.summ.S, chr==i)$BP.wg, subset(mel.rec.summ.S, chr==i)$meanS, spar = 0.5)
  chr.splines.pred.mel[[i]] <- predict(chr.splines.mel[[i]], subset(mel.rec.summ.S, chr==i)[,"BP.wg"])$y 
}
mel.rec.summ.S$meanS.spline <- unlist(rbindlist(lapply(chr.splines.pred.mel, as.data.table)))


#

#################################### erato PBS plotting ####################################
# prep chr midpoints for plotting
axisdf.era = era.all.pbs %>% dplyr::group_by(CHR) %>% dplyr::summarize(center=( max(BP.wg) + min(BP.wg) ) / 2 ); axisdf.era
axisdf.era$CHR[21]<-"Z";axisdf.era$CHR[21]
# create start ends of chr for plotting
head(ref.scaff)
ref.scaff.era<-ref.scaff.era[,1:2]
names(ref.scaff.era)<-c("scaff","length")
ref.scaff.era$add<-c(0,cumsum(ref.scaff.era$length)[-length(ref.scaff.era$scaff)])
ref.scaff.era<-ref.scaff.era[ref.scaff.era$scaff!="Hera_complete_mtDNA",]
ref.scaff.era$CHR<-as.integer(substr(ref.scaff.era$scaff,7,8))
#ref.scaff.era<-ref.scaff.era[!grepl(ref.scaff.era$scaff,pattern = "Hera200"),]
era.ref.chr.pos <- ref.scaff.era; era.ref.chr.pos$start <- era.ref.chr.pos$add
era.ref.chr.pos$end <-  era.ref.chr.pos$start +era.ref.chr.pos$length; head(era.ref.chr.pos)
era.ref.chr.pos <- summarise(group_by(era.ref.chr.pos, CHR),start=min(start),end=max(end)); era.ref.chr.pos
era.ref.chr.pos <- subset(era.ref.chr.pos, era.ref.chr.pos$CHR!="")

# margins #top, right, bottom, left
names(era.all.pbs)
co.w.era.p <- ggplot(subset(era.all.pbs), aes(x=BP.wg, y=zfst.2)) +
  geom_rect(inherit.aes = F, data=era.ref.chr.pos, aes(xmin=start, xmax=end, ymin=0,ymax=30), colour="transparent", fill=c(rep(c("white", "grey40"), 10 ), "white"), alpha=.1 )+
  geom_hline(yintercept = 4 , colour="grey40", lty="dashed",alpha=.99)+
  geom_point(color="darkgrey", alpha=0.9, size=0.02) +
  #geom_rect(inherit.aes = F,data=zoom.sections.spread.era, aes(xmin=start, xmax=end, ymin=0.2,ymax=3.1), colour="transparent", fill="yellow", alpha=0.8 )+
  geom_point(inherit.aes = T, aes(x=BP.wg, y=zfst.2),   data=subset(era.all.pbs, !(is.na(hdr.CoW)) & is.hdr=="hdr"), alpha=0.7, size=.1, colour="grey40") +
  geom_point(inherit.aes = T, aes(x=BP.wg, y=zfst.2),   data=subset(era.all.pbs, sharing.type=="within.west"), alpha=0.7, size=.4, colour="#0372B2") +
  geom_point(inherit.aes = T, aes(x=BP.wg, y=zfst.2),   data=subset(era.all.pbs, sharing.type=="across"), alpha=0.7, size=.4, colour="#D65D00") +
  geom_text(data = axisdf.era, aes(x=center, y=27, label=CHR))+
  scale_x_continuous( expand = c(0, 0),label = axisdf.era$CHR, breaks= axisdf.era$center ) +
  scale_y_continuous( expand = c(0, 0), limits = c(0,30), labels = c(0, 10, 20, 30) ) +    # remove space between plot area and x axis
  ylab("Col. W\nzFst")+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=.5),legend.position="none",
        plot.margin = unit(c(0.1,0.1,0,0.1), "lines"),axis.text = element_text(size=10), axis.text.x = element_blank(),axis.ticks.x = element_blank(),
        axis.text.y = element_text(size=10),axis.title.x = element_blank(), panel.background = element_rect(fill = NA), plot.background = element_rect(fill = NA, color = NA),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.background = element_rect(fill = "transparent"), legend.box.background = element_rect(fill = "transparent")); co.w.era.p

ec.w.era.p <- ggplot(subset(era.all.pbs), aes(x=BP.wg, y=zPBS0.4)) +
  geom_rect(inherit.aes = F, data=era.ref.chr.pos, aes(xmin=start, xmax=end, ymin=0,ymax=30), colour="transparent", fill=c(rep(c("white", "grey40"), 10 ), "white"), alpha=.1 )+
  geom_hline(yintercept = 4 , colour="grey40", lty="dashed",alpha=.99)+
  geom_point(color="darkgrey", alpha=0.9, size=0.02) +
  #geom_rect(inherit.aes = F,data=zoom.sections.spread.era, aes(xmin=start, xmax=end, ymin=0.2,ymax=3.1), colour="transparent", fill="yellow", alpha=0.8 )+
  geom_point(inherit.aes = T, aes(x=BP.wg, y=zPBS0.4),   data=subset(era.all.pbs, !(is.na(hdr.EcW)) & is.hdr=="hdr"), alpha=0.7, size=.1, colour="grey40") +
  
  geom_point(inherit.aes = T, aes(x=BP.wg, y=zPBS0.4),   data=subset(era.all.pbs, sharing.type=="within.west"), alpha=0.7, size=.4, colour="#0372B2") +
  geom_point(inherit.aes = T, aes(x=BP.wg, y=zPBS0.4),   data=subset(era.all.pbs, sharing.type=="across"), alpha=0.7, size=.4, colour="#D65D00") +
  scale_x_continuous( expand = c(0, 0),label = axisdf.era$CHR, breaks= axisdf.era$center ) +
  scale_y_continuous( expand = c(0, 0), limits = c(0,30), labels = c(0, 10, 20, "") ) +    # remove space between plot area and x axis
  ylab("Ecu. W\nzPBS high")+
  #theme_bw()+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=.5),legend.position="none",
        plot.margin = unit(c(0,0.1,0,0.1), "lines"),axis.text = element_text(size=10), axis.text.x = element_blank(),axis.ticks.x = element_blank(),
        axis.text.y = element_text(size=10),axis.title.x = element_blank(), panel.background = element_rect(fill = NA), plot.background = element_rect(fill = NA, color = NA),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.background = element_rect(fill = "transparent"), legend.box.background = element_rect(fill = "transparent")); ec.w.era.p

co.e.era.p <- ggplot(subset(era.all.pbs), aes(x=BP.wg, y=zPBS0.1)) +
  geom_rect(inherit.aes = F, data=era.ref.chr.pos, aes(xmin=start, xmax=end, ymin=0,ymax=30), colour="transparent", fill=c(rep(c("white", "grey40"), 10 ), "white"), alpha=.1 )+
  geom_hline(yintercept = 4 , colour="grey40", lty="dashed",alpha=.99)+
  geom_point(color="darkgrey", alpha=0.9, size=0.02) +
  #geom_rect(inherit.aes = F,data=zoom.sections.spread.era, aes(xmin=start, xmax=end, ymin=0.2,ymax=3.1), colour="transparent", fill="yellow", alpha=0.8 )+
  geom_point(inherit.aes = T, aes(x=BP.wg, y=zPBS0.1),   data=subset(era.all.pbs, !(is.na(hdr.CoE)) & is.hdr=="hdr"), alpha=0.7, size=.1, colour="grey40") +
  geom_point(inherit.aes = T, aes(x=BP.wg, y=zPBS0.1),   data=subset(era.all.pbs, sharing.type=="within.east"), alpha=0.7, size=.4, colour="#049E73") +
  geom_point(inherit.aes = T, aes(x=BP.wg, y=zPBS0.1),   data=subset(era.all.pbs, sharing.type=="across"), alpha=0.7, size=.4, colour="#D65D00") +
  scale_x_continuous( expand = c(0, 0),label = axisdf.era$CHR, breaks= axisdf.era$center ) +
  scale_y_continuous( expand = c(0, 0), limits = c(0,30), labels = c(0, 10, 20, "") ) +    # remove space between plot area and x axis
  ylab("Col. E\nzPBS high")+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=.5),legend.position="none",
        plot.margin = unit(c(0,0.1,0,0.1), "lines"),axis.text = element_text(size=10), axis.text.x = element_blank(),axis.ticks.x = element_blank(),
        axis.text.y = element_text(size=10),axis.title.x = element_blank(), panel.background = element_rect(fill = "transparent"), plot.background = element_rect(fill = "transparent", color = NA),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.background = element_rect(fill = "transparent"), legend.box.background = element_rect(fill = "transparent")); co.e.era.p

ec.e.era.p <- ggplot(subset(era.all.pbs), aes(x=BP.wg, y=zPBS0.3)) +
  geom_rect(inherit.aes = F, data=era.ref.chr.pos, aes(xmin=start, xmax=end, ymin=0,ymax=30), colour="transparent", fill=c(rep(c("white", "grey40"), 10 ), "white"), alpha=.1 )+
  geom_hline(yintercept = 4 , colour="grey40", lty="dashed",alpha=.99)+
  geom_point(color="darkgrey", alpha=0.9, size=0.02) +
  geom_point(inherit.aes = T, aes(x=BP.wg, y=zPBS0.3),   data=subset(era.all.pbs, !(is.na(hdr.EcE)) & is.hdr=="hdr"), alpha=0.7, size=.1, colour="grey40") +
  geom_point(inherit.aes = T, aes(x=BP.wg, y=zPBS0.3),   data=subset(era.all.pbs, sharing.type=="within.east"), alpha=0.7, size=.4, colour="#049E73") +
  geom_point(inherit.aes = T, aes(x=BP.wg, y=zPBS0.3),   data=subset(era.all.pbs, sharing.type=="across"), alpha=0.7, size=.4, colour="#D65D00") +
  scale_x_continuous( expand = c(0, 0),label = axisdf.era$CHR, breaks= axisdf.era$center ) +
  scale_y_continuous( expand = c(0, 0), limits = c(0,30), labels = c(0, 10, 20, "") ) +    # remove space between plot area and x axis
  ylab("Ecu. E\nzPBS high")+
  #theme_bw()+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=.5),legend.position="none",
        plot.margin = unit(c(0,0.1,0,0.1), "lines"),axis.text = element_text(size=10), axis.text.x = element_blank(),axis.ticks.x = element_blank(),
        axis.text.y = element_text(size=10),axis.title.x = element_blank(), panel.background = element_rect(fill = "transparent"), plot.background = element_rect(fill = "transparent", color = NA),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.background = element_rect(fill = "transparent"), legend.box.background = element_rect(fill = "transparent")); ec.e.era.p


# use patchwork

co.w.era.p / ec.w.era.p / co.e.era.p / ec.e.era.p

#################################### melpomene PBS plotting ####################################

# prep chr midpoints for plotting
axisdf.mel = mel.all.pbs %>% dplyr::group_by(CHR) %>% dplyr::summarize(center=( max(BP.wg) + min(BP.wg) ) / 2 ); axisdf.mel
axisdf.mel$CHR[21]<-"Z";axisdf.mel$CHR[21]
# create start ends of chr for plotting
head(ref.scaff)
ref.scaff.mel<-ref.scaff.mel[,1:2]
names(ref.scaff.mel)<-c("scaff","length")
ref.scaff.mel$add<-c(0,cumsum(ref.scaff.mel$length)[-length(ref.scaff.mel$scaff)])
ref.scaff.mel<-ref.scaff.mel[ref.scaff.mel$scaff!="Hmel_complete_mtDNA",]
ref.scaff.mel$CHR<-as.integer(substr(ref.scaff.mel$scaff,6,7))
ref.scaff.mel<-ref.scaff.mel[!grepl(ref.scaff.mel$scaff,pattern = "Hmel200"),]
mel.ref.chr.pos <- ref.scaff.mel; mel.ref.chr.pos$start <- mel.ref.chr.pos$add
mel.ref.chr.pos$end <-  mel.ref.chr.pos$start +mel.ref.chr.pos$length; head(mel.ref.chr.pos)
mel.ref.chr.pos <- summarise(group_by(mel.ref.chr.pos, CHR),start=min(start),end=max(end)); mel.ref.chr.pos
mel.ref.chr.pos <- subset(mel.ref.chr.pos, mel.ref.chr.pos$CHR!="")

axisdf.mel = mel.all.pbs %>% dplyr::group_by(CHR) %>% dplyr::summarize(center=( max(BP.wg) + min(BP.wg) ) / 2 ); axisdf.mel
axisdf.mel$CHR[21]<-"Z";axisdf.mel$CHR[21]

# create liftedover start ends of chr for plotting
mel.lift.ref.chr.pos <- summarise(group_by(era.ref.chr.pos, CHR),start=min(start),end=max(end)); mel.lift.ref.chr.pos
mel.lift.ref.chr.pos <- subset(mel.lift.ref.chr.pos, mel.lift.ref.chr.pos$CHR!="")
axisdf.mel.lift = mel.all.pbs %>% dplyr::group_by(era.CHR) %>% dplyr::summarize(center=( max(era.BP.wg) + min(era.BP.wg) ) / 2 ); axisdf.mel.lift
axisdf.mel.lift$era.CHR[21]<-"Z";axisdf.mel.lift$era.CHR[21]

#### melpomene not lift #####
# margins #top, right, bottom, left
names(mel.all.pbs)
co.w.mel.p <- ggplot(subset(mel.all.pbs), aes(x=BP.wg, y=zPBS0.2)) +
  geom_rect(inherit.aes = F, data=mel.ref.chr.pos, aes(xmin=start, xmax=end, ymin=0,ymax=20), colour="transparent", fill=c(rep(c("white", "grey40"), 10 ), "white"), alpha=.1 )+
  geom_hline(yintercept = 4 , colour="grey40", lty="dashed",alpha=.99)+
  geom_point(color="darkgrey", alpha=0.9, size=0.02) +
  geom_point(inherit.aes = T, aes(x=BP.wg, y=zPBS0.2),   data=subset(mel.all.pbs, !(is.na(hdr.CoW)) & is.hdr=="hdr"), alpha=0.7, size=.1, colour="grey40") +
  geom_point(inherit.aes = T, aes(x=BP.wg, y=zPBS0.2),   data=subset(mel.all.pbs, sharing.type=="within.west"), alpha=0.7, size=.4, colour="#0372B2") +
  geom_point(inherit.aes = T, aes(x=BP.wg, y=zPBS0.2),   data=subset(mel.all.pbs, sharing.type=="across"), alpha=0.7, size=.4, colour="#D65D00") +
  geom_text(data = axisdf.mel, aes(x=center, y=18, label=CHR))+
  scale_x_continuous( expand = c(0, 0),label = axisdf.mel$CHR, breaks= axisdf.mel$center ) +
  scale_y_continuous( expand = c(0, 0), limits = c(0,20), labels = c(0, 10, 20), breaks = c(0,10,20) ) +    # remove space between plot area and x axis
  ylab("Col. W\nzPBS high")+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=.5),legend.position="none",
        plot.margin = unit(c(0.5,0.1,0,0.1), "lines"),axis.text = element_text(size=10), axis.text.x = element_blank(),axis.ticks.x = element_blank(),
        axis.text.y = element_text(size=10),axis.title.x = element_blank(), panel.background = element_rect(fill = NA), plot.background = element_rect(fill = NA, color = NA),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.background = element_rect(fill = "transparent"), legend.box.background = element_rect(fill = "transparent")); co.w.mel.p

ec.w.mel.p <- ggplot(subset(mel.all.pbs), aes(x=BP.wg, y=zfst.4)) +
  geom_rect(inherit.aes = F, data=mel.ref.chr.pos, aes(xmin=start, xmax=end, ymin=0,ymax=20), colour="transparent", fill=c(rep(c("white", "grey40"), 10 ), "white"), alpha=.1 )+
  geom_hline(yintercept = 4 , colour="grey40", lty="dashed",alpha=.99)+
  geom_point(color="darkgrey", alpha=0.9, size=0.02) +
  geom_point(inherit.aes = T, aes(x=BP.wg, y=zfst.4),   data=subset(mel.all.pbs, !(is.na(hdr.EcW)) & is.hdr=="hdr"), alpha=0.7, size=.1, colour="grey40") +
  geom_point(inherit.aes = T, aes(x=BP.wg, y=zfst.4),   data=subset(mel.all.pbs, sharing.type=="within.west"), alpha=0.7, size=.4, colour="#0372B2") +
  geom_point(inherit.aes = T, aes(x=BP.wg, y=zfst.4),   data=subset(mel.all.pbs, sharing.type=="across"), alpha=0.7, size=.4, colour="#D65D00") +
  scale_x_continuous( expand = c(0, 0),label = axisdf.mel$CHR, breaks= axisdf.mel$center ) +
  scale_y_continuous( expand = c(0, 0), limits = c(0,20), labels = c(0, 10, ""), breaks = c(0,10,20) ) +    # remove space between plot area and x axis
  ylab("Ecu. W\nzFst")+
  #theme_bw()+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=.5),legend.position="none",
        plot.margin = unit(c(0,0.1,0,0.1), "lines"),axis.text = element_text(size=10), axis.text.x = element_blank(),axis.ticks.x = element_blank(),
        axis.text.y = element_text(size=10),axis.title.x = element_blank(), panel.background = element_rect(fill = NA), plot.background = element_rect(fill = NA, color = NA),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.background = element_rect(fill = "transparent"), legend.box.background = element_rect(fill = "transparent")); ec.w.mel.p

co.e.mel.p <- ggplot(subset(mel.all.pbs), aes(x=BP.wg, y=zPBS0.1)) +
  geom_rect(inherit.aes = F, data=mel.ref.chr.pos, aes(xmin=start, xmax=end, ymin=0,ymax=20), colour="transparent", fill=c(rep(c("white", "grey40"), 10 ), "white"), alpha=.1 )+
  geom_hline(yintercept = 4 , colour="grey40", lty="dashed",alpha=.99)+
  geom_point(color="darkgrey", alpha=0.9, size=0.02) +
  geom_point(inherit.aes = T, aes(x=BP.wg, y=zPBS0.1),   data=subset(mel.all.pbs, !(is.na(hdr.CoE)) & is.hdr=="hdr"), alpha=0.7, size=.1, colour="grey40") +
  geom_point(inherit.aes = T, aes(x=BP.wg, y=zPBS0.1),   data=subset(mel.all.pbs, sharing.type=="within.east"), alpha=0.7, size=.4, colour="#049E73") +
  geom_point(inherit.aes = T, aes(x=BP.wg, y=zPBS0.1),   data=subset(mel.all.pbs, sharing.type=="across"), alpha=0.7, size=.4, colour="#D65D00") +
  scale_x_continuous( expand = c(0, 0),label = axisdf.mel$CHR, breaks= axisdf.mel$center ) +
  scale_y_continuous( expand = c(0, 0), limits = c(0,20), labels = c(0, 10, ""), breaks = c(0,10,20) ) +    # remove space between plot area and x axis
  ylab("Col. E\nzPBS high")+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=.5),legend.position="none",
        plot.margin = unit(c(0,0.1,0,0.1), "lines"),axis.text = element_text(size=10), axis.text.x = element_blank(),axis.ticks.x = element_blank(),
        axis.text.y = element_text(size=10),axis.title.x = element_blank(), panel.background = element_rect(fill = "transparent"), plot.background = element_rect(fill = "transparent", color = NA),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.background = element_rect(fill = "transparent"), legend.box.background = element_rect(fill = "transparent")); co.e.mel.p

ec.e.mel.p <- ggplot(subset(mel.all.pbs), aes(x=BP.wg, y=zPBS0.3)) +
  geom_rect(inherit.aes = F, data=mel.ref.chr.pos, aes(xmin=start, xmax=end, ymin=0,ymax=20), colour="transparent", fill=c(rep(c("white", "grey40"), 10 ), "white"), alpha=.1 )+
  geom_hline(yintercept = 4 , colour="grey40", lty="dashed",alpha=.99)+
  geom_point(color="darkgrey", alpha=0.9, size=0.02) +
  #geom_rect(inherit.aes = F,data=zoom.sections.spread.mel, aes(xmin=start, xmax=end, ymin=0.2,ymax=3.1), colour="transparent", fill="yellow", alpha=0.8 )+
  geom_point(inherit.aes = T, aes(x=BP.wg, y=zPBS0.3),   data=subset(mel.all.pbs, !(is.na(hdr.EcE)) & is.hdr=="hdr"), alpha=0.7, size=.1, colour="grey40") +
  geom_point(inherit.aes = T, aes(x=BP.wg, y=zPBS0.3),   data=subset(mel.all.pbs, sharing.type=="within.east"), alpha=0.7, size=.4, colour="#049E73") +
  geom_point(inherit.aes = T, aes(x=BP.wg, y=zPBS0.3),   data=subset(mel.all.pbs, sharing.type=="across"), alpha=0.7, size=.4, colour="#D65D00") +
  scale_x_continuous( expand = c(0, 0),label = axisdf.mel$CHR, breaks= axisdf.mel$center ) +
  scale_y_continuous( expand = c(0, 0), limits = c(0,20), labels = c(0, 10, ""), breaks = c(0,10,20) ) +    # remove space between plot area and x axis
  ylab("Ecu. E\nzPBS high")+
  #theme_bw()+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=.5),legend.position="none",
        plot.margin = unit(c(0,0.1,0,0.1), "lines"),axis.text = element_text(size=10), axis.text.x = element_blank(),axis.ticks.x = element_blank(),
        axis.text.y = element_text(size=10),axis.title.x = element_blank(), panel.background = element_rect(fill = "transparent"), plot.background = element_rect(fill = "transparent", color = NA),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.background = element_rect(fill = "transparent"), legend.box.background = element_rect(fill = "transparent")); ec.e.mel.p

#### melpomene liftover #####
co.w.mel.p.lift <- ggplot(subset(mel.all.pbs), aes(x=era.BP.wg, y=zPBS0.2)) +
  geom_rect(inherit.aes = F, data=mel.lift.ref.chr.pos, aes(xmin=start, xmax=end, ymin=0,ymax=20), colour="transparent", fill=c(rep(c("white", "grey40"), 10 ), "white"), alpha=.1 )+
  geom_hline(yintercept = 4 , colour="grey40", lty="dashed",alpha=.99)+
  geom_point(color="darkgrey", alpha=0.9, size=0.02) +
  geom_point(inherit.aes = T, aes(x=era.BP.wg, y=zPBS0.2),   data=subset(mel.all.pbs, !(is.na(hdr.CoW)) & is.hdr=="hdr"), alpha=0.7, size=.1, colour="grey40") +
  geom_point(inherit.aes = T, aes(x=era.BP.wg, y=zPBS0.2),   data=subset(mel.all.pbs, sharing.type=="within.west"), alpha=0.7, size=.4, colour="#0372B2") +
  geom_point(inherit.aes = T, aes(x=era.BP.wg, y=zPBS0.2),   data=subset(mel.all.pbs, sharing.type=="across"), alpha=0.7, size=.4, colour="#D65D00") +
  geom_text(data = axisdf.mel.lift, aes(x=center, y=18, label=era.CHR))+
  scale_x_continuous( expand = c(0, 0),label = axisdf.mel.lift$era.CHR, breaks= axisdf.mel.lift$center ) +
  scale_y_continuous( expand = c(0, 0), limits = c(0,20), labels = c(0, 10, 20), breaks = c(0,10,20) ) +    # remove space between plot area and x axis
  ylab("zPBS high")+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=.5),legend.position="none",
        plot.margin = unit(c(0.5,0.1,0,0.1), "lines"),axis.text = element_text(size=10), axis.text.x = element_blank(),axis.ticks.x = element_blank(),
        axis.text.y = element_text(size=10),axis.title.x = element_blank(), panel.background = element_rect(fill = NA), plot.background = element_rect(fill = NA, color = NA),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.background = element_rect(fill = "transparent"), legend.box.background = element_rect(fill = "transparent")); co.w.mel.p.lift


ec.w.mel.p.lift <- ggplot(subset(mel.all.pbs), aes(x=era.BP.wg, y=zfst.4)) +
  geom_rect(inherit.aes = F, data=mel.lift.ref.chr.pos, aes(xmin=start, xmax=end, ymin=0,ymax=20), colour="transparent", fill=c(rep(c("white", "grey40"), 10 ), "white"), alpha=.1 )+
  geom_hline(yintercept = 4 , colour="grey40", lty="dashed",alpha=.99)+
  geom_point(color="darkgrey", alpha=0.9, size=0.02) +
  #geom_rect(inherit.aes = F,data=zoom.sections.spread.mel, aes(xmin=start, xmax=end, ymin=0.2,ymax=3.1), colour="transparent", fill="yellow", alpha=0.8 )+
  geom_point(inherit.aes = T, aes(x=era.BP.wg, y=zfst.4),   data=subset(mel.all.pbs, !(is.na(hdr.EcW)) & is.hdr=="hdr"), alpha=0.7, size=.1, colour="grey40") +
  geom_point(inherit.aes = T, aes(x=era.BP.wg, y=zfst.4),   data=subset(mel.all.pbs, sharing.type=="within.west"), alpha=0.7, size=.4, colour="#0372B2") +
  geom_point(inherit.aes = T, aes(x=era.BP.wg, y=zfst.4),   data=subset(mel.all.pbs, sharing.type=="across"), alpha=0.7, size=.4, colour="#D65D00") +
  scale_x_continuous( expand = c(0, 0),label = axisdf.mel.lift$era.CHR, breaks= axisdf.mel.lift$center ) +
  scale_y_continuous( expand = c(0, 0), limits = c(0,20), labels = c(0, 10, ""), breaks = c(0,10,20) ) +    # remove space between plot area and x axis
  ylab("zFst")+
  #theme_bw()+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=.5),legend.position="none",
        plot.margin = unit(c(0,0.1,0,0.1), "lines"),axis.text = element_text(size=10), axis.text.x = element_blank(),axis.ticks.x = element_blank(),
        axis.text.y = element_text(size=10),axis.title.x = element_blank(), panel.background = element_rect(fill = NA), plot.background = element_rect(fill = NA, color = NA),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.background = element_rect(fill = "transparent"), legend.box.background = element_rect(fill = "transparent")); ec.w.mel.p.lift

co.e.mel.p.lift <- ggplot(subset(mel.all.pbs), aes(x=era.BP.wg, y=zPBS0.1)) +
  geom_rect(inherit.aes = F, data=mel.lift.ref.chr.pos, aes(xmin=start, xmax=end, ymin=0,ymax=20), colour="transparent", fill=c(rep(c("white", "grey40"), 10 ), "white"), alpha=.1 )+
  geom_hline(yintercept = 4 , colour="grey40", lty="dashed",alpha=.99)+
  geom_point(color="darkgrey", alpha=0.9, size=0.02) +
  geom_point(inherit.aes = T, aes(x=era.BP.wg, y=zPBS0.1),   data=subset(mel.all.pbs, !(is.na(hdr.CoE)) & is.hdr=="hdr"), alpha=0.7, size=.1, colour="grey40") +
  geom_point(inherit.aes = T, aes(x=era.BP.wg, y=zPBS0.1),   data=subset(mel.all.pbs, sharing.type=="within.east"), alpha=0.7, size=.4, colour="#049E73") +
  geom_point(inherit.aes = T, aes(x=era.BP.wg, y=zPBS0.1),   data=subset(mel.all.pbs, sharing.type=="across"), alpha=0.7, size=.4, colour="#D65D00") +
  scale_x_continuous( expand = c(0, 0),label = axisdf.mel.lift$era.CHR, breaks= axisdf.mel.lift$center ) +
  scale_y_continuous( expand = c(0, 0), limits = c(0,20), labels = c(0, 10, ""), breaks = c(0,10,20) ) +    # remove space between plot area and x axis
  ylab("zPBS high")+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=.5),legend.position="none",
        plot.margin = unit(c(0,0.1,0,0.1), "lines"),axis.text = element_text(size=10), axis.text.x = element_blank(),axis.ticks.x = element_blank(),
        axis.text.y = element_text(size=10),axis.title.x = element_blank(), panel.background = element_rect(fill = "transparent"), plot.background = element_rect(fill = "transparent", color = NA),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.background = element_rect(fill = "transparent"), legend.box.background = element_rect(fill = "transparent")); co.e.mel.p.lift

ec.e.mel.p.lift <- ggplot(subset(mel.all.pbs), aes(x=era.BP.wg, y=zPBS0.3)) +
  geom_rect(inherit.aes = F, data=mel.lift.ref.chr.pos, aes(xmin=start, xmax=end, ymin=0,ymax=20), colour="transparent", fill=c(rep(c("white", "grey40"), 10 ), "white"), alpha=.1 )+
  geom_hline(yintercept = 4 , colour="grey40", lty="dashed",alpha=.99)+
  geom_point(color="darkgrey", alpha=0.9, size=0.02) +
  geom_point(inherit.aes = T, aes(x=era.BP.wg, y=zPBS0.3),   data=subset(mel.all.pbs, !(is.na(hdr.EcE)) & is.hdr=="hdr"), alpha=0.7, size=.1, colour="grey40") +
  geom_point(inherit.aes = T, aes(x=era.BP.wg, y=zPBS0.3),   data=subset(mel.all.pbs, sharing.type=="within.east"), alpha=0.7, size=.4, colour="#049E73") +
  geom_point(inherit.aes = T, aes(x=era.BP.wg, y=zPBS0.3),   data=subset(mel.all.pbs, sharing.type=="across"), alpha=0.7, size=.4, colour="#D65D00") +
  scale_x_continuous( expand = c(0, 0),label = axisdf.mel.lift$era.CHR, breaks= axisdf.mel.lift$center ) +
  scale_y_continuous( expand = c(0, 0), limits = c(0,20), labels = c(0, 10, ""), breaks = c(0,10,20) ) +    # remove space between plot area and x axis
  ylab("zPBS high")+
  #theme_bw()+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=.5),legend.position="none",
        plot.margin = unit(c(0,0.1,0,0.1), "lines"),axis.text = element_text(size=10), axis.text.x = element_blank(),axis.ticks.x = element_blank(),
        axis.text.y = element_text(size=10),axis.title.x = element_blank(), panel.background = element_rect(fill = "transparent"), plot.background = element_rect(fill = "transparent", color = NA),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.background = element_rect(fill = "transparent"), legend.box.background = element_rect(fill = "transparent")); ec.e.mel.p.lift




#################################### recombination rate plotting ####################################
### plot ###
era.rec.p <- ggplot(data=era.rec.summ.S, aes(x=BP.wg, y=meanS)) +
  # facet_grid(~chromosome, scales = 'free_x', space = 'free_x', switch = 'x')+ # to emulate sunin's plots
  geom_rect(inherit.aes = F, data=era.ref.chr.pos, aes(xmin=start, xmax=end, ymin=0,ymax=7), colour="transparent", fill=c(rep(c("white", "grey40"), 10 ), "white"), alpha=.2 )+
  geom_rect(inherit.aes = F, data=subset(era.shared.regions, sharing.type=="within.west"), aes(xmin=start.BP.wg-100000, xmax=end.BP.wg+100000, ymin=0,ymax=7), colour="transparent", fill="#0372B2", alpha=0.7 )+
  geom_rect(inherit.aes = F, data=subset(era.shared.regions, sharing.type=="within.east"), aes(xmin=start.BP.wg-100000, xmax=end.BP.wg+100000, ymin=0,ymax=7), colour="transparent", fill="#049E73", alpha=0.7 )+
  geom_rect(inherit.aes = F, data=subset(era.shared.regions, sharing.type=="across"), aes(xmin=start.BP.wg-100000, xmax=end.BP.wg+100000, ymin=0,ymax=7), colour="transparent", fill="#D65D00", alpha=1 )+
  geom_line(aes(y=meanS.spline, group=chromosome))+
  #geom_line(aes(y=rollmean(meanS, 10, na.pad=TRUE, ))) +
  #geom_point(color="darkgrey", alpha=0.9, size=0.02) +
  scale_x_continuous( expand = c(0, 0),label = axisdf.era$CHR, breaks= axisdf.era$center ) +
  scale_y_continuous( expand = c(0, 0), limits = c(0,7), labels = c(0, 5, ""), breaks = c(0,5,20) ) +    # remove space between plot area and x axis
  ylab("Rec. rate\n(cM/Mb)")+
  #theme_bw()+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=.5),legend.position="none",
        plot.margin = unit(c(0,0.1,0.5,0.1), "lines"),axis.text = element_text(size=10), axis.text.x = element_blank(),axis.ticks.x = element_blank(),
        axis.text.y = element_text(size=10),axis.title.x = element_blank(), panel.background = element_rect(fill = "transparent"), plot.background = element_rect(fill = "transparent", color = NA),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.background = element_rect(fill = "transparent"), legend.box.background = element_rect(fill = "transparent"));era.rec.p

mel.rec.p <- ggplot(data=mel.rec.summ.S, aes(x=BP.wg, y=meanS)) +
  geom_rect(inherit.aes = F, data=mel.ref.chr.pos, aes(xmin=start, xmax=end, ymin=0,ymax=15), colour="transparent", fill=c(rep(c("white", "grey40"), 10 ), "white"), alpha=.2 )+
  geom_rect(inherit.aes = F, data=subset(mel.shared.regions, sharing.type=="within.west"), aes(xmin=start.BP.wg-100000, xmax=end.BP.wg+100000, ymin=0,ymax=15), colour="transparent", fill="#0372B2", alpha=0.7 )+
  geom_rect(inherit.aes = F, data=subset(mel.shared.regions, sharing.type=="within.east"), aes(xmin=start.BP.wg-100000, xmax=end.BP.wg+100000, ymin=0,ymax=15), colour="transparent", fill="#049E73", alpha=0.7 )+
  geom_rect(inherit.aes = F, data=subset(mel.shared.regions, sharing.type=="across"), aes(xmin=start.BP.wg-100000, xmax=end.BP.wg+100000, ymin=0,ymax=15), colour="transparent", fill="#D65D00", alpha=1 )+
  geom_line(aes(y=meanS.spline, group=chr))+
  #geom_line(aes(y=rollmean(meanS, 10, na.pad=TRUE))) + # best use splines instead of rolling means
  #geom_point(color="darkgrey", alpha=0.9, size=0.02) +
   scale_x_continuous( expand = c(0, 0),label = axisdf.mel$CHR, breaks= axisdf.mel$center ) +
   scale_y_continuous( expand = c(0, 0), limits = c(0,15), labels = c(0, 10, ""), breaks = c(0,10,20) ) +    # remove space between plot area and x axis
  ylab("Rec. rate\n(cM/Mb)")+
  #theme_bw()+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=.5),legend.position="none",
        plot.margin = unit(c(0,0.1,0,0.1), "lines"),axis.text = element_text(size=10), axis.text.x = element_blank(),axis.ticks.x = element_blank(),
        axis.text.y = element_text(size=10),axis.title.x = element_blank(), panel.background = element_rect(fill = "transparent"), plot.background = element_rect(fill = "transparent", color = NA),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.background = element_rect(fill = "transparent"), legend.box.background = element_rect(fill = "transparent"));mel.rec.p



#################################### Fig 2 A- PBS clines use patchwork ###############
# with positions lifted over in melpomene (i.e. missing some windows)
co.w.era.p / ec.w.era.p / co.e.era.p / ec.e.era.p /era.rec.p / co.w.mel.p.lift / ec.w.mel.p.lift / co.e.mel.p.lift / ec.e.mel.p.lift / mel.rec.p
ggsave("~/Dropbox (Cambridge University)/PhD/22_pop.gen.paper/figures/fig2/fig2A.mel.lift.png", width=10, height = 10)

# with diff positions for erato/melpomene
co.w.era.p / ec.w.era.p / co.e.era.p / ec.e.era.p /era.rec.p / co.w.mel.p / ec.w.mel.p / co.e.mel.p / ec.e.mel.p / mel.rec.p
ggsave("~/Dropbox (Cambridge University)/PhD/22_pop.gen.paper/figures/fig2/fig2A.png", width=10, height = 10)


setwd("/Users/gabrielamontejokovacevich/Dropbox/tmp/")
packages <- installed.packages()[,"Package"]
save(packages, file="Rpackages")

#################################### Fig 2 B - recomb rate distribution ###############
head(era.all.pbs)


###### plot distribution of scaled recomb rate, plus lines on top of mean values of recomb rate at HDRs ######
head(era.rec.summ.S)
# less skewed if you remove colour pattern loci &CHR!=18&CHR!=15)
rec.dens.era <-  ggplot(subset(era.rec.summ.S), aes(x=meanS)) + 
  geom_density(fill="grey", colour="grey", alpha=0.99)+
  geom_density(inherit.aes = F, aes(x=r), data=subset(rec.hdr.era.summ.para.west, sharing.type=="within.west"), fill="transparent", colour="#0372B2", alpha=0, size=1)+
  geom_density(inherit.aes = F, aes(x=r), data= subset(rec.hdr.era.summ.para.east, sharing.type=="within.east"), fill="transparent", colour="#049E73", alpha=0, size=1)+
  geom_density(inherit.aes = F, aes(x=r), data=subset(rec.hdr.era.summ.allo, sharing.type=="across"), fill="transparent", colour="#D65D00", alpha=0, size=1)+
  geom_vline(xintercept = mean(era.all.pbs$r), colour="darkgrey", size=1, lty="dashed")+
  geom_vline(xintercept = mean(subset(rec.hdr.era.summ.para.west, sharing.type=="within.west")$r), colour="#0372B2", size=1, lty="solid")+
  geom_vline(xintercept = mean(subset(rec.hdr.era.summ.para.east, sharing.type=="within.east")$r), colour="#049E73", size=1, lty="solid")+
  geom_vline(xintercept = mean(subset(rec.hdr.era.summ.allo, sharing.type=="across")$r), colour="#D65D00", size=1, lty="solid")+
  geom_rect(inherit.aes = F, data=subset(rec.hdr.era.summ.para.west, sharing.type=="within.west"), aes(xmin=r-.005, xmax=r+.005, ymin=-0.05,ymax=0), colour="transparent", fill="#0372B2", alpha=0.7 )+
    geom_rect(inherit.aes = F, data=subset(rec.hdr.era.summ.para.east, sharing.type=="within.east"), aes(xmin=r-.005, xmax=r+.005, ymin=-0.05,ymax=0), colour="transparent", fill="#049E73", alpha=0.7 )+
    geom_rect(inherit.aes = F, data=subset(rec.hdr.era.summ.allo, sharing.type=="across"), aes(xmin=r-.005, xmax=r+.005, ymin=-0.05,ymax=0), colour="transparent", fill="#D65D00", alpha=1 )+
    xlab(expression(paste("Recombination rate cM/Mb",  italic(" H. erato"))))+
  theme_classic()+
  theme(axis.line.x = element_line(color="black", size = 0.5),axis.line.y = element_line(color="black", size = 0.5),axis.title.y = element_text(size=14),
        plot.margin = unit(c(0.5,0.2,0.2,0.2), "lines"),axis.text = element_text(size=10), #axis.text.x = element_blank(),
        axis.text.y = element_text(size=10),axis.title.x = element_blank(), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_rect(fill = "transparent"), plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
        legend.background = element_rect(fill = "transparent"), legend.box.background = element_rect(fill = "transparent"))+
  scale_y_continuous(limits = c(-0.05,0.8), expand = c(0,0))+scale_x_continuous(limits = c(0,7), expand = c(0,0)); rec.dens.era

rec.dens.mel <-  ggplot(mel.rec.summ.S, aes(x=meanS)) + 
  geom_density(fill="grey", colour="grey", alpha=0.99)+
  geom_density(inherit.aes = F, aes(x=r), data=subset(rec.hdr.mel.summ.para.west, sharing.type=="within.west"), fill="transparent", colour="#0372B2", alpha=0, size=1)+
  geom_density(inherit.aes = F, aes(x=r), data= subset(rec.hdr.mel.summ.para.east, sharing.type=="within.east"), fill="transparent", colour="#049E73", alpha=0, size=1)+
  geom_density(inherit.aes = F, aes(x=r), data=subset(rec.hdr.mel.summ.allo, sharing.type=="across"), fill="transparent", colour="#D65D00", alpha=0, size=1)+
  geom_vline(xintercept = mean(mel.all.pbs$r), colour="darkgrey", size=1, lty="dashed")+
  geom_vline(xintercept = mean(subset(rec.hdr.mel.summ.para.west, sharing.type=="within.west")$r), colour="#0372B2", size=1, lty="solid")+
  geom_vline(xintercept = mean(subset(rec.hdr.mel.summ.para.east, sharing.type=="within.east")$r), colour="#049E73", size=1, lty="solid")+
  geom_vline(xintercept = mean(subset(rec.hdr.mel.summ.allo, sharing.type=="across")$r), colour="#D65D00", size=1, lty="solid")+
  #geom_rect(inherit.aes = F, data=rec.hdr.mel.summ, aes(xmin=r-.005, xmax=r+.005, ymin=0,ymax=0.5), colour="transparent", fill="red", alpha=1 )+
  geom_rect(inherit.aes = F, data=subset(rec.hdr.mel.summ.para.west, sharing.type=="within.west"), aes(xmin=r-.005, xmax=r+.005, ymin=-0.05,ymax=0), colour="transparent", fill="#0372B2", alpha=0.7 )+
  geom_rect(inherit.aes = F, data=subset(rec.hdr.mel.summ.para.east, sharing.type=="within.east"), aes(xmin=r-.005, xmax=r+.005, ymin=-0.05,ymax=0), colour="transparent", fill="#049E73", alpha=0.7 )+
  geom_rect(inherit.aes = F, data=subset(rec.hdr.mel.summ.allo, sharing.type=="across"), aes(xmin=r-.005, xmax=r+.005, ymin=-0.05,ymax=0), colour="transparent", fill="#D65D00", alpha=1 )+
  xlab(expression(paste("Recombination rate cM/Mb",  italic(" H. melpomene"))))+
  theme_classic()+
  theme(axis.line.x = element_line(color="black", size = 0.5),axis.line.y = element_line(color="black", size = 0.5),axis.title.y = element_text(size=14),
        plot.margin = unit(c(0.5,0.2,0.2,0.2), "lines"),axis.text = element_text(size=10), #axis.text.x = element_blank(),
        axis.text.y = element_text(size=10),axis.title.x = element_blank(), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_rect(fill = "transparent"), plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
        legend.background = element_rect(fill = "transparent"), legend.box.background = element_rect(fill = "transparent"))+
  scale_y_continuous(limits = c(-0.05,0.6), expand = c(0,0)) +scale_x_continuous(limits = c(0,15.5), expand = c(0,0)); rec.dens.mel


# save
rec.dens.era /rec.dens.mel
ggsave("~/Dropbox (Cambridge University)/PhD/22_pop.gen.paper/figures/fig2/fig2B.png", width=5, height = 5)

rec.dens.era
ggsave("~/Dropbox (Cambridge University)/PhD/22_pop.gen.paper/figures/fig2/fig2Bera.png", width=3, height = 2.5)

rec.dens.mel
ggsave("~/Dropbox (Cambridge University)/PhD/22_pop.gen.paper/figures/fig2/fig2Bmel.png", width=3, height = 2.5)



#################################### Fig 2 C - sharing sims densities ####################################
density.10ksim.50kb <- read.csv("local/data/sharing/density.sharing.percs/density.10ksim.50kb.csv") # simulation densities
sharing.percs.50kb <- read.csv("local/data/sharing/density.sharing.percs/sharing.percs.50kb.csv") # percentages
names(density.10ksim.50kb)

#### across sides data prep ####
# to obtain perc of total, multiply %shared within(e&w)*across, then average (for both sims and observed)
# and average densities in within(e&w) and across
density.10ksim.50kb.era.e.within <- subset(density.10ksim.50kb, type.comp=="within"&species=="era"&type.side=="e")
density.10ksim.50kb.era.w.within <- subset(density.10ksim.50kb, type.comp=="within"&species=="era"&type.side=="w")
density.10ksim.50kb.era.across <- subset(density.10ksim.50kb, type.comp=="across"&species=="era")
density.10ksim.50kb.era.across.of.total <- data.frame(perc=(((density.10ksim.50kb.era.e.within$perc/100)*(density.10ksim.50kb.era.across$perc/100)*100)+
                                                              (density.10ksim.50kb.era.w.within$perc/100)*(density.10ksim.50kb.era.across$perc/100)*100)/2,
                                                      density=((((density.10ksim.50kb.era.e.within$density)+(density.10ksim.50kb.era.across$density))/2)+
                                                                 ((density.10ksim.50kb.era.w.within$density)+(density.10ksim.50kb.era.across$density))/2)/2);density.10ksim.50kb.era.across.of.total

sharing.percs.50kb.era.e.within <- subset(sharing.percs.50kb, type.comp=="within"&species=="era"&type.side=="e")
sharing.percs.50kb.era.w.within <- subset(sharing.percs.50kb, type.comp=="within"&species=="era"&type.side=="w")
sharing.percs.50kb.era.across <- subset(sharing.percs.50kb, type.comp=="across"&species=="era")

sharing.percs.50kb.era.across.of.total <- data.frame(overlap.int.prop=(((sharing.percs.50kb.era.e.within$overlap.int.prop/100)*(sharing.percs.50kb.era.across$overlap.int.prop/100)*100)+
                                                                         (sharing.percs.50kb.era.w.within$overlap.int.prop/100)*(sharing.percs.50kb.era.across$overlap.int.prop/100)*100)/2,
                                                     overlap.prop.pseudo.std_err=((((sharing.percs.50kb.era.e.within$overlap.prop.pseudo.std_err)+(sharing.percs.50kb.era.across$overlap.prop.pseudo.std_err))/2)+
                                                                                    ((sharing.percs.50kb.era.w.within$overlap.prop.pseudo.std_err)+(sharing.percs.50kb.era.across$overlap.prop.pseudo.std_err))/2)/2);sharing.percs.50kb.era.across.of.total
# and avmelge densities in within(e&w) and across
density.10ksim.50kb.mel.e.within <- subset(density.10ksim.50kb, type.comp=="within"&species=="mel"&type.side=="e")
density.10ksim.50kb.mel.w.within <- subset(density.10ksim.50kb, type.comp=="within"&species=="mel"&type.side=="w")
density.10ksim.50kb.mel.across <- subset(density.10ksim.50kb, type.comp=="across"&species=="mel")
density.10ksim.50kb.mel.across.of.total <- data.frame(perc=(((density.10ksim.50kb.mel.e.within$perc/100)*(density.10ksim.50kb.mel.across$perc/100)*100)+
                                                              (density.10ksim.50kb.mel.w.within$perc/100)*(density.10ksim.50kb.mel.across$perc/100)*100)/2,
                                                      density=((((density.10ksim.50kb.mel.e.within$density)+(density.10ksim.50kb.mel.across$density))/2)+
                                                                 ((density.10ksim.50kb.mel.w.within$density)+(density.10ksim.50kb.mel.across$density))/2)/2);density.10ksim.50kb.mel.across.of.total

sharing.percs.50kb.mel.e.within <- subset(sharing.percs.50kb, type.comp=="within"&species=="mel"&type.side=="e")
sharing.percs.50kb.mel.w.within <- subset(sharing.percs.50kb, type.comp=="within"&species=="mel"&type.side=="w")
sharing.percs.50kb.mel.across <- subset(sharing.percs.50kb, type.comp=="across"&species=="mel")

sharing.percs.50kb.mel.across.of.total <- data.frame(overlap.int.prop=(((sharing.percs.50kb.mel.e.within$overlap.int.prop/100)*(sharing.percs.50kb.mel.across$overlap.int.prop/100)*100)+
                                                                         (sharing.percs.50kb.mel.w.within$overlap.int.prop/100)*(sharing.percs.50kb.mel.across$overlap.int.prop/100)*100)/2,
                                                     overlap.prop.pseudo.std_err=((((sharing.percs.50kb.mel.e.within$overlap.prop.pseudo.std_err)+(sharing.percs.50kb.mel.across$overlap.prop.pseudo.std_err))/2)+
                                                                                    ((sharing.percs.50kb.mel.w.within$overlap.prop.pseudo.std_err)+(sharing.percs.50kb.mel.across$overlap.prop.pseudo.std_err))/2)/2);sharing.percs.50kb.mel.across.of.total

#### 6 plots separately  ####

# era.w.within
era.w.within <- ggplot(data = density.10ksim.50kb.era.w.within, aes(x=perc, y=density)) + 
  geom_line() +
  geom_ribbon(aes(ymin=0, ymax=pmax(density,0)), fill="grey", col="#0372B2", alpha=0.5)+
  scale_x_continuous(limits = c(0,55), expand = c(0, 0)) +
  scale_y_continuous( expand = c(0, 0)) +
  theme_classic()+
  theme(axis.line.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.title.y = element_blank(), 
        axis.title.x = element_blank(), axis.text.x = element_text(size=14))+
  geom_vline(aes(xintercept =sharing.percs.50kb.era.w.within$overlap.int.prop), col="#0372B2", lwd=2)+
  geom_vline(aes(xintercept =sharing.percs.50kb.era.w.within$overlap.int.prop+sharing.percs.50kb.era.w.within$overlap.prop.pseudo.std_err*1.96), 
             lty = 2, col="#0372B2", lwd=.7)+
  geom_vline(aes(xintercept =sharing.percs.50kb.era.w.within$overlap.int.prop-sharing.percs.50kb.era.w.within$overlap.prop.pseudo.std_err*1.96), 
             lty = 2, col="#0372B2", lwd=.7)+ 
  theme(
    panel.background = element_rect(fill = "transparent"), # bg of the panel
    plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
    panel.grid.major = element_blank(), # get rid of major grid
    panel.grid.minor = element_blank(), # get rid of minor grid
    legend.background = element_rect(fill = "transparent"), # get rid of legend bg
    legend.box.background = element_rect(fill = "transparent") # get rid of legend panel bg
  ); era.w.within

# era.e.within
era.e.within <- ggplot(data = density.10ksim.50kb.era.e.within, aes(x=perc, y=density)) + 
  geom_line() +
  geom_ribbon(aes(ymin=0, ymax=pmax(density,0)), fill="grey", col="#049E73", alpha=0.5)+
  scale_x_continuous(limits = c(0,55), expand = c(0, 0)) +
  scale_y_continuous( expand = c(0, 0)) +
  theme_classic()+
  theme(axis.line.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.title.y = element_blank(), 
        axis.title.x = element_blank(), axis.text.x = element_text(size=14))+
  geom_vline(aes(xintercept =sharing.percs.50kb.era.e.within$overlap.int.prop), col="#049E73", lwd=2)+
  geom_vline(aes(xintercept =sharing.percs.50kb.era.e.within$overlap.int.prop+sharing.percs.50kb.era.e.within$overlap.prop.pseudo.std_err*1.96),
             lty = 2, col="#049E73", lwd=.7)+
  geom_vline(aes(xintercept =sharing.percs.50kb.era.e.within$overlap.int.prop-sharing.percs.50kb.era.e.within$overlap.prop.pseudo.std_err*1.96),
             lty = 2,  col="#049E73", lwd=.7)+ 
  theme(
    panel.background = element_rect(fill = "transparent"), # bg of the panel
    plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
    panel.grid.major = element_blank(), # get rid of major grid
    panel.grid.minor = element_blank(), # get rid of minor grid
    legend.background = element_rect(fill = "transparent"), # get rid of legend bg
    legend.box.background = element_rect(fill = "transparent") # get rid of legend panel bg
  ); era.e.within

# mel.w.within
mel.w.within <- ggplot(data = density.10ksim.50kb.mel.w.within, aes(x=perc, y=density)) + 
  geom_line() +
  geom_ribbon(aes(ymin=0, ymax=pmax(density,0)), fill="grey", col="#0372B2", alpha=0.5)+
  scale_x_continuous(limits = c(0,55), expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme_classic()+
  theme(axis.line.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.title.y = element_blank(), 
        axis.title.x = element_blank(), axis.text.x = element_text(size=14))+
  geom_vline(aes(xintercept =sharing.percs.50kb.mel.w.within$overlap.int.prop), col="#0372B2", lwd=2)+
  geom_vline(aes(xintercept =sharing.percs.50kb.mel.w.within$overlap.int.prop+sharing.percs.50kb.mel.w.within$overlap.prop.pseudo.std_err*1.96),
             lty = 2, col="#0372B2", lwd=.7)+
  geom_vline(aes(xintercept =sharing.percs.50kb.mel.w.within$overlap.int.prop-sharing.percs.50kb.mel.w.within$overlap.prop.pseudo.std_err*1.96), 
             lty = 2, col="#0372B2", lwd=.7)+ 
  theme(
    panel.background = element_rect(fill = "transparent"), # bg of the panel
    plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
    panel.grid.major = element_blank(), # get rid of major grid
    panel.grid.minor = element_blank(), # get rid of minor grid
    legend.background = element_rect(fill = "transparent"), # get rid of legend bg
    legend.box.background = element_rect(fill = "transparent") # get rid of legend panel bg
  ); mel.w.within



# mel.e.within
mel.e.within <- ggplot(data = density.10ksim.50kb.mel.e.within, aes(x=perc, y=density)) + 
  geom_line() +
  geom_ribbon(aes(ymin=0, ymax=pmax(density,0)), fill="grey", col="#049E73", alpha=0.5)+
  scale_x_continuous(limits = c(0,55), expand = c(0, 0)) +
  scale_y_continuous( expand = c(0, 0)) +
  theme_classic()+
  theme(axis.line.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.title.y = element_blank(), 
        axis.title.x = element_blank(), axis.text.x = element_text(size=14))+
  geom_vline(aes(xintercept =sharing.percs.50kb.mel.e.within$overlap.int.prop), col="#049E73", lwd=2)+
  geom_vline(aes(xintercept =sharing.percs.50kb.mel.e.within$overlap.int.prop+sharing.percs.50kb.mel.e.within$overlap.prop.pseudo.std_err*1.96), 
             lty = 2, col="#049E73", lwd=.7)+
  geom_vline(aes(xintercept =sharing.percs.50kb.mel.e.within$overlap.int.prop-sharing.percs.50kb.mel.e.within$overlap.prop.pseudo.std_err*1.96), 
             lty = 2, col="#049E73", lwd=.7)+ 
  theme(
    panel.background = element_rect(fill = "transparent"), # bg of the panel
    plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
    panel.grid.major = element_blank(), # get rid of major grid
    panel.grid.minor = element_blank(), # get rid of minor grid
    legend.background = element_rect(fill = "transparent"), # get rid of legend bg
    legend.box.background = element_rect(fill = "transparent") # get rid of legend panel bg
  ); mel.e.within

# era.across

era.across<- ggplot(data = density.10ksim.50kb.era.across.of.total, aes(x=perc, y=density)) + 
  geom_line() +
  geom_ribbon(aes(ymin=0, ymax=pmax(density,0)), fill="grey", col="#D65D00", alpha=0.5)+
  scale_x_continuous(limits = c(0,55), expand = c(0, 0)) +
  scale_y_continuous( expand = c(0, 0)) +
  theme_classic()+
  theme(axis.line.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.title.y = element_blank(), 
        axis.title.x = element_blank(), axis.text.x = element_text(size=14))+
  geom_vline(aes(xintercept =sharing.percs.50kb.era.across.of.total$overlap.int.prop), color="#D65D00", lwd=2)+
  geom_vline(aes(xintercept =sharing.percs.50kb.era.across.of.total$overlap.int.prop + sharing.percs.50kb.era.across.of.total$overlap.prop.pseudo.std_err*1.96),
             lty = 2, color="#D65D00", lwd=.7) +
  geom_vline(aes(xintercept =sharing.percs.50kb.era.across.of.total$overlap.int.prop - sharing.percs.50kb.era.across.of.total$overlap.prop.pseudo.std_err*1.96), 
             lty = 2, color="#D65D00", lwd=.7)+ 
  theme(
    panel.background = element_rect(fill = "transparent"), # bg of the panel
    plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
    panel.grid.major = element_blank(), # get rid of major grid
    panel.grid.minor = element_blank(), # get rid of minor grid
    legend.background = element_rect(fill = "transparent"), # get rid of legend bg
    legend.box.background = element_rect(fill = "transparent") # get rid of legend panel bg
  ) ; era.across


# mel.across
mel.across<- ggplot(data = density.10ksim.50kb.mel.across.of.total, aes(x=perc, y=density)) + 
  geom_line() +
  geom_ribbon(aes(ymin=0, ymax=pmax(density,0)), fill="grey", color="#D65D00", alpha=0.5)+
  scale_x_continuous(limits = c(0,55), expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme_classic()+
  theme(axis.line.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.title.y = element_blank(), 
        axis.title.x = element_blank(), axis.text.x = element_text(size=14))+
  geom_vline(aes(xintercept =sharing.percs.50kb.mel.across.of.total$overlap.int.prop),  color="#D65D00", lwd=2)+
  geom_vline(aes(xintercept =sharing.percs.50kb.mel.across.of.total$overlap.int.prop + sharing.percs.50kb.mel.across.of.total$overlap.prop.pseudo.std_err*1.96), 
             lty = 2, color="#D65D00", lwd=.7) +
  geom_vline(aes(xintercept =sharing.percs.50kb.mel.across.of.total$overlap.int.prop - sharing.percs.50kb.mel.across.of.total$overlap.prop.pseudo.std_err*1.96), 
             lty = 2, color="#D65D00", lwd=.7) + 
  theme(
    panel.background = element_rect(fill = "transparent"), # bg of the panel
    plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
    panel.grid.major = element_blank(), # get rid of major grid
    panel.grid.minor = element_blank(), # get rid of minor grid
    legend.background = element_rect(fill = "transparent"), # get rid of legend bg
    legend.box.background = element_rect(fill = "transparent") # get rid of legend panel bg
  ); mel.across


sharing.all.50kb <- plot_grid(era.w.within, era.across, era.e.within, mel.w.within, mel.across, mel.e.within, nrow = 2, scale=.95); sharing.all.50kb
ggsave("local/plots/fig.2/b.sharing.sims.xlim.50kb.png",  bg = "transparent" )



####################################  SI Fig - block size vs % genome covered ###########
tables <- list()
tables.listpaths <- dir(pattern = ".*block.summ.*txt", path = "local/plots/sharing.sims", full.names = TRUE); tables.listpaths
for (i in 1:8) {
  tables[[i]] <- read.table(tables.listpaths[i], header = TRUE)
  #tables[[i]] <- tables[[i]][,c(1:9)]
  tables[[i]]$buffer.kb <- as.numeric(as.character(gsub("([0-9]+).*$", "\\1", word(tables.listpaths[i],9,sep = fixed(".")))))
  tables[[i]]$species <- word(tables.listpaths[i],8,sep = fixed("."))
}

tables.df <-rbindlist(tables)
tables.df
ggplot(aes(x=mean.block.size.kb, y=perc.genome.cov.block, colour=pop), data=tables.df)+
  geom_point()+geom_line()+
  facet_grid(~species)
ggsave("local/plots/sharing.sims/perc.genome.cov.vs.block.size.png")

ggplot(aes(x=mean.block.size.kb, y=n.block, colour=pop), data=tables.df)+
  geom_point()+geom_line()+
  facet_grid(~species)
ggsave("local/plots/sharing.sims/n.block.vs.block.size.png")

ggplot(aes(x=buffer.kb , y=mean.block.size.kb-(buffer.kb*2), colour=pop), data=tables.df)+
  geom_point()+geom_line()+
  facet_grid(~species)
ggsave("local/plots/sharing.sims/n.block.vs.block.size.png")

## sharing plots
ggplot(aes(x=buffer.kb , y=overlap.int.prop, colour=pop), data=tables.df)+
  geom_point()+geom_line()+
  facet_grid(~species)
ggsave("local/plots/sharing.sims/overlap.int.prop.vs.buffer.kb.png")

ggplot(aes(x=mean.block.size.kb , y=overlap.int.prop, colour=pop), data=tables.df)+
  geom_point()+geom_line()+
  facet_grid(~species)
ggsave("local/plots/sharing.sims/overlap.int.prop.vs.mean.block.size.kb.png")

ggplot(aes(x=perc.genome.cov.block , y=overlap.int.prop, colour=pop, shape=pop), data=tables.df)+
  geom_point()+geom_line()+
  geom_line(aes(y=sim10k.overlap.prop, x=perc.genome.cov.block),linetype = "dashed")+
  geom_point(aes(y=sim10k.overlap.prop, x=perc.genome.cov.block))+
  facet_grid(~species)
ggsave("local/plots/sharing.sims/overlap.int.prop.and.10ksim.vs.perc.genome.cov.block.png")

