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
library(RcppCNPy)
library(ggfortify)
library(ggrepel)
library(zoo)
library(ggExtra)
options(scipen = 999)
#devtools::install_github("petrelharp/local_pca/lostruct")
library(lostruct)
library(ggpubr)
#install_github("kassambara/factoextra")
library("factoextra")
library(data.table)
library(relaimpo)


##### functions ##### 
setwd("/Users/gabrielamontejokovacevich/Dropbox (Cambridge University)/PhD/7_Assocation_studies/9_ANGSD/12.era.mel.altitude/local/")

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



# ### regions ###
# era.e <- read.csv("local/data/sharing/shared.e.regions.era.75.csv")
# head(era.e)
# 
# unique(era.e$scaff)
# ggplot(aes(x=BP.wg, y=zPBS0.1), data=subset(era.e, scaff=="Herato0209"))+
#   geom_point()


# Function that removes windows with too low number of ABBA and BABA counts
# sets fd to 0 if D is negative
corrNegfd<-function(dataset){
  dataset<-dataset[!is.na(dataset$fd),]
  dataset<-dataset[dataset$ABBA+dataset$BABA>1,]
  dataset[dataset$D<0 & !is.na(dataset$D),"fd"]<-0
  dataset[dataset$fd<0 & !is.na(dataset$fd),"fd"]<-0
  return(dataset)
}

# Reformat and correct Negative fd

getFd.era<-function(combination){
  prefix="data/joana.fd/output/eratoRelatives.withSRA.chr1-21.max0.5N.minDP3.biallSNPs.mac2."
  
  # Read in the fd.df
  fd.df<-read.csv(paste0(prefix,combination,"_50kb"))
  
  # Correct negative values
  fd.df<-corrNegfd(fd.df)
  
  # reformat to datatable
  fd.df<-data.table(fd.df)
  
  # Add chromosome name
  fd.df$chr<-substring(fd.df$scaffold,first = 7,last=8)
  
  # add window information
  fd.df<-data.table(fd.df)
  fd.df$start.WG <-fd.df$start+ref.scaff.era[match(fd.df$scaffold,ref.scaff.era$scaffold),"add"]
  fd.df$end.WG <-fd.df$end+ref.scaff.era[match(fd.df$scaffold,ref.scaff.era$scaffold),"add"]
  fd.df$mid.WG <-(fd.df$end.WG +  fd.df$start.WG)/2
  
  # add if east (only) PBS outlier (para SHDR)
  shdr.east.ranges<-data.table(era.outliersEast[,c("scaff","start.bp","end.bp")])
  names(shdr.east.ranges)<-c("scaffold","start","end") # to match fd.df names
  setkey(shdr.east.ranges)
  # by.x concatenates scaffold and start/end
  oveast<-foverlaps(x = fd.df,y = shdr.east.ranges, which=T, by.x = names(shdr.east.ranges), type="any",nomatch=0L)
  fd.df$shdr.para.east.id<- NA
  # use keys to select corresponding shdr names
  fd.df$shdr.para.east.id[oveast$xid]<- paste(era.outliersEast[ oveast$yid,]$shdr.para.east.id)
  
  # add if west (only) PBS outlier (para SHDR)
  shdr.west.ranges<-data.table(era.outliersWest[,c("scaff","start.bp","end.bp")])
  names(shdr.west.ranges)<-c("scaffold","start","end") # to match fd.df names
  setkey(shdr.west.ranges)
  # by.x concatenates scaffold and start/end
  ovWest<-foverlaps(x = fd.df,y = shdr.west.ranges, which=T, by.x = names(shdr.west.ranges), type="any",nomatch=0L)
  fd.df$shdr.para.west.id<- NA
  # use keys to select corresponding shdr names
  fd.df$shdr.para.west.id[ovWest$xid]<-  paste(era.outliersWest[ ovWest$yid,]$shdr.para.west.id)
  
  # add if across (only) PBS outlier (para SHDR)
  shdr.across.ranges<-data.table(era.outliersAcross[,c("scaff","start.bp","end.bp")])
  names(shdr.across.ranges)<-c("scaffold","start","end") # to match fd.df names
  setkey(shdr.across.ranges)
  # by.x concatenates scaffold and start/end
  ovacross<-foverlaps(x = fd.df,y = shdr.across.ranges, which=T, by.x = names(shdr.across.ranges), type="any",nomatch=0L)
  fd.df$shdr.allo.id<- NA
  # use keys to select corresponding shdr names
  fd.df$shdr.allo.id[ovacross$xid]<-  paste(era.outliersAcross[ ovacross$yid,]$shdr.allo.id)
  
  
  return(fd.df)
}

# personalised fill layer https://stackoverflow.com/questions/43440068/ggplot2-fix-colors-to-factor-levels
scale_color_shdr.type <- function(...){
  ggplot2:::manual_scale(
    'color', 
    values = setNames(c("#D45D00", "#0372B2","#049E73" ), c("shdr.allo", "shdr.west", "shdr.east")), 
    ...
  )
}

# barplot(c(2,5,10), col=c("#D45D00", "#0372B2","#049E73"))

#################################### 0. data prep #####
setwd("/Users/gabrielamontejokovacevich/Dropbox (Cambridge University)/PhD/7_Assocation_studies/9_ANGSD/12.era.mel.altitude/local/")
ref.scaff.era <- read.table("local/data/local/data/0.ref/Heliconius_erato_demophoon_v1_-_scaffolds.fa.fai", row.names = NULL)
all.bams.list.info <- read.csv("local/data/02.info/all.bams.list.info.csv"); names(all.bams.list.info)
era.e.list.info <- subset(all.bams.list.info,to.use.pbs=="yes" &species=="erato"&side.short=="e" ); nrow(era.e.list.info)
era.w.list.info <- subset(all.bams.list.info,to.use.pbs=="yes" &species=="erato"&side.short=="w" ); nrow(era.w.list.info)
era.a.list.info <- subset(all.bams.list.info,to.use.pbs=="yes" &species=="erato" ); nrow(era.a.list.info)

mel.e.list.info <- subset(all.bams.list.info,to.use.pbs=="yes" &species=="melpomene"&side.short=="e" )
mel.w.list.info <- subset(all.bams.list.info,to.use.pbs=="yes" &species=="melpomene"&side.short=="w" )
mel.a.list.info <- subset(all.bams.list.info,to.use.pbs=="yes" &species=="melpomene" )

#era.all.pbs <- read.csv("local/data/pbs.out/era.all.pbs.recomb.csv")
era.all.pbs$CHR <- as.numeric(as.character(substr(era.all.pbs$scaff,7,8))); unique(era.all.pbs$CHR)
era.all.pbs$Fst02.1 <- if_else(era.all.pbs$Fst02.1<0, 0, era.all.pbs$Fst02.1)
era.all.pbs$Fst02.3 <- if_else(era.all.pbs$Fst02.3<0, 0, era.all.pbs$Fst02.3)
era.all.pbs$zFst02.east.col <- zscore(era.all.pbs$Fst02.1)
era.all.pbs[era.all.pbs$Fst02.3!="" & !(is.na(era.all.pbs$Fst02.3)),]$zFst02.east.ecu <- zscore(era.all.pbs[era.all.pbs$Fst02.3!=""& !(is.na(era.all.pbs$Fst02.3)),]$Fst02.3); era.all.pbs$zFst02.east.ecu
era.all.pbs$Fst02.east.col <- era.all.pbs$Fst02.1
era.all.pbs$Fst01.east.col <- era.all.pbs$Fst01.1
era.all.pbs$Fst12.east.col <- era.all.pbs$Fst12.1

era.all.pbs$Fst02.east.ecu <- era.all.pbs$Fst02.3
era.all.pbs$Fst01.east.ecu <- era.all.pbs$Fst01.3
era.all.pbs$Fst12.east.ecu <- era.all.pbs$Fst12.3

era.all.pbs$zPBS0.east.col <- era.all.pbs$zPBS0.1
era.all.pbs$zPBS0.east.ecu <- era.all.pbs$zPBS0.3


era.chr2.pbs <-  subset(era.all.pbs, CHR=="2") # we only need chr 2

# read in outlier shdr datasets
era.e.shdr <- read.csv("data/shdr.summ/shdr.era.east.outlier.df.csv"); names(era.e.shdr.outlier)
era.w.shdr <- read.csv("data/shdr.summ/shdr.era.west.outlier.df.csv"); head(era.w.shdr.outlier)

era.e.shdr$chr <- as.numeric(as.character(substr(era.e.shdr$scaff, 7,8)))
era.w.shdr$chr <- as.numeric(as.character(substr(era.w.shdr$scaff, 7,8)))


## recomb rates
# load cyrbia and lativitta chr2
era.rec.summ.S <-  read.csv("data/svb.shm/recomb/era.rec.summ.scaled.csv")
era.rec.cyr.S <-  read.csv("data/svb.shm/recomb/era.rec.cyr.scaled.csv")
era.rec.lat.S <-  read.csv("data/svb.shm/recomb/era.rec.lat.scaled.csv")
era.rec.emm.S <-  read.csv("data/svb.shm/recomb/era.rec.emm.scaled.csv")
era.rec.fav.S <-  read.csv("data/svb.shm/recomb/era.rec.fav.scaled.csv")
era.rec.not.S <-  read.csv("data/svb.shm/recomb/era.rec.not.scaled.csv")


############### prep gff ######
# load gff, get genes, add offsets, get intervals
ref.scaff.era <- read.table("local/data/local/data/local/data/local/data/7_Assocation_studies/9_ANGSD/0.ref/Heliconius_erato_demophoon_v1_-_scaffolds.fa.fai", row.names = NULL)
era.scafEnds <- cumsum(ref.scaff.era[,2]); era.offset <- era.scafEnds - ref.scaff.era[,2]; names(era.offset) <- ref.scaff.era[,1]
era.gff <- read.table(pipe("cut -f 1,3,4,5 local/data/local/data/local/data/local/data/7_Assocation_studies/9_ANGSD/0.ref/heliconius_erato_demophoon_v1_core_32_85_1.gff"),
                      sep = "\t", as.is = T, col.names = c("scaffold", "feature", "start", "end"))
era.gff_gene <- subset(era.gff, feature == "gene")
era.gff_gene$start.offset <- era.gff_gene$start + era.offset[era.gff_gene$scaffold]; era.gff_gene$end.offset <- era.gff_gene$end + era.offset[era.gff_gene$scaffold]
era.gff_gene$chr <- as.integer(substr(era.gff_gene$scaffold, 7,8))



############### selection stats outliers tajima delta pi dxy ###############
era.e.shdr.outlier <- read.csv("data/shdr.summ/shdr.era.east.outlier.df.csv"); names(era.e.shdr.outlier)
era.w.shdr.outlier <- read.csv("data/shdr.summ/shdr.era.west.outlier.df.csv"); head(era.w.shdr.outlier)
mel.e.shdr.outlier <- read.csv("data/shdr.summ/shdr.mel.east.outlier.df.csv"); head(mel.e.shdr.outlier)
mel.w.shdr.outlier <- read.csv("data/shdr.summ/shdr.mel.west.outlier.df.csv"); head(mel.w.shdr.outlier)

#################################### 1. FIG5 PARTS ####################################
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
era.ref.chr.pos <- ref.scaff.era; era.ref.chr.pos$start <- era.ref.chr.pos$add
era.ref.chr.pos$end <-  era.ref.chr.pos$start +era.ref.chr.pos$length; head(era.ref.chr.pos)
era.ref.chr.pos <- summarise(group_by(era.ref.chr.pos, scaff),start=min(start),end=max(end)); era.ref.chr.pos # to highlight scaffolds within chr2

era.ref.chr.pos$chr <- as.integer(substr(era.ref.chr.pos$scaff,7,8))
##################### PANEL A 1. example pipeline east shdr.east.77 pc1 vs altitude ##############
####### zPBS #######
names(era.all.pbs)
id <- c("shdr.east.077")
min(subset(era.all.pbs, shdr.para.east.id==id &zPBS0.3>4)$midPos)
max(subset(era.all.pbs, shdr.para.east.id==id &zPBS0.3>4)$midPos)

zpbs.shdr77 <- ggplot(subset(era.all.pbs, shdr.para.east.id==id), aes(x=midPos/10^6, y=zPBS0.1)) +
  # geom_rect(inherit.aes = F, data=subset(era.e.shdr, shdr.para.east.id==id ),
  #           aes(xmin=start, xmax=end,ymin=-Inf,ymax=Inf, fill=is.allo.hdr), colour="transparent", alpha=.4 )+
  # scale_fill_manual(values=c("#D45D00"))+
  annotate("rect", xmin = 9351500/10^6, xmax = 9576500/10^6, ymin = 4, ymax = 40, colour = "transparent", fill="#EEEEEE", size=1) +
  geom_hline(yintercept = 4) +
  geom_line(inherit.aes = F, color="black", alpha=0.9, size=1, data=subset(era.all.pbs, shdr.para.east.id==id), aes(x=midPos/10^6, y=rollmean(zPBS0.1,1, na.pad=TRUE)), lty="solid") +
  geom_line(inherit.aes = F, color="grey40", alpha=0.9, size=1, data=subset(era.all.pbs,  shdr.para.east.id==id), aes(x=midPos/10^6, y=rollmean(zPBS0.3, 1, na.pad=TRUE)), lty=c("11")) +
  scale_y_continuous(expand = c(0, 2))+
  scale_x_continuous(expand = c(0, 0))+
  ylab("East zPBS high")+ xlab("Position (MB), Herato2101")+ 
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1),#axis.ticks.x=element_blank(),
        axis.ticks.length =unit(-0.15, "cm"), #axis.ticks.margin=unit(0.5, "cm"),
        axis.title = element_text(size=16),axis.text.y = element_text(size=10),axis.text.x = element_text(size=10),
        #axis.text.y = element_blank(),
        #axis.text.x = element_blank(), axis.title = element_blank(), axis.ticks=element_blank(),
        #axis.text.y.right = element_text(size = 10, margin = unit(c(t = 0, r = 0, b = 0, l = -8), "mm"), colour="springgreen", face="bold"),
        panel.background = element_blank(), #legend.text = element_text(size=8), legend.title =  element_text(size=12),
        panel.grid.major = element_blank(),# axis.text.x=element_blank(), #axis.text.y=element_blank(), #
        plot.margin = unit(c(0.1, 0.5, 0.1,0.5), "cm"),
        panel.grid.minor = element_blank(),# axis.title.x=element_blank(),#axis.title.y=element_blank(),
        legend.text.align = 0, legend.position = "none",strip.text = element_blank(),rect = element_rect(fill = "transparent"),plot.background = element_blank(),
        strip.background =element_rect(fill="transparent", color="transparent")); zpbs.shdr77 


####### local PCA #######
### indiv data ####
era.e.bam_names<-read.table("local/data/02.info/pop.list/era.e.txt",header=F)
era.e.bam_names$V1<-substr(era.e.bam_names$V1, 76,84); era.e.bam_names$V1
era.e.bam_names$V2<-gsub('[.rm]', '', era.e.bam_names$V1); era.e.bam_names$V2
era.e.bam_names$V2<- if_else(substr(era.e.bam_names$V2,0,2)=="CS" , gsub('[_]', '', era.e.bam_names$V2), era.e.bam_names$V2); era.e.bam_names$V2
era.e.bam_names$V2[122] <-  "BC0411_ena"
era.e.bam_names$V2

# important !!!!!!!!!!!!!!!!!!!@!!!!!! - order 
era.e.list.info <- era.e.list.info[match(era.e.bam_names$V2, era.e.list.info$id),]
era.e.bam_names$V2 ==era.e.list.info$id

### pca plot ####
era.e.77.npy <-npyLoad(file="output/mds.shared/era.e.stdv4/era.e.region.shdr.east.077.cov.npy")
names(era.e.77.npy) <- "shdr.east.077"
region <-"shdr.east.077"
era.e.77.pca <- region
era.e.77.pca<-prcomp(era.e.77.npy)
rownames(era.e.77.pca$x)<-era.e.bam_names$V2

# store pca info for plotting/analyses
era.e.77.pca.df <- as.data.frame(era.e.77.pca$x[,c(1:2)]); era.e.77.pca.df$id <- rownames(era.e.77.pca.df)
era.e.77.pca.df$altitude <- subset(era.e.list.info, id %in% era.e.bam_names$V2)$altitude[match(era.e.77.pca.df$id, subset(era.e.list.info, id %in% era.e.bam_names$V2)$id)]
era.e.77.pca.df$alt.type <- subset(era.e.list.info, id %in% era.e.bam_names$V2)$alt.type[match(era.e.77.pca.df$id, subset(era.e.list.info, id %in% era.e.bam_names$V2)$id)]
era.e.77.pca.df$country<- subset(era.e.list.info, id %in% era.e.bam_names$V2)$country[match(era.e.77.pca.df$id, subset(era.e.list.info, id %in% era.e.bam_names$V2)$id)]
era.e.77.pca.df$alt.type.country <- paste(era.e.77.pca.df$alt.type, era.e.77.pca.df$country, sep=".")

era.e.77.pca.df$region <- region

# plot pcas
# loadings
summary(era.e.77.pca)$importance[2,][1]
summary(era.e.77.pca)$importance[2,][2]

era.e.77.pca.plot <- ggplot(data=era.e.77.pca.df, aes(x=PC2, y=-PC1, color=altitude))+
  scale_color_gradient( low = "#00FF00B3", high = "#0000FFF2")+
  xlab("PC2 (20%)")+ylab("PC1 (66%)")+
  geom_point(aes(shape=alt.type.country), size=4)+
  scale_shape_manual(values=c(17,24,15,0,19,1))+
  #annotate(geom = 'text', label = era.e.regions.list.sig[[6]] , x =  0.2*(min(era.e.pca.list.sig.df$altitude) + max(era.e.pca.list.sig.df$altitude)), y = -Inf, hjust = 0, vjust = -1, size=4)+
  #annotate(geom = 'text', label = paste("PC1 ", round(summary(era.e.pca.list.sig[[i]])$importance[2,][1], digits = 2)*100, "%", sep="") ,  x =  0*(min(era.e.pca.list.sig.df$altitude) + max(era.e.pca.list.sig.df$altitude)), y = Inf, hjust = 0, vjust = 3, size=4)+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1),#axis.ticks.x=element_blank(),
        axis.ticks.length =unit(-0.15, "cm"), #axis.ticks.margin=unit(0.5, "cm"),
        axis.title = element_text(size=16),axis.text.y = element_text(size=10),axis.text.x = element_text(size=10),
        #axis.text.y = element_blank(),
        #axis.text.x = element_blank(), axis.title = element_blank(), axis.ticks=element_blank(),
        #axis.text.y.right = element_text(size = 10, margin = unit(c(t = 0, r = 0, b = 0, l = -8), "mm"), colour="springgreen", face="bold"),
        panel.background = element_blank(), #legend.text = element_text(size=8), legend.title =  element_text(size=12),
        panel.grid.major = element_blank(),# axis.text.x=element_blank(), #axis.text.y=element_blank(), #
        plot.margin = unit(c(0.1, 0.5, 0.1,0.5), "cm"),
        panel.grid.minor = element_blank(),# axis.title.x=element_blank(),#axis.title.y=element_blank(),
        legend.text.align = 0, legend.position = "none",strip.text = element_blank(),rect = element_rect(fill = "transparent"),plot.background = element_blank(),
        strip.background =element_rect(fill="transparent", color="transparent")); era.e.77.pca.plot


####### PC1 vs altitude #######
era.e.77.pca.alt.plot <- ggplot(data=era.e.77.pca.df, aes(y=-PC1, x=altitude, color=altitude))+
  scale_color_gradient( low = "#00FF00B3", high = "#0000FFF2")+
  ylab("PC1 (66%)")+xlab("Altitude")+
  geom_point(aes(shape=alt.type.country), size=4)+
  geom_smooth(method="lm", color="black")+
  scale_shape_manual(values=c(17,24,15,0,19,1))+
  #annotate(geom = 'text', label = era.e.regions.list.sig[[6]] , x =  0.2*(min(era.e.pca.list.sig.df$altitude) + max(era.e.pca.list.sig.df$altitude)), y = -Inf, hjust = 0, vjust = -1, size=4)+
  #annotate(geom = 'text', label = paste("PC1 ", round(summary(era.e.pca.list.sig[[i]])$importance[2,][1], digits = 2)*100, "%", sep="") ,  x =  0*(min(era.e.pca.list.sig.df$altitude) + max(era.e.pca.list.sig.df$altitude)), y = Inf, hjust = 0, vjust = 3, size=4)+
  stat_cor( hjust = -0.02, vjust = 0.01, size=4, aes(label =paste(..r.label.., cut(..p..,  breaks = c(-Inf, 0.0001, 0.001, 0.01, 0.05, Inf), labels = c("'****'", "'***'", "'**'", "'*'", "'ns'")), sep = "~"))) +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1),#axis.ticks.x=element_blank(),
        axis.ticks.length =unit(-0.15, "cm"), #axis.ticks.margin=unit(0.5, "cm"),
        axis.title = element_text(size=16),axis.text.y = element_text(size=10),axis.text.x = element_text(size=10),
        #axis.text.y = element_blank(),
        #axis.text.x = element_blank(), axis.title = element_blank(), axis.ticks=element_blank(),
        #axis.text.y.right = element_text(size = 10, margin = unit(c(t = 0, r = 0, b = 0, l = -8), "mm"), colour="springgreen", face="bold"),
        panel.background = element_blank(), #legend.text = element_text(size=8), legend.title =  element_text(size=12),
        panel.grid.major = element_blank(),# axis.text.x=element_blank(), #axis.text.y=element_blank(), #
        plot.margin = unit(c(0.1, 0.5, 0.1,0.5), "cm"),
        panel.grid.minor = element_blank(),# axis.title.x=element_blank(),#axis.title.y=element_blank(),
        legend.text.align = 0, legend.position = "none",strip.text = element_blank(),rect = element_rect(fill = "transparent"),plot.background = element_blank(),
        strip.background =element_rect(fill="transparent", color="transparent")); era.e.77.pca.alt.plot



####### concatenate 3, panel A #######

zpbs.shdr77 | era.e.77.pca.plot | era.e.77.pca.alt.plot
cowplot::plot_grid(zpbs.shdr77 , era.e.77.pca.plot, era.e.77.pca.alt.plot, nrow = 1)
ggsave("local/data/local/data/local/data/local/data/22_pop.gen.paper/figures/fig4.inv.local.pca/v2/local.pca.example.png", width = 27, height = 8, units = 'cm', bg="transparent")


##################### PANEL B- BARPLOT NUMBER OF SIGNIFICANT CORR WITH ALTITUDE #####################
head(era.e.shdr.outlier)
era.e.shdr.outlier$pc1.alt.altitude.pval.is.signif <- if_else(era.e.shdr.outlier$pc1.alt.altitude.pval <=0.05, "yes", "no"); era.e.shdr.outlier$pc1.alt.altitude.pval.is.signif
era.e.shdr.outlier<-subset(era.e.shdr.outlier, !(is.na(era.e.shdr.outlier$pc1.alt.altitude.pval)))
era.w.shdr.outlier$pc1.alt.altitude.pval.is.signif <- if_else(era.w.shdr.outlier$pc1.alt.altitude.pval <=0.05, "yes", "no"); era.w.shdr.outlier$pc1.alt.altitude.pval.is.signif
era.w.shdr.outlier<-subset(era.w.shdr.outlier, !(is.na(era.w.shdr.outlier$pc1.alt.altitude.pval)))

era.e.shdr.outlier.summ <- summarise(group_by(era.e.shdr.outlier, is.allo.hdr,pc1.alt.altitude.pval.is.signif),
          n=n()); era.e.shdr.outlier.summ 
era.w.shdr.outlier.summ <- summarise(group_by(era.w.shdr.outlier, is.allo.hdr,pc1.alt.altitude.pval.is.signif),
                                     n=n()); era.w.shdr.outlier.summ 
era.e.shdr.outlier.summ$side <- c("East"); era.w.shdr.outlier.summ$side <- c("West")
era.e.shdr.outlier.summ$is.allo.hdr.side <- paste( era.e.shdr.outlier.summ$is.allo.hdr, era.e.shdr.outlier.summ$side)
era.w.shdr.outlier.summ$is.allo.hdr.side <- paste( era.w.shdr.outlier.summ$is.allo.hdr, era.w.shdr.outlier.summ$side)
no.pca.outlier.era.summ <- rbind(era.e.shdr.outlier.summ, era.w.shdr.outlier.summ)
no.pca.outlier.era.summ$side <- factor(no.pca.outlier.era.summ$side, levels=c("West", "East"))

era.stats.summ.p <- ggplot(data=no.pca.outlier.era.summ , aes(fill=pc1.alt.altitude.pval.is.signif , x=is.allo.hdr, y=n)) + 
  geom_col(aes(color=is.allo.hdr.side), size=2, width = 0.75)+ 
  ylab("PC1 correlates with altitude (No. SHDR)") + 
  scale_color_manual(values=c( "#049E73","#0372B2","#D65D00","#D65D00" )) +
  scale_y_continuous( expand = c(0,0), limits = c(0,100))+
  scale_fill_manual(values=c("white", "black")) +theme_bw()+facet_wrap(~side)+
  theme(#axis.title.y = element_text(size=14,vjust=-1),
    plot.margin = unit(c(1,0.2,1.5,0.2), "lines"),axis.text.x = element_blank(), #axis.text.x = element_blank(),
    axis.text.y = element_text(size=10),axis.title = element_blank(),  axis.ticks.x = element_blank(),
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_rect(fill = "transparent",  size = 1),
    panel.background = element_rect(fill = "transparent"), plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
    legend.background = element_rect(fill = "transparent"), legend.box.background = element_rect(fill = "transparent"), strip.background.y = element_blank(),
    strip.text = element_text(size=17, face = "bold"), strip.background = element_blank(), legend.position = "none"); era.stats.summ.p
ggsave("local/data/local/data/local/data/local/data/22_pop.gen.paper/figures/fig4.inv.local.pca/v3/era.no.sig.pca.alt.corr.png", width =8, height = 7, units = 'cm', bg="transparent")

mel.e.shdr.outlier$pc1.alt.altitude.pval.is.signif <- if_else(mel.e.shdr.outlier$pc1.alt.altitude.pval <=0.05, "yes", "no"); mel.e.shdr.outlier$pc1.alt.altitude.pval.is.signif
mel.e.shdr.outlier<-subset(mel.e.shdr.outlier, !(is.na(mel.e.shdr.outlier$pc1.alt.altitude.pval)))
mel.w.shdr.outlier$pc1.alt.altitude.pval.is.signif <- if_else(mel.w.shdr.outlier$pc1.alt.altitude.pval <=0.05, "yes", "no"); mel.w.shdr.outlier$pc1.alt.altitude.pval.is.signif
mel.w.shdr.outlier<-subset(mel.w.shdr.outlier, !(is.na(mel.w.shdr.outlier$pc1.alt.altitude.pval)))

mel.e.shdr.outlier.summ <- summarise(group_by(mel.e.shdr.outlier, is.allo.hdr,pc1.alt.altitude.pval.is.signif),
                                     n=n()); mel.e.shdr.outlier.summ 
mel.w.shdr.outlier.summ <- summarise(group_by(mel.w.shdr.outlier, is.allo.hdr,pc1.alt.altitude.pval.is.signif),
                                     n=n()); mel.w.shdr.outlier.summ 
mel.e.shdr.outlier.summ$side <- c("East"); mel.w.shdr.outlier.summ$side <- c("West")
mel.e.shdr.outlier.summ$is.allo.hdr.side <- paste( mel.e.shdr.outlier.summ$is.allo.hdr, mel.e.shdr.outlier.summ$side)
mel.w.shdr.outlier.summ$is.allo.hdr.side <- paste( mel.w.shdr.outlier.summ$is.allo.hdr, mel.w.shdr.outlier.summ$side)
no.pca.outlier.mel.summ <- rbind(mel.e.shdr.outlier.summ, mel.w.shdr.outlier.summ)
no.pca.outlier.mel.summ$side <- factor(no.pca.outlier.mel.summ$side, levels=c("West", "East"))

mel.stats.summ.p <- ggplot(data=no.pca.outlier.mel.summ , aes(fill=pc1.alt.altitude.pval.is.signif , x=is.allo.hdr, y=n)) + 
  geom_col(aes(color=is.allo.hdr.side), size=2, width = 0.75)+ 
  ylab("PC1 correlates with altitude (No. SHDR)") + 
  scale_color_manual(values=c( "#049E73","#0372B2","#D65D00","#D65D00" )) +
  scale_y_continuous( expand = c(0,0), limits = c(0,100))+
  #scale_x_discrete(limits=c("Parapatric", "Allopatric"))+
  scale_fill_manual(values=c("white", "black")) +theme_bw()+facet_wrap(~side)+
  theme(#axis.title.y = element_text(size=14,vjust=-1),
        plot.margin = unit(c(1,0.2,1.5,0.2), "lines"),axis.text.x = element_blank(), #axis.text.x = element_blank(),
        axis.text.y = element_text(size=10),axis.title = element_blank(),  axis.ticks.x = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_rect(fill = "transparent",  size = 1),
        panel.background = element_rect(fill = "transparent"), plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
        legend.background = element_rect(fill = "transparent"), legend.box.background = element_rect(fill = "transparent"), strip.background.y = element_blank(),
        strip.text = element_text(size=17, face = "bold"), strip.background = element_blank(), legend.position = "none"); mel.stats.summ.p
ggsave("local/data/local/data/local/data/local/data/22_pop.gen.paper/figures/fig4.inv.local.pca/v3/mel.no.sig.pca.alt.corr.png", width =8, height = 7, units = 'cm', bg="transparent")

plot_grid(era.stats.summ.p, mel.stats.summ.p, ncol = 1)
ggsave("local/data/local/data/local/data/local/data/22_pop.gen.paper/figures/fig4.inv.local.pca/v3/no.sig.pca.alt.corr.png", width =7, height = 16, units = 'cm', bg="transparent")



##################### PANEL C 3. inv chr 2 ##############

min(subset(era.chr2.pbs, scaff!="Herato0201" &scaff!="Herato0202"&scaff!="Herato0203"&BP.wg>24220000 &BP.wg<34300000)$scaff)
min(subset(era.chr2.pbs, scaff!="Herato0201" &scaff!="Herato0202"&scaff!="Herato0203"& scaff=="Herato0204"&BP.wg>24220000 &BP.wg<34300000)$midPos)

max(subset(era.chr2.pbs, scaff!="Herato0201" &scaff!="Herato0202"&scaff!="Herato0203"&BP.wg>24220000 &BP.wg<34300000)$scaff)
max(subset(era.chr2.pbs, scaff!="Herato0201" &scaff!="Herato0202"&scaff!="Herato0203"& scaff=="Herato0215"&BP.wg>24220000 &BP.wg<34300000)$midPos)


zpbs.co.ec.e <- ggplot(subset(era.chr2.pbs, scaff!="Herato0201" &scaff!="Herato0202"&scaff!="Herato0203"&scaff!="Herato0204"&scaff!="Herato0203"), aes(x=BP.wg, y=zPBS0.east.col)) +
  annotate("rect", xmin = 25750000, xmax = 32400000, ymin = 0, ymax = Inf, colour = "transparent", fill="#EEEEEE", size=1) +
  
  geom_rect(inherit.aes = F, data=subset(era.e.shdr, chr==2 & scaff!="Herato0201" &scaff!="Herato0202"&scaff!="Herato0203"&scaff!="Herato0204"& end<33800000),
           aes(xmin=start, xmax=end,ymin=-Inf,ymax=Inf, fill=is.allo.hdr), colour="transparent", alpha=.4 )+
  scale_fill_manual(values=c("#049E73", "#D45D00"))+
  
  # # inverted region
  # annotate("rect", xmin = 25750000, xmax = 32400000, ymin = -1, ymax = -0.5, colour = "transparent", fill="#EEEEEE", size=1) +
  # annotate("segment", x = 25750000, xend = 25750000, y = -1, yend = -0.5, colour = "black", lty=1,size=0.7) +
  # annotate("segment", x = 32400000, xend = 32400000, y = -1, yend = -0.5, colour = "black", lty=1,size=0.7) +
  # # non inverted region
  # annotate("rect", xmin = 33100000, xmax = 34100000, ymin = -1, ymax = -0.5, colour = "transparent", fill="#EEEEEE", size=1) +
  # annotate("segment", x = 33100000, xend = 33100000, y = -1, yend = -0.5, colour = "black", lty=1,size=0.7) +
  # annotate("segment", x = 34100000, xend = 34100000, y = -1, yend = -0.5, colour = "black", lty=1,size=0.7) +
  #geom_rect(inherit.aes = F,data= subset(era.gff_gene, chr==2 &scaffold!="Herato0201" &scaffold!="Herato0202"&scaffold!="Herato0203"&scaffold!="Herato0204"&scaffold!="Herato0203"&scaffold!="Herato0204" & end.offset<33800000), aes(xmin=start.offset, xmax=end.offset,ymin=7.5,ymax=8) )+
  geom_line(inherit.aes = F, color="black", alpha=0.9, size=1, data=subset(era.chr2.pbs, scaff!="Herato0201" &scaff!="Herato0202"&scaff!="Herato0203"&BP.wg>24220000 &BP.wg<34300000), aes(x=BP.wg, y=rollmean(zPBS0.east.col,100, na.pad=TRUE)), lty="solid") +
  geom_line(inherit.aes = F, color="grey40", alpha=0.9, size=1, data=subset(era.chr2.pbs, scaff!="Herato0201" &scaff!="Herato0202"&scaff!="Herato0203"&BP.wg>24220000  & BP.wg<34300000), aes(x=BP.wg, y=rollmean(zPBS0.east.ecu, 100, na.pad=TRUE)), lty=c("11")) +
  #geom_label_repel(inherit.aes = F, data=subset(era.e.shdr, chr==2 & scaff!="Herato0201" &scaff!="Herato0202"&scaff!="Herato0203"&scaff!="Herato0204"& end<33800000 ), aes(x=(start+end)/2, y=4, label=substr(shdr.para.east.id,11,13 ), color=is.allo.hdr), min.segment.length = 0.05) +
  #scale_colour_manual(values=c("#049E73", "#D45D00"))+
  # annotate("rect", xmin = -Inf, xmax = 34300000, ymin = 3.8, ymax = 4.1, colour = "black", fill="white", size=0.5) +
  # geom_rect(inherit.aes = F, data=subset(era.ref.chr.pos, chr==2 & scaff!="Herato0201" & scaff!="Herato0202"&scaff!="Herato0203"&scaff!="Herato0204"& end<34300000), 
  #           aes(xmin=start, xmax=end, ymin=3.8,ymax=4.1),  fill=c(rep(c("white","grey70"), 5 )),color="black", alpha=1, size=0.5 )+
  # coord_cartesian(ylim = c(-1,4.1), clip = 'off',expand = c(0, 0)) + 
  scale_y_continuous(expand = c(0, 0))+
  scale_x_continuous(expand = c(0, 0), labels=function(x)x/10^6)+
  ylab("East zPBS high")+  xlab("Position (MB)")+ 
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1),#axis.ticks.x=element_blank(),
        axis.ticks.length =unit(-0.15, "cm"), #axis.ticks.margin=unit(0.5, "cm"),
        axis.title.y = element_text(size=16),axis.text.y = element_text(size=10),#axis.text.x = element_text(size=10),
        #axis.text.y = element_blank(),
        axis.text.x = element_blank(), axis.title = element_blank(), axis.ticks.x=element_blank(),
        #axis.text.y.right = element_text(size = 10, margin = unit(c(t = 0, r = 0, b = 0, l = -8), "mm"), colour="springgreen", face="bold"),
        panel.background = element_blank(), #legend.text = element_text(size=8), legend.title =  element_text(size=12),
        panel.grid.major = element_blank(),# axis.text.x=element_blank(), #axis.text.y=element_blank(), #
        plot.margin = unit(c(0.1, 0.5, 0.1,0.5), "cm"),
        panel.grid.minor = element_blank(),# axis.title.x=element_blank(),#axis.title.y=element_blank(),
        legend.text.align = 0, legend.position = "none",strip.text = element_blank(),rect = element_rect(fill = "transparent"),
        strip.background =element_rect(fill="transparent", color="transparent")); zpbs.co.ec.e 

ggsave("local/data/local/data/local/data/local/data/22_pop.gen.paper/figures/fig4.inv.local.pca/v2/chr.2.simple1.png", width = 18, height = 7, units = 'cm')




############## 4. PC1 vs alt inv/non-inv chr 2 ##############
era.e.inv.npy <-npyLoad(file="output/mds.shared/era.e.inv/era.e.region.inv2.full.cov.npy")
names(era.e.inv.npy) <- "shdr.east.inv"
region <-"shdr.east.inv"
era.e.inv.pca <- region
era.e.inv.pca<-prcomp(era.e.inv.npy)
rownames(era.e.inv.pca$x)<-era.e.bam_names$V2

# store pca info for plotting/analyses
era.e.inv.pca.df <- as.data.frame(era.e.inv.pca$x[,c(1:2)]); era.e.inv.pca.df$id <- rownames(era.e.inv.pca.df)
era.e.inv.pca.df$altitude <- subset(era.e.list.info, id %in% era.e.bam_names$V2)$altitude[match(era.e.inv.pca.df$id, subset(era.e.list.info, id %in% era.e.bam_names$V2)$id)]
era.e.inv.pca.df$alt.type <- subset(era.e.list.info, id %in% era.e.bam_names$V2)$alt.type[match(era.e.inv.pca.df$id, subset(era.e.list.info, id %in% era.e.bam_names$V2)$id)]
era.e.inv.pca.df$country<- subset(era.e.list.info, id %in% era.e.bam_names$V2)$country[match(era.e.inv.pca.df$id, subset(era.e.list.info, id %in% era.e.bam_names$V2)$id)]
era.e.inv.pca.df$subsp<- subset(era.e.list.info, id %in% era.e.bam_names$V2)$subsp[match(era.e.inv.pca.df$id, subset(era.e.list.info, id %in% era.e.bam_names$V2)$id)]
era.e.inv.pca.df$col.pattern<- if_else(era.e.inv.pca.df$subsp=="dignus", "black", "rayed")
era.e.inv.pca.df$alt.type.country <- paste(era.e.inv.pca.df$alt.type, era.e.inv.pca.df$country, sep=".")
era.e.inv.pca.df$region <- region

summary(era.e.inv.pca)$importance[2,][1]
era.e.inv.pca.alt.plot <- ggplot(data=era.e.inv.pca.df, aes(y=-PC1, x=altitude, color=altitude))+
  scale_color_gradient( low = "#00FF00B3", high = "#0000FFF2")+
  ylab("Local PCA PC1 (87%)")+xlab("Altitude")+
  geom_smooth(method="lm", color="black")+
  geom_point(aes(shape=alt.type.country), size=4)+
  scale_shape_manual(values=c(17,24,15,0,19,1))+
  #annotate(geom = 'text', label = era.e.regions.list.sig[[6]] , x =  0.2*(min(era.e.pca.list.sig.df$altitude) + max(era.e.pca.list.sig.df$altitude)), y = -Inf, hjust = 0, vjust = -1, size=4)+
  #annotate(geom = 'text', label = paste("PC1 ", round(summary(era.e.pca.list.sig[[i]])$importance[2,][1], digits = 2)*100, "%", sep="") ,  x =  0*(min(era.e.pca.list.sig.df$altitude) + max(era.e.pca.list.sig.df$altitude)), y = Inf, hjust = 0, vjust = 3, size=4)+
  stat_cor( hjust = -0.6, vjust = 0.8, size=6, aes(label =paste(..r.label.., cut(..p..,  breaks = c(-Inf, 0.0001, 0.001, 0.01, 0.05, Inf), labels = c("'****'", "'***'", "'**'", "'*'", "'ns'")), sep = "~"))) +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1),#axis.ticks.x=element_blank(),
        axis.ticks.length =unit(-0.15, "cm"), #axis.ticks.margin=unit(0.5, "cm"),
        axis.title = element_text(size=16),axis.text.y = element_text(size=10),axis.text.x = element_text(size=10),
        #axis.text.y = element_blank(),
        #axis.text.x = element_blank(), axis.title = element_blank(), axis.ticks=element_blank(),
        #axis.text.y.right = element_text(size = 10, margin = unit(c(t = 0, r = 0, b = 0, l = -8), "mm"), colour="springgreen", face="bold"),
        panel.background = element_blank(), #legend.text = element_text(size=8), legend.title =  element_text(size=12),
        panel.grid.major = element_blank(),# axis.text.x=element_blank(), #axis.text.y=element_blank(), #
        plot.margin = unit(c(0.1, 0.5, 0.1,0.5), "cm"),
        panel.grid.minor = element_blank(),# axis.title.x=element_blank(),#axis.title.y=element_blank(),
        legend.text.align = 0, legend.position = "none",strip.text = element_blank(),rect = element_rect(fill = "transparent"),
        strip.background =element_rect(fill="transparent", color="transparent")); era.e.inv.pca.alt.plot

ggsave("local/data/local/data/local/data/local/data/22_pop.gen.paper/figures/fig4.inv.local.pca/v2/era.e.chr2.inv.pca.alt.plot.png", width = 9.5, height = 8.5, units = 'cm', bg="transparent")

##### non inv ###
era.e.non.inv.npy <-npyLoad(file="output/mds.shared/era.e.inv/era.e.region.no.inv2.Herato0215.cov.npy")
era.e.non.inv.npy <-npyLoad(file="output/mds.shared/era.e.inv/era.e.region.non.inv2.full.cov.npy")
era.e.non.inv.npy <-npyLoad(file="output/mds.shared/era.mel.chr1.neutral/era.e.chr1.cov.npy")
era.e.non.inv.npy <-npyLoad(file="output/mds.shared/era.mel.wg.neutral/era.e.wg.cov.npy")

names(era.e.non.inv.npy) <- "shdr.east.non.inv"
region <-"shdr.east.non.inv"
era.e.non.inv.pca <- region
era.e.non.inv.pca<-prcomp(era.e.non.inv.npy)
rownames(era.e.non.inv.pca$x)<-era.e.bam_names$V2

# store pca info for plotting/analyses
era.e.non.inv.pca.df <- as.data.frame(era.e.non.inv.pca$x[,c(1:2)]); era.e.non.inv.pca.df$id <- rownames(era.e.non.inv.pca.df)
era.e.non.inv.pca.df$altitude <- subset(era.e.list.info, id %in% era.e.bam_names$V2)$altitude[match(era.e.non.inv.pca.df$id, subset(era.e.list.info, id %in% era.e.bam_names$V2)$id)]
era.e.non.inv.pca.df$alt.type <- subset(era.e.list.info, id %in% era.e.bam_names$V2)$alt.type[match(era.e.non.inv.pca.df$id, subset(era.e.list.info, id %in% era.e.bam_names$V2)$id)]
era.e.non.inv.pca.df$country<- subset(era.e.list.info, id %in% era.e.bam_names$V2)$country[match(era.e.non.inv.pca.df$id, subset(era.e.list.info, id %in% era.e.bam_names$V2)$id)]
era.e.non.inv.pca.df$subsp<- subset(era.e.list.info, id %in% era.e.bam_names$V2)$subsp[match(era.e.non.inv.pca.df$id, subset(era.e.list.info, id %in% era.e.bam_names$V2)$id)]
era.e.non.inv.pca.df$col.pattern<- if_else(era.e.non.inv.pca.df$subsp=="dignus", "black", "rayed")
era.e.non.inv.pca.df$alt.type.country <- paste(era.e.non.inv.pca.df$alt.type, era.e.non.inv.pca.df$country, sep=".")
era.e.non.inv.pca.df$region <- region

summary(lm(PC1 ~ altitude + col.pattern, data=era.e.inv.pca.df ))
summary(lm(PC1 ~ altitude + col.pattern, data=era.e.non.inv.pca.df ))

summary(era.e.non.inv.pca)$importance[2,][1]
era.e.non.inv.pca.alt.plot <- ggplot(data=era.e.non.inv.pca.df, aes(y=PC1, x=altitude, color=altitude))+
  scale_color_gradient( low = "#00FF00B3", high = "#0000FFF2")+
  ylab("Local PCA PC1 (9%)")+xlab("Altitude")+
  geom_point(aes(shape=alt.type.country), size=4)+
  geom_smooth(method="lm", color="black")+
  scale_shape_manual(values=c(17,24,15,0,19,1))+
  #annotate(geom = 'text', label = era.e.regions.list.sig[[6]] , x =  0.2*(min(era.e.pca.list.sig.df$altitude) + max(era.e.pca.list.sig.df$altitude)), y = -Inf, hjust = 0, vjust = -1, size=4)+
  #annotate(geom = 'text', label = paste("PC1 ", round(summary(era.e.pca.list.sig[[i]])$importance[2,][1], digits = 2)*100, "%", sep="") ,  x =  0*(min(era.e.pca.list.sig.df$altitude) + max(era.e.pca.list.sig.df$altitude)), y = Inf, hjust = 0, vjust = 3, size=4)+
  stat_cor( hjust = -0.6, vjust = 0.8, size=6, aes(label =paste(..r.label.., cut(..p..,  breaks = c(-Inf, 0.0001, 0.001, 0.01, 0.05, Inf), labels = c("'****'", "'***'", "'**'", "'*'", "'ns'")), sep = "~"))) +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1),#axis.ticks.x=element_blank(),
        axis.ticks.length =unit(-0.15, "cm"), #axis.ticks.margin=unit(0.5, "cm"),
        axis.title = element_text(size=16),axis.text.y = element_text(size=10),axis.text.x = element_text(size=10),
        #axis.text.y = element_blank(),
        #axis.text.x = element_blank(), axis.title = element_blank(), axis.ticks=element_blank(),
        #axis.text.y.right = element_text(size = 10, margin = unit(c(t = 0, r = 0, b = 0, l = -8), "mm"), colour="springgreen", face="bold"),
        panel.background = element_blank(), #legend.text = element_text(size=8), legend.title =  element_text(size=12),
        panel.grid.major = element_blank(),# axis.text.x=element_blank(), #axis.text.y=element_blank(), #
        plot.margin = unit(c(0.1, 0.5, 0.1,0.5), "cm"),
        panel.grid.minor = element_blank(),# axis.title.x=element_blank(),#axis.title.y=element_blank(),
        legend.text.align = 0, legend.position = "none",strip.text = element_blank(),rect = element_rect(fill = "transparent"),
        strip.background =element_rect(fill="transparent", color="transparent")); era.e.non.inv.pca.alt.plot

ggsave("local/data/local/data/local/data/local/data/22_pop.gen.paper/figures/fig4.inv.local.pca/v2/era.e.chr2.non.inv.pca.alt.plot.png", width = 10, height = 8.5, units = 'cm')
ggsave("local/data/local/data/local/data/local/data/22_pop.gen.paper/figures/fig4.inv.local.pca/v2/era.e.chr2.no.inv2.Herato0215.pca.alt.plot.png", width = 10, height = 8.5, units = 'cm')
ggsave("local/data/local/data/local/data/local/data/22_pop.gen.paper/figures/fig4.inv.local.pca/v2/era.e.wg.neutral.pca.alt.plot.png", width = 10, height = 8.5, units = 'cm')


############## 5. PCA inv/non-inv chr 2 ##############
era.e.inv.npy <-npyLoad(file="output/mds.shared/era.e.inv/era.e.region.inv2.full.cov.npy")
names(era.e.inv.npy) <- "shdr.east.inv"
region <-"shdr.east.inv"
era.e.inv.pca <- region
era.e.inv.pca<-prcomp(era.e.inv.npy)
rownames(era.e.inv.pca$x)<-era.e.bam_names$V2

# store pca info for plotting/analyses
era.e.inv.pca.df <- as.data.frame(era.e.inv.pca$x[,c(1:2)]); era.e.inv.pca.df$id <- rownames(era.e.inv.pca.df)
era.e.inv.pca.df$altitude <- subset(era.e.list.info, id %in% era.e.bam_names$V2)$altitude[match(era.e.inv.pca.df$id, subset(era.e.list.info, id %in% era.e.bam_names$V2)$id)]
era.e.inv.pca.df$alt.type <- subset(era.e.list.info, id %in% era.e.bam_names$V2)$alt.type[match(era.e.inv.pca.df$id, subset(era.e.list.info, id %in% era.e.bam_names$V2)$id)]
era.e.inv.pca.df$country<- subset(era.e.list.info, id %in% era.e.bam_names$V2)$country[match(era.e.inv.pca.df$id, subset(era.e.list.info, id %in% era.e.bam_names$V2)$id)]
era.e.inv.pca.df$alt.type.country <- paste(era.e.inv.pca.df$alt.type, era.e.inv.pca.df$country, sep=".")
era.e.inv.pca.df$region <- region

summary(era.e.inv.pca)$importance[2,][1]
summary(era.e.inv.pca)$importance[2,][2]
era.e.inv.pca.plot <- ggplot(data=era.e.inv.pca.df, aes(x=PC1, y=PC2, color=altitude))+
  scale_color_gradient( low = "#00FF00B3", high = "#0000FFF2")+
  xlab("PC1 (87%)")+ylab("PC2 (20%)")+
  geom_point(aes(shape=alt.type.country), size=4)+
  scale_shape_manual(values=c(17,24,15,0,19,1))+
  #annotate(geom = 'text', label = era.e.regions.list.sig[[6]] , x =  0.2*(min(era.e.pca.list.sig.df$altitude) + max(era.e.pca.list.sig.df$altitude)), y = -Inf, hjust = 0, vjust = -1, size=4)+
  #annotate(geom = 'text', label = paste("PC1 ", round(summary(era.e.pca.list.sig[[i]])$importance[2,][1], digits = 2)*100, "%", sep="") ,  x =  0*(min(era.e.pca.list.sig.df$altitude) + max(era.e.pca.list.sig.df$altitude)), y = Inf, hjust = 0, vjust = 3, size=4)+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1),#axis.ticks.x=element_blank(),
        axis.ticks.length =unit(-0.15, "cm"), #axis.ticks.margin=unit(0.5, "cm"),
        axis.title = element_text(size=16),axis.text.y = element_text(size=10),axis.text.x = element_text(size=10),
        #axis.text.y = element_blank(),
        #axis.text.x = element_blank(), axis.title = element_blank(), axis.ticks=element_blank(),
        #axis.text.y.right = element_text(size = 10, margin = unit(c(t = 0, r = 0, b = 0, l = -8), "mm"), colour="springgreen", face="bold"),
        panel.background = element_blank(), #legend.text = element_text(size=8), legend.title =  element_text(size=12),
        panel.grid.major = element_blank(),# axis.text.x=element_blank(), #axis.text.y=element_blank(), #
        plot.margin = unit(c(0.1, 0.5, 0.1,0.5), "cm"),
        panel.grid.minor = element_blank(),# axis.title.x=element_blank(),#axis.title.y=element_blank(),
        legend.text.align = 0, legend.position = "none",strip.text = element_blank(),rect = element_rect(fill = "transparent"),
        strip.background =element_rect(fill="transparent", color="transparent")); era.e.inv.pca.plot

ggsave("local/data/local/data/local/data/local/data/22_pop.gen.paper/figures/fig4.inv.local.pca/v2/era.e.chr2.inv.pca.plot.png", width = 10, height = 8.5, units = 'cm')

##### non .inv ###
era.e.non.inv.npy <-npyLoad(file="output/mds.shared/era.e.inv/era.e.region.non.inv2.full.cov.npy")
era.e.non.inv.npy <-npyLoad(file="output/mds.shared/era.mel.chr1.neutral/era.e.chr1.cov.npy")
#era.e.non.inv.npy <-npyLoad(file="output/mds.shared/era.mel.wg.neutral/era.e.wg.cov.npy")

names(era.e.non.inv.npy) <- "shdr.east.non.inv"
region <-"shdr.east.non.inv"
era.e.non.inv.pca <- region
era.e.non.inv.pca<-prcomp(era.e.non.inv.npy)
rownames(era.e.non.inv.pca$x)<-era.e.bam_names$V2

# store pca info for plotting/analyses
era.e.non.inv.pca.df <- as.data.frame(era.e.non.inv.pca$x[,c(1:2)]); era.e.non.inv.pca.df$id <- rownames(era.e.non.inv.pca.df)
era.e.non.inv.pca.df$altitude <- subset(era.e.list.info, id %in% era.e.bam_names$V2)$altitude[match(era.e.non.inv.pca.df$id, subset(era.e.list.info, id %in% era.e.bam_names$V2)$id)]
era.e.non.inv.pca.df$alt.type <- subset(era.e.list.info, id %in% era.e.bam_names$V2)$alt.type[match(era.e.non.inv.pca.df$id, subset(era.e.list.info, id %in% era.e.bam_names$V2)$id)]
era.e.non.inv.pca.df$country<- subset(era.e.list.info, id %in% era.e.bam_names$V2)$country[match(era.e.non.inv.pca.df$id, subset(era.e.list.info, id %in% era.e.bam_names$V2)$id)]
era.e.non.inv.pca.df$alt.type.country <- paste(era.e.non.inv.pca.df$alt.type, era.e.non.inv.pca.df$country, sep=".")
era.e.non.inv.pca.df$region <- region

subset(era.e.non.inv.pca.df, PC2>0.5)

summary(era.e.non.inv.pca)$importance[2,][1]
summary(era.e.non.inv.pca)$importance[2,][2]
era.e.non.inv.pca.plot <- ggplot(data=era.e.non.inv.pca.df, aes(x=PC1, y=PC2, color=altitude))+
  scale_color_gradient( low = "#00FF00B3", high = "#0000FFF2")+
  xlab("PC1 (11%)")+ylab("PC2 (2%)")+
  geom_point(aes(shape=alt.type.country), size=4)+
  scale_shape_manual(values=c(17,24,15,0,19,1))+
  #annotate(geom = 'text', label = era.e.regions.list.sig[[6]] , x =  0.2*(min(era.e.pca.list.sig.df$altitude) + max(era.e.pca.list.sig.df$altitude)), y = -Inf, hjust = 0, vjust = -1, size=4)+
  #annotate(geom = 'text', label = paste("PC1 ", round(summary(era.e.pca.list.sig[[i]])$importance[2,][1], digits = 2)*100, "%", sep="") ,  x =  0*(min(era.e.pca.list.sig.df$altitude) + max(era.e.pca.list.sig.df$altitude)), y = Inf, hjust = 0, vjust = 3, size=4)+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1),#axis.ticks.x=element_blank(),
        axis.ticks.length =unit(-0.15, "cm"), #axis.ticks.margin=unit(0.5, "cm"),
        axis.title = element_text(size=16),axis.text.y = element_text(size=10),axis.text.x = element_text(size=10),
        #axis.text.y = element_blank(),
        #axis.text.x = element_blank(), axis.title = element_blank(), axis.ticks=element_blank(),
        #axis.text.y.right = element_text(size = 10, margin = unit(c(t = 0, r = 0, b = 0, l = -8), "mm"), colour="springgreen", face="bold"),
        panel.background = element_blank(), #legend.text = element_text(size=8), legend.title =  element_text(size=12),
        panel.grid.major = element_blank(),# axis.text.x=element_blank(), #axis.text.y=element_blank(), #
        plot.margin = unit(c(0.1, 0.5, 0.1,0.5), "cm"),
        panel.grid.minor = element_blank(),# axis.title.x=element_blank(),#axis.title.y=element_blank(),
        legend.text.align = 0, legend.position = "none",strip.text = element_blank(),rect = element_rect(fill = "transparent"),
        strip.background =element_rect(fill="transparent", color="transparent")); era.e.non.inv.pca.plot

ggsave("local/data/local/data/local/data/local/data/22_pop.gen.paper/figures/fig4.inv.local.pca/v2/era.e.chr1.neutral.pca.plot.png", width = 9.5, height = 8.5, units = 'cm')

####### concatenate 2 PCAs, panel B ####### 

plot_grid( era.e.non.inv.pca.plot, era.e.inv.pca.plot )
ggsave("local/data/local/data/local/data/local/data/22_pop.gen.paper/figures/fig4.inv.local.pca/v2/era.e.chr1.neutral.and.inv.pca.plot.png", width = 19, height = 8.5, units = 'cm')




#################################### 2. SI CHR 2 INV FDM and TREE  ##############
############## fdm ##############
zpbs.co.ec.e.fullchr2 <- ggplot(subset(era.chr2.pbs), aes(x=BP.wg, y=zPBS0.east.col)) +
  geom_rect(inherit.aes = F, data=subset(era.e.shdr, chr==2),
           aes(xmin=start, xmax=end,ymin=-Inf,ymax=Inf, fill=is.allo.hdr), colour="grey50", alpha=.4 )+
  scale_fill_manual(values=c("#049E73", "#D45D00"))+
  #geom_rect(inherit.aes = F,data= subset(era.gff_gene, chr==2 &scaffold!="Herato0201" &scaffold!="Herato0202"&scaffold!="Herato0203"&scaffold!="Herato0204"&scaffold!="Herato0203"&scaffold!="Herato0204" & end.offset<33800000), aes(xmin=start.offset, xmax=end.offset,ymin=7.5,ymax=8) )+
  geom_line(inherit.aes = F, color="black", alpha=0.9, size=1, data=subset(era.chr2.pbs), aes(x=BP.wg, y=rollmean(zPBS0.east.col,50, na.pad=TRUE)), lty="solid") +
  geom_line(inherit.aes = F, color="grey40", alpha=0.9, size=1, data=subset(era.chr2.pbs), aes(x=BP.wg, y=rollmean(zPBS0.east.ecu, 50, na.pad=TRUE)), lty=c("11")) +
  geom_label_repel(inherit.aes = F, data=subset(era.e.shdr, chr==2  ), aes(x=(start+end)/2, y=12, label=substr(shdr.para.east.id,11,13 ), color=is.allo.hdr), min.segment.length = 0.05) +
  scale_colour_manual(values=c("#049E73", "#D45D00"))+
  annotate("rect", xmin = -Inf, xmax = 34300000, ymin = 14.2, ymax = 15.1, colour = "black", fill="white", size=0.5) +
  geom_rect(inherit.aes = F, data=subset(era.ref.chr.pos, chr==2 ),
            aes(xmin=start, xmax=end, ymin=14.2,ymax=15.1),  fill=c(rep(c("white","grey70"), 7 ), "white"),color="black", alpha=1, size=0.5 )+
  coord_cartesian(ylim = c(-1,15.1), clip = 'off',expand = c(0, 0)) + 
  scale_x_continuous(expand = c(0, 0))+
  ylab("East zPBS high")+ 
  theme(panel.border = element_rect(colour = "black", fill=NA, size=.5),legend.position="none",
        plot.margin = unit(c(0.75,0.1,0.1,0.1), "lines"),axis.text = element_text(size=10), axis.text.x = element_blank(),axis.ticks.x = element_blank(),
        axis.text.y = element_text(size=10),axis.title.y = element_text(size=12),axis.title.x = element_blank(), panel.background = element_rect(fill = NA), plot.background = element_rect(fill = NA, color = NA),
        axis.text.y.right = element_text(size = 10, margin = unit(c(t = 0, r = 0, b = 0, l = -4), "mm"), colour="#eb8055ff"),
        axis.title.y.right  = element_text(size = 12, margin = unit(c(t = 0, r = 0, b = 0, l = 2), "mm"),  colour="#eb8055ff", face="bold"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.background = element_rect(fill = "transparent"), legend.box.background = element_rect(fill = "transparent")); zpbs.co.ec.e.fullchr2


  
fdm.him.era <- ggplot(subset(fd_eraCvlowE_eraChighE_himHE_eleHW, chr=="02"), aes(x=mid.WG, y=fdM)) +
  geom_rect(inherit.aes = F, data=subset(era.e.shdr, chr==2 ), 
            aes(xmin=start, xmax=end,ymin=-Inf,ymax=Inf, fill=is.allo.hdr), colour="grey50", alpha=.4 )+
  scale_fill_manual(values=c("#049E73", "#D65D00"))+
  geom_hline(yintercept = 0, lty="dashed")+
    geom_line(color="black", alpha=0.9, size=1, lty="solid") +
  ylab("H. himera fdM")+
    geom_line(inherit.aes = F,data=subset(fd_eraEvlowE_eraEhighE_himHE_eleHW, chr=="02"), aes(x=mid.WG, y=fdM), color="grey30", alpha=0.9, size=1, lty=c("11")) +
  #geom_label_repel(inherit.aes = F, data=subset(era.e.shdr, chr==2  ), aes(x=(start+end)/2, y=0.28, label=substr(shdr.para.east.id,11,13 ), color=is.allo.hdr), min.segment.length = 0.05) +
  scale_colour_manual(values=c("#049E73", "#D65D00"))+
  coord_cartesian(ylim = c(-0.12,0.35), clip = 'off',expand = c(0, 0)) + 
    theme(panel.border = element_rect(colour = "#eb8055ff", fill=NA, size=1),legend.position="none",
          plot.margin = unit(c(0.1,0.1,0.1,0.1), "lines"),axis.text = element_text(size=10), axis.text.x = element_blank(),axis.ticks.x = element_blank(),
          axis.text.y = element_text(size=10),axis.title.y = element_text(size=12),axis.title.x = element_blank(), panel.background = element_rect(fill = NA), plot.background = element_rect(fill = NA, color = NA),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.background = element_rect(fill = "transparent"), legend.box.background = element_rect(fill = "transparent"));   fdm.him.era

fdm.clyEhighE.era <- ggplot(subset(fd_eraCvlowE_eraChighE_clyEhighE_eleHW, chr=="02"), aes(x=mid.WG, y=fdM)) +
  geom_rect(inherit.aes = F, data=subset(era.e.shdr, chr==2 ), 
            aes(xmin=start, xmax=end,ymin=-Inf,ymax=Inf, fill=is.allo.hdr), colour="grey50", alpha=.4 )+
  scale_fill_manual(values=c("#049E73", "#D65D00"))+
  geom_hline(yintercept = 0, lty="dashed")+
  
  geom_line(color="black", alpha=0.9, size=1, lty="solid") +
  ylab("H. clysonymus fdM")+
  geom_line(inherit.aes = F,data=subset(fd_eraEvlowE_eraEhighE_clyEhighE_eleHW, chr=="02"), aes(x=mid.WG, y=fdM), color="grey30", alpha=0.9, size=1, lty=c("11")) +
  coord_cartesian(ylim = c(-0.1,0.1), clip = 'off',expand = c(0, 0)) + 
  theme(panel.border = element_rect(colour = "#efe350ff", fill=NA, size=1),legend.position="none",
        plot.margin = unit(c(0.1,0.1,0.1,0.1), "lines"),axis.text = element_text(size=10), axis.text.x = element_blank(),axis.ticks.x = element_blank(),
        axis.text.y = element_text(size=10),axis.title.y = element_text(size=12),axis.title.x = element_blank(), panel.background = element_rect(fill = NA), plot.background = element_rect(fill = NA, color = NA),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.background = element_rect(fill = "transparent"), legend.box.background = element_rect(fill = "transparent"));   fdm.clyEhighE.era


fdm.telEhighE.era <- ggplot(subset(fd_eraCvlowE_eraChighE_telhighE_eleHW, chr=="02"), aes(x=mid.WG, y=fdM)) +
  geom_rect(inherit.aes = F, data=subset(era.e.shdr, chr==2 ), 
            aes(xmin=start, xmax=end,ymin=-Inf,ymax=Inf, fill=is.allo.hdr), colour="grey50", alpha=.4 )+
  scale_fill_manual(values=c("#049E73", "#D65D00"))+
  geom_hline(yintercept = 0, lty="dashed")+
  geom_line(color="black", alpha=0.9, size=1, lty="solid") +
  ylab("H. telesiphe fdM")+xlab("Genome-wide adjusted position (bp)")+
  geom_line(inherit.aes = F,data=subset(fd_eraEvlowE_eraEhighE_telhighE_eleHW, chr=="02"), aes(x=mid.WG, y=fdM), color="grey30", alpha=0.9, size=1, lty=c("11")) +
  coord_cartesian(ylim = c(-0.1,0.1), clip = 'off',expand = c(0, 0)) + 
  theme(panel.border = element_rect(colour = "#efe350ff", fill=NA, size=1),legend.position="none",
        plot.margin = unit(c(0.1,0.1,0.1,0.1), "lines"),axis.text = element_text(size=10), #axis.text.x = element_blank(),axis.ticks.x = element_blank(),
        axis.text.y = element_text(size=10),axis.title.y = element_text(size=12), panel.background = element_rect(fill = NA), plot.background = element_rect(fill = NA, color = NA),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.background = element_rect(fill = "transparent"), legend.box.background = element_rect(fill = "transparent"));   fdm.telEhighE.era

(zpbs.co.ec.e.fullchr2/ fdm.him.era / fdm.clyEhighE.era/ fdm.telEhighE.era)+ plot_annotation(tag_levels = 'A') 
ggsave("local/data/local/data/local/data/local/data/22_pop.gen.paper/figures/supp.fdm.chr2/fdm.chr2.png", width = 18, height = 23, units = 'cm', bg="transparent")


############## chr 2 inv tree ##############
library(ggtree)
tree.era <- read.tree(file = "data/local.pca.tree/RAxML_bipartitions.eratoRelatives.withSRA.chr2.max0.5N.minDP3.chr2inv.subset.bs100.renamed")
plot(tree.era)

tree.era.df <- tibble(orig= tree.era$tip.label, sp=substr(tree.era$tip.label, 0, 3), 
       country=if_else(str_extract(tree.era$tip.label, "([:upper:]){1,}")=="E", "Ecuador", if_else(str_extract(tree.era$tip.label, "([:upper:]){1,}")=="C", "Colombia", "")),
       side=c( "East", "East","East","East","East", "West", "West","West","West", "East","East", "", "", "", "East", "West", "" ),
       alt.type=c("", "low distant", "high", "low distant", "high","low distant","high","high", "low distant","low distant","low distant","","","","","",""),
       haplotype= c("", "hom. WT", "hom. WT", "hom. WT", "hom. WT", "", "", "", "", "hom. INV", "hom. INV","","","","","" ,""),
       species=c("H. himera", "H. erato", "H. erato", "H. erato", "H. erato", "H. erato", "H. erato", "H. erato", "H. erato", "H. erato", "H. erato", 
                 "H. hermathena", "H. telesiphe", "H. hortense",  "H. clysnymus","H. clysonymus", "H. eleuchia"),
       tip.label=paste( species, country, side, haplotype)); tree.era.df 

tree.era.df[is.na(tree.era.df$country),]$country <- ""
tree.era.df$tip.label <- paste( tree.era.df$species, tree.era.df$country, tree.era.df$side, tree.era.df$haplotype)
tree.era.df[tree.era.df$tip.label=="H. clysnymus Ecuador East ",]$tip.label<- "H. clysnymus East "

era.full.tree <- ggtree(tree.era, ladderize =F ) %<+% tree.era.df  + 
  geom_tiplab(aes(label=tip.label), align = T, hjust = -0.05) +
  #geom_nodelab(geom = "text",aes(label = node) )+
  #geom_label2(aes(label=label, subset = !is.na(as.numeric(label)) & as.numeric(label) > 20), alpha=0.7, size=2)+
  geom_tippoint(aes(shape = alt.type, fill=alt.type),   size=3) + 
  scale_shape_manual(values=c(NA,24,21))+
  scale_fill_manual(values = c( NA, 'blue','green')) + 
  
  guides(fill = guide_legend(override.aes = list(shape = 21), title = "Type" ),
         shape = guide_legend(override.aes = list(fill = "black"), title = "Altitude"  ) ) +
  xlim(0,0.08) +
  theme(axis.title.y = element_text(size=14,vjust=-2),
        plot.margin = unit(c(0,0.2,0,0.2), "lines"),axis.text.x = element_blank(), rect = element_rect(fill = "transparent") ,#axis.text.x = element_blank(),
        axis.text.y = element_blank(),axis.title.x = element_blank(),  axis.ticks.length=unit(-0.2, "cm"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_rect(fill = "transparent", color = NA),
        panel.background = element_rect(fill = "transparent", color = NA), plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
        legend.background = element_rect(fill = "transparent", color = NA), legend.box.background = element_rect(fill = "transparent", color = NA), legend.key = element_rect(fill = "transparent", color = NA),strip.background.y = element_blank(),
        strip.text = element_text(size=14, face = "bold"), strip.background = element_blank(), legend.position = "right") ; era.full.tree



ggsave("plots/joana.tree/ch2.inv.haps.full.png", width = 7, height = 4, bg="transparent")




################################################################## PCAS all, east, west ##################################################################

############################ erato ############################
############## neutral era.E ##############
### indiv data ####
era.e.bam_names<-read.table("local/data/02.info/pop.list/era.e.txt",header=F)
era.e.bam_names$V1<-substr(era.e.bam_names$V1, 76,84); era.e.bam_names$V1
era.e.bam_names$V2<-gsub('[.rm]', '', era.e.bam_names$V1); era.e.bam_names$V2
era.e.bam_names$V2<- if_else(substr(era.e.bam_names$V2,0,2)=="CS" , gsub('[_]', '', era.e.bam_names$V2), era.e.bam_names$V2); era.e.bam_names$V2
era.e.bam_names$V2[122] <-  "BC0411_ena"
era.e.bam_names$V2

# important !!!!!!!!!!!!!!!!!!!@!!!!!! - order 
era.e.list.info <- era.e.list.info[match(era.e.bam_names$V2, era.e.list.info$id),]
era.e.bam_names$V2 ==era.e.list.info$id

### pca ####
era.e.neutral.wg.npy <-npyLoad(file="output/mds.shared/era.mel.wg.neutral/era.e.wg.cov.npy")
#era.e.neutral.wg.npy<-npyLoad(file="output/mds.shared/era.mel.chr1.neutral/era.e.chr1.cov.npy")

era.e.neutral.wg.pca<-prcomp(era.e.neutral.wg.npy)
rownames(era.e.neutral.wg.pca$x)<-era.e.bam_names$V2

# store pca info for plotting/analyses
era.e.neutral.wg.pca.df <- as.data.frame(era.e.neutral.wg.pca$x[,c(1:2)]); era.e.neutral.wg.pca.df$id <- rownames(era.e.neutral.wg.pca.df)
era.e.neutral.wg.pca.df$altitude <- subset(era.e.list.info, id %in% era.e.bam_names$V2)$altitude[match(era.e.neutral.wg.pca.df$id, subset(era.e.list.info, id %in% era.e.bam_names$V2)$id)]; era.e.neutral.wg.pca.df$altitude
era.e.neutral.wg.pca.df$alt.type <- subset(era.e.list.info, id %in% era.e.bam_names$V2)$alt.type[match(era.e.neutral.wg.pca.df$id, subset(era.e.list.info, id %in% era.e.bam_names$V2)$id)]
era.e.neutral.wg.pca.df$country<- subset(era.e.list.info, id %in% era.e.bam_names$V2)$country[match(era.e.neutral.wg.pca.df$id, subset(era.e.list.info, id %in% era.e.bam_names$V2)$id)]
era.e.neutral.wg.pca.df$alt.type.country <- paste(era.e.neutral.wg.pca.df$alt.type, era.e.neutral.wg.pca.df$country, sep=".")
era.e.neutral.wg.pca.df$region <- region

# careful when plotting, we are missing lowland colombia so must remove filled squares
summary(era.e.neutral.wg.pca)$importance[2,][1]
summary(era.e.neutral.wg.pca)$importance[2,][2]
era.e.neutral.wg.pca.plot <- ggplot(data=era.e.neutral.wg.pca.df, aes(x=PC1, y=PC2, color=altitude))+
  scale_color_gradient( low = "#00FF00B3", high = "#0000FFF2")+
  xlab("PC1 (8.8%)")+ylab("PC2 (1.8%)")+
  geom_point(aes(shape=alt.type.country), size=4)+
  scale_shape_manual(values=c(17,24,15,0,19,1))+
  #annotate(geom = 'text', label = era.e.regions.list.sig[[6]] , x =  0.2*(min(era.e.pca.list.sig.df$altitude) + max(era.e.pca.list.sig.df$altitude)), y = -Inf, hjust = 0, vjust = -1, size=4)+
  #annotate(geom = 'text', label = paste("PC1 ", round(summary(era.e.pca.list.sig[[i]])$importance[2,][1], digits = 2)*100, "%", sep="") ,  x =  0*(min(era.e.pca.list.sig.df$altitude) + max(era.e.pca.list.sig.df$altitude)), y = Inf, hjust = 0, vjust = 3, size=4)+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1),#axis.ticks.x=element_blank(),
        axis.ticks.length =unit(-0.15, "cm"), #axis.ticks.margin=unit(0.5, "cm"),
        axis.title = element_text(size=16),axis.text.y = element_text(size=10),axis.text.x = element_text(size=10),
        #axis.text.y = element_blank(),
        #axis.text.x = element_blank(), axis.title = element_blank(), axis.ticks=element_blank(),
        #axis.text.y.right = element_text(size = 10, margin = unit(c(t = 0, r = 0, b = 0, l = -8), "mm"), colour="springgreen", face="bold"),
        panel.background = element_blank(), #legend.text = element_text(size=8), legend.title =  element_text(size=12),
        panel.grid.major = element_blank(),# axis.text.x=element_blank(), #axis.text.y=element_blank(), #
        plot.margin = unit(c(0.1, 0.5, 0.1,0.5), "cm"),
        panel.grid.minor = element_blank(),# axis.title.x=element_blank(),#axis.title.y=element_blank(),
        legend.text.align = 0, legend.position = "none",strip.text = element_blank(),plot.background = element_blank(),
        strip.background =element_rect(fill="transparent", color="transparent")); era.e.neutral.wg.pca.plot

############## neutral era.w ##############
### indiv data ####
era.w.bam_names<-read.table("local/data/02.info/pop.list/era.w.txt",header=F)
era.w.bam_names$V1<-substr(era.w.bam_names$V1, 76,84); era.w.bam_names$V1
era.w.bam_names$V2<-gsub('[.rm]', '', era.w.bam_names$V1)
era.w.bam_names$V2<- if_else(substr(era.w.bam_names$V2,9, 9)=="_" , paste(era.w.bam_names$V2, "ena", sep = ""), era.w.bam_names$V2); era.w.bam_names$V2; era.w.list.info$id
era.w.bam_names$V2

# important !!!!!!!!!!!!!!!!!!!@!!!!!! - order 
era.w.list.info <- era.w.list.info[match(era.w.bam_names$V2, era.w.list.info$id),]
era.w.bam_names$V2 ==era.w.list.info$id

### pca ####
era.w.neutral.wg.npy <-npyLoad(file="output/mds.shared/era.mel.wg.neutral/era.w.wg.cov.npy")

era.w.neutral.wg.pca<-prcomp(era.w.neutral.wg.npy)
rownames(era.w.neutral.wg.pca$x)<-era.w.bam_names$V2

# store pca info for plotting/analyses
era.w.neutral.wg.pca.df <- as.data.frame(era.w.neutral.wg.pca$x[,c(1:2)]); era.w.neutral.wg.pca.df$id <- rownames(era.w.neutral.wg.pca.df)
era.w.neutral.wg.pca.df$altitude <- subset(era.w.list.info, id %in% era.w.bam_names$V2)$altitude[match(era.w.neutral.wg.pca.df$id, subset(era.w.list.info, id %in% era.w.bam_names$V2)$id)]
era.w.neutral.wg.pca.df$alt.type <- subset(era.w.list.info, id %in% era.w.bam_names$V2)$alt.type[match(era.w.neutral.wg.pca.df$id, subset(era.w.list.info, id %in% era.w.bam_names$V2)$id)]
era.w.neutral.wg.pca.df$country<- subset(era.w.list.info, id %in% era.w.bam_names$V2)$country[match(era.w.neutral.wg.pca.df$id, subset(era.w.list.info, id %in% era.w.bam_names$V2)$id)]
era.w.neutral.wg.pca.df$alt.type.country <- paste(era.w.neutral.wg.pca.df$alt.type, era.w.neutral.wg.pca.df$country, sep=".")
era.w.neutral.wg.pca.df$region <- region

# careful when plotting, we are missing lowland colombia so must remove filled squares
summary(era.w.neutral.wg.pca)$importance[2,][1]
summary(era.w.neutral.wg.pca)$importance[2,][2]
era.w.neutral.wg.pca.plot <- ggplot(data=era.w.neutral.wg.pca.df, aes(x=PC1, y=PC2, color=altitude))+
  scale_color_gradient( low = "#00FF00B3", high = "#0000FFF2")+
  xlab("PC1 (21%)")+ylab("PC2 (2%)")+
  geom_point(aes(shape=alt.type.country), size=4)+
  scale_shape_manual(values=c(17,24,0,19,1))+
  #annotate(geom = 'text', label = era.e.regions.list.sig[[6]] , x =  0.2*(min(era.e.pca.list.sig.df$altitude) + max(era.e.pca.list.sig.df$altitude)), y = -Inf, hjust = 0, vjust = -1, size=4)+
  #annotate(geom = 'text', label = paste("PC1 ", round(summary(era.e.pca.list.sig[[i]])$importance[2,][1], digits = 2)*100, "%", sep="") ,  x =  0*(min(era.e.pca.list.sig.df$altitude) + max(era.e.pca.list.sig.df$altitude)), y = Inf, hjust = 0, vjust = 3, size=4)+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1),#axis.ticks.x=element_blank(),
        axis.ticks.length =unit(-0.15, "cm"), #axis.ticks.margin=unit(0.5, "cm"),
        axis.title = element_text(size=16),axis.text.y = element_text(size=10),axis.text.x = element_text(size=10),
        #axis.text.y = element_blank(),
        #axis.text.x = element_blank(), axis.title = element_blank(), axis.ticks=element_blank(),
        #axis.text.y.right = element_text(size = 10, margin = unit(c(t = 0, r = 0, b = 0, l = -8), "mm"), colour="springgreen", face="bold"),
        panel.background = element_blank(), #legend.text = element_text(size=8), legend.title =  element_text(size=12),
        panel.grid.major = element_blank(),# axis.text.x=element_blank(), #axis.text.y=element_blank(), #
        plot.margin = unit(c(0.1, 0.5, 0.1,0.5), "cm"),plot.background = element_blank(),
        panel.grid.minor = element_blank(),# axis.title.x=element_blank(),#axis.title.y=element_blank(),
        legend.text.align = 0, legend.position = "none",strip.text = element_blank(),rect = element_rect(fill = "transparent"),
        strip.background =element_rect(fill="transparent", color="transparent")); era.w.neutral.wg.pca.plot

############## neutral era.a ##############
### indiv data ####
era.a.list.info <-subset(era.a.list.info, id!="")
era.a.bam_names<-read.table("local/data/02.info/pop.list/era.a.txt",header=F)
era.a.bam_names$V1<-substr(era.a.bam_names$V1, 76,84); era.a.bam_names$V1
era.a.bam_names$V2<- if_else(substr(era.a.bam_names$V1,  nchar(era.a.bam_names$V1), nchar(era.a.bam_names$V1))=="_" , paste(era.a.bam_names$V1, "ena", sep = ""), era.a.bam_names$V1); era.a.bam_names$V2

era.a.bam_names$V2<-gsub('[.rm]', '', era.a.bam_names$V2); era.a.bam_names$V2
era.a.bam_names$V2<- if_else(substr(era.a.bam_names$V2,  nchar(era.a.bam_names$V2), nchar(era.a.bam_names$V2))=="n" , paste(era.a.bam_names$V2, "a", sep = ""), era.a.bam_names$V2); era.a.bam_names$V2

# important !!!!!!!!!!!!!!!!!!!@!!!!!! - order 
era.a.list.info <- era.a.list.info[match(era.a.bam_names$V2, era.a.list.info$id),]
era.a.bam_names$V2 ==era.a.list.info$id
subset(era.a.bam_names, !(era.a.bam_names$V2 %in% era.a.list.info$id))


### pca ####
era.a.neutral.wg.npy <-npyLoad(file="output/mds.shared/era.mel.wg.neutral/era.a.wg.cov.npy")

era.a.neutral.wg.pca<-prcomp(era.a.neutral.wg.npy)
rownames(era.a.neutral.wg.pca$x)<-era.a.bam_names$V2

# store pca info for plotting/analyses
era.a.neutral.wg.pca.df <- as.data.frame(era.a.neutral.wg.pca$x[,c(1:2)]); era.a.neutral.wg.pca.df$id <- rownames(era.a.neutral.wg.pca.df)
era.a.neutral.wg.pca.df$altitude <- subset(era.a.list.info, id %in% era.a.bam_names$V2)$altitude[match(era.a.neutral.wg.pca.df$id, subset(era.a.list.info, id %in% era.a.bam_names$V2)$id)]
era.a.neutral.wg.pca.df$alt.type <- subset(era.a.list.info, id %in% era.a.bam_names$V2)$alt.type[match(era.a.neutral.wg.pca.df$id, subset(era.a.list.info, id %in% era.a.bam_names$V2)$id)]
era.a.neutral.wg.pca.df$country<- subset(era.a.list.info, id %in% era.a.bam_names$V2)$country[match(era.a.neutral.wg.pca.df$id, subset(era.a.list.info, id %in% era.a.bam_names$V2)$id)]
era.a.neutral.wg.pca.df$subsp<- subset(era.a.list.info, id %in% era.a.bam_names$V2)$subsp[match(era.a.neutral.wg.pca.df$id, subset(era.a.list.info, id %in% era.a.bam_names$V2)$id)]

era.a.neutral.wg.pca.df$alt.type.country <- paste(era.a.neutral.wg.pca.df$alt.type, era.a.neutral.wg.pca.df$country, sep=".")
era.a.neutral.wg.pca.df$region <- region

# careful when plotting, we are missing lowland colombia so must remove filled squares
summary(era.a.neutral.wg.pca)$importance[2,][1]
summary(era.a.neutral.wg.pca)$importance[2,][2]
era.a.neutral.wg.pca.plot <- ggplot(data=era.a.neutral.wg.pca.df, aes(x=PC1, y=PC2, color=altitude))+
  scale_color_gradient( low = "#00FF00B3", high = "#0000FFF2")+
  xlab("PC1 (94%)")+ylab("PC2 (0.6%)")+
  geom_point(aes(shape=alt.type.country), size=4)+
  scale_shape_manual(values=c(17,24,15,0,19,1))+
  scale_x_reverse()+
  #annotate(geom = 'text', label = era.e.regions.list.sig[[6]] , x =  0.2*(min(era.e.pca.list.sig.df$altitude) + max(era.e.pca.list.sig.df$altitude)), y = -Inf, hjust = 0, vjust = -1, size=4)+
  #annotate(geom = 'text', label = paste("PC1 ", round(summary(era.e.pca.list.sig[[i]])$importance[2,][1], digits = 2)*100, "%", sep="") ,  x =  0*(min(era.e.pca.list.sig.df$altitude) + max(era.e.pca.list.sig.df$altitude)), y = Inf, hjust = 0, vjust = 3, size=4)+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1),#axis.ticks.x=element_blank(),
        axis.ticks.length =unit(-0.15, "cm"), #axis.ticks.margin=unit(0.5, "cm"),
        axis.title = element_text(size=16),axis.text.y = element_text(size=10),axis.text.x = element_text(size=10),
        #axis.text.y = element_blank(),
        #axis.text.x = element_blank(), axis.title = element_blank(), axis.ticks=element_blank(),
        #axis.text.y.right = element_text(size = 10, margin = unit(c(t = 0, r = 0, b = 0, l = -8), "mm"), colour="springgreen", face="bold"),
        panel.background = element_blank(), #legend.text = element_text(size=8), legend.title =  element_text(size=12),
        panel.grid.major = element_blank(),# axis.text.x=element_blank(), #axis.text.y=element_blank(), #
        plot.margin = unit(c(0.1, 0.5, 0.1,0.5), "cm"),plot.background = element_blank(),
        panel.grid.minor = element_blank(),# axis.title.x=element_blank(),#axis.title.y=element_blank(),
        legend.text.align = 0, legend.position = "none",strip.text = element_blank(),rect = element_rect(fill = "transparent"),
        strip.background =element_rect(fill="transparent", color="transparent")); era.a.neutral.wg.pca.plot




############################ melpomene ############################
############## neutral mel.E ##############
### indiv data ####
mel.e.bam_names<-read.table("local/data/02.info/pop.list/mel.e.txt",header=F)
mel.e.bam_names$V1<-substr(mel.e.bam_names$V1, 76,84); mel.e.bam_names$V1

mel.e.bam_names$V2<- if_else(substr(mel.e.bam_names$V1,  nchar(mel.e.bam_names$V1), nchar(mel.e.bam_names$V1))=="_" , paste(mel.e.bam_names$V1, "ena", sep = ""), mel.e.bam_names$V1); mel.e.bam_names$V2
mel.e.bam_names$V2<-gsub('[.rm]', '', mel.e.bam_names$V2); mel.e.bam_names$V2
mel.e.bam_names$V2<- if_else(substr(mel.e.bam_names$V2,  0,5)=="CAM00" , paste(mel.e.bam_names$V2, "ena", sep = "_"), mel.e.bam_names$V2); mel.e.bam_names$V2
mel.e.bam_names$V2<- if_else(substr(mel.e.bam_names$V2,  0,5)=="CAM01" , paste(mel.e.bam_names$V2, "ena", sep = "_"), mel.e.bam_names$V2); mel.e.bam_names$V2


# important !!!!!!!!!!!!!!!!!!!@!!!!!! - order 
mel.e.list.info <- mel.e.list.info[match(mel.e.bam_names$V2, mel.e.list.info$id),]
mel.e.bam_names$V2 ==mel.e.list.info$id

### pca ####
mel.e.neutral.wg.npy <-npyLoad(file="output/mds.shared/era.mel.wg.neutral/mel.e.wg.cov.npy")

mel.e.neutral.wg.pca<-prcomp(mel.e.neutral.wg.npy)
rownames(mel.e.neutral.wg.pca$x)<-mel.e.bam_names$V2

# store pca info for plotting/analyses
mel.e.neutral.wg.pca.df <- as.data.frame(mel.e.neutral.wg.pca$x[,c(1:2)]); mel.e.neutral.wg.pca.df$id <- rownames(mel.e.neutral.wg.pca.df)
mel.e.neutral.wg.pca.df$altitude <- subset(mel.e.list.info, id %in% mel.e.bam_names$V2)$altitude[match(mel.e.neutral.wg.pca.df$id, subset(mel.e.list.info, id %in% mel.e.bam_names$V2)$id)]; mel.e.neutral.wg.pca.df$altitude
mel.e.neutral.wg.pca.df$alt.type <- subset(mel.e.list.info, id %in% mel.e.bam_names$V2)$alt.type[match(mel.e.neutral.wg.pca.df$id, subset(mel.e.list.info, id %in% mel.e.bam_names$V2)$id)]
mel.e.neutral.wg.pca.df$country<- subset(mel.e.list.info, id %in% mel.e.bam_names$V2)$country[match(mel.e.neutral.wg.pca.df$id, subset(mel.e.list.info, id %in% mel.e.bam_names$V2)$id)]
mel.e.neutral.wg.pca.df$alt.type.country <- paste(mel.e.neutral.wg.pca.df$alt.type, mel.e.neutral.wg.pca.df$country, sep=".")
mel.e.neutral.wg.pca.df$region <- region

# careful when plotting, we are missing lowland colombia so must remove filled squares
summary(mel.e.neutral.wg.pca)$importance[2,][1]
summary(mel.e.neutral.wg.pca)$importance[2,][2]
mel.e.neutral.wg.pca.plot <- ggplot(data=mel.e.neutral.wg.pca.df, aes(x=PC1, y=PC2, color=altitude))+
  scale_color_gradient( low = "#00FF00B3", high = "#0000FFF2")+
  xlab("PC1 (11%)")+ylab("PC2 (3.9%)")+
  geom_point(aes(shape=alt.type.country), size=4)+
  scale_shape_manual(values=c(17,24,15,0,19,1))+
  #annotate(geom = 'text', label = mel.e.regions.list.sig[[6]] , x =  0.2*(min(mel.e.pca.list.sig.df$altitude) + max(mel.e.pca.list.sig.df$altitude)), y = -Inf, hjust = 0, vjust = -1, size=4)+
  #annotate(geom = 'text', label = paste("PC1 ", round(summary(mel.e.pca.list.sig[[i]])$importance[2,][1], digits = 2)*100, "%", sep="") ,  x =  0*(min(mel.e.pca.list.sig.df$altitude) + max(mel.e.pca.list.sig.df$altitude)), y = Inf, hjust = 0, vjust = 3, size=4)+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1),#axis.ticks.x=element_blank(),
        axis.ticks.length =unit(-0.15, "cm"), #axis.ticks.margin=unit(0.5, "cm"),
        axis.title = element_text(size=16),axis.text.y = element_text(size=10),axis.text.x = element_text(size=10),
        #axis.text.y = element_blank(),
        #axis.text.x = element_blank(), axis.title = element_blank(), axis.ticks=element_blank(),
        #axis.text.y.right = element_text(size = 10, margin = unit(c(t = 0, r = 0, b = 0, l = -8), "mm"), colour="springgreen", face="bold"),
        panel.background = element_blank(), #legend.text = element_text(size=8), legend.title =  element_text(size=12),
        panel.grid.major = element_blank(),# axis.text.x=element_blank(), #axis.text.y=element_blank(), #
        plot.margin = unit(c(0.1, 0.5, 0.1,0.5), "cm"),
        panel.grid.minor = element_blank(),# axis.title.x=element_blank(),#axis.title.y=element_blank(),
        legend.text.align = 0, legend.position = "none",strip.text = element_blank(),plot.background = element_blank(),
        strip.background =element_rect(fill="transparent", color="transparent")); mel.e.neutral.wg.pca.plot

############## neutral mel.w ##############
### indiv data ####
mel.w.bam_names<-read.table("local/data/02.info/pop.list/mel.w.txt",header=F)
mel.w.bam_names$V1<-substr(mel.w.bam_names$V1, 76,84); mel.w.bam_names$V1

mel.w.bam_names$V2<- if_else(substr(mel.w.bam_names$V1,  nchar(mel.w.bam_names$V1), nchar(mel.w.bam_names$V1))=="_" , paste(mel.w.bam_names$V1, "ena", sep = ""), mel.w.bam_names$V1); mel.w.bam_names$V2
mel.w.bam_names$V2<-gsub('[.rm]', '', mel.w.bam_names$V2); mel.w.bam_names$V2
mel.w.bam_names$V2<- if_else(substr(mel.w.bam_names$V2,  0,5)=="CAM00" , paste(mel.w.bam_names$V2, "ena", sep = "_"), mel.w.bam_names$V2); mel.w.bam_names$V2
mel.w.bam_names$V2<- if_else(substr(mel.w.bam_names$V2,  0,5)=="CAM01" , paste(mel.w.bam_names$V2, "ena", sep = "_"), mel.w.bam_names$V2); mel.w.bam_names$V2


# important !!!!!!!!!!!!!!!!!!!@!!!!!! - order 
mel.w.list.info <- mel.w.list.info[match(mel.w.bam_names$V2, mel.w.list.info$id),]
mel.w.bam_names$V2 ==mel.w.list.info$id

### pca ####
mel.w.neutral.wg.npy <-npyLoad(file="output/mds.shared/era.mel.wg.neutral/mel.w.wg.cov.npy")

mel.w.neutral.wg.pca<-prcomp(mel.w.neutral.wg.npy)
rownames(mel.w.neutral.wg.pca$x)<-mel.w.bam_names$V2

# store pca info for plotting/analyses
mel.w.neutral.wg.pca.df <- as.data.frame(mel.w.neutral.wg.pca$x[,c(1:2)]); mel.w.neutral.wg.pca.df$id <- rownames(mel.w.neutral.wg.pca.df)
mel.w.neutral.wg.pca.df$altitude <- subset(mel.w.list.info, id %in% mel.w.bam_names$V2)$altitude[match(mel.w.neutral.wg.pca.df$id, subset(mel.w.list.info, id %in% mel.w.bam_names$V2)$id)]
mel.w.neutral.wg.pca.df$alt.type <- subset(mel.w.list.info, id %in% mel.w.bam_names$V2)$alt.type[match(mel.w.neutral.wg.pca.df$id, subset(mel.w.list.info, id %in% mel.w.bam_names$V2)$id)]
mel.w.neutral.wg.pca.df$country<- subset(mel.w.list.info, id %in% mel.w.bam_names$V2)$country[match(mel.w.neutral.wg.pca.df$id, subset(mel.w.list.info, id %in% mel.w.bam_names$V2)$id)]
mel.w.neutral.wg.pca.df$alt.type.country <- paste(mel.w.neutral.wg.pca.df$alt.type, mel.w.neutral.wg.pca.df$country, sep=".")
mel.w.neutral.wg.pca.df$region <- region

# careful when plotting, we are missing lowland colombia so must remove filled squares
summary(mel.w.neutral.wg.pca)$importance[2,][1]
summary(mel.w.neutral.wg.pca)$importance[2,][2]
mel.w.neutral.wg.pca.plot <- ggplot(data=mel.w.neutral.wg.pca.df, aes(x=PC1, y=PC2, color=altitude))+
  scale_color_gradient( low = "#00FF00B3", high = "#0000FFF2")+
  xlab("PC1 (36%)")+ylab("PC2 (7%)")+
  geom_point(aes(shape=alt.type.country), size=4)+
  scale_shape_manual(values=c(17,24,15,0,19) )+
  #annotate(geom = 'text', label = mel.e.regions.list.sig[[6]] , x =  0.2*(min(mel.e.pca.list.sig.df$altitude) + max(mel.e.pca.list.sig.df$altitude)), y = -Inf, hjust = 0, vjust = -1, size=4)+
  #annotate(geom = 'text', label = paste("PC1 ", round(summary(mel.e.pca.list.sig[[i]])$importance[2,][1], digits = 2)*100, "%", sep="") ,  x =  0*(min(mel.e.pca.list.sig.df$altitude) + max(mel.e.pca.list.sig.df$altitude)), y = Inf, hjust = 0, vjust = 3, size=4)+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1),#axis.ticks.x=element_blank(),
        axis.ticks.length =unit(-0.15, "cm"), #axis.ticks.margin=unit(0.5, "cm"),
        axis.title = element_text(size=16),axis.text.y = element_text(size=10),axis.text.x = element_text(size=10),
        #axis.text.y = element_blank(),
        #axis.text.x = element_blank(), axis.title = element_blank(), axis.ticks=element_blank(),
        #axis.text.y.right = element_text(size = 10, margin = unit(c(t = 0, r = 0, b = 0, l = -8), "mm"), colour="springgreen", face="bold"),
        panel.background = element_blank(), #legend.text = element_text(size=8), legend.title =  element_text(size=12),
        panel.grid.major = element_blank(),# axis.text.x=element_blank(), #axis.text.y=element_blank(), #
        plot.margin = unit(c(0.1, 0.5, 0.1,0.5), "cm"),plot.background = element_blank(),
        panel.grid.minor = element_blank(),# axis.title.x=element_blank(),#axis.title.y=element_blank(),
        legend.text.align = 0, legend.position = "none",strip.text = element_blank(),rect = element_rect(fill = "transparent"),
        strip.background =element_rect(fill="transparent", color="transparent")); mel.w.neutral.wg.pca.plot

############## neutral mel.a ##############
### indiv data ####
mel.a.list.info <- subset(all.bams.list.info,to.use.pbs=="yes" &species=="melpomene" )
mel.a.list.info <-subset(mel.a.list.info, id!=""); mel.a.list.info$id
mel.a.bam_names<-read.table("local/data/02.info/pop.list/mel.a.txt",header=F)
mel.a.bam_names$V1<-substr(mel.a.bam_names$V1, 76,84); mel.a.bam_names$V1

mel.a.bam_names$V2<- if_else(substr(mel.a.bam_names$V1,  nchar(mel.a.bam_names$V1), nchar(mel.a.bam_names$V1))=="_" , paste(mel.a.bam_names$V1, "ena", sep = ""), mel.a.bam_names$V1); mel.a.bam_names$V2
mel.a.bam_names$V2<-gsub('[.rm]', '', mel.a.bam_names$V2); mel.a.bam_names$V2
mel.a.bam_names$V2<- if_else(substr(mel.a.bam_names$V2,  0,5)=="CAM00" , paste(mel.a.bam_names$V2, "ena", sep = "_"), mel.a.bam_names$V2); mel.a.bam_names$V2
mel.a.bam_names$V2<- if_else(substr(mel.a.bam_names$V2,  0,5)=="CAM01" , paste(mel.a.bam_names$V2, "ena", sep = "_"), mel.a.bam_names$V2); mel.a.bam_names$V2

# important !!!!!!!!!!!!!!!!!!!@!!!!!! - order 
mel.a.list.info <- mel.a.list.info[match(mel.a.bam_names$V2, mel.a.list.info$id),]
mel.a.bam_names$V2 ==mel.a.list.info$id
subset(mel.a.bam_names, !(mel.a.bam_names$V2 %in% mel.a.list.info$id))


### pca ####
mel.a.neutral.wg.npy <-npyLoad(file="output/mds.shared/era.mel.wg.neutral/mel.a.wg.cov.npy")

mel.a.neutral.wg.pca<-prcomp(mel.a.neutral.wg.npy)
rownames(mel.a.neutral.wg.pca$x)<-mel.a.bam_names$V2

# store pca info for plotting/analyses
mel.a.neutral.wg.pca.df <- as.data.frame(mel.a.neutral.wg.pca$x[,c(1:2)]); mel.a.neutral.wg.pca.df$id <- rownames(mel.a.neutral.wg.pca.df)
mel.a.neutral.wg.pca.df$altitude <- subset(mel.a.list.info, id %in% mel.a.bam_names$V2)$altitude[match(mel.a.neutral.wg.pca.df$id, subset(mel.a.list.info, id %in% mel.a.bam_names$V2)$id)]
mel.a.neutral.wg.pca.df$alt.type <- subset(mel.a.list.info, id %in% mel.a.bam_names$V2)$alt.type[match(mel.a.neutral.wg.pca.df$id, subset(mel.a.list.info, id %in% mel.a.bam_names$V2)$id)]
mel.a.neutral.wg.pca.df$country<- subset(mel.a.list.info, id %in% mel.a.bam_names$V2)$country[match(mel.a.neutral.wg.pca.df$id, subset(mel.a.list.info, id %in% mel.a.bam_names$V2)$id)]
mel.a.neutral.wg.pca.df$subsp<- subset(mel.a.list.info, id %in% mel.a.bam_names$V2)$subsp[match(mel.a.neutral.wg.pca.df$id, subset(mel.a.list.info, id %in% mel.a.bam_names$V2)$id)]
mel.a.neutral.wg.pca.df$side<- subset(mel.a.list.info, id %in% mel.a.bam_names$V2)$side[match(mel.a.neutral.wg.pca.df$id, subset(mel.a.list.info, id %in% mel.a.bam_names$V2)$id)]
mel.a.neutral.wg.pca.df$alt.type.country <- paste(mel.a.neutral.wg.pca.df$alt.type, mel.a.neutral.wg.pca.df$country, sep=".")
mel.a.neutral.wg.pca.df$region <- region

# careful when plotting, we are missing lowland colombia so must remove filled squares
summary(mel.a.neutral.wg.pca)$importance[2,][1]
summary(mel.a.neutral.wg.pca)$importance[2,][2]
mel.a.neutral.wg.pca.plot <- ggplot(data=mel.a.neutral.wg.pca.df, aes(x=PC1, y=PC2, color=altitude))+
  scale_color_gradient( low = "#00FF00B3", high = "#0000FFF2")+
  xlab("PC1 (88%)")+ylab("PC2 (2.9%)")+
  geom_point(aes(shape=alt.type.country), size=4)+
  scale_shape_manual(values=c(17,24,15,0,19,1))+
  #scale_x_reverse()+
  #annotate(geom = 'text', label = mel.e.regions.list.sig[[6]] , x =  0.2*(min(mel.e.pca.list.sig.df$altitude) + max(mel.e.pca.list.sig.df$altitude)), y = -Inf, hjust = 0, vjust = -1, size=4)+
  #annotate(geom = 'text', label = paste("PC1 ", round(summary(mel.e.pca.list.sig[[i]])$importance[2,][1], digits = 2)*100, "%", sep="") ,  x =  0*(min(mel.e.pca.list.sig.df$altitude) + max(mel.e.pca.list.sig.df$altitude)), y = Inf, hjust = 0, vjust = 3, size=4)+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1),#axis.ticks.x=element_blank(),
        axis.ticks.length =unit(-0.15, "cm"), #axis.ticks.margin=unit(0.5, "cm"),
        axis.title = element_text(size=16),axis.text.y = element_text(size=10),axis.text.x = element_text(size=10),
        #axis.text.y = element_blank(),
        #axis.text.x = element_blank(), axis.title = element_blank(), axis.ticks=element_blank(),
        #axis.text.y.right = element_text(size = 10, margin = unit(c(t = 0, r = 0, b = 0, l = -8), "mm"), colour="springgreen", face="bold"),
        panel.background = element_blank(), #legend.text = element_text(size=8), legend.title =  element_text(size=12),
        panel.grid.major = element_blank(),# axis.text.x=element_blank(), #axis.text.y=element_blank(), #
        plot.margin = unit(c(0.1, 0.5, 0.1,0.5), "cm"),plot.background = element_blank(),
        panel.grid.minor = element_blank(),# axis.title.x=element_blank(),#axis.title.y=element_blank(),
        legend.text.align = 0, legend.position = "none",strip.text = element_blank(),rect = element_rect(fill = "transparent"),
        strip.background =element_rect(fill="transparent", color="transparent")); mel.a.neutral.wg.pca.plot



############## save all ##############

summary(lm(PC1~ altitude, data=era.e.neutral.wg.pca.df))
summary(lm(PC1~ altitude, data=era.w.neutral.wg.pca.df))
summary(lm(PC1~ altitude, data=mel.e.neutral.wg.pca.df))
summary(lm(PC1~ altitude, data=mel.w.neutral.wg.pca.df))


ggsave2(plot = era.e.neutral.wg.pca.plot, filename = "local/data/local/data/local/data/local/data/22_pop.gen.paper/figures/supp.pca/era.e.neutral.wg.pca.plot.png", width = 11.5, height = 10, units = 'cm', bg="transparent")
ggsave2(plot = era.w.neutral.wg.pca.plot, filename = "local/data/local/data/local/data/local/data/22_pop.gen.paper/figures/supp.pca/era.w.neutral.wg.pca.plot.png", width = 11.5, height = 10, units = 'cm', bg="transparent")
ggsave2(plot = era.a.neutral.wg.pca.plot, filename = "local/data/local/data/local/data/local/data/22_pop.gen.paper/figures/supp.pca/era.a.neutral.wg.pca.plot.png", width = 11.5, height = 10, units = 'cm', bg="transparent")

ggsave2(plot = mel.e.neutral.wg.pca.plot, filename = "local/data/local/data/local/data/local/data/22_pop.gen.paper/figures/supp.pca/mel.e.neutral.wg.pca.plot.png", width = 11.5, height = 10, units = 'cm', bg="transparent")
ggsave2(plot = mel.w.neutral.wg.pca.plot, filename = "local/data/local/data/local/data/local/data/22_pop.gen.paper/figures/supp.pca/mel.w.neutral.wg.pca.plot.png", width = 11.5, height = 10, units = 'cm', bg="transparent")
ggsave2(plot = mel.a.neutral.wg.pca.plot, filename = "local/data/local/data/local/data/local/data/22_pop.gen.paper/figures/supp.pca/mel.a.neutral.wg.pca.plot.png", width = 11.1, height = 10, units = 'cm', bg="transparent")


################################################################## LOCAL PCA PLOTTING, stdv4 all, SI ##############################################################################
#################################### 2. erato east local PCAs ####################################
##### load indiv data #####
#add rownames
era.e.bam_names<-read.table("local/data/02.info/pop.list/era.e.txt",header=F)
era.e.bam_names$V1<-substr(era.e.bam_names$V1, 76,84); era.e.bam_names$V1
era.e.bam_names$V2<-gsub('[.rm]', '', era.e.bam_names$V1); era.e.bam_names$V2
era.e.bam_names$V2<- if_else(substr(era.e.bam_names$V2,0,2)=="CS" , gsub('[_]', '', era.e.bam_names$V2), era.e.bam_names$V2); era.e.bam_names$V2
era.e.bam_names$V2[122] <-  "BC0411_ena"
era.e.bam_names$V2

# important !!!!!!!!!!!!!!!!!!!@!!!!!! - order 
era.e.list.info <- era.e.list.info[match(era.e.bam_names$V2, era.e.list.info$id),]
era.e.bam_names$V2 ==era.e.list.info$id
# info on regions
era.e.region_info <-read.csv("local/data/02.info/sites/shared.exp50kb.era.e.std4.new.era.e.summ.csv")

##### loop #####
era.e.file.list <- list.files("output/mds.shared/era.e.stdv4", full.names = T); era.e.file.list

era.e.npy.list <- list(); era.e.pca.list <- list(); era.e.pca.plot.list <- list()
era.e.pc1.vs.alt.plot.list <- list(); era.e.regions.list <- list()
era.e.pca.list.df<- list()

for (i in 1:length(era.e.file.list)) {
  file <- era.e.file.list[[i]]
  era.e.npy.list[[i]] <-npyLoad(file=paste0(era.e.file.list[[i]]))
  names(era.e.npy.list)[[i]] <- file
  region <- names(era.e.npy.list)[[i]] %>% strsplit( "/" ) %>% sapply( tail, 1 )  %>% strsplit( "[.]" )  %>%   sapply( "[", 4:6 )  %>% paste(collapse=".")
  era.e.regions.list[[i]] <- region
  era.e.pca.list[[i]]<-prcomp(era.e.npy.list[[i]])
  rownames(era.e.pca.list[[i]]$x)<-era.e.bam_names$V2

  # store pca info for plotting/analyses
  era.e.pca.list.df[[i]] <- as.data.frame(era.e.pca.list[[i]]$x[,c(1:2)]); era.e.pca.list.df[[i]]$id <- rownames(era.e.pca.list.df[[i]])
  era.e.pca.list.df[[i]]$altitude <- subset(era.e.list.info, id %in% era.e.bam_names$V2)$altitude[match(era.e.pca.list.df[[i]]$id, subset(era.e.list.info, id %in% era.e.bam_names$V2)$id)]
  era.e.pca.list.df[[i]]$alt.type <- subset(era.e.list.info, id %in% era.e.bam_names$V2)$alt.type[match(era.e.pca.list.df[[i]]$id, subset(era.e.list.info, id %in% era.e.bam_names$V2)$id)]
  era.e.pca.list.df[[i]]$country <- subset(era.e.list.info, id %in% era.e.bam_names$V2)$country[match(era.e.pca.list.df[[i]]$id, subset(era.e.list.info, id %in% era.e.bam_names$V2)$id)]
  era.e.pca.list.df[[i]]$alt.type.country <- paste(era.e.pca.list.df[[i]]$alt.type, era.e.pca.list.df[[i]]$country, sep=".")
  era.e.pca.list.df[[i]]$region <- region
  era.e.pca.list.df[[i]]$shdr.type <- if_else(era.e.shdr.outlier[era.e.shdr.outlier$shdr.para.east.id==region,]$is.allo.hdr=="yes", "shdr.allo", "shdr.east")
  
  # save PCs 
  era.e.list.info[,c(paste(region, ".PC1", sep = ""))] <-NA; era.e.list.info[,c(paste(region, ".PC2", sep = ""))] <-NA
  era.e.list.info[,c(paste(region, ".PC1", sep = ""))] <- era.e.pca.list.df[[i]]$PC1; era.e.list.info[,c(paste(region, ".PC2", sep = ""))] <- era.e.pca.list.df[[i]]$PC2
  
  # plot pcas
  era.e.pca.plot.list[[i]] <- ggplot(data=era.e.pca.list.df[[i]], aes(x=PC1, y=PC2, colour=altitude))+ 
    geom_point(size=1, alpha=0.6 , aes(shape=alt.type.country))+
    scale_color_gradient( low = "#00FF00B3", high = "#0000FFF2")+
    # colour of border according to shdr type and personalised color scale
    ggnewscale::new_scale_color()+
    geom_rect(aes(color = shdr.type), fill="transparent", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, size=2) +
    scale_color_shdr.type()+
    scale_shape_manual(values=c(17,24,15,0,19,1))+
    #scale_color_manual(values = c("#0332FA", "#6ECE6F", "#0CFA00"))+
    annotate(geom = 'text', label = region , x =  -Inf, y = Inf, hjust = -0.1, vjust = 2, size=3)+
    annotate(geom = 'text', label = paste("PC1 ", round(summary(era.e.pca.list[[i]])$importance[2,][1], digits = 2)*100, "%", sep="") , 
             x =  -Inf, y = Inf, hjust = -0.1, vjust = 3.5, size=3)+
    annotate(geom = 'text', label = paste("PC2 ", round(summary(era.e.pca.list[[i]])$importance[2,][2], digits = 2)*100, "%", sep="") , 
             x =  -Inf, y = Inf, hjust = -0.1, vjust = 5, size=3)+
    theme_classic()+ theme(panel.border = element_rect(colour = NA, fill=NA, size=1),#axis.ticks.x=element_blank(),
                           axis.ticks.length.y =unit(-0.15, "cm"), #axis.ticks.margin=unit(0.5, "cm"),
                           axis.text.y = element_blank(),
                           axis.text.x = element_blank(), axis.title = element_blank(), axis.ticks=element_blank(),
                           #axis.text.y.right = element_text(size = 10, margin = unit(c(t = 0, r = 0, b = 0, l = -8), "mm"), colour="springgreen", face="bold"),
                           panel.background = element_blank(), #legend.text = element_text(size=8), legend.title =  element_text(size=12),
                           panel.grid.major = element_blank(),# axis.text.x=element_blank(), #axis.text.y=element_blank(), #
                           plot.margin = unit(c(0, 0, -.2,0), "cm"),
                           panel.grid.minor = element_blank(),# axis.title.x=element_blank(),#axis.title.y=element_blank(),
                           legend.text.align = 0, legend.position = "none",strip.text = element_blank(),
                           strip.background =element_rect(fill="transparent", color="transparent")); era.e.pca.plot.list[[i]]
 # era.e.pca.plot.list[[i]] <- ggExtra::ggMarginal(era.e.pca.plot.list[[i]], type = "histogram", margins = "x", groupFill = T, bins=4); era.e.pca.plot.list[[i]]
  

  era.e.pc1.vs.alt.plot.list[[i]] <- ggplot(data=era.e.pca.list.df[[i]], aes(y=PC1, x=altitude, color=altitude))+
    scale_color_gradient( low = "#00FF00B3", high = "#0000FFF2")+
    geom_point(size=1, aes(shape=alt.type.country))+
    scale_shape_manual(values=c(17,24,15,0,19,1))+
    geom_smooth(method="lm", color="black" )+
    # colour of border according to shdr type and personalised color scale
    ggnewscale::new_scale_color()+
    geom_rect(aes(color = shdr.type), fill="transparent", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, size=2) +
    scale_color_shdr.type()+
    annotate(geom = 'text', label = region , x =  0.15*(min(era.e.pca.list.df[[i]]$altitude) + max(era.e.pca.list.df[[i]]$altitude)), y = -Inf, hjust = 0, vjust = -1, size=3)+
    annotate(geom = 'text', label = paste("PC1 ", round(summary(era.e.pca.list[[i]])$importance[2,][1], digits = 2)*100, "%", sep="") , 
             x =  0*(min(era.e.pca.list.df[[i]]$altitude) + max(era.e.pca.list.df[[i]]$altitude)), y = Inf, hjust = 0, vjust = 4.5, size=3)+
    stat_cor( hjust = -0.0, vjust = 0.6, size=3, aes(label =paste(..r.label.., cut(..p..,  breaks = c(-Inf, 0.0001, 0.001, 0.01, 0.05, Inf), labels = c("'****'", "'***'", "'**'", "'*'", "'ns'")), sep = "~"))) +
    theme_classic()+ theme(panel.border = element_rect(colour = "black", fill=NA, size=1),#axis.ticks.x=element_blank(),
                           axis.ticks.length.y =unit(-0.15, "cm"), #axis.ticks.margin=unit(0.5, "cm"),
                           axis.text.y = element_blank(),
                           axis.text.x = element_blank(), axis.title = element_blank(), axis.ticks=element_blank(),
                           #axis.text.y.right = element_text(size = 10, margin = unit(c(t = 0, r = 0, b = 0, l = -8), "mm"), colour="springgreen", face="bold"),
                           panel.background = element_blank(), #legend.text = element_text(size=8), legend.title =  element_text(size=12),
                           panel.grid.major = element_blank(),# axis.text.x=element_blank(), #axis.text.y=element_blank(), #
                           plot.margin = unit(c(0, 0, -.2,0), "cm"),
                           panel.grid.minor = element_blank(),# axis.title.x=element_blank(),#axis.title.y=element_blank(),
                           legend.text.align = 0, legend.position = "none",strip.text = element_blank(),
                           strip.background =element_rect(fill="transparent", color="transparent")); era.e.pc1.vs.alt.plot.list[[i]]
  

  } 

# plot_grid(plotlist=era.e.pca.plot.list, ncol=6, greedy =T)
# ggsave2("plots/local.pca/era.e.stdv4.no.intersect.all.png",  width = 210, height = 297, units = "mm")
# 
# plot_grid(plotlist=era.e.pc1.vs.alt.plot.list, ncol=7, greedy = T)
# ggsave2("plots/local.pca/era.e.pc1.vs.alt.stdv4.no.intersect.png",  width = 210, height = 297, units = "mm")


##### check whther individuals have the same haplotyoe in first half of inversion vs the second #####
era.e.pca.plot.list[[2]]+theme_bw()+geom_label_repel(aes(label=id)) | era.e.pca.plot.list[[3]] +theme_bw()+geom_label_repel(aes(label=id))
era.e.pca.list.df[[2]]$cluster <- if_else(era.e.pca.list.df[[2]]$PC1 < 1.5, "hom.wt", if_else(era.e.pca.list.df[[2]]$PC1 > 5, "hom.inv", "het" ) )
era.e.pca.list.df[[3]]$cluster <- if_else(era.e.pca.list.df[[3]]$PC1 < 1.5, "hom.wt", if_else(era.e.pca.list.df[[3]]$PC1 > 5, "hom.inv", "het" ) )
subset(era.e.pca.list.df[[2]], !(era.e.pca.list.df[[2]]$cluster == era.e.pca.list.df[[3]]$cluster ))
subset(era.e.pca.list.df[[3]], !(era.e.pca.list.df[[3]]$cluster == era.e.pca.list.df[[2]]$cluster )) # same
era.e.pca.plot.list[[2]]+theme_bw()+geom_label_repel(aes(label=id)) | era.e.pca.plot.list[[3]] +theme_bw()+geom_label_repel(aes(label=id))



###### PC1 vs altitude coefficients ###### 
### add neutral pc1 coefficients
era.e.list.info$neutral.wg.PC1 <- era.e.neutral.wg.pca.df$PC1[match(era.e.list.info$id, era.e.neutral.wg.pca.df$id)]; era.e.list.info$neutral.wg.PC1

era.e.region_info <-read.csv("local/data/02.info/sites/shared.exp50kb.era.e.std4.new.era.e.summ.csv")

names( era.e.list.info)
era.e.region_info$pc1.alt.altitude.pval<- NA; era.e.region_info$pc1.alt.neutral.wg.PC1.pval<- NA; era.e.region_info$pc1.alt.total.r2<- NA
era.e.region_info$pc1.alt.altitude.r2 <- NA; era.e.region_info$pc1.alt.neutral.wg.PC1.r2 <- NA
era.e.region_info$PC1.var.explained <- NA; era.e.region_info$PC2.var.explained <- NA
era.e.region_info$region.r.name <- unlist(era.e.regions.list)[pmatch(era.e.region_info$shdr.para.east.id ,unlist(era.e.regions.list))]; era.e.region_info$region.r.name
era.e.region_info$is.allo.hdr <- era.e.shdr.outlier$is.allo.hdr[match(era.e.region_info$shdr.para.east.id, era.e.shdr.outlier$shdr.para.east.id)]

# check vif not >4, any shdr
car::vif(lm(era.e.list.info[,c(paste("shdr.east.001", ".PC1", sep = ""))] ~ altitude +neutral.wg.PC1, data=era.e.list.info))

for (i in 1:length(era.e.file.list)) {
  region <- names(era.e.npy.list)[[i]] %>% strsplit( "/" ) %>% sapply( tail, 1 )  %>% strsplit( "[.]" )  %>%   sapply( "[", 4:6 )  %>% paste(collapse=".")
  era.e.regions.list[[i]] <- region
  # get p.values for each predictor 
  altitude.pval <- round(summary(lm(era.e.list.info[,c(paste(region, ".PC1", sep = ""))] ~ altitude +neutral.wg.PC1, data=era.e.list.info))$coefficients[,4][2], digits = 3)
  neutral.wg.PC1.pval <- round(summary(lm(era.e.list.info[,c(paste(region, ".PC1", sep = ""))] ~ altitude +neutral.wg.PC1, data=era.e.list.info))$coefficients[,4][3], digits = 3)
  
  # get relative importance of each predictor
  total.model.r2 <- calc.relimp(lm(era.e.list.info[,c(paste(region, ".PC1", sep = ""))] ~ altitude +neutral.wg.PC1, data=era.e.list.info),rela=F)@R2 
  altitude.r2 <- calc.relimp(lm(era.e.list.info[,c(paste(region, ".PC1", sep = ""))] ~ altitude +neutral.wg.PC1, data=era.e.list.info),rela=F)@lmg["altitude"] 
  neutral.wg.PC1.r2 <- calc.relimp(lm(era.e.list.info[,c(paste(region, ".PC1", sep = ""))] ~ altitude +neutral.wg.PC1, data=era.e.list.info),rela=F)@lmg["neutral.wg.PC1"] 

  # cor.result <- cor.test(era.e.list.info[,c(paste(region, ".PC1", sep = ""))] , era.e.list.info$altitude , method = "pearson")
  # estimate <- cor.result$estimate; pval <- cor.result$p.value
  era.e.region_info[subset(era.e.region_info, region.r.name!="")$region.r.name == region,]$pc1.alt.altitude.pval<- altitude.pval
  era.e.region_info[subset(era.e.region_info, region.r.name!="")$region.r.name == region,]$pc1.alt.neutral.wg.PC1.pval<- neutral.wg.PC1.pval
  era.e.region_info[subset(era.e.region_info, region.r.name!="")$region.r.name == region,]$pc1.alt.total.r2<- total.model.r2
  era.e.region_info[subset(era.e.region_info, region.r.name!="")$region.r.name == region,]$pc1.alt.altitude.r2 <- altitude.r2
  era.e.region_info[subset(era.e.region_info, region.r.name!="")$region.r.name == region,]$pc1.alt.neutral.wg.PC1.r2 <- neutral.wg.PC1.r2
  era.e.region_info[subset(era.e.region_info, region.r.name!="")$region.r.name == region,]$PC1.var.explained <- summary(era.e.pca.list[[i]])$importance[2,][1]
  era.e.region_info[subset(era.e.region_info, region.r.name!="")$region.r.name == region,]$PC2.var.explained  <- summary(era.e.pca.list[[i]])$importance[2,][2]
  }

era.e.region_info$is.alt.sig.predictor <- if_else(era.e.region_info$pc1.alt.altitude.pval<=0.05, "yes", "no"); nrow(subset(era.e.region_info, is.alt.sig.predictor=="yes"))/nrow(era.e.region_info)*100
era.e.region_info$is.neutral.wg.sig.predictor <- if_else(era.e.region_info$pc1.alt.neutral.wg.PC1.r2<=0.05, "yes", "no"); nrow(subset(era.e.region_info, is.neutral.wg.sig.predictor=="yes"))/nrow(era.e.region_info)*100
mean(subset(era.e.region_info, is.alt.sig.predictor=="yes")$pc1.alt.altitude.r2)
mean(subset(era.e.region_info, is.neutral.wg.sig.predictor =="no")$pc1.alt.neutral.wg.PC1.r2)

###### save coefficients to era.e.shdr.outlier ###### 
names(era.e.region_info); names(era.e.shdr.outlier)
# era.e.shdr.outlier$pc1.alt.cor.coef.nointersect.sites.4stdv <- era.e.region_info$pc1.alt.cor[match(era.e.shdr.outlier$shdr.para.east.id, era.e.region_info$shdr.para.east.id)]
# era.e.shdr.outlier$pc1.alt.cor.pval.nointersect.sites.4stdv <- era.e.region_info$pc1.alt.cor.pval[match(era.e.shdr.outlier$shdr.para.east.id, era.e.region_info$shdr.para.east.id)]

era.e.shdr.outlier$pc1.alt.altitude.pval <- era.e.region_info$pc1.alt.altitude.pval[match(era.e.shdr.outlier$shdr.para.east.id, era.e.region_info$shdr.para.east.id)]
era.e.shdr.outlier$pc1.alt.neutral.wg.PC1.pval <- era.e.region_info$pc1.alt.neutral.wg.PC1.pval[match(era.e.shdr.outlier$shdr.para.east.id, era.e.region_info$shdr.para.east.id)]
era.e.shdr.outlier$pc1.alt.altitude.r2 <- era.e.region_info$pc1.alt.altitude.r2[match(era.e.shdr.outlier$shdr.para.east.id, era.e.region_info$shdr.para.east.id)]; head(era.e.shdr.outlier)
era.e.shdr.outlier$pc1.alt.neutral.wg.PC1.r2 <- era.e.region_info$pc1.alt.neutral.wg.PC1.r2[match(era.e.shdr.outlier$shdr.para.east.id, era.e.region_info$shdr.para.east.id)]; head(era.e.shdr.outlier)

era.e.shdr.outlier$PC1.var.explained.nointersect.sites.4stdv <- era.e.region_info$PC1.var.explained[match(era.e.shdr.outlier$shdr.para.east.id, era.e.region_info$shdr.para.east.id)]
era.e.shdr.outlier$PC2.var.explained.nointersect.sites.4stdv <- era.e.region_info$PC2.var.explained[match(era.e.shdr.outlier$shdr.para.east.id, era.e.region_info$shdr.para.east.id)]; era.e.shdr.outlier$PC2.var.explained.nointersect.sites.4stdv

write.csv(era.e.shdr.outlier, "data/shdr.summ/shdr.era.east.outlier.df.csv", row.names = F); names(era.e.shdr.outlier)


###### how similar are pcas across shdrs ###### 
# plot altitde r2 instead of pearson cor
N_PC=2
window_eigs <-sapply(era.e.npy.list, function (cm) cov_pca(k=N_PC, covmat=cm))
window_dist<- pc_dist(t(window_eigs), npc=N_PC)
mds_axe<-cmdscale(window_dist)
rownames(mds_axe) <- sapply(era.e.file.list, substr,44,56  ) 

## add info
era.e.mds_axes.df <- as.data.frame(mds_axe)
era.e.mds_axes.df$region <- sapply(era.e.file.list, substr, 44,56  ) ; era.e.mds_axes.df$region
era.e.mds_axes.df$region.short <- sapply(era.e.file.list, substr, 44,56  ) ; era.e.mds_axes.df$region.short 
#era.e.mds_axes.df$region.short[1:5] <- c("inv1.Herato0204", "inv2.2.Herato0211" ,"inv2.Herato0211" ,"non.inv2.Herato0215" ,"non.inv.Herato0215")
era.e.mds_axes.df$region.size.bp <- era.e.region_info$region.size[pmatch(era.e.mds_axes.df$region, era.e.region_info$shdr.para.east.id )]; era.e.mds_axes.df$region.size.bp
era.e.mds_axes.df$PC1.var.explained <- era.e.region_info$PC1.var.explained[pmatch(era.e.mds_axes.df$region, era.e.region_info$shdr.para.east.id )]; era.e.mds_axes.df$region.size.bp
#era.e.mds_axes.df$pc1.alt.cor <- era.e.region_info$pc1.alt.cor[pmatch(era.e.mds_axes.df$region, era.e.region_info$shdr.para.east.id )]; era.e.mds_axes.df$region.size.bp
era.e.mds_axes.df$pc1.alt.altitude.r2 <- era.e.region_info$pc1.alt.altitude.r2[pmatch(era.e.mds_axes.df$region, era.e.region_info$shdr.para.east.id )]; era.e.mds_axes.df$pc1.alt.altitude.r2
era.e.mds_axes.df$is.non.inv <- if_else(substr(era.e.mds_axes.df$region.short, 0,1)=="s", "no", "yes")
era.e.mds_axes.df$is.allo.hdr <- era.e.region_info$is.allo.hdr[pmatch(era.e.mds_axes.df$region, era.e.region_info$shdr.para.east.id )]; era.e.mds_axes.df$region.size.bp

## plot
era.e.all.regions.mds <- ggplot(data=subset(era.e.mds_axes.df, is.non.inv!="yes"), aes(x=V1, y=V2 , fill=pc1.alt.altitude.r2))+
  geom_label_repel(aes(label=region.short ), max.overlaps = 30, fill="transparent", color="black", size=3 ) +
  scale_color_manual(values=c("#049E73", "#D65D00"))+
  geom_point(aes(size=region.size.bp/1000, colour=is.allo.hdr), alpha=.9,  shape=21, stroke=1.5)+
  ylab("MDS2") + xlab("MDS1")+ 
  scale_size_continuous(range=c(3,9), limits = c(0,950) ,name="Shared region\nsize KB (>4stdv)")+
  scale_fill_continuous(name= expression(paste(  "Altitutde R"^{2}*"")), limits = c(0,0.5), low="white", high="black" )+
  theme_bw()+ theme(legend.position = "none"); era.e.all.regions.mds


###### save pc1 vs alt plots for significant shdrs only ######
names(era.e.region_info)
era.e.region_info.sig <- subset(era.e.region_info, pc1.alt.altitude.pval<=0.05)

# loop only with significant regions and shdrs, for SI
era.e.file.list.df <- as.data.frame(era.e.file.list)
era.e.region_info.sig$file.name <- paste("output/mds.shared/era.e.stdv4/era.e.region.", era.e.region_info.sig$shdr.para.east.id, ".cov.npy", sep=""); era.e.region_info.sig$file.name
era.e.region_info.sig$is.shdr <- if_else(substr(era.e.region_info.sig$shdr.para.east.id, 0,4)=="shdr", "yes", "no"); era.e.region_info.sig$is.shdr 
era.e.region_info.sig <- subset(era.e.region_info.sig, is.shdr=="yes"); era.e.region_info.sig 

era.e.npy.list.sig <- list(); era.e.pca.list.sig <- list()
era.e.regions.list.sig <- list(); era.e.pc1.vs.alt.plot.list.sig <- list()

for (i in 1:length(era.e.region_info.sig$file.name)) {
  file <- era.e.region_info.sig$file.name[[i]]
  era.e.npy.list.sig[[i]] <-npyLoad(file=paste0(era.e.region_info.sig$file.name[i]))
  names(era.e.npy.list.sig)[[i]] <- file
  region <- names(era.e.npy.list.sig)[[i]] %>% strsplit( "/" ) %>% sapply( tail, 1 )  %>% strsplit( "[.]" )  %>%   sapply( "[", 4:6 )  %>% paste(collapse=".")
  era.e.regions.list.sig[[i]] <- region
  era.e.pca.list.sig[[i]]<-prcomp(era.e.npy.list.sig[[i]])
  rownames(era.e.pca.list.sig[[i]]$x)<-era.e.bam_names$V2
  
  # store pca info for plotting/analyses
  era.e.pca.list.sig.df <- as.data.frame(era.e.pca.list.sig[[i]]$x[,c(1:2)]); era.e.pca.list.sig.df$id <- rownames(era.e.pca.list.sig.df)
  era.e.pca.list.sig.df$altitude <- subset(era.e.list.info, id %in% era.e.bam_names$V2)$altitude[match(era.e.pca.list.sig.df$id, subset(era.e.list.info, id %in% era.e.bam_names$V2)$id)]
  era.e.pca.list.sig.df$alt.type <- subset(era.e.list.info, id %in% era.e.bam_names$V2)$alt.type[match(era.e.pca.list.sig.df$id, subset(era.e.list.info, id %in% era.e.bam_names$V2)$id)]
  era.e.pca.list.sig.df$country <- subset(era.e.list.info, id %in% era.e.bam_names$V2)$country[match(era.e.pca.list.sig.df$id, subset(era.e.list.info, id %in% era.e.bam_names$V2)$id)]
  era.e.pca.list.sig.df$region <- region
  era.e.pca.list.sig.df$alt.type.country <- paste(era.e.pca.list.sig.df$alt.type, era.e.pca.list.sig.df$country, sep=".")
  era.e.pca.list.sig.df$is.allo <- era.e.shdr.outlier$is.allo.hdr[match(era.e.pca.list.sig.df$region, era.e.shdr.outlier$shdr.para.east.id)]
  era.e.pca.list.sig.df$shdr.type <- if_else(era.e.shdr.outlier[era.e.shdr.outlier$shdr.para.east.id==region,]$is.allo.hdr=="yes", "shdr.allo", "shdr.east")
  r2.alt.label <-  as.expression(bquote(Alt.~R^2==.( round(era.e.region_info.sig[i,]$pc1.alt.altitude.r2, digits = 2)  )))
  
  # plot
  era.e.pc1.vs.alt.plot.list.sig[[i]] <- ggplot(data=era.e.pca.list.sig.df, aes(y=PC1, x=altitude, color=altitude))+
  scale_color_gradient( low = "#00FF00B3", high = "#0000FFF2")+
  geom_point(size=1 , aes(shape=alt.type.country ))+geom_smooth(method="lm", color="black")+
  scale_shape_manual(values=c(17,24,15,0,19,1))+
    # colour of border according to shdr type and personalised color scale
    ggnewscale::new_scale_color()+
    geom_rect(aes(color = shdr.type), fill="transparent", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, size=2) +
    scale_color_shdr.type()+
  annotate(geom = 'text', label = region , x =  0.2*(min(era.e.pca.list.sig.df$altitude) + max(era.e.pca.list.sig.df$altitude)), y = -Inf, hjust = 0, vjust = -0.25, size=4)+
  annotate(geom = "text", label = paste("PC1 ", round(summary(era.e.pca.list.sig[[i]])$importance[2,][1], digits = 2)*100, "%", sep="") , 
           x =  0*(min(era.e.pca.list.sig.df$altitude) + max(era.e.pca.list.sig.df$altitude)), y = Inf, hjust = 0, vjust = 2.5, size=4)+
  stat_cor( hjust = -0.0, vjust = 0.25, size=4, aes(label =paste(..r.label.., cut(..p..,  breaks = c(-Inf, 0.0001, 0.001, 0.01, 0.05, Inf), labels = c("'****'", "'***'", "'**'", "'*'", "'ns'")), sep = "~"))) +
    annotate(geom = 'text', label =r2.alt.label , parse=TRUE, x =  0*(min(era.e.pca.list.sig.df$altitude) + max(era.e.pca.list.sig.df$altitude)), y = Inf, hjust = 0, vjust = 3, size=4)+
    scale_y_continuous(expand = c(.1,0))+
    theme_classic()+ theme(panel.border = element_rect(colour = "transparent", fill=NA, size=1),#axis.ticks.x=element_blank(),
                           axis.ticks.length.y =unit(-0.15, "cm"), #axis.ticks.margin=unit(0.5, "cm"),
                           axis.text.y = element_blank(),
                           axis.text.x = element_blank(), axis.title = element_blank(), axis.ticks=element_blank(),
                           #axis.text.y.right = element_text(size = 10, margin = unit(c(t = 0, r = 0, b = 0, l = -8), "mm"), colour="springgreen", face="bold"),
                           panel.background = element_blank(), #legend.text = element_text(size=8), legend.title =  element_text(size=12),
                           panel.grid.major = element_blank(),# axis.text.x=element_blank(), #axis.text.y=element_blank(), #
                           plot.margin = unit(c(0, 0, -.05,0), "cm"),
                           panel.grid.minor = element_blank(),# axis.title.x=element_blank(),#axis.title.y=element_blank(),
                           legend.text.align = 0, legend.position = "none",strip.text = element_blank(),
                           strip.background =element_rect(fill="transparent", color="transparent")); era.e.pc1.vs.alt.plot.list.sig[[i]]
}

plot_grid(plotlist=era.e.pc1.vs.alt.plot.list.sig, ncol=6, greedy = T)
ggsave2("plots/local.pca/era.e.pc1.vs.alt.significant.no.intersect.png",  width = 210, height = 227, units = "mm")

## get legends for SI

ggplot(data=era.e.pca.list.sig.df, aes(y=PC1, x=altitude, color=altitude))+
  scale_color_gradient( low = "#00FF00B3", high = "#0000FFF2")+
  geom_point(size=3 , aes(shape=country))+geom_smooth(method="lm", color="black")+
  annotate(geom = 'text', label = region , x =  0.2*(min(era.e.pca.list.sig.df$altitude) + max(era.e.pca.list.sig.df$altitude)), y = -Inf, hjust = 0, vjust = -1, size=4)+
  annotate(geom = 'text', label = paste("PC1 ", round(summary(era.e.pca.list.sig[[i]])$importance[2,][1], digits = 2)*100, "%", sep="") , 
           x =  0*(min(era.e.pca.list.sig.df$altitude) + max(era.e.pca.list.sig.df$altitude)), y = Inf, hjust = 0, vjust = 3, size=4)+
  stat_cor( hjust = -0.0, vjust = 0.8, size=4, aes(label =paste(..r.label.., cut(..p..,  breaks = c(-Inf, 0.0001, 0.001, 0.01, 0.05, Inf), labels = c("'****'", "'***'", "'**'", "'*'", "'ns'")), sep = "~")))+theme_classic()+
  theme(legend.text = element_text(size=18), legend.title =  element_text(size=18))

############ get plots for fig 4 ############
### shdr 12, number 6 ####
era.e.pc1.vs.alt.plot.list.sig[[6]]
era.e.pca.list.sig.df <- as.data.frame(era.e.pca.list.sig[[6]]$x[,c(1:2)]); era.e.pca.list.sig.df$id <- rownames(era.e.pca.list.sig.df)
era.e.pca.list.sig.df$altitude <- subset(era.e.list.info, id %in% era.e.bam_names$V2)$altitude[match(era.e.pca.list.sig.df$id, subset(era.e.list.info, id %in% era.e.bam_names$V2)$id)]
era.e.pca.list.sig.df$alt.type <- subset(era.e.list.info, id %in% era.e.bam_names$V2)$alt.type[match(era.e.pca.list.sig.df$id, subset(era.e.list.info, id %in% era.e.bam_names$V2)$id)]
era.e.pca.list.sig.df$country <- subset(era.e.list.info, id %in% era.e.bam_names$V2)$country[match(era.e.pca.list.sig.df$id, subset(era.e.list.info, id %in% era.e.bam_names$V2)$id)]
era.e.pca.list.sig.df$region <- era.e.regions.list.sig[[6]]
era.e.pca.list.sig.df$alt.type.country <- paste(era.e.pca.list.sig.df$alt.type, era.e.pca.list.sig.df$country, sep=".")
unique(era.e.pca.list.sig.df$alt.type.country)

ggplot(data=era.e.pca.list.sig.df, aes(y=PC1, x=altitude, color=altitude))+
  scale_color_gradient( low = "#00FF00B3", high = "#0000FFF2")+
  ylab("Local PCA- PC1 (93%)")+xlab("Altitude")+
  geom_point(aes(shape=alt.type.country), size=4)+
  geom_smooth(method="lm", color="black")+
  scale_shape_manual(values=c(17,24,15,0,19,1))+
  #annotate(geom = 'text', label = era.e.regions.list.sig[[6]] , x =  0.2*(min(era.e.pca.list.sig.df$altitude) + max(era.e.pca.list.sig.df$altitude)), y = -Inf, hjust = 0, vjust = -1, size=4)+
  #annotate(geom = 'text', label = paste("PC1 ", round(summary(era.e.pca.list.sig[[i]])$importance[2,][1], digits = 2)*100, "%", sep="") ,  x =  0*(min(era.e.pca.list.sig.df$altitude) + max(era.e.pca.list.sig.df$altitude)), y = Inf, hjust = 0, vjust = 3, size=4)+
  stat_cor( hjust = -0.75, vjust = 0.8, size=6, aes(label =paste(..r.label.., cut(..p..,  breaks = c(-Inf, 0.0001, 0.001, 0.01, 0.05, Inf), labels = c("'****'", "'***'", "'**'", "'*'", "'ns'")), sep = "~"))) +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1),#axis.ticks.x=element_blank(),
        axis.ticks.length.y =unit(-0.15, "cm"), #axis.ticks.margin=unit(0.5, "cm"),
        axis.title = element_text(size=16),axis.text.y = element_text(size=10),axis.text.x = element_text(size=12),
        #axis.text.y = element_blank(),
        #axis.text.x = element_blank(), axis.title = element_blank(), axis.ticks=element_blank(),
        #axis.text.y.right = element_text(size = 10, margin = unit(c(t = 0, r = 0, b = 0, l = -8), "mm"), colour="springgreen", face="bold"),
        panel.background = element_blank(), #legend.text = element_text(size=8), legend.title =  element_text(size=12),
        panel.grid.major = element_blank(),# axis.text.x=element_blank(), #axis.text.y=element_blank(), #
        plot.margin = unit(c(0.01, 0.1, 0.1,0.1), "cm"),
        panel.grid.minor = element_blank(),# axis.title.x=element_blank(),#axis.title.y=element_blank(),
        legend.text.align = 0, legend.position = "none",strip.text = element_blank(),
        strip.background =element_rect(fill="transparent", color="transparent"))
ggsave("local/data/local/data/local/data/local/data/22_pop.gen.paper/figures/fig4.inv.local.pca/shdr.12.e.png", width = 9, height = 9, units = 'cm')

ggplot(data=era.e.pca.list.sig.df, aes(y=PC1, x=altitude, color=altitude))+
  scale_color_gradient( low = "#00FF00B3", high = "#0000FFF2")+
  ylab("Local PCA- PC1 (93%)")+xlab("Altitude")+
  geom_point(aes(shape=alt.type.country), size=4)+
  geom_smooth(method="lm", color="black")+
  scale_shape_manual(values=c(17,24,15,0,19,1))+theme_classic()
ggsave("local/data/local/data/local/data/local/data/22_pop.gen.paper/figures/fig4.inv.local.pca/legend.shdr.12.e.png", width = 9, height = 9, units = 'cm')


### inv2.2, number 2 ####
era.e.pc1.vs.alt.plot.list[[2]]
era.e.pca.list.df <- as.data.frame(era.e.pca.list[[2]]$x[,c(1:2)]); era.e.pca.list.df$id <- rownames(era.e.pca.list.df)
era.e.pca.list.df$altitude <- subset(era.e.list.info, id %in% era.e.bam_names$V2)$altitude[match(era.e.pca.list.df$id, subset(era.e.list.info, id %in% era.e.bam_names$V2)$id)]
era.e.pca.list.df$alt.type <- subset(era.e.list.info, id %in% era.e.bam_names$V2)$alt.type[match(era.e.pca.list.df$id, subset(era.e.list.info, id %in% era.e.bam_names$V2)$id)]
era.e.pca.list.df$country <- subset(era.e.list.info, id %in% era.e.bam_names$V2)$country[match(era.e.pca.list.df$id, subset(era.e.list.info, id %in% era.e.bam_names$V2)$id)]
era.e.pca.list.df$region <- era.e.regions.list[[2]]
era.e.pca.list.df$alt.type.country <- paste(era.e.pca.list.df$alt.type, era.e.pca.list.df$country, sep=".")
unique(era.e.pca.list.df$alt.type.country)
round(summary(era.e.pca.list[[2]])$importance[2,][1], digits = 2)*100
ggplot(data=era.e.pca.list.df, aes(y=PC1, x=altitude, color=altitude))+
  scale_color_gradient( low = "#00FF00B3", high = "#0000FFF2")+
  ylab("Local PCA- PC1 (88%)")+xlab("Altitude")+
  geom_point(aes(shape=alt.type.country), size=4)+
  geom_smooth(method="lm", color="black")+
  scale_shape_manual(values=c(17,24,15,0,19,1))+
  #annotate(geom = 'text', label = era.e.regions.list[[6]] , x =  0.2*(min(era.e.pca.list.df$altitude) + max(era.e.pca.list.df$altitude)), y = -Inf, hjust = 0, vjust = -1, size=4)+
  #annotate(geom = 'text', label = paste("PC1 ", round(summary(era.e.pca.list[[i]])$importance[2,][1], digits = 2)*100, "%", sep="") ,  x =  0*(min(era.e.pca.list.df$altitude) + max(era.e.pca.list.df$altitude)), y = Inf, hjust = 0, vjust = 3, size=4)+
  stat_cor( hjust = -0.75, vjust = 0.8, size=6, aes(label =paste(..r.label.., cut(..p..,  breaks = c(-Inf, 0.0001, 0.001, 0.01, 0.05, Inf), labels = c("'****'", "'***'", "'**'", "'*'", "'ns'")), sep = "~"))) +
   theme(panel.border = element_rect(colour = "black", fill=NA, size=1),#axis.ticks.x=element_blank(),
                         axis.ticks.length.y =unit(-0.15, "cm"), #axis.ticks.margin=unit(0.5, "cm"),
                         axis.title = element_text(size=16),axis.text = element_text(size=12),
                         #axis.text.y = element_blank(),
                         #axis.text.x = element_blank(), axis.title = element_blank(), axis.ticks=element_blank(),
                         #axis.text.y.right = element_text(size = 10, margin = unit(c(t = 0, r = 0, b = 0, l = -8), "mm"), colour="springgreen", face="bold"),
                         panel.background = element_blank(), #legend.text = element_text(size=8), legend.title =  element_text(size=12),
                         panel.grid.major = element_blank(),# axis.text.x=element_blank(), #axis.text.y=element_blank(), #
                         plot.margin = unit(c(0.01, 0.1, 0.1,0.1), "cm"),
                         panel.grid.minor = element_blank(),# axis.title.x=element_blank(),#axis.title.y=element_blank(),
                         legend.text.align = 0, legend.position = "none",strip.text = element_blank(),
                         strip.background =element_rect(fill="transparent", color="transparent"))

ggsave("local/data/local/data/local/data/local/data/22_pop.gen.paper/figures/fig4.inv.local.pca/inv2.2.e.png", width = 9, height = 9, units = 'cm')

 ### non inverted region, number 4 ####
era.e.pc1.vs.alt.plot.list[[4]]
era.e.pca.list.df <- as.data.frame(era.e.pca.list[[4]]$x[,c(1:2)]); era.e.pca.list.df$id <- rownames(era.e.pca.list.df)
era.e.pca.list.df$altitude <- subset(era.e.list.info, id %in% era.e.bam_names$V2)$altitude[match(era.e.pca.list.df$id, subset(era.e.list.info, id %in% era.e.bam_names$V2)$id)]
era.e.pca.list.df$alt.type <- subset(era.e.list.info, id %in% era.e.bam_names$V2)$alt.type[match(era.e.pca.list.df$id, subset(era.e.list.info, id %in% era.e.bam_names$V2)$id)]
era.e.pca.list.df$country <- subset(era.e.list.info, id %in% era.e.bam_names$V2)$country[match(era.e.pca.list.df$id, subset(era.e.list.info, id %in% era.e.bam_names$V2)$id)]
era.e.pca.list.df$region <- era.e.regions.list[[4]]
era.e.pca.list.df$alt.type.country <- paste(era.e.pca.list.df$alt.type, era.e.pca.list.df$country, sep=".")
unique(era.e.pca.list.df$alt.type.country)
round(summary(era.e.pca.list[[4]])$importance[4,][1], digits = 2)*100

ggplot(data=era.e.pca.list.df, aes(y=PC1, x=altitude, color=altitude))+
  scale_color_gradient( low = "#00FF00B3", high = "#0000FFF2")+
  ylab("Local PCA- PC1 (7%)")+xlab("Altitude")+
  geom_point(aes(shape=alt.type.country), size=4)+
  geom_smooth(method="lm", color="black")+
  scale_shape_manual(values=c(17,24,15,0,19,1))+
  #annotate(geom = 'text', label = era.e.regions.list[[6]] , x =  0.2*(min(era.e.pca.list.df$altitude) + max(era.e.pca.list.df$altitude)), y = -Inf, hjust = 0, vjust = -1, size=4)+
  #annotate(geom = 'text', label = paste("PC1 ", round(summary(era.e.pca.list[[i]])$importance[2,][1], digits = 2)*100, "%", sep="") ,  x =  0*(min(era.e.pca.list.df$altitude) + max(era.e.pca.list.df$altitude)), y = Inf, hjust = 0, vjust = 3, size=4)+
  stat_cor( hjust = -0.75, vjust = 0.8, size=6, aes(label =paste(..r.label.., cut(..p..,  breaks = c(-Inf, 0.0001, 0.001, 0.01, 0.05, Inf), labels = c("'****'", "'***'", "'**'", "'*'", "'ns'")), sep = "~"))) +
 scale_y_continuous(breaks = c(0,1))+
   theme(panel.border = element_rect(colour = "black", fill=NA, size=1),#axis.ticks.x=element_blank(),
        axis.ticks.length.y =unit(-0.15, "cm"), #axis.ticks.margin=unit(0.5, "cm"),
        axis.title = element_text(size=16),axis.text = element_text(size=12),
        #axis.text.y = element_blank(),
        #axis.text.x = element_blank(), axis.title = element_blank(), axis.ticks=element_blank(),
        #axis.text.y.right = element_text(size = 10, margin = unit(c(t = 0, r = 0, b = 0, l = -8), "mm"), colour="springgreen", face="bold"),
        panel.background = element_blank(), #legend.text = element_text(size=8), legend.title =  element_text(size=12),
        panel.grid.major = element_blank(),# axis.text.x=element_blank(), #axis.text.y=element_blank(), #
        plot.margin = unit(c(0.01, 0.1, 0.1,0.1), "cm"),
        panel.grid.minor = element_blank(),# axis.title.x=element_blank(),#axis.title.y=element_blank(),
        legend.text.align = 0, legend.position = "none",strip.text = element_blank(),
        strip.background =element_rect(fill="transparent", color="transparent"))
ggsave("local/data/local/data/local/data/local/data/22_pop.gen.paper/figures/fig4.inv.local.pca/non.inv2.e.png", width = 9, height = 9, units = 'cm')



###### save PCA plots  shdrs within putative inversions only ######



###### local PCA trees os inversion ############## 
library(ggtree)
tree.era <- read.tree(file = "data/local.pca.tree/SHDR/RAxML_bipartitions.eratoRelatives.withSRA.chr1.max0.5N.minDP3.balSubset.shdr.east.004.bs100.renamed.nwk")
tree.era <- read.tree(file = "data/local.pca.tree/SHDR/RAxML_bipartitions.eratoRelatives.withSRA.chr1.max0.5N.minDP3.balSubset.shdr.east.005.bs100.renamed.nwk")
tree.era <- read.tree(file = "data/local.pca.tree/SHDR/RAxML_bipartitions.eratoRelatives.withSRA.chr2.max0.5N.minDP3.balSubset.shdr.east.009.bs100.renamed.nwk")
tree.era <- read.tree(file = "data/local.pca.tree/SHDR/RAxML_bipartitions.eratoRelatives.withSRA.chr2.max0.5N.minDP3.max0.5N.shdr10.renamed.nwk")

tree.era <- read.tree(file = "data/local.pca.tree/SHDR/RAxML_bipartitions.eratoRelatives.withSRA.chr17.max0.5N.minDP3.balSubset.shdr.east.056.bs100.renamed.renamed.nwk")
tree.era <- read.tree(file = "data/local.pca.tree/SHDR/RAxML_bipartitions.eratoRelatives.withSRA.chr21.max0.5N.minDP3.balSubset.shdr.east.077.bs100.renamed.nwk")
tree.era <- read.tree(file = "data/local.pca.tree/SHDR/RAxML_bipartitions.eratoRelatives.withSRA.chr3.max0.5N.minDP3.balSubset.shdr.east.017.bs100.renamed.nwk")

tree.era <- read.tree(file = "data/local.pca.tree/SHDR/RAxML_bipartitions.eratoRelatives.withSRA.chr3.max0.5N.minDP3.balSubset.shdr.east.017.bs100.renamed.nwk")

tree.era <- read.tree(file = "data/local.pca.tree/RAxML_bipartitions.eratoRelatives.withSRA.chr2.max0.5N.minDP3.chr2inv.subset.bs100.renamed")
plot(tree.era)

#################################### 3 erato west local PCA new ############## 
##### load indiv data #####
#add rownames
era.w.bam_names<-read.table("local/data/02.info/pop.list/era.w.txt",header=F)
era.w.bam_names$V1<-substr(era.w.bam_names$V1, 76,84); era.w.bam_names$V1
era.w.bam_names$V2<-gsub('[.rm]', '', era.w.bam_names$V1)
era.w.bam_names$V2<- if_else(substr(era.w.bam_names$V2,9, 9)=="_" , paste(era.w.bam_names$V2, "ena", sep = ""), era.w.bam_names$V2); era.w.bam_names$V2; era.w.list.info$id


# important !!!!!!!!!!!!!!!!!!!@!!!!!! - order 
era.w.list.info <- era.w.list.info[match(era.w.bam_names$V2, era.w.list.info$id),]
era.w.bam_names$V2 ==era.w.list.info$id
# info on regions
era.w.region_info <-read.csv("local/data/02.info/sites/shared.exp50kb.era.w.std4.new.era.w.summ.csv")

##### loop #####

era.w.file.list <- list.files("output/mds.shared/era.w.stdv4/", full.names = T); era.w.file.list

era.w.npy.list <- list(); era.w.pca.list <- list(); era.w.pca.plot.list <- list()
era.w.pc1.vs.alt.plot.list <- list(); era.w.regions.list <- list()
era.w.pca.list.df<- list()
for (i in 1:length(era.w.file.list)) {
  file <- era.w.file.list[[i]]
  era.w.npy.list[[i]] <-npyLoad(file=paste0(era.w.file.list[[i]]))
  names(era.w.npy.list)[[i]] <- file
  region <- names(era.w.npy.list)[[i]] %>% strsplit( "/" ) %>% sapply( tail, 1 )  %>% strsplit( "[.]" )  %>%   sapply( "[", 4:6 )  %>% paste(collapse=".")
  era.w.regions.list[[i]] <- region
  era.w.pca.list[[i]]<-prcomp(era.w.npy.list[[i]])
  rownames(era.w.pca.list[[i]]$x)<-era.w.bam_names$V2
  
  # store pca info for plotting/analyses
  era.w.pca.list.df[[i]] <- as.data.frame(era.w.pca.list[[i]]$x[,c(1:2)]); era.w.pca.list.df[[i]]$id <- rownames(era.w.pca.list.df[[i]])
  era.w.pca.list.df[[i]]$altitude <- subset(era.w.list.info, id %in% era.w.bam_names$V2)$altitude[match(era.w.pca.list.df[[i]]$id, subset(era.w.list.info, id %in% era.w.bam_names$V2)$id)]
  era.w.pca.list.df[[i]]$alt.type <- subset(era.w.list.info, id %in% era.w.bam_names$V2)$alt.type[match(era.w.pca.list.df[[i]]$id, subset(era.w.list.info, id %in% era.w.bam_names$V2)$id)]
  era.w.pca.list.df[[i]]$country <- subset(era.w.list.info, id %in% era.w.bam_names$V2)$country[match(era.w.pca.list.df[[i]]$id, subset(era.w.list.info, id %in% era.w.bam_names$V2)$id)]
  era.w.pca.list.df[[i]]$alt.type.country <- paste(era.w.pca.list.df[[i]]$alt.type, era.w.pca.list.df[[i]]$country, sep=".")
  era.w.pca.list.df[[i]]$region <- region
  era.w.pca.list.df[[i]]$shdr.type <- if_else(era.w.shdr.outlier[era.w.shdr.outlier$shdr.para.west.id==region,]$is.allo.hdr=="yes", "shdr.allo", "shdr.west")
  
  # save PCs 
  era.w.list.info[,c(paste(region, ".PC1", sep = ""))] <-NA; era.w.list.info[,c(paste(region, ".PC2", sep = ""))] <-NA
  era.w.list.info[,c(paste(region, ".PC1", sep = ""))] <- era.w.pca.list.df[[i]]$PC1; era.w.list.info[,c(paste(region, ".PC2", sep = ""))] <- era.w.pca.list.df[[i]]$PC2
  
  # plot pcas
  era.w.pca.plot.list[[i]] <- ggplot(data=era.w.pca.list.df[[i]], aes(x=PC1, y=PC2, colour=altitude))+ 
    geom_point(size=1, alpha=0.6 , aes(shape=alt.type.country))+
    scale_color_gradient( low = "#00FF00B3", high = "#0000FFF2")+
    # colour of border according to shdr type and personalised color scale
    ggnewscale::new_scale_color()+
    geom_rect(aes(color = shdr.type), fill="transparent", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, size=2) +
    scale_color_shdr.type()+
    scale_shape_manual(values=c(17,24,0,19,1))+
    #scale_color_manual(values = c("#0332FA", "#6ECE6F", "#0CFA00"))+
    annotate(geom = 'text', label = region , x =  -Inf, y = Inf, hjust = -0.1, vjust = 2, size=3)+
    annotate(geom = 'text', label = paste("PC1 ", round(summary(era.w.pca.list[[i]])$importance[2,][1], digits = 2)*100, "%", sep="") , 
             x =  -Inf, y = Inf, hjust = -0.1, vjust = 3.5, size=3)+
    annotate(geom = 'text', label = paste("PC2 ", round(summary(era.w.pca.list[[i]])$importance[2,][2], digits = 2)*100, "%", sep="") , 
             x =  -Inf, y = Inf, hjust = -0.1, vjust = 5, size=3)+
    theme_classic()+ theme(panel.border = element_rect(colour = NA, fill=NA, size=1),#axis.ticks.x=element_blank(),
                           axis.ticks.length.y =unit(-0.15, "cm"), #axis.ticks.margin=unit(0.5, "cm"),
                           axis.text.y = element_blank(),
                           axis.text.x = element_blank(), axis.title = element_blank(), axis.ticks=element_blank(),
                           #axis.text.y.right = element_text(size = 10, margin = unit(c(t = 0, r = 0, b = 0, l = -8), "mm"), colour="springgreen", face="bold"),
                           panel.background = element_blank(), #legend.text = element_text(size=8), legend.title =  element_text(size=12),
                           panel.grid.major = element_blank(),# axis.text.x=element_blank(), #axis.text.y=element_blank(), #
                           plot.margin = unit(c(0, 0, -.2,0), "cm"),
                           panel.grid.minor = element_blank(),# axis.title.x=element_blank(),#axis.title.y=element_blank(),
                           legend.text.align = 0, legend.position = "none",strip.text = element_blank(),
                           strip.background =element_rect(fill="transparent", color="transparent")); era.w.pca.plot.list[[i]]
  # era.w.pca.plot.list[[i]] <- ggExtra::ggMarginal(era.w.pca.plot.list[[i]], type = "histogram", margins = "x", groupFill = T, bins=4); era.w.pca.plot.list[[i]]
  
  
  era.w.pc1.vs.alt.plot.list[[i]] <- ggplot(data=era.w.pca.list.df[[i]], aes(y=PC1, x=altitude, color=altitude))+
    scale_color_gradient( low = "#00FF00B3", high = "#0000FFF2")+
    geom_point(size=1, aes(shape=alt.type.country))+
    scale_shape_manual(values=c(17,24,0,19,1))+
    geom_smooth(method="lm", color="black" )+
    # colour of border according to shdr type and personalised color scale
    ggnewscale::new_scale_color()+
    geom_rect(aes(color = shdr.type), fill="transparent", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, size=2) +
    scale_color_shdr.type()+
    annotate(geom = 'text', label = region , x =  0.15*(min(era.w.pca.list.df[[i]]$altitude) + max(era.w.pca.list.df[[i]]$altitude)), y = -Inf, hjust = 0, vjust = -1, size=3)+
    annotate(geom = 'text', label = paste("PC1 ", round(summary(era.w.pca.list[[i]])$importance[2,][1], digits = 2)*100, "%", sep="") , 
             x =  0*(min(era.w.pca.list.df[[i]]$altitude) + max(era.w.pca.list.df[[i]]$altitude)), y = Inf, hjust = 0, vjust = 4.5, size=3)+
    stat_cor( hjust = -0.0, vjust = 0.6, size=3, aes(label =paste(..r.label.., cut(..p..,  breaks = c(-Inf, 0.0001, 0.001, 0.01, 0.05, Inf), labels = c("'****'", "'***'", "'**'", "'*'", "'ns'")), sep = "~"))) +
    theme_classic()+ theme(panel.border = element_rect(colour = "black", fill=NA, size=1),#axis.ticks.x=element_blank(),
                           axis.ticks.length.y =unit(-0.15, "cm"), #axis.ticks.margin=unit(0.5, "cm"),
                           axis.text.y = element_blank(),
                           axis.text.x = element_blank(), axis.title = element_blank(), axis.ticks=element_blank(),
                           #axis.text.y.right = element_text(size = 10, margin = unit(c(t = 0, r = 0, b = 0, l = -8), "mm"), colour="springgreen", face="bold"),
                           panel.background = element_blank(), #legend.text = element_text(size=8), legend.title =  element_text(size=12),
                           panel.grid.major = element_blank(),# axis.text.x=element_blank(), #axis.text.y=element_blank(), #
                           plot.margin = unit(c(0, 0, -.2,0), "cm"),
                           panel.grid.minor = element_blank(),# axis.title.x=element_blank(),#axis.title.y=element_blank(),
                           legend.text.align = 0, legend.position = "none",strip.text = element_blank(),
                           strip.background =element_rect(fill="transparent", color="transparent")); era.w.pc1.vs.alt.plot.list[[i]]
  
  
} 

# plot_grid(plotlist=era.w.pca.plot.list, ncol=6, greedy =T)
# ggsave2("plots/local.pca/era.w.stdv4.no.intersect.all.png",  width = 210, height = 550, units = "mm")
# 
# plot_grid(plotlist=era.w.pc1.vs.alt.plot.list, ncol=7, greedy = T)
# ggsave2("plots/local.pca/era.w.pc1.vs.alt.stdv4.no.intersect.png",  width = 210, height = 550, units = "mm")
# 


###### PC1 vs altitude coefficients ###### 
### add neutral pc1 coefficients
era.w.list.info$neutral.wg.PC1 <- era.w.neutral.wg.pca.df$PC1[match(era.w.list.info$id, era.w.neutral.wg.pca.df$id)]; era.w.list.info$neutral.wg.PC1

era.w.region_info <-read.csv("local/data/02.info/sites/shared.exp50kb.era.w.std4.new.era.w.summ.csv")

names( era.w.list.info)
era.w.region_info$pc1.alt.altitude.pval<- NA; era.w.region_info$pc1.alt.neutral.wg.PC1.pval<- NA; era.w.region_info$pc1.alt.total.r2<- NA
era.w.region_info$pc1.alt.altitude.r2 <- NA; era.w.region_info$pc1.alt.neutral.wg.PC1.r2 <- NA
era.w.region_info$PC1.var.explained <- NA; era.w.region_info$PC2.var.explained <- NA
era.w.region_info$region.r.name <- unlist(era.w.regions.list)[pmatch(era.w.region_info$shdr.para.west.id ,unlist(era.w.regions.list))]; era.w.region_info$region.r.name
era.w.region_info$is.allo.hdr <- era.w.shdr.outlier$is.allo.hdr[match(era.w.region_info$shdr.para.west.id, era.w.shdr.outlier$shdr.para.west.id)]

# check vif not >4, any shdr
car::vif(lm(era.w.list.info[,c(paste("shdr.west.001", ".PC1", sep = ""))] ~ altitude +neutral.wg.PC1, data=era.w.list.info))

for (i in 1:length(era.w.file.list)) {
  region <- names(era.w.npy.list)[[i]] %>% strsplit( "/" ) %>% sapply( tail, 1 )  %>% strsplit( "[.]" )  %>%   sapply( "[", 4:6 )  %>% paste(collapse=".")
  era.w.regions.list[[i]] <- region
  # get p.values for each predictor 
  altitude.pval <- round(summary(lm(era.w.list.info[,c(paste(region, ".PC1", sep = ""))] ~ altitude +neutral.wg.PC1, data=era.w.list.info))$coefficients[,4][2], digits = 3)
  neutral.wg.PC1.pval <- round(summary(lm(era.w.list.info[,c(paste(region, ".PC1", sep = ""))] ~ altitude +neutral.wg.PC1, data=era.w.list.info))$coefficients[,4][3], digits = 3)
  
  # get relative importance of each predictor
  total.model.r2 <- calc.relimp(lm(era.w.list.info[,c(paste(region, ".PC1", sep = ""))] ~ altitude +neutral.wg.PC1, data=era.w.list.info),rela=F)@R2 
  altitude.r2 <- calc.relimp(lm(era.w.list.info[,c(paste(region, ".PC1", sep = ""))] ~ altitude +neutral.wg.PC1, data=era.w.list.info),rela=F)@lmg["altitude"] 
  neutral.wg.PC1.r2 <- calc.relimp(lm(era.w.list.info[,c(paste(region, ".PC1", sep = ""))] ~ altitude +neutral.wg.PC1, data=era.w.list.info),rela=F)@lmg["neutral.wg.PC1"] 
  
  # cor.result <- cor.test(era.w.list.info[,c(paste(region, ".PC1", sep = ""))] , era.w.list.info$altitude , method = "pearson")
  # estimate <- cor.result$estimate; pval <- cor.result$p.value
  era.w.region_info[subset(era.w.region_info, region.r.name!="")$region.r.name == region,]$pc1.alt.altitude.pval<- altitude.pval
  era.w.region_info[subset(era.w.region_info, region.r.name!="")$region.r.name == region,]$pc1.alt.neutral.wg.PC1.pval<- neutral.wg.PC1.pval
  era.w.region_info[subset(era.w.region_info, region.r.name!="")$region.r.name == region,]$pc1.alt.total.r2<- total.model.r2
  era.w.region_info[subset(era.w.region_info, region.r.name!="")$region.r.name == region,]$pc1.alt.altitude.r2 <- altitude.r2
  era.w.region_info[subset(era.w.region_info, region.r.name!="")$region.r.name == region,]$pc1.alt.neutral.wg.PC1.r2 <- neutral.wg.PC1.r2
  era.w.region_info[subset(era.w.region_info, region.r.name!="")$region.r.name == region,]$PC1.var.explained <- summary(era.w.pca.list[[i]])$importance[2,][1]
  era.w.region_info[subset(era.w.region_info, region.r.name!="")$region.r.name == region,]$PC2.var.explained  <- summary(era.w.pca.list[[i]])$importance[2,][2]
}

era.w.region_info$is.alt.sig.predictor <- if_else(era.w.region_info$pc1.alt.altitude.pval<=0.05, "yes", "no"); nrow(subset(era.w.region_info, is.alt.sig.predictor=="yes"))/nrow(era.w.region_info)*100
era.w.region_info$is.neutral.wg.sig.predictor <- if_else(era.w.region_info$pc1.alt.neutral.wg.PC1.r2<=0.05, "yes", "no"); nrow(subset(era.w.region_info, is.neutral.wg.sig.predictor=="yes"))/nrow(era.w.region_info)*100
mean(subset(era.w.region_info, is.alt.sig.predictor=="yes")$pc1.alt.altitude.r2)
mean(subset(era.w.region_info, is.neutral.wg.sig.predictor =="no")$pc1.alt.neutral.wg.PC1.r2)

###### save coefficients to era.w.shdr.outlier ###### 
names(era.w.region_info); names(era.w.shdr.outlier)
# era.w.shdr.outlier$pc1.alt.cor.coef.nointersect.sites.4stdv <- era.w.region_info$pc1.alt.cor[match(era.w.shdr.outlier$shdr.para.west.id, era.w.region_info$shdr.para.west.id)]
# era.w.shdr.outlier$pc1.alt.cor.pval.nointersect.sites.4stdv <- era.w.region_info$pc1.alt.cor.pval[match(era.w.shdr.outlier$shdr.para.west.id, era.w.region_info$shdr.para.west.id)]

era.w.shdr.outlier$pc1.alt.altitude.pval <- era.w.region_info$pc1.alt.altitude.pval[match(era.w.shdr.outlier$shdr.para.west.id, era.w.region_info$shdr.para.west.id)]
era.w.shdr.outlier$pc1.alt.neutral.wg.PC1.pval <- era.w.region_info$pc1.alt.neutral.wg.PC1.pval[match(era.w.shdr.outlier$shdr.para.west.id, era.w.region_info$shdr.para.west.id)]
era.w.shdr.outlier$pc1.alt.altitude.r2 <- era.w.region_info$pc1.alt.altitude.r2[match(era.w.shdr.outlier$shdr.para.west.id, era.w.region_info$shdr.para.west.id)]; head(era.w.shdr.outlier)
era.w.shdr.outlier$pc1.alt.neutral.wg.PC1.r2 <- era.w.region_info$pc1.alt.neutral.wg.PC1.r2[match(era.w.shdr.outlier$shdr.para.west.id, era.w.region_info$shdr.para.west.id)]; head(era.w.shdr.outlier)

era.w.shdr.outlier$PC1.var.explained.nointersect.sites.4stdv <- era.w.region_info$PC1.var.explained[match(era.w.shdr.outlier$shdr.para.west.id, era.w.region_info$shdr.para.west.id)]
era.w.shdr.outlier$PC2.var.explained.nointersect.sites.4stdv <- era.w.region_info$PC2.var.explained[match(era.w.shdr.outlier$shdr.para.west.id, era.w.region_info$shdr.para.west.id)]; era.w.shdr.outlier$PC2.var.explained.nointersect.sites.4stdv

#write.csv(era.w.shdr.outlier, "data/shdr.summ/shdr.era.west.outlier.df.csv", row.names = F); names(era.w.shdr.outlier)


###### how similar are pcas across shdrs ###### 
# plot altitde r2 instead of pearson cor
N_PC=2
window_eigs <-sapply(era.w.npy.list, function (cm) cov_pca(k=N_PC, covmat=cm))
window_dist<- pc_dist(t(window_eigs), npc=N_PC)
mds_axe<-cmdscale(window_dist)
rownames(mds_axe) <- sapply(era.w.file.list, substr,45,57  ) 

## add info
era.w.mds_axes.df <- as.data.frame(mds_axe)
era.w.mds_axes.df$region <- sapply(era.w.file.list, substr, 45,57  ) ; era.w.mds_axes.df$region
era.w.mds_axes.df$region.short <- sapply(era.w.file.list, substr, 45,57  ) ; era.w.mds_axes.df$region.short 
#era.w.mds_axes.df$region.short[1:5] <- c("inv1.Herato0204", "inv2.2.Herato0211" ,"inv2.Herato0211" ,"non.inv2.Herato0215" ,"non.inv.Herato0215")
era.w.mds_axes.df$region.size.bp <- era.w.region_info$region.size[pmatch(era.w.mds_axes.df$region, era.w.region_info$shdr.para.west.id )]; era.w.mds_axes.df$region.size.bp
era.w.mds_axes.df$PC1.var.explained <- era.w.region_info$PC1.var.explained[pmatch(era.w.mds_axes.df$region, era.w.region_info$shdr.para.west.id )]; era.w.mds_axes.df$region.size.bp
#era.w.mds_axes.df$pc1.alt.cor <- era.w.region_info$pc1.alt.cor[pmatch(era.w.mds_axes.df$region, era.w.region_info$shdr.para.west.id )]; era.w.mds_axes.df$region.size.bp
era.w.mds_axes.df$pc1.alt.altitude.r2 <- era.w.region_info$pc1.alt.altitude.r2[pmatch(era.w.mds_axes.df$region, era.w.region_info$shdr.para.west.id )]; era.w.mds_axes.df$pc1.alt.altitude.r2
era.w.mds_axes.df$is.non.inv <- if_else(substr(era.w.mds_axes.df$region.short, 0,1)=="s", "no", "yes")
era.w.mds_axes.df$is.allo.hdr <- era.w.region_info$is.allo.hdr[pmatch(era.w.mds_axes.df$region, era.w.region_info$shdr.para.west.id )]; era.w.mds_axes.df$region.size.bp

## plot
era.w.all.regions.mds <- ggplot(data=subset(era.w.mds_axes.df, is.non.inv!="yes"), aes(x=V1, y=V2 , fill=pc1.alt.altitude.r2))+
  geom_label_repel(aes(label=region.short ), max.overlaps = 30, fill="transparent", color="black", size=3 ) +
  scale_color_manual(values=c( "#0372B2", "#D65D00"))+
  geom_point(aes(size=region.size.bp/1000, colour=is.allo.hdr), alpha=.9,  shape=21, stroke=1.5)+
  ylab("MDS2") + xlab("MDS1")+ 
  scale_size_continuous(range=c(3,9), limits = c(0,950) ,name="Shared region\nsize KB (>4stdv)")+
  scale_fill_continuous(name= expression(paste(  "Altitutde R"^{2}*"")), limits = c(0,0.5), low="white", high="black" )+
  theme_bw()+ theme(legend.position = "none"); era.w.all.regions.mds


###### save pc1 vs alt plots for significant shdrs only ######
names(era.w.region_info)
era.w.region_info.sig <- subset(era.w.region_info, pc1.alt.altitude.pval<=0.05)

# loop only with significant regions and shdrs, for SI
era.w.file.list.df <- as.data.frame(era.w.file.list)
era.w.region_info.sig$file.name <- paste("output/mds.shared/era.w.stdv4/era.w.region.", era.w.region_info.sig$shdr.para.west.id, ".cov.npy", sep=""); era.w.region_info.sig$file.name
era.w.region_info.sig$is.shdr <- if_else(substr(era.w.region_info.sig$shdr.para.west.id, 0,4)=="shdr", "yes", "no"); era.w.region_info.sig$is.shdr 
era.w.region_info.sig <- subset(era.w.region_info.sig, is.shdr=="yes"); era.w.region_info.sig 

era.w.npy.list.sig <- list(); era.w.pca.list.sig <- list()
era.w.regions.list.sig <- list(); era.w.pc1.vs.alt.plot.list.sig <- list()

for (i in 1:length(era.w.region_info.sig$file.name)) {
  file <- era.w.region_info.sig$file.name[[i]]
  era.w.npy.list.sig[[i]] <-npyLoad(file=paste0(era.w.region_info.sig$file.name[i]))
  names(era.w.npy.list.sig)[[i]] <- file
  region <- names(era.w.npy.list.sig)[[i]] %>% strsplit( "/" ) %>% sapply( tail, 1 )  %>% strsplit( "[.]" )  %>%   sapply( "[", 4:6 )  %>% paste(collapse=".")
  era.w.regions.list.sig[[i]] <- region
  era.w.pca.list.sig[[i]]<-prcomp(era.w.npy.list.sig[[i]])
  rownames(era.w.pca.list.sig[[i]]$x)<-era.w.bam_names$V2
  
  # store pca info for plotting/analyses
  era.w.pca.list.sig.df <- as.data.frame(era.w.pca.list.sig[[i]]$x[,c(1:2)]); era.w.pca.list.sig.df$id <- rownames(era.w.pca.list.sig.df)
  era.w.pca.list.sig.df$altitude <- subset(era.w.list.info, id %in% era.w.bam_names$V2)$altitude[match(era.w.pca.list.sig.df$id, subset(era.w.list.info, id %in% era.w.bam_names$V2)$id)]
  era.w.pca.list.sig.df$alt.type <- subset(era.w.list.info, id %in% era.w.bam_names$V2)$alt.type[match(era.w.pca.list.sig.df$id, subset(era.w.list.info, id %in% era.w.bam_names$V2)$id)]
  era.w.pca.list.sig.df$country <- subset(era.w.list.info, id %in% era.w.bam_names$V2)$country[match(era.w.pca.list.sig.df$id, subset(era.w.list.info, id %in% era.w.bam_names$V2)$id)]
  era.w.pca.list.sig.df$region <- region
  era.w.pca.list.sig.df$alt.type.country <- paste(era.w.pca.list.sig.df$alt.type, era.w.pca.list.sig.df$country, sep=".")
  era.w.pca.list.sig.df$is.allo <- era.w.shdr.outlier$is.allo.hdr[match(era.w.pca.list.sig.df$region, era.w.shdr.outlier$shdr.para.west.id)]
  era.w.pca.list.sig.df$shdr.type <- if_else(era.w.shdr.outlier[era.w.shdr.outlier$shdr.para.west.id==region,]$is.allo.hdr=="yes", "shdr.allo", "shdr.west")
  r2.alt.label <-  as.expression(bquote(Alt.~R^2==.( round(era.w.region_info.sig[i,]$pc1.alt.altitude.r2, digits = 2)  )))
  
  # plot
  era.w.pc1.vs.alt.plot.list.sig[[i]] <- ggplot(data=era.w.pca.list.sig.df, aes(y=PC1, x=altitude, color=altitude))+
    scale_color_gradient( low = "#00FF00B3", high = "#0000FFF2")+
    geom_point(size=1 , aes(shape=alt.type.country ))+geom_smooth(method="lm", color="black")+
    scale_shape_manual(values=c(17,24,0,19,1))+
    # colour of border according to shdr type and personalised color scale
    ggnewscale::new_scale_color()+
    geom_rect(aes(color = shdr.type), fill="transparent", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, size=2) +
    scale_color_shdr.type()+
    annotate(geom = 'text', label = region , x =  0.2*(min(era.w.pca.list.sig.df$altitude) + max(era.w.pca.list.sig.df$altitude)), y = -Inf, hjust = 0, vjust = -0.25, size=4)+
    annotate(geom = "text", label = paste("PC1 ", round(summary(era.w.pca.list.sig[[i]])$importance[2,][1], digits = 2)*100, "%", sep="") , 
             x =  0*(min(era.w.pca.list.sig.df$altitude) + max(era.w.pca.list.sig.df$altitude)), y = Inf, hjust = 0, vjust = 2.5, size=4)+
    stat_cor( hjust = -0.0, vjust = 0.25, size=4, aes(label =paste(..r.label.., cut(..p..,  breaks = c(-Inf, 0.0001, 0.001, 0.01, 0.05, Inf), labels = c("'****'", "'***'", "'**'", "'*'", "'ns'")), sep = "~"))) +
    annotate(geom = 'text', label =r2.alt.label , parse=TRUE, x =  0*(min(era.w.pca.list.sig.df$altitude) + max(era.w.pca.list.sig.df$altitude)), y = Inf, hjust = 0, vjust = 3, size=4)+
    scale_y_continuous(expand = c(.1,0))+
    theme_classic()+ theme(panel.border = element_rect(colour = "transparent", fill=NA, size=1),#axis.ticks.x=element_blank(),
                           axis.ticks.length.y =unit(-0.15, "cm"), #axis.ticks.margin=unit(0.5, "cm"),
                           axis.text.y = element_blank(),
                           axis.text.x = element_blank(), axis.title = element_blank(), axis.ticks=element_blank(),
                           #axis.text.y.right = element_text(size = 10, margin = unit(c(t = 0, r = 0, b = 0, l = -8), "mm"), colour="springgreen", face="bold"),
                           panel.background = element_blank(), #legend.text = element_text(size=8), legend.title =  element_text(size=12),
                           panel.grid.major = element_blank(),# axis.text.x=element_blank(), #axis.text.y=element_blank(), #
                           plot.margin = unit(c(0, 0, -.05,0), "cm"),
                           panel.grid.minor = element_blank(),# axis.title.x=element_blank(),#axis.title.y=element_blank(),
                           legend.text.align = 0, legend.position = "none",strip.text = element_blank(),
                           strip.background =element_rect(fill="transparent", color="transparent")); era.w.pc1.vs.alt.plot.list.sig[[i]]
}

plot_grid(plotlist=era.w.pc1.vs.alt.plot.list.sig, ncol=6, greedy = T)
ggsave2("plots/local.pca/era.w.pc1.vs.alt.significant.no.intersect.png",  width = 210, height = 227, units = "mm")


#################################### 4. melpomene east local PCAs ####################################
##### load indiv data #####
#add rownames
mel.e.bam_names<-read.table("local/data/02.info/pop.list/mel.e.txt",header=F)
mel.e.bam_names$V1<-substr(mel.e.bam_names$V1, 76,84); mel.e.bam_names$V1
mel.e.bam_names$V2<-gsub('[.rm]', '', mel.e.bam_names$V1)
mel.e.bam_names$V2<- if_else(sapply(mel.e.bam_names$V2, str_detect, "_")==TRUE , paste(mel.e.bam_names$V2, "ena", sep = ""), mel.e.bam_names$V2); mel.e.bam_names$V2; mel.e.list.info$id

mel.e.bam_names[mel.e.bam_names$V2=="CAM016550",]$V2 <- "CAM016550_ena"
mel.e.bam_names[mel.e.bam_names$V2=="CAM017162",]$V2 <- "CAM017162_ena"

mel.e.bam_names$V2; mel.e.list.info$id


# important !!!!!!!!!!!!!!!!!!!@!!!!!! - order 
mel.e.list.info <- mel.e.list.info[match(mel.e.bam_names$V2, mel.e.list.info$id),]
mel.e.bam_names$V2 ==mel.e.list.info$id
# info on regions
mel.e.region_info <-read.csv("local/data/02.info/sites/shared.exp50kb.mel.e.std4.new.mel.e.summ.csv")

##### loop #####
mel.e.file.list <- list.files("output/mds.shared/mel.e.stdv4/", full.names = T); mel.e.file.list

mel.e.npy.list <- list(); mel.e.pca.list <- list(); mel.e.pca.plot.list <- list()
mel.e.pc1.vs.alt.plot.list <- list(); mel.e.regions.list <- list()
mel.e.pca.list.df<- list()

for (i in 1:length(mel.e.file.list)) {
  file <- mel.e.file.list[[i]]
  mel.e.npy.list[[i]] <-npyLoad(file=paste0(mel.e.file.list[[i]]))
  names(mel.e.npy.list)[[i]] <- file
  region <- names(mel.e.npy.list)[[i]] %>% strsplit( "/" ) %>% sapply( tail, 1 )  %>% strsplit( "[.]" )  %>%   sapply( "[", 4:6 )  %>% paste(collapse=".")
  mel.e.regions.list[[i]] <- region
  mel.e.pca.list[[i]]<-prcomp(mel.e.npy.list[[i]])
  rownames(mel.e.pca.list[[i]]$x)<-mel.e.bam_names$V2
  
  # store pca info for plotting/analyses
  mel.e.pca.list.df[[i]] <- as.data.frame(mel.e.pca.list[[i]]$x[,c(1:2)]); mel.e.pca.list.df[[i]]$id <- rownames(mel.e.pca.list.df[[i]])
  mel.e.pca.list.df[[i]]$altitude <- subset(mel.e.list.info, id %in% mel.e.bam_names$V2)$altitude[match(mel.e.pca.list.df[[i]]$id, subset(mel.e.list.info, id %in% mel.e.bam_names$V2)$id)]
  mel.e.pca.list.df[[i]]$alt.type <- subset(mel.e.list.info, id %in% mel.e.bam_names$V2)$alt.type[match(mel.e.pca.list.df[[i]]$id, subset(mel.e.list.info, id %in% mel.e.bam_names$V2)$id)]
  mel.e.pca.list.df[[i]]$country <- subset(mel.e.list.info, id %in% mel.e.bam_names$V2)$country[match(mel.e.pca.list.df[[i]]$id, subset(mel.e.list.info, id %in% mel.e.bam_names$V2)$id)]
  mel.e.pca.list.df[[i]]$alt.type.country <- paste(mel.e.pca.list.df[[i]]$alt.type, mel.e.pca.list.df[[i]]$country, sep=".")
  mel.e.pca.list.df[[i]]$region <- region
  mel.e.pca.list.df[[i]]$shdr.type <- if_else(mel.e.shdr.outlier[mel.e.shdr.outlier$shdr.para.east.id==region,]$is.allo.hdr=="yes", "shdr.allo", "shdr.east")
  
  # save PCs 
  mel.e.list.info[,c(paste(region, ".PC1", sep = ""))] <-NA; mel.e.list.info[,c(paste(region, ".PC2", sep = ""))] <-NA
  mel.e.list.info[,c(paste(region, ".PC1", sep = ""))] <- mel.e.pca.list.df[[i]]$PC1; mel.e.list.info[,c(paste(region, ".PC2", sep = ""))] <- mel.e.pca.list.df[[i]]$PC2
  
  # plot pcas
  mel.e.pca.plot.list[[i]] <- ggplot(data=mel.e.pca.list.df[[i]], aes(x=PC1, y=PC2, colour=altitude))+ 
    geom_point(size=1, alpha=0.6 , aes(shape=alt.type.country))+
    scale_color_gradient( low = "#00FF00B3", high = "#0000FFF2")+
    # colour of border according to shdr type and personalised color scale
    ggnewscale::new_scale_color()+
    geom_rect(aes(color = shdr.type), fill="transparent", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, size=2) +
    scale_color_shdr.type()+
    scale_shape_manual(values=c(17,24,15,0,19,1))+
    #scale_color_manual(values = c("#0332FA", "#6ECE6F", "#0CFA00"))+
    annotate(geom = 'text', label = region , x =  -Inf, y = Inf, hjust = -0.1, vjust = 2, size=3)+
    annotate(geom = 'text', label = paste("PC1 ", round(summary(mel.e.pca.list[[i]])$importance[2,][1], digits = 2)*100, "%", sep="") , 
             x =  -Inf, y = Inf, hjust = -0.1, vjust = 3.5, size=3)+
    annotate(geom = 'text', label = paste("PC2 ", round(summary(mel.e.pca.list[[i]])$importance[2,][2], digits = 2)*100, "%", sep="") , 
             x =  -Inf, y = Inf, hjust = -0.1, vjust = 5, size=3)+
    theme_classic()+ theme(panel.border = element_rect(colour = "transparent", fill=NA, size=1),#axis.ticks.x=element_blank(),
                           axis.ticks.length.y =unit(-0.15, "cm"), #axis.ticks.margin=unit(0.5, "cm"),
                           axis.text.y = element_blank(),
                           axis.text.x = element_blank(), axis.title = element_blank(), axis.ticks=element_blank(),
                           #axis.text.y.right = element_text(size = 10, margin = unit(c(t = 0, r = 0, b = 0, l = -8), "mm"), colour="springgreen", face="bold"),
                           panel.background = element_blank(), #legend.text = element_text(size=8), legend.title =  element_text(size=12),
                           panel.grid.major = element_blank(),# axis.text.x=element_blank(), #axis.text.y=element_blank(), #
                           plot.margin = unit(c(0, 0, -.2,0), "cm"),
                           panel.grid.minor = element_blank(),# axis.title.x=element_blank(),#axis.title.y=element_blank(),
                           legend.text.align = 0, legend.position = "none",strip.text = element_blank(),
                           strip.background =element_rect(fill="transparent", color="transparent")); mel.e.pca.plot.list[[i]]
  # mel.e.pca.plot.list[[i]] <- ggExtra::ggMarginal(mel.e.pca.plot.list[[i]], type = "histogram", margins = "x", groupFill = T, bins=4); mel.e.pca.plot.list[[i]]
  
  
  mel.e.pc1.vs.alt.plot.list[[i]] <- ggplot(data=mel.e.pca.list.df[[i]], aes(y=PC1, x=altitude, color=altitude))+
    scale_color_gradient( low = "#00FF00B3", high = "#0000FFF2")+
    geom_point(size=1, aes(shape=alt.type.country))+
    scale_shape_manual(values=c(17,24,15,0,19,1))+
    geom_smooth(method="lm", color="black" )+
    # colour of border according to shdr type and personalised color scale
    ggnewscale::new_scale_color()+
    geom_rect(aes(color = shdr.type), fill="transparent", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, size=2) +
    scale_color_shdr.type()+
    annotate(geom = 'text', label = region , x =  0.15*(min(mel.e.pca.list.df[[i]]$altitude) + max(mel.e.pca.list.df[[i]]$altitude)), y = -Inf, hjust = 0, vjust = -1, size=3)+
    annotate(geom = 'text', label = paste("PC1 ", round(summary(mel.e.pca.list[[i]])$importance[2,][1], digits = 2)*100, "%", sep="") , 
             x =  0*(min(mel.e.pca.list.df[[i]]$altitude) + max(mel.e.pca.list.df[[i]]$altitude)), y = Inf, hjust = 0, vjust = 4.5, size=3)+
    stat_cor( hjust = -0.0, vjust = 0.6, size=3, aes(label =paste(..r.label.., cut(..p..,  breaks = c(-Inf, 0.0001, 0.001, 0.01, 0.05, Inf), labels = c("'****'", "'***'", "'**'", "'*'", "'ns'")), sep = "~"))) +
    theme_classic()+ theme(panel.border = element_rect(colour = "transparent", fill=NA, size=1),#axis.ticks.x=element_blank(),
                           axis.ticks.length.y =unit(-0.15, "cm"), #axis.ticks.margin=unit(0.5, "cm"),
                           axis.text.y = element_blank(),
                           axis.text.x = element_blank(), axis.title = element_blank(), axis.ticks=element_blank(),
                           #axis.text.y.right = element_text(size = 10, margin = unit(c(t = 0, r = 0, b = 0, l = -8), "mm"), colour="springgreen", face="bold"),
                           panel.background = element_blank(), #legend.text = element_text(size=8), legend.title =  element_text(size=12),
                           panel.grid.major = element_blank(),# axis.text.x=element_blank(), #axis.text.y=element_blank(), #
                           plot.margin = unit(c(0, 0, -.2,0), "cm"),
                           panel.grid.minor = element_blank(),# axis.title.x=element_blank(),#axis.title.y=element_blank(),
                           legend.text.align = 0, legend.position = "none",strip.text = element_blank(),
                           strip.background =element_rect(fill="transparent", color="transparent")); mel.e.pc1.vs.alt.plot.list[[i]]
  
  
} 



plot_grid(plotlist=mel.e.pca.plot.list, ncol=6, greedy =T)
ggsave2("plots/local.pca/mel.e.stdv4.no.intersect.all.png",  width = 210, height = 297, units = "mm")

plot_grid(plotlist=mel.e.pc1.vs.alt.plot.list, ncol=7, greedy = T)
ggsave2("plots/local.pca/mel.e.pc1.vs.alt.stdv4.no.intersect.png",  width = 210, height = 247, units = "mm")


###### PC1 vs altitude coefficients ###### 
### add neutral pc1 coefficients
mel.e.list.info$neutral.wg.PC1 <- mel.e.neutral.wg.pca.df$PC1[match(mel.e.list.info$id, mel.e.neutral.wg.pca.df$id)]; mel.e.list.info$neutral.wg.PC1

mel.e.region_info <-read.csv("local/data/02.info/sites/shared.exp50kb.mel.e.std4.new.mel.e.summ.csv")

names( mel.e.list.info)
mel.e.region_info$pc1.alt.altitude.pval<- NA; mel.e.region_info$pc1.alt.neutral.wg.PC1.pval<- NA; mel.e.region_info$pc1.alt.total.r2<- NA
mel.e.region_info$pc1.alt.altitude.r2 <- NA; mel.e.region_info$pc1.alt.neutral.wg.PC1.r2 <- NA
mel.e.region_info$PC1.var.explained <- NA; mel.e.region_info$PC2.var.explained <- NA
mel.e.region_info$region.r.name <- unlist(mel.e.regions.list)[pmatch(mel.e.region_info$shdr.para.east.id ,unlist(mel.e.regions.list))]; mel.e.region_info$region.r.name
mel.e.region_info$is.allo.hdr <- mel.e.shdr.outlier$is.allo.hdr[match(mel.e.region_info$shdr.para.east.id, mel.e.shdr.outlier$shdr.para.east.id)]

# check vif not >4, any shdr
car::vif(lm(mel.e.list.info[,c(paste("shdr.east.001", ".PC1", sep = ""))] ~ altitude +neutral.wg.PC1, data=mel.e.list.info))

for (i in 1:length(mel.e.file.list)) {
  region <- names(mel.e.npy.list)[[i]] %>% strsplit( "/" ) %>% sapply( tail, 1 )  %>% strsplit( "[.]" )  %>%   sapply( "[", 4:6 )  %>% paste(collapse=".")
  mel.e.regions.list[[i]] <- region
  # get p.values for each predictor 
  altitude.pval <- round(summary(lm(mel.e.list.info[,c(paste(region, ".PC1", sep = ""))] ~ altitude +neutral.wg.PC1, data=mel.e.list.info))$coefficients[,4][2], digits = 3)
  neutral.wg.PC1.pval <- round(summary(lm(mel.e.list.info[,c(paste(region, ".PC1", sep = ""))] ~ altitude +neutral.wg.PC1, data=mel.e.list.info))$coefficients[,4][3], digits = 3)
  
  # get relative importance of each predictor
  total.model.r2 <- calc.relimp(lm(mel.e.list.info[,c(paste(region, ".PC1", sep = ""))] ~ altitude +neutral.wg.PC1, data=mel.e.list.info),rela=F)@R2 
  altitude.r2 <- calc.relimp(lm(mel.e.list.info[,c(paste(region, ".PC1", sep = ""))] ~ altitude +neutral.wg.PC1, data=mel.e.list.info),rela=F)@lmg["altitude"] 
  neutral.wg.PC1.r2 <- calc.relimp(lm(mel.e.list.info[,c(paste(region, ".PC1", sep = ""))] ~ altitude +neutral.wg.PC1, data=mel.e.list.info),rela=F)@lmg["neutral.wg.PC1"] 
  
  # cor.result <- cor.test(mel.e.list.info[,c(paste(region, ".PC1", sep = ""))] , mel.e.list.info$altitude , method = "pearson")
  # estimate <- cor.result$estimate; pval <- cor.result$p.value
  mel.e.region_info[subset(mel.e.region_info, region.r.name!="")$region.r.name == region,]$pc1.alt.altitude.pval<- altitude.pval
  mel.e.region_info[subset(mel.e.region_info, region.r.name!="")$region.r.name == region,]$pc1.alt.neutral.wg.PC1.pval<- neutral.wg.PC1.pval
  mel.e.region_info[subset(mel.e.region_info, region.r.name!="")$region.r.name == region,]$pc1.alt.total.r2<- total.model.r2
  mel.e.region_info[subset(mel.e.region_info, region.r.name!="")$region.r.name == region,]$pc1.alt.altitude.r2 <- altitude.r2
  mel.e.region_info[subset(mel.e.region_info, region.r.name!="")$region.r.name == region,]$pc1.alt.neutral.wg.PC1.r2 <- neutral.wg.PC1.r2
  mel.e.region_info[subset(mel.e.region_info, region.r.name!="")$region.r.name == region,]$PC1.var.explained <- summary(mel.e.pca.list[[i]])$importance[2,][1]
  mel.e.region_info[subset(mel.e.region_info, region.r.name!="")$region.r.name == region,]$PC2.var.explained  <- summary(mel.e.pca.list[[i]])$importance[2,][2]
}

mel.e.region_info$is.alt.sig.predictor <- if_else(mel.e.region_info$pc1.alt.altitude.pval<=0.05, "yes", "no"); nrow(subset(mel.e.region_info, is.alt.sig.predictor=="yes"))/nrow(mel.e.region_info)*100
mel.e.region_info$is.neutral.wg.sig.predictor <- if_else(mel.e.region_info$pc1.alt.neutral.wg.PC1.r2<=0.05, "yes", "no"); nrow(subset(mel.e.region_info, is.neutral.wg.sig.predictor=="yes"))/nrow(mel.e.region_info)*100
mean(subset(mel.e.region_info, is.alt.sig.predictor=="yes")$pc1.alt.altitude.r2)
mean(subset(mel.e.region_info, is.neutral.wg.sig.predictor =="no")$pc1.alt.neutral.wg.PC1.r2)

###### save coefficients to mel.e.shdr.outlier ###### 
names(mel.e.region_info); names(mel.e.shdr.outlier)
# mel.e.shdr.outlier$pc1.alt.cor.coef.nointersect.sites.4stdv <- mel.e.region_info$pc1.alt.cor[match(mel.e.shdr.outlier$shdr.para.east.id, mel.e.region_info$shdr.para.east.id)]
# mel.e.shdr.outlier$pc1.alt.cor.pval.nointersect.sites.4stdv <- mel.e.region_info$pc1.alt.cor.pval[match(mel.e.shdr.outlier$shdr.para.east.id, mel.e.region_info$shdr.para.east.id)]

mel.e.shdr.outlier$pc1.alt.altitude.pval <- mel.e.region_info$pc1.alt.altitude.pval[match(mel.e.shdr.outlier$shdr.para.east.id, mel.e.region_info$shdr.para.east.id)]
mel.e.shdr.outlier$pc1.alt.neutral.wg.PC1.pval <- mel.e.region_info$pc1.alt.neutral.wg.PC1.pval[match(mel.e.shdr.outlier$shdr.para.east.id, mel.e.region_info$shdr.para.east.id)]
mel.e.shdr.outlier$pc1.alt.altitude.r2 <- mel.e.region_info$pc1.alt.altitude.r2[match(mel.e.shdr.outlier$shdr.para.east.id, mel.e.region_info$shdr.para.east.id)]; head(mel.e.shdr.outlier)
mel.e.shdr.outlier$pc1.alt.neutral.wg.PC1.r2 <- mel.e.region_info$pc1.alt.neutral.wg.PC1.r2[match(mel.e.shdr.outlier$shdr.para.east.id, mel.e.region_info$shdr.para.east.id)]; head(mel.e.shdr.outlier)

mel.e.shdr.outlier$PC1.var.explained.nointersect.sites.4stdv <- mel.e.region_info$PC1.var.explained[match(mel.e.shdr.outlier$shdr.para.east.id, mel.e.region_info$shdr.para.east.id)]
mel.e.shdr.outlier$PC2.var.explained.nointersect.sites.4stdv <- mel.e.region_info$PC2.var.explained[match(mel.e.shdr.outlier$shdr.para.east.id, mel.e.region_info$shdr.para.east.id)]; mel.e.shdr.outlier$PC2.var.explained.nointersect.sites.4stdv

#write.csv(mel.e.shdr.outlier, "data/shdr.summ/shdr.mel.east.outlier.df.csv", row.names = F); names(mel.e.shdr.outlier)


###### how similar are pcas across shdrs ###### 
# plot altitde r2 instead of pearson cor
N_PC=2
window_eigs <-sapply(mel.e.npy.list, function (cm) cov_pca(k=N_PC, covmat=cm))
window_dist<- pc_dist(t(window_eigs), npc=N_PC)
mds_axe<-cmdscale(window_dist)
rownames(mds_axe) <- sapply(mel.e.file.list, substr,45,57  ) 

## add info
mel.e.mds_axes.df <- as.data.frame(mds_axe)
mel.e.mds_axes.df$region <- sapply(mel.e.file.list, substr, 45,57 ) ; mel.e.mds_axes.df$region
mel.e.mds_axes.df$region.short <- sapply(mel.e.file.list, substr, 45,57  ) ; mel.e.mds_axes.df$region.short 
#mel.e.mds_axes.df$region.short[1:5] <- c("inv1.Hmelto0204", "inv2.2.Hmelto0211" ,"inv2.Hmelto0211" ,"non.inv2.Hmelto0215" ,"non.inv.Hmelto0215")
mel.e.mds_axes.df$region.size.bp <- mel.e.region_info$region.size[pmatch(mel.e.mds_axes.df$region, mel.e.region_info$shdr.para.east.id )]; mel.e.mds_axes.df$region.size.bp
mel.e.mds_axes.df$PC1.var.explained <- mel.e.region_info$PC1.var.explained[pmatch(mel.e.mds_axes.df$region, mel.e.region_info$shdr.para.east.id )]; mel.e.mds_axes.df$region.size.bp
#mel.e.mds_axes.df$pc1.alt.cor <- mel.e.region_info$pc1.alt.cor[pmatch(mel.e.mds_axes.df$region, mel.e.region_info$shdr.para.east.id )]; mel.e.mds_axes.df$region.size.bp
mel.e.mds_axes.df$pc1.alt.altitude.r2 <- mel.e.region_info$pc1.alt.altitude.r2[pmatch(mel.e.mds_axes.df$region, mel.e.region_info$shdr.para.east.id )]; mel.e.mds_axes.df$pc1.alt.altitude.r2
mel.e.mds_axes.df$is.non.inv <- if_else(substr(mel.e.mds_axes.df$region.short, 0,1)=="s", "no", "yes")
mel.e.mds_axes.df$is.allo.hdr <- mel.e.region_info$is.allo.hdr[pmatch(mel.e.mds_axes.df$region, mel.e.region_info$shdr.para.east.id )]; mel.e.mds_axes.df$region.size.bp

## plot
mel.e.all.regions.mds <- ggplot(data=subset(mel.e.mds_axes.df, is.non.inv!="yes"), aes(x=V1, y=V2 , fill=pc1.alt.altitude.r2))+
  geom_label_repel(aes(label=region.short ), max.overlaps = 30, fill="transparent", color="black", size=3 ) +
  scale_color_manual(values=c("#049E73", "#D65D00"))+
  geom_point(aes(size=region.size.bp/1000, colour=is.allo.hdr), alpha=.9,  shape=21, stroke=1.5)+
  ylab("MDS2") + xlab("MDS1")+ 
  scale_size_continuous(range=c(3,9), limits = c(0,950) ,name="Shared region\nsize KB (>4stdv)")+
  scale_fill_continuous(name= expression(paste(  "Altitutde R"^{2}*"")), limits = c(0,0.5), low="white", high="black" )+
  theme_bw()+ theme(legend.position = "none"); mel.e.all.regions.mds


###### save pc1 vs alt plots for significant shdrs only ######
names(mel.e.region_info)
mel.e.region_info.sig <- subset(mel.e.region_info, pc1.alt.altitude.pval<=0.05)

# loop only with significant regions and shdrs, for SI
mel.e.file.list.df <- as.data.frame(mel.e.file.list)
mel.e.region_info.sig$file.name <- paste("output/mds.shared/mel.e.stdv4/mel.e.region.", mel.e.region_info.sig$shdr.para.east.id, ".cov.npy", sep=""); mel.e.region_info.sig$file.name
mel.e.region_info.sig$is.shdr <- if_else(substr(mel.e.region_info.sig$shdr.para.east.id, 0,4)=="shdr", "yes", "no"); mel.e.region_info.sig$is.shdr 
mel.e.region_info.sig <- subset(mel.e.region_info.sig, is.shdr=="yes"); mel.e.region_info.sig 

mel.e.npy.list.sig <- list(); mel.e.pca.list.sig <- list()
mel.e.regions.list.sig <- list(); mel.e.pc1.vs.alt.plot.list.sig <- list()

for (i in 1:length(mel.e.region_info.sig$file.name)) {
  file <- mel.e.region_info.sig$file.name[[i]]
  mel.e.npy.list.sig[[i]] <-npyLoad(file=paste0(mel.e.region_info.sig$file.name[i]))
  names(mel.e.npy.list.sig)[[i]] <- file
  region <- names(mel.e.npy.list.sig)[[i]] %>% strsplit( "/" ) %>% sapply( tail, 1 )  %>% strsplit( "[.]" )  %>%   sapply( "[", 4:6 )  %>% paste(collapse=".")
  mel.e.regions.list.sig[[i]] <- region
  mel.e.pca.list.sig[[i]]<-prcomp(mel.e.npy.list.sig[[i]])
  rownames(mel.e.pca.list.sig[[i]]$x)<-mel.e.bam_names$V2
  
  # store pca info for plotting/analyses
  mel.e.pca.list.sig.df <- as.data.frame(mel.e.pca.list.sig[[i]]$x[,c(1:2)]); mel.e.pca.list.sig.df$id <- rownames(mel.e.pca.list.sig.df)
  mel.e.pca.list.sig.df$altitude <- subset(mel.e.list.info, id %in% mel.e.bam_names$V2)$altitude[match(mel.e.pca.list.sig.df$id, subset(mel.e.list.info, id %in% mel.e.bam_names$V2)$id)]
  mel.e.pca.list.sig.df$alt.type <- subset(mel.e.list.info, id %in% mel.e.bam_names$V2)$alt.type[match(mel.e.pca.list.sig.df$id, subset(mel.e.list.info, id %in% mel.e.bam_names$V2)$id)]
  mel.e.pca.list.sig.df$country <- subset(mel.e.list.info, id %in% mel.e.bam_names$V2)$country[match(mel.e.pca.list.sig.df$id, subset(mel.e.list.info, id %in% mel.e.bam_names$V2)$id)]
  mel.e.pca.list.sig.df$region <- region
  mel.e.pca.list.sig.df$alt.type.country <- paste(mel.e.pca.list.sig.df$alt.type, mel.e.pca.list.sig.df$country, sep=".")
  mel.e.pca.list.sig.df$is.allo <- mel.e.shdr.outlier$is.allo.hdr[match(mel.e.pca.list.sig.df$region, mel.e.shdr.outlier$shdr.para.east.id)]
  mel.e.pca.list.sig.df$shdr.type <- if_else(mel.e.shdr.outlier[mel.e.shdr.outlier$shdr.para.east.id==region,]$is.allo.hdr=="yes", "shdr.allo", "shdr.east")
  r2.alt.label <-  as.expression(bquote(Alt.~R^2==.( round(mel.e.region_info.sig[i,]$pc1.alt.altitude.r2, digits = 2)  )))
  
  # plot
  mel.e.pc1.vs.alt.plot.list.sig[[i]] <- ggplot(data=mel.e.pca.list.sig.df, aes(y=PC1, x=altitude, color=altitude))+
    scale_color_gradient( low = "#00FF00B3", high = "#0000FFF2")+
    geom_point(size=1 , aes(shape=alt.type.country ))+geom_smooth(method="lm", color="black")+
    scale_shape_manual(values=c(17,24,15,0,19,1))+
    # colour of border according to shdr type and personalised color scale
    ggnewscale::new_scale_color()+
    geom_rect(aes(color = shdr.type), fill="transparent", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, size=2) +
    scale_color_shdr.type()+
    annotate(geom = 'text', label = region , x =  0.2*(min(mel.e.pca.list.sig.df$altitude) + max(mel.e.pca.list.sig.df$altitude)), y = -Inf, hjust = 0, vjust = -0.25, size=4)+
    annotate(geom = "text", label = paste("PC1 ", round(summary(mel.e.pca.list.sig[[i]])$importance[2,][1], digits = 2)*100, "%", sep="") , 
             x =  0*(min(mel.e.pca.list.sig.df$altitude) + max(mel.e.pca.list.sig.df$altitude)), y = Inf, hjust = 0, vjust = 2.5, size=4)+
    stat_cor( hjust = -0.0, vjust = 0.25, size=4, aes(label =paste(..r.label.., cut(..p..,  breaks = c(-Inf, 0.0001, 0.001, 0.01, 0.05, Inf), labels = c("'****'", "'***'", "'**'", "'*'", "'ns'")), sep = "~"))) +
    annotate(geom = 'text', label =r2.alt.label , parse=TRUE, x =  0*(min(mel.e.pca.list.sig.df$altitude) + max(mel.e.pca.list.sig.df$altitude)), y = Inf, hjust = 0, vjust = 3, size=4)+
    scale_y_continuous(expand = c(.1,0))+
    theme_classic()+ theme(panel.border = element_rect(colour = "transparent", fill=NA, size=1),#axis.ticks.x=element_blank(),
                           axis.ticks.length.y =unit(-0.15, "cm"), #axis.ticks.margin=unit(0.5, "cm"),
                           axis.text.y = element_blank(),
                           axis.text.x = element_blank(), axis.title = element_blank(), axis.ticks=element_blank(),
                           #axis.text.y.right = element_text(size = 10, margin = unit(c(t = 0, r = 0, b = 0, l = -8), "mm"), colour="springgreen", face="bold"),
                           panel.background = element_blank(), #legend.text = element_text(size=8), legend.title =  element_text(size=12),
                           panel.grid.major = element_blank(),# axis.text.x=element_blank(), #axis.text.y=element_blank(), #
                           plot.margin = unit(c(0, 0, -.05,0), "cm"),
                           panel.grid.minor = element_blank(),# axis.title.x=element_blank(),#axis.title.y=element_blank(),
                           legend.text.align = 0, legend.position = "none",strip.text = element_blank(),
                           strip.background =element_rect(fill="transparent", color="transparent")); mel.e.pc1.vs.alt.plot.list.sig[[i]]
}

plot_grid(plotlist=mel.e.pc1.vs.alt.plot.list.sig, ncol=6, greedy = T)
ggsave2("plots/local.pca/mel.e.pc1.vs.alt.significant.no.intersect.png",  width = 210, height = 227, units = "mm")

## get legends for SI

ggplot(data=mel.e.pca.list.sig.df, aes(y=PC1, x=altitude, color=altitude))+
  scale_color_gradient( low = "#00FF00B3", high = "#0000FFF2")+
  geom_point(size=1 , aes(shape=alt.type.country ))+geom_smooth(method="lm", color="black")+
  geom_smooth(method="lm", color="black")+
  scale_shape_manual(values=c(17,24,15,0,19,1))+
  annotate(geom = 'text', label = region , x =  0.2*(min(mel.e.pca.list.sig.df$altitude) + max(mel.e.pca.list.sig.df$altitude)), y = -Inf, hjust = 0, vjust = -1, size=4)+
  annotate(geom = 'text', label = paste("PC1 ", round(summary(mel.e.pca.list.sig[[i]])$importance[2,][1], digits = 2)*100, "%", sep="") , 
           x =  0*(min(mel.e.pca.list.sig.df$altitude) + max(mel.e.pca.list.sig.df$altitude)), y = Inf, hjust = 0, vjust = 3, size=4)+
  stat_cor( hjust = -0.0, vjust = 0.8, size=4, aes(label =paste(..r.label.., cut(..p..,  breaks = c(-Inf, 0.0001, 0.001, 0.01, 0.05, Inf), labels = c("'****'", "'***'", "'**'", "'*'", "'ns'")), sep = "~")))+theme_classic()+
  theme(legend.text = element_text(size=18), legend.title =  element_text(size=18))

#################################### 5 melpomene west local PCA new ############## 
##### load indiv data #####
#add rownames
mel.w.list.info <- subset(all.bams.list.info,to.use.pbs=="yes" &species=="melpomene"&side.short=="w" )
mel.w.bam_names<-read.table("local/data/02.info/pop.list/mel.w.txt",header=F)
mel.w.bam_names$V1<-substr(mel.w.bam_names$V1, 76,84); mel.w.bam_names$V1
mel.w.bam_names$V2<-gsub('[.rm]', '', mel.w.bam_names$V1);mel.w.bam_names$V2
mel.w.bam_names$V2<- if_else(substr(mel.w.bam_names$V2,9, 9)=="_" , paste(mel.w.bam_names$V2, "ena", sep = ""), mel.w.bam_names$V2); mel.w.bam_names$V2; mel.w.list.info$id
mel.w.bam_names[mel.w.bam_names$V2=="CAM002856",]$V2 <- "CAM002856_ena"
mel.w.bam_names[mel.w.bam_names$V2=="CAM002857",]$V2 <- "CAM002857_ena"

# important !!!!!!!!!!!!!!!!!!!@!!!!!! - order 
mel.w.list.info <- mel.w.list.info[match(mel.w.bam_names$V2, mel.w.list.info$id),]
mel.w.bam_names$V2 ==mel.w.list.info$id
# info on regions
mel.w.region_info <-read.csv("local/data/02.info/sites/shared.exp50kb.mel.w.std4.new.mel.w.summ.csv")

##### loop #####

mel.w.file.list <- list.files("output/mds.shared/mel.w.stdv4/", full.names = T); mel.w.file.list

mel.w.npy.list <- list(); mel.w.pca.list <- list(); mel.w.pca.plot.list <- list()
mel.w.pc1.vs.alt.plot.list <- list(); mel.w.regions.list <- list()
mel.w.pca.list.df<- list()


for (i in 1:length(mel.w.file.list)) {
  file <- mel.w.file.list[[i]]
  mel.w.npy.list[[i]] <-npyLoad(file=paste0(mel.w.file.list[[i]]))
  names(mel.w.npy.list)[[i]] <- file
  region <- names(mel.w.npy.list)[[i]] %>% strsplit( "/" ) %>% sapply( tail, 1 )  %>% strsplit( "[.]" )  %>%   sapply( "[", 4:6 )  %>% paste(collapse=".")
  mel.w.regions.list[[i]] <- region
  mel.w.pca.list[[i]]<-prcomp(mel.w.npy.list[[i]])
  rownames(mel.w.pca.list[[i]]$x)<-mel.w.bam_names$V2
  
  # store pca info for plotting/analyses
  mel.w.pca.list.df[[i]] <- as.data.frame(mel.w.pca.list[[i]]$x[,c(1:2)]); mel.w.pca.list.df[[i]]$id <- rownames(mel.w.pca.list.df[[i]])
  mel.w.pca.list.df[[i]]$altitude <- subset(mel.w.list.info, id %in% mel.w.bam_names$V2)$altitude[match(mel.w.pca.list.df[[i]]$id, subset(mel.w.list.info, id %in% mel.w.bam_names$V2)$id)]
  mel.w.pca.list.df[[i]]$alt.type <- subset(mel.w.list.info, id %in% mel.w.bam_names$V2)$alt.type[match(mel.w.pca.list.df[[i]]$id, subset(mel.w.list.info, id %in% mel.w.bam_names$V2)$id)]
  mel.w.pca.list.df[[i]]$country <- subset(mel.w.list.info, id %in% mel.w.bam_names$V2)$country[match(mel.w.pca.list.df[[i]]$id, subset(mel.w.list.info, id %in% mel.w.bam_names$V2)$id)]
  mel.w.pca.list.df[[i]]$alt.type.country <- paste(mel.w.pca.list.df[[i]]$alt.type, mel.w.pca.list.df[[i]]$country, sep=".")
  mel.w.pca.list.df[[i]]$region <- region
  mel.w.pca.list.df[[i]]$shdr.type <- if_else(mel.w.shdr.outlier[mel.w.shdr.outlier$shdr.para.west.id==region,]$is.allo.hdr=="yes", "shdr.allo", "shdr.west")
  
  # save PCs 
  mel.w.list.info[,c(paste(region, ".PC1", sep = ""))] <-NA; mel.w.list.info[,c(paste(region, ".PC2", sep = ""))] <-NA
  mel.w.list.info[,c(paste(region, ".PC1", sep = ""))] <- mel.w.pca.list.df[[i]]$PC1; mel.w.list.info[,c(paste(region, ".PC2", sep = ""))] <- mel.w.pca.list.df[[i]]$PC2
  
  # plot pcas
  mel.w.pca.plot.list[[i]] <- ggplot(data=mel.w.pca.list.df[[i]], aes(x=PC1, y=PC2, colour=altitude))+ 
    geom_point(size=1, alpha=0.6 , aes(shape=alt.type.country))+
    scale_color_gradient( low = "#00FF00B3", high = "#0000FFF2")+
    # colour of border according to shdr type and personalised color scale
    ggnewscale::new_scale_color()+
    geom_rect(aes(color = shdr.type), fill="transparent", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, size=2) +
    scale_color_shdr.type()+
    scale_shape_manual(values=c(17,24,15,0,19))+
    #scale_color_manual(values = c("#0332FA", "#6ECE6F", "#0CFA00"))+
    annotate(geom = 'text', label = region , x =  -Inf, y = Inf, hjust = -0.1, vjust = 2, size=3)+
    annotate(geom = 'text', label = paste("PC1 ", round(summary(mel.w.pca.list[[i]])$importance[2,][1], digits = 2)*100, "%", sep="") , 
             x =  -Inf, y = Inf, hjust = -0.1, vjust = 3.5, size=3)+
    annotate(geom = 'text', label = paste("PC2 ", round(summary(mel.w.pca.list[[i]])$importance[2,][2], digits = 2)*100, "%", sep="") , 
             x =  -Inf, y = Inf, hjust = -0.1, vjust = 5, size=3)+
    theme_classic()+ theme(panel.border = element_rect(colour = "transparent", fill=NA, size=1),#axis.ticks.x=element_blank(),
                           axis.ticks.length.y =unit(-0.15, "cm"), #axis.ticks.margin=unit(0.5, "cm"),
                           axis.text.y = element_blank(),
                           axis.text.x = element_blank(), axis.title = element_blank(), axis.ticks=element_blank(),
                           #axis.text.y.right = element_text(size = 10, margin = unit(c(t = 0, r = 0, b = 0, l = -8), "mm"), colour="springgreen", face="bold"),
                           panel.background = element_blank(), #legend.text = element_text(size=8), legend.title =  element_text(size=12),
                           panel.grid.major = element_blank(),# axis.text.x=element_blank(), #axis.text.y=element_blank(), #
                           plot.margin = unit(c(0, 0, -.2,0), "cm"),
                           panel.grid.minor = element_blank(),# axis.title.x=element_blank(),#axis.title.y=element_blank(),
                           legend.text.align = 0, legend.position = "none",strip.text = element_blank(),
                           strip.background =element_rect(fill="transparent", color="transparent")); mel.w.pca.plot.list[[i]]
  # mel.w.pca.plot.list[[i]] <- ggExtra::ggMarginal(mel.w.pca.plot.list[[i]], type = "histogram", margins = "x", groupFill = T, bins=4); mel.w.pca.plot.list[[i]]
  
  
  mel.w.pc1.vs.alt.plot.list[[i]] <- ggplot(data=mel.w.pca.list.df[[i]], aes(y=PC1, x=altitude, color=altitude))+
    scale_color_gradient( low = "#00FF00B3", high = "#0000FFF2")+
    geom_point(size=1, aes(shape=alt.type.country))+
    scale_shape_manual(values=c(17,24,15,0,19))+
    geom_smooth(method="lm", color="black" )+
    # colour of border according to shdr type and personalised color scale
    ggnewscale::new_scale_color()+
    geom_rect(aes(color = shdr.type), fill="transparent", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, size=2) +
    scale_color_shdr.type()+
    annotate(geom = 'text', label = region , x =  0.15*(min(mel.w.pca.list.df[[i]]$altitude) + max(mel.w.pca.list.df[[i]]$altitude)), y = -Inf, hjust = 0, vjust = -1, size=3)+
    annotate(geom = 'text', label = paste("PC1 ", round(summary(mel.w.pca.list[[i]])$importance[2,][1], digits = 2)*100, "%", sep="") , 
             x =  0*(min(mel.w.pca.list.df[[i]]$altitude) + max(mel.w.pca.list.df[[i]]$altitude)), y = Inf, hjust = 0, vjust = 4.5, size=3)+
    stat_cor( hjust = -0.0, vjust = 0.6, size=3, aes(label =paste(..r.label.., cut(..p..,  breaks = c(-Inf, 0.0001, 0.001, 0.01, 0.05, Inf), labels = c("'****'", "'***'", "'**'", "'*'", "'ns'")), sep = "~"))) +
    theme_classic()+ theme(panel.border = element_rect(colour = "transparent", fill=NA, size=1),#axis.ticks.x=element_blank(),
                           axis.ticks.length.y =unit(-0.15, "cm"), #axis.ticks.margin=unit(0.5, "cm"),
                           axis.text.y = element_blank(),
                           axis.text.x = element_blank(), axis.title = element_blank(), axis.ticks=element_blank(),
                           #axis.text.y.right = element_text(size = 10, margin = unit(c(t = 0, r = 0, b = 0, l = -8), "mm"), colour="springgreen", face="bold"),
                           panel.background = element_blank(), #legend.text = element_text(size=8), legend.title =  element_text(size=12),
                           panel.grid.major = element_blank(),# axis.text.x=element_blank(), #axis.text.y=element_blank(), #
                           plot.margin = unit(c(0, 0, -.2,0), "cm"),
                           panel.grid.minor = element_blank(),# axis.title.x=element_blank(),#axis.title.y=element_blank(),
                           legend.text.align = 0, legend.position = "none",strip.text = element_blank(),
                           strip.background =element_rect(fill="transparent", color="transparent")); mel.w.pc1.vs.alt.plot.list[[i]]
  
  
} 

# plot_grid(plotlist=mel.w.pca.plot.list, ncol=6, greedy =T)
# ggsave2("plots/local.pca/mel.w.stdv4.no.intersect.all.png",  width = 190, height = 297, units = "mm")
# 
# plot_grid(plotlist=mel.w.pc1.vs.alt.plot.list, ncol=7, greedy = T)
# ggsave2("plots/local.pca/mel.w.pc1.vs.alt.stdv4.no.intersect.png",  width = 210, height = 327, units = "mm")


###### PC1 vs altitude coefficients ###### 
### add neutral pc1 coefficients
mel.w.list.info$neutral.wg.PC1 <- mel.w.neutral.wg.pca.df$PC1[match(mel.w.list.info$id, mel.w.neutral.wg.pca.df$id)]; mel.w.list.info$neutral.wg.PC1

mel.w.region_info <-read.csv("local/data/02.info/sites/shared.exp50kb.mel.w.std4.new.mel.w.summ.csv")

names( mel.w.list.info)
mel.w.region_info$pc1.alt.altitude.pval<- NA; mel.w.region_info$pc1.alt.neutral.wg.PC1.pval<- NA; mel.w.region_info$pc1.alt.total.r2<- NA
mel.w.region_info$pc1.alt.altitude.r2 <- NA; mel.w.region_info$pc1.alt.neutral.wg.PC1.r2 <- NA
mel.w.region_info$PC1.var.explained <- NA; mel.w.region_info$PC2.var.explained <- NA
mel.w.region_info$region.r.name <- unlist(mel.w.regions.list)[pmatch(mel.w.region_info$shdr.para.west.id ,unlist(mel.w.regions.list))]; mel.w.region_info$region.r.name
mel.w.region_info$is.allo.hdr <- mel.w.shdr.outlier$is.allo.hdr[match(mel.w.region_info$shdr.para.west.id, mel.w.shdr.outlier$shdr.para.west.id)]

# check vif not >4, any shdr
car::vif(lm(mel.w.list.info[,c(paste("shdr.west.001", ".PC1", sep = ""))] ~ altitude +neutral.wg.PC1, data=mel.w.list.info))

for (i in 1:length(mel.w.file.list)) {
  region <- names(mel.w.npy.list)[[i]] %>% strsplit( "/" ) %>% sapply( tail, 1 )  %>% strsplit( "[.]" )  %>%   sapply( "[", 4:6 )  %>% paste(collapse=".")
  mel.w.regions.list[[i]] <- region
  # get p.values for each predictor 
  altitude.pval <- round(summary(lm(mel.w.list.info[,c(paste(region, ".PC1", sep = ""))] ~ altitude +neutral.wg.PC1, data=mel.w.list.info))$coefficients[,4][2], digits = 3)
  neutral.wg.PC1.pval <- round(summary(lm(mel.w.list.info[,c(paste(region, ".PC1", sep = ""))] ~ altitude +neutral.wg.PC1, data=mel.w.list.info))$coefficients[,4][3], digits = 3)
  
  # get relative importance of each predictor
  total.model.r2 <- calc.relimp(lm(mel.w.list.info[,c(paste(region, ".PC1", sep = ""))] ~ altitude +neutral.wg.PC1, data=mel.w.list.info),rela=F)@R2 
  altitude.r2 <- calc.relimp(lm(mel.w.list.info[,c(paste(region, ".PC1", sep = ""))] ~ altitude +neutral.wg.PC1, data=mel.w.list.info),rela=F)@lmg["altitude"] 
  neutral.wg.PC1.r2 <- calc.relimp(lm(mel.w.list.info[,c(paste(region, ".PC1", sep = ""))] ~ altitude +neutral.wg.PC1, data=mel.w.list.info),rela=F)@lmg["neutral.wg.PC1"] 
  
  # cor.result <- cor.test(mel.w.list.info[,c(paste(region, ".PC1", sep = ""))] , mel.w.list.info$altitude , method = "pearson")
  # estimate <- cor.result$estimate; pval <- cor.result$p.value
  mel.w.region_info[subset(mel.w.region_info, region.r.name!="")$region.r.name == region,]$pc1.alt.altitude.pval<- altitude.pval
  mel.w.region_info[subset(mel.w.region_info, region.r.name!="")$region.r.name == region,]$pc1.alt.neutral.wg.PC1.pval<- neutral.wg.PC1.pval
  mel.w.region_info[subset(mel.w.region_info, region.r.name!="")$region.r.name == region,]$pc1.alt.total.r2<- total.model.r2
  mel.w.region_info[subset(mel.w.region_info, region.r.name!="")$region.r.name == region,]$pc1.alt.altitude.r2 <- altitude.r2
  mel.w.region_info[subset(mel.w.region_info, region.r.name!="")$region.r.name == region,]$pc1.alt.neutral.wg.PC1.r2 <- neutral.wg.PC1.r2
  mel.w.region_info[subset(mel.w.region_info, region.r.name!="")$region.r.name == region,]$PC1.var.explained <- summary(mel.w.pca.list[[i]])$importance[2,][1]
  mel.w.region_info[subset(mel.w.region_info, region.r.name!="")$region.r.name == region,]$PC2.var.explained  <- summary(mel.w.pca.list[[i]])$importance[2,][2]
}

mel.w.region_info$is.alt.sig.predictor <- if_else(mel.w.region_info$pc1.alt.altitude.pval<=0.05, "yes", "no"); nrow(subset(mel.w.region_info, is.alt.sig.predictor=="yes"))/nrow(mel.w.region_info)*100
mel.w.region_info$is.neutral.wg.sig.predictor <- if_else(mel.w.region_info$pc1.alt.neutral.wg.PC1.r2<=0.05, "yes", "no"); nrow(subset(mel.w.region_info, is.neutral.wg.sig.predictor=="yes"))/nrow(mel.w.region_info)*100
mean(subset(mel.w.region_info, is.alt.sig.predictor=="yes")$pc1.alt.altitude.r2)
mean(subset(mel.w.region_info, is.neutral.wg.sig.predictor =="no")$pc1.alt.neutral.wg.PC1.r2)

###### save coefficients to mel.w.shdr.outlier ###### 
names(mel.w.region_info); names(mel.w.shdr.outlier)
# mel.w.shdr.outlier$pc1.alt.cor.coef.nointersect.sites.4stdv <- mel.w.region_info$pc1.alt.cor[match(mel.w.shdr.outlier$shdr.para.west.id, mel.w.region_info$shdr.para.west.id)]
# mel.w.shdr.outlier$pc1.alt.cor.pval.nointersect.sites.4stdv <- mel.w.region_info$pc1.alt.cor.pval[match(mel.w.shdr.outlier$shdr.para.west.id, mel.w.region_info$shdr.para.west.id)]

mel.w.shdr.outlier$pc1.alt.altitude.pval <- mel.w.region_info$pc1.alt.altitude.pval[match(mel.w.shdr.outlier$shdr.para.west.id, mel.w.region_info$shdr.para.west.id)]
mel.w.shdr.outlier$pc1.alt.neutral.wg.PC1.pval <- mel.w.region_info$pc1.alt.neutral.wg.PC1.pval[match(mel.w.shdr.outlier$shdr.para.west.id, mel.w.region_info$shdr.para.west.id)]
mel.w.shdr.outlier$pc1.alt.altitude.r2 <- mel.w.region_info$pc1.alt.altitude.r2[match(mel.w.shdr.outlier$shdr.para.west.id, mel.w.region_info$shdr.para.west.id)]; head(mel.w.shdr.outlier)
mel.w.shdr.outlier$pc1.alt.neutral.wg.PC1.r2 <- mel.w.region_info$pc1.alt.neutral.wg.PC1.r2[match(mel.w.shdr.outlier$shdr.para.west.id, mel.w.region_info$shdr.para.west.id)]; head(mel.w.shdr.outlier)

mel.w.shdr.outlier$PC1.var.explained.nointersect.sites.4stdv <- mel.w.region_info$PC1.var.explained[match(mel.w.shdr.outlier$shdr.para.west.id, mel.w.region_info$shdr.para.west.id)]
mel.w.shdr.outlier$PC2.var.explained.nointersect.sites.4stdv <- mel.w.region_info$PC2.var.explained[match(mel.w.shdr.outlier$shdr.para.west.id, mel.w.region_info$shdr.para.west.id)]; mel.w.shdr.outlier$PC2.var.explained.nointersect.sites.4stdv

#write.csv(mel.w.shdr.outlier, "data/shdr.summ/shdr.mel.west.outlier.df.csv", row.names = F); names(mel.w.shdr.outlier)


###### how similar are pcas across shdrs ###### 
# plot altitde r2 instead of pearson cor
N_PC=2
window_eigs <-sapply(mel.w.npy.list, function (cm) cov_pca(k=N_PC, covmat=cm))
window_dist<- pc_dist(t(window_eigs), npc=N_PC)
mds_axe<-cmdscale(window_dist)
rownames(mds_axe) <- sapply(mel.w.file.list, substr,45,57  ) 

## add info
mel.w.mds_axes.df <- as.data.frame(mds_axe)
mel.w.mds_axes.df$region <- sapply(mel.w.file.list, substr, 45,57  ) ; mel.w.mds_axes.df$region
mel.w.mds_axes.df$region.short <- sapply(mel.w.file.list, substr, 45,57  ) ; mel.w.mds_axes.df$region.short 
#mel.w.mds_axes.df$region.short[1:5] <- c("inv1.Hmelto0204", "inv2.2.Hmelto0211" ,"inv2.Hmelto0211" ,"non.inv2.Hmelto0215" ,"non.inv.Hmelto0215")
mel.w.mds_axes.df$region.size.bp <- mel.w.region_info$region.size[pmatch(mel.w.mds_axes.df$region, mel.w.region_info$shdr.para.west.id )]; mel.w.mds_axes.df$region.size.bp
mel.w.mds_axes.df$PC1.var.explained <- mel.w.region_info$PC1.var.explained[pmatch(mel.w.mds_axes.df$region, mel.w.region_info$shdr.para.west.id )]; mel.w.mds_axes.df$region.size.bp
#mel.w.mds_axes.df$pc1.alt.cor <- mel.w.region_info$pc1.alt.cor[pmatch(mel.w.mds_axes.df$region, mel.w.region_info$shdr.para.west.id )]; mel.w.mds_axes.df$region.size.bp
mel.w.mds_axes.df$pc1.alt.altitude.r2 <- mel.w.region_info$pc1.alt.altitude.r2[pmatch(mel.w.mds_axes.df$region, mel.w.region_info$shdr.para.west.id )]; mel.w.mds_axes.df$pc1.alt.altitude.r2
mel.w.mds_axes.df$is.non.inv <- if_else(substr(mel.w.mds_axes.df$region.short, 0,1)=="s", "no", "yes")
mel.w.mds_axes.df$is.allo.hdr <- mel.w.region_info$is.allo.hdr[pmatch(mel.w.mds_axes.df$region, mel.w.region_info$shdr.para.west.id )]; mel.w.mds_axes.df$region.size.bp

## plot
mel.w.all.regions.mds <- ggplot(data=subset(mel.w.mds_axes.df, is.non.inv!="yes"), aes(x=V1, y=V2 , fill=pc1.alt.altitude.r2))+
  geom_label_repel(aes(label=region.short ), max.overlaps = 30, fill="transparent", color="black", size=3 ) +
  scale_color_manual(values=c( "#0372B2", "#D65D00"))+
  geom_point(aes(size=region.size.bp/1000, colour=is.allo.hdr), alpha=.9,  shape=21, stroke=1.5)+
  ylab("MDS2") + xlab("MDS1")+ 
  scale_size_continuous(range=c(3,9), limits = c(0,950) ,name="Shared region\nsize KB (>4stdv)")+
  scale_fill_continuous(name= expression(paste(  "Altitutde R"^{2}*"")), limits = c(0,0.5), low="white", high="black" )+
  theme_bw()+ theme(legend.position = "none"); mel.w.all.regions.mds


###### save pc1 vs alt plots for significant shdrs only ######
names(mel.w.region_info)
mel.w.region_info.sig <- subset(mel.w.region_info, pc1.alt.altitude.pval<=0.05)

# loop only with significant regions and shdrs, for SI
mel.w.file.list.df <- as.data.frame(mel.w.file.list)
mel.w.region_info.sig$file.name <- paste("output/mds.shared/mel.w.stdv4/mel.w.region.", mel.w.region_info.sig$shdr.para.west.id, ".cov.npy", sep=""); mel.w.region_info.sig$file.name
mel.w.region_info.sig$is.shdr <- if_else(substr(mel.w.region_info.sig$shdr.para.west.id, 0,4)=="shdr", "yes", "no"); mel.w.region_info.sig$is.shdr 
mel.w.region_info.sig <- subset(mel.w.region_info.sig, is.shdr=="yes"); mel.w.region_info.sig 

mel.w.npy.list.sig <- list(); mel.w.pca.list.sig <- list()
mel.w.regions.list.sig <- list(); mel.w.pc1.vs.alt.plot.list.sig <- list()

for (i in 1:length(mel.w.region_info.sig$file.name)) {
  file <- mel.w.region_info.sig$file.name[[i]]
  mel.w.npy.list.sig[[i]] <-npyLoad(file=paste0(mel.w.region_info.sig$file.name[i]))
  names(mel.w.npy.list.sig)[[i]] <- file
  region <- names(mel.w.npy.list.sig)[[i]] %>% strsplit( "/" ) %>% sapply( tail, 1 )  %>% strsplit( "[.]" )  %>%   sapply( "[", 4:6 )  %>% paste(collapse=".")
  mel.w.regions.list.sig[[i]] <- region
  mel.w.pca.list.sig[[i]]<-prcomp(mel.w.npy.list.sig[[i]])
  rownames(mel.w.pca.list.sig[[i]]$x)<-mel.w.bam_names$V2
  
  # store pca info for plotting/analyses
  mel.w.pca.list.sig.df <- as.data.frame(mel.w.pca.list.sig[[i]]$x[,c(1:2)]); mel.w.pca.list.sig.df$id <- rownames(mel.w.pca.list.sig.df)
  mel.w.pca.list.sig.df$altitude <- subset(mel.w.list.info, id %in% mel.w.bam_names$V2)$altitude[match(mel.w.pca.list.sig.df$id, subset(mel.w.list.info, id %in% mel.w.bam_names$V2)$id)]
  mel.w.pca.list.sig.df$alt.type <- subset(mel.w.list.info, id %in% mel.w.bam_names$V2)$alt.type[match(mel.w.pca.list.sig.df$id, subset(mel.w.list.info, id %in% mel.w.bam_names$V2)$id)]
  mel.w.pca.list.sig.df$country <- subset(mel.w.list.info, id %in% mel.w.bam_names$V2)$country[match(mel.w.pca.list.sig.df$id, subset(mel.w.list.info, id %in% mel.w.bam_names$V2)$id)]
  mel.w.pca.list.sig.df$region <- region
  mel.w.pca.list.sig.df$alt.type.country <- paste(mel.w.pca.list.sig.df$alt.type, mel.w.pca.list.sig.df$country, sep=".")
  mel.w.pca.list.sig.df$is.allo <- mel.w.shdr.outlier$is.allo.hdr[match(mel.w.pca.list.sig.df$region, mel.w.shdr.outlier$shdr.para.west.id)]
  mel.w.pca.list.sig.df$shdr.type <- if_else(mel.w.shdr.outlier[mel.w.shdr.outlier$shdr.para.west.id==region,]$is.allo.hdr=="yes", "shdr.allo", "shdr.west")
  r2.alt.label <-  as.expression(bquote(Alt.~R^2==.( round(mel.w.region_info.sig[i,]$pc1.alt.altitude.r2, digits = 2)  )))
  
  # plot
  mel.w.pc1.vs.alt.plot.list.sig[[i]] <- ggplot(data=mel.w.pca.list.sig.df, aes(y=PC1, x=altitude, color=altitude))+
    scale_color_gradient( low = "#00FF00B3", high = "#0000FFF2")+
    geom_point(size=1 , aes(shape=alt.type.country ))+geom_smooth(method="lm", color="black")+
    scale_shape_manual(values=c(17,24,15,0,19))+
    # colour of border according to shdr type and personalised color scale
    ggnewscale::new_scale_color()+
    geom_rect(aes(color = shdr.type), fill="transparent", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, size=2) +
    scale_color_shdr.type()+
    annotate(geom = 'text', label = region , x =  0.2*(min(mel.w.pca.list.sig.df$altitude) + max(mel.w.pca.list.sig.df$altitude)), y = -Inf, hjust = 0, vjust = -0.25, size=4)+
    annotate(geom = "text", label = paste("PC1 ", round(summary(mel.w.pca.list.sig[[i]])$importance[2,][1], digits = 2)*100, "%", sep="") , 
             x =  0*(min(mel.w.pca.list.sig.df$altitude) + max(mel.w.pca.list.sig.df$altitude)), y = Inf, hjust = 0, vjust = 2.5, size=4)+
    stat_cor( hjust = -0.0, vjust = 0.25, size=4, aes(label =paste(..r.label.., cut(..p..,  breaks = c(-Inf, 0.0001, 0.001, 0.01, 0.05, Inf), labels = c("'****'", "'***'", "'**'", "'*'", "'ns'")), sep = "~"))) +
    annotate(geom = 'text', label =r2.alt.label , parse=TRUE, x =  0*(min(mel.w.pca.list.sig.df$altitude) + max(mel.w.pca.list.sig.df$altitude)), y = Inf, hjust = 0, vjust = 3, size=4)+
    scale_y_continuous(expand = c(.1,0))+
    theme_classic()+ theme(panel.border = element_rect(colour = "transparent", fill=NA, size=1),#axis.ticks.x=element_blank(),
                           axis.ticks.length.y =unit(-0.15, "cm"), #axis.ticks.margin=unit(0.5, "cm"),
                           axis.text.y = element_blank(),
                           axis.text.x = element_blank(), axis.title = element_blank(), axis.ticks=element_blank(),
                           #axis.text.y.right = element_text(size = 10, margin = unit(c(t = 0, r = 0, b = 0, l = -8), "mm"), colour="springgreen", face="bold"),
                           panel.background = element_blank(), #legend.text = element_text(size=8), legend.title =  element_text(size=12),
                           panel.grid.major = element_blank(),# axis.text.x=element_blank(), #axis.text.y=element_blank(), #
                           plot.margin = unit(c(0, 0, -.05,0), "cm"),
                           panel.grid.minor = element_blank(),# axis.title.x=element_blank(),#axis.title.y=element_blank(),
                           legend.text.align = 0, legend.position = "none",strip.text = element_blank(),
                           strip.background =element_rect(fill="transparent", color="transparent")); mel.w.pc1.vs.alt.plot.list.sig[[i]]
}

plot_grid(plotlist=mel.w.pc1.vs.alt.plot.list.sig, ncol=7, greedy = T)
ggsave2("plots/local.pca/mel.w.pc1.vs.alt.significant.no.intersect.png",  width = 210, height = 320, units = "mm")


#################################### 6. plot mds of PCAs ####################################

era.e.all.regions.mds <- ggplot(data=subset(era.e.mds_axes.df, is.non.inv!="yes"), aes(x=V1, y=V2 , fill=pc1.alt.altitude.r2))+
  theme_bw()+ theme(legend.position = "none")+
  geom_point(aes(size=region.size.bp/1000, colour=is.allo.hdr), alpha=.9,  shape=21, stroke=1.5)+
  geom_label_repel(aes(label=substr(region.short, 11,13 )), max.overlaps = 30, fill="transparent", color="black", size=3 ) +
  scale_color_manual(values=c("#049E73", "#D65D00"))+
  ylab("MDS2") + xlab("MDS1")+ 
  scale_size_continuous(range=c(3,13), limits = c(0,1150) ,name="Shared region\nsize KB (>4stdv)")+
  scale_fill_continuous(name=expression(paste(  "Altitutde R"^{2}*"")), limits = c(0,0.75), low="white", high="black" ); era.e.all.regions.mds

era.w.all.regions.mds <- ggplot(data=subset(era.w.mds_axes.df), aes(x=V1, y=V2 , fill=pc1.alt.altitude.r2))+
  theme_bw()+ theme(legend.position = "none")+
  geom_label_repel(aes(label=substr(region.short, 11,13 )), max.overlaps = 30, fill="transparent", color="black", size=3 ) +
  geom_point(aes(size=region.size.bp/1000, colour=is.allo.hdr), alpha=.9,  shape=21, stroke=1.5)+
  scale_color_manual(values=c("#0372B2", "#D65D00"))+
  ylab("MDS2") + xlab("MDS1")+ 
  scale_size_continuous(range=c(3,13), limits = c(0,1150) ,name="Shared region\nsize KB (>4stdv)")+
  scale_fill_continuous(name=expression(paste(  "Altitutde R"^{2}*"")), limits = c(0,0.75), low="white", high="black" ); era.w.all.regions.mds


mel.e.all.regions.mds <- ggplot(data=subset(mel.e.mds_axes.df), aes(x=V1, y=V2 , fill=pc1.alt.altitude.r2))+
  theme_bw()+ theme(legend.position = "none")+
  geom_point(aes(size=region.size.bp/1000, colour=is.allo.hdr), alpha=.9,  shape=21, stroke=1.5)+
  geom_label_repel(aes(label=substr(region.short, 11,13 )), max.overlaps = 30, fill="transparent", color="black", size=3 ) +
  scale_color_manual(values=c("#049E73", "#D65D00"))+
  ylab("MDS2") + xlab("MDS1")+ 
  scale_size_continuous(range=c(3,13), limits = c(0,1150) ,name="Shared region\nsize KB (>4stdv)")+
  scale_fill_continuous(name=expression(paste(  "Altitutde R"^{2}*"")), limits = c(0,0.75), low="white", high="black" ); mel.e.all.regions.mds

mel.w.all.regions.mds <- ggplot(data=subset(mel.w.mds_axes.df), aes(x=V1, y=V2 , fill=pc1.alt.altitude.r2))+
  theme_bw()+ theme(legend.position = "none")+
  geom_point(aes(size=region.size.bp/1000, colour=is.allo.hdr), alpha=.9,  shape=21, stroke=1.5)+
  geom_label_repel(aes(label=substr(region.short, 11,13 )), max.overlaps = 30, fill="transparent", color="black", size=3 ) +
  scale_color_manual(values=c("#0372B2", "#D65D00"))+
  ylab("MDS2") + xlab("MDS1")+ 
  scale_size_continuous(range=c(3,13), limits = c(0,1150) ,name="Shared region\nsize KB (>4stdv)")+
  scale_fill_continuous(name=expression(paste(  "Altitutde R"^{2}*"")), limits = c(0,0.75), low="white", high="black" ); mel.w.all.regions.mds


(( era.w.all.regions.mds |era.e.all.regions.mds ) / ( mel.w.all.regions.mds| mel.e.all.regions.mds ))+
  plot_layout(guides = "collect")+ plot_annotation(tag_levels = 'A') & theme(legend.position = 'right')

ggsave2("plots/local.pca/era.mel.mds.of.local.pca.no.intersect.png",  width = 300, height = 260, units = "mm")

#################################### 7. analyses of shdrs stats x pca ####################################
names(era.e.shdr.outlier)
write.csv()

# write.csv(era.e.shdr.outlier, "data/shdr.summ/shdr.era.east.outlier.df.csv", row.names = F); names(era.e.shdr.outlier)
# write.csv(era.w.shdr.outlier,"data/shdr.summ/shdr.era.west.outlier.df.csv", row.names = F); head(era.w.shdr.outlier)
# write.csv(mel.e.shdr.outlier, "data/shdr.summ/shdr.mel.east.outlier.df.csv", row.names = F); head(mel.e.shdr.outlier)
# write.csv(mel.w.shdr.outlier, "data/shdr.summ/shdr.mel.west.outlier.df.csv", row.names = F); head(mel.w.shdr.outlier)

era.e.shdr.outlier <-subset(era.e.shdr.outlier, shdr.para.east.id   !="" & is.max.pbs.hig.above4=="yes")

era.e.shdr.outlier<-subset(era.e.shdr.outlier, is.max.pbs.hig.above4=="yes")
era.w.shdr.outlier<-subset(era.w.shdr.outlier, is.max.pbs.hig.above4=="yes")

nrow(subset(era.e.shdr.outlier, is.allo.hdr=="yes"))
nrow(subset(era.w.shdr.outlier, is.allo.hdr=="yes"))
nrow(subset(mel.e.shdr.outlier, is.allo.hdr=="yes")); nrow(subset(mel.e.shdr.outlier))
nrow(subset(mel.w.shdr.outlier, is.allo.hdr=="yes")); nrow(subset(mel.w.shdr.outlier))


subset(era.e.shdr.outlier, pc1.alt.cor.pval.nointersect.sites.4stdv<0.05 & is.allo.hdr=="yes")$shdr.para.east.id
subset(era.w.shdr.outlier, pc1.alt.cor.pval.nointersect.sites.4stdv<0.05 & is.allo.hdr=="yes")$shdr.para.west.id

subset(mel.e.shdr.outlier, pc1.alt.cor.pval.nointersect.sites.4stdv<0.05 & is.allo.hdr=="yes")$shdr.para.east.id
subset(mel.w.shdr.outlier, pc1.alt.cor.pval.nointersect.sites.4stdv<0.05 & is.allo.hdr=="yes")$shdr.para.west.id

shdr.list <-list()
shdr.list[[1]]<- era.e.shdr.outlier; shdr.list[[2]]<- era.w.shdr.outlier
shdr.list[[3]]<- mel.e.shdr.outlier; shdr.list[[4]]<- mel.w.shdr.outlier
shdr.list.summs <-list()
i=1
names(shdr.list[[i]])
for (i in 1:length(shdr.list)) {
  shdr.list[[i]]$is.signif <- if_else( shdr.list[[i]]$pc1.alt.altitude.pval<0.05, "yes", "no")
  shdr.list.summs[[i]] <-summarise(group_by(subset(shdr.list[[i]] ),  is.signif),
          n=n(),
          total.shdr = nrow(shdr.list[[i]]),
          perc.of.total = c((n/total.shdr)*100),
          mean.r2= mean(pc1.alt.altitude.r2),
          mean.pc1= mean(PC1.var.explained*100, na.rm=T),
          perc.of.sig.correlations.of.those.with.pca = ( n/ nrow(subset(shdr.list[[i]], pc1.alt.altitude.pval!="") ))*100  )
}; 
shdr.list.summs[[1]] ; shdr.list.summs[[2]] 

shdr.list.summs[[3]] ; shdr.list.summs[[4]]

100-27.5 # west shdrs used
100-39.8
100-84.6
100-74.1
  i=1
for (i in 1:length(shdr.list)) {
  shdr.list[[i]]$is.signif <- if_else( shdr.list[[i]]$pc1.alt.altitude.pval<0.05, "yes", "no")
  shdr.list.summs[[i]] <-summarise(group_by(subset(shdr.list[[i]] ),  is.signif, is.allo.hdr),
                                   n=n(),
                                   total.shdr=nrow(shdr.list[[i]]),
                                   total.allo.shdr=nrow(subset(shdr.list[[i]], is.allo.hdr=="yes" ) ),
                                   perc.of.total = (n/total.shdr)*100,
                                   perc.of.allo.hdr = (n/total.allo.shdr)*100,
                                   mean.alt.r2= mean(abs(pc1.alt.altitude.r2)),
                                   no.shdr.with.pca = nrow(subset(shdr.list[[i]], pc1.alt.altitude.pval!="") ),
                                   perc.of.sig.correlations.of.those.with.pca= ( n/ nrow(subset(shdr.list[[i]], pc1.alt.altitude.pval!="") ))*100 ,
                                   )
}; shdr.list.summs[[i]] 

shdr.list.summs[[1]] ; shdr.list.summs[[2]] 
shdr.list.summs[[3]] ; shdr.list.summs[[4]]

((10/27*100)+ (4/12*100))/2
((5/12*100)+ (1/10*100))/2

((10/27*100)+ (4/12*100))/2
((5/12*100)+ (1/10*100))/2

nrow(subset(shdr.list[[1]], is.allo.hdr=="yes"))
((38/133*100)+ (38/75*100))/2
((11/104*100)+ (11/58*100))/2


# mean % variation explained by PC1 across all pcas
(mean(shdr.list[[1]]$PC1.var.explained, na.rm = T) + mean(shdr.list[[2]]$PC1.var.explained, na.rm = T) + 
    mean(shdr.list[[3]]$PC1.var.explained, na.rm = T) + mean(shdr.list[[4]]$PC1.var.explained, na.rm = T) )/4

mean( c(shdr.list[[1]]$PC1.var.explained, shdr.list[[2]]$PC1.var.explained, shdr.list[[3]]$PC1.var.explained, shdr.list[[4]]$PC1.var.explained), na.rm = T)
sd( c(shdr.list[[1]]$PC1.var.explained, shdr.list[[2]]$PC1.var.explained, shdr.list[[3]]$PC1.var.explained, shdr.list[[4]]$PC1.var.explained), na.rm = T)



( )/2 ; shdr.list.summs[[2]] 

shdr.list.summs[[3]] ; shdr.list.summs[[4]]

########## partial r2 ##########

names(era.e.shdr.outlier)
ggplot(era.e.shdr.outlier, aes(x=1, y=pc1.alt.altitude.r2) )+
  geom_point()+
  geom_point(inherit.aes = F, data=era.e.shdr.outlier, aes(x=2, y=pc1.alt.neutral.wg.PC1.r2))

gather()


r2.era.e.shdr.outlier.summ <- gather(era.e.shdr.outlier[,c("shdr.para.east.id", "pc1.alt.altitude.r2", "pc1.alt.neutral.wg.PC1.r2")], type, r2, -shdr.para.east.id)

ggplot(r2.era.e.shdr.outlier.summ, aes(x=type, y=r2) )+
  geom_point()+geom_boxplot()+ geom_line(aes(group=shdr.para.east.id))+stat_compare_means()

r2.mel.e.shdr.outlier.summ <- gather(mel.e.shdr.outlier[,c("shdr.para.east.id", "pc1.alt.altitude.r2", "pc1.alt.neutral.wg.PC1.r2")], type, r2, -shdr.para.east.id)
ggplot(r2.mel.e.shdr.outlier.summ, aes(x=type, y=r2) )+
  geom_point()+geom_boxplot()+ geom_line(aes(group=shdr.para.east.id))+stat_compare_means()




