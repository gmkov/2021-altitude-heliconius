##### GMK April 2021 #####
#### outlier/snp/gene densities- use gwas results 10snp 10snp ####

rm(list=ls())
dev.off()
######################### packages ####


library(devtools)
library("windowscanr")
#install_github("drveera/ggman")
library(ggplot2)
library(dplyr)
library(tidyr)
library(qqman)
library(cowplot)
library(data.table)
library(stringr)
library(devtools)
library(ggman)
library(cowplot)
library(EnvStats)
library(zoo)
library(ggpubr)
options(scipen = 999)
#devtools::install_github("thomasp85/patchwork")
library(patchwork)
library(RIdeogram)
library(EnvStats)
library(ggnewscale)
library(ggrepel)
setwd("git/2021-altitude-heliconius/")
pop.info <- read.csv("02.info/pop.short.info.csv")

rand_non_overlapping_intervals <- function(total, interval_lengths){
  n_intervals <- length(interval_lengths) # number of intervals
  #subtract all interval lengths and add inserts in the remaining positions
  insert_sites <- sort(sample(0:(total-sum(interval_lengths)), n_intervals, replace = T)) 
  reordered_lengths <- sample(interval_lengths) # randomly order insert lengths
  lengths_cum_sum <- cumsum(reordered_lengths) #get the cumulative sum of insert lengths
  insert_sites <- insert_sites + lengths_cum_sum #adjust insert sites by cumulative insert lengths
  starts <- insert_sites - reordered_lengths + 1 #starts are insert sites minus lengths
  ends <- insert_sites
  intervals = Intervals(as.matrix(cbind(starts,ends)), closed = c(TRUE,TRUE),type="R")
  intervals
}


######################### 0. prep data ######
# shdr
# era.e.shdr <- read.csv("data/sharing/era.shdr.para.east.csv"); head(era.e.shdr)
# era.w.shdr <- read.csv("data/sharing/era.shdr.para.west.csv"); head(era.w.shdr)
# era.allo.shdr <- read.csv("data/sharing/era.shdr.allo.summ.csv"); head(era.allo.shdr)
# mel.e.shdr <- read.csv("data/sharing/mel.shdr.para.east.csv"); head(mel.e.shdr)
# mel.w.shdr <- read.csv("data/sharing/mel.shdr.para.west.csv"); head(mel.w.shdr)
# mel.allo.shdr <- read.csv("data/sharing/mel.shdr.allo.summ.csv"); head(mel.allo.shdr)

# selection stats outliers tajima delta pi dxy
era.e.shdr.outlier <- read.csv("data/shdr.summ/shdr.era.east.outlier.df.csv"); names(era.e.shdr.outlier)
era.w.shdr.outlier <- read.csv("data/shdr.summ/shdr.era.west.outlier.df.csv"); head(era.w.shdr.outlier)
mel.e.shdr.outlier <- read.csv("data/shdr.summ/shdr.mel.east.outlier.df.csv"); head(mel.e.shdr.outlier)
mel.w.shdr.outlier <- read.csv("data/shdr.summ/shdr.mel.west.outlier.df.csv"); head(mel.w.shdr.outlier)

era.e.shdr.outlier<- subset(era.e.shdr.outlier, start!="")

# altitude gwas
era.alt.gwas <- read.table("~/Dropbox/PhD/29.shape.paper/data/data.clean/era.n479/PC062_merged_Herato.PL.AD.HAPCUT2.inf0.05.corr.gwas.alt.admix.lrt0.gz.50snp.10snp.pos.ranks.short")
mel.alt.gwas <- read.table("~/Dropbox/PhD/29.shape.paper/data/data.clean/mel.187/run164_merged_Hmel.PL.AD.HAPCUT2.inf0.05.corr.gwas.alt.admix.doasso6.50snp.10snp.pos.ranks.short")

mel.alt.gwas$CHR <- as.numeric(as.character(substr(mel.alt.gwas$scaff, 6, 7)))
  
# fst
era.fst <- read.table("~/Dropbox/PhD/29.shape.paper/data/data.clean/era.n479/era.fst.5kb.1kb.idx", row.names = NULL, header = TRUE ); str(era.fst)
mel.fst <- read.table("~/Dropbox/PhD/29.shape.paper/data/data.clean/mel.187/mel.fst.5kb.1kb.idx", row.names = NULL, header = TRUE ); str(era.fst)

setwd("/Users/gabrielamontejokovacevich/Dropbox (Cambridge University)/PhD/29.shape.paper/data/data.clean/selectionStatsHaplotaggingData/")

# genome sizes
ref.scaff.era <- read.table("~/Dropbox (Cambridge University)/PhD/7_Assocation_studies/9_ANGSD/0.ref/Heliconius_erato_demophoon_v1_-_scaffolds.fa.fai", row.names = NULL)
ref.scaff.mel <- read.table("~/Dropbox (Cambridge University)/PhD/7_Assocation_studies/9_ANGSD/0.ref/Hmel2.5.scaffolds.fa.fai", row.names = NULL)
scafEnds.era <- cumsum(ref.scaff.era[,2]); offset.era <- scafEnds.era - ref.scaff.era[,2]
scafEnds.mel <- cumsum(ref.scaff.mel[,2]); offset.mel <- scafEnds.mel - ref.scaff.mel[,2]
era.genome_size <- tail(scafEnds.era,1)
mel.genome_size <- tail(scafEnds.mel,1)

############################################### 1. plotting sel. scans ERA  ##################################################################
################# 0. prep data ##################################
# Read in erato scaffold info:
ref.scaff.era <- read.table('Heliconius_erato_demophoon_v1_-_scaffolds.fa.fai')
ref.scaff.era<-ref.scaff.era[,1:2]
names(ref.scaff.era)<-c("scaff","length")
ref.scaff.era$add<-c(0,cumsum(ref.scaff.era$length)[-length(ref.scaff.era$scaff)])
ref.scaff.era$CHR<-as.integer(substr(ref.scaff.era$scaff,7,8))
ref.scaff.era<-ref.scaff.era[ref.scaff.era$scaff!="Herato_mt",]

scafEnds.era <- cumsum(ref.scaff.era[,2]); offset.era <- scafEnds.era - ref.scaff.era[,2]
era.genome_size <- tail(scafEnds.era,1)

# Generate chromosome info:
chr.era<-aggregate(ref.scaff.era[,3],list(ref.scaff.era$CHR),min)
names(chr.era)<-c("CHR","add")
chr.era$length<-aggregate(ref.scaff.era[,2],list(ref.scaff.era$CHR),sum)[,2]
chr.era$mid<-chr.era$add+chr.era$length/2

# Read in Fst values of erato
erato<-read.table("erato.fst.10kb.idx.withHMM",header=T,sep="\t")
erato$win_mid_add<-erato$midPos+ref.scaff.era[match(as.character(erato$chr),as.character(ref.scaff.era$scaff)),"add"]
erato<-erato[!is.na(erato$win_mid_add),]
erato$CHR = as.factor(as.integer(str_sub(erato$chr,7,8)))
fst<-erato[!is.na(erato$CHR),]

# Read in SWEED and CLR
sweed.era<-read.table("PC062_merged_all.sweeD.notabilis_lativitta.out",header=T)
sweed.era$win_mid_add<-sweed.era$POS_GRID+ref.scaff.era[match(as.character(sweed.era$CHROM),as.character(ref.scaff.era$scaff)),"add"]
sweed.era$CHR=as.factor(as.integer(str_sub(sweed.era$CHROM,7,8)))

# Read in pi values
piLativ<-read.table("PC062_merged.PL.AD.HAPCUT2.inf0.05.lativitta.10kb.windowed.pi",header=T)
piNota<-read.table("PC062_merged.PL.AD.HAPCUT2.inf0.05.notabilis.10kb.windowed.pi",header=T)
names(piLativ)<-c("Chr","start","end","n","pi")
names(piNota)<-c("Chr","start","end","n","pi")
piLativ$WinCenter<-piLativ$start+5000
piNota$WinCenter<-piNota$start+5000
piLativ$win_mid_add<-piLativ$WinCenter+ref.scaff.era[match(as.character(piLativ$Chr),as.character(ref.scaff.era$scaff)),"add"]
piNota$win_mid_add<-piNota$WinCenter+ref.scaff.era[match(as.character(piNota$Chr),as.character(ref.scaff.era$scaff)),"add"]
piLativ$CHR=as.factor(as.integer(str_sub(piLativ$Chr,7,8)))
piNota$CHR=as.factor(as.integer(str_sub(piNota$Chr,7,8)))

piLativ$piNota <- piNota$pi[match(piLativ$win_mid_add, piNota$win_mid_add)]
names(piLativ)[5] <- "piLativ"; names(piLativ)
piLativ$piNota.minus.piLativ <- piLativ$piNota -piLativ$piLativ
piLativ.Nota <- piLativ

# iHS,  CLR
iHSera<-read.table("PC062_merged_all.erato.notabilis_lativitta.iHS.out",header=T)
iHSera$win_mid_add<-iHSera$START+ref.scaff.era[match(as.character(iHSera$CHROM),as.character(ref.scaff.era$scaff)),"add"]
iHSera$CHR=as.factor(as.integer(str_sub(iHSera$CHROM,7,8)))

omega.era<-read.table("PC062_merged_all.omega.notabilis_lativitta.wHeader.out",header=T)
omega.era$win_mid_add<-omega.era$POS_NOTABILIS+ref.scaff.era[match(as.character(omega.era$CHROM),as.character(ref.scaff.era$scaff)),"add"]
omega.era$CHR=as.factor(as.integer(str_sub(omega.era$CHROM,7,8)))

### prep chr for plotting ##
# prep chr midpoints for plotting
axisdf.era = era.fst %>% dplyr::group_by(CHR) %>% dplyr::summarize(center=( max(BP.wg) + min(BP.wg) ) / 2 ); axisdf.era
axisdf.era$CHR[21]<-"Z";axisdf.era$CHR[21]
# create start ends of chr for plotting
head(ref.scaff.era)
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

################# 1. plot selection stats with ROI ##################################
########## 1.1 deltapi, fst, ROI ##########
era.east.west.hdr.piLativ.Nota.plot <- list()
for ( i in 1:21) {
  # so that the loop carries on despite errors in some of the chr (to do with there not being any shdrs)
  # but problems with colouring those chromosomes where only one level present, allo/not, like chr8
  tryCatch({
    era.east.west.hdr.piLativ.Nota.plot[[i]] <- ggplot(subset(piLativ.Nota, CHR==i), aes(x=win_mid_add, y=piNota.minus.piLativ)) +
      geom_rect(inherit.aes = F, data=subset(era.e.shdr, chr==i ), aes(xmin=start, xmax=end,ymin=-0.02,ymax=0.009, fill=is.allopatric), colour="transparent", alpha=.5 )+
      scale_fill_manual(values=c("#049E73", "#D65D00"))+
      ggnewscale::new_scale_fill()+
      geom_rect(inherit.aes = F, data=subset(era.w.shdr, chr==i ), aes(xmin=start, xmax=end, ymin=-0.02,ymax=0.009, fill=is.allopatric), colour="transparent", alpha=.5 )+
      scale_fill_manual(values=c("#0372B2", "#D65D00"))+
      #geom_rect(inherit.aes = F, data=subset(era.w.shdr, chr==i &is.allopatric=="no"), aes(xmin=start, xmax=end, ymin=0,ymax=1), colour="transparent", fill="#0372B2", alpha=.3 )+
      annotate("text",  x=Inf, y = Inf, label = i, vjust=1.5, hjust=1.1)  +
      geom_label_repel(inherit.aes = F, data=subset(era.e.shdr, chr==i ), aes(x=(start+end)/2, y=0.009, label=substr(shdr.para.east.id,11,13 ), color=is.allopatric), min.segment.length = 0.05) +
      scale_colour_manual(values=c("#049E73", "#D65D00"))+
      #geom_label_repel(inherit.aes = F, data=subset(era.w.shdr, chr==i &is.allopatric=="no"), aes(x=(start+end)/2, y=0.89, label=substr(shdr.para.west.id,11,13 )), color="#0372B2", min.segment.length = 0.05) +
      geom_point(aes(y=piNota.minus.piLativ), colour="black",size=.5,alpha=.8) +
      geom_line(aes(y=rollmean(piNota.minus.piLativ, 10, na.pad=TRUE)), colour="black",size=.8,alpha=.8) +
      scale_x_continuous( expand = c(0, 0)) +
      scale_y_continuous( expand = c(0, 0), limits=c(-0.02,0.01)  ) +    # remove space between plot area and x axis
      ylab("|piLativ.Nota|")+
      theme_classic()+
      theme(legend.position="none")+ # Remove legend
      theme(axis.line.x = element_line(color="black", size = 0.5),axis.line.y = element_line(color="black", size = 0.5),axis.title.y = element_text(size=14),
            plot.margin = unit(c(0.5,0.2,0.2,0.2), "lines"),axis.text = element_text(size=10), #axis.text.x = element_blank(),
            axis.text.y = element_text(size=10),axis.title.x = element_blank(),       panel.background = element_rect(fill = "transparent"), plot.background = element_rect(fill = "transparent", color = NA),
            panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.background = element_rect(fill = "transparent"), legend.box.background = element_rect(fill = "transparent")) 
    
    ggsave2(plot = era.east.west.hdr.piLativ.Nota.plot[[i]],width = 7, height = 3, filename = paste0("git/2021-altitude-heliconius/plots/haplo.stats/deltapi/era/chr", i, ".png"))
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}

########## 1.1 omega (black, orange), fst, ROI ##########
era.east.west.hdr.omega.era.plot <- list()
for ( i in 1:21) {
  # so that the loop carries on despite errors in some of the chr (to do with there not being any shdrs)
  # but problems with colouring those chromosomes where only one level present, allo/not, like chr8
  tryCatch({
    era.east.west.hdr.omega.era.plot[[i]] <- ggplot(subset(omega.era, CHR==i), aes(x=win_mid_add, y=Omega_MAX_NOTABILIS)) +
      geom_rect(inherit.aes = F, data=subset(era.e.shdr, chr==i ), aes(xmin=start, xmax=end, ymin=0,ymax=1750, fill=is.allopatric), colour="transparent", alpha=.5 )+
      scale_fill_manual(values=c("#049E73", "#D65D00"))+
      ggnewscale::new_scale_fill()+
      geom_rect(inherit.aes = F, data=subset(era.w.shdr, chr==i ), aes(xmin=start, xmax=end, ymin=0,ymax=1750, fill=is.allopatric), colour="transparent", alpha=.5 )+
      scale_fill_manual(values=c("#0372B2", "#D65D00"))+
      #geom_rect(inherit.aes = F, data=subset(era.w.shdr, chr==i &is.allopatric=="no"), aes(xmin=start, xmax=end, ymin=0,ymax=1), colour="transparent", fill="#0372B2", alpha=.3 )+
      annotate("text",  x=Inf, y = Inf, label = i, vjust=1.5, hjust=1.1)  +
      geom_label_repel(inherit.aes = F, data=subset(era.e.shdr, chr==i ), aes(x=(start+end)/2, y=1750, label=substr(shdr.para.east.id,11,13 ), color=is.allopatric), min.segment.length = 0.05) +
      scale_colour_manual(values=c("#049E73", "#D65D00"))+
      #geom_label_repel(inherit.aes = F, data=subset(era.w.shdr, chr==i &is.allopatric=="no"), aes(x=(start+end)/2, y=0.89, label=substr(shdr.para.west.id,11,13 )), color="#0372B2", min.segment.length = 0.05) +
      geom_point(aes(y=Omega_MAX_NOTABILIS), colour="black",size=.5,alpha=.8) +
      geom_point(aes(y=Omega_MAX_LATIVITTA), colour="orange",size=.5,alpha=.9) +
      geom_line(aes(y=rollmean(Omega_MAX_NOTABILIS, 10, na.pad=TRUE)), colour="black",size=.8,alpha=.8) +
      geom_line(aes(y=rollmean(Omega_MAX_LATIVITTA, 10, na.pad=TRUE)), colour="orange",size=.8,alpha=.8) +
      scale_x_continuous( expand = c(0, 0)) +
      scale_y_continuous( expand = c(0, 0), limits=c(0,2000)  ) +    # remove space between plot area and x axis
      ylab("|omega.era|")+
      theme_classic()+
      theme(legend.position="none")+ # Remove legend
      theme(axis.line.x = element_line(color="black", size = 0.5),axis.line.y = element_line(color="black", size = 0.5),axis.title.y = element_text(size=14),
            plot.margin = unit(c(0.5,0.2,0.2,0.2), "lines"),axis.text = element_text(size=10), #axis.text.x = element_blank(),
            axis.text.y = element_text(size=10),axis.title.x = element_blank(),       panel.background = element_rect(fill = "transparent"), plot.background = element_rect(fill = "transparent", color = NA),
            panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.background = element_rect(fill = "transparent"), legend.box.background = element_rect(fill = "transparent")) 
    
    ggsave2(plot = era.east.west.hdr.omega.era.plot[[i]],width = 7, height = 3, filename = paste0("git/2021-altitude-heliconius/plots/haplo.stats/omega/era/chr", i, ".png"))
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}



########## 1.1 iHS (black, orange), fst, ROI ##########
era.east.west.hdr.iHS.plot <- list()
for ( i in 1:21) {
  # so that the loop carries on despite errors in some of the chr (to do with there not being any shdrs)
  # but problems with colouring those chromosomes where only one level present, allo/not, like chr8
  tryCatch({
  era.east.west.hdr.iHS.plot[[i]] <-  ggplot(subset(iHSera, CHR==i), aes(x=win_mid_add, y=iHS_NOTABILIS)) +
    geom_rect(inherit.aes = F, data=subset(era.e.shdr, chr==i ), aes(xmin=start, xmax=end, ymin=0,ymax=1, fill=is.allopatric), colour="transparent", alpha=.5 )+
    scale_fill_manual(values=c("#049E73", "#D65D00"))+
    ggnewscale::new_scale_fill()+
    geom_rect(inherit.aes = F, data=subset(era.w.shdr, chr==i ), aes(xmin=start, xmax=end, ymin=0,ymax=1, fill=is.allopatric), colour="transparent", alpha=.5 )+
    scale_fill_manual(values=c("#0372B2", "#D65D00"))+
    #geom_rect(inherit.aes = F, data=subset(era.w.shdr, chr==i &is.allopatric=="no"), aes(xmin=start, xmax=end, ymin=0,ymax=1), colour="transparent", fill="#0372B2", alpha=.3 )+
    annotate("text",  x=Inf, y = Inf, label = i, vjust=1.5, hjust=1.1)  +
    geom_label_repel(inherit.aes = F, data=subset(era.e.shdr, chr==i ), aes(x=(start+end)/2, y=0.99, label=substr(shdr.para.east.id,11,13 ), color=is.allopatric), min.segment.length = 0.05) +
    scale_colour_manual(values=c("#049E73", "#D65D00"))+
    #geom_label_repel(inherit.aes = F, data=subset(era.w.shdr, chr==i &is.allopatric=="no"), aes(x=(start+end)/2, y=0.89, label=substr(shdr.para.west.id,11,13 )), color="#0372B2", min.segment.length = 0.05) +
    geom_point(aes(y=iHS_NOTABILIS), colour="black",size=.5,alpha=.8) +
    geom_point(aes(y=iHS_LATIVITTA), colour="orange",size=.5,alpha=.9) +
    geom_line(aes(y=rollmean(iHS_NOTABILIS, 10, na.pad=TRUE)), colour="black",size=.8,alpha=.8) +
    geom_line(aes(y=rollmean(iHS_LATIVITTA, 10, na.pad=TRUE)), colour="orange",size=.8,alpha=.8) +
    scale_x_continuous( expand = c(0, 0)) +
    scale_y_continuous( expand = c(0, 0), limits=c(0,1)  ) +    # remove space between plot area and x axis
    ylab("|iHS|")+
    theme_classic()+
    theme(legend.position="none")+ # Remove legend
    theme(axis.line.x = element_line(color="black", size = 0.5),axis.line.y = element_line(color="black", size = 0.5),axis.title.y = element_text(size=14),
          plot.margin = unit(c(0.5,0.2,0.2,0.2), "lines"),axis.text = element_text(size=10), #axis.text.x = element_blank(),
          axis.text.y = element_text(size=10),axis.title.x = element_blank(),       panel.background = element_rect(fill = "transparent"), plot.background = element_rect(fill = "transparent", color = NA),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.background = element_rect(fill = "transparent"), legend.box.background = element_rect(fill = "transparent")) 
  ggsave2(plot = era.east.west.hdr.iHS.plot[[i]],width = 7, height = 3, filename = paste0("git/2021-altitude-heliconius/plots/haplo.stats/iHS/era/chr", i, ".png"))
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}

########## 1.4 sweed (black, orange), fst, ROI ##########
head(sweed.era)
era.east.west.hdr.CLR.plot <- list()
for ( i in 1:21) {
  # so that the loop carries on despite errors in some of the chr (to do with there not being any shdrs)
  # but problems with colouring those chromosomes where only one level present, allo/not, like chr8
  tryCatch({
    era.east.west.hdr.CLR.plot[[i]] <-  ggplot(subset(sweed.era, CHR==i), aes(x=win_mid_add, y=CLR_NOTABILIS)) +
      geom_rect(inherit.aes = F, data=subset(era.e.shdr, chr==i ), aes(xmin=start, xmax=end, ymin=0,ymax=200, fill=is.allopatric), colour="transparent", alpha=.5 )+
      scale_fill_manual(values=c("#049E73", "#D65D00"))+
      ggnewscale::new_scale_fill()+
      geom_rect(inherit.aes = F, data=subset(era.w.shdr, chr==i ), aes(xmin=start, xmax=end, ymin=0,ymax=200, fill=is.allopatric), colour="transparent", alpha=.5 )+
      scale_fill_manual(values=c("#0372B2", "#D65D00"))+
      #geom_rect(inherit.aes = F, data=subset(era.w.shdr, chr==i &is.allopatric=="no"), aes(xmin=start, xmax=end, ymin=0,ymax=1), colour="transparent", fill="#0372B2", alpha=.3 )+
      annotate("text",  x=Inf, y = Inf, label = i, vjust=1.5, hjust=1.1)  +
      geom_label_repel(inherit.aes = F, data=subset(era.e.shdr, chr==i ), aes(x=(start+end)/2, y=200, label=substr(shdr.para.east.id,11,13 ), color=is.allopatric), min.segment.length = 0.05) +
      scale_colour_manual(values=c("#049E73", "#D65D00"))+
      #geom_label_repel(inherit.aes = F, data=subset(era.w.shdr, chr==i &is.allopatric=="no"), aes(x=(start+end)/2, y=0.89, label=substr(shdr.para.west.id,11,13 )), color="#0372B2", min.segment.length = 0.05) +
      geom_point(aes(y=CLR_NOTABILIS), colour="black",size=.75,alpha=.8) +
      geom_point(aes(y=CLR_LATIVITTA), colour="orange",size=.75,alpha=.9) +
      #geom_line(aes(y=rollmean(CLR_NOTABILIS, 10, na.pad=TRUE)), colour="black",size=.8,alpha=.8) +
      # geom_line(aes(y=rollmean(CLR_LATIVITTA, 10, na.pad=TRUE)), colour="orange",size=.8,alpha=.8) +
      scale_x_continuous( expand = c(0, 0)) +
      scale_y_continuous( expand = c(0, 0), limits=c(0,200)  ) +    # remove space between plot area and x axis
      ylab("|CLR|")+
      theme_classic()+
      theme(legend.position="none")+ # Remove legend
      theme(axis.line.x = element_line(color="black", size = 0.5),axis.line.y = element_line(color="black", size = 0.5),axis.title.y = element_text(size=14),
            plot.margin = unit(c(0.5,0.2,0.2,0.2), "lines"),axis.text = element_text(size=10), #axis.text.x = element_blank(),
            axis.text.y = element_text(size=10),axis.title.x = element_blank(),       panel.background = element_rect(fill = "transparent"), plot.background = element_rect(fill = "transparent", color = NA),
            panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.background = element_rect(fill = "transparent"), legend.box.background = element_rect(fill = "transparent")) 
    ggsave2(plot = era.east.west.hdr.CLR.plot[[i]],width = 7, height = 3, filename = paste0("git/2021-altitude-heliconius/plots/haplo.stats/sweed/era/chr", i, ".png"))
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}


############################################### 2. plotting sel.scans  MEL ##################################################################
################# 0. prep data ##################################
# Read in scaffold info:
ref.scaff.mel <- read.table('Hmel2.5.fa.fai')
ref.scaff.mel<-ref.scaff.mel[,1:2]
names(ref.scaff.mel)<-c("scaff","length")
ref.scaff.mel$add<-c(0,cumsum(ref.scaff.mel$length)[-length(ref.scaff.mel$scaff)])
ref.scaff.mel<-ref.scaff.mel[ref.scaff.mel$scaff!="Hmel_complete_mtDNA",]
ref.scaff.mel$CHR<-as.integer(substr(ref.scaff.mel$scaff,6,7))
ref.scaff.mel<-ref.scaff.mel[!grepl(ref.scaff.mel$scaff,pattern = "Hmel200"),]
scafEnds.mel <- cumsum(ref.scaff.mel[,2]); offset.mel <- scafEnds.mel - ref.scaff.mel[,2]
mel.genome_size <- tail(scafEnds.mel,1)

# Genera chr.melomosome info:
chr.mel<-aggregate(ref.scaff.mel[,3],list(ref.scaff.mel$CHR),min)
names(chr.mel)<-c("CHR","add")
chr.mel$length<-aggregate(ref.scaff.mel[,2],list(ref.scaff.mel$CHR),sum)[,2]
chr.mel$mid<-chr.mel$add+chr.mel$length/2

# Read in STITCH pi
piMalle<-read.table("run164_merged.PL.AD.HAPCUT2.inf0.05.malleti.10kb.windowed.pi",header=T)
piPless<-read.table("run164_merged.PL.AD.HAPCUT2.inf0.05.plesseni.10kb.windowed.pi",header=T)
names(piMalle)<-c("Chr","start","end","n","pi")
names(piPless)<-c("Chr","start","end","n","pi")
piMalle$WinCenter<-piMalle$start+5000
piPless$WinCenter<-piPless$start+5000
piPless$win_mid_add<-piPless$WinCenter+ref.scaff.mel[match(as.character(piPless$Chr),as.character(ref.scaff.mel$scaff)),"add"]
piMalle$win_mid_add<-piMalle$WinCenter+ref.scaff.mel[match(as.character(piMalle$Chr),as.character(ref.scaff.mel$scaff)),"add"]
piMalle$CHR=as.factor(as.integer(str_sub(piMalle$Chr,6,7)))
piPless$CHR=as.factor(as.integer(str_sub(piPless$Chr,6,7)))

piMalle$piPless <- piPless$pi[match(piMalle$win_mid_add, piPless$win_mid_add)]
names(piMalle)[5] <- "piMalle"; names(piMalle)
piMalle$piPless.minus.piMalle <- piMalle$piPless -piMalle$piMalle
piMalle.Pless <- piMalle
head(piMalle.Pless)


# Read in SWEED and CLR
sweed.mel<-read.table("run164_merged_all.sweeD.plesseni_malleti.out",header=T)
sweed.mel$win_mid_add<-sweed.mel$POS_GRID+ref.scaff.mel[match(as.character(sweed.mel$CHROM),as.character(ref.scaff.mel$scaff)),"add"]
sweed.mel$CHR=as.factor(as.integer(str_sub(sweed.mel$CHROM,6,7)))

iHSmel<-read.table("run164_merged_all.melpomene.plesseni_malleti.iHS.out",header=T)
iHSmel$win_mid_add<-iHSmel$START+ref.scaff.mel[match(as.character(iHSmel$CHROM),as.character(ref.scaff.mel$scaff)),"add"]
iHSmel$CHR=as.factor(as.integer(str_sub(iHSmel$CHROM,6,7)))

omega.mel<-read.table("run164_merged_all.omega.plesseni_malleti.wHeader.out",header=T)
omega.mel$win_mid_add<-omega.mel$POS_PLESSENI+ref.scaff.mel[match(as.character(omega.mel$CHROM),as.character(ref.scaff.mel$scaff)),"add"]
omega.mel$CHR=as.factor(as.integer(str_sub(omega.mel$CHROM,6,7)))

################# 1. plot selection stats with ROI ##################################
head(piMalle); head(sweed.mel); head(omega.mel)
# create midpos per chr for plotting ggplot
axisdf.mel = piMalle %>% dplyr::group_by(CHR) %>% dplyr::summarize(center=( max(win_mid_add ) + min(win_mid_add ) ) / 2 ); axisdf.mel
str(axisdf.mel)
axisdf.mel$CHR <- as.character(axisdf.mel$CHR)
axisdf.mel$CHR[21]<-"Z"; axisdf.mel$CHR[21]

# create start ends of chr for plotting
head(ref.scaff)
mel.ref.chr.pos <- ref.scaff.mel
mel.ref.chr.pos$start <- mel.ref.chr.pos$add
mel.ref.chr.pos$end <-  mel.ref.chr.pos$start +mel.ref.chr.pos$length; head(mel.ref.chr.pos)
mel.ref.chr.pos <- summarise(group_by(mel.ref.chr.pos, CHR),
                             start=min(start),
                             end=max(end)); mel.ref.chr.pos

########## 1.1 deltapi, fst, ROI ##########
mel.east.west.hdr.era.alt.gwas.plot <- list()
for ( i in 1:21) {
  # so that the loop carries on despite errors in some of the chr (to do with there not being any shdrs)
  # but problems with colouring those chromosomes where only one level present, allo/not, like chr8
  tryCatch({
    mel.east.west.hdr.era.alt.gwas.plot[[i]] <- ggplot(subset(era.alt.gwas, CHR==i), aes(x=win_mid_add, y=piPless.minus.piMalle)) +
      geom_rect(inherit.aes = F, data=subset(mel.e.shdr, chr==i ), aes(xmin=start, xmax=end,ymin=-0.02,ymax=0.009, fill=is.allopatric), colour="transparent", alpha=.5 )+
      scale_fill_manual(values=c("#049E73", "#D65D00"))+
      ggnewscale::new_scale_fill()+
      geom_rect(inherit.aes = F, data=subset(mel.w.shdr, chr==i ), aes(xmin=start, xmax=end, ymin=-0.02,ymax=0.009, fill=is.allopatric), colour="transparent", alpha=.5 )+
      scale_fill_manual(values=c("#0372B2", "#D65D00"))+
      #geom_rect(inherit.aes = F, data=subset(mel.w.shdr, chr==i &is.allopatric=="no"), aes(xmin=start, xmax=end, ymin=0,ymax=1), colour="transparent", fill="#0372B2", alpha=.3 )+
      annotate("text",  x=Inf, y = Inf, label = i, vjust=1.5, hjust=1.1)  +
      geom_label_repel(inherit.aes = F, data=subset(mel.e.shdr, chr==i ), aes(x=(start+end)/2, y=0.009, label=substr(shdr.para.east.id,11,13 ), color=is.allopatric), min.segment.length = 0.05) +
      scale_colour_manual(values=c("#049E73", "#D65D00"))+
      #geom_label_repel(inherit.aes = F, data=subset(mel.w.shdr, chr==i &is.allopatric=="no"), aes(x=(start+end)/2, y=0.89, label=substr(shdr.para.west.id,11,13 )), color="#0372B2", min.segment.length = 0.05) +
      geom_point(aes(y=piPless.minus.piMalle), colour="black",size=.5,alpha=.8) +
      geom_line(aes(y=rollmean(piPless.minus.piMalle, 10, na.pad=TRUE)), colour="black",size=.8,alpha=.8) +
      scale_x_continuous( expand = c(0, 0)) +
      scale_y_continuous( expand = c(0, 0), limits=c(-0.02,0.01)  ) +    # remove space between plot area and x axis
      ylab("|era.alt.gwas|")+
      theme_classic()+
      theme(legend.position="none")+ # Remove legend
      theme(axis.line.x = element_line(color="black", size = 0.5),axis.line.y = element_line(color="black", size = 0.5),axis.title.y = element_text(size=14),
            plot.margin = unit(c(0.5,0.2,0.2,0.2), "lines"),axis.text = element_text(size=10), #axis.text.x = element_blank(),
            axis.text.y = element_text(size=10),axis.title.x = element_blank(),       panel.background = element_rect(fill = "transparent"), plot.background = element_rect(fill = "transparent", color = NA),
            panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.background = element_rect(fill = "transparent"), legend.box.background = element_rect(fill = "transparent")) 
    
    ggsave2(plot = mel.east.west.hdr.era.alt.gwas.plot[[i]],width = 7, height = 3, filename = paste0("git/2021-altitude-heliconius/plots/haplo.stats/deltapi/mel/chr", i, ".png"))
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}

########## 1.1 omega (black, orange), fst, ROI ##########

mel.east.west.hdr.omega.mel.plot <- list()
for ( i in 1:21) {
  # so that the loop carries on despite errors in some of the chr (to do with there not being any shdrs)
  # but problems with colouring those chromosomes where only one level present, allo/not, like chr8
  tryCatch({
    mel.east.west.hdr.omega.mel.plot[[i]] <- ggplot(subset(omega.mel, CHR==i), aes(x=win_mid_add, y=Omega_MAX_PLESSENI)) +
      geom_rect(inherit.aes = F, data=subset(mel.e.shdr, chr==i ), aes(xmin=start, xmax=end, ymin=0,ymax=2750, fill=is.allopatric), colour="transparent", alpha=.5 )+
      scale_fill_manual(values=c("#049E73", "#D65D00"))+
      ggnewscale::new_scale_fill()+
      geom_rect(inherit.aes = F, data=subset(mel.w.shdr, chr==i ), aes(xmin=start, xmax=end, ymin=0,ymax=2750, fill=is.allopatric), colour="transparent", alpha=.5 )+
      scale_fill_manual(values=c("#0372B2", "#D65D00"))+
      #geom_rect(inherit.aes = F, data=subset(mel.w.shdr, chr==i &is.allopatric=="no"), aes(xmin=start, xmax=end, ymin=0,ymax=1), colour="transparent", fill="#0372B2", alpha=.3 )+
      annotate("text",  x=Inf, y = Inf, label = i, vjust=1.5, hjust=1.1)  +
      geom_label_repel(inherit.aes = F, data=subset(mel.e.shdr, chr==i ), aes(x=(start+end)/2, y=2750, label=substr(shdr.para.east.id,11,13 ), color=is.allopatric), min.segment.length = 0.05) +
      scale_colour_manual(values=c("#049E73", "#D65D00"))+
      #geom_label_repel(inherit.aes = F, data=subset(mel.w.shdr, chr==i &is.allopatric=="no"), aes(x=(start+end)/2, y=0.89, label=substr(shdr.para.west.id,11,13 )), color="#0372B2", min.segment.length = 0.05) +
      geom_point(aes(y=Omega_MAX_PLESSENI), colour="black",size=.5,alpha=.8) +
      geom_point(aes(y=Omega_MAX_MALLETI), colour="orange",size=.5,alpha=.9) +
      geom_line(aes(y=rollmean(Omega_MAX_PLESSENI, 10, na.pad=TRUE)), colour="black",size=.8,alpha=.8) +
      geom_line(aes(y=rollmean(Omega_MAX_MALLETI, 10, na.pad=TRUE)), colour="orange",size=.8,alpha=.8) +
      scale_x_continuous( expand = c(0, 0)) +
      scale_y_continuous( expand = c(0, 0), limits=c(0,3000)  ) +    # remove space between plot area and x axis
      ylab("|omega.mel|")+
      theme_classic()+
      theme(legend.position="none")+ # Remove legend
      theme(axis.line.x = element_line(color="black", size = 0.5),axis.line.y = element_line(color="black", size = 0.5),axis.title.y = element_text(size=14),
            plot.margin = unit(c(0.5,0.2,0.2,0.2), "lines"),axis.text = element_text(size=10), #axis.text.x = element_blank(),
            axis.text.y = element_text(size=10),axis.title.x = element_blank(),       panel.background = element_rect(fill = "transparent"), plot.background = element_rect(fill = "transparent", color = NA),
            panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.background = element_rect(fill = "transparent"), legend.box.background = element_rect(fill = "transparent")) 
    
    ggsave2(plot = mel.east.west.hdr.omega.mel.plot[[i]],width = 7, height = 3, filename = paste0("git/2021-altitude-heliconius/plots/haplo.stats/omega/mel/chr", i, ".png"))
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}



########## 1.1 iHS (black, orange), fst, ROI ##########
mel.east.west.hdr.iHS.plot <- list()
for ( i in 1:21) {
  # so that the loop carries on despite errors in some of the chr (to do with there not being any shdrs)
  # but problems with colouring those chromosomes where only one level present, allo/not, like chr8
  tryCatch({
    mel.east.west.hdr.iHS.plot[[i]] <-  ggplot(subset(iHSmel, CHR==i), aes(x=win_mid_add, y=iHS_PLESSENI)) +
      geom_rect(inherit.aes = F, data=subset(mel.e.shdr, chr==i ), aes(xmin=start, xmax=end, ymin=0,ymax=1, fill=is.allopatric), colour="transparent", alpha=.5 )+
      scale_fill_manual(values=c("#049E73", "#D65D00"))+
      ggnewscale::new_scale_fill()+
      geom_rect(inherit.aes = F, data=subset(mel.w.shdr, chr==i ), aes(xmin=start, xmax=end, ymin=0,ymax=1, fill=is.allopatric), colour="transparent", alpha=.5 )+
      scale_fill_manual(values=c("#0372B2", "#D65D00"))+
      #geom_rect(inherit.aes = F, data=subset(mel.w.shdr, chr==i &is.allopatric=="no"), aes(xmin=start, xmax=end, ymin=0,ymax=1), colour="transparent", fill="#0372B2", alpha=.3 )+
      annotate("text",  x=Inf, y = Inf, label = i, vjust=1.5, hjust=1.1)  +
      geom_label_repel(inherit.aes = F, data=subset(mel.e.shdr, chr==i ), aes(x=(start+end)/2, y=0.99, label=substr(shdr.para.east.id,11,13 ), color=is.allopatric), min.segment.length = 0.05) +
      scale_colour_manual(values=c("#049E73", "#D65D00"))+
      #geom_label_repel(inherit.aes = F, data=subset(mel.w.shdr, chr==i &is.allopatric=="no"), aes(x=(start+end)/2, y=0.89, label=substr(shdr.para.west.id,11,13 )), color="#0372B2", min.segment.length = 0.05) +
      geom_point(aes(y=iHS_PLESSENI), colour="black",size=.5,alpha=.8) +
      geom_point(aes(y=iHS_MALLETI), colour="orange",size=.5,alpha=.9) +
      geom_line(aes(y=rollmean(iHS_PLESSENI, 10, na.pad=TRUE)), colour="black",size=.8,alpha=.8) +
      geom_line(aes(y=rollmean(iHS_MALLETI, 10, na.pad=TRUE)), colour="orange",size=.8,alpha=.8) +
      scale_x_continuous( expand = c(0, 0)) +
      scale_y_continuous( expand = c(0, 0), limits=c(0,1)  ) +    # remove space between plot area and x axis
      ylab("|iHS|")+
      theme_classic()+
      theme(legend.position="none")+ # Remove legend
      theme(axis.line.x = element_line(color="black", size = 0.5),axis.line.y = element_line(color="black", size = 0.5),axis.title.y = element_text(size=14),
            plot.margin = unit(c(0.5,0.2,0.2,0.2), "lines"),axis.text = element_text(size=10), #axis.text.x = element_blank(),
            axis.text.y = element_text(size=10),axis.title.x = element_blank(),       panel.background = element_rect(fill = "transparent"), plot.background = element_rect(fill = "transparent", color = NA),
            panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.background = element_rect(fill = "transparent"), legend.box.background = element_rect(fill = "transparent")) 
    ggsave2(plot = mel.east.west.hdr.iHS.plot[[i]],width = 7, height = 3, filename = paste0("git/2021-altitude-heliconius/plots/haplo.stats/iHS/mel/chr", i, ".png"))
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}

########## 1.4 sweed (black, orange), fst, ROI ##########
head(sweed.mel)
mel.east.west.hdr.CLR.plot <- list()
for ( i in 1:21) {
  # so that the loop carries on despite errors in some of the chr (to do with there not being any shdrs)
  # but problems with colouring those chromosomes where only one level present, allo/not, like chr8
  tryCatch({
    mel.east.west.hdr.CLR.plot[[i]] <-  ggplot(subset(sweed.mel, CHR==i), aes(x=win_mid_add, y=CLR_PLESSENI)) +
      geom_rect(inherit.aes = F, data=subset(mel.e.shdr, chr==i ), aes(xmin=start, xmax=end, ymin=0,ymax=200, fill=is.allopatric), colour="transparent", alpha=.5 )+
      scale_fill_manual(values=c("#049E73", "#D65D00"))+
      ggnewscale::new_scale_fill()+
      geom_rect(inherit.aes = F, data=subset(mel.w.shdr, chr==i ), aes(xmin=start, xmax=end, ymin=0,ymax=200, fill=is.allopatric), colour="transparent", alpha=.5 )+
      scale_fill_manual(values=c("#0372B2", "#D65D00"))+
      #geom_rect(inherit.aes = F, data=subset(mel.w.shdr, chr==i &is.allopatric=="no"), aes(xmin=start, xmax=end, ymin=0,ymax=1), colour="transparent", fill="#0372B2", alpha=.3 )+
      annotate("text",  x=Inf, y = Inf, label = i, vjust=1.5, hjust=1.1)  +
      geom_label_repel(inherit.aes = F, data=subset(mel.e.shdr, chr==i ), aes(x=(start+end)/2, y=200, label=substr(shdr.para.east.id,11,13 ), color=is.allopatric), min.segment.length = 0.05) +
      scale_colour_manual(values=c("#049E73", "#D65D00"))+
      #geom_label_repel(inherit.aes = F, data=subset(mel.w.shdr, chr==i &is.allopatric=="no"), aes(x=(start+end)/2, y=0.89, label=substr(shdr.para.west.id,11,13 )), color="#0372B2", min.segment.length = 0.05) +
      geom_point(aes(y=CLR_PLESSENI), colour="black",size=.75,alpha=.8) +
      geom_point(aes(y=CLR_MALLETI), colour="orange",size=.75,alpha=.9) +
      #geom_line(aes(y=rollmean(CLR_PLESSENI, 10, na.pad=TRUE)), colour="black",size=.8,alpha=.8) +
      # geom_line(aes(y=rollmean(CLR_MALLETI, 10, na.pad=TRUE)), colour="orange",size=.8,alpha=.8) +
      scale_x_continuous( expand = c(0, 0)) +
      scale_y_continuous( expand = c(0, 0), limits=c(0,200)  ) +    # remove space between plot area and x axis
      ylab("|CLR|")+
      theme_classic()+
      theme(legend.position="none")+ # Remove legend
      theme(axis.line.x = element_line(color="black", size = 0.5),axis.line.y = element_line(color="black", size = 0.5),axis.title.y = element_text(size=14),
            plot.margin = unit(c(0.5,0.2,0.2,0.2), "lines"),axis.text = element_text(size=10), #axis.text.x = element_blank(),
            axis.text.y = element_text(size=10),axis.title.x = element_blank(),       panel.background = element_rect(fill = "transparent"), plot.background = element_rect(fill = "transparent", color = NA),
            panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.background = element_rect(fill = "transparent"), legend.box.background = element_rect(fill = "transparent")) 
    ggsave2(plot = mel.east.west.hdr.CLR.plot[[i]],width = 7, height = 3, filename = paste0("git/2021-altitude-heliconius/plots/haplo.stats/sweed/mel/chr", i, ".png"))
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}



############################################### 3. plotting ALTITUDE GWAS #################################
names(era.alt.gwas)
era.east.west.hdr.era.alt.gwas.plot <- list()
dev.off()
for ( i in 1:21) {
  # so that the loop carries on despite errors in some of the chr (to do with there not being any shdrs)
  # but problems with colouring those chromosomes where only one level present, allo/not, like chr8
  tryCatch({
     era.east.west.hdr.era.alt.gwas.plot[[i]] <- ggplot(subset(era.alt.gwas, CHR==i), aes(x=BP.wg, y=-log10(LRT.pval_median))) +
      geom_rect(inherit.aes = F, data=subset(era.e.shdr, chr==i ), aes(xmin=start, xmax=end,ymin=-0.02,ymax=25, fill=is.allopatric), colour="transparent", alpha=.5 )+
      scale_fill_manual(values=c("#049E73", "#D65D00"))+
      ggnewscale::new_scale_fill()+
      geom_rect(inherit.aes = F, data=subset(era.w.shdr, chr==i ), aes(xmin=start, xmax=end, ymin=-0.02,ymax=25, fill=is.allopatric), colour="transparent", alpha=.5 )+
      scale_fill_manual(values=c("#0372B2", "#D65D00"))+
      #geom_rect(inherit.aes = F, data=subset(era.w.shdr, chr==i &is.allopatric=="no"), aes(xmin=start, xmax=end, ymin=0,ymax=1), colour="transparent", fill="#0372B2", alpha=.3 )+
      annotate("text",  x=Inf, y = Inf, label = i, vjust=1.5, hjust=1.1)  +
      geom_label_repel(inherit.aes = F, data=subset(era.e.shdr, chr==i ), aes(x=(start+end)/2, y=23, label=substr(shdr.para.east.id,11,13 ), color=is.allopatric), min.segment.length = 0.05) +
      scale_colour_manual(values=c("#049E73", "#D65D00"))+
      #geom_label_repel(inherit.aes = F, data=subset(era.w.shdr, chr==i &is.allopatric=="no"), aes(x=(start+end)/2, y=0.89, label=substr(shdr.para.west.id,11,13 )), color="#0372B2", min.segment.length = 0.05) +
      geom_point(aes(y=-log10(LRT.pval_median)), colour="black",size=.5,alpha=.8) +
      #geom_line(aes(y=rollmean(-log10(LRT.pval_median), 10, na.pad=TRUE)), colour="black",size=.8,alpha=.8) +
      scale_x_continuous( expand = c(0, 0)) +
      scale_y_continuous( expand = c(0, 0), limits=c(-0.02,25)  ) +    # remove space between plot area and x axis
      ylab("era.alt.gwas")+
      theme_classic()+
      theme(legend.position="none")+ # Remove legend
      theme(axis.line.x = element_line(color="black", size = 0.5),axis.line.y = element_line(color="black", size = 0.5),axis.title.y = element_text(size=14),
            plot.margin = unit(c(0.5,0.2,0.2,0.2), "lines"),axis.text = element_text(size=10), #axis.text.x = element_blank(),
            axis.text.y = element_text(size=10),axis.title.x = element_blank(),       panel.background = element_rect(fill = "transparent"), plot.background = element_rect(fill = "transparent", color = NA),
            panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.background = element_rect(fill = "transparent"), legend.box.background = element_rect(fill = "transparent")) 
    
     cowplot::ggsave2(plot = era.east.west.hdr.era.alt.gwas.plot[[i]], width = 19, height = 8, units = "cm", dpi = 100, filename = paste0("git/2021-altitude-heliconius/plots/haplo.stats/alt.gwas/era/chr", i, ".png"))
    dev.off()
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}

mel.east.west.hdr.mel.alt.gwas.plot <- list()
dev.off()
for ( i in 1:21) {
  # so that the loop carries on despite errors in some of the chr (to do with there not being any shdrs)
  # but problems with colouring those chromosomes where only one level present, allo/not, like chr8
  tryCatch({
    mel.east.west.hdr.mel.alt.gwas.plot[[i]] <- ggplot(subset(mel.alt.gwas, CHR==i), aes(x=BP.wg, y=-log10(LRT.pval_median))) +
      geom_rect(inherit.aes = F, data=subset(mel.e.shdr, chr==i ), aes(xmin=start, xmax=end,ymin=-0.02,ymax=25, fill=is.allopatric), colour="transparent", alpha=.5 )+
      scale_fill_manual(values=c("#049E73", "#D65D00"))+
      ggnewscale::new_scale_fill()+
      geom_rect(inherit.aes = F, data=subset(mel.w.shdr, chr==i ), aes(xmin=start, xmax=end, ymin=-0.02,ymax=25, fill=is.allopatric), colour="transparent", alpha=.5 )+
      scale_fill_manual(values=c("#0372B2", "#D65D00"))+
      #geom_rect(inherit.aes = F, data=subset(mel.w.shdr, chr==i &is.allopatric=="no"), aes(xmin=start, xmax=end, ymin=0,ymax=1), colour="transparent", fill="#0372B2", alpha=.3 )+
      annotate("text",  x=Inf, y = Inf, label = i, vjust=1.5, hjust=1.1)  +
      geom_label_repel(inherit.aes = F, data=subset(mel.e.shdr, chr==i ), aes(x=(start+end)/2, y=23, label=substr(shdr.para.east.id,11,13 ), color=is.allopatric), min.segment.length = 0.05) +
      scale_colour_manual(values=c("#049E73", "#D65D00"))+
      #geom_label_repel(inherit.aes = F, data=subset(mel.w.shdr, chr==i &is.allopatric=="no"), aes(x=(start+end)/2, y=0.89, label=substr(shdr.para.west.id,11,13 )), color="#0372B2", min.segment.length = 0.05) +
      geom_point(aes(y=-log10(LRT.pval_median)), colour="black",size=.5,alpha=.8) +
      #geom_line(aes(y=rollmean(-log10(LRT.pval_median), 10, na.pad=TRUE)), colour="black",size=.8,alpha=.8) +
      scale_x_continuous( expand = c(0, 0)) +
      scale_y_continuous( expand = c(0, 0), limits=c(-0.02,25)  ) +    # remove space between plot area and x axis
      ylab("mel.alt.gwas")+
      theme_classic()+
      theme(legend.position="none")+ # Remove legend
      theme(axis.line.x = element_line(color="black", size = 0.5),axis.line.y = element_line(color="black", size = 0.5),axis.title.y = element_text(size=14),
            plot.margin = unit(c(0.5,0.2,0.2,0.2), "lines"),axis.text = element_text(size=10), #axis.text.x = element_blank(),
            axis.text.y = element_text(size=10),axis.title.x = element_blank(),       panel.background = element_rect(fill = "transparent"), plot.background = element_rect(fill = "transparent", color = NA),
            panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.background = element_rect(fill = "transparent"), legend.box.background = element_rect(fill = "transparent")) 
    
    cowplot::ggsave2(plot = mel.east.west.hdr.mel.alt.gwas.plot[[i]], width = 19, height = 8, units = "cm", dpi = 100, filename = paste0("git/2021-altitude-heliconius/plots/haplo.stats/alt.gwas/mel/chr", i, ".png"))
    dev.off()
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}


############################################### 4. era SIMS ALL STATS against east shdr ###############################################
################# 0. prep data ####################
########## add hdr.id to datasets; east only ##########
head(piLativ.Nota); head(omega.era); head(iHSera); head(sweed.era); head(era.alt.gwas)
match.by.range <- function(vec) {
  x <- era.e.shdr
  if(length(.x <- which(vec >= x$start & vec <= x$end ))) .x else NA }

piLativ.Nota$shdr.para.east.id  <-era.e.shdr$shdr.para.east.id[mapply(match.by.range, piLativ.Nota$win_mid_add)]
omega.era$shdr.para.east.id  <-era.e.shdr$shdr.para.east.id[mapply(match.by.range, omega.era$win_mid_add)]
iHSera$shdr.para.east.id  <-era.e.shdr$shdr.para.east.id[mapply(match.by.range, iHSera$win_mid_add)]
sweed.era$shdr.para.east.id  <-era.e.shdr$shdr.para.east.id[mapply(match.by.range, sweed.era$win_mid_add)]
era.alt.gwas$shdr.para.east.id  <-era.e.shdr$shdr.para.east.id[mapply(match.by.range, era.alt.gwas$BP.wg)]

########## summarise per hdr and concatenate ##########
head(piLativ.Nota)
piLativ.Nota.summ <- summarise(group_by(subset(piLativ.Nota), shdr.para.east.id),
                               piNota.minus.piLativ.min=min(piNota.minus.piLativ)); piLativ.Nota.summ <- subset(piLativ.Nota.summ , shdr.para.east.id!=""); piLativ.Nota.summ

omega.era.summ <- summarise(group_by(subset(omega.era), shdr.para.east.id),
                            omega.era.max=max(Omega_MAX_NOTABILIS)); omega.era.summ <- subset(omega.era.summ , shdr.para.east.id!=""); omega.era.summ

iHSera.summ <- summarise(group_by(subset(iHSera), shdr.para.east.id),
                         iHS_NOTABILIS.max=max(iHS_NOTABILIS)); iHSera.summ <- subset(iHSera.summ , shdr.para.east.id!=""); iHSera.summ

sweed.era.summ <- summarise(group_by(subset(sweed.era), shdr.para.east.id),
                            CLR_NOTABILIS.max=max(CLR_NOTABILIS)); sweed.era.summ <- subset(sweed.era.summ , shdr.para.east.id!=""); sweed.era.summ

era.alt.gwas.summ <- summarise(group_by(subset(era.alt.gwas), shdr.para.east.id),
                               log10.LRT.pval_median.max=max(-log10(LRT.pval_median))); era.alt.gwas.summ <- subset(era.alt.gwas.summ , shdr.para.east.id!=""); era.alt.gwas.summ
library(purrr)
era.stat.shdr.summ <- list(piLativ.Nota.summ, omega.era.summ,iHSera.summ, sweed.era.summ ,era.alt.gwas.summ) %>% reduce(left_join, by = "shdr.para.east.id"); era.stat.shdr.summ


########## obtain 90th/9th percentile per stat with perms ##########
match.by.range.1k <- function(vec) {
  x <- sim1k.era.df
  if(length(.x <- which(vec >= x$starts & vec <= x$ends ))) .x else NA }

# first create intervals with real shdr
era.shdr.e.int <- intervals::reduce(Intervals(as.matrix(era.e.shdr[,c("start","end")]), closed=c(TRUE,TRUE), type="R"))

# random intervals
sim1k.era.east <- lapply(1:1000, function(x){rand_non_overlapping_intervals(era.genome_size, size(era.shdr.e.int))})

sim.piLativ.Nota.summ <-list()
sim.omega.era.summ <-list()
sim.iHSera.summ <-list()
sim.sweed.era.summ <-list()
sim.era.alt.gwas.summ <-list()


for (i in 1:1000) {
  # prep random intervals for action
  sim1k.era.df <- as.data.frame(sim1k.era.east[i])
  sim1k.era.df$ran.hdr.name <- paste("ran", i, ".hdr.", seq(1:nrow(sim1k.era.df )), sep = "")
 
  # create tmps
  piLativ.Nota.tmp <- piLativ.Nota
  omega.era.tmp <- omega.era
  iHSera.tmp <- iHSera
  sweed.era.tmp <- sweed.era
  era.alt.gwas.tmp <- era.alt.gwas
  
  # add ran hdr names to theta.tmp
  piLativ.Nota.tmp$ran.hdr.name <- sim1k.era.df$ran.hdr.name[mapply(match.by.range.1k, piLativ.Nota.tmp$win_mid_add)]
  omega.era.tmp$ran.hdr.name <- sim1k.era.df$ran.hdr.name[mapply(match.by.range.1k, omega.era.tmp$win_mid_add)]
  iHSera.tmp$ran.hdr.name <- sim1k.era.df$ran.hdr.name[mapply(match.by.range.1k, iHSera.tmp$win_mid_add)]
  sweed.era.tmp$ran.hdr.name <- sim1k.era.df$ran.hdr.name[mapply(match.by.range.1k, sweed.era.tmp$win_mid_add)]
  era.alt.gwas.tmp$ran.hdr.name <- sim1k.era.df$ran.hdr.name[mapply(match.by.range.1k, era.alt.gwas.tmp$BP.wg)]
  
  # use observed genome-wide quantile thresholds
  # get min td in col and ecu
  # subset by hdr being present to avoid getting an NA category
  sim.piLativ.Nota.summ[[i]] <- summarise(group_by(subset(piLativ.Nota.tmp), ran.hdr.name),
                               piNota.minus.piLativ.min=min(piNota.minus.piLativ)); sim.piLativ.Nota.summ <- subset(sim.piLativ.Nota.summ[[i]] , ran.hdr.name!="")
  
  sim.omega.era.summ[[i]] <- summarise(group_by(subset(omega.era.tmp), ran.hdr.name),
                            omega.era.max=max(Omega_MAX_NOTABILIS)); sim.omega.era.summ <- subset(sim.omega.era.summ[[i]] , ran.hdr.name!="")

  sim.iHSera.summ[[i]] <- summarise(group_by(subset(iHSera.tmp), ran.hdr.name),
                         iHS_NOTABILIS.max=max(iHS_NOTABILIS)); sim.iHSera.summ <- subset(sim.iHSera.summ[[i]] , ran.hdr.name!="")

  sim.sweed.era.summ[[i]] <- summarise(group_by(subset(sweed.era.tmp), ran.hdr.name),
                            CLR_NOTABILIS.max=max(CLR_NOTABILIS)); sim.sweed.era.summ <- subset(sim.sweed.era.summ[[i]] , ran.hdr.name!="");

  sim.era.alt.gwas.summ[[i]] <- summarise(group_by(subset(era.alt.gwas.tmp), ran.hdr.name),
                               log10.LRT.pval_median.max=max(-log10(LRT.pval_median))); sim.era.alt.gwas.summ <- subset(sim.era.alt.gwas.summ[[i]] , ran.hdr.name!="")
  }

sim.piLativ.Nota.summ.df <- bind_rows(sim.piLativ.Nota.summ, .id = "column_label")
sim.omega.era.summ.df <- bind_rows(sim.omega.era.summ, .id = "column_label")
sim.iHSera.summ.df <- bind_rows(sim.iHSera.summ, .id = "column_label")
sim.sweed.era.summ.df <- bind_rows(sim.sweed.era.summ, .id = "column_label")
sim.era.alt.gwas.summ.df <- bind_rows(sim.era.alt.gwas.summ, .id = "column_label")

# sim.piLativ.Nota.summ.df <- subset(sim.piLativ.Nota.summ.df, ran.hdr.name!="")
# sim.omega.era.summ.df  <- subset(sim.omega.era.summ.df , ran.hdr.name!="")
# sim.iHSera.summ.df <- subset(sim.iHSera.summ.df, ran.hdr.name!="")
# sim.sweed.era.summ.df  <- subset(sim.sweed.era.summ.df , ran.hdr.name!="")
# sim.era.alt.gwas.summ.df <- subset(sim.era.alt.gwas.summ.df, ran.hdr.name!="")

write.csv(sim.piLativ.Nota.summ.df, "git/2021-altitude-heliconius/data/haplo.stats/sim.1k.piLativ.Nota.summ.df.csv", row.names = F)
write.csv(sim.omega.era.summ.df, "git/2021-altitude-heliconius/data/haplo.stats/sim.1k.omega.era.summ.df.csv", row.names = F)
write.csv(sim.iHSera.summ.df, "git/2021-altitude-heliconius/data/haplo.stats/sim.1k.iHSera.summ.df.csv", row.names = F)
write.csv(sim.sweed.era.summ.df, "git/2021-altitude-heliconius/data/haplo.stats/sim.1k.sweed.era.summ.df.csv", row.names = F)
write.csv(sim.era.alt.gwas.summ.df , "git/2021-altitude-heliconius/data/haplo.stats/sim.1k.era.alt.gwas.summ.df.csv", row.names = F)

# what percentage of the genome is covered by >90th gwas outliers
sim.era.alt.gwas.summ.df  <- read.csv("git/2021-altitude-heliconius/data/haplo.stats/sim.1k.era.alt.gwas.summ.df.csv")
sim.mel.alt.gwas.summ.df  <- read.csv("git/2021-altitude-heliconius/data/haplo.stats/sim.1k.mel.alt.gwas.summ.df.csv")

quantile(sim.era.alt.gwas.summ.df$log10.LRT.pval_median.max, 0.90, na.rm = T); head(era.alt.gwas)
quantile(sim.mel.alt.gwas.summ.df$log10.LRT.pval_median.max, 0.90, na.rm = T)

# careful windows overlap! divide by 5, as 
era.alt.gwas$bp.window <- era.alt.gwas$Position_max- era.alt.gwas$Position_min
(sum(subset(era.alt.gwas, -log10( LRT.pval_median) > quantile(sim.era.alt.gwas.summ.df$log10.LRT.pval_median.max, 0.90, na.rm = T))$bp.window)/5) / era.genome_size *100

era.genome_size; sum(subset(era.alt.gwas)$bp.window)/5

## add whether SHDRs are outliers or not
head(era.stat.shdr.summ)

era.stat.shdr.summ$piNota.minus.piLativ.min.less.5th.perc.sims <- if_else(era.stat.shdr.summ$piNota.minus.piLativ.min < quantile(sim.piLativ.Nota.summ.df$piNota.minus.piLativ.min, 0.05, na.rm = T), "yes", "no"); head(era.stat.shdr.summ)
era.stat.shdr.summ$piNota.minus.piLativ.min.less.10th.perc.sims <- if_else(era.stat.shdr.summ$piNota.minus.piLativ.min < quantile(sim.piLativ.Nota.summ.df$piNota.minus.piLativ.min, 0.1, na.rm = T), "yes", "no"); head(era.stat.shdr.summ)

era.stat.shdr.summ$omega.era.max.more.95th.perc.sims <- if_else(era.stat.shdr.summ$omega.era.max > quantile(sim.omega.era.summ.df$omega.era.max, 0.95, na.rm = T), "yes", "no"); head(era.stat.shdr.summ)
era.stat.shdr.summ$omega.era.max.more.90th.perc.sims <- if_else(era.stat.shdr.summ$omega.era.max > quantile(sim.omega.era.summ.df$omega.era.max, 0.9, na.rm = T), "yes", "no"); head(era.stat.shdr.summ)

era.stat.shdr.summ$iHS_NOTABILIS.max.more.95th.perc.sims <- if_else(era.stat.shdr.summ$iHS_NOTABILIS.max > quantile(sim.iHSera.summ.df$iHS_NOTABILIS.max, 0.95, na.rm = T), "yes", "no"); head(era.stat.shdr.summ)
era.stat.shdr.summ$iHS_NOTABILIS.max.more.90th.perc.sims <- if_else(era.stat.shdr.summ$iHS_NOTABILIS.max > quantile(sim.iHSera.summ.df$iHS_NOTABILIS.max, 0.9, na.rm = T), "yes", "no"); head(era.stat.shdr.summ)

era.stat.shdr.summ$CLR_NOTABILIS.max.more.95th.perc.sims <- if_else(era.stat.shdr.summ$CLR_NOTABILIS.max > quantile(sim.sweed.era.summ.df$CLR_NOTABILIS.max, 0.95, na.rm = T), "yes", "no"); head(era.stat.shdr.summ)
era.stat.shdr.summ$CLR_NOTABILIS.max.more.90th.perc.sims <- if_else(era.stat.shdr.summ$CLR_NOTABILIS.max > quantile(sim.sweed.era.summ.df$CLR_NOTABILIS.max, 0.90, na.rm = T), "yes", "no"); head(era.stat.shdr.summ)

era.stat.shdr.summ$log10.LRT.pval_median.max.95th.perc.sims <- if_else(era.stat.shdr.summ$log10.LRT.pval_median.max > quantile(sim.era.alt.gwas.summ.df$log10.LRT.pval_median.max, 0.95, na.rm = T), "yes", "no"); head(era.stat.shdr.summ)
era.stat.shdr.summ$log10.LRT.pval_median.max.90th.perc.sims <- if_else(era.stat.shdr.summ$log10.LRT.pval_median.max > quantile(sim.era.alt.gwas.summ.df$log10.LRT.pval_median.max, 0.90, na.rm = T), "yes", "no"); head(era.stat.shdr.summ)

######### how many stats (/5) outliers  ########
era.stat.shdr.summ[is.na(era.stat.shdr.summ$piNota.minus.piLativ.min.less.10th.perc.sims),]$piNota.minus.piLativ.min.less.10th.perc.sims <-"no"
era.stat.shdr.summ[is.na(era.stat.shdr.summ$piNota.minus.piLativ.min.less.5th.perc.sims),]$piNota.minus.piLativ.min.less.5th.perc.sims <-"no"

era.stat.shdr.summ$no.stats.outlier.min.max.10th.90th.sims <- if_else(era.stat.shdr.summ$piNota.minus.piLativ.min.less.10th.perc.sims=="yes", 1, 0) +
  if_else(era.stat.shdr.summ$omega.era.max.more.90th.perc.sims=="yes", 1, 0) +
  if_else(era.stat.shdr.summ$iHS_NOTABILIS.max.more.90th.perc.sims=="yes", 1, 0)+
  if_else(era.stat.shdr.summ$CLR_NOTABILIS.max.more.90th.perc.sims=="yes", 1, 0)+
  if_else(era.stat.shdr.summ$log10.LRT.pval_median.max.90th.perc.sims=="yes", 1, 0) ; era.stat.shdr.summ$no.stats.outlier.min.max.10th.90th.sims

era.stat.shdr.summ$no.stats.outlier.min.max.5th.95th.sims <- if_else(era.stat.shdr.summ$piNota.minus.piLativ.min.less.5th.perc.sims=="yes", 1, 0) +
  if_else(era.stat.shdr.summ$omega.era.max.more.95th.perc.sims=="yes", 1, 0) +
  if_else(era.stat.shdr.summ$iHS_NOTABILIS.max.more.95th.perc.sims=="yes", 1, 0)+
  if_else(era.stat.shdr.summ$CLR_NOTABILIS.max.more.95th.perc.sims=="yes", 1, 0)+
  if_else(era.stat.shdr.summ$log10.LRT.pval_median.max.95th.perc.sims=="yes", 1, 0) ; era.stat.shdr.summ$no.stats.outlier.min.max.5th.95th.sims
names(era.stat.shdr.summ)

era.stat.shdr.summ$is.allo.hdr <- era.e.shdr$is.allopatric[match(era.stat.shdr.summ$shdr.para.east.id, era.e.shdr$shdr.para.east.id)]

# save
era.stat.shdr.summ; write.csv(era.stat.shdr.summ, "git/2021-altitude-heliconius/data/haplo.stats/era.stat.shdr.summ.csv", row.names = F)
era.stat.shdr.summ <- read.csv( "git/2021-altitude-heliconius/data/haplo.stats/era.stat.shdr.summ.csv")

# summarise by number of shdr with 1,2,3,4,5
era.stat.shdr.summ.10th.90th.summ <- summarise(group_by(era.stat.shdr.summ,is.allo.hdr, no.stats.outlier.min.max.10th.90th.sims), n=n()); era.stat.shdr.summ.10th.90th.summ
era.stat.shdr.summ.5th.95th.summ <- summarise(group_by(era.stat.shdr.summ,is.allo.hdr, no.stats.outlier.min.max.5th.95th.sims), n=n()); era.stat.shdr.summ.5th.95th.summ


era.stat.shdr.summ
era.e.outlier.number.min.max.10th.90th.plot <- ggplot(data=era.stat.shdr.summ.10th.90th.summ, aes(fill=no.stats.outlier.min.max.10th.90th.sims, x=is.allo.hdr, y=n)) + 
  geom_col(aes(color=is.allo.hdr), size=1)+ 
  ylab("Number of SHDR with haplotagging selection/gwas outliers") + 
  scale_color_manual(values=c( "#049E73","#D65D00" )) +
  scale_y_continuous( expand = c(0,0), limits = c(0,50))+
  scale_fill_continuous(low="white", high="black") +theme_bw()+
  theme(axis.line.x = element_line(color="black", size = 0.5),axis.line.y = element_line(color="black", size = 0.5),axis.title.y = element_text(size=14),
        plot.margin = unit(c(0.5,0.2,0.2,0.2), "lines"),axis.text = element_text(size=14), #axis.text.x = element_blank(),
        axis.text.y = element_text(size=10),axis.title.x = element_blank(), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_rect(fill = "transparent"), plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
        legend.background = element_rect(fill = "transparent"), legend.box.background = element_rect(fill = "transparent"), strip.background.y = element_blank(),
        strip.text = element_text(size=14, face = "bold"), strip.background = element_blank(), legend.position = "none"); era.e.outlier.number.min.max.10th.90th.plot
cowplot::ggsave2(filename  = "git/2021-altitude-heliconius/plots/haplo.stats/era.e.outlier.number.min.max.10th.90th.sims.deltapi.omega.ihs.sweed.altgwas.barplot.png", height = 5, width = 3,dpi=72)


dev.off()
ggplot(data=era.stat.shdr.summ.5th.95th.summ, aes(fill=no.stats.outlier.min.max.5th.95th.sims, x=is.allo.hdr, y=n)) + 
  geom_col(aes(color=is.allo.hdr), size=1)+ 
  ylab("Number of SHDR with haplotagging selection/gwas outliers") + 
  scale_color_manual(values=c( "#049E73","#D65D00" )) +
  scale_y_continuous( expand = c(0,0), limits = c(0,50))+
  scale_fill_continuous(low="white", high="black") +theme_bw()+
  theme(axis.line.x = element_line(color="black", size = 0.5),axis.line.y = element_line(color="black", size = 0.5),axis.title.y = element_text(size=14),
        plot.margin = unit(c(0.5,0.2,0.2,0.2), "lines"),axis.text = element_text(size=14), #axis.text.x = element_blank(),
        axis.text.y = element_text(size=10),axis.title.x = element_blank(), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_rect(fill = "transparent"), plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
        legend.background = element_rect(fill = "transparent"), legend.box.background = element_rect(fill = "transparent"), strip.background.y = element_blank(),
        strip.text = element_text(size=14, face = "bold"), strip.background = element_blank(), legend.position = "none")
cowplot::ggsave2(filename  = "git/2021-altitude-heliconius/plots/haplo.stats/era.e.outlier.number.min.max.5th.95th.sims.deltapi.omega.ihs.sweed.altgwas.barplot.png", height = 5, width = 3,dpi=72)


### do peaks with higher pbs correlate with more stats significant #####
era.stat.shdr.summ <-  read.csv("git/2021-altitude-heliconius/data/haplo.stats/era.stat.shdr.summ.csv")

era.stat.shdr.summ_era.e.shdr <- merge(era.stat.shdr.summ, era.e.shdr, by="shdr.para.east.id"  )
ggplot(data=era.stat.shdr.summ_era.e.shdr, aes(x=as.factor(no.stats.outlier.min.max.10th.90th.sims), y=zPBS.hig.ec.e.max.value, fill=no.stats.outlier.min.max.10th.90th.sims))+
  geom_boxplot()+theme_bw()+theme(legend.position = "none")+ scale_fill_continuous(low="white", high="black") + stat_n_text()
ggplot(data=era.stat.shdr.summ_era.e.shdr, aes(x=as.factor(no.stats.outlier.min.max.10th.90th.sims), y=zPBS.hig.co.e.max.value, fill=no.stats.outlier.min.max.10th.90th.sims))+
  geom_boxplot()+theme_bw()+theme(legend.position = "none")+ scale_fill_continuous(low="white", high="black") + stat_n_text()

ggplot(data=era.stat.shdr.summ_era.e.shdr, aes(x=as.factor(no.stats.outlier.min.max.5th.95th.sims), y=zPBS.hig.ec.e.max.value, fill=no.stats.outlier.min.max.5th.95th.sims))+
  geom_boxplot()+theme_bw()+theme(legend.position = "none")+ scale_fill_continuous(low="white", high="black") + stat_n_text()

ggplot(data=era.stat.shdr.summ_era.e.shdr, aes(x=as.factor(no.stats.outlier.min.max.5th.95th.sims), y=zPBS.hig.co.e.max.value, fill=no.stats.outlier.min.max.5th.95th.sims))+
  geom_boxplot()+theme_bw()+theme(legend.position = "none")+ scale_fill_continuous(low="white", high="black") + stat_n_text()


############################################### 5. mel SIMS ALL STATS against east shdr ###############################################
################# 0. prep data ####################
########## add hdr.id to datasets; east only ##########
head(piMalle.Pless); head(omega.mel); head(iHSmel); head(sweed.mel); head(mel.alt.gwas)
match.by.range <- function(vec) {
  x <- mel.e.shdr
  if(length(.x <- which(vec >= x$start & vec <= x$end ))) .x else NA }

piMalle.Pless$shdr.para.east.id  <-mel.e.shdr$shdr.para.east.id[mapply(match.by.range, piMalle.Pless$win_mid_add)]
omega.mel$shdr.para.east.id  <-mel.e.shdr$shdr.para.east.id[mapply(match.by.range, omega.mel$win_mid_add)]
iHSmel$shdr.para.east.id  <-mel.e.shdr$shdr.para.east.id[mapply(match.by.range, iHSmel$win_mid_add)]
sweed.mel$shdr.para.east.id  <-mel.e.shdr$shdr.para.east.id[mapply(match.by.range, sweed.mel$win_mid_add)]
mel.alt.gwas$shdr.para.east.id  <-mel.e.shdr$shdr.para.east.id[mapply(match.by.range, mel.alt.gwas$BP.wg)]

########## summarise per hdr and concatenate ##########
head(piMalle.Pless)
piMalle.Pless.summ <- summarise(group_by(subset(piMalle.Pless), shdr.para.east.id),
                               piPless.minus.piMalle.min=min(piPless.minus.piMalle)); piMalle.Pless.summ <- subset(piMalle.Pless.summ , shdr.para.east.id!=""); piMalle.Pless.summ

omega.mel.summ <- summarise(group_by(subset(omega.mel), shdr.para.east.id),
                            omega.mel.max=max(Omega_MAX_PLESSENI)); omega.mel.summ <- subset(omega.mel.summ , shdr.para.east.id!=""); omega.mel.summ

iHSmel.summ <- summarise(group_by(subset(iHSmel), shdr.para.east.id),
                         iHS_PLESSENI.max=max(iHS_PLESSENI)); iHSmel.summ <- subset(iHSmel.summ , shdr.para.east.id!=""); iHSmel.summ

sweed.mel.summ <- summarise(group_by(subset(sweed.mel), shdr.para.east.id),
                            CLR_PLESSENI.max=max(CLR_PLESSENI)); sweed.mel.summ <- subset(sweed.mel.summ , shdr.para.east.id!=""); sweed.mel.summ

mel.alt.gwas.summ <- summarise(group_by(subset(mel.alt.gwas), shdr.para.east.id),
                               log10.LRT.pval_median.max=max(-log10(LRT.pval_median))); mel.alt.gwas.summ <- subset(mel.alt.gwas.summ , shdr.para.east.id!=""); mel.alt.gwas.summ
library(purrr)
mel.stat.shdr.summ <- list(piMalle.Pless.summ, omega.mel.summ,iHSmel.summ, sweed.mel.summ ,mel.alt.gwas.summ) %>% reduce(left_join, by = "shdr.para.east.id"); mel.stat.shdr.summ


########## obtain 90th/9th percentile per stat with perms ##########
match.by.range.1k <- function(vec) {
  x <- sim1k.mel.df
  if(length(.x <- which(vec >= x$starts & vec <= x$ends ))) .x else NA }

# first create intervals with real shdr
mel.shdr.e.int <- intervals::reduce(Intervals(as.matrix(mel.e.shdr[,c("start","end")]), closed=c(TRUE,TRUE), type="R"))

# random intervals
sim1k.mel.east <- lapply(1:1000, function(x){rand_non_overlapping_intervals(mel.genome_size, size(mel.shdr.e.int))})

sim.piMalle.Pless.summ <-list()
sim.omega.mel.summ <-list()
sim.iHSmel.summ <-list()
sim.sweed.mel.summ <-list()
sim.mel.alt.gwas.summ <-list()


for (i in 1:1000) {
  # prep random intervals for action
  sim1k.mel.df <- as.data.frame(sim1k.mel.east[i])
  sim1k.mel.df$ran.hdr.name <- paste("ran", i, ".hdr.", seq(1:nrow(sim1k.mel.df )), sep = "")
  
  # create tmps
  piMalle.Pless.tmp <- piMalle.Pless
  omega.mel.tmp <- omega.mel
  iHSmel.tmp <- iHSmel
  sweed.mel.tmp <- sweed.mel
  mel.alt.gwas.tmp <- mel.alt.gwas
  
  # add ran hdr names to theta.tmp
  piMalle.Pless.tmp$ran.hdr.name <- sim1k.mel.df$ran.hdr.name[mapply(match.by.range.1k, piMalle.Pless.tmp$win_mid_add)]
  omega.mel.tmp$ran.hdr.name <- sim1k.mel.df$ran.hdr.name[mapply(match.by.range.1k, omega.mel.tmp$win_mid_add)]
  iHSmel.tmp$ran.hdr.name <- sim1k.mel.df$ran.hdr.name[mapply(match.by.range.1k, iHSmel.tmp$win_mid_add)]
  sweed.mel.tmp$ran.hdr.name <- sim1k.mel.df$ran.hdr.name[mapply(match.by.range.1k, sweed.mel.tmp$win_mid_add)]
  mel.alt.gwas.tmp$ran.hdr.name <- sim1k.mel.df$ran.hdr.name[mapply(match.by.range.1k, mel.alt.gwas.tmp$BP.wg)]
  
  # use observed genome-wide quantile thresholds
  # get min td in col and ecu
  # subset by hdr being present to avoid getting an NA category
  sim.piMalle.Pless.summ[[i]] <- summarise(group_by(subset(piMalle.Pless.tmp), ran.hdr.name),
                                          piPless.minus.piMalle.min=min(piPless.minus.piMalle)); sim.piMalle.Pless.summ[[i]] <- subset(sim.piMalle.Pless.summ[[i]] , ran.hdr.name!="")
  
  sim.omega.mel.summ[[i]] <- summarise(group_by(subset(omega.mel.tmp), ran.hdr.name),
                                       omega.mel.max=max(Omega_MAX_PLESSENI)); sim.omega.mel.summ[[i]] <- subset(sim.omega.mel.summ[[i]] , ran.hdr.name!="")
  
  sim.iHSmel.summ[[i]] <- summarise(group_by(subset(iHSmel.tmp), ran.hdr.name),
                                    iHS_PLESSENI.max=max(iHS_PLESSENI)); sim.iHSmel.summ[[i]] <- subset(sim.iHSmel.summ[[i]] , ran.hdr.name!="")
  
  sim.sweed.mel.summ[[i]] <- summarise(group_by(subset(sweed.mel.tmp), ran.hdr.name),
                                       CLR_PLESSENI.max=max(CLR_PLESSENI)); sim.sweed.mel.summ[[i]] <- subset(sim.sweed.mel.summ[[i]] , ran.hdr.name!="");
  
  sim.mel.alt.gwas.summ[[i]] <- summarise(group_by(subset(mel.alt.gwas.tmp), ran.hdr.name),
                                          log10.LRT.pval_median.max=max(-log10(LRT.pval_median))); sim.mel.alt.gwas.summ[[i]] <- subset(sim.mel.alt.gwas.summ[[i]] , ran.hdr.name!="")
}

sim.piMalle.Pless.summ.df <- bind_rows(sim.piMalle.Pless.summ, .id = "column_label")
sim.omega.mel.summ.df <- bind_rows(sim.omega.mel.summ, .id = "column_label")
sim.iHSmel.summ.df <- bind_rows(sim.iHSmel.summ, .id = "column_label")
sim.sweed.mel.summ.df <- bind_rows(sim.sweed.mel.summ, .id = "column_label")
sim.mel.alt.gwas.summ.df <- bind_rows(sim.mel.alt.gwas.summ, .id = "column_label")

# sim.piMalle.Pless.summ.df <- subset(sim.piMalle.Pless.summ.df, ran.hdr.name!="")
# sim.omega.mel.summ.df  <- subset(sim.omega.mel.summ.df , ran.hdr.name!="")
# sim.iHSmel.summ.df <- subset(sim.iHSmel.summ.df, ran.hdr.name!="")
# sim.sweed.mel.summ.df  <- subset(sim.sweed.mel.summ.df , ran.hdr.name!="")
# sim.mel.alt.gwas.summ.df <- subset(sim.mel.alt.gwas.summ.df, ran.hdr.name!="")

write.csv(sim.piMalle.Pless.summ.df, "git/2021-altitude-heliconius/data/haplo.stats/sim.1k.piMalle.Pless.summ.df.csv", row.names = F)
write.csv(sim.omega.mel.summ.df, "git/2021-altitude-heliconius/data/haplo.stats/sim.1k.omega.mel.summ.df.csv", row.names = F)
write.csv(sim.iHSmel.summ.df, "git/2021-altitude-heliconius/data/haplo.stats/sim.1k.iHSmel.summ.df.csv", row.names = F)
write.csv(sim.sweed.mel.summ.df, "git/2021-altitude-heliconius/data/haplo.stats/sim.1k.sweed.mel.summ.df.csv", row.names = F)
write.csv(sim.mel.alt.gwas.summ.df , "git/2021-altitude-heliconius/data/haplo.stats/sim.1k.mel.alt.gwas.summ.df.csv", row.names = F)


## add whether SHDRs are outliers or not
head(mel.stat.shdr.summ)

mel.stat.shdr.summ$piPless.minus.piMalle.min.less.5th.perc.sims <- if_else(mel.stat.shdr.summ$piPless.minus.piMalle.min < quantile(sim.piMalle.Pless.summ.df$piPless.minus.piMalle.min, 0.05, na.rm = T), "yes", "no"); head(mel.stat.shdr.summ)
mel.stat.shdr.summ$piPless.minus.piMalle.min.less.10th.perc.sims <- if_else(mel.stat.shdr.summ$piPless.minus.piMalle.min < quantile(sim.piMalle.Pless.summ.df$piPless.minus.piMalle.min, 0.1, na.rm = T), "yes", "no"); head(mel.stat.shdr.summ)

mel.stat.shdr.summ$omega.mel.max.more.95th.perc.sims <- if_else(mel.stat.shdr.summ$omega.mel.max > quantile(sim.omega.mel.summ.df$omega.mel.max, 0.95, na.rm = T), "yes", "no"); head(mel.stat.shdr.summ)
mel.stat.shdr.summ$omega.mel.max.more.90th.perc.sims <- if_else(mel.stat.shdr.summ$omega.mel.max > quantile(sim.omega.mel.summ.df$omega.mel.max, 0.9, na.rm = T), "yes", "no"); head(mel.stat.shdr.summ)

mel.stat.shdr.summ$iHS_PLESSENI.max.more.95th.perc.sims <- if_else(mel.stat.shdr.summ$iHS_PLESSENI.max > quantile(sim.iHSmel.summ.df$iHS_PLESSENI.max, 0.95, na.rm = T), "yes", "no"); head(mel.stat.shdr.summ)
mel.stat.shdr.summ$iHS_PLESSENI.max.more.90th.perc.sims <- if_else(mel.stat.shdr.summ$iHS_PLESSENI.max > quantile(sim.iHSmel.summ.df$iHS_PLESSENI.max, 0.9, na.rm = T), "yes", "no"); head(mel.stat.shdr.summ)

mel.stat.shdr.summ$CLR_PLESSENI.max.more.95th.perc.sims <- if_else(mel.stat.shdr.summ$CLR_PLESSENI.max > quantile(sim.sweed.mel.summ.df$CLR_PLESSENI.max, 0.95, na.rm = T), "yes", "no"); head(mel.stat.shdr.summ)
mel.stat.shdr.summ$CLR_PLESSENI.max.more.90th.perc.sims <- if_else(mel.stat.shdr.summ$CLR_PLESSENI.max > quantile(sim.sweed.mel.summ.df$CLR_PLESSENI.max, 0.90, na.rm = T), "yes", "no"); head(mel.stat.shdr.summ)

mel.stat.shdr.summ$log10.LRT.pval_median.max.95th.perc.sims <- if_else(mel.stat.shdr.summ$log10.LRT.pval_median.max > quantile(sim.mel.alt.gwas.summ.df$log10.LRT.pval_median.max, 0.95, na.rm = T), "yes", "no"); head(mel.stat.shdr.summ)
mel.stat.shdr.summ$log10.LRT.pval_median.max.90th.perc.sims <- if_else(mel.stat.shdr.summ$log10.LRT.pval_median.max > quantile(sim.mel.alt.gwas.summ.df$log10.LRT.pval_median.max, 0.90, na.rm = T), "yes", "no"); head(mel.stat.shdr.summ)

######### how many stats (/5) outliers  ########
mel.stat.shdr.summ[is.na(mel.stat.shdr.summ$piPless.minus.piMalle.min.less.10th.perc.sims),]$piPless.minus.piMalle.min.less.10th.perc.sims <-"no"
mel.stat.shdr.summ[is.na(mel.stat.shdr.summ$piPless.minus.piMalle.min.less.5th.perc.sims),]$piPless.minus.piMalle.min.less.5th.perc.sims <-"no"

mel.stat.shdr.summ$no.stats.outlier.min.max.10th.90th.sims <- if_else(mel.stat.shdr.summ$piPless.minus.piMalle.min.less.10th.perc.sims=="yes", 1, 0) +
  if_else(mel.stat.shdr.summ$omega.mel.max.more.90th.perc.sims=="yes", 1, 0) +
  if_else(mel.stat.shdr.summ$iHS_PLESSENI.max.more.90th.perc.sims=="yes", 1, 0)+
  if_else(mel.stat.shdr.summ$CLR_PLESSENI.max.more.90th.perc.sims=="yes", 1, 0)+
  if_else(mel.stat.shdr.summ$log10.LRT.pval_median.max.90th.perc.sims=="yes", 1, 0) ; mel.stat.shdr.summ$no.stats.outlier.min.max.10th.90th.sims

mel.stat.shdr.summ$no.stats.outlier.min.max.5th.95th.sims <- if_else(mel.stat.shdr.summ$piPless.minus.piMalle.min.less.5th.perc.sims=="yes", 1, 0) +
  if_else(mel.stat.shdr.summ$omega.mel.max.more.95th.perc.sims=="yes", 1, 0) +
  if_else(mel.stat.shdr.summ$iHS_PLESSENI.max.more.95th.perc.sims=="yes", 1, 0)+
  if_else(mel.stat.shdr.summ$CLR_PLESSENI.max.more.95th.perc.sims=="yes", 1, 0)+
  if_else(mel.stat.shdr.summ$log10.LRT.pval_median.max.95th.perc.sims=="yes", 1, 0) ; mel.stat.shdr.summ$no.stats.outlier.min.max.5th.95th.sims
names(mel.stat.shdr.summ)

mel.stat.shdr.summ$is.allo.hdr <- mel.e.shdr$is.allopatric[match(mel.stat.shdr.summ$shdr.para.east.id, mel.e.shdr$shdr.para.east.id)]

# save
mel.stat.shdr.summ; write.csv(mel.stat.shdr.summ, "git/2021-altitude-heliconius/data/haplo.stats/mel.stat.shdr.summ.csv", row.names = F)

# summarise by number of shdr with 1,2,3,4,5
mel.stat.shdr.summ.10th.90th.summ <- summarise(group_by(mel.stat.shdr.summ,is.allo.hdr, no.stats.outlier.min.max.10th.90th.sims), n=n()); mel.stat.shdr.summ.10th.90th.summ
mel.stat.shdr.summ.5th.95th.summ <- summarise(group_by(mel.stat.shdr.summ,is.allo.hdr, no.stats.outlier.min.max.5th.95th.sims), n=n()); mel.stat.shdr.summ.5th.95th.summ


mel.stat.shdr.summ
dev.off()
ggplot(data=mel.stat.shdr.summ.5th.95th.summ, aes(fill=no.stats.outlier.min.max.5th.95th.sims, x=is.allo.hdr, y=n)) + 
  geom_col(aes(color=is.allo.hdr), size=1)+ 
  ylab("Number of SHDR with haplotagging selection/gwas outliers") + 
  scale_color_manual(values=c( "#049E73","#D65D00" )) +
  scale_y_continuous( expand = c(0,0), limits = c(0,50))+
  scale_fill_continuous(low="white", high="black") +theme_bw()+
  theme(axis.line.x = element_line(color="black", size = 0.5),axis.line.y = element_line(color="black", size = 0.5),axis.title.y = element_text(size=14),
        plot.margin = unit(c(0.5,0.2,0.2,0.2), "lines"),axis.text = element_text(size=14), #axis.text.x = element_blank(),
        axis.text.y = element_text(size=10),axis.title.x = element_blank(), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_rect(fill = "transparent"), plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
        legend.background = element_rect(fill = "transparent"), legend.box.background = element_rect(fill = "transparent"), strip.background.y = element_blank(),
        strip.text = element_text(size=14, face = "bold"), strip.background = element_blank(), legend.position = "none")
cowplot::ggsave2(filename  = "git/2021-altitude-heliconius/plots/haplo.stats/mel.e.outlier.number.min.max.5th.95th.sims.deltapi.omega.ihs.sweed.altgwas.barplot.png", height = 5, width = 3,dpi=72)


### do peaks with higher pbs correlate with more stats significant #####
mel.stat.shdr.summ <-  read.csv("git/2021-altitude-heliconius/data/haplo.stats/mel.stat.shdr.summ.csv")

mel.stat.shdr.summ_mel.e.shdr <- merge(mel.stat.shdr.summ, mel.e.shdr, by="shdr.para.east.id"  )
head(mel.stat.shdr.summ_mel.e.shdr)
ggplot(data=mel.stat.shdr.summ_mel.e.shdr, aes(x=as.factor(no.stats.outlier.min.max.10th.90th.sims), y=zPBS.hig.ec.e.max.value.x, fill=no.stats.outlier.min.max.10th.90th.sims))+
  geom_boxplot()+theme_bw()+theme(legend.position = "none")+ scale_fill_continuous(low="white", high="black") + stat_n_text()
ggplot(data=mel.stat.shdr.summ_mel.e.shdr, aes(x=as.factor(no.stats.outlier.min.max.10th.90th.sims), y=zPBS.hig.co.e.max.value.x, fill=no.stats.outlier.min.max.10th.90th.sims))+
  geom_boxplot()+theme_bw()+theme(legend.position = "none")+ scale_fill_continuous(low="white", high="black") + stat_n_text()

ggplot(data=mel.stat.shdr.summ_mel.e.shdr, aes(x=as.factor(no.stats.outlier.min.max.5th.95th.sims), y=zPBS.hig.ec.e.max.value.x, fill=no.stats.outlier.min.max.5th.95th.sims))+
  geom_boxplot()+theme_bw()+theme(legend.position = "none")+ scale_fill_continuous(low="white", high="black") + stat_n_text()

ggplot(data=mel.stat.shdr.summ_mel.e.shdr, aes(x=as.factor(no.stats.outlier.min.max.5th.95th.sims), y=zPBS.hig.co.e.max.value.x, fill=no.stats.outlier.min.max.5th.95th.sims))+
  geom_boxplot()+theme_bw()+theme(legend.position = "none")+ scale_fill_continuous(low="white", high="black") + stat_n_text()


############################################### 6. era SIMS ALL STATS against west shdr ###############################################
################# 0. prep data ####################
########## add hdr.id to datasets; west only ##########
head(piLativ.Nota); head(omega.era); head(iHSera); head(sweed.era); head(era.alt.gwas)
match.by.range <- function(vec) {
  x <- era.w.shdr
  if(length(.x <- which(vec >= x$start & vec <= x$end ))) .x else NA }

piLativ.Nota$shdr.para.west.id  <-era.w.shdr$shdr.para.west.id[mapply(match.by.range, piLativ.Nota$win_mid_add)]
omega.era$shdr.para.west.id  <-era.w.shdr$shdr.para.west.id[mapply(match.by.range, omega.era$win_mid_add)]
iHSera$shdr.para.west.id  <-era.w.shdr$shdr.para.west.id[mapply(match.by.range, iHSera$win_mid_add)]
sweed.era$shdr.para.west.id  <-era.w.shdr$shdr.para.west.id[mapply(match.by.range, sweed.era$win_mid_add)]
era.alt.gwas$shdr.para.west.id  <-era.w.shdr$shdr.para.west.id[mapply(match.by.range, era.alt.gwas$BP.wg)]; unique(era.alt.gwas$shdr.para.west.id)

########## summarise per hdr and concatenate ##########
head(piLativ.Nota)
piLativ.Nota.summ <- summarise(group_by(subset(piLativ.Nota), shdr.para.west.id),
                               piNota.minus.piLativ.min=min(piNota.minus.piLativ)); piLativ.Nota.summ <- subset(piLativ.Nota.summ , shdr.para.west.id!=""); piLativ.Nota.summ

omega.era.summ <- summarise(group_by(subset(omega.era), shdr.para.west.id),
                            omega.era.max=max(Omega_MAX_NOTABILIS)); omega.era.summ <- subset(omega.era.summ , shdr.para.west.id!=""); omega.era.summ

iHSera.summ <- summarise(group_by(subset(iHSera), shdr.para.west.id),
                         iHS_NOTABILIS.max=max(iHS_NOTABILIS)); iHSera.summ <- subset(iHSera.summ , shdr.para.west.id!=""); iHSera.summ

sweed.era.summ <- summarise(group_by(subset(sweed.era), shdr.para.west.id),
                            CLR_NOTABILIS.max=max(CLR_NOTABILIS)); sweed.era.summ <- subset(sweed.era.summ , shdr.para.west.id!=""); sweed.era.summ

era.alt.gwas.summ <- summarise(group_by(subset(era.alt.gwas), shdr.para.west.id),
                               log10.LRT.pval_median.max=max(-log10(LRT.pval_median))); era.alt.gwas.summ <- subset(era.alt.gwas.summ , shdr.para.west.id!=""); era.alt.gwas.summ
library(purrr)
era.stat.shdr.summ <- list(piLativ.Nota.summ, omega.era.summ,iHSera.summ, sweed.era.summ ,era.alt.gwas.summ) %>% purrr::reduce(left_join, by = "shdr.para.west.id"); era.stat.shdr.summ


########## obtain 90th/9th percentile per stat with perms ##########
match.by.range.1k <- function(vec) {
  x <- sim1k.era.df
  if(length(.x <- which(vec >= x$starts & vec <= x$ends ))) .x else NA }

# first create intervals with real shdr
era.shdr.e.int <- intervals::reduce(Intervals(as.matrix(era.w.shdr[,c("start","end")]), closed=c(TRUE,TRUE), type="R"))

# random intervals
sim1k.era.west <- lapply(1:1000, function(x){rand_non_overlapping_intervals(era.genome_size, size(era.shdr.e.int))})

sim.piLativ.Nota.summ <-list()
sim.omega.era.summ <-list()
sim.iHSera.summ <-list()
sim.sweed.era.summ <-list()
sim.era.alt.gwas.summ <-list()


for (i in 1:1000) {
  # prep random intervals for action
  sim1k.era.df <- as.data.frame(sim1k.era.west[i])
  sim1k.era.df$ran.hdr.name <- paste("ran", i, ".hdr.", seq(1:nrow(sim1k.era.df )), sep = "")
  
  # create tmps
  piLativ.Nota.tmp <- piLativ.Nota
  omega.era.tmp <- omega.era
  iHSera.tmp <- iHSera
  sweed.era.tmp <- sweed.era
  era.alt.gwas.tmp <- era.alt.gwas
  
  # add ran hdr names to theta.tmp
  piLativ.Nota.tmp$ran.hdr.name <- sim1k.era.df$ran.hdr.name[mapply(match.by.range.1k, piLativ.Nota.tmp$win_mid_add)]
  omega.era.tmp$ran.hdr.name <- sim1k.era.df$ran.hdr.name[mapply(match.by.range.1k, omega.era.tmp$win_mid_add)]
  iHSera.tmp$ran.hdr.name <- sim1k.era.df$ran.hdr.name[mapply(match.by.range.1k, iHSera.tmp$win_mid_add)]
  sweed.era.tmp$ran.hdr.name <- sim1k.era.df$ran.hdr.name[mapply(match.by.range.1k, sweed.era.tmp$win_mid_add)]
  era.alt.gwas.tmp$ran.hdr.name <- sim1k.era.df$ran.hdr.name[mapply(match.by.range.1k, era.alt.gwas.tmp$BP.wg)]
  
  # use observed genome-wide quantile thresholds
  # get min td in col and ecu
  # subset by hdr being present to avoid getting an NA category
  sim.piLativ.Nota.summ[[i]] <- summarise(group_by(subset(piLativ.Nota.tmp), ran.hdr.name),
                                          piNota.minus.piLativ.min=min(piNota.minus.piLativ)); sim.piLativ.Nota.summ[[i]] <- subset(sim.piLativ.Nota.summ[[i]] , ran.hdr.name!="")
  
  sim.omega.era.summ[[i]] <- summarise(group_by(subset(omega.era.tmp), ran.hdr.name),
                                       omega.era.max=max(Omega_MAX_NOTABILIS)); sim.omega.era.summ[[i]] <- subset(sim.omega.era.summ[[i]] , ran.hdr.name!="")
  
  sim.iHSera.summ[[i]] <- summarise(group_by(subset(iHSera.tmp), ran.hdr.name),
                                    iHS_NOTABILIS.max=max(iHS_NOTABILIS)); sim.iHSera.summ[[i]] <- subset(sim.iHSera.summ[[i]] , ran.hdr.name!="")
  
  sim.sweed.era.summ[[i]] <- summarise(group_by(subset(sweed.era.tmp), ran.hdr.name),
                                       CLR_NOTABILIS.max=max(CLR_NOTABILIS)); sim.sweed.era.summ[[i]] <- subset(sim.sweed.era.summ[[i]] , ran.hdr.name!="");
  
  sim.era.alt.gwas.summ[[i]] <- summarise(group_by(subset(era.alt.gwas.tmp), ran.hdr.name),
                                          log10.LRT.pval_median.max=max(-log10(LRT.pval_median))); sim.era.alt.gwas.summ[[i]] <- subset(sim.era.alt.gwas.summ[[i]] , ran.hdr.name!="")
}

sim.piLativ.Nota.summ.df <- bind_rows(sim.piLativ.Nota.summ, .id = "column_label")
sim.omega.era.summ.df <- bind_rows(sim.omega.era.summ, .id = "column_label")
sim.iHSera.summ.df <- bind_rows(sim.iHSera.summ, .id = "column_label")
sim.sweed.era.summ.df <- bind_rows(sim.sweed.era.summ, .id = "column_label")
sim.era.alt.gwas.summ.df <- bind_rows(sim.era.alt.gwas.summ, .id = "column_label")

# sim.piLativ.Nota.summ.df <- subset(sim.piLativ.Nota.summ.df, ran.hdr.name!="")
# sim.omega.era.summ.df  <- subset(sim.omega.era.summ.df , ran.hdr.name!="")
# sim.iHSera.summ.df <- subset(sim.iHSera.summ.df, ran.hdr.name!="")
# sim.sweed.era.summ.df  <- subset(sim.sweed.era.summ.df , ran.hdr.name!="")
# sim.era.alt.gwas.summ.df <- subset(sim.era.alt.gwas.summ.df, ran.hdr.name!="")

write.csv(sim.piLativ.Nota.summ.df, "git/2021-altitude-heliconius/data/haplo.stats/west.hdr/sim.1k.piLativ.Nota.summ.df.csv", row.names = F)
write.csv(sim.omega.era.summ.df, "git/2021-altitude-heliconius/data/haplo.stats/west.hdr/sim.1k.omega.era.summ.df.csv", row.names = F)
write.csv(sim.iHSera.summ.df, "git/2021-altitude-heliconius/data/haplo.stats/west.hdr/sim.1k.iHSera.summ.df.csv", row.names = F)
write.csv(sim.sweed.era.summ.df, "git/2021-altitude-heliconius/data/haplo.stats/west.hdr/sim.1k.sweed.era.summ.df.csv", row.names = F)
write.csv(sim.era.alt.gwas.summ.df , "git/2021-altitude-heliconius/data/haplo.stats/west.hdr/sim.1k.era.alt.gwas.summ.df.csv", row.names = F)


## add whether SHDRs are outliers or not
head(era.stat.shdr.summ)

era.stat.shdr.summ$piNota.minus.piLativ.min.less.5th.perc.sims <- if_else(era.stat.shdr.summ$piNota.minus.piLativ.min < quantile(sim.piLativ.Nota.summ.df$piNota.minus.piLativ.min, 0.05, na.rm = T), "yes", "no"); head(era.stat.shdr.summ)
era.stat.shdr.summ$piNota.minus.piLativ.min.less.10th.perc.sims <- if_else(era.stat.shdr.summ$piNota.minus.piLativ.min < quantile(sim.piLativ.Nota.summ.df$piNota.minus.piLativ.min, 0.1, na.rm = T), "yes", "no"); head(era.stat.shdr.summ)

era.stat.shdr.summ$omega.era.max.more.95th.perc.sims <- if_else(era.stat.shdr.summ$omega.era.max > quantile(sim.omega.era.summ.df$omega.era.max, 0.95, na.rm = T), "yes", "no"); head(era.stat.shdr.summ)
era.stat.shdr.summ$omega.era.max.more.90th.perc.sims <- if_else(era.stat.shdr.summ$omega.era.max > quantile(sim.omega.era.summ.df$omega.era.max, 0.9, na.rm = T), "yes", "no"); head(era.stat.shdr.summ)

era.stat.shdr.summ$iHS_NOTABILIS.max.more.95th.perc.sims <- if_else(era.stat.shdr.summ$iHS_NOTABILIS.max > quantile(sim.iHSera.summ.df$iHS_NOTABILIS.max, 0.95, na.rm = T), "yes", "no"); head(era.stat.shdr.summ)
era.stat.shdr.summ$iHS_NOTABILIS.max.more.90th.perc.sims <- if_else(era.stat.shdr.summ$iHS_NOTABILIS.max > quantile(sim.iHSera.summ.df$iHS_NOTABILIS.max, 0.9, na.rm = T), "yes", "no"); head(era.stat.shdr.summ)

era.stat.shdr.summ$CLR_NOTABILIS.max.more.95th.perc.sims <- if_else(era.stat.shdr.summ$CLR_NOTABILIS.max > quantile(sim.sweed.era.summ.df$CLR_NOTABILIS.max, 0.95, na.rm = T), "yes", "no"); head(era.stat.shdr.summ)
era.stat.shdr.summ$CLR_NOTABILIS.max.more.90th.perc.sims <- if_else(era.stat.shdr.summ$CLR_NOTABILIS.max > quantile(sim.sweed.era.summ.df$CLR_NOTABILIS.max, 0.90, na.rm = T), "yes", "no"); head(era.stat.shdr.summ)

era.stat.shdr.summ$log10.LRT.pval_median.max.95th.perc.sims <- if_else(era.stat.shdr.summ$log10.LRT.pval_median.max > quantile(sim.era.alt.gwas.summ.df$log10.LRT.pval_median.max, 0.95, na.rm = T), "yes", "no"); head(era.stat.shdr.summ)
era.stat.shdr.summ$log10.LRT.pval_median.max.90th.perc.sims <- if_else(era.stat.shdr.summ$log10.LRT.pval_median.max > quantile(sim.era.alt.gwas.summ.df$log10.LRT.pval_median.max, 0.90, na.rm = T), "yes", "no"); head(era.stat.shdr.summ)

######### how many stats (/5) outliers  ########
era.stat.shdr.summ[is.na(era.stat.shdr.summ$piNota.minus.piLativ.min.less.10th.perc.sims),]$piNota.minus.piLativ.min.less.10th.perc.sims <-"no"
era.stat.shdr.summ[is.na(era.stat.shdr.summ$piNota.minus.piLativ.min.less.5th.perc.sims),]$piNota.minus.piLativ.min.less.5th.perc.sims <-"no"

era.stat.shdr.summ[is.na(era.stat.shdr.summ$iHS_NOTABILIS.max.more.90th.perc.sims),]$iHS_NOTABILIS.max.more.90th.perc.sims <-"no"
era.stat.shdr.summ[is.na(era.stat.shdr.summ$iHS_NOTABILIS.max.more.95th.perc.sims),]$iHS_NOTABILIS.max.more.95th.perc.sims <-"no"

era.stat.shdr.summ[is.na(era.stat.shdr.summ$omega.era.max.more.95th.perc.sims),]$omega.era.max.more.95th.perc.sims <-"no"
era.stat.shdr.summ[is.na(era.stat.shdr.summ$omega.era.max.more.90th.perc.sims),]$omega.era.max.more.90th.perc.sims <-"no"

era.stat.shdr.summ$no.stats.outlier.min.max.10th.90th.sims <- if_else(era.stat.shdr.summ$piNota.minus.piLativ.min.less.10th.perc.sims=="yes", 1, 0) +
  if_else(era.stat.shdr.summ$omega.era.max.more.90th.perc.sims=="yes", 1, 0) +
  if_else(era.stat.shdr.summ$iHS_NOTABILIS.max.more.90th.perc.sims=="yes", 1, 0)+
  if_else(era.stat.shdr.summ$CLR_NOTABILIS.max.more.90th.perc.sims=="yes", 1, 0)+
  if_else(era.stat.shdr.summ$log10.LRT.pval_median.max.90th.perc.sims=="yes", 1, 0) ; era.stat.shdr.summ$no.stats.outlier.min.max.10th.90th.sims

era.stat.shdr.summ$no.stats.outlier.min.max.5th.95th.sims <- if_else(era.stat.shdr.summ$piNota.minus.piLativ.min.less.5th.perc.sims=="yes", 1, 0) +
  if_else(era.stat.shdr.summ$omega.era.max.more.95th.perc.sims=="yes", 1, 0) +
  if_else(era.stat.shdr.summ$iHS_NOTABILIS.max.more.95th.perc.sims=="yes", 1, 0)+
  if_else(era.stat.shdr.summ$CLR_NOTABILIS.max.more.95th.perc.sims=="yes", 1, 0)+
  if_else(era.stat.shdr.summ$log10.LRT.pval_median.max.95th.perc.sims=="yes", 1, 0) ; era.stat.shdr.summ$no.stats.outlier.min.max.5th.95th.sims
names(era.stat.shdr.summ)

era.stat.shdr.summ$is.allo.hdr <- era.w.shdr$is.allopatric[match(era.stat.shdr.summ$shdr.para.west.id, era.w.shdr$shdr.para.west.id)]

# save
era.stat.shdr.summ; write.csv(era.stat.shdr.summ, "git/2021-altitude-heliconius/data/haplo.stats/west.hdr/era.stat.shdr.summ.csv", row.names = F)

# summarise by number of shdr with 1,2,3,4,5
era.stat.shdr.summ.10th.90th.summ <- summarise(group_by(era.stat.shdr.summ,is.allo.hdr, no.stats.outlier.min.max.10th.90th.sims), n=n()); era.stat.shdr.summ.10th.90th.summ
era.stat.shdr.summ.5th.95th.summ <- summarise(group_by(era.stat.shdr.summ,is.allo.hdr, no.stats.outlier.min.max.5th.95th.sims), n=n()); era.stat.shdr.summ.5th.95th.summ


era.stat.shdr.summ
dev.off()
ggplot(data=era.stat.shdr.summ.5th.95th.summ, aes(fill=no.stats.outlier.min.max.5th.95th.sims, x=is.allo.hdr, y=n)) + 
  geom_col(aes(color=is.allo.hdr), size=1)+ 
  ylab("Number of SHDR with haplotagging selection/gwas outliers") + 
  scale_color_manual(values=c( "#049E73","#D65D00" )) +
  scale_y_continuous( expand = c(0,0), limits = c(0,100))+
  scale_fill_continuous(low="white", high="black") +theme_bw()+
  theme(axis.line.x = element_line(color="black", size = 0.5),axis.line.y = element_line(color="black", size = 0.5),axis.title.y = element_text(size=14),
        plot.margin = unit(c(0.5,0.2,0.2,0.2), "lines"),axis.text = element_text(size=14), #axis.text.x = element_blank(),
        axis.text.y = element_text(size=10),axis.title.x = element_blank(), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_rect(fill = "transparent"), plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
        legend.background = element_rect(fill = "transparent"), legend.box.background = element_rect(fill = "transparent"), strip.background.y = element_blank(),
        strip.text = element_text(size=14, face = "bold"), strip.background = element_blank(), legend.position = "none")
cowplot::ggsave2(filename  = "git/2021-altitude-heliconius/plots/haplo.stats/west.hdr/era.e.outlier.number.min.max.5th.95th.sims.deltapi.omega.ihs.sweed.altgwas.barplot.png", height = 5, width = 3,dpi=72)

### do peaks with higher pbs correlate with more stats significant #####
era.stat.shdr.summ <-  read.csv("git/2021-altitude-heliconius/data/haplo.stats/west.hdr/era.stat.shdr.summ.csv", row.names = F)

era.stat.shdr.summ_era.w.shdr <- merge(era.stat.shdr.summ, era.w.shdr, by="shdr.para.west.id"  )
ggplot(data=era.stat.shdr.summ_era.w.shdr, aes(x=as.factor(no.stats.outlier.min.max.10th.90th.sims), y=zPBS.hig.ec.w.max.value, fill=no.stats.outlier.min.max.10th.90th.sims))+
  geom_boxplot()+theme_bw()+theme(legend.position = "none")+ scale_fill_continuous(low="white", high="black") 

ggplot(data=era.stat.shdr.summ_era.w.shdr, aes(x=as.factor(no.stats.outlier.min.max.10th.90th.sims), y=zPBS.hig.ec.w.max.value, fill=no.stats.outlier.min.max.10th.90th.sims))+
  geom_boxplot()+theme_bw()+theme(legend.position = "none")+ scale_fill_continuous(low="white", high="black") 


############################################### 7. mel SIMS ALL STATS against west shdr ###############################################
################# 0. prep data ####################
########## add hdr.id to datasets; west only ##########
head(piMalle.Pless); head(omega.mel); head(iHSmel); head(sweed.mel); head(mel.alt.gwas)
match.by.range <- function(vec) {
  x <- mel.w.shdr
  if(length(.x <- which(vec >= x$start & vec <= x$end ))) .x else NA }

piMalle.Pless$shdr.para.west.id  <-mel.w.shdr$shdr.para.west.id[mapply(match.by.range, piMalle.Pless$win_mid_add)]
omega.mel$shdr.para.west.id  <-mel.w.shdr$shdr.para.west.id[mapply(match.by.range, omega.mel$win_mid_add)]
iHSmel$shdr.para.west.id  <-mel.w.shdr$shdr.para.west.id[mapply(match.by.range, iHSmel$win_mid_add)]
sweed.mel$shdr.para.west.id  <-mel.w.shdr$shdr.para.west.id[mapply(match.by.range, sweed.mel$win_mid_add)]
mel.alt.gwas$shdr.para.west.id  <-mel.w.shdr$shdr.para.west.id[mapply(match.by.range, mel.alt.gwas$BP.wg)]

########## summarise per hdr and concatenate ##########
head(piMalle.Pless)
piMalle.Pless.summ <- summarise(group_by(subset(piMalle.Pless), shdr.para.west.id),
                                piPless.minus.piMalle.min=min(piPless.minus.piMalle)); piMalle.Pless.summ <- subset(piMalle.Pless.summ , shdr.para.west.id!=""); piMalle.Pless.summ

omega.mel.summ <- summarise(group_by(subset(omega.mel), shdr.para.west.id),
                            omega.mel.max=max(Omega_MAX_PLESSENI)); omega.mel.summ <- subset(omega.mel.summ , shdr.para.west.id!=""); omega.mel.summ

iHSmel.summ <- summarise(group_by(subset(iHSmel), shdr.para.west.id),
                         iHS_PLESSENI.max=max(iHS_PLESSENI)); iHSmel.summ <- subset(iHSmel.summ , shdr.para.west.id!=""); iHSmel.summ

sweed.mel.summ <- summarise(group_by(subset(sweed.mel), shdr.para.west.id),
                            CLR_PLESSENI.max=max(CLR_PLESSENI)); sweed.mel.summ <- subset(sweed.mel.summ , shdr.para.west.id!=""); sweed.mel.summ

mel.alt.gwas.summ <- summarise(group_by(subset(mel.alt.gwas), shdr.para.west.id),
                               log10.LRT.pval_median.max=max(-log10(LRT.pval_median))); mel.alt.gwas.summ <- subset(mel.alt.gwas.summ , shdr.para.west.id!=""); mel.alt.gwas.summ
library(purrr)
mel.stat.shdr.summ <- list(piMalle.Pless.summ, omega.mel.summ,iHSmel.summ, sweed.mel.summ ,mel.alt.gwas.summ) %>% reduce(left_join, by = "shdr.para.west.id"); mel.stat.shdr.summ


########## obtain 90th/9th percentile per stat with perms ##########
match.by.range.1k <- function(vec) {
  x <- sim1k.mel.df
  if(length(.x <- which(vec >= x$starts & vec <= x$ends ))) .x else NA }

# first create intervals with real shdr
mel.shdr.e.int <- intervals::reduce(Intervals(as.matrix(mel.w.shdr[,c("start","end")]), closed=c(TRUE,TRUE), type="R"))

# random intervals
sim1k.mel.west <- lapply(1:1000, function(x){rand_non_overlapping_intervals(mel.genome_size, size(mel.shdr.e.int))})

sim.piMalle.Pless.summ <-list()
sim.omega.mel.summ <-list()
sim.iHSmel.summ <-list()
sim.sweed.mel.summ <-list()
sim.mel.alt.gwas.summ <-list()


for (i in 1:1000) {
  # prep random intervals for action
  sim1k.mel.df <- as.data.frame(sim1k.mel.west[i])
  sim1k.mel.df$ran.hdr.name <- paste("ran", i, ".hdr.", seq(1:nrow(sim1k.mel.df )), sep = "")
  
  # create tmps
  piMalle.Pless.tmp <- piMalle.Pless
  omega.mel.tmp <- omega.mel
  iHSmel.tmp <- iHSmel
  sweed.mel.tmp <- sweed.mel
  mel.alt.gwas.tmp <- mel.alt.gwas
  
  # add ran hdr names to theta.tmp
  piMalle.Pless.tmp$ran.hdr.name <- sim1k.mel.df$ran.hdr.name[mapply(match.by.range.1k, piMalle.Pless.tmp$win_mid_add)]
  omega.mel.tmp$ran.hdr.name <- sim1k.mel.df$ran.hdr.name[mapply(match.by.range.1k, omega.mel.tmp$win_mid_add)]
  iHSmel.tmp$ran.hdr.name <- sim1k.mel.df$ran.hdr.name[mapply(match.by.range.1k, iHSmel.tmp$win_mid_add)]
  sweed.mel.tmp$ran.hdr.name <- sim1k.mel.df$ran.hdr.name[mapply(match.by.range.1k, sweed.mel.tmp$win_mid_add)]
  mel.alt.gwas.tmp$ran.hdr.name <- sim1k.mel.df$ran.hdr.name[mapply(match.by.range.1k, mel.alt.gwas.tmp$BP.wg)]
  
  # use observed genome-wide quantile thresholds
  # get min td in col and ecu
  # subset by hdr being present to avoid getting an NA category
  sim.piMalle.Pless.summ[[i]] <- summarise(group_by(subset(piMalle.Pless.tmp), ran.hdr.name),
                                           piPless.minus.piMalle.min=min(piPless.minus.piMalle)); sim.piMalle.Pless.summ[[i]] <- subset(sim.piMalle.Pless.summ[[i]] , ran.hdr.name!="")
  
  sim.omega.mel.summ[[i]] <- summarise(group_by(subset(omega.mel.tmp), ran.hdr.name),
                                       omega.mel.max=max(Omega_MAX_PLESSENI)); sim.omega.mel.summ[[i]] <- subset(sim.omega.mel.summ[[i]] , ran.hdr.name!="")
  
  sim.iHSmel.summ[[i]] <- summarise(group_by(subset(iHSmel.tmp), ran.hdr.name),
                                    iHS_PLESSENI.max=max(iHS_PLESSENI)); sim.iHSmel.summ[[i]] <- subset(sim.iHSmel.summ[[i]] , ran.hdr.name!="")
  
  sim.sweed.mel.summ[[i]] <- summarise(group_by(subset(sweed.mel.tmp), ran.hdr.name),
                                       CLR_PLESSENI.max=max(CLR_PLESSENI)); sim.sweed.mel.summ[[i]] <- subset(sim.sweed.mel.summ[[i]] , ran.hdr.name!="");
  
  sim.mel.alt.gwas.summ[[i]] <- summarise(group_by(subset(mel.alt.gwas.tmp), ran.hdr.name),
                                          log10.LRT.pval_median.max=max(-log10(LRT.pval_median))); sim.mel.alt.gwas.summ[[i]] <- subset(sim.mel.alt.gwas.summ[[i]] , ran.hdr.name!="")
}

sim.piMalle.Pless.summ.df <- bind_rows(sim.piMalle.Pless.summ, .id = "column_label")
sim.omega.mel.summ.df <- bind_rows(sim.omega.mel.summ, .id = "column_label")
sim.iHSmel.summ.df <- bind_rows(sim.iHSmel.summ, .id = "column_label")
sim.sweed.mel.summ.df <- bind_rows(sim.sweed.mel.summ, .id = "column_label")
sim.mel.alt.gwas.summ.df <- bind_rows(sim.mel.alt.gwas.summ, .id = "column_label")

# sim.piMalle.Pless.summ.df <- subset(sim.piMalle.Pless.summ.df, ran.hdr.name!="")
# sim.omega.mel.summ.df  <- subset(sim.omega.mel.summ.df , ran.hdr.name!="")
# sim.iHSmel.summ.df <- subset(sim.iHSmel.summ.df, ran.hdr.name!="")
# sim.sweed.mel.summ.df  <- subset(sim.sweed.mel.summ.df , ran.hdr.name!="")
# sim.mel.alt.gwas.summ.df <- subset(sim.mel.alt.gwas.summ.df, ran.hdr.name!="")

write.csv(sim.piMalle.Pless.summ.df, "git/2021-altitude-heliconius/data/haplo.stats/west.hdr/sim.1k.piMalle.Pless.summ.df.csv", row.names = F)
write.csv(sim.omega.mel.summ.df, "git/2021-altitude-heliconius/data/haplo.stats/west.hdr/sim.1k.omega.mel.summ.df.csv", row.names = F)
write.csv(sim.iHSmel.summ.df, "git/2021-altitude-heliconius/data/haplo.stats/west.hdr/sim.1k.iHSmel.summ.df.csv", row.names = F)
write.csv(sim.sweed.mel.summ.df, "git/2021-altitude-heliconius/data/haplo.stats/west.hdr/sim.1k.sweed.mel.summ.df.csv", row.names = F)
write.csv(sim.mel.alt.gwas.summ.df , "git/2021-altitude-heliconius/data/haplo.stats/west.hdr/sim.1k.mel.alt.gwas.summ.df.csv", row.names = F)


## add whether SHDRs are outliers or not
head(mel.stat.shdr.summ)

mel.stat.shdr.summ$piPless.minus.piMalle.min.less.5th.perc.sims <- if_else(mel.stat.shdr.summ$piPless.minus.piMalle.min < quantile(sim.piMalle.Pless.summ.df$piPless.minus.piMalle.min, 0.05, na.rm = T), "yes", "no"); head(mel.stat.shdr.summ)
mel.stat.shdr.summ$piPless.minus.piMalle.min.less.10th.perc.sims <- if_else(mel.stat.shdr.summ$piPless.minus.piMalle.min < quantile(sim.piMalle.Pless.summ.df$piPless.minus.piMalle.min, 0.1, na.rm = T), "yes", "no"); head(mel.stat.shdr.summ)

mel.stat.shdr.summ$omega.mel.max.more.95th.perc.sims <- if_else(mel.stat.shdr.summ$omega.mel.max > quantile(sim.omega.mel.summ.df$omega.mel.max, 0.95, na.rm = T), "yes", "no"); head(mel.stat.shdr.summ)
mel.stat.shdr.summ$omega.mel.max.more.90th.perc.sims <- if_else(mel.stat.shdr.summ$omega.mel.max > quantile(sim.omega.mel.summ.df$omega.mel.max, 0.9, na.rm = T), "yes", "no"); head(mel.stat.shdr.summ)

mel.stat.shdr.summ$iHS_PLESSENI.max.more.95th.perc.sims <- if_else(mel.stat.shdr.summ$iHS_PLESSENI.max > quantile(sim.iHSmel.summ.df$iHS_PLESSENI.max, 0.95, na.rm = T), "yes", "no"); head(mel.stat.shdr.summ)
mel.stat.shdr.summ$iHS_PLESSENI.max.more.90th.perc.sims <- if_else(mel.stat.shdr.summ$iHS_PLESSENI.max > quantile(sim.iHSmel.summ.df$iHS_PLESSENI.max, 0.9, na.rm = T), "yes", "no"); head(mel.stat.shdr.summ)

mel.stat.shdr.summ$CLR_PLESSENI.max.more.95th.perc.sims <- if_else(mel.stat.shdr.summ$CLR_PLESSENI.max > quantile(sim.sweed.mel.summ.df$CLR_PLESSENI.max, 0.95, na.rm = T), "yes", "no"); head(mel.stat.shdr.summ)
mel.stat.shdr.summ$CLR_PLESSENI.max.more.90th.perc.sims <- if_else(mel.stat.shdr.summ$CLR_PLESSENI.max > quantile(sim.sweed.mel.summ.df$CLR_PLESSENI.max, 0.90, na.rm = T), "yes", "no"); head(mel.stat.shdr.summ)

mel.stat.shdr.summ$log10.LRT.pval_median.max.95th.perc.sims <- if_else(mel.stat.shdr.summ$log10.LRT.pval_median.max > quantile(sim.mel.alt.gwas.summ.df$log10.LRT.pval_median.max, 0.95, na.rm = T), "yes", "no"); head(mel.stat.shdr.summ)
mel.stat.shdr.summ$log10.LRT.pval_median.max.90th.perc.sims <- if_else(mel.stat.shdr.summ$log10.LRT.pval_median.max > quantile(sim.mel.alt.gwas.summ.df$log10.LRT.pval_median.max, 0.90, na.rm = T), "yes", "no"); head(mel.stat.shdr.summ)

######### how many stats (/5) outliers  ########
mel.stat.shdr.summ[is.na(mel.stat.shdr.summ$piPless.minus.piMalle.min.less.10th.perc.sims),]$piPless.minus.piMalle.min.less.10th.perc.sims <-"no"
mel.stat.shdr.summ[is.na(mel.stat.shdr.summ$piPless.minus.piMalle.min.less.5th.perc.sims),]$piPless.minus.piMalle.min.less.5th.perc.sims <-"no"

mel.stat.shdr.summ$no.stats.outlier.min.max.10th.90th.sims <- if_else(mel.stat.shdr.summ$piPless.minus.piMalle.min.less.10th.perc.sims=="yes", 1, 0) +
  if_else(mel.stat.shdr.summ$omega.mel.max.more.90th.perc.sims=="yes", 1, 0) +
  if_else(mel.stat.shdr.summ$iHS_PLESSENI.max.more.90th.perc.sims=="yes", 1, 0)+
  if_else(mel.stat.shdr.summ$CLR_PLESSENI.max.more.90th.perc.sims=="yes", 1, 0)+
  if_else(mel.stat.shdr.summ$log10.LRT.pval_median.max.90th.perc.sims=="yes", 1, 0) ; mel.stat.shdr.summ$no.stats.outlier.min.max.10th.90th.sims

mel.stat.shdr.summ$no.stats.outlier.min.max.5th.95th.sims <- if_else(mel.stat.shdr.summ$piPless.minus.piMalle.min.less.5th.perc.sims=="yes", 1, 0) +
  if_else(mel.stat.shdr.summ$omega.mel.max.more.95th.perc.sims=="yes", 1, 0) +
  if_else(mel.stat.shdr.summ$iHS_PLESSENI.max.more.95th.perc.sims=="yes", 1, 0)+
  if_else(mel.stat.shdr.summ$CLR_PLESSENI.max.more.95th.perc.sims=="yes", 1, 0)+
  if_else(mel.stat.shdr.summ$log10.LRT.pval_median.max.95th.perc.sims=="yes", 1, 0) ; mel.stat.shdr.summ$no.stats.outlier.min.max.5th.95th.sims
names(mel.stat.shdr.summ)

mel.stat.shdr.summ$is.allo.hdr <- mel.w.shdr$is.allopatric[match(mel.stat.shdr.summ$shdr.para.west.id, mel.w.shdr$shdr.para.west.id)]

# save
mel.stat.shdr.summ; write.csv(mel.stat.shdr.summ, "git/2021-altitude-heliconius/data/haplo.stats/west.hdr/mel.stat.shdr.summ.csv", row.names = F)
mel.stat.shdr.summ <- read.csv( "git/2021-altitude-heliconius/data/haplo.stats/mel.stat.shdr.summ.csv")

# summarise by number of shdr with 1,2,3,4,5
mel.stat.shdr.summ.10th.90th.summ <- summarise(group_by(mel.stat.shdr.summ,is.allo.hdr, no.stats.outlier.min.max.10th.90th.sims), n=n()); mel.stat.shdr.summ.10th.90th.summ


mel.stat.shdr.summ
mel.e.outlier.number.min.max.10th.90th.plot <- ggplot(data=mel.stat.shdr.summ.10th.90th.summ, aes(fill=no.stats.outlier.min.max.10th.90th.sims, x=is.allo.hdr, y=n)) + 
  geom_col(aes(color=is.allo.hdr), size=1)+ 
  ylab("Number of SHDR with haplotagging selection/gwas outliers") + 
  scale_color_manual(values=c( "#049E73","#D65D00" )) +
  scale_y_continuous( expand = c(0,0), limits = c(0,50))+
  scale_fill_continuous(low="white", high="black") +theme_bw()+
  theme(axis.line.x = element_line(color="black", size = 0.5),axis.line.y = element_line(color="black", size = 0.5),axis.title.y = element_text(size=14),
        plot.margin = unit(c(0.5,0.2,0.2,0.2), "lines"),axis.text = element_text(size=14), #axis.text.x = element_blank(),
        axis.text.y = element_text(size=10),axis.title.x = element_blank(), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_rect(fill = "transparent"), plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
        legend.background = element_rect(fill = "transparent"), legend.box.background = element_rect(fill = "transparent"), strip.background.y = element_blank(),
        strip.text = element_text(size=14, face = "bold"), strip.background = element_blank(), legend.position = "none"); mel.e.outlier.number.min.max.10th.90th.plot

cowplot::ggsave2(filename  = "git/2021-altitude-heliconius/plots/haplo.stats/mel.e.outlier.number.min.max.10th.90th.sims.deltapi.omega.ihs.sweed.altgwas.barplot.png", height = 5, width = 3,dpi=72)


dev.off()
ggplot(data=mel.stat.shdr.summ.5th.95th.summ, aes(fill=no.stats.outlier.min.max.5th.95th.sims, x=is.allo.hdr, y=n)) + 
  geom_col(aes(color=is.allo.hdr), size=1)+ 
  ylab("Number of SHDR with haplotagging selection/gwas outliers") + 
  scale_color_manual(values=c( "#049E73","#D65D00" )) +
  scale_y_continuous( expand = c(0,0), limits = c(0,50))+
  scale_fill_continuous(low="white", high="black") +theme_bw()+
  theme(axis.line.x = element_line(color="black", size = 0.5),axis.line.y = element_line(color="black", size = 0.5),axis.title.y = element_text(size=14),
        plot.margin = unit(c(0.5,0.2,0.2,0.2), "lines"),axis.text = element_text(size=14), #axis.text.x = element_blank(),
        axis.text.y = element_text(size=10),axis.title.x = element_blank(), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_rect(fill = "transparent"), plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
        legend.background = element_rect(fill = "transparent"), legend.box.background = element_rect(fill = "transparent"), strip.background.y = element_blank(),
        strip.text = element_text(size=14, face = "bold"), strip.background = element_blank(), legend.position = "none")
cowplot::ggsave2(filename  = "git/2021-altitude-heliconius/plots/haplo.stats/west.hdr/mel.e.outlier.number.min.max.5th.95th.sims.deltapi.omega.ihs.sweed.altgwas.barplot.png", height = 5, width = 3,dpi=72)

############################################### 8. plot era mel outlier barplot for S.I. ###############################################
### load all data ####
era.stat.shdr.east.summ <- read.csv( "git/2021-altitude-heliconius/data/haplo.stats/era.stat.shdr.summ.csv")
era.stat.shdr.west.summ <- read.csv( "git/2021-altitude-heliconius/data/haplo.stats/west.hdr/era.stat.shdr.summ.csv")

mel.stat.shdr.east.summ <- read.csv( "git/2021-altitude-heliconius/data/haplo.stats/mel.stat.shdr.summ.csv")
mel.stat.shdr.west.summ <- read.csv( "git/2021-altitude-heliconius/data/haplo.stats/west.hdr/mel.stat.shdr.summ.csv")

### summarise by number of shdr with 1,2,3,4,5 ####
era.stat.shdr.east.summ.10th.90th.summ <- summarise(group_by(era.stat.shdr.east.summ,is.allo.hdr, no.stats.outlier.min.max.10th.90th.sims), n=n()); era.stat.shdr.east.summ.10th.90th.summ
era.stat.shdr.west.summ.10th.90th.summ <- summarise(group_by(era.stat.shdr.west.summ,is.allo.hdr, no.stats.outlier.min.max.10th.90th.sims), n=n()); era.stat.shdr.west.summ.10th.90th.summ

mel.stat.shdr.east.summ.10th.90th.summ <- summarise(group_by(mel.stat.shdr.east.summ,is.allo.hdr, no.stats.outlier.min.max.10th.90th.sims), n=n()); mel.stat.shdr.east.summ.10th.90th.summ
mel.stat.shdr.west.summ.10th.90th.summ <- summarise(group_by(mel.stat.shdr.west.summ,is.allo.hdr, no.stats.outlier.min.max.10th.90th.sims), n=n()); mel.stat.shdr.west.summ.10th.90th.summ


(sum(subset(era.stat.shdr.east.summ.10th.90th.summ, no.stats.outlier.min.max.10th.90th.sims>0 &is.allo.hdr=="no")$n)/ sum(subset(era.stat.shdr.east.summ.10th.90th.summ,is.allo.hdr=="no")$n))*100
(sum(subset(era.stat.shdr.west.summ.10th.90th.summ, no.stats.outlier.min.max.10th.90th.sims>0 &is.allo.hdr=="no")$n)/ sum(subset(era.stat.shdr.west.summ.10th.90th.summ,is.allo.hdr=="no")$n))*100

(sum(subset(mel.stat.shdr.east.summ.10th.90th.summ, no.stats.outlier.min.max.10th.90th.sims>0 &is.allo.hdr=="no")$n)/ sum(subset(mel.stat.shdr.east.summ.10th.90th.summ,is.allo.hdr=="no")$n))*100
(sum(subset(mel.stat.shdr.west.summ.10th.90th.summ, no.stats.outlier.min.max.10th.90th.sims>0 &is.allo.hdr=="no")$n)/ sum(subset(mel.stat.shdr.west.summ.10th.90th.summ,is.allo.hdr=="no")$n))*100

(51.1+65.8)/2; sd(c(51.1,65.8))
(26.80412+37.63441)/2; sd(c(26.80412,37.63441))

(sum(subset(era.stat.shdr.east.summ.10th.90th.summ, no.stats.outlier.min.max.10th.90th.sims>0 &is.allo.hdr=="yes")$n)/ sum(subset(era.stat.shdr.east.summ.10th.90th.summ,is.allo.hdr=="yes")$n))*100
(sum(subset(era.stat.shdr.west.summ.10th.90th.summ, no.stats.outlier.min.max.10th.90th.sims>0 &is.allo.hdr=="yes")$n)/ sum(subset(era.stat.shdr.west.summ.10th.90th.summ,is.allo.hdr=="yes")$n))*100

(sum(subset(mel.stat.shdr.east.summ.10th.90th.summ, no.stats.outlier.min.max.10th.90th.sims>0 &is.allo.hdr=="yes")$n)/ sum(subset(mel.stat.shdr.east.summ.10th.90th.summ,is.allo.hdr=="yes")$n))*100
(sum(subset(mel.stat.shdr.west.summ.10th.90th.summ, no.stats.outlier.min.max.10th.90th.sims>0 &is.allo.hdr=="yes")$n)/ sum(subset(mel.stat.shdr.west.summ.10th.90th.summ,is.allo.hdr=="yes")$n))*100
(42.5+40+54.5 + 72.7)/4; sd(c(42.5,40,54.5 , 72.7))


(sum(subset(era.stat.shdr.west.summ.10th.90th.summ, no.stats.outlier.min.max.10th.90th.sims>0)$n)/ sum(era.stat.shdr.west.summ.10th.90th.summ$n))*100
(sum(subset(mel.stat.shdr.east.summ.10th.90th.summ, no.stats.outlier.min.max.10th.90th.sims>0)$n)/ sum(mel.stat.shdr.east.summ.10th.90th.summ$n))*100
(sum(subset(mel.stat.shdr.west.summ.10th.90th.summ, no.stats.outlier.min.max.10th.90th.sims>0)$n)/ sum(mel.stat.shdr.west.summ.10th.90th.summ$n))*100
# percentage of total (lumped allo para)
era.stat.shdr.east.summ.10th.90th.summ1 <- summarise(group_by(era.stat.shdr.east.summ, no.stats.outlier.min.max.10th.90th.sims), n=n()); era.stat.shdr.east.summ.10th.90th.summ1
era.stat.shdr.west.summ.10th.90th.summ1 <- summarise(group_by(era.stat.shdr.west.summ, no.stats.outlier.min.max.10th.90th.sims), n=n()); era.stat.shdr.west.summ.10th.90th.summ1

mel.stat.shdr.east.summ.10th.90th.summ1 <- summarise(group_by(mel.stat.shdr.east.summ, no.stats.outlier.min.max.10th.90th.sims), n=n()); mel.stat.shdr.east.summ.10th.90th.summ1
mel.stat.shdr.west.summ.10th.90th.summ1 <- summarise(group_by(mel.stat.shdr.west.summ, no.stats.outlier.min.max.10th.90th.sims), n=n()); mel.stat.shdr.west.summ.10th.90th.summ1

(sum(subset(era.stat.shdr.east.summ.10th.90th.summ1, no.stats.outlier.min.max.10th.90th.sims>0)$n)/ sum(era.stat.shdr.east.summ.10th.90th.summ1$n))*100
(sum(subset(era.stat.shdr.west.summ.10th.90th.summ1, no.stats.outlier.min.max.10th.90th.sims>0)$n)/ sum(era.stat.shdr.west.summ.10th.90th.summ1$n))*100

(sum(subset(mel.stat.shdr.east.summ.10th.90th.summ1, no.stats.outlier.min.max.10th.90th.sims>0)$n)/ sum(mel.stat.shdr.east.summ.10th.90th.summ1$n))*100
(sum(subset(mel.stat.shdr.west.summ.10th.90th.summ1, no.stats.outlier.min.max.10th.90th.sims>0)$n)/ sum(mel.stat.shdr.west.summ.10th.90th.summ1$n))*100


### plot ####
library(patchwork)
era.e.outlier.number.min.max.10th.90th.plot <- ggplot(data=era.stat.shdr.east.summ.10th.90th.summ, aes(fill=no.stats.outlier.min.max.10th.90th.sims, x=is.allo.hdr, y=n)) + 
  geom_col(aes(color=is.allo.hdr), size=1)+ 
  ylab("") + 
  scale_color_manual(values=c( "#049E73","#D65D00" )) +
  scale_y_continuous( expand = c(0,0), limits = c(0,50))+
  scale_x_discrete(labels= c("Parapatric", "Allopatric")) +
  scale_fill_continuous(low="white", high="black")  +
  theme(axis.title.y = element_text(size=13),
        plot.margin = unit(c(0.5,0.2,0.2,0.2), "lines"),axis.text = element_text(size=14), #axis.text.x = element_blank(),
        axis.text.y = element_text(size=10),axis.title.x = element_blank(), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_rect(fill = "transparent"), plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        legend.background = element_rect(fill = "transparent"), legend.box.background = element_rect(fill = "transparent"), strip.background.y = element_blank(),
        strip.text = element_text(size=14, face = "bold"), strip.background = element_blank(), legend.position = "none"); era.e.outlier.number.min.max.10th.90th.plot

era.w.outlier.number.min.max.10th.90th.plot <- ggplot(data=era.stat.shdr.west.summ.10th.90th.summ, aes(fill=no.stats.outlier.min.max.10th.90th.sims, x=is.allo.hdr, y=n)) + 
  geom_col(aes(color=is.allo.hdr), size=1)+ 
  ylab("") + 
  scale_color_manual(values=c( "#0372B2","#D65D00" )) +
  scale_x_discrete(labels= c("Parapatric", "Allopatric")) +
  scale_y_continuous( expand = c(0,0), limits = c(0,100))+
  scale_fill_continuous(low="white", high="black")  +
  theme(axis.title.y = element_text(size=13),
        plot.margin = unit(c(0.5,0.2,0.2,0.2), "lines"),axis.text = element_text(size=14), #axis.text.x = element_blank(),
        axis.text.y = element_text(size=10),axis.title.x = element_blank(), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_rect(fill = "transparent"), plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        legend.background = element_rect(fill = "transparent"), legend.box.background = element_rect(fill = "transparent"), strip.background.y = element_blank(),
        strip.text = element_text(size=14, face = "bold"), strip.background = element_blank(), legend.position = "none"); era.w.outlier.number.min.max.10th.90th.plot

mel.e.outlier.number.min.max.10th.90th.plot <- ggplot(data=mel.stat.shdr.east.summ.10th.90th.summ, aes(fill=no.stats.outlier.min.max.10th.90th.sims, x=is.allo.hdr, y=n)) + 
  geom_col(aes(color=is.allo.hdr), size=1)+ 
  ylab("") + 
  scale_color_manual(values=c( "#049E73","#D65D00" )) +
  scale_y_continuous( expand = c(0,0), limits = c(0,50))+  
  scale_x_discrete(labels= c("Parapatric", "Allopatric")) +
  scale_fill_continuous(low="white", high="black")  +
  theme(axis.title.y = element_text(size=13),
        plot.margin = unit(c(0.5,0.2,0.2,0.2), "lines"),axis.text = element_text(size=14), #axis.text.x = element_blank(),
        axis.text.y = element_text(size=10),axis.title.x = element_blank(), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_rect(fill = "transparent"), plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        legend.background = element_rect(fill = "transparent"), legend.box.background = element_rect(fill = "transparent"), strip.background.y = element_blank(),
        strip.text = element_text(size=14, face = "bold"), strip.background = element_blank(), legend.position = "none"); mel.e.outlier.number.min.max.10th.90th.plot

mel.w.outlier.number.min.max.10th.90th.plot <- ggplot(data=mel.stat.shdr.west.summ.10th.90th.summ, aes(fill=no.stats.outlier.min.max.10th.90th.sims, x=is.allo.hdr, y=n)) + 
  geom_col(aes(color=is.allo.hdr), size=1)+ 
  ylab("Number of SHDR with haplotagging selection/gwas outliers") + 
  scale_color_manual(values=c( "#0372B2","#D65D00" )) +
  scale_x_discrete(labels= c("Parapatric", "Allopatric")) +
  scale_y_continuous( expand = c(0,0), limits = c(0,100))+
  scale_fill_continuous(low="white", high="black") +
  theme(axis.title.y = element_text(size=13),
        plot.margin = unit(c(0.5,0.2,0.2,0.2), "lines"),axis.text = element_text(size=14), #axis.text.x = element_blank(),
        axis.text.y = element_text(size=10),axis.title.x = element_blank(), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_rect(fill = "transparent"), plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        legend.background = element_rect(fill = "transparent"), legend.box.background = element_rect(fill = "transparent"), strip.background.y = element_blank(),
        strip.text = element_text(size=14, face = "bold"), strip.background = element_blank(), legend.position = "none"); mel.w.outlier.number.min.max.10th.90th.plot



### concatenate ####
(era.w.outlier.number.min.max.10th.90th.plot | era.e.outlier.number.min.max.10th.90th.plot )/
  ( mel.w.outlier.number.min.max.10th.90th.plot | mel.e.outlier.number.min.max.10th.90th.plot)

cowplot::ggsave2(filename  = "git/2021-altitude-heliconius/plots/haplo.stats/era.mel.outlier.number.min.max.10th.90th.sims.deltapi.omega.ihs.sweed.altgwas.barplot.png", height = 5.4, width = 6,dpi=300)

############################################### 9. plot altitude  GWAS ###############################################
# shdr
setwd("git/2021-altitude-heliconius/")
era.e.shdr <- read.csv("data/sharing/era.shdr.para.east.csv"); head(era.e.shdr)
era.w.shdr <- read.csv("data/sharing/era.shdr.para.west.csv"); head(era.w.shdr)
mel.e.shdr <- read.csv("data/sharing/mel.shdr.para.east.csv"); head(mel.e.shdr)
mel.w.shdr <- read.csv("data/sharing/mel.shdr.para.west.csv"); head(mel.w.shdr)

era.stat.shdr.summ <-  read.csv("git/2021-altitude-heliconius/data/haplo.stats/era.stat.shdr.summ.csv")
mel.stat.shdr.summ <-  read.csv("git/2021-altitude-heliconius/data/haplo.stats/mel.stat.shdr.summ.csv")

unique(era.stat.shdr.summ$shdr.para.east.id) %in% unique(era.e.shdr$shdr.para.east.id)
head(era.stat.shdr.summ)
era.e.shdr$is.outlier.alt.gwas <- if_else(era.e.shdr$shdr.para.east.id %in% subset(era.stat.shdr.summ, log10.LRT.pval_median.max.90th.perc.sims =="yes")$shdr.para.east.id, "yes", "no" ); era.e.shdr$is.outlier.alt.gwas
mel.e.shdr$is.outlier.alt.gwas <- if_else(mel.e.shdr$shdr.para.east.id %in% subset(mel.stat.shdr.summ, log10.LRT.pval_median.max.90th.perc.sims =="yes")$shdr.para.east.id, "yes", "no" ); mel.e.shdr$is.outlier.alt.gwas

nrow(subset(era.stat.shdr.summ, log10.LRT.pval_median.max.90th.perc.sims =="yes"))/nrow(era.stat.shdr.summ) *100
nrow(subset(mel.stat.shdr.summ, log10.LRT.pval_median.max.90th.perc.sims =="yes"))/nrow(mel.stat.shdr.summ) *100

nrow(subset(era.stat.shdr.summ, log10.LRT.pval_median.max.90th.perc.sims =="yes"&is.allo.hdr=="no"))/nrow(subset(era.stat.shdr.summ,is.allo.hdr=="no")) *100
nrow(subset(mel.stat.shdr.summ, log10.LRT.pval_median.max.90th.perc.sims =="yes"&is.allo.hdr=="no"))/nrow(subset(mel.stat.shdr.summ,is.allo.hdr=="no")) *100





# west
era.stat.shdr.summ.w <-  read.csv("git/2021-altitude-heliconius/data/haplo.stats/west.hdr/era.stat.shdr.summ.csv")
mel.stat.shdr.summ.w <-  read.csv("git/2021-altitude-heliconius/data/haplo.stats/west.hdr/mel.stat.shdr.summ.csv")

head(era.stat.shdr.summ)
era.w.shdr$is.outlier.alt.gwas <- if_else(era.w.shdr$shdr.para.west.id %in% subset(era.stat.shdr.summ.w, log10.LRT.pval_median.max.90th.perc.sims =="yes")$shdr.para.west.id, "yes", "no" ); era.w.shdr$is.outlier.alt.gwas
mel.w.shdr$is.outlier.alt.gwas <- if_else(mel.w.shdr$shdr.para.west.id %in% subset(mel.stat.shdr.summ.w, log10.LRT.pval_median.max.90th.perc.sims =="yes")$shdr.para.west.id, "yes", "no" ); mel.w.shdr$is.outlier.alt.gwas

nrow(subset(era.stat.shdr.summ.w, log10.LRT.pval_median.max.90th.perc.sims =="yes" &is.allo.hdr=="no"))/nrow(subset(era.stat.shdr.summ.w, is.allo.hdr=="no")) *100
nrow(subset(mel.stat.shdr.summ.w, log10.LRT.pval_median.max.90th.perc.sims =="yes"&is.allo.hdr=="no"))/nrow(subset(mel.stat.shdr.summ.w, is.allo.hdr=="no")) *100



# altitude gwas
era.alt.gwas <- read.table("~/Dropbox/PhD/29.shape.paper/data/data.clean/era.n479/PC062_merged_Herato.PL.AD.HAPCUT2.inf0.05.corr.gwas.alt.admix.lrt0.gz.50snp.10snp.pos.ranks.short")
mel.alt.gwas <- read.table("~/Dropbox/PhD/29.shape.paper/data/data.clean/mel.187/run164_merged_Hmel.PL.AD.HAPCUT2.inf0.05.corr.gwas.alt.admix.doasso6.50snp.10snp.pos.ranks.short")
mel.alt.gwas$CHR <- as.numeric(as.character(substr(mel.alt.gwas$scaff, 6,7))); unique(mel.alt.gwas$CHR)
axisdf.era = era.alt.gwas %>% dplyr::group_by(CHR) %>% dplyr::summarize(center=( max(BP.wg) + min(BP.wg) ) / 2 ); axisdf.era
axisdf.era$CHR[21]<-"Z";axisdf.era$CHR[21]

axisdf.mel = mel.alt.gwas %>% dplyr::group_by(CHR) %>% dplyr::summarize(center=( max(BP.wg) + min(BP.wg) ) / 2 ); axisdf.mel
axisdf.mel$CHR[21]<-"Z";axisdf.mel$CHR[21]


dev.off()
names(era.alt.gwas)
ggplot(subset(era.alt.gwas, -log10(LRT.pval_median)>0), aes(x=BP.wg, y=-log10(LRT.pval_median))) +
  geom_rect(inherit.aes = F, data=subset(era.e.shdr), aes(xmin=start-100000, xmax=end+100000,ymin=-0.02,ymax=25, fill=is.allopatric), colour="transparent", alpha=.5 )+
  scale_fill_manual(values=c("#049E73", "#D65D00"))+
  ggnewscale::new_scale_fill()+
  geom_rect(inherit.aes = F, data=subset(era.w.shdr), aes(xmin=start-100000, xmax=end+100000, ymin=-0.02,ymax=25, fill=is.allopatric), colour="transparent", alpha=.5 )+
  scale_fill_manual(values=c("#0372B2", "#D65D00"))+
  #geom_rect(inherit.aes = F, data=subset(era.w.shdr&is.allopatric=="no"), aes(xmin=start, xmax=end, ymin=0,ymax=1), colour="transparent", fill="#0372B2", alpha=.3 )+
  #annotate("text",  x=Inf, y = Inf, label = i, vjust=1.5, hjust=1.1)  +
  #geom_label_repel(inherit.aes = F, data=subset(era.e.shdr), aes(x=(start+end)/2, y=23, label=substr(shdr.para.east.id,11,13 ), color=is.allopatric), min.segment.length = 0.05) +
  scale_colour_manual(values=c("#049E73", "#D65D00"))+
  #geom_label_repel(inherit.aes = F, data=subset(era.w.shdr&is.allopatric=="no"), aes(x=(start+end)/2, y=0.89, label=substr(shdr.para.west.id,11,13 )), color="#0372B2", min.segment.length = 0.05) +
  ggnewscale::new_scale_color()+
  geom_point(aes(y=-log10(LRT.pval_median),  color=as.factor(CHR)), size=.5) +
  scale_color_manual(values = rep(c("black", "grey36"), 22 )) +
  
  geom_segment(inherit.aes = F, data=(subset(era.e.shdr, is.outlier.alt.gwas=="yes"&is.allopatric=="no")), aes(x=(start+end)/2, y=24.98, xend=((start+end)/2)+.1, yend=22.82),
  size = .5, arrow = arrow(length = unit(.3, "inches"), type="closed"),  arrow.fill="#049E73", alpha=0.7) +
  geom_segment(inherit.aes = F, data=(subset(era.e.shdr, is.outlier.alt.gwas=="yes"&is.allopatric=="yes")), aes(x=(start+end)/2, y=24.98, xend=((start+end)/2)+.1, yend=22.82),
               size = .5, arrow = arrow(length = unit(.3, "inches"), type="closed"),  arrow.fill="#D65D00", alpha=0.7) +
    scale_x_continuous( expand = c(0, 0),label = axisdf.era$CHR, breaks= axisdf.era$center ) +
  scale_y_continuous( expand = c(0, 0), limits=c(-0.02,25.1)  ) +    # remove space between plot area and x axis
  ylab("-log10(median p-value) altitude GWAS")+
  theme_classic()+
  theme(legend.position="none")+ # Remove legend
  theme(axis.line.x = element_line(color="black", size = 0.5),axis.line.y = element_line(color="black", size = 0.5),axis.title.y = element_text(size=14),
        plot.margin = unit(c(0.5,0.2,0.2,0.2), "lines"),axis.text = element_text(size=10), #axis.text.x = element_blank(),
        axis.text.y = element_text(size=10),axis.title.x = element_blank(),       panel.background = element_rect(fill = "transparent"), plot.background = element_rect(fill = "transparent", color = NA),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.background = element_rect(fill = "transparent"), legend.box.background = element_rect(fill = "transparent")) 

cowplot::ggsave2(width = 36, height = 9, units = "cm", dpi = 300, filename = paste0("git/2021-altitude-heliconius/plots/haplo.stats/alt.gwas/era.alt.gwas.png"))


names(mel.alt.gwas)
ggplot(subset(mel.alt.gwas, -log10(LRT.pval_median)>0), aes(x=BP.wg, y=-log10(LRT.pval_median))) +
  geom_rect(inherit.aes = F, data=subset(mel.e.shdr), aes(xmin=start-10000, xmax=end+10000,ymin=-0.02,ymax=6.1, fill=is.allopatric), colour="transparent", alpha=.5 )+
  scale_fill_manual(values=c("#049E73", "#D65D00"))+
  ggnewscale::new_scale_fill()+
  geom_rect(inherit.aes = F, data=subset(mel.w.shdr), aes(xmin=start-10000, xmax=end+10000, ymin=-0.02,ymax=6.1, fill=is.allopatric), colour="transparent", alpha=.5 )+
  scale_fill_manual(values=c("#0372B2", "#D65D00"))+
  #geom_rect(inherit.aes = F, data=subset(mel.w.shdr&is.allopatric=="no"), aes(xmin=start, xmax=end, ymin=0,ymax=1), colour="transparent", fill="#0372B2", alpha=.3 )+
  #annotate("text",  x=Inf, y = Inf, label = i, vjust=1.5, hjust=1.1)  +
  #geom_label_repel(inherit.aes = F, data=subset(mel.e.shdr), aes(x=(start+end)/2, y=23, label=substr(shdr.para.east.id,11,13 ), color=is.allopatric), min.segment.length = 0.05) +
  scale_colour_manual(values=c("#049E73", "#D65D00"))+
  #geom_label_repel(inherit.aes = F, data=subset(mel.w.shdr&is.allopatric=="no"), aes(x=(start+end)/2, y=0.89, label=substr(shdr.para.west.id,11,13 )), color="#0372B2", min.segment.length = 0.05) +
  ggnewscale::new_scale_color()+
  geom_point(aes(y=-log10(LRT.pval_median),  color=as.factor(CHR)), size=.5) +
  scale_color_manual(values = rep(c("black", "grey36"), 22 )) +
  
  geom_segment(inherit.aes = F, data=(subset(mel.e.shdr, is.outlier.alt.gwas=="yes"&is.allopatric=="no")), aes(x=(start+end)/2, y=5.98, xend=((start+end)/2)+.1, yend=5.42),
               size = .5, arrow = arrow(length = unit(.3, "inches"), type="closed"),  arrow.fill="#049E73", alpha=0.7) +
  geom_segment(inherit.aes = F, data=(subset(mel.e.shdr, is.outlier.alt.gwas=="yes"&is.allopatric=="yes")), aes(x=(start+end)/2, y=5.98, xend=((start+end)/2)+.1, yend=5.42),
               size = .5, arrow = arrow(length = unit(.3, "inches"), type="closed"),  arrow.fill="#D65D00", alpha=0.7) +
  scale_x_continuous( expand = c(0, 0),label = axisdf.mel$CHR, breaks= axisdf.mel$center ) +
  scale_y_continuous( expand = c(0, 0), limits=c(-0.02,6.1)  ) +    # remove space between plot area and x axis
  ylab("-log10(median p-value) altitude GWAS")+
  theme_classic()+
  theme(legend.position="none")+ # Remove legend
  theme(axis.line.x = element_line(color="black", size = 0.5),axis.line.y = element_line(color="black", size = 0.5),axis.title.y = element_text(size=14),
        plot.margin = unit(c(0.5,0.2,0.2,0.2), "lines"),axis.text = element_text(size=10), #axis.text.x = element_blank(),
        axis.text.y = element_text(size=10),axis.title.x = element_blank(),       panel.background = element_rect(fill = "transparent"), plot.background = element_rect(fill = "transparent", color = NA),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.background = element_rect(fill = "transparent"), legend.box.background = element_rect(fill = "transparent")) 

cowplot::ggsave2(width = 36, height = 9, units = "cm", dpi = 300, filename = paste0("git/2021-altitude-heliconius/plots/haplo.stats/alt.gwas/mel.alt.gwas.png"))



############################################### INVERSIONS OVERLAPS, era mel #################################
setwd("/Users/gabrielamontejokovacevich/Dropbox (Cambridge University)/PhD/29.shape.paper/data/data.clean/selectionStatsHaplotaggingData/")
inv <-  read.csv("../../data.clean/inversions/meier.inversions.csv")
inv.era <- subset(inv, species=="erato")
inv.mel <- subset(inv, species=="melpomene")

# add whole genome offsets
inv.era$junction1.left.in.wg <- inv.era$junction1.left.in +ref.scaff.era[match(as.character(inv.era$Scaffold),as.character(ref.scaff.era$scaff)),"add"]
inv.era$junction2.right.in.wg <- inv.era$junction2.right.in +ref.scaff.era[match(as.character(inv.era$Scaffold),as.character(ref.scaff.era$scaff)),"add"]
inv.era$junction1.left.in.wg_junction2.right.in.wg_mean <- (inv.era$junction1.left.in.wg+inv.era$junction2.right.in.wg)/2

inv.mel$junction1.left.in.wg <- inv.mel$junction1.left.in +ref.scaff.mel[match(as.character(inv.mel$Scaffold),as.character(ref.scaff.mel$scaff)),"add"]
inv.mel$junction2.right.in.wg <- inv.mel$junction2.right.in +ref.scaff.mel[match(as.character(inv.mel$Scaffold),as.character(ref.scaff.mel$scaff)),"add"]
inv.mel$junction1.left.in.wg_junction2.right.in.wg_mean <- (inv.mel$junction1.left.in.wg+inv.mel$junction2.right.in.wg)/2

# check if the midpos, start or end of the inversion falls within any hdr
#function to check whether each of ints1 overlaps any of ints2
library(intervals)
overlaps_any <- function(ints1, ints2){
  overlaps <- interval_overlap(ints1,ints2)
  sapply(overlaps, length, simplify = T) != 0
}

## create hdr and inversion intervals 
era.e.shdr.int <- intervals::reduce(Intervals(as.matrix(era.e.shdr[,c("start","end")]), closed=c(TRUE,TRUE), type="R")) 
era.w.shdr.int <- intervals::reduce(Intervals(as.matrix(era.w.shdr[,c("start","end")]), closed=c(TRUE,TRUE), type="R")) 

inv.era.int <- intervals::reduce(Intervals(as.matrix(inv.era[,c("junction1.left.in.wg","junction2.right.in.wg")]), closed=c(TRUE,TRUE), type="R")) 


era.e.shdr$overlaps.with.inversion <- if_else(overlaps_any(era.e.shdr.int , inv.era.int)=="TRUE", "yes", "no")
era.stat.shdr.summ$overlaps.with.inversion <- if_else(overlaps_any(era.e.shdr.int , inv.era.int)=="TRUE", "yes", "no"); head(era.stat.shdr.summ)

era.w.shdr$overlaps.with.inversion <- if_else(overlaps_any(era.w.shdr.int , inv.era.int)=="TRUE", "yes", "no")

# save
era.e.shdr; write.csv(era.e.shdr, "git/2021-altitude-heliconius/data/sharing/era.shdr.para.east.csv", row.names = F)
era.stat.shdr.summ; write.csv(era.stat.shdr.summ, "git/2021-altitude-heliconius/data/haplo.stats/era.stat.shdr.summ.csv", row.names = F)
era.w.shdr; write.csv(era.w.shdr, "git/2021-altitude-heliconius/data/sharing/era.shdr.para.west.csv", row.names = F)

#### mel ######
## create hdr and inversion intervals 
mel.e.shdr.int <- intervals::reduce(Intervals(as.matrix(mel.e.shdr[,c("start","end")]), closed=c(TRUE,TRUE), type="R")) 
inv.mel.int <- intervals::reduce(Intervals(as.matrix(inv.mel[,c("junction1.left.in.wg","junction2.right.in.wg")]), closed=c(TRUE,TRUE), type="R")) 
mel.w.shdr.int <- intervals::reduce(Intervals(as.matrix(mel.w.shdr[,c("start","end")]), closed=c(TRUE,TRUE), type="R")) 

mel.e.shdr$overlaps.with.inversion <- if_else(overlaps_any(mel.e.shdr.int , inv.mel.int)=="TRUE", "yes", "no")
mel.stat.shdr.summ$overlaps.with.inversion <- if_else(overlaps_any(mel.e.shdr.int , inv.mel.int)=="TRUE", "yes", "no"); head(mel.stat.shdr.summ)
mel.w.shdr$overlaps.with.inversion <- if_else(overlaps_any(mel.w.shdr.int , inv.mel.int)=="TRUE", "yes", "no")

# save
mel.e.shdr; write.csv(mel.e.shdr, "git/2021-altitude-heliconius/data/sharing/mel.shdr.para.east.csv", row.names = F)
#mel.stat.shdr.summ; write.csv(mel.stat.shdr.summ, "git/2021-altitude-heliconius/data/haplo.stats/mel.stat.shdr.summ.csv", row.names = F)
mel.w.shdr; write.csv(mel.w.shdr, "git/2021-altitude-heliconius/data/sharing/mel.shdr.para.west.csv", row.names = F)





############################################### shape.roiS OVERLAPS, era mel #################################
shape.roi.era <-  read.csv("data/shape.roi/ROI.zoom.sections.spread.era.csv")
shape.roi.mel <-  read.csv("data/shape.roi/ROI.zoom.sections.spread.mel.csv")

shape.roi.era$size <- shape.roi.era$end - shape.roi.era$start
sum(shape.roi.era$size)/ era.genome_size *100

shape.roi.mel$size <- shape.roi.mel$end - shape.roi.mel$start
sum(shape.roi.mel$size)/ mel.genome_size *100

# check if the midpos, start or end of the shape.roi falls within any hdr
#function to check whether each of ints1 overlaps any of ints2
library(intervals)
overlaps_any <- function(ints1, ints2){
  overlaps <- interval_overlap(ints1,ints2)
  sapply(overlaps, length, simplify = T) != 0
}


##### era ######


## create hdr and shape.roi intervals 
era.e.shdr.outlier.int <- intervals::reduce(Intervals(as.matrix(era.e.shdr.outlier[,c("start","end")]), closed=c(TRUE,TRUE), type="R")) 
era.w.shdr.outlier.int <- intervals::reduce(Intervals(as.matrix(era.w.shdr.outlier[,c("start","end")]), closed=c(TRUE,TRUE), type="R")) 
head(shape.roi.era)
shape.roi.era.int <- intervals::reduce(Intervals(as.matrix(shape.roi.era[,c("start","end")]), closed=c(TRUE,TRUE), type="R")) 


era.e.shdr.outlier$overlaps.with.shape.roi <- if_else(overlaps_any(era.e.shdr.outlier.int , shape.roi.era.int)=="TRUE", "yes", "no")
era.w.shdr.outlier$overlaps.with.shape.roi <- if_else(overlaps_any(era.w.shdr.outlier.int , shape.roi.era.int)=="TRUE", "yes", "no")

era.e.shdr.outlier[era.e.shdr.outlier$overlaps.with.shape.roi=="yes",]$overlaps.with.shape.roi.ID <- shape.roi.era[overlaps_any(shape.roi.era.int, era.e.shdr.outlier.int ),]$zoom.no

era.e.shdr.outlier[overlaps_any(era.e.shdr.outlier.int, shape.roi.era.int ),]

shape.roi.era[overlaps_any(shape.roi.era.int , era.e.shdr.outlier.int ),]
shape.roi.era[overlaps_any(shape.roi.era.int , era.w.shdr.outlier.int ),]

shape.roi.era[overlaps_any(shape.roi.era.int , era.shdr.outlier.int ),]


# what number of ROIs overlapped
nrow(shape.roi.era[overlaps_any(shape.roi.era.int , era.e.shdr.outlier.int ),])
nrow(shape.roi.era[overlaps_any(shape.roi.era.int , era.w.shdr.outlier.int ),])



# get the roi ID

if_else(overlaps_any(era.e.shdr.outlier.int , shape.roi.era.int)=="TRUE", "yes", "no")

# save
# era.e.shdr.outlier; write.csv(era.e.shdr.outlier, "git/2021-altitude-heliconius/data/shdr.summ/shdr.era.east.outlier.df.csv", row.names = F)
# era.w.shdr.outlier; write.csv(era.w.shdr.outlier, "git/2021-altitude-heliconius/data/shdr.summ/shdr.era.west.outlier.df.csv", row.names = F)

##### mel ######


## create hdr and shape.roi intervals 
mel.e.shdr.outlier.int <- intervals::reduce(Intervals(as.matrix(mel.e.shdr.outlier[,c("start","end")]), closed=c(TRUE,TRUE), type="R")) 
mel.w.shdr.outlier.int <- intervals::reduce(Intervals(as.matrix(mel.w.shdr.outlier[,c("start","end")]), closed=c(TRUE,TRUE), type="R")) 
head(shape.roi.mel)
shape.roi.mel.int <- intervals::reduce(Intervals(as.matrix(shape.roi.mel[,c("start","end")]), closed=c(TRUE,TRUE), type="R")) 

mel.e.shdr.outlier$overlaps.with.shape.roi <- if_else(overlaps_any(mel.e.shdr.outlier.int , shape.roi.mel.int)=="TRUE", "yes", "no")
mel.w.shdr.outlier$overlaps.with.shape.roi <- if_else(overlaps_any(mel.w.shdr.outlier.int , shape.roi.mel.int)=="TRUE", "yes", "no")

shape.roi.mel[overlaps_any(shape.roi.mel.int , mel.e.shdr.outlier.int ),]
shape.roi.mel[overlaps_any(shape.roi.mel.int , mel.w.shdr.outlier.int ),]


# save
# mel.e.shdr.outlier; write.csv(mel.e.shdr.outlier, "git/2021-altitude-heliconius/data/shdr.summ/shdr.mel.east.outlier.df.csv", row.names = F)
# mel.w.shdr.outlier; write.csv(mel.w.shdr.outlier, "git/2021-altitude-heliconius/data/shdr.summ/shdr.mel.west.outlier.df.csv", row.names = F)





##### era simulations with all shdr in one dataset ######




era.shdr.outlier <- rbind(era.e.shdr.outlier[,c("start","end")], era.w.shdr.outlier[,c("start","end")])
era.shdr.outlier.int <- intervals::reduce(Intervals(as.matrix(era.shdr.outlier), closed=c(TRUE,TRUE), type="R")) 
nrow(era.shdr.outlier.int)

if_else(overlaps_any(era.shdr.outlier.int , shape.roi.era.int)=="TRUE", "yes", "no")

# get numbers for later
era.overlaps.intA <- overlaps_any(shape.roi.era.int , era.shdr.outlier.int )
era.overlaps.int.countA <- sum(era.overlaps.intA ); era.overlaps.int.countA



# random intervals
sim1k.shape.roi.era.int <- lapply(1:10000, function(x){rand_non_overlapping_intervals(era.genome_size, size(shape.roi.era.int))})

era.shdr.overlap.ran.ROI.summ <-list()

for (i in 1:10000) {
  # prep random intervals for action
  sim1k.shape.roi.era.int.tmp <- sim1k.shape.roi.era.int[[i]]
  
  # what number of ROIs overlapped
  era.shdr.overlap.ran.ROI.summ[[i]] <-nrow(sim1k.shape.roi.era.int[[i]][overlaps_any(sim1k.shape.roi.era.int[[i]] , era.shdr.outlier.int ),])

}

era.shdr.overlap.ran.ROI.summ.df <- as.data.frame(unlist(era.shdr.overlap.ran.ROI.summ))
names(era.shdr.overlap.ran.ROI.summ.df)<- c("no.roi.overlaps")

era.ran.roi.overlap.shdr <- ggplot(era.shdr.overlap.ran.ROI.summ.df, aes(x=no.roi.overlaps)) +
  geom_histogram(binwidth=1, colour="black", fill="white", position="dodge") +
  geom_vline(xintercept = 5, color="Red", size=2)+
  geom_vline(xintercept = quantile(era.shdr.overlap.ran.ROI.summ.df$no.roi.overlaps, 0.90), size=1, color="black", lty="dashed")+
  xlab("Number of wing shape ROI\noverlaps with SHDR")+
  scale_y_continuous(expand = c(0,0))+
  theme_classic()+
  theme(
    panel.background = element_rect(fill = "transparent"), # bg of the panel
    plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
    panel.grid.major = element_blank(), # get rid of major grid
    panel.grid.minor = element_blank(), # get rid of minor grid
    legend.background = element_rect(fill = "transparent"), # get rid of legend bg
  );era.ran.roi.overlap.shdr




##### mel simulations with all shdr in one dataset ######

mel.shdr.outlier <- rbind(mel.e.shdr.outlier[,c("start","end")], mel.w.shdr.outlier[,c("start","end")])
mel.shdr.outlier.int <- intervals::reduce(Intervals(as.matrix(mel.shdr.outlier), closed=c(TRUE,TRUE), type="R")) 
nrow(mel.shdr.outlier.int)

# random intervals
sim1k.shape.roi.mel.int <- lapply(1:10000, function(x){rand_non_overlapping_intervals(mel.genome_size, size(shape.roi.mel.int))})

mel.shdr.overlap.ran.ROI.summ <-list()

for (i in 1:10000) {
  # prep random intervals for action
  sim1k.shape.roi.mel.int.tmp <- sim1k.shape.roi.mel.int[[i]]
  
  # what number of ROIs overlapped
  mel.shdr.overlap.ran.ROI.summ[[i]] <-nrow(sim1k.shape.roi.mel.int[[i]][overlaps_any(sim1k.shape.roi.mel.int[[i]] , mel.shdr.outlier.int ),])
  
}

mel.shdr.overlap.ran.ROI.summ.df <- as.data.frame(unlist(mel.shdr.overlap.ran.ROI.summ))
names(mel.shdr.overlap.ran.ROI.summ.df)<- c("no.roi.overlaps")

mel.ran.roi.overlap.shdr <- ggplot(mel.shdr.overlap.ran.ROI.summ.df, aes(x=no.roi.overlaps)) +
  geom_histogram(binwidth=1, colour="black", fill="white", position="dodge") +
  geom_density()+
  geom_vline(xintercept = 2, color="Red", size=2)+
  geom_vline(xintercept = quantile(mel.shdr.overlap.ran.ROI.summ.df$no.roi.overlaps, 0.9), size=1, color="black", lty="dashed")+
  xlab("Number of wing shape ROI\noverlaps with SHDR")+
  scale_y_continuous(expand = c(0,0))+
  theme_classic()+
  theme(
    panel.background = element_rect(fill = "transparent"), # bg of the panel
    plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
    panel.grid.major = element_blank(), # get rid of major grid
    panel.grid.minor = element_blank(), # get rid of minor grid
    legend.background = element_rect(fill = "transparent"), # get rid of legend bg
  );mel.ran.roi.overlap.shdr


plot_grid(era.ran.roi.overlap.shdr, mel.ran.roi.overlap.shdr, labels = c("A", "B"))
ggsave("plots/shape.roi.overlaps/shape.roi.overlap.shdr.10ksims.png", width = 5, height = 2)

##### simulations ######
# the shape ROIs tend to be bigger than the SHDR, 
# so permutate position of ROI intervals and check overlaps with real SHDRs

# random intervals
sim1k.shape.roi.era.int <- lapply(1:10000, function(x){rand_non_overlapping_intervals(era.genome_size, size(shape.roi.era.int))})

era.e.shdr.overlap.ran.ROI.summ <-list()
era.w.shdr.overlap.ran.ROI.summ <-list()

for (i in 1:10000) {
  # prep random intervals for action
  sim1k.shape.roi.era.int.tmp <- sim1k.shape.roi.era.int[[i]]

  # what number of ROIs overlapped
  era.e.shdr.overlap.ran.ROI.summ[[i]] <-nrow(sim1k.shape.roi.era.int[[i]][overlaps_any(sim1k.shape.roi.era.int[[i]] , era.e.shdr.outlier.int ),])
  era.w.shdr.overlap.ran.ROI.summ[[i]] <-nrow(sim1k.shape.roi.era.int[[i]][overlaps_any(sim1k.shape.roi.era.int[[i]] , era.w.shdr.outlier.int ),])
  
  }

era.e.shdr.overlap.ran.ROI.summ.df <- as.data.frame(unlist(era.e.shdr.overlap.ran.ROI.summ))
era.w.shdr.overlap.ran.ROI.summ.df <- as.data.frame(unlist(era.w.shdr.overlap.ran.ROI.summ))

names(era.e.shdr.overlap.ran.ROI.summ.df)<- c("no.roi.overlaps")
names(era.w.shdr.overlap.ran.ROI.summ.df)<- c("no.roi.overlaps")

era.e.ran.roi.overlap.shdr <- ggplot(era.e.shdr.overlap.ran.ROI.summ.df, aes(x=no.roi.overlaps)) +
  geom_histogram(binwidth=1, colour="black", fill="white", position="dodge") +
  geom_vline(xintercept = 3, color="Red", size=2)+
  geom_vline(xintercept = quantile(era.e.shdr.overlap.ran.ROI.summ.df$no.roi.overlaps, 0.90), size=1, color="blue", lty="dashed");era.e.ran.roi.overlap.shdr

era.w.ran.roi.overlap.shdr <-ggplot(era.w.shdr.overlap.ran.ROI.summ.df, aes(x=no.roi.overlaps)) +
  geom_histogram(binwidth=1, colour="black", fill="white", position="dodge") +
  geom_vline(xintercept = 5, color="Red", size=2)+
  geom_vline(xintercept = quantile(era.w.shdr.overlap.ran.ROI.summ.df$no.roi.overlaps, 0.90), size=1, color="blue", lty="dashed");era.w.ran.roi.overlap.shdr

plot_grid(era.e.ran.roi.overlap.shdr, era.w.ran.roi.overlap.shdr, labels = c("era.e.shdr", "era.w.shdr"))

##### mel ####
# random intervals
sim1k.shape.roi.mel.int <- lapply(1:10000, function(x){rand_non_overlapping_intervals(mel.genome_size, size(shape.roi.mel.int))})

mel.e.shdr.overlap.ran.ROI.summ <-list()
mel.w.shdr.overlap.ran.ROI.summ <-list()

for (i in 1:10000) {
  # prep random intervals for action
  sim1k.shape.roi.mel.int.tmp <- sim1k.shape.roi.mel.int[[i]]
  
  # what number of ROIs overlapped
  mel.e.shdr.overlap.ran.ROI.summ[[i]] <-nrow(sim1k.shape.roi.mel.int[[i]][overlaps_any(sim1k.shape.roi.mel.int[[i]] , mel.e.shdr.outlier.int ),])
  mel.w.shdr.overlap.ran.ROI.summ[[i]] <-nrow(sim1k.shape.roi.mel.int[[i]][overlaps_any(sim1k.shape.roi.mel.int[[i]] , mel.w.shdr.outlier.int ),])
  
  
}

mel.e.shdr.overlap.ran.ROI.summ.df <- as.data.frame(unlist(mel.e.shdr.overlap.ran.ROI.summ))
mel.w.shdr.overlap.ran.ROI.summ.df <- as.data.frame(unlist(mel.w.shdr.overlap.ran.ROI.summ))

names(mel.e.shdr.overlap.ran.ROI.summ.df)<- c("no.roi.overlaps")
names(mel.w.shdr.overlap.ran.ROI.summ.df)<- c("no.roi.overlaps")

mel.e.ran.roi.overlap.shdr <- ggplot(mel.e.shdr.overlap.ran.ROI.summ.df, aes(x=no.roi.overlaps)) +
  geom_histogram(binwidth=1, colour="black", fill="white", position="dodge") +
  geom_vline(xintercept = 3, color="Red", size=2)+
  geom_vline(xintercept = quantile(mel.e.shdr.overlap.ran.ROI.summ.df$no.roi.overlaps, 0.90), size=1, color="blue", lty="dashed");mel.e.ran.roi.overlap.shdr

mel.w.ran.roi.overlap.shdr <-ggplot(mel.w.shdr.overlap.ran.ROI.summ.df, aes(x=no.roi.overlaps)) +
  geom_histogram(binwidth=1, colour="black", fill="white", position="dodge") +
  geom_vline(xintercept = 5, color="Red", size=2)+
  geom_vline(xintercept = quantile(mel.w.shdr.overlap.ran.ROI.summ.df$no.roi.overlaps, 0.90), size=1, color="blue", lty="dashed");mel.w.ran.roi.overlap.shdr

plot_grid(mel.e.ran.roi.overlap.shdr, mel.w.ran.roi.overlap.shdr, labels = c("mel.e.shdr", "mel.w.shdr"))


######### jackknife ############
### random jacknife blocks across genome
genome_size <- tail(scafEnds.era,1)
#genome_size <- genome_size-bp.removed.era.mean 
block_size <- 1000000
block_starts <- seq(1,genome_size, block_size)
n_blocks <- length(block_starts)
block_Intervals <- Intervals(as.matrix(cbind(block_starts,block_starts+block_size-1)), closed = c(TRUE,TRUE), type="R")

era.overlaps.int.countA #5

## across sides
#jackknifed trio intervals (remove 1 blcok in each case)
era.int.jackknifed.shape.roi <- lapply(1:n_blocks, function(x){interval_difference(shape.roi.era.int, block_Intervals[x,])})

era.overlap.counts.jackA <- unlist(lapply(1:n_blocks, function(x){sum(overlaps_any(era.int.jackknifed.shape.roi[[x]],era.shdr.outlier.int ))}))
era.overlaps.int.prop.jackknifedA <-  (era.overlap.counts.jackA / (nrow(era.e.overlap.union)))*100 

era.overlap.counts.pseudoA <- era.overlaps.int.countA*n_blocks - era.overlap.counts.jackA*(n_blocks-1)
era.overlap.counts.sdA <- sd(era.overlap.counts.pseudoA)
era.overlap.counts.pseudo.std_errA <- sqrt(var(era.overlap.counts.pseudoA)/length(era.overlap.counts.pseudoA))


, 
