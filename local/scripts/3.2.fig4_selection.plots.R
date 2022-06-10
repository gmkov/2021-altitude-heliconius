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
library(ggpubr)
options(scipen = 999)
library(intervals)
library(ggtext)
library(zoo)
library(patchwork)
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

match.by.range.1k <- function(vec) {
  x <- sim1k.era.df
  if(length(.x <- which(vec >= x$starts & vec <= x$ends ))) .x else NA }


#################### 0. data prep #####
setwd("git/2021-altitude-heliconius/data/pbs.out/")
pop.info <- read.csv("local/data/local/data/local/data/02.info/pop.short.info.csv")

# ## all windows
# era.all.pbs <- read.csv("local/data/pbs.out/era.all.pbs.recomb.csv")
# mel.all.pbs <- read.csv("local/data/pbs.out/mel.all.pbs.recomb.csv")
# 
####### shdr with oulier info from sims #######
shdr.era.east.outlier.df <- read.csv("local/data/shdr.summ/shdr.era.east.outlier.df.csv"); head(shdr.era.east.outlier.df)
shdr.era.west.outlier.df <- read.csv("local/data/shdr.summ/shdr.era.west.outlier.df.csv"); head(shdr.era.west.outlier.df)
shdr.era.outlier.df.summ <- read.csv( "local/data/shdr.summ/shdr.era.outlier.df.summ.csv")
shdr.era.outlier.df.summ$side <- factor(shdr.era.outlier.df.summ$side, levels=c("West", "East"))
shdr.era.outlier.df.summ$is.allo.hdr <- factor(shdr.era.outlier.df.summ$is.allo.hdr, levels=c("Parapatric", "Allopatric")); shdr.era.outlier.df.summ$is.allo.hdr

shdr.mel.east.outlier.df <- read.csv("local/data/shdr.summ/shdr.mel.east.outlier.df.csv"); head(shdr.mel.east.outlier.df)
shdr.mel.west.outlier.df <- read.csv("local/data/shdr.summ/shdr.mel.west.outlier.df.csv"); head(shdr.mel.west.outlier.df)
shdr.mel.outlier.df.summ <- read.csv( "local/data/shdr.summ/shdr.mel.outlier.df.summ.csv")
shdr.mel.outlier.df.summ$side <- factor(shdr.mel.outlier.df.summ$side, levels=c("West", "East"))
shdr.mel.outlier.df.summ$is.allo.hdr <- factor(shdr.mel.outlier.df.summ$is.allo.hdr, levels=c("Parapatric", "Allopatric")); shdr.mel.outlier.df.summ$is.allo.hdr

names(shdr.era.east.outlier.df)
########## max erato ######
## parapatric
clyEhighW.max <- nrow(subset(shdr.era.west.outlier.df, is.allo.hdr=="no" & (fdm.max.more.90th.perc.sims_fd_eraEvlowW_eraEhighW_clyEhighW_eleHW=="yes" | fdm.max.more.90th.perc.sims_fd_eraCvlowW_eraChighW_clyEhighW_eleHW=="yes")))/ nrow(subset(shdr.era.west.outlier.df, is.allo.hdr=="no" )) *100
nrow(subset(shdr.era.west.outlier.df, is.allo.hdr=="no" & (fdm.max.more.90th.perc.sims_fd_eraEvlowW_eraEhighW_clyEhighW_eleHW=="yes"  & fdm.max.more.90th.perc.sims_fd_eraCvlowW_eraChighW_clyEhighW_eleHW=="yes")))/ nrow(subset(shdr.era.west.outlier.df, is.allo.hdr=="no" )) *100

nrow(subset(shdr.era.west.outlier.df, is.allo.hdr=="yes" & (fdm.max.more.90th.perc.sims_fd_eraEvlowW_eraEhighW_clyEhighW_eleHW=="yes"  | fdm.max.more.90th.perc.sims_fd_eraCvlowW_eraChighW_clyEhighW_eleHW=="yes")))/ nrow(subset(shdr.era.west.outlier.df, is.allo.hdr=="yes" )) *100
nrow(subset(shdr.era.west.outlier.df, is.allo.hdr=="yes" & (fdm.max.more.90th.perc.sims_fd_eraEvlowW_eraEhighW_clyEhighW_eleHW=="yes"  & fdm.max.more.90th.perc.sims_fd_eraCvlowW_eraChighW_clyEhighW_eleHW=="yes")))/ nrow(subset(shdr.era.west.outlier.df, is.allo.hdr=="yes" )) *100


nrow(subset(shdr.era.east.outlier.df, is.allo.hdr=="no" & (fdm.max.more.90th.perc.sims_fd_eraEvlowE_eraEhighE_clyEhighE_eleHW=="yes" | fdm.max.more.90th.perc.sims_fd_eraCvlowE_eraChighE_clyEhighE_eleHW=="yes")))/ nrow(subset(shdr.era.east.outlier.df, is.allo.hdr=="no" )) *100
nrow(subset(shdr.era.east.outlier.df, is.allo.hdr=="no" & (fdm.max.more.90th.perc.sims_fd_eraEvlowE_eraEhighE_clyEhighE_eleHW=="yes"  & fdm.max.more.90th.perc.sims_fd_eraCvlowE_eraChighE_clyEhighE_eleHW=="yes")))/ nrow(subset(shdr.era.east.outlier.df, is.allo.hdr=="no" )) *100

nrow(subset(shdr.era.east.outlier.df, is.allo.hdr=="no" & (fdm.max.more.90th.perc.sims_fd_eraEvlowE_eraEhighE_telhighE_eleHW=="yes" | fdm.max.more.90th.perc.sims_fd_eraCvlowE_eraChighE_telhighE_eleHW=="yes")))/ nrow(subset(shdr.era.east.outlier.df, is.allo.hdr=="no" )) *100
nrow(subset(shdr.era.east.outlier.df, is.allo.hdr=="no" & (fdm.max.more.90th.perc.sims_fd_eraEvlowE_eraEhighE_telhighE_eleHW=="yes"  & fdm.max.more.90th.perc.sims_fd_eraCvlowE_eraChighE_telhighE_eleHW=="yes")))/ nrow(subset(shdr.era.east.outlier.df, is.allo.hdr=="no" )) *100

nrow(subset(shdr.era.east.outlier.df, is.allo.hdr=="no" & (fdm.max.more.90th.perc.sims_fd_eraEvlowE_eraEhighE_himHE_eleHW=="yes" | fdm.max.more.90th.perc.sims_fd_eraCvlowE_eraChighE_himHE_eleHW=="yes")))/ nrow(subset(shdr.era.east.outlier.df, is.allo.hdr=="no" )) *100
nrow(subset(shdr.era.east.outlier.df, is.allo.hdr=="no" & (fdm.max.more.90th.perc.sims_fd_eraEvlowE_eraEhighE_himHE_eleHW=="yes"  & fdm.max.more.90th.perc.sims_fd_eraCvlowE_eraChighE_himHE_eleHW=="yes")))/ nrow(subset(shdr.era.east.outlier.df, is.allo.hdr=="no" )) *100

## allopatric
nrow(subset(shdr.era.east.outlier.df, is.allo.hdr=="yes" & (fdm.max.more.90th.perc.sims_fd_eraEvlowE_eraEhighE_clyEhighE_eleHW=="yes" | fdm.max.more.90th.perc.sims_fd_eraCvlowE_eraChighE_clyEhighE_eleHW=="yes")))/ nrow(subset(shdr.era.east.outlier.df, is.allo.hdr=="yes" )) *100
nrow(subset(shdr.era.east.outlier.df, is.allo.hdr=="yes" & (fdm.max.more.90th.perc.sims_fd_eraEvlowE_eraEhighE_clyEhighE_eleHW=="yes"  & fdm.max.more.90th.perc.sims_fd_eraCvlowE_eraChighE_clyEhighE_eleHW=="yes")))/ nrow(subset(shdr.era.east.outlier.df, is.allo.hdr=="yes" )) *100

nrow(subset(shdr.era.east.outlier.df, is.allo.hdr=="yes" & (fdm.max.more.90th.perc.sims_fd_eraEvlowE_eraEhighE_telhighE_eleHW=="yes" | fdm.max.more.90th.perc.sims_fd_eraCvlowE_eraChighE_telhighE_eleHW=="yes")))/ nrow(subset(shdr.era.east.outlier.df, is.allo.hdr=="yes" )) *100
nrow(subset(shdr.era.east.outlier.df, is.allo.hdr=="yes" & (fdm.max.more.90th.perc.sims_fd_eraEvlowE_eraEhighE_telhighE_eleHW=="yes"  & fdm.max.more.90th.perc.sims_fd_eraCvlowE_eraChighE_telhighE_eleHW=="yes")))/ nrow(subset(shdr.era.east.outlier.df, is.allo.hdr=="yes" )) *100

nrow(subset(shdr.era.east.outlier.df, is.allo.hdr=="yes" & (fdm.max.more.90th.perc.sims_fd_eraEvlowE_eraEhighE_himHE_eleHW=="yes" | fdm.max.more.90th.perc.sims_fd_eraCvlowE_eraChighE_himHE_eleHW=="yes")))/ nrow(subset(shdr.era.east.outlier.df, is.allo.hdr=="yes" )) *100
nrow(subset(shdr.era.east.outlier.df, is.allo.hdr=="yes" & (fdm.max.more.90th.perc.sims_fd_eraEvlowE_eraEhighE_himHE_eleHW=="yes"  & fdm.max.more.90th.perc.sims_fd_eraCvlowE_eraChighE_himHE_eleHW=="yes")))/ nrow(subset(shdr.era.east.outlier.df, is.allo.hdr=="yes" )) *100

########## min erato ######
## parapatric
clyEhighW.para.min <- nrow(subset(shdr.era.west.outlier.df, is.allo.hdr=="no" & (fdm.min.less.10th.perc.sims_fd_eraEvlowW_eraEhighW_clyEhighW_eleHW=="yes" | fdm.min.less.10th.perc.sims_fd_eraCvlowW_eraChighW_clyEhighW_eleHW=="yes")))/ nrow(subset(shdr.era.west.outlier.df, is.allo.hdr=="no" )) *100; clyEhighW.para.min
clyEhighW.allo.min <-nrow(subset(shdr.era.west.outlier.df, is.allo.hdr=="no" & (fdm.min.less.10th.perc.sims_fd_eraEvlowW_eraEhighW_clyEhighW_eleHW=="yes"  & fdm.min.less.10th.perc.sims_fd_eraCvlowW_eraChighW_clyEhighW_eleHW=="yes")))/ nrow(subset(shdr.era.west.outlier.df, is.allo.hdr=="no" )) *100; clyEhighW.allo.min

nrow(subset(shdr.era.west.outlier.df, is.allo.hdr=="yes" & (fdm.min.less.10th.perc.sims_fd_eraEvlowW_eraEhighW_clyEhighW_eleHW=="yes"  | fdm.min.less.10th.perc.sims_fd_eraCvlowW_eraChighW_clyEhighW_eleHW=="yes")))/ nrow(subset(shdr.era.west.outlier.df, is.allo.hdr=="yes" )) *100
nrow(subset(shdr.era.west.outlier.df, is.allo.hdr=="yes" & (fdm.min.less.10th.perc.sims_fd_eraEvlowW_eraEhighW_clyEhighW_eleHW=="yes"  & fdm.min.less.10th.perc.sims_fd_eraCvlowW_eraChighW_clyEhighW_eleHW=="yes")))/ nrow(subset(shdr.era.west.outlier.df, is.allo.hdr=="yes" )) *100

nrow(subset(shdr.era.east.outlier.df, is.allo.hdr=="no" & (fdm.min.less.10th.perc.sims_fd_eraEvlowE_eraEhighE_clyEhighE_eleHW=="yes" | fdm.min.less.10th.perc.sims_fd_eraCvlowE_eraChighE_clyEhighE_eleHW=="yes")))/ nrow(subset(shdr.era.east.outlier.df, is.allo.hdr=="no" )) *100
nrow(subset(shdr.era.east.outlier.df, is.allo.hdr=="no" & (fdm.min.less.10th.perc.sims_fd_eraEvlowE_eraEhighE_clyEhighE_eleHW=="yes"  & fdm.min.less.10th.perc.sims_fd_eraCvlowE_eraChighE_clyEhighE_eleHW=="yes")))/ nrow(subset(shdr.era.east.outlier.df, is.allo.hdr=="no" )) *100

nrow(subset(shdr.era.east.outlier.df, is.allo.hdr=="no" & (fdm.min.less.10th.perc.sims_fd_eraEvlowE_eraEhighE_telhighE_eleHW=="yes" | fdm.min.less.10th.perc.sims_fd_eraCvlowE_eraChighE_telhighE_eleHW=="yes")))/ nrow(subset(shdr.era.east.outlier.df, is.allo.hdr=="no" )) *100
nrow(subset(shdr.era.east.outlier.df, is.allo.hdr=="no" & (fdm.min.less.10th.perc.sims_fd_eraEvlowE_eraEhighE_telhighE_eleHW=="yes"  & fdm.min.less.10th.perc.sims_fd_eraCvlowE_eraChighE_telhighE_eleHW=="yes")))/ nrow(subset(shdr.era.east.outlier.df, is.allo.hdr=="no" )) *100

nrow(subset(shdr.era.east.outlier.df, is.allo.hdr=="no" & (fdm.min.less.10th.perc.sims_fd_eraEvlowE_eraEhighE_himHE_eleHW=="yes" | fdm.min.less.10th.perc.sims_fd_eraCvlowE_eraChighE_himHE_eleHW=="yes")))/ nrow(subset(shdr.era.east.outlier.df, is.allo.hdr=="no" )) *100
nrow(subset(shdr.era.east.outlier.df, is.allo.hdr=="no" & (fdm.min.less.10th.perc.sims_fd_eraEvlowE_eraEhighE_himHE_eleHW=="yes"  & fdm.min.less.10th.perc.sims_fd_eraCvlowE_eraChighE_himHE_eleHW=="yes")))/ nrow(subset(shdr.era.east.outlier.df, is.allo.hdr=="no" )) *100

## allopatric
nrow(subset(shdr.era.east.outlier.df, is.allo.hdr=="yes" & (fdm.min.less.10th.perc.sims_fd_eraEvlowE_eraEhighE_clyEhighE_eleHW=="yes" | fdm.min.less.10th.perc.sims_fd_eraCvlowE_eraChighE_clyEhighE_eleHW=="yes")))/ nrow(subset(shdr.era.east.outlier.df, is.allo.hdr=="yes" )) *100
nrow(subset(shdr.era.east.outlier.df, is.allo.hdr=="yes" & (fdm.min.less.10th.perc.sims_fd_eraEvlowE_eraEhighE_clyEhighE_eleHW=="yes"  & fdm.min.less.10th.perc.sims_fd_eraCvlowE_eraChighE_clyEhighE_eleHW=="yes")))/ nrow(subset(shdr.era.east.outlier.df, is.allo.hdr=="yes" )) *100

nrow(subset(shdr.era.east.outlier.df, is.allo.hdr=="yes" & (fdm.min.less.10th.perc.sims_fd_eraEvlowE_eraEhighE_telhighE_eleHW=="yes" | fdm.min.less.10th.perc.sims_fd_eraCvlowE_eraChighE_telhighE_eleHW=="yes")))/ nrow(subset(shdr.era.east.outlier.df, is.allo.hdr=="yes" )) *100
nrow(subset(shdr.era.east.outlier.df, is.allo.hdr=="yes" & (fdm.min.less.10th.perc.sims_fd_eraEvlowE_eraEhighE_telhighE_eleHW=="yes"  & fdm.min.less.10th.perc.sims_fd_eraCvlowE_eraChighE_telhighE_eleHW=="yes")))/ nrow(subset(shdr.era.east.outlier.df, is.allo.hdr=="yes" )) *100

nrow(subset(shdr.era.east.outlier.df, is.allo.hdr=="yes" & (fdm.min.less.10th.perc.sims_fd_eraEvlowE_eraEhighE_himHE_eleHW=="yes" | fdm.min.less.10th.perc.sims_fd_eraCvlowE_eraChighE_himHE_eleHW=="yes")))/ nrow(subset(shdr.era.east.outlier.df, is.allo.hdr=="yes" )) *100
nrow(subset(shdr.era.east.outlier.df, is.allo.hdr=="yes" & (fdm.min.less.10th.perc.sims_fd_eraEvlowE_eraEhighE_himHE_eleHW=="yes"  & fdm.min.less.10th.perc.sims_fd_eraCvlowE_eraChighE_himHE_eleHW=="yes")))/ nrow(subset(shdr.era.east.outlier.df, is.allo.hdr=="yes" )) *100


## parapatric
names(shdr.mel.east.outlier.df)

nrow(subset(shdr.mel.west.outlier.df, is.allo.hdr=="no" & (fdm.max.more.90th.perc.sims_fd_melElowW_melEhighW_cydEhighW_hecale=="yes" | fdm.max.more.90th.perc.sims_fd_melCvlowW_melChighW_cydEhighW_hecale=="yes")))/ nrow(subset(shdr.mel.west.outlier.df, is.allo.hdr=="no" )) *100
nrow(subset(shdr.mel.west.outlier.df, is.allo.hdr=="no" & (fdm.max.more.90th.perc.sims_fd_melElowW_melEhighW_cydEhighW_hecale=="yes"  & fdm.max.more.90th.perc.sims_fd_melCvlowW_melChighW_cydEhighW_hecale=="yes")))/ nrow(subset(shdr.mel.west.outlier.df, is.allo.hdr=="no" )) *100

nrow(subset(shdr.mel.east.outlier.df, is.allo.hdr=="no" & (fdm.max.more.90th.perc.sims_fd_melEvlowE_melEhighE_timEhighE_hecale=="yes" | fdm.max.more.90th.perc.sims_fd_melCvlowE_melChighE_timEhighE_hecale=="yes")))/ nrow(subset(shdr.mel.east.outlier.df, is.allo.hdr=="no" )) *100
nrow(subset(shdr.mel.east.outlier.df, is.allo.hdr=="no" & (fdm.max.more.90th.perc.sims_fd_melEvlowE_melEhighE_timEhighE_hecale=="yes"  & fdm.max.more.90th.perc.sims_fd_melCvlowE_melChighE_timEhighE_hecale=="yes")))/ nrow(subset(shdr.mel.east.outlier.df, is.allo.hdr=="no" )) *100

## allopatric
nrow(subset(shdr.mel.west.outlier.df, is.allo.hdr=="yes" & (fdm.max.more.90th.perc.sims_fd_melElowW_melEhighW_cydEhighW_hecale=="yes" | fdm.max.more.90th.perc.sims_fd_melCvlowW_melChighW_cydEhighW_hecale=="yes")))/ nrow(subset(shdr.mel.west.outlier.df, is.allo.hdr=="yes" )) *100
nrow(subset(shdr.mel.west.outlier.df, is.allo.hdr=="yes" & (fdm.max.more.90th.perc.sims_fd_melElowW_melEhighW_cydEhighW_hecale=="yes"  & fdm.max.more.90th.perc.sims_fd_melCvlowW_melChighW_cydEhighW_hecale=="yes")))/ nrow(subset(shdr.mel.west.outlier.df, is.allo.hdr=="yes" )) *100

nrow(subset(shdr.mel.east.outlier.df, is.allo.hdr=="yes" & (fdm.max.more.90th.perc.sims_fd_melEvlowE_melEhighE_timEhighE_hecale=="yes" | fdm.max.more.90th.perc.sims_fd_melCvlowE_melChighE_timEhighE_hecale=="yes")))/ nrow(subset(shdr.mel.east.outlier.df, is.allo.hdr=="yes" )) *100
nrow(subset(shdr.mel.east.outlier.df, is.allo.hdr=="yes" & (fdm.max.more.90th.perc.sims_fd_melEvlowE_melEhighE_timEhighE_hecale=="yes"  & fdm.max.more.90th.perc.sims_fd_melCvlowE_melChighE_timEhighE_hecale=="yes")))/ nrow(subset(shdr.mel.east.outlier.df, is.allo.hdr=="yes" )) *100


nrow(subset(shdr.era.west.outlier.df, is.allo.hdr=="no" & (fdm.max.more.90th.perc.sims_fd_eraEvlowW_eraEhighW_clyEhighW_eleHW=="yes" | fdm.max.more.90th.perc.sims_fd_eraCvlowW_eraChighW_clyEhighW_eleHW=="yes" |
                                                             fdm.max.more.90th.perc.sims_fd_eraEvlowE_eraEhighE_clyEhighE_eleHW=="yes" | fdm.max.more.90th.perc.sims_fd_eraCvlowE_eraChighE_clyEhighE_eleHW=="yes"|
                                                             )))

####### shdr para, allo summs #######
era.shdr.para.east <- read.csv("local/data/sharing/era.shdr.para.east.csv")
era.shdr.para.west <- read.csv("local/data/sharing/era.shdr.para.west.csv")
era.shdr.allo.all  <- read.csv("local/data/sharing/era.shdr.allo.all.csv")

era.shdr.para.east <- read.csv("local/data/sharing/era.shdr.para.east.csv")
era.shdr.para.west <- read.csv("local/data/sharing/era.shdr.para.west.csv")
era.shdr.allo.all  <- read.csv("local/data/sharing/era.shdr.allo.all.csv")

# what para shdr make up the allo shdr
era.shdr.allo.summ <- read.csv("local/data/sharing/era.shdr.allo.summ.csv")
mel.shdr.allo.summ <- read.csv("local/data/sharing/mel.shdr.allo.summ.csv")

####### load genomes, for simulations #######
ref.scaff.era <- read.table("local/data/local/data/local/data/local/data/0.ref/Heliconius_erato_demophoon_v1_-_scaffolds.fa.fai", row.names = NULL)
ref.scaff.mel <- read.table("local/data/local/data/local/data/local/data/0.ref/Hmel2.5.scaffolds.fa.fai", row.names = NULL)
scafEnds.era <- cumsum(ref.scaff.era[,2]); offset.era <- scafEnds.era - ref.scaff.era[,2]
scafEnds.mel <- cumsum(ref.scaff.mel[,2]); offset.mel <- scafEnds.mel - ref.scaff.mel[,2]
era.genome_size <- tail(scafEnds.era,1)
mel.genome_size <- tail(scafEnds.mel,1)

####### load genes gff positions for plotting #######


library(intervals)
ref.scaff.era <- read.table("local/data/local/data/local/data/local/data/0.ref/Heliconius_erato_demophoon_v1_-_scaffolds.fa.fai", row.names = NULL)
ref.scaff.mel <- read.table("local/data/local/data/local/data/local/data/0.ref/Hmel2.5.scaffolds.fa.fai", row.names = NULL)

#get whole genome scaffold offsets from gff
era.scafEnds <- cumsum(ref.scaff.era[,2]); era.offset <- era.scafEnds - ref.scaff.era[,2]; names(era.offset) <- ref.scaff.era[,1]
mel.scafEnds <- cumsum(ref.scaff.mel[,2]); mel.offset <- mel.scafEnds - ref.scaff.mel[,2]; names(mel.offset) <- ref.scaff.mel[,1]

# load gff, get genes, add offsets, get intervals
era.gff <- read.table(pipe("cut -f 1,3,4,5,9 local/data/local/data/local/data/local/data/0.ref/heliconius_erato_demophoon_v1_core_32_85_1.gff"),
                      sep = "\t", as.is = T, col.names = c("scaffold", "feature", "start", "end", "info")); head(era.gff)
era.gff_gene <- subset(era.gff, feature == "gene"); head(era.gff_gene)
era.gff_gene$gene_name <- str_replace(str_remove(sub("\\;.*", "", era.gff_gene$info), "ID="), "TU", "model"); era.gff_gene$gene_name     


era.gff_gene$start.offset <- era.gff_gene$start + era.offset[era.gff_gene$scaffold]; era.gff_gene$end.offset <- era.gff_gene$end + era.offset[era.gff_gene$scaffold]

#era.gene_Intervals <- reduce(Intervals(as.matrix(era.gff_gene[,c("start","end")]), closed=c(TRUE,TRUE), type="Z"))
#write.table(era.gff_gene, "local/data/local/data/data/gff/era.gff_gene.txt")
mel.gff <- read.table(pipe("cut -f 1,3,4,5,9 local/data/local/data/local/data/local/data/0.ref/Hmel2.5.gff3"),
                      sep = "\t", as.is = T, col.names = c("scaffold", "feature", "start", "end", "info"))
mel.gff_gene <- subset(mel.gff, feature == "gene")
mel.gff_gene$gene_name <- sub(".*stable_id=", "", mel.gff_gene$info); mel.gff_gene$gene_name     
mel.gff_gene$start.offset <- mel.gff_gene$start + mel.offset[mel.gff_gene$scaffold]; mel.gff_gene$end.offset <- mel.gff_gene$end + mel.offset[mel.gff_gene$scaffold]
#mel.gene_Intervals <- reduce(Intervals(as.matrix(mel.gff_gene[,c("start","end")]), closed=c(TRUE,TRUE), type="Z"))


####### load pbs with rec #######
era.all.pbs <-  read.csv("era.all.pbs.recomb.csv"); head(era.all.pbs)
####### load thetas sorted #####
thetas.era.hig.vlo.pop <- read.csv("local/data/thetas/thetas.era.hig.vlo.pop.csv"); head(thetas.era.hig.vlo.pop)
thetas.era.hig.vlo.pop$tp.mean <- thetas.era.hig.vlo.pop$tP / thetas.era.hig.vlo.pop$nSites


### summarise (so mean pi of pops) by alt.type per pos, then spread by altitude, then plot deltapi
pi.era.all.alt <- summarise(group_by(thetas.era.hig.vlo.pop, chr, WinCenter,BP.wg, shdr.para.east.id, shdr.para.west.id, alt.type, cline,  background),
                            #Tajima=mean(Tajima),
                            pi=mean(tP/nSites)); head(pi.era.all.alt)
# spread
pi.era.all.alt.wide <-spread(pi.era.all.alt, alt.type, pi)
pi.era.all.alt.wide$side <- substr(pi.era.all.alt.wide$cline, 8,9); pi.era.all.alt.wide$side

# calculate delta pi (if not subsetting hig low, can do the others)
#pi.era.all.alt.wide$delta.pi.hig.low <- pi.era.all.alt.wide$hig-pi.era.all.alt.wide$low
#pi.era.all.alt.wide$delta.pi.low.vlo <- pi.era.all.alt.wide$low-pi.era.all.alt.wide$vlo
pi.era.all.alt.wide$delta.pi.hig.vlo <- pi.era.all.alt.wide$hig-pi.era.all.alt.wide$vlo
pi.era.all.alt.wide$shdr.allo.id <- era.shdr.allo.summ$shdr.allo.id[match(pi.era.all.alt.wide$shdr.para.west.id, era.shdr.allo.summ$shdr.para.west)]
pi.era.all.alt.wide[is.na(pi.era.all.alt.wide$shdr.allo.id),]$shdr.allo.id <- era.shdr.allo.summ$shdr.allo.id[match(pi.era.all.alt.wide[is.na(pi.era.all.alt.wide$shdr.allo.id),]$shdr.para.west.id, era.shdr.allo.summ$shdr.para.west)]
head(pi.era.all.alt.wide); unique(pi.era.all.alt.wide$shdr.allo.id)


####### load dxy ####
dxy.era.all <- read.csv("local/data/dxy/dxy.era.all.csv")

######################################## 1. BARPLOT NUMBER OUTLIERS ##########


# what percentage of hdrs have 1 2 or 3

era.stats.summ.p <- ggplot(data=shdr.era.outlier.df.summ, aes(fill=no.stats.outlier.min.max.10th.90th.sims, x=is.allo.hdr, y=n, label=total.shdr)) + 
  geom_col(aes(color=is.allo.hdr.side), size=1)+ 
  ylab("Number of SHDR") + 
  scale_color_manual(values=c( "#049E73","#0372B2","#D65D00","#D65D00" )) +
  scale_y_continuous( expand = c(0,0), limits = c(0,100))+
  scale_fill_continuous(low="white", high="black") +theme_bw()+facet_wrap(~side)+
  theme(axis.title.y = element_text(size=14,vjust=-1),
        plot.margin = unit(c(1,0.2,1.5,0.2), "lines"),axis.text.x = element_blank(), #axis.text.x = element_blank(),
        axis.text.y = element_text(size=10),axis.title.x = element_blank(),  axis.ticks.x = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_rect(fill = "transparent",  size = 1),
        panel.background = element_rect(fill = "transparent"), plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
        legend.background = element_rect(fill = "transparent"), legend.box.background = element_rect(fill = "transparent"), strip.background.y = element_blank(),
        strip.text = element_text(size=14, face = "bold"), strip.background = element_blank(), legend.position = "none"); era.stats.summ.p
ggsave("local/plots/shdr.summ/era.e.w.outlier.number.min.max.10th.90th.sims.tajima.deltapi.dxy.barplot.png", height = 3, width = 3)


mel.stats.summ.p <- ggplot(data=shdr.mel.outlier.df.summ, aes(fill=no.stats.outlier.min.max.10th.90th.sims, x=is.allo.hdr, y=n, label=total.shdr)) + 
  geom_col(aes(color=is.allo.hdr.side), size=1)+ 
  ylab("Number of SHDR") + 
  scale_color_manual(values=c( "#049E73","#0372B2","#D65D00","#D65D00" )) +
  scale_y_continuous( expand = c(0,0), limits = c(0,100))+
  scale_fill_continuous(low="white", high="black") +theme_bw()+facet_wrap(~side)+
  theme(axis.title.y = element_text(size=14,vjust=-1),
        plot.margin = unit(c(1.5,0.2,1,0.2), "lines"),axis.text.x = element_blank(), #axis.text.x = element_blank(),
        axis.text.y = element_text(size=10),axis.title.x = element_blank(),  axis.ticks.x = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_rect(fill = "transparent", , size = 1),
        panel.background = element_rect(fill = "transparent"), plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
        legend.background = element_rect(fill = "transparent"), legend.box.background = element_rect(fill = "transparent"), strip.background.y = element_blank(),
        strip.text = element_text(size=14, face = "bold"), strip.background = element_blank(), legend.position = "none"); mel.stats.summ.p
ggsave("local/plots/shdr.summ/mel.e.w.outlier.number.min.max.10th.90th.sims.tajima.deltapi.dxy.barplot.png", height = 3, width = 3)


era.stats.summ.p
######################################## 2. ZOOM IN  ##########
# arrange no space between https://stackoverflow.com/questions/55151531/ggplot2-arrange-multiple-plots-all-the-same-size-no-gaps-in-between

####################### 2.1 era para/allo region #######################
########### shdr east 005 ###########
# plot lines for east col and east ecu
#### PBS #### 
# start and end of region above 4 zpbs
start <- min(subset(era.all.pbs, shdr.para.east.id=="shdr.east.005"&zPBS0.1>4)$BP.wg)
end <- max(subset(era.all.pbs, shdr.para.east.id=="shdr.east.005"&zPBS0.1>4)$BP.wg)

head(era.gff_gene)
era.para.east.p <- ggplot(data=subset(era.all.pbs, shdr.para.east.id=="shdr.east.005"), aes(x=midPos, y=zPBS0.1))+
  annotate(geom = "rect", xmin = start, xmax = end, ymin = -Inf, ymax=Inf, fill="grey90")+
  geom_rect(inherit.aes = F,data=subset(era.gff_gene, start.offset>min(subset(era.all.pbs, shdr.para.east.id=="shdr.east.005")$BP.wg) & 
                                          end.offset<max(subset(era.all.pbs, shdr.para.east.id=="shdr.east.005")$BP.wg)), aes(xmin=start.offset, xmax=end.offset, ymin=13.3,ymax=14.4) )+
  geom_rect(inherit.aes = F,data=subset(era.gff_gene, gene_name=="evm.model.Herato0101.551"), color=NA, fill="purple", aes(xmin=start.offset, xmax=end.offset, ymin=13.3,ymax=14.4) )+
  geom_rect(inherit.aes = F,data=subset(era.gff_gene, gene_name=="evm.model.Herato0101.552"), color=NA, fill="purple", aes(xmin=start.offset, xmax=end.offset, ymin=13.3,ymax=14.4) )+
  annotate(geom = "text", y=15.5, x=(subset(era.gff_gene, gene_name=="evm.model.Herato0101.551")$start.offset +subset(era.gff_gene, gene_name=="evm.model.Herato0101.551")$end.offset)/2,
            label=expression(italic("Col4α1, vkg")) , size=5, colour="purple") +
  annotate(geom = "text", y=15.5, x=(subset(era.gff_gene, gene_name=="evm.model.Herato0101.552")$start.offset +subset(era.gff_gene, gene_name=="evm.model.Herato0101.552")$end.offset)/2,
           label=expression(italic("PITH")) , size=5, colour="purple") +
  #geom_point( shape=24, fill="#049E73", size=3)+
  geom_line(aes(y=rollmean(zPBS0.1, 1, na.pad=TRUE )), color="#049E73", size=2, alpha=.8) +
  #geom_point(inherit.aes = F, data=subset(era.all.pbs, shdr.para.east.id=="shdr.east.005"), aes(x=midPos, y=zPBS0.3), fill="#049E73", shape=25, size=3, color="grey")+
  geom_line(inherit.aes = F, data=subset(era.all.pbs, shdr.para.east.id=="shdr.east.005"), aes(x=midPos, y=rollmean(zPBS0.3, 3, na.pad=TRUE )), color="#049E73", size=2, alpha=.8, lty=c("11")) +
  ylab("zPBS/zFst highs")+
  theme_bw()+ ylim(-1.05,16)+
  scale_x_continuous( expand = c(0, 0)) +
  theme(axis.title.y = element_text(size=14,vjust=-2),rect = element_rect(fill = "transparent") ,
        plot.margin = unit(c(0,0.2,0,0.2), "lines"),axis.text =element_blank(), #axis.text.x = element_blank(),
        axis.text.y = element_text(size=10),axis.title.x = element_blank(),  axis.ticks.length=unit(-0.2, "cm"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_rect(fill = "transparent",  color = "#049E73", size = 2),
        panel.background = element_rect(fill = "transparent",color="#049E73", size = 1), plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
        legend.background = element_rect(fill = "transparent"), legend.box.background = element_rect(fill = "transparent",color="#049E73", size = 1), strip.background.y = element_blank(),
        strip.text = element_text(size=14, face = "bold"), strip.background = element_blank(), legend.position = "none"); era.para.east.p

####  delta.pi #### 
head(pi.era.all.alt.wide)
era.para.east.delta.pi.p <- ggplot(data=subset(pi.era.all.alt.wide, shdr.para.east.id=="shdr.east.005" &(cline=="era.co.e"|cline=="era.ec.e")), 
                             aes(x=BP.wg, y=delta.pi.hig.vlo, color=cline, fill=cline, group=cline))+
  annotate(geom = "rect", xmin = start, xmax = end, ymin = -Inf, ymax=Inf, fill="grey90")+
  geom_line(aes(y=rollmean(delta.pi.hig.vlo, 3, na.pad=TRUE ), linetype=cline), size=2, alpha=.8) +
  scale_color_manual(values=c("#049E73", "#049E73"))+
  scale_linetype_manual(values=c('solid','11'))+
  ylab(c(expression(Delta*pi*' high - low')))+
  theme_bw()+
  scale_y_continuous(n.breaks = 4)+
  scale_x_continuous( expand = c(0, 0)) +
  theme(axis.title.y = element_text(size=14,vjust=-2),
        plot.margin = unit(c(0,0.2,0,0.2), "lines"),axis.text.x = element_blank(), rect = element_rect(fill = "transparent") ,#axis.text.x = element_blank(),
        axis.text.y = element_text(size=10),axis.title.x = element_blank(),  axis.ticks.length=unit(-0.2, "cm"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_rect(fill = "transparent",  color = "#049E73", size = 2),
        panel.background = element_rect(fill = "transparent",color="#049E73", size = 1), plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
        legend.background = element_rect(fill = "transparent"), legend.box.background = element_rect(fill = "transparent",color="#049E73", size = 1), strip.background.y = element_blank(),
        strip.text = element_text(size=14, face = "bold"), strip.background = element_blank(), legend.position = "none"); era.para.east.delta.pi.p

####  high tajimas d #### 
head(thetas.era.hig.vlo.pop)
era.para.east.tajima.high.p <- ggplot(data=subset(thetas.era.hig.vlo.pop, shdr.para.east.id=="shdr.east.005" &(cline=="era.co.e"|cline=="era.ec.e") &alt.type=="hig"), 
                                   aes(x=BP.wg, y=Tajima, color=cline, fill=cline, group=cline))+
  annotate(geom = "rect", xmin = start, xmax = end, ymin = -Inf, ymax=Inf, fill="grey90")+
  #geom_point( shape=24, size=3)+
  geom_line(aes(y=rollmean(Tajima, 3, na.pad=TRUE ), linetype=cline), size=2, alpha=.8) +
  scale_color_manual(values=c("#049E73", "#049E73"))+
  scale_linetype_manual(values=c('solid','11'))+
  ylab('Tajima\'s D')+
  theme_bw()+
  scale_x_continuous( expand = c(0, 0)) +  scale_y_continuous(n.breaks = 4)+
  theme(axis.title.y = element_text(size=14,vjust=-2),
        plot.margin = unit(c(0,0.2,0,0.2), "lines"),axis.text.x = element_blank(), #axis.text.x = element_blank(),
        axis.text.y = element_text(size=10),axis.title.x = element_blank(),  axis.ticks.length=unit(-0.2, "cm"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_rect(fill = "transparent",  color = "#049E73", size = 2),
        panel.background = element_rect(fill = "transparent",color="#049E73", size = 1), plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
        legend.background = element_rect(fill = "transparent"), legend.box.background = element_rect(fill = "transparent",color="#049E73", size = 1), strip.background.y = element_blank(),
        strip.text = element_text(size=14, face = "bold"), strip.background = element_blank(), legend.position = "none"); era.para.east.tajima.high.p


####  dxy #### 
head(dxy.era.all )
unique(dxy.era.all$pop1.pop2.alt.type)
era.para.east.dxy.p <- ggplot(data=subset(dxy.era.all, shdr.para.east.id=="shdr.east.005" &(cline=="era.co.e"|cline=="era.ec.e") & pop1.pop2.alt.type=="hig_vlo"), 
                                      aes(x=BP.wg/1000000, y=dxy_mean , color=cline, fill=cline, group=cline))+
  annotate(geom = "rect", xmin = start/1000000, xmax = end/1000000, ymin = -Inf, ymax=Inf, fill="grey90")+
  geom_line(aes(y=rollmean(dxy_mean ,5 , na.pad=TRUE ), linetype=cline), size=2, alpha=.8) +
  scale_color_manual(values=c("#049E73", "#049E73"))+
  scale_linetype_manual(values=c('solid','11'))+
  ylab('Dxy high - low')+ xlab('Position (MB), Herato0101') +
  theme_bw()+
  scale_x_continuous( expand = c(0, 0)) +
  theme(axis.title.y = element_text(size=14,vjust=-2),
        plot.margin = unit(c(0,0.2,0,0.2), "lines"),axis.text.x = element_text(size=10), #axis.text.x = element_blank(),
        axis.text.y = element_text(size=10),axis.title.x = element_text(size=14), axis.ticks.length=unit(-0.2, "cm"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_rect(fill = "transparent",  color = "#049E73", size = 2),
        panel.background = element_rect(fill = "transparent",color="#049E73", size = 1), plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
        legend.background = element_rect(fill = "transparent"), legend.box.background = element_rect(fill = "transparent",color="#049E73", size = 1), strip.background.y = element_blank(),
        strip.text = element_text(size=14, face = "bold"), strip.background = element_blank(), legend.position = "none"); era.para.east.dxy.p


########### shdr.all.028  ###########
#### PBS #### 
names(era.all.pbs)
start <- min(subset(era.all.pbs, shdr.allo.id=="shdr.all.028"&(zPBS0.1>4 | zfst.2>4 | zPBS0.3>4 | zPBS0.4>4 ))$BP.wg); start
end <- max(subset(era.all.pbs, shdr.allo.id=="shdr.all.028"&(zPBS0.1>4 | zfst.2>4 | zPBS0.3>4 | zPBS0.4>4 ))$BP.wg); end

era.allo.p <- ggplot(data=subset(era.all.pbs, shdr.allo.id=="shdr.all.028"), aes(x=BP.wg, y=zPBS0.1))+
  annotate(geom = "rect", xmin = start, xmax = end, ymin = -Inf, ymax=Inf, fill="grey90")+
  geom_rect(inherit.aes = F,data=subset(era.gff_gene, start.offset>min(subset(era.all.pbs, shdr.allo.id=="shdr.all.028")$BP.wg) & 
                                          end.offset<max(subset(era.all.pbs, shdr.allo.id=="shdr.all.028")$BP.wg)), aes(xmin=start.offset, xmax=end.offset, ymin=10.8,ymax=11.6) )+
  geom_rect(inherit.aes = F,data=subset(era.gff_gene, gene_name=="evm.model.Herato1507.182"), color=NA, fill="purple", aes(xmin=start.offset, xmax=end.offset, ymin=10.8,ymax=11.6) )+
  annotate(geom = "text", y=12.5, x=(subset(era.gff_gene, gene_name=="evm.model.Herato1507.182")$start.offset +subset(era.gff_gene, gene_name=="evm.model.Herato1507.182")$end.offset)/2,
           label=expression(italic("stem cell tumor")) , size=5, colour="purple") +
  #geom_point( shape=24, fill="#049E73", size=3)+
  geom_line(aes(y=rollmean(zPBS0.1, 1, na.pad=TRUE )), color="#049E73", size=2, alpha=.8) +
  #geom_point(inherit.aes = F, data=subset(era.all.pbs, shdr.allo.id=="shdr.all.028"), aes(x=BP.wg, y=zPBS0.3), fill="#049E73", shape=25, size=3, color="grey")+
  geom_line(inherit.aes = F, data=subset(era.all.pbs, shdr.allo.id=="shdr.all.028"), aes(x=BP.wg,y=rollmean(zPBS0.3, 1, na.pad=TRUE )), color="#049E73", size=2, alpha=.8, lty=c("11")) +
  #geom_point(inherit.aes = F, data=subset(era.all.pbs, shdr.allo.id=="shdr.all.028"), aes(x=BP.wg, y=zfst.2), fill="#0372B2", shape=24, size=3, color="black")+
  geom_line(inherit.aes = F, data=subset(era.all.pbs, shdr.allo.id=="shdr.all.028"), aes(x=BP.wg,y=rollmean(zfst.2, 1, na.pad=TRUE )), color="#0372B2", size=2, alpha=.8, lty=1) +
  #geom_point(inherit.aes = F, data=subset(era.all.pbs, shdr.allo.id=="shdr.all.028"), aes(x=BP.wg, y=zPBS0.4), fill="#0372B2", shape=25, size=3, color="grey")+
  geom_line(inherit.aes = F, data=subset(era.all.pbs, shdr.allo.id=="shdr.all.028"), aes(x=BP.wg,y=rollmean(zPBS0.4, 1, na.pad=TRUE )), color="#0372B2", size=2, alpha=.8, lty=c("11")) +
  ylab("zPBS highs")+
  theme_bw()+ylim(-1.05,13)+
  scale_x_continuous( expand = c(0, 0)) +
  theme(axis.title.y = element_text(size=14,vjust=-2),
        plot.margin = unit(c(0,0.2,0,0.2), "lines"),axis.text =element_blank(), axis.text.x = element_blank(),
        axis.text.y = element_text(size=10),axis.title.x = element_blank(),  axis.ticks.length=unit(-0.2, "cm"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_rect(fill = "transparent",  color = "#D65D00", size = 2),
        panel.background = element_rect(fill = "transparent"), plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
        legend.background = element_rect(fill = "transparent"), legend.box.background = element_rect(fill = "transparent"), strip.background.y = element_blank(),
        strip.text = element_text(size=14, face = "bold"), strip.background = element_blank(), legend.position = "none"); era.allo.p

####  delta.pi #### 
head(pi.era.all.alt.wide)
era.allo.delta.pi.p <- ggplot(data=subset(pi.era.all.alt.wide, shdr.allo.id=="shdr.all.028" ), 
                                   aes(x=BP.wg, y=delta.pi.hig.vlo, color=cline, fill=cline, group=cline))+
  annotate(geom = "rect", xmin = start, xmax = end, ymin = -Inf, ymax=Inf, fill="grey90")+
  #geom_hline(yintercept = 0) +
  geom_line(aes(y=rollmean(delta.pi.hig.vlo, 1, na.pad=TRUE, group=cline ), linetype=cline), size=2, alpha=.8) +
  scale_color_manual(values=c("#049E73", "#049E73", "#0372B2", "#0372B2"))+
  scale_linetype_manual(values=c('solid','11','solid','11'))+
  ylab(c(expression(Delta*pi*' high - low')))+
  theme_bw()+
  scale_y_continuous(n.breaks = 4)+
  scale_x_continuous( expand = c(0, 0)) +
  theme(axis.title.y = element_text(size=14,vjust=-2),
        plot.margin = unit(c(0,0.2,0,0.2), "lines"),axis.text.x = element_blank(), #axis.text.x = element_blank(),
        axis.text.y = element_text(size=10),axis.title.x = element_blank(),  axis.ticks.length=unit(-0.2, "cm"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_rect(fill = "transparent",  color = "#D65D00", size = 2),
        panel.background = element_rect(fill = "transparent"), plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
        legend.background = element_rect(fill = "transparent"), legend.box.background = element_rect(fill = "transparent"), strip.background.y = element_blank(),
        strip.text = element_text(size=14, face = "bold"), strip.background = element_blank(), legend.position = "none"); era.allo.delta.pi.p

####  high tajimas d #### 
head(thetas.era.hig.pop )
era.allo.tajima.high.p <- ggplot(data=subset(thetas.era.hig.vlo.pop, shdr.allo.id=="shdr.all.028"&alt.type=="hig" ), 
                                      aes(x=BP.wg, y=Tajima, color=cline, fill=cline, group=cline))+
  annotate(geom = "rect", xmin = start, xmax = end, ymin = -Inf, ymax=Inf, fill="grey90")+
  geom_line(aes(y=rollmean(Tajima, 1, na.pad=TRUE ), linetype=cline), size=2, alpha=.8) +
  #geom_point( shape=24, size=3)+
  scale_color_manual(values=c("#049E73", "#049E73", "#0372B2", "#0372B2"))+
  scale_linetype_manual(values=c('solid','11','solid','11'))+
  ylab('Tajima\'s D')+
  theme_bw()+
  scale_x_continuous( expand = c(0, 0)) +  scale_y_continuous(n.breaks = 4)+
  theme(axis.title.y = element_text(size=14,vjust=-2),
        plot.margin = unit(c(0,0.2,0,0.2), "lines"),axis.text.x = element_blank(), #axis.text.x = element_blank(),
        axis.text.y = element_text(size=10),axis.title.x = element_blank(),  axis.ticks.length=unit(-0.2, "cm"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_rect(fill = "transparent",  color = "#D65D00", size = 2),
        panel.background = element_rect(fill = "transparent"), plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
        legend.background = element_rect(fill = "transparent"), legend.box.background = element_rect(fill = "transparent"), strip.background.y = element_blank(),
        strip.text = element_text(size=14, face = "bold"), strip.background = element_blank(), legend.position = "none"); era.allo.tajima.high.p


####  dxy #### 
head(dxy.era.all )
unique(dxy.era.all$pop1.pop2.alt.type)
start.bp <- min(subset(era.all.pbs, shdr.allo.id=="shdr.all.028"&(zPBS0.1>4 | zfst.2>4 | zPBS0.3>4 | zPBS0.4>4 ))$midPos); start.bp
end.bp <- max(subset(era.all.pbs, shdr.allo.id=="shdr.all.028"&(zPBS0.1>4 | zfst.2>4 | zPBS0.3>4 | zPBS0.4>4 ))$midPos); end.bp

era.allo.dxy.p <- ggplot(data=subset(dxy.era.all, shdr.allo.id=="shdr.all.028" & pop1.pop2.alt.type=="hig_vlo"), 
                              aes(x=win_mid/1000000, y=dxy_mean , color=cline,shape=cline, fill=cline, group=cline))+
  annotate(geom = "rect", xmin = start.bp/1000000, xmax = end.bp/1000000, ymin = -Inf, ymax=Inf, fill="grey90")+
  geom_line(aes(y=rollmean(dxy_mean ,1 , na.pad=TRUE ), linetype=cline), size=2, alpha=.8) +
  #geom_point()+
  scale_color_manual(values=c("#049E73", "#049E73", "#0372B2", "#0372B2"))+
  scale_linetype_manual(values=c('solid','11','solid','11'))+
  ylab('Dxy high - low')+ xlab('Position (MB), Herato1507') +
  theme_bw()+
  scale_x_continuous( expand = c(0, 0)) +
  theme(axis.title.y = element_text(size=14,vjust=-2),
        plot.margin = unit(c(0,0.2,0,0.2), "lines"),axis.text.x = element_text(size=10), #axis.text.x = element_blank(),
        axis.text.y = element_text(size=10),axis.title.x = element_text(size=14,margin = margin(t = 0, r = 0, b = 0, l = 0)), axis.ticks.length=unit(-0.2, "cm"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_rect(fill = "transparent",  color = "#D65D00", size = 2),
        panel.background = element_rect(fill = "transparent"), plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
        legend.background = element_rect(fill = "transparent"), legend.box.background = element_rect(fill = "transparent"), strip.background.y = element_blank(),
        strip.text = element_text(size=14, face = "bold"), strip.background = element_blank(), legend.position = "none"); era.allo.dxy.p

#### plot all #### 

left1 <- (era.para.east.p/ era.para.east.delta.pi.p /era.para.east.tajima.high.p /era.para.east.dxy.p ) 
left2 <- (era.allo.p/ era.allo.delta.pi.p /era.allo.tajima.high.p /era.allo.dxy.p)

((era.stats.summ.p / mel.stats.summ.p)  | 
  left1 | left2 + 
       plot_layout(widths = c(0.35, 1, 1), heights =c(1,1,1,1) ) )&
  plot_annotation(theme = theme(plot.background = element_rect(color  = 'transparent', size = 2,linetype = 'dotted', fill ="transparent")))

ggsave("~/Dropbox/PhD/22_pop.gen.paper/figures/fig3/fig3.era.mel.selection.png", width = 11, height = 7,  bg = "transparent")







####################### 2.1 era para/allo region no green lines #######################
########### shdr east 005 ###########
# plot lines for east col and east ecu
#### PBS #### 
# start and end of region above 4 zpbs
start <- min(subset(era.all.pbs, shdr.para.east.id=="shdr.east.005"&zPBS0.1>4)$BP.wg)
end <- max(subset(era.all.pbs, shdr.para.east.id=="shdr.east.005"&zPBS0.1>4)$BP.wg)

head(era.gff_gene)
era.para.east.p <- ggplot(data=subset(era.all.pbs, shdr.para.east.id=="shdr.east.005"), aes(x=midPos, y=zPBS0.1))+
  annotate(geom = "rect", xmin = start, xmax = end, ymin = -Inf, ymax=Inf, fill="#049E73", alpha=.3)+
  geom_rect(inherit.aes = F,data=subset(era.gff_gene, start.offset>min(subset(era.all.pbs, shdr.para.east.id=="shdr.east.005")$BP.wg) & 
                                          end.offset<max(subset(era.all.pbs, shdr.para.east.id=="shdr.east.005")$BP.wg)), aes(xmin=start.offset, xmax=end.offset, ymin=13.3,ymax=14.4) )+
  geom_rect(inherit.aes = F,data=subset(era.gff_gene, gene_name=="evm.model.Herato0101.551"), color=NA, fill="purple", aes(xmin=start.offset, xmax=end.offset, ymin=13.3,ymax=14.4) )+
  geom_rect(inherit.aes = F,data=subset(era.gff_gene, gene_name=="evm.model.Herato0101.552"), color=NA, fill="purple", aes(xmin=start.offset, xmax=end.offset, ymin=13.3,ymax=14.4) )+
  annotate(geom = "text", y=15.5, x=(subset(era.gff_gene, gene_name=="evm.model.Herato0101.551")$start.offset +subset(era.gff_gene, gene_name=="evm.model.Herato0101.551")$end.offset)/2,
           label=expression(italic("Col4α1, vkg")) , size=5, colour="purple") +
  annotate(geom = "text", y=15.5, x=(subset(era.gff_gene, gene_name=="evm.model.Herato0101.552")$start.offset +subset(era.gff_gene, gene_name=="evm.model.Herato0101.552")$end.offset)/2,
           label=expression(italic("PITH")) , size=5, colour="purple") +
  #geom_point( shape=24, fill="black", size=3)+
  geom_line(aes(y=rollmean(zPBS0.1, 1, na.pad=TRUE )), color="black", size=1, alpha=.8) +
  #geom_point(inherit.aes = F, data=subset(era.all.pbs, shdr.para.east.id=="shdr.east.005"), aes(x=midPos, y=zPBS0.3), fill="black", shape=25, size=3, color="grey")+
  geom_line(inherit.aes = F, data=subset(era.all.pbs, shdr.para.east.id=="shdr.east.005"), aes(x=midPos, y=rollmean(zPBS0.3, 1, na.pad=TRUE )), color="black", size=1, alpha=.8, lty=c("11")) +
  ylab("zPBS high / zFst")+
  theme_bw()+ ylim(-1.05,16)+
  scale_x_continuous( expand = c(0, 0)) +
  theme(axis.title.y = element_text(size=14,vjust=-2),rect = element_rect(fill = "transparent") ,
        plot.margin = unit(c(0,0.2,0,0.2), "lines"),axis.text =element_blank(), #axis.text.x = element_blank(),
        axis.text.y = element_text(size=10),axis.title.x = element_blank(),  axis.ticks.length=unit(-0.2, "cm"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_rect(fill = "transparent",  color = "black", size = 1),
        panel.background = element_rect(fill = "transparent",color="black", size = 1), plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
        legend.background = element_rect(fill = "transparent"), legend.box.background = element_rect(fill = "transparent",color="black", size = 1), strip.background.y = element_blank(),
        strip.text = element_text(size=14, face = "bold"), strip.background = element_blank(), legend.position = "none"); era.para.east.p

####  delta.pi #### 
head(pi.era.all.alt.wide)
era.para.east.delta.pi.p <- ggplot(data=subset(pi.era.all.alt.wide, shdr.para.east.id=="shdr.east.005" &(cline=="era.co.e"|cline=="era.ec.e")), 
                                   aes(x=BP.wg, y=delta.pi.hig.vlo, color=cline, fill=cline, group=cline))+
  annotate(geom = "rect", xmin = start, xmax = end, ymin = -Inf, ymax=Inf, fill="#049E73", alpha=.3)+
  geom_line(aes(y=rollmean(delta.pi.hig.vlo, 1, na.pad=TRUE ), linetype=cline), size=1, alpha=.8) +
  #geom_point()+
  scale_color_manual(values=c("black", "black"))+
  scale_linetype_manual(values=c('solid','11'))+
  ylab(c(expression(Delta*pi*' high - low')))+
  theme_bw()+
  scale_y_continuous(n.breaks = 4)+
  scale_x_continuous( expand = c(0, 0)) +
  theme(axis.title.y = element_text(size=14,vjust=-1),
        plot.margin = unit(c(0,0.2,0,0.2), "lines"),axis.text.x = element_blank(), rect = element_rect(fill = "transparent") ,#axis.text.x = element_blank(),
        axis.text.y = element_text(size=9),axis.title.x = element_blank(),  axis.ticks.length=unit(-0.2, "cm"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_rect(fill = "transparent",  color = "black", size = 1),
        panel.background = element_rect(fill = "transparent",color="black", size = 1), plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
        legend.background = element_rect(fill = "transparent"), legend.box.background = element_rect(fill = "transparent",color="black", size = 1), strip.background.y = element_blank(),
        strip.text = element_text(size=14, face = "bold"), strip.background = element_blank(), legend.position = "none"); era.para.east.delta.pi.p

####  high tajimas d #### 
head(thetas.era.hig.vlo.pop)
era.para.east.tajima.high.p <- ggplot(data=subset(thetas.era.hig.vlo.pop, shdr.para.east.id=="shdr.east.005" &(cline=="era.co.e"|cline=="era.ec.e") &alt.type=="hig"), 
                                      aes(x=BP.wg, y=Tajima, color=cline, fill=cline, group=cline))+
  annotate(geom = "rect", xmin = start, xmax = end, ymin = -Inf, ymax=Inf, fill="#049E73", alpha=.3)+
  #geom_point( shape=24, size=3)+
  geom_line(aes(y=rollmean(Tajima, 1, na.pad=TRUE ), linetype=cline), size=1, alpha=.8) +
  scale_color_manual(values=c("black", "black"))+
  scale_linetype_manual(values=c('solid','11'))+
  ylab('Tajima\'s D')+
  theme_bw()+
  scale_x_continuous( expand = c(0, 0)) +  scale_y_continuous(n.breaks = 4)+
  theme(axis.title.y = element_text(size=14,vjust=-2),
        plot.margin = unit(c(0,0.2,0,0.2), "lines"),axis.text.x = element_blank(), #axis.text.x = element_blank(),
        axis.text.y = element_text(size=10),axis.title.x = element_blank(),  axis.ticks.length=unit(-0.2, "cm"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_rect(fill = "transparent",  color = "black", size = 1),
        panel.background = element_rect(fill = "transparent",color="black", size = 1), plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
        legend.background = element_rect(fill = "transparent"), legend.box.background = element_rect(fill = "transparent",color="black", size = 1), strip.background.y = element_blank(),
        strip.text = element_text(size=14, face = "bold"), strip.background = element_blank(), legend.position = "none"); era.para.east.tajima.high.p


####  dxy #### 
head(dxy.era.all )
unique(dxy.era.all$pop1.pop2.alt.type)
era.para.east.dxy.p <- ggplot(data=subset(dxy.era.all, shdr.para.east.id=="shdr.east.005" &(cline=="era.co.e"|cline=="era.ec.e") & pop1.pop2.alt.type=="hig_vlo"), 
                              aes(x=BP.wg/1000000, y=dxy_mean , color=cline, fill=cline, group=cline))+
  annotate(geom = "rect", xmin = start/1000000, xmax = end/1000000, ymin = -Inf, ymax=Inf, fill="#049E73", alpha=.3)+
  geom_line(aes(y=rollmean(dxy_mean ,1 , na.pad=TRUE ), linetype=cline), size=1, alpha=.8) +
  #geom_point( shape=24, size=3)+
  scale_color_manual(values=c("black", "black"))+
  scale_linetype_manual(values=c('solid','11'))+
  ylab('Dxy')+ xlab('Position (MB), Herato0101') +
  theme_bw()+
  scale_x_continuous( expand = c(0, 0)) +
  theme(axis.title.y = element_text(size=14,vjust=-2),
        plot.margin = unit(c(0,0.6,0,0.2), "lines"),axis.text.x = element_text(size=10), #axis.text.x = element_blank(),
        axis.text.y = element_text(size=10),axis.title.x = element_text(size=14), axis.ticks.length=unit(-0.2, "cm"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_rect(fill = "transparent",  color = "black", size = 1),
        panel.background = element_rect(fill = "transparent",color="black", size = 1), plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
        legend.background = element_rect(fill = "transparent"), legend.box.background = element_rect(fill = "transparent",color="black", size = 1), strip.background.y = element_blank(),
        strip.text = element_text(size=14, face = "bold"), strip.background = element_blank(), legend.position = "none"); era.para.east.dxy.p

#### plot all #### 

left1 <- (era.para.east.p/ era.para.east.delta.pi.p /era.para.east.tajima.high.p /era.para.east.dxy.p ) ; left1

ggsave("~/Dropbox/PhD/22_pop.gen.paper/figures/fig3/fig3.era.mel.selection.thin.png", width = 4, height = 7,  bg = "transparent")







######################################## 3. pbs value with number outlier ######################################## 
names(shdr.era.east.outlier.df)
library(ggpubr); library(ggbeeswarm); library(EnvStats)
co.pbs.east.p.era <- ggplot(data=subset(shdr.era.east.outlier.df, is.max.pbs.hig.above4=="yes"), 
                        aes(y=zPBS.hig.co.e.max.value, x=no.stats.outlier.min.max.10th.90th.sims))+
  geom_boxplot(aes(group=as.factor(no.stats.outlier.min.max.10th.90th.sims)))+
  geom_beeswarm(aes(colour=is.allo.hdr), size=2, alpha=0.7)+ 
  scale_color_manual(values=c("#049E73", "#D65D00"))+
  #stat_cor()+   geom_smooth(method="lm", color="black")+
  stat_n_text() +
  ylab("zPBS highland Colombia East") +
  theme(axis.title.y = element_text(size=10),
        plot.margin = unit(c(0,0.2,0,0.2), "lines"),axis.text.x = element_text(size=10),
        axis.text.y = element_text(size=10), axis.title.x= element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_rect(fill = "transparent"),
        panel.background = element_rect(fill = "transparent"), plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
        legend.background = element_rect(fill = "transparent"), legend.box.background = element_rect(fill = "transparent"), strip.background.y = element_blank(),
        strip.text = element_text(size=14, face = "bold"), strip.background = element_blank(), legend.position = "none"); co.pbs.east.p.era

ec.pbs.east.p.era <- ggplot(data=subset(shdr.era.east.outlier.df, is.max.pbs.hig.above4=="yes"), 
       aes(y=zPBS.hig.ec.e.max.value, x=no.stats.outlier.min.max.10th.90th.sims))+
  geom_boxplot(aes(group=as.factor(no.stats.outlier.min.max.10th.90th.sims)))+
  geom_beeswarm(aes(colour=is.allo.hdr), size=2, alpha=0.7)+ 
  scale_color_manual(values=c("#049E73", "#D65D00"))+
  #stat_cor()+   geom_smooth(method="lm", color="black")+
  stat_n_text() +
  ylab("zPBS highland Ecuador East") +
  theme(axis.title.y = element_text(size=10),
        plot.margin = unit(c(0,0.2,0,0.2), "lines"),axis.text.x = element_text(size=10),
        axis.text.y = element_text(size=10), axis.title.x= element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_rect(fill = "transparent"),
        panel.background = element_rect(fill = "transparent"), plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
        legend.background = element_rect(fill = "transparent"), legend.box.background = element_rect(fill = "transparent"), strip.background.y = element_blank(),
        strip.text = element_text(size=14, face = "bold"), strip.background = element_blank(), legend.position = "none"); ec.pbs.east.p.era

names(shdr.era.west.outlier.df)
co.pbs.west.p.era <- ggplot(data=subset(shdr.era.west.outlier.df, is.max.pbs.hig.above4=="yes"), 
                        aes(y=zFst.co.w.max.value, x=no.stats.outlier.min.max.10th.90th.sims))+
  geom_boxplot(aes(group=as.factor(no.stats.outlier.min.max.10th.90th.sims)))+
  geom_beeswarm(aes(colour=is.allo.hdr), size=2, alpha=0.7)+ 
  scale_color_manual(values=c("#0372B2", "#D65D00"))+
  #stat_cor()+   geom_smooth(method="lm", color="black")+
  stat_n_text() +
  ylab("zFst Colombia West") +
  theme(axis.title.y = element_text(size=10),
        plot.margin = unit(c(0,0.2,0,0.2), "lines"),axis.text.x = element_text(size=10),
        axis.text.y = element_text(size=10), axis.title.x= element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_rect(fill = "transparent"),
        panel.background = element_rect(fill = "transparent"), plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
        legend.background = element_rect(fill = "transparent"), legend.box.background = element_rect(fill = "transparent"), strip.background.y = element_blank(),
        strip.text = element_text(size=14, face = "bold"), strip.background = element_blank(), legend.position = "none"); co.pbs.west.p.era

ec.pbs.west.p.era <- ggplot(data=subset(shdr.era.west.outlier.df, is.max.pbs.hig.above4=="yes"), 
                        aes(y=zPBS.hig.ec.w.max.value, x=no.stats.outlier.min.max.10th.90th.sims))+
  geom_boxplot(aes(group=as.factor(no.stats.outlier.min.max.10th.90th.sims)))+
  geom_beeswarm(aes(colour=is.allo.hdr), size=2, alpha=0.7)+ 
  scale_color_manual(values=c("#0372B2", "#D65D00"))+
  #stat_cor()+   geom_smooth(method="lm", color="black")+
  stat_n_text() +
  ylab("zPBS highland Ecuador West") +
  theme(axis.title.y = element_text(size=10),
        plot.margin = unit(c(0,0.2,0,0.2), "lines"),axis.text.x = element_text(size=10),
        axis.text.y = element_text(size=10), axis.title.x= element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_rect(fill = "transparent"),
        panel.background = element_rect(fill = "transparent"), plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
        legend.background = element_rect(fill = "transparent"), legend.box.background = element_rect(fill = "transparent"), strip.background.y = element_blank(),
        strip.text = element_text(size=14, face = "bold"), strip.background = element_blank(), legend.position = "none"); ec.pbs.west.p.era



co.pbs.east.p.mel <- ggplot(data=subset(shdr.mel.east.outlier.df, is.max.pbs.hig.above4=="yes"), 
                            aes(y=zPBS.hig.co.e.max.value, x=no.stats.outlier.min.max.10th.90th.sims))+
  geom_boxplot(aes(group=as.factor(no.stats.outlier.min.max.10th.90th.sims)))+
  geom_beeswarm(aes(colour=is.allo.hdr), size=2, alpha=0.7)+ 
  scale_color_manual(values=c("#049E73", "#D65D00"))+
  #stat_cor()+   geom_smooth(method="lm", color="black")+
  stat_n_text() +
  ylab("zPBS highland Colombia East") +
  theme(axis.title.y = element_text(size=10),
        plot.margin = unit(c(0,0.2,0,0.2), "lines"),axis.text.x = element_text(size=10),
        axis.text.y = element_text(size=10), axis.title.x= element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_rect(fill = "transparent"),
        panel.background = element_rect(fill = "transparent"), plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
        legend.background = element_rect(fill = "transparent"), legend.box.background = element_rect(fill = "transparent"), strip.background.y = element_blank(),
        strip.text = element_text(size=14, face = "bold"), strip.background = element_blank(), legend.position = "none"); co.pbs.east.p.mel

ec.pbs.east.p.mel <- ggplot(data=subset(shdr.mel.east.outlier.df, is.max.pbs.hig.above4=="yes"), 
                            aes(y=zPBS.hig.ec.e.max.value, x=no.stats.outlier.min.max.10th.90th.sims))+
  geom_boxplot(aes(group=as.factor(no.stats.outlier.min.max.10th.90th.sims)))+
  geom_beeswarm(aes(colour=is.allo.hdr), size=2, alpha=0.7)+ 
  scale_color_manual(values=c("#049E73", "#D65D00"))+
  #stat_cor()+   geom_smooth(method="lm", color="black")+
  stat_n_text() +
  ylab("zPBS highland Ecuador East") +
  theme(axis.title.y = element_text(size=10),
        plot.margin = unit(c(0,0.2,0,0.2), "lines"),axis.text.x = element_text(size=10),
        axis.text.y = element_text(size=10), axis.title.x= element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_rect(fill = "transparent"),
        panel.background = element_rect(fill = "transparent"), plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
        legend.background = element_rect(fill = "transparent"), legend.box.background = element_rect(fill = "transparent"), strip.background.y = element_blank(),
        strip.text = element_text(size=14, face = "bold"), strip.background = element_blank(), legend.position = "none"); ec.pbs.east.p.mel

names(shdr.mel.west.outlier.df)
co.pbs.west.p.mel <- ggplot(data=subset(shdr.mel.west.outlier.df, is.max.pbs.hig.above4=="yes"), 
                            aes(y=zPBS.hig.co.w.max.value, x=no.stats.outlier.min.max.10th.90th.sims))+
  geom_boxplot(aes(group=as.factor(no.stats.outlier.min.max.10th.90th.sims)))+
  geom_beeswarm(aes(colour=is.allo.hdr), size=2, alpha=0.7)+ 
  scale_color_manual(values=c("#0372B2", "#D65D00"))+
  #stat_cor()+   geom_smooth(method="lm", color="black")+
  stat_n_text() +
  xlim(-0.5,3.5)+
  ylab("zPBS highland Colombia West") +
  theme(axis.title.y = element_text(size=10),
        plot.margin = unit(c(0,0.2,0,0.2), "lines"),axis.text.x = element_text(size=10),
        axis.text.y = element_text(size=10), axis.title.x= element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_rect(fill = "transparent"),
        panel.background = element_rect(fill = "transparent"), plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
        legend.background = element_rect(fill = "transparent"), legend.box.background = element_rect(fill = "transparent"), strip.background.y = element_blank(),
        strip.text = element_text(size=14, face = "bold"), strip.background = element_blank(), legend.position = "none"); co.pbs.west.p.mel

ec.pbs.west.p.mel <- ggplot(data=subset(shdr.mel.west.outlier.df, is.max.pbs.hig.above4=="yes"), 
                            aes(y=zFst.ec.w.max.value, x=no.stats.outlier.min.max.10th.90th.sims))+
  geom_boxplot(aes(group=as.factor(no.stats.outlier.min.max.10th.90th.sims)))+
  geom_beeswarm(aes(colour=is.allo.hdr), size=2, alpha=0.7)+ 
  scale_color_manual(values=c("#0372B2", "#D65D00"))+
  #stat_cor()+   geom_smooth(method="lm", color="black")+
  stat_n_text() +
  ylab("zFst Ecuador West") +
  xlim(-0.5,3.5)+
  theme(axis.title.y = element_text(size=10),
        plot.margin = unit(c(0,0.2,0,0.2), "lines"),axis.text.x = element_text(size=10),
        axis.text.y = element_text(size=10), axis.title.x= element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_rect(fill = "transparent"),
        panel.background = element_rect(fill = "transparent"), plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
        legend.background = element_rect(fill = "transparent"), legend.box.background = element_rect(fill = "transparent"), strip.background.y = element_blank(),
        strip.text = element_text(size=14, face = "bold"), strip.background = element_blank(), legend.position = "none"); ec.pbs.west.p.mel

(((co.pbs.west.p.era / ec.pbs.west.p.era) | (co.pbs.east.p.era / ec.pbs.east.p.era )) /
  ((co.pbs.west.p.mel / ec.pbs.west.p.mel) | (co.pbs.east.p.mel / ec.pbs.east.p.mel )) ) + plot_annotation(tag_levels = 'A' )
ggsave("local/plots/shdr.summ/era.mel.pbs.vs.no.stats.outlier.png", height = 10, width = 7)





