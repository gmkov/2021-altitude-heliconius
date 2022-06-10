########### packages ########### 
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
options(scipen = 999)
library(intervals)


########### functions ###########
# expands dataset +-x windows
bp.list.expanded <- function(dat, x){
  BP.wg <- dat[["BP.wg"]]
  # x= windows (1kb), convert to bp
  y <- x*1000
  seq(from=as.numeric(as.character(BP.wg))-y, to=as.numeric(as.character(BP.wg))+y, by=1000)}

# ztransform anything, equivalent to scale(centre=TRUE, scale=TRUE)
zscore <- function(x) { 
  x<- as.numeric(as.character(x))
  (x-mean(x))/sd(x) }

## shift one row (either dir)
rowShift <- function(x, shiftLen = 1L) {
  r <- (1L + shiftLen):(length(x) + shiftLen)
  r[r<1] <- NA
  return(x[r])
}
## find start ends of blocks obtained with bp.list.expanded
findBlocks <- function(x, windowSize){
  y <- as.numeric(as.character(x$BP.wg))
  #import windowsize as KB, so make B
  windowSize <- windowSize*1000
  # subset so that the diff between row above and below
  starts <- x[((is.na(y - rowShift( y,-1)) | y - rowShift( y,-1)>windowSize)),]; starts$pos.type <- c("start")
  ends <- x[((is.na(rowShift( y,+1)- y) | rowShift( y,+1)- y>windowSize)),]; ends$pos.type <- c("end")
  arrange(rbind(starts,ends), BP.wg) # concatenate, then order by BP.wg increasing
}

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

total_intersection <- function(ints1, ints2){
  #get intersection
  intersec <- interval_intersection(ints1,ints2)
  #close ends (make inclusive)
  type(intersec) <- "R"
  intersec <- close_intervals(intersec)
  #sum all intersection lengths
  sum(size(intersec))
}

# function to check whether each of ints1 overlaps any of ints2
overlaps_anys <- function(ints1, ints2, ints3, ints4){
  overlaps <- interval_overlap(ints1,ints2, ints3, ints4)
  sapply(overlaps, length, simplify = T) != 0
}

#function to check whether each of ints1 overlaps any of ints2
overlaps_any <- function(ints1, ints2){
  overlaps <- interval_overlap(ints1,ints2)
  sapply(overlaps, length, simplify = T) != 0
}


# function to check whether each of ints1 totally includes any of ints2 
includes_any <- function(ints1,ints2){
  inclusions <- interval_included(ints1,ints2)
  sapply(inclusions, length, simplify = T) != 0
}

# shift one row (either direction)
rowShift <- function(x, shiftLen = 1L) {
  r <- (1L + shiftLen):(length(x) + shiftLen)
  r[r<1] <- NA
  return(x[r])
}

# find start ends of blocks
findBlocks <- function(x, windowSize){
  y <- as.numeric(as.character(x$BP.wg))
  #import windowsize as KB, so make B
  windowSize <- windowSize*1000
  # subset so that the diff between row above and below
  starts <- x[((is.na(y - rowShift( y,-1)) | y - rowShift( y,-1)>windowSize)),]; starts$pos.type <- c("start")
  ends <- x[((is.na(rowShift( y,+1)- y) | rowShift( y,+1)- y>windowSize)),]; ends$pos.type <- c("end")
  arrange(rbind(starts,ends), BP.wg) # concatenate, then order by BP.wg increasing
}

# store density values from simulations
density.df <- function(x) {
  temp <- density(x)
  return(data.frame(perc = temp$x, density = temp$y))}


########### data ########### 
setwd("/Users/gabrielamontejokovacevich/Dropbox (Cambridge University)/git/2021-altitude-heliconius/")

# shdr
era.e.shdr <- read.csv("local/data/sharing/era.shdr.para.east.csv"); head(era.e.shdr)
era.w.shdr <- read.csv("local/data/sharing/era.shdr.para.west.csv"); head(era.w.shdr)
era.allo.shdr <- read.csv("local/data/sharing/era.shdr.allo.summ.csv"); head(era.allo.shdr)
mel.e.shdr <- read.csv("local/data/sharing/mel.shdr.para.east.csv"); head(mel.e.shdr)
mel.w.shdr <- read.csv("local/data/sharing/mel.shdr.para.west.csv"); head(mel.w.shdr)
mel.allo.shdr <- read.csv("local/data/sharing/mel.shdr.allo.summ.csv"); head(mel.allo.shdr)
#mel.allo.shdr <- read.csv("local/data/sharing/mel.shdr.allo.all.csv"); head(mel.allo.shdr)

ref.scaff.era <- read.table("local/data/ref/Heliconius_erato_demophoon_v1_-_scaffolds.fa.fai", row.names = NULL)
ref.scaff.mel <- read.table("local/data/ref/Hmel2.5.scaffolds.fa.fai", row.names = NULL)
scafEnds.era <- cumsum(ref.scaff.era[,2]); offset.era <- scafEnds.era - ref.scaff.era[,2]
scafEnds.mel <- cumsum(ref.scaff.mel[,2]); offset.mel <- scafEnds.mel - ref.scaff.mel[,2]

########### 2. INTERVALS OVERLAPS - OBSERVED ###########
# are there overlaps between intervals of erato and intervals of lifted over melpomene shdr
# bare in mind era.e./w. include allo hdrs, but with start ends relevant to each side
# best use maximal allo hdrs FOR THE SIMULATIONS  - so merge intervals of east west

era.e.shdr.int <- reduce(Intervals(as.matrix(era.e.shdr[,c("start","end")]), closed=c(TRUE,TRUE), type="R")) 
era.w.shdr.int <- reduce(Intervals(as.matrix(era.w.shdr[,c("start","end")]), closed=c(TRUE,TRUE), type="R")) 
era.e.w.int <- c(era.e.shdr.int, era.w.shdr.int); nrow(era.e.w.int)
era.e.w.int <-reduce(era.e.w.int); nrow(era.e.w.int)

# use lifted over postiions
mel.e.shdr.int <- reduce(Intervals(as.matrix(mel.e.shdr[,c("start.era.BP.wg","end.era.BP.wg" )]), closed=c(TRUE,TRUE), type="R")) 
mel.w.shdr.int <- reduce(Intervals(as.matrix(mel.w.shdr[,c("start.era.BP.wg","end.era.BP.wg" )]), closed=c(TRUE,TRUE), type="R")) 
mel.e.w.int <- c(mel.e.shdr.int, mel.w.shdr.int); length(mel.e.w.int)
mel.e.w.int <-reduce(mel.e.w.int); nrow(mel.e.w.int)

## check if any overlap between era and mel east/west dfs that have allopatric
era.mel.overlaps.A <- overlaps_any(era.e.w.int, mel.e.w.int); era.mel.overlaps.A
era.mel.overlaps.intA <- sum(era.mel.overlaps.A); era.mel.overlaps.intA #112
era.mel.overlaps.int.propA <-  (era.mel.overlaps.intA / (nrow(era.e.w.int)))*100; era.mel.overlaps.int.propA # 17% of hdrs overlap

era.mel.overlaps.B <- overlaps_any( mel.e.w.int, era.e.w.int); era.mel.overlaps.B
era.mel.overlaps.intB <- sum(era.mel.overlaps.B); era.mel.overlaps.intB #112
era.mel.overlaps.int.propB <-  (era.mel.overlaps.intB / (nrow(era.e.w.int)))*100; era.mel.overlaps.int.propB # 17% of hdrs overlap

(era.mel.overlaps.int.propA + era.mel.overlaps.int.propB)/2

########### 3. INTERVALS OVERLAPS - SIMULATED ###########
#### random blocks overlap sim era ####
# one sim per trio (unlike Ana, only had one sim for )
era.genome_size <- tail(scafEnds.era,1)
mel.genome_size <- tail(scafEnds.mel,1)

#genome_size <- genome_size-bp.removed.era.mean 
sim10k.era <-  list(); sim10k.era <- lapply(1:10000, function(x){rand_non_overlapping_intervals(era.genome_size,size(era.e.w.int))})
sim10k.mel <-  list(); sim10k.mel <- lapply(1:10000, function(x){rand_non_overlapping_intervals(mel.genome_size,size(mel.e.w.int))})

# do overlapsAB (both directions), counts, prop on the fly
sim10k.era.mel.overlap.prop <- (((sapply(lapply(sim10k.era, function(x){overlaps_any(x,mel.e.w.int)}), sum, simplify = T))+
                                 (sapply(lapply(sim10k.mel, function(x){overlaps_any(x,era.e.w.int)}), sum, simplify = T))) / 
                                (nrow(era.e.w.int)+nrow(mel.e.w.int)))*100 

mean(sim10k.era.mel.overlap.prop); sd(sim10k.era.mel.overlap.prop); ( era.mel.overlaps.int.propA+ era.mel.overlaps.int.propB)/2


