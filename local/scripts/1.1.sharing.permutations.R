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
setwd("git/2021-altitude-heliconius/")
ref.scaff.era <- read.table("local/data/ref/Heliconius_erato_demophoon_v1_-_scaffolds.fa.fai", row.names = NULL)
ref.scaff.mel <- read.table("local/data/ref/Hmel2.5.scaffolds.fa.fai", row.names = NULL)

# read in files and modify accordingly 
era.tests <- list()
era.listpaths <- dir(pattern = ".*era.*txt", path = "local/data/pbs.out/"); era.listpaths
for (i in 1:4) {
  era.tests[[i]] <- read.table(era.listpaths[i])
  # fix headers
  names(era.tests[[i]]) <- lapply(era.tests[[i]][1,], as.character)
  era.tests[[i]] <- era.tests[[i]][-1,]
  era.tests[[i]]$SNP <- as.character(era.tests[[i]]$SNP)
  if(names(era.tests[[i]])[8]=="PBS0"){
    # make variables that will be used numeric
    era.tests[[i]]$PBS0 <- as.numeric(as.character(era.tests[[i]]$PBS0))
    #era.tests[[i]]$PBE0 <- as.numeric(as.character(era.tests[[i]]$PBE0))
    era.tests[[i]]$zPBS0 <- zscore(era.tests[[i]]$PBS0)
  } else if(names(era.tests[[i]])[5]=="fst"){
    era.tests[[i]]$fst <- as.numeric(as.character(era.tests[[i]]$fst))
    era.tests[[i]]$zfst <- zscore(era.tests[[i]]$fst)
  }
}

mel.tests <- list()
mel.listpaths <- dir(pattern = ".*mel.*txt", path = "local/data/pbs.out/"); mel.listpaths
for (i in 1:4) {
  mel.tests[[i]] <- read.table(mel.listpaths[i])
  names(mel.tests[[i]]) <- lapply(mel.tests[[i]][1,], as.character)
  mel.tests[[i]] <- mel.tests[[i]][-1,]
  mel.tests[[i]]$SNP <- as.character(mel.tests[[i]]$SNP)
  if(names(mel.tests[[i]])[8]=="PBS0"){
    mel.tests[[i]]$PBS0 <- as.numeric(as.character(mel.tests[[i]]$PBS0))
    #mel.tests[[i]]$PBE0 <- as.numeric(as.character(mel.tests[[i]]$PBE0))
    mel.tests[[i]]$zPBS0 <- zscore(mel.tests[[i]]$PBS0)
  } else if(names(mel.tests[[i]])[5]=="fst"){
    mel.tests[[i]]$fst <- as.numeric(as.character(mel.tests[[i]]$fst))
    mel.tests[[i]]$zfst <- zscore(mel.tests[[i]]$fst)
  }
}

### subset ref.scaffs that we have in our data, and drop those we dont ###
# especially important for melpomene where many unassigned scaffolds
ref.scaff.era <- subset(ref.scaff.era, V1 %in% unique(era.tests[[1]]$scaff))
ref.scaff.mel <- subset(ref.scaff.mel, V1 %in% unique(mel.tests[[1]]$scaff))

########### data prep - top1, exp100, BP.wg ########### 
###### BP whole genome offsets ######
# create BP.wg first, use all along- STICHING scaff together.
# get scaffold offsets from fai
scafEnds.era <- cumsum(ref.scaff.era[,2]); offset.era <- scafEnds.era - ref.scaff.era[,2]
plot(scafEnds.era,offset.era); names(offset.era) <- ref.scaff.era[,1]
scafEnds.mel <- cumsum(ref.scaff.mel[,2]); offset.mel <- scafEnds.mel - ref.scaff.mel[,2]
plot(scafEnds.mel,offset.mel); names(offset.mel) <- ref.scaff.mel[,1]
# check levels present in data
levels(era.tests[[1]]$scaff); levels(mel.tests[[1]]$scaff)

#  sub offset by scaffolds that are present (done already)
offset.era.sub <- offset.era[names(offset.era) %in% unique(era.tests[[1]]$scaff)]
offset.mel.sub <- offset.mel[names(offset.mel) %in% unique(mel.tests[[1]]$scaff)]

# adjust midpos of window by adding whole genome offset, to make whole genome BP - to all tests
# before we had added scaff offset within chromosome (here no need for BP)
for (i in 1:4) {
  era.tests[[i]]$BP.wg <- as.numeric(as.character(era.tests[[i]]$midPos)) + offset.era.sub[era.tests[[i]]$scaff]}
for (i in 1:4) {
  mel.tests[[i]]$BP.wg <- as.numeric(as.character(mel.tests[[i]]$midPos)) + offset.mel.sub[mel.tests[[i]]$scaff]}

###### top 4std ######

### era zPBS0
era.tests.top1.pbs0 <- list()
for (i in seq_along(era.tests)) {
  n<-1
  if(names(era.tests[[i]])[8]=="PBS0"){
    # get sites above 3std from mean
    era.tests.top1.pbs0[[i]] <- subset(era.tests[[i]], zPBS0 > 4)
  } else if(names(era.tests[[i]])[5]=="fst"){
    # get sites above 3std from mean
    era.tests.top1.pbs0[[i]] <- subset(era.tests[[i]], zfst > 4)
  }
}

### mel zPBS0
mel.tests.top1.pbs0 <- list()
for (i in seq_along(mel.tests)) {
  n<-1
  if(names(mel.tests[[i]])[8]=="PBS0"){
    # get sites within top 1% PBE values
    mel.tests.top1.pbs0[[i]] <- subset(mel.tests[[i]], zPBS0 > 4)
  } else if(names(mel.tests[[i]])[5]=="fst"){
    # get sites within top 1% fst values
    mel.tests.top1.pbs0[[i]] <- subset(mel.tests[[i]], zfst > 4)
  }
}

###### top 4std ######
## test to compare sharing
### era zPBS1
era.tests.top1.PBS1 <- list()
for (i in seq_along(era.tests)) {
  n<-1
  if(names(era.tests[[i]])[8]=="PBS1"){
    # get sites above 3std from mean
    era.tests.top1.PBS1[[i]] <- subset(era.tests[[i]], zPBS1 > 4)
  } else if(names(era.tests[[i]])[5]=="fst"){
    # get sites above 3std from mean
    era.tests.top1.PBS1[[i]] <- subset(era.tests[[i]], zfst > 4)
  }
}

### mel zPBS1
mel.tests.top1.PBS1 <- list()
for (i in seq_along(mel.tests)) {
  n<-1
  if(names(mel.tests[[i]])[8]=="PBS1"){
    # get sites within top 1% PBE values
    mel.tests.top1.PBS1[[i]] <- subset(mel.tests[[i]], zPBS1 > 4)
  } else if(names(mel.tests[[i]])[5]=="fst"){
    # get sites within top 1% fst values
    mel.tests.top1.PBS1[[i]] <- subset(mel.tests[[i]], zfst > 4)
  }
}

###### expand BP.wg 50kb ######
### start with 1050kb, work your way down (when re-running, testing kb)
## replace expanded50for expanded50 , 50kb for 50kb , (!!!!)
# era-  buffer50 top1
list.pos.to.include <- list()
era.tests.top1.pbs0.expanded50<- list()
for (i in 1:4) {
  list.pos.to.include[[i]] <- unique((apply(era.tests.top1.pbs0[[i]], 1 , bp.list.expanded,50))[TRUE])
  era.tests.top1.pbs0.expanded50[[i]] <-subset(era.tests[[i]], BP.wg %in% list.pos.to.include[[i]] )
}


# mel- - buffer50 top1
list.pos.to.include <- list()
mel.tests.top1.pbs0.expanded50<- list()
for (i in 1:4) {
  list.pos.to.include[[i]] <- unique((apply(mel.tests.top1.pbs0[[i]], 1 , bp.list.expanded, 50))[TRUE])
  mel.tests.top1.pbs0.expanded50[[i]] <-subset(mel.tests[[i]], BP.wg %in% list.pos.to.include[[i]] )
}

# #check
# plot(era.tests.top1.pbs0.expanded50[[2]]$BP.wg, as.numeric(as.character(era.tests.top1.pbs0.expanded50[[2]]$CHR)))
# plot(mel.tests.top1.pbs0.expanded50[[2]]$BP.wg, as.numeric(as.character(mel.tests.top1.pbs0.expanded50[[2]]$CHR)))

###################### 1. FIND BLOCKS ###########
# find blocks in selected dataset, WITHIN SCAFFOLDS, with whole-genome BP start-end
# spread struggles with repeated keys, so create new key/index by row
# find 1kb blocks (windows)
blocks.era <- list()
for (i in 1:4) {
  blocks.era[[i]] <- era.tests.top1.pbs0.expanded50[[i]] %>% 
    findBlocks(1) %>% dplyr::select(scaff, BP.wg, pos.type)%>% 
    group_by_at(vars(-BP.wg)) %>%  # group by everything other than the value column. 
    mutate(row_id=1:n()) %>% ungroup() %>%  # build group index
    spread(key=pos.type, value=BP.wg) %>%    # spread
    dplyr::select(-row_id) %>%  # drop the index
    subset(end-start!=0) # remove weird block length 0
}

blocks.mel <- list()
for (i in 1:4) {
  blocks.mel[[i]] <- mel.tests.top1.pbs0.expanded50[[i]] %>% 
    findBlocks(1) %>% dplyr::select(scaff, BP.wg, pos.type)%>% 
    group_by_at(vars(-BP.wg)) %>%  # group by everything other than the value column. 
    mutate(row_id=1:n()) %>% ungroup() %>%  # build group index
    spread(key=pos.type, value=BP.wg) %>%    # spread
    dplyr::select(-row_id) %>% # drop the index
    subset(end-start!=0) # remove weird block length 0
}

## test blocks, look
ggplot(data=subset( blocks.era[[1]], scaff=="Herato0101"), aes(x=end)) +
  geom_rect(inherit.aes = F, data=subset(blocks.era[[1]], scaff=="Herato0101"), aes(xmin=start, xmax=end, ymin=3,ymax=4.5), colour="transparent", fill="#049E73", alpha=0.7 )+
  geom_rect(inherit.aes = F, data=subset(blocks.era[[2]], scaff=="Herato0101"), aes(xmin=start, xmax=end, ymin=0,ymax=1.5), colour="transparent", fill="#0372B2", alpha=0.7 )+ 
  geom_rect(inherit.aes = F, data=subset(blocks.era[[3]], scaff=="Herato0101"), aes(xmin=start, xmax=end, ymin=4.5,ymax=6), colour="transparent", fill="#049E73", alpha=0.7 )+
  geom_rect(inherit.aes = F, data=subset(blocks.era[[4]], scaff=="Herato0101"), aes(xmin=start, xmax=end, ymin=1.5,ymax=3), colour="transparent", fill="#0372B2", alpha=0.7 )+
  scale_y_continuous( expand = c(0, 0) )

ggplot(data=subset( blocks.mel[[1]], scaff=="Hmel201001o"), aes(x=end)) +
  geom_rect(inherit.aes = F, data=subset(blocks.mel[[1]], scaff=="Hmel201001o"), aes(xmin=start, xmax=end, ymin=3,ymax=4.5), colour="transparent", fill="#049E73", alpha=0.7 )+
  geom_rect(inherit.aes = F, data=subset(blocks.mel[[2]], scaff=="Hmel201001o"), aes(xmin=start, xmax=end, ymin=0,ymax=1.5), colour="transparent", fill="#0372B2", alpha=0.7 )+ 
  geom_rect(inherit.aes = F, data=subset(blocks.mel[[3]], scaff=="Hmel201001o"), aes(xmin=start, xmax=end, ymin=4.5,ymax=6), colour="transparent", fill="#049E73", alpha=0.7 )+
  geom_rect(inherit.aes = F, data=subset(blocks.mel[[4]], scaff=="Hmel201001o"), aes(xmin=start, xmax=end, ymin=1.5,ymax=3), colour="transparent", fill="#0372B2", alpha=0.7 )+
  scale_y_continuous( expand = c(0, 0) )



#### summary table blocks ####
block.summ.era <-list()
for (i in 1:4) {
  block.summ.era[[i]] <- data_frame(pop=i,
           mean.block.size.kb=(mean(blocks.era[[i]]$end- blocks.era[[i]]$start))/1000,
           n.block=length(blocks.era[[i]]$end),
           perc.genome.cov.block=(1-(sum(ref.scaff.era[,2]) - sum(blocks.era[[i]]$end- blocks.era[[i]]$start))/sum(ref.scaff.era[,2]))*100)
}

block.summ.era <- rbind(block.summ.era[[1]], block.summ.era[[2]], block.summ.era[[3]], block.summ.era[[4]])
block.summ.era$pop <- c("co.e", "co.w", "ec.e", "ec.w");block.summ.era
write.csv(block.summ.era, "local/plots/sharing.sims/pbs.out/block.summ.era.50kb.csv", row.names = FALSE)

block.summ.mel <-list()
for (i in 1:4) {
  block.summ.mel[[i]] <- data_frame(pop=i,
                                    mean.block.size.kb=(mean(blocks.mel[[i]]$end- blocks.mel[[i]]$start))/1000,
                                    n.block=length(blocks.mel[[i]]$end),
                                    perc.genome.cov.block=(1-(sum(ref.scaff.mel[,2]) - sum(blocks.mel[[i]]$end- blocks.mel[[i]]$start))/sum(ref.scaff.mel[,2]))*100)
}

block.summ.mel <- rbind(block.summ.mel[[1]], block.summ.mel[[2]], block.summ.mel[[3]], block.summ.mel[[4]])
block.summ.mel$pop <- c("co.e", "co.w", "ec.e", "ec.w");block.summ.mel
write.csv(block.summ.mel, "local/plots/sharing.sims/pbs.out/block.summ.mel.50kb.csv", row.names = FALSE)


###################### 2. SELECTED BLOCKS INTERVALS ###########
## create "selected blocks" intervals, within species
# use type R so that start and end are kept
# use liftover mel positions when across species
# scaff start end
blocks.era.int <-list()
for (i in 1:4) {
  blocks.era.int[[i]] <- reduce(Intervals(as.matrix(blocks.era[[i]][,c("start","end")]), closed=c(TRUE,TRUE), type="R")) }

blocks.mel.int <-list()
for (i in 1:4) {
  blocks.mel.int[[i]] <- reduce(Intervals(as.matrix(blocks.mel[[i]][,c("start","end")]), closed=c(TRUE,TRUE), type="R")) }

#plot(blocks.era.int[[1]][20:30])
#plot(blocks.mel.int[[1]][1:10])
Intervals(as.matrix(blocks.era[[i]][,c("start","end")]), closed=c(TRUE,TRUE), type="R")[18]
blocks.era[[i]][,c("start","end")][18,]
blocks.era.int[[i]][18,]


###################### 3. INTERVALS OVERLAPS - OBSERVED ###########
###### era ######
## overlaps
# check whether each of ints1 overlaps any of ints2- answer YES/NO- obtain indices
# era.e, in both directions, then means for percentage shared
era.e.overlaps.intA <- overlaps_any(blocks.era.int[[1]], blocks.era.int[[3]])
era.e.overlaps.int.countA <- sum(era.e.overlaps.intA); era.e.overlaps.int.countA #112
era.e.overlaps.int.propA <-  (era.e.overlaps.int.countA / (nrow(blocks.era.int[[1]])))*100 # 53% of the top, 242 blocks total
era.e.overlaps.intB <- overlaps_any(blocks.era.int[[3]], blocks.era.int[[1]])
era.e.overlaps.int.countB <- sum(era.e.overlaps.intB);era.e.overlaps.int.countB #10
era.e.overlaps.int.propB <-  (era.e.overlaps.int.countB / (nrow(blocks.era.int[[3]])))*100 # 25% of the top, 242 blocks total
era.e.overlaps.int.prop <-  ((era.e.overlaps.int.countA +era.e.overlaps.int.countB) / (nrow(blocks.era.int[[1]])+nrow(blocks.era.int[[3]])))*100 # 46% of the top, 242 blocks total

# era.w
era.w.overlaps.intA <- overlaps_any(blocks.era.int[[2]], blocks.era.int[[4]])
era.w.overlaps.int.countA <- sum(era.w.overlaps.intA) #112
era.w.overlaps.int.propA <-  (era.w.overlaps.int.countA / (nrow(blocks.era.int[[2]])))*100 # 46% of the top, 242 blocks total
era.w.overlaps.intB <- overlaps_any(blocks.era.int[[4]], blocks.era.int[[2]])
era.w.overlaps.int.countB <- sum(era.w.overlaps.intB) #10
era.w.overlaps.int.propB <-  (era.w.overlaps.int.countB / (nrow(blocks.era.int[[4]])))*100 # 46% of the top, 242 blocks total
era.w.overlaps.int.prop <-  ((era.w.overlaps.int.countA +era.w.overlaps.int.countB) / (nrow(blocks.era.int[[2]])+nrow(blocks.era.int[[4]])))*100 # 46% of the top, 242 blocks total

# sharing across sides option1 
# create new intervals for overlaps (not interesections, hence 'union') on each side of Andes
# grab unions (c, reduce) of overlapping ints within sides
# 80 intervals east
era.e.overlap.union <- interval_union(reduce(Intervals(as.matrix(subset(blocks.era.int[[1]],  era.e.overlaps.intA)),type = "R")),
                                      reduce(Intervals(as.matrix(subset(blocks.era.int[[3]],  era.e.overlaps.intB)),type = "R"))); head(era.e.overlap.union )
plot(era.e.overlap.union[1:2])

# 143 intervals east
era.w.overlap.union <- interval_union(reduce(Intervals(as.matrix(subset(blocks.era.int[[2]],  era.w.overlaps.intA)),type = "R")),
                                      reduce(Intervals(as.matrix(subset(blocks.era.int[[4]],  era.w.overlaps.intB)),type = "R")))

# check overlaps across sides
era.overlaps.intA <- overlaps_any(era.e.overlap.union , era.w.overlap.union )
era.overlaps.int.countA <- sum(era.overlaps.intA); era.overlaps.int.countA #112
era.overlaps.int.propA <-  (era.overlaps.int.countA / (nrow(era.e.overlap.union )))*100 #
era.overlaps.intB <- overlaps_any(era.w.overlap.union, era.e.overlap.union )
era.overlaps.int.countB <- sum(era.overlaps.intB);era.overlaps.int.countB #10
era.overlaps.int.propB <-  (era.overlaps.int.countB / (nrow(era.w.overlap.union)))*100 # 
era.overlaps.int.prop <-  ((era.overlaps.int.countA +era.overlaps.int.countB) / (nrow(era.e.overlap.union)+nrow(era.w.overlap.union)))*100 ; era.overlaps.int.prop # 

# 41.25% of the within sides shared blocks are shared across sides 

###### mel ######
## overlaps
# check whether each of ints1 overlaps any of ints2
# mel.e, in both directions, then means for percentage shared
mel.e.overlaps.intA <- overlaps_any(blocks.mel.int[[1]], blocks.mel.int[[3]])
mel.e.overlaps.int.countA <- sum(mel.e.overlaps.intA); mel.e.overlaps.int.countA #112
mel.e.overlaps.int.propA <-  (mel.e.overlaps.int.countA / (nrow(blocks.mel.int[[1]])))*100 # 53% of the top, 242 blocks total
mel.e.overlaps.intB <- overlaps_any(blocks.mel.int[[3]], blocks.mel.int[[1]])
mel.e.overlaps.int.countB <- sum(mel.e.overlaps.intB);mel.e.overlaps.int.countB #10
mel.e.overlaps.int.propB <-  (mel.e.overlaps.int.countB / (nrow(blocks.mel.int[[3]])))*100 # 25% of the top, 242 blocks total
mel.e.overlaps.int.prop <-  ((mel.e.overlaps.int.countA +mel.e.overlaps.int.countB) / (nrow(blocks.mel.int[[1]])+nrow(blocks.mel.int[[3]])))*100 # 46% of the top, 242 blocks total

# mel.w
mel.w.overlaps.intA <- overlaps_any(blocks.mel.int[[2]], blocks.mel.int[[4]])
mel.w.overlaps.int.countA <- sum(mel.w.overlaps.intA) #112
mel.w.overlaps.int.propA <-  (mel.w.overlaps.int.countA / (nrow(blocks.mel.int[[2]])))*100 # 46% of the top, 242 blocks total
mel.w.overlaps.intB <- overlaps_any(blocks.mel.int[[4]], blocks.mel.int[[2]])
mel.w.overlaps.int.countB <- sum(mel.w.overlaps.intB) #10
mel.w.overlaps.int.propB <-  (mel.w.overlaps.int.countB / (nrow(blocks.mel.int[[4]])))*100 # 46% of the top, 242 blocks total
mel.w.overlaps.int.prop <-  ((mel.w.overlaps.int.countA +mel.w.overlaps.int.countB) / (nrow(blocks.mel.int[[2]])+nrow(blocks.mel.int[[4]])))*100 # 46% of the top, 242 blocks total


# create new intervals for overlaps on each side of Andes
# grab unions (c, reduce) of overlapping ints within sides
# 18 intervals east
mel.e.overlap.union <- interval_union(reduce(Intervals(as.matrix(subset(blocks.mel.int[[1]],  mel.e.overlaps.intA)),type = "R")),
                                      reduce(Intervals(as.matrix(subset(blocks.mel.int[[3]],  mel.e.overlaps.intB)),type = "R"))); head(mel.e.overlap.union )
plot(reduce(Intervals(as.matrix(subset(blocks.mel.int[[1]],  mel.e.overlaps.intA)),type = "R"))[1])
plot(mel.e.overlap.union[1:2])

# 104 intervals east
mel.w.overlap.union <- interval_union(reduce(Intervals(as.matrix(subset(blocks.mel.int[[2]],  mel.w.overlaps.intA)),type = "R")),
                                      reduce(Intervals(as.matrix(subset(blocks.mel.int[[4]],  mel.w.overlaps.intB)),type = "R")))

mel.overlaps.intA <- overlaps_any(mel.e.overlap.union , mel.w.overlap.union )
mel.overlaps.int.countA <- sum(mel.overlaps.intA); mel.overlaps.int.countA #112
mel.overlaps.int.propA <-  (mel.overlaps.int.countA / (nrow(mel.e.overlap.union )))*100 
mel.overlaps.intB <- overlaps_any(mel.w.overlap.union, mel.e.overlap.union )
mel.overlaps.int.countB <- sum(mel.overlaps.intB);mel.overlaps.int.countB #10
mel.overlaps.int.propB <-  (mel.overlaps.int.countB / (nrow(mel.w.overlap.union)))*100 
mel.overlaps.int.prop <-  ((mel.overlaps.int.countA +mel.overlaps.int.countB) / (nrow(mel.e.overlap.union)+nrow(mel.w.overlap.union)))*100 

# creates new interval with intersections across all pops- NO
mel.inters.int <- interval_intersection(blocks.mel.int[[1]], blocks.mel.int[[3]], blocks.mel.int[[2]], blocks.mel.int[[4]])
plot(mel.inters.int[1])

###################### 4. ADD TO ALL.PBS for analyses/plotting ######
era.all.pbs <- read.csv("local/data/pbs.out/era.all.pbs.out.csv"); head(era.all.pbs)
mel.all.pbs <- read.csv("local/data/pbs.out/mel.all.pbs.out.csv"); head(mel.all.pbs)

###### era all hdrs per cline ######
# create list of private hdrs, with start ends names
blocks.era.int.df <-list()
for (i in 1:4) {blocks.era.int.df[[i]] <- as.data.frame(blocks.era.int[[i]])
  names(blocks.era.int.df[[i]]) <- c("start", "end")
  blocks.era.int.df[[i]]$block.no <- sprintf('%0.3d', 1:nrow(blocks.era.int.df[[i]])) }; head(blocks.era.int.df[[i]])

# add cline to each dataset, then concatenate with number
blocks.era.int.df[[1]]$cline <- "CoE"; blocks.era.int.df[[2]]$cline <- "CoW"; blocks.era.int.df[[3]]$cline <- "EcE"; blocks.era.int.df[[4]]$cline <- "EcW"
for (i in 1:4) {blocks.era.int.df[[i]]$hdr.name <- paste(blocks.era.int.df[[i]]$cline, blocks.era.int.df[[i]]$block.no, sep = ".") } ; head(blocks.era.int.df[[i]])

# add private hdr names matching start ends
blocks.era.int.df.all <-  rbind(blocks.era.int.df[[1]], blocks.era.int.df[[2]],  blocks.era.int.df[[3]],  blocks.era.int.df[[4]]); head(blocks.era.int.df.all)
blocks.era.int.df.all.wide <-  spread(blocks.era.int.df.all, cline, hdr.name); head(blocks.era.int.df.all.wide)

# make long so that start ends in one BP.wg column (and both get matched)
blocks.era.int.df.all.wide.longbp <- gather(blocks.era.int.df.all.wide, type, BP.wg, -c(block.no, CoE ,  CoW , EcE  , EcW)); head(blocks.era.int.df.all.wide.longbp)

# add start or ends bp.wg
# prevent NAs from trying to replace values, as some blocks may have same start/ends across different clines
era.all.pbs$hdr.CoE <-  blocks.era.int.df.all.wide.longbp[!is.na(blocks.era.int.df.all.wide.longbp$CoE),]$CoE[
  match(era.all.pbs$BP.wg, blocks.era.int.df.all.wide.longbp[!is.na(blocks.era.int.df.all.wide.longbp$CoE),]$BP.wg)] 
era.all.pbs$hdr.CoW <-  blocks.era.int.df.all.wide.longbp[!is.na(blocks.era.int.df.all.wide.longbp$CoW),]$CoW[
  match(era.all.pbs$BP.wg, blocks.era.int.df.all.wide.longbp[!is.na(blocks.era.int.df.all.wide.longbp$CoW),]$BP.wg)] 
era.all.pbs$hdr.EcE <-  blocks.era.int.df.all.wide.longbp[!is.na(blocks.era.int.df.all.wide.longbp$EcE),]$EcE[
  match(era.all.pbs$BP.wg, blocks.era.int.df.all.wide.longbp[!is.na(blocks.era.int.df.all.wide.longbp$EcE),]$BP.wg)] 
era.all.pbs$hdr.EcW <-  blocks.era.int.df.all.wide.longbp[!is.na(blocks.era.int.df.all.wide.longbp$EcW),]$EcW[
  match(era.all.pbs$BP.wg, blocks.era.int.df.all.wide.longbp[!is.na(blocks.era.int.df.all.wide.longbp$EcW),]$BP.wg)] 

# check that starts and ends of each private hdr are there
test <- subset(era.all.pbs, hdr.CoE!=""| hdr.CoW!="" | hdr.EcW!=""| hdr.EcE!="")

# fill in hdrs with their names, but make sure to only use the subset of the positions within the relevant (eg CoE) block ranges
# use data table inrange https://stackoverflow.com/questions/43641874/subset-by-multiple-ranges
# then fill from last value upwards  with zoo (fill between start and end)
library(zoo)

era.all.pbs[ era.all.pbs$BP.wg %in% as.data.frame(as.data.table(era.all.pbs)[era.all.pbs$BP.wg %inrange% subset(blocks.era.int.df.all, cline=="CoE" )[,1:2]])$BP.wg ,]$hdr.CoE <- 
  c(na.locf(as.data.table(era.all.pbs)[era.all.pbs$BP.wg %inrange% subset(blocks.era.int.df.all, cline=="CoE" )[,1:2]]$hdr.CoE, fromLast = TRUE))
era.all.pbs[ era.all.pbs$BP.wg %in% as.data.frame(as.data.table(era.all.pbs)[era.all.pbs$BP.wg %inrange% subset(blocks.era.int.df.all, cline=="CoW" )[,1:2]])$BP.wg ,]$hdr.CoW <- 
  c(na.locf(as.data.table(era.all.pbs)[era.all.pbs$BP.wg %inrange% subset(blocks.era.int.df.all, cline=="CoW" )[,1:2]]$hdr.CoW, fromLast = TRUE))
era.all.pbs[ era.all.pbs$BP.wg %in% as.data.frame(as.data.table(era.all.pbs)[era.all.pbs$BP.wg %inrange% subset(blocks.era.int.df.all, cline=="EcE" )[,1:2]])$BP.wg ,]$hdr.EcE <- 
  c(na.locf(as.data.table(era.all.pbs)[era.all.pbs$BP.wg %inrange% subset(blocks.era.int.df.all, cline=="EcE" )[,1:2]]$hdr.EcE, fromLast = TRUE))
era.all.pbs[ era.all.pbs$BP.wg %in% as.data.frame(as.data.table(era.all.pbs)[era.all.pbs$BP.wg %inrange% subset(blocks.era.int.df.all, cline=="EcW" )[,1:2]])$BP.wg ,]$hdr.EcW <- 
  c(na.locf(as.data.table(era.all.pbs)[era.all.pbs$BP.wg %inrange% subset(blocks.era.int.df.all, cline=="EcW" )[,1:2]]$hdr.EcW, fromLast = TRUE))

test <- subset(era.all.pbs, hdr.CoE!=""| hdr.CoW!="" | hdr.EcW!=""| hdr.EcE!="")


###### era parapatric Shdrs ######
# now we have all pop-specific HDRs in the master sheet. 
# add overlaps within sides of Andes, parapatric, use the Intervals output (as a check)
# subset intervals of CoE (1) that overlap with EcE OR viceversa, then add hdr CoE.no_EcE.no to column in dataset

# ranges of overlapping blocks (reduced)
era.e.overlap.union.df <- data.frame(era.e.overlap.union); names(era.e.overlap.union.df) <- c("start", "end"); head(era.e.overlap.union.df); nrow(era.e.overlap.union.df)
era.w.overlap.union.df <- data.frame(era.w.overlap.union); names(era.w.overlap.union.df) <- c("start", "end"); head(era.w.overlap.union.df); nrow(era.w.overlap.union.df)

# give them new para.east sharing IDs
era.e.overlap.union.df$shdr.para.east.id <- paste("shdr.east", sprintf('%0.3d', 1:nrow(era.e.overlap.union.df)), sep="."); head(era.e.overlap.union.df)
era.w.overlap.union.df$shdr.para.west.id <- paste("shdr.west", sprintf('%0.3d', 1:nrow(era.w.overlap.union.df)), sep="."); head(era.w.overlap.union.df)

# add new shdr names to pbs.all- add the starts of each region, then fill in
era.all.pbs$shdr.para.east.id <-NA
era.all.pbs[era.all.pbs$BP.wg %in%  as.data.frame(as.data.table(era.all.pbs)[era.all.pbs$BP.wg %inrange% era.e.overlap.union.df[,1:2]])$BP.wg ,]$shdr.para.east.id <- 
  era.e.overlap.union.df$shdr.para.east.id[match(era.all.pbs[era.all.pbs$BP.wg %in%  as.data.frame(as.data.table(era.all.pbs)[era.all.pbs$BP.wg %inrange% era.e.overlap.union.df[,1:2]])$BP.wg ,]$BP.wg , era.e.overlap.union.df[,1])  ]  
# fill
era.all.pbs[era.all.pbs$BP.wg %in%  as.data.frame(as.data.table(era.all.pbs)[era.all.pbs$BP.wg %inrange% era.e.overlap.union.df[,1:2]])$BP.wg ,]$shdr.para.east.id <-
  c(na.locf(era.all.pbs[era.all.pbs$BP.wg %in%  as.data.frame(as.data.table(era.all.pbs)[era.all.pbs$BP.wg %inrange% era.e.overlap.union.df[,1:2]])$BP.wg ,]$shdr.para.east.id, fromLast = F))

era.all.pbs$shdr.para.west.id <-NA
era.all.pbs[era.all.pbs$BP.wg %in%  as.data.frame(as.data.table(era.all.pbs)[era.all.pbs$BP.wg %inrange% era.w.overlap.union.df[,1:2]])$BP.wg ,]$shdr.para.west.id <- 
  era.w.overlap.union.df$shdr.para.west.id[match(era.all.pbs[era.all.pbs$BP.wg %in%  as.data.frame(as.data.table(era.all.pbs)[era.all.pbs$BP.wg %inrange% era.w.overlap.union.df[,1:2]])$BP.wg ,]$BP.wg , era.w.overlap.union.df[,1])  ]  
# fill
era.all.pbs[era.all.pbs$BP.wg %in%  as.data.frame(as.data.table(era.all.pbs)[era.all.pbs$BP.wg %inrange% era.w.overlap.union.df[,1:2]])$BP.wg ,]$shdr.para.west.id <-
  c(na.locf(era.all.pbs[era.all.pbs$BP.wg %in%  as.data.frame(as.data.table(era.all.pbs)[era.all.pbs$BP.wg %inrange% era.w.overlap.union.df[,1:2]])$BP.wg ,]$shdr.para.west.id, fromLast = F))


# obtain shdr concatenated hdr names (only use subset of data within overlaps)
era.all.pbs$shdr.para.east.id.concat <-NA
era.all.pbs[era.all.pbs$BP.wg %in%  as.data.frame(as.data.table(era.all.pbs)[era.all.pbs$BP.wg %inrange% era.e.overlap.union.df[,1:2]])$BP.wg ,]$shdr.para.east.id.concat <-
  paste(era.all.pbs[era.all.pbs$BP.wg %in%  as.data.frame(as.data.table(era.all.pbs)[era.all.pbs$BP.wg %inrange% era.e.overlap.union.df[,1:2]])$BP.wg ,]$hdr.CoE,
        era.all.pbs[era.all.pbs$BP.wg %in%  as.data.frame(as.data.table(era.all.pbs)[era.all.pbs$BP.wg %inrange% era.e.overlap.union.df[,1:2]])$BP.wg ,]$hdr.EcE, sep = "_")
era.all.pbs$shdr.para.west.id.concat <-NA
era.all.pbs[era.all.pbs$BP.wg %in%  as.data.frame(as.data.table(era.all.pbs)[era.all.pbs$BP.wg %inrange% era.w.overlap.union.df[,1:2]])$BP.wg ,]$shdr.para.west.id.concat <-
  paste(era.all.pbs[era.all.pbs$BP.wg %in%  as.data.frame(as.data.table(era.all.pbs)[era.all.pbs$BP.wg %inrange% era.w.overlap.union.df[,1:2]])$BP.wg ,]$hdr.CoW,
        era.all.pbs[era.all.pbs$BP.wg %in%  as.data.frame(as.data.table(era.all.pbs)[era.all.pbs$BP.wg %inrange% era.w.overlap.union.df[,1:2]])$BP.wg ,]$hdr.EcW, sep = "_")

test <- subset(era.all.pbs, shdr.para.east.id!="")
  
# is it an overlap or an intersect, add for each side
era.all.pbs$shdr.para.east.overlap.intersect <- if_else((substr(era.all.pbs$shdr.para.east.id.concat,1, 3)=="NA_" | substr(era.all.pbs$shdr.para.east.id.concat,8,10) =="_NA"), "overlap", "intersect")
era.all.pbs$shdr.para.west.overlap.intersect <- if_else((substr(era.all.pbs$shdr.para.west.id.concat,1, 3)=="NA_" | substr(era.all.pbs$shdr.para.west.id.concat,8,10) =="_NA"), "overlap", "intersect")
summarise(group_by(era.all.pbs,hdr.para.east.overlap.intersect,shdr.para.east.id.concat, hdr.CoE, hdr.EcE, hdr.para.west.overlap.intersect,shdr.para.west.id.concat, hdr.CoW, hdr.EcW ))

test <- subset(era.all.pbs, hdr.para.east!="")

###### era allopatric Shdrs ######

# grab overlaps of unions in both directions, then reduce
era.overlap.union.df <- data.frame(interval_union(reduce(Intervals(as.matrix(subset(era.e.overlap.union,  era.overlaps.intA)),type = "R")),
                                      reduce(Intervals(as.matrix(subset(era.w.overlap.union,  era.overlaps.intB)),type = "R")))); names(era.overlap.union.df) <- c("start", "end"); head(era.overlap.union.df ); nrow(era.overlap.union.df)

# give them new allo.east sharing IDs
era.overlap.union.df$shdr.allo.id <- paste("shdr.all", sprintf('%0.3d', 1:nrow(era.overlap.union.df)), sep="."); head(era.overlap.union.df)

# add new shdr names to pbs.all- add the starts of each region, then fill in
era.all.pbs$shdr.allo.id <-NA
era.all.pbs[era.all.pbs$BP.wg %in%  as.data.frame(as.data.table(era.all.pbs)[era.all.pbs$BP.wg %inrange% era.overlap.union.df[,1:2]])$BP.wg ,]$shdr.allo.id <- 
  era.overlap.union.df$shdr.allo.id[match(era.all.pbs[era.all.pbs$BP.wg %in%  as.data.frame(as.data.table(era.all.pbs)[era.all.pbs$BP.wg %inrange% era.overlap.union.df[,1:2]])$BP.wg ,]$BP.wg , era.overlap.union.df[,1])  ]  
# fill
era.all.pbs[era.all.pbs$BP.wg %in%  as.data.frame(as.data.table(era.all.pbs)[era.all.pbs$BP.wg %inrange% era.overlap.union.df[,1:2]])$BP.wg ,]$shdr.allo.id <-
  c(na.locf(era.all.pbs[era.all.pbs$BP.wg %in%  as.data.frame(as.data.table(era.all.pbs)[era.all.pbs$BP.wg %inrange% era.overlap.union.df[,1:2]])$BP.wg ,]$shdr.allo.id, fromLast = F))


# obtain shdr concatenated hdr names (only use subset of data within overlaps)
era.all.pbs$shdr.allo.id.concat <-NA
era.all.pbs[era.all.pbs$BP.wg %in%  as.data.frame(as.data.table(era.all.pbs)[era.all.pbs$BP.wg %inrange% era.overlap.union.df[,1:2]])$BP.wg ,]$shdr.allo.id.concat <-
  paste(era.all.pbs[era.all.pbs$BP.wg %in%  as.data.frame(as.data.table(era.all.pbs)[era.all.pbs$BP.wg %inrange% era.overlap.union.df[,1:2]])$BP.wg ,]$hdr.CoE,
        era.all.pbs[era.all.pbs$BP.wg %in%  as.data.frame(as.data.table(era.all.pbs)[era.all.pbs$BP.wg %inrange% era.overlap.union.df[,1:2]])$BP.wg ,]$hdr.EcE,
        era.all.pbs[era.all.pbs$BP.wg %in%  as.data.frame(as.data.table(era.all.pbs)[era.all.pbs$BP.wg %inrange% era.overlap.union.df[,1:2]])$BP.wg ,]$hdr.CoW,
        era.all.pbs[era.all.pbs$BP.wg %in%  as.data.frame(as.data.table(era.all.pbs)[era.all.pbs$BP.wg %inrange% era.overlap.union.df[,1:2]])$BP.wg ,]$hdr.EcW, sep = "_")

# is it an overlap or an intersect, if it's not overlapping all of the original per cline HDRs then it's an overlap, otherwise an intersection (they all had it)
era.all.pbs$shdr.allo.overlap.intersect <- NA
era.all.pbs[!(is.na(era.all.pbs$shdr.allo.id)),]$shdr.allo.overlap.intersect <- if_else(grepl("NA", subset(era.all.pbs, !(is.na(shdr.allo.id)))$shdr.allo.id.concat, fixed = TRUE)==T, "overlap", "intersect")

# add sharing type for plotting
era.all.pbs$sharing.type <- if_else(!(is.na(era.all.pbs$shdr.allo.id)), "across", if_else(!(is.na(era.all.pbs$shdr.para.east.id)), "within.east", if_else(!(is.na(era.all.pbs$shdr.para.west.id)), "within.west", "NA" )))

###### store info of hdrs making up allo shdr ######
unique(era.all.pbs$sharing.type)
era.shdr.allo.summ <- summarise(group_by(subset(era.all.pbs,sharing.type=="across") ,scaff, sharing.type,shdr.allo.id),
                            shdr.para.west = paste(unique( na.omit(shdr.para.west.id))),
                            shdr.para.east = unique( na.omit(shdr.para.east.id )),
                            start=min(midPos), end=max(midPos ),
                            start.BP.wg=min(BP.wg), end.BP.wg=max(BP.wg ));era.shdr.allo.summ
write.csv(era.shdr.allo.summ, "local/data/sharing/era.shdr.allo.summ.csv", row.names = F)

####### era old naming shdr by concatenating #####
era.all.pbs$hdr.para.east <-NA
era.all.pbs[era.all.pbs$BP.wg %in% 
              as.data.frame(as.data.table(era.all.pbs)[
                (era.all.pbs$BP.wg %inrange% data.frame(blocks.era.int[[1]][
                  overlaps_any(blocks.era.int[[1]], blocks.era.int[[3]])]) )|(era.all.pbs$BP.wg %inrange% data.frame(blocks.era.int[[3]][overlaps_any(blocks.era.int[[3]], blocks.era.int[[1]])]) ) ])
            $BP.wg ,]$hdr.para.east <- c(paste(era.all.pbs[era.all.pbs$BP.wg %in% 
                                                             as.data.frame(as.data.table(era.all.pbs)[
                                                               (era.all.pbs$BP.wg %inrange% data.frame(blocks.era.int[[1]][
                                                                 overlaps_any(blocks.era.int[[1]], blocks.era.int[[3]])]) )|(era.all.pbs$BP.wg %inrange% data.frame(blocks.era.int[[3]][overlaps_any(blocks.era.int[[3]], blocks.era.int[[1]])]) ) ])
                                                           $BP.wg ,]$hdr.CoE, 
                                               era.all.pbs[era.all.pbs$BP.wg %in% 
                                                             as.data.frame(as.data.table(era.all.pbs)[
                                                               (era.all.pbs$BP.wg %inrange% data.frame(blocks.era.int[[1]][
                                                                 overlaps_any(blocks.era.int[[1]], blocks.era.int[[3]])]) )|(era.all.pbs$BP.wg %inrange% data.frame(blocks.era.int[[3]][overlaps_any(blocks.era.int[[3]], blocks.era.int[[1]])]) ) ])
                                                           $BP.wg ,]$hdr.EcE, sep = "_"))
###### mel all hdrs per cline ######
# create list of private hdrs, with start ends names
blocks.mel.int.df <-list()
for (i in 1:4) {blocks.mel.int.df[[i]] <- as.data.frame(blocks.mel.int[[i]])
names(blocks.mel.int.df[[i]]) <- c("start", "end")
blocks.mel.int.df[[i]]$block.no <- sprintf('%0.3d', 1:nrow(blocks.mel.int.df[[i]])) }; head(blocks.mel.int.df[[i]])

# add cline to each dataset, then concatenate with number
blocks.mel.int.df[[1]]$cline <- "CoE"; blocks.mel.int.df[[2]]$cline <- "CoW"; blocks.mel.int.df[[3]]$cline <- "EcE"; blocks.mel.int.df[[4]]$cline <- "EcW"
for (i in 1:4) {blocks.mel.int.df[[i]]$hdr.name <- paste(blocks.mel.int.df[[i]]$cline, blocks.mel.int.df[[i]]$block.no, sep = ".") } ; head(blocks.mel.int.df[[i]])

# add private hdr names matching start ends
blocks.mel.int.df.all <-  rbind(blocks.mel.int.df[[1]], blocks.mel.int.df[[2]],  blocks.mel.int.df[[3]],  blocks.mel.int.df[[4]]); head(blocks.mel.int.df.all)
blocks.mel.int.df.all.wide <-  spread(blocks.mel.int.df.all, cline, hdr.name); head(blocks.mel.int.df.all.wide)

# make long so that start ends in one BP.wg column (and both get matched)
blocks.mel.int.df.all.wide.longbp <- gather(blocks.mel.int.df.all.wide, type, BP.wg, -c(block.no, CoE ,  CoW , EcE  , EcW)); head(blocks.mel.int.df.all.wide.longbp)

# add start or ends bp.wg
# prevent NAs from trying to replace values, as some blocks may have same start/ends across different clines
mel.all.pbs$hdr.CoE <-  blocks.mel.int.df.all.wide.longbp[!is.na(blocks.mel.int.df.all.wide.longbp$CoE),]$CoE[
  match(mel.all.pbs$BP.wg, blocks.mel.int.df.all.wide.longbp[!is.na(blocks.mel.int.df.all.wide.longbp$CoE),]$BP.wg)] 
mel.all.pbs$hdr.CoW <-  blocks.mel.int.df.all.wide.longbp[!is.na(blocks.mel.int.df.all.wide.longbp$CoW),]$CoW[
  match(mel.all.pbs$BP.wg, blocks.mel.int.df.all.wide.longbp[!is.na(blocks.mel.int.df.all.wide.longbp$CoW),]$BP.wg)] 
mel.all.pbs$hdr.EcE <-  blocks.mel.int.df.all.wide.longbp[!is.na(blocks.mel.int.df.all.wide.longbp$EcE),]$EcE[
  match(mel.all.pbs$BP.wg, blocks.mel.int.df.all.wide.longbp[!is.na(blocks.mel.int.df.all.wide.longbp$EcE),]$BP.wg)] 
mel.all.pbs$hdr.EcW <-  blocks.mel.int.df.all.wide.longbp[!is.na(blocks.mel.int.df.all.wide.longbp$EcW),]$EcW[
  match(mel.all.pbs$BP.wg, blocks.mel.int.df.all.wide.longbp[!is.na(blocks.mel.int.df.all.wide.longbp$EcW),]$BP.wg)] 

# check that starts and ends of each private hdr are there
test <- subset(mel.all.pbs, hdr.CoE!=""| hdr.CoW!="" | hdr.EcW!=""| hdr.EcE!="")

# fill in hdrs with their names, but make sure to only use the subset of the positions within the relevant (eg CoE) block ranges
# use data table inrange https://stackoverflow.com/questions/43641874/subset-by-multiple-ranges
# then fill from last value upwards  with zoo (fill between start and end)
library(zoo)

mel.all.pbs[ mel.all.pbs$BP.wg %in% as.data.frame(as.data.table(mel.all.pbs)[mel.all.pbs$BP.wg %inrange% subset(blocks.mel.int.df.all, cline=="CoE" )[,1:2]])$BP.wg ,]$hdr.CoE <- 
  c(na.locf(as.data.table(mel.all.pbs)[mel.all.pbs$BP.wg %inrange% subset(blocks.mel.int.df.all, cline=="CoE" )[,1:2]]$hdr.CoE, fromLast = TRUE))
mel.all.pbs[ mel.all.pbs$BP.wg %in% as.data.frame(as.data.table(mel.all.pbs)[mel.all.pbs$BP.wg %inrange% subset(blocks.mel.int.df.all, cline=="CoW" )[,1:2]])$BP.wg ,]$hdr.CoW <- 
  c(na.locf(as.data.table(mel.all.pbs)[mel.all.pbs$BP.wg %inrange% subset(blocks.mel.int.df.all, cline=="CoW" )[,1:2]]$hdr.CoW, fromLast = TRUE))
mel.all.pbs[ mel.all.pbs$BP.wg %in% as.data.frame(as.data.table(mel.all.pbs)[mel.all.pbs$BP.wg %inrange% subset(blocks.mel.int.df.all, cline=="EcE" )[,1:2]])$BP.wg ,]$hdr.EcE <- 
  c(na.locf(as.data.table(mel.all.pbs)[mel.all.pbs$BP.wg %inrange% subset(blocks.mel.int.df.all, cline=="EcE" )[,1:2]]$hdr.EcE, fromLast = TRUE))
mel.all.pbs[ mel.all.pbs$BP.wg %in% as.data.frame(as.data.table(mel.all.pbs)[mel.all.pbs$BP.wg %inrange% subset(blocks.mel.int.df.all, cline=="EcW" )[,1:2]])$BP.wg ,]$hdr.EcW <- 
  c(na.locf(as.data.table(mel.all.pbs)[mel.all.pbs$BP.wg %inrange% subset(blocks.mel.int.df.all, cline=="EcW" )[,1:2]]$hdr.EcW, fromLast = TRUE))

test <- subset(mel.all.pbs, hdr.CoE!=""| hdr.CoW!="" | hdr.EcW!=""| hdr.EcE!="")


###### mel parapatric Shdrs ######
# now we have all pop-specific HDRs in the master sheet. 
# add overlaps within sides of Andes, parapatric, use the Intervals output (as a check)
# subset intervals of CoE (1) that overlap with EcE OR viceversa, then add hdr CoE.no_EcE.no to column in dataset

# ranges of overlapping blocks (reduced)
mel.e.overlap.union.df <- data.frame(mel.e.overlap.union); names(mel.e.overlap.union.df) <- c("start", "end"); head(mel.e.overlap.union.df); nrow(mel.e.overlap.union.df)
mel.w.overlap.union.df <- data.frame(mel.w.overlap.union); names(mel.w.overlap.union.df) <- c("start", "end"); head(mel.w.overlap.union.df); nrow(mel.w.overlap.union.df)

# give them new para.east sharing IDs
mel.e.overlap.union.df$shdr.para.east.id <- paste("shdr.east", sprintf('%0.3d', 1:nrow(mel.e.overlap.union.df)), sep="."); head(mel.e.overlap.union.df)
mel.w.overlap.union.df$shdr.para.west.id <- paste("shdr.west", sprintf('%0.3d', 1:nrow(mel.w.overlap.union.df)), sep="."); head(mel.w.overlap.union.df)

# add new shdr names to pbs.all- add the starts of each region, then fill in
mel.all.pbs$shdr.para.east.id <-NA
mel.all.pbs[mel.all.pbs$BP.wg %in%  as.data.frame(as.data.table(mel.all.pbs)[mel.all.pbs$BP.wg %inrange% mel.e.overlap.union.df[,1:2]])$BP.wg ,]$shdr.para.east.id <- 
  mel.e.overlap.union.df$shdr.para.east.id[match(mel.all.pbs[mel.all.pbs$BP.wg %in%  as.data.frame(as.data.table(mel.all.pbs)[mel.all.pbs$BP.wg %inrange% mel.e.overlap.union.df[,1:2]])$BP.wg ,]$BP.wg , mel.e.overlap.union.df[,1])  ]  
# fill
mel.all.pbs[mel.all.pbs$BP.wg %in%  as.data.frame(as.data.table(mel.all.pbs)[mel.all.pbs$BP.wg %inrange% mel.e.overlap.union.df[,1:2]])$BP.wg ,]$shdr.para.east.id <-
  c(na.locf(mel.all.pbs[mel.all.pbs$BP.wg %in%  as.data.frame(as.data.table(mel.all.pbs)[mel.all.pbs$BP.wg %inrange% mel.e.overlap.union.df[,1:2]])$BP.wg ,]$shdr.para.east.id, fromLast = F))

mel.all.pbs$shdr.para.west.id <-NA
mel.all.pbs[mel.all.pbs$BP.wg %in%  as.data.frame(as.data.table(mel.all.pbs)[mel.all.pbs$BP.wg %inrange% mel.w.overlap.union.df[,1:2]])$BP.wg ,]$shdr.para.west.id <- 
  mel.w.overlap.union.df$shdr.para.west.id[match(mel.all.pbs[mel.all.pbs$BP.wg %in%  as.data.frame(as.data.table(mel.all.pbs)[mel.all.pbs$BP.wg %inrange% mel.w.overlap.union.df[,1:2]])$BP.wg ,]$BP.wg , mel.w.overlap.union.df[,1])  ]  
# fill
mel.all.pbs[mel.all.pbs$BP.wg %in%  as.data.frame(as.data.table(mel.all.pbs)[mel.all.pbs$BP.wg %inrange% mel.w.overlap.union.df[,1:2]])$BP.wg ,]$shdr.para.west.id <-
  c(na.locf(mel.all.pbs[mel.all.pbs$BP.wg %in%  as.data.frame(as.data.table(mel.all.pbs)[mel.all.pbs$BP.wg %inrange% mel.w.overlap.union.df[,1:2]])$BP.wg ,]$shdr.para.west.id, fromLast = F))


# obtain shdr concatenated hdr names (only use subset of data within overlaps)
mel.all.pbs$shdr.para.east.id.concat <-NA
mel.all.pbs[mel.all.pbs$BP.wg %in%  as.data.frame(as.data.table(mel.all.pbs)[mel.all.pbs$BP.wg %inrange% mel.e.overlap.union.df[,1:2]])$BP.wg ,]$shdr.para.east.id.concat <-
  paste(mel.all.pbs[mel.all.pbs$BP.wg %in%  as.data.frame(as.data.table(mel.all.pbs)[mel.all.pbs$BP.wg %inrange% mel.e.overlap.union.df[,1:2]])$BP.wg ,]$hdr.CoE,
        mel.all.pbs[mel.all.pbs$BP.wg %in%  as.data.frame(as.data.table(mel.all.pbs)[mel.all.pbs$BP.wg %inrange% mel.e.overlap.union.df[,1:2]])$BP.wg ,]$hdr.EcE, sep = "_")
mel.all.pbs$shdr.para.west.id.concat <-NA
mel.all.pbs[mel.all.pbs$BP.wg %in%  as.data.frame(as.data.table(mel.all.pbs)[mel.all.pbs$BP.wg %inrange% mel.w.overlap.union.df[,1:2]])$BP.wg ,]$shdr.para.west.id.concat <-
  paste(mel.all.pbs[mel.all.pbs$BP.wg %in%  as.data.frame(as.data.table(mel.all.pbs)[mel.all.pbs$BP.wg %inrange% mel.w.overlap.union.df[,1:2]])$BP.wg ,]$hdr.CoW,
        mel.all.pbs[mel.all.pbs$BP.wg %in%  as.data.frame(as.data.table(mel.all.pbs)[mel.all.pbs$BP.wg %inrange% mel.w.overlap.union.df[,1:2]])$BP.wg ,]$hdr.EcW, sep = "_")

test <- subset(mel.all.pbs, shdr.para.east.id!="")

# is it an overlap or an intersect, add for each side
mel.all.pbs$shdr.para.east.overlap.intersect <- if_else((substr(mel.all.pbs$shdr.para.east.id.concat,1, 3)=="NA_" | substr(mel.all.pbs$shdr.para.east.id.concat,8,10) =="_NA"), "overlap", "intersect")
mel.all.pbs$shdr.para.west.overlap.intersect <- if_else((substr(mel.all.pbs$shdr.para.west.id.concat,1, 3)=="NA_" | substr(mel.all.pbs$shdr.para.west.id.concat,8,10) =="_NA"), "overlap", "intersect")
summarise(group_by(mel.all.pbs,shdr.para.east.overlap.intersect,shdr.para.east.id.concat, hdr.CoE, hdr.EcE, shdr.para.west.overlap.intersect, shdr.para.west.id.concat, hdr.CoW, hdr.EcW ))

###### mel allopatric Shdrs ######

# grab overlaps of unions in both directions, then reduce
mel.overlap.union.df <- data.frame(interval_union(reduce(Intervals(as.matrix(subset(mel.e.overlap.union,  mel.overlaps.intA)),type = "R")),
                                                  reduce(Intervals(as.matrix(subset(mel.w.overlap.union,  mel.overlaps.intB)),type = "R")))); names(mel.overlap.union.df) <- c("start", "end"); head(mel.overlap.union.df ); nrow(mel.overlap.union.df)

# give them new allo.east sharing IDs
mel.overlap.union.df$shdr.allo.id <- paste("shdr.all", sprintf('%0.3d', 1:nrow(mel.overlap.union.df)), sep="."); head(mel.overlap.union.df)

# add new shdr names to pbs.all- add the starts of each region, then fill in
mel.all.pbs$shdr.allo.id <-NA
mel.all.pbs[mel.all.pbs$BP.wg %in%  as.data.frame(as.data.table(mel.all.pbs)[mel.all.pbs$BP.wg %inrange% mel.overlap.union.df[,1:2]])$BP.wg ,]$shdr.allo.id <- 
  mel.overlap.union.df$shdr.allo.id[match(mel.all.pbs[mel.all.pbs$BP.wg %in%  as.data.frame(as.data.table(mel.all.pbs)[mel.all.pbs$BP.wg %inrange% mel.overlap.union.df[,1:2]])$BP.wg ,]$BP.wg , mel.overlap.union.df[,1])  ]  
# fill
mel.all.pbs[mel.all.pbs$BP.wg %in%  as.data.frame(as.data.table(mel.all.pbs)[mel.all.pbs$BP.wg %inrange% mel.overlap.union.df[,1:2]])$BP.wg ,]$shdr.allo.id <-
  c(na.locf(mel.all.pbs[mel.all.pbs$BP.wg %in%  as.data.frame(as.data.table(mel.all.pbs)[mel.all.pbs$BP.wg %inrange% mel.overlap.union.df[,1:2]])$BP.wg ,]$shdr.allo.id, fromLast = F))


# obtain shdr concatenated hdr names (only use subset of data within overlaps)
mel.all.pbs$shdr.allo.id.concat <-NA
mel.all.pbs[mel.all.pbs$BP.wg %in%  as.data.frame(as.data.table(mel.all.pbs)[mel.all.pbs$BP.wg %inrange% mel.overlap.union.df[,1:2]])$BP.wg ,]$shdr.allo.id.concat <-
  paste(mel.all.pbs[mel.all.pbs$BP.wg %in%  as.data.frame(as.data.table(mel.all.pbs)[mel.all.pbs$BP.wg %inrange% mel.overlap.union.df[,1:2]])$BP.wg ,]$hdr.CoE,
        mel.all.pbs[mel.all.pbs$BP.wg %in%  as.data.frame(as.data.table(mel.all.pbs)[mel.all.pbs$BP.wg %inrange% mel.overlap.union.df[,1:2]])$BP.wg ,]$hdr.EcE,
        mel.all.pbs[mel.all.pbs$BP.wg %in%  as.data.frame(as.data.table(mel.all.pbs)[mel.all.pbs$BP.wg %inrange% mel.overlap.union.df[,1:2]])$BP.wg ,]$hdr.CoW,
        mel.all.pbs[mel.all.pbs$BP.wg %in%  as.data.frame(as.data.table(mel.all.pbs)[mel.all.pbs$BP.wg %inrange% mel.overlap.union.df[,1:2]])$BP.wg ,]$hdr.EcW, sep = "_")

# is it an overlap or an intersect, if it's not overlapping all of the original per cline HDRs then it's an overlap, otherwise an intersection (they all had it)
mel.all.pbs$shdr.allo.overlap.intersect <- if_else(grepl("NA", mel.all.pbs$shdr.allo.id.concat, fixed = TRUE)==T, "overlap", "intersect")
test <- subset(mel.all.pbs, shdr.allo.id!="")

# add sharing type for plotting
mel.all.pbs$sharing.type <- if_else(!(is.na(mel.all.pbs$shdr.allo.id)), "across", if_else(!(is.na(mel.all.pbs$shdr.para.east.id)), "within.east", if_else(!(is.na(mel.all.pbs$shdr.para.west.id)), "within.west", "NA" )))


###### store info of hdrs making up allo shdr ######
mel.shdr.allo.summ <- summarise(group_by(subset(mel.all.pbs,sharing.type=="across") ,scaff, sharing.type,shdr.allo.id),
                                shdr.para.west = paste(unique( na.omit(shdr.para.west.id))),
                                shdr.para.east = unique( na.omit(shdr.para.east.id )),
                                start=min(midPos), end=max(midPos ),
                                start.BP.wg=min(BP.wg), end.BP.wg=max(BP.wg ));mel.shdr.allo.summ
write.csv(mel.shdr.allo.summ, "local/data/sharing/mel.shdr.allo.summ.csv", row.names = F)

#### store all.pbs ######
head(era.all.pbs); head(mel.all.pbs)
write.csv(era.all.pbs, "local/data/pbs.out/era.all.pbs.out.csv", row.names = F)
write.csv(mel.all.pbs, "local/data/pbs.out/mel.all.pbs.out.csv", row.names = F)

#### store hdr ranges? 1.all per cline, 2. shdr.para ,  3.shdr.allo ####
write.csv(blocks.era.int.df.all, "local/data/sharing/era.cline.hdr.csv", row.names = F)
write.csv(era.e.overlap.union.df, "local/data/sharing/era.shdr.para.east.csv", row.names = F)
write.csv(era.w.overlap.union.df, "local/data/sharing/era.shdr.para.west.csv", row.names = F)
write.csv(era.overlap.union.df, "local/data/sharing/era.shdr.allo.all.csv", row.names = F)

write.csv(blocks.mel.int.df.all, "local/data/sharing/mel.cline.hdr.csv", row.names = F)
write.csv(mel.e.overlap.union.df, "local/data/sharing/mel.shdr.para.east.csv", row.names = F)
write.csv(mel.w.overlap.union.df, "local/data/sharing/mel.shdr.para.west.csv", row.names = F)
write.csv(mel.overlap.union.df, "local/data/sharing/mel.shdr.allo.all.csv", row.names = F)


era.hdr.summ <- data_frame(sharing.type=c("hdr.CoE", "hdr.CoW","hdr.EcE","hdr.EcW","shdr.para.east", "shdr.para.west", "shdr.allo.all"),
           no.blocks=c(nrow(blocks.era.int.df[[1]]), nrow(blocks.era.int.df[[2]]), nrow(blocks.era.int.df[[3]]), nrow(blocks.era.int.df[[4]]), nrow(era.e.overlap.union.df), nrow(era.w.overlap.union.df), nrow(era.overlap.union.df) ),
           block.length.mean= c(mean(blocks.era.int.df[[1]]$end - blocks.era.int.df[[1]]$start), mean(blocks.era.int.df[[2]]$end - blocks.era.int.df[[2]]$start), mean(blocks.era.int.df[[3]]$end - blocks.era.int.df[[3]]$start), mean(blocks.era.int.df[[4]]$end - blocks.era.int.df[[4]]$start), 
                                mean(era.e.overlap.union.df$end-era.e.overlap.union.df$start), mean(era.w.overlap.union.df$end-era.w.overlap.union.df$start), mean(era.overlap.union.df$end - era.overlap.union.df$start) ),
           block.length.median= c(median(blocks.era.int.df[[1]]$end - blocks.era.int.df[[1]]$start), median(blocks.era.int.df[[2]]$end - blocks.era.int.df[[2]]$start), median(blocks.era.int.df[[3]]$end - blocks.era.int.df[[3]]$start), median(blocks.era.int.df[[4]]$end - blocks.era.int.df[[4]]$start), 
                                median(era.e.overlap.union.df$end-era.e.overlap.union.df$start), median(era.w.overlap.union.df$end-era.w.overlap.union.df$start), median(era.overlap.union.df$end - era.overlap.union.df$start) ),
           block.length.sd= c(sd(blocks.era.int.df[[1]]$end - blocks.era.int.df[[1]]$start), sd(blocks.era.int.df[[2]]$end - blocks.era.int.df[[2]]$start), sd(blocks.era.int.df[[3]]$end - blocks.era.int.df[[3]]$start), sd(blocks.era.int.df[[4]]$end - blocks.era.int.df[[4]]$start), 
                                sd(era.e.overlap.union.df$end-era.e.overlap.union.df$start), sd(era.w.overlap.union.df$end-era.w.overlap.union.df$start), sd(era.overlap.union.df$end - era.overlap.union.df$start) ) ); era.hdr.summ

mel.hdr.summ <- data_frame(sharing.type=c("hdr.CoE", "hdr.CoW","hdr.EcE","hdr.EcW","shdr.para.east", "shdr.para.west", "shdr.allo.all"),
                           no.blocks=c(nrow(blocks.mel.int.df[[1]]), nrow(blocks.mel.int.df[[2]]), nrow(blocks.mel.int.df[[3]]), nrow(blocks.mel.int.df[[4]]), nrow(mel.e.overlap.union.df), nrow(mel.w.overlap.union.df), nrow(mel.overlap.union.df) ),
                           block.length.mean= c(mean(blocks.mel.int.df[[1]]$end - blocks.mel.int.df[[1]]$start), mean(blocks.mel.int.df[[2]]$end - blocks.mel.int.df[[2]]$start), mean(blocks.mel.int.df[[3]]$end - blocks.mel.int.df[[3]]$start), mean(blocks.mel.int.df[[4]]$end - blocks.mel.int.df[[4]]$start), 
                                                mean(mel.e.overlap.union.df$end-mel.e.overlap.union.df$start), mean(mel.w.overlap.union.df$end-mel.w.overlap.union.df$start), mean(mel.overlap.union.df$end - mel.overlap.union.df$start) ),
                           block.length.median= c(median(blocks.mel.int.df[[1]]$end - blocks.mel.int.df[[1]]$start), median(blocks.mel.int.df[[2]]$end - blocks.mel.int.df[[2]]$start), median(blocks.mel.int.df[[3]]$end - blocks.mel.int.df[[3]]$start), median(blocks.mel.int.df[[4]]$end - blocks.mel.int.df[[4]]$start), 
                                                  median(mel.e.overlap.union.df$end-mel.e.overlap.union.df$start), median(mel.w.overlap.union.df$end-mel.w.overlap.union.df$start), median(mel.overlap.union.df$end - mel.overlap.union.df$start) ),
                           block.length.sd= c(sd(blocks.mel.int.df[[1]]$end - blocks.mel.int.df[[1]]$start), sd(blocks.mel.int.df[[2]]$end - blocks.mel.int.df[[2]]$start), sd(blocks.mel.int.df[[3]]$end - blocks.mel.int.df[[3]]$start), sd(blocks.mel.int.df[[4]]$end - blocks.mel.int.df[[4]]$start), 
                                              sd(mel.e.overlap.union.df$end-mel.e.overlap.union.df$start), sd(mel.w.overlap.union.df$end-mel.w.overlap.union.df$start), sd(mel.overlap.union.df$end - mel.overlap.union.df$start) ) ); mel.hdr.summ

write.csv(era.hdr.summ, "local/data/sharing/era.hdr.summ.csv", row.names = F)
write.csv(mel.hdr.summ, "local/data/sharing/mel.hdr.summ.csv", row.names = F)


### create updated era.summ.regions.shared
head(era.all.pbs)
summarise(group_by(era.all.pbs, shdr.allo.id,shdr.para.east.id, shdr.para.west.id ),
          shdr.para.west.id.a=c(unique(shdr.para.west.id )),
          scaff=unique(scaff),
          start.BP.wg=min(BP.wg), eng.BP.wg=max(BP.wg))

## summ hdr numbers
summarise(group_by(era.all.pbs),
          hdr.CoE.no = length(unique( hdr.CoE)),
          hdr.EcE.no = length(unique( hdr.EcE)),
          hdr.CoW.no = length(unique( hdr.CoW)),
          hdr.EcW.no = length(unique( hdr.EcW)),
          shdr.para.east.no = length(unique( shdr.para.east.id )),
          shdr.para.west.no = length(unique( shdr.para.west.id )),
          shdr.allo.no= length(unique( shdr.allo.id )))

summarise(group_by(mel.all.pbs),
          hdr.CoE.no = length(unique( hdr.CoE)),
          hdr.EcE.no = length(unique( hdr.EcE)),
          hdr.CoW.no = length(unique( hdr.CoW)),
          hdr.EcW.no = length(unique( hdr.EcW)),
          shdr.para.east.no = length(unique( shdr.para.east.id )),
          shdr.para.west.no = length(unique( shdr.para.west.id )),
          shdr.allo.no= length(unique( shdr.allo.id )))

########### 4. INTERVALS OVERLAPS - SIMULATED ###########
#### random blocks overlap sim era ####
# one sim per trio (unlike Ana, only had one sim for )
genome_size <- tail(scafEnds.era,1)
#genome_size <- genome_size-bp.removed.era.mean 
sim10k.era <-  list()
for (i in 1:4) {
  sim10k.era[[i]] <- lapply(1:10000, function(x){rand_non_overlapping_intervals(genome_size,size(blocks.era.int[[i]]))}) }

# overlap between sim10k (for each trio) and the trio we want to compare to (as with the observed), 
# must have same closing (z/r) to be able to check overlaps across intervals
sim10k.era.e.overlap.propA <- ((sapply(lapply(sim10k.era[[1]], function(x){overlaps_any(x,blocks.era.int[[3]])}), sum, simplify = T))  / 
                                (nrow(blocks.era.int[[1]])))*100 
sim10k.era.e.overlap.propB <- ((sapply(lapply(sim10k.era[[3]], function(x){overlaps_any(x,blocks.era.int[[1]])}), sum, simplify = T))  / 
                                 (nrow(blocks.era.int[[3]])))*100 

sim10k.era.w.overlap.propA <- ((sapply(lapply(sim10k.era[[2]], function(x){overlaps_any(x,blocks.era.int[[4]])}), sum, simplify = T))  / 
                                 (nrow(blocks.era.int[[2]])))*100 
sim10k.era.w.overlap.propB <- ((sapply(lapply(sim10k.era[[4]], function(x){overlaps_any(x,blocks.era.int[[2]])}), sum, simplify = T))  / 
                                 (nrow(blocks.era.int[[4]])))*100 

# do overlapsAB (both directions), counts, prop on the fly
sim10k.era.e.overlap.prop <- (((sapply(lapply(sim10k.era[[1]], function(x){overlaps_any(x,blocks.era.int[[3]])}), sum, simplify = T))+
                                            (sapply(lapply(sim10k.era[[3]], function(x){overlaps_any(x,blocks.era.int[[1]])}), sum, simplify = T))) / 
                                (nrow(blocks.era.int[[1]])+nrow(blocks.era.int[[3]])))*100 

sim10k.era.w.overlap.prop <- (((sapply(lapply(sim10k.era[[2]], function(x){overlaps_any(x,blocks.era.int[[4]])}), sum, simplify = T))+
                                 (sapply(lapply(sim10k.era[[4]], function(x){overlaps_any(x,blocks.era.int[[2]])}), sum, simplify = T))) / 
                                (nrow(blocks.era.int[[2]])+nrow(blocks.era.int[[4]])))*100 

### simulations for era sharing across sides, only two backgrounds
# use overlaps unions for number and size of intervals
sim10k.era.e.across <- lapply(1:10000, function(x){rand_non_overlapping_intervals(genome_size,size(era.e.overlap.union))})
sim10k.era.w.across <- lapply(1:10000, function(x){rand_non_overlapping_intervals(genome_size,size(era.w.overlap.union))})

# how many overlaps in the simulations
# the union intervals needs to be Z type
sim10k.era.overlap.across.propA <- ((sapply(lapply(sim10k.era.e.across, function(x){overlaps_any(x,era.w.overlap.union)}), sum, simplify = T))  / 
                                 (nrow(era.e.overlap.union)))*100 

sim10k.era.overlap.across.propB <- ((sapply(lapply(sim10k.era.w.across, function(x){overlaps_any(x,era.e.overlap.union)}), sum, simplify = T))  / 
                                 (nrow(era.w.overlap.union)))*100 


sim10k.era.overlap.across.prop <- (((sapply(lapply(sim10k.era.e.across, function(x){overlaps_any(x,era.w.overlap.union)}), sum, simplify = T))+
                                 (sapply(lapply(sim10k.era.w.across, function(x){overlaps_any(x,era.e.overlap.union)}), sum, simplify = T))) / 
                                (nrow(era.e.overlap.union)+nrow(era.w.overlap.union)))*100 


#### random blocks overlap sim mel ####
# one sim per trio (unlike Ana, only had one sim for )
genome_size <- tail(scafEnds.mel,1)
#genome_size <- genome_size-bp.removed.era.mean ### CHANGE TO MEL
sim10k.mel <-  list()
for (i in 1:4) {
  sim10k.mel[[i]] <- lapply(1:10000, function(x){rand_non_overlapping_intervals(genome_size,size(blocks.mel.int[[i]]))}) }

# overlap between sim10k (for each trio) and the trio we want to compare to (as with the observed), 
sim10k.mel.e.overlap.propA <- ((sapply(lapply(sim10k.mel[[1]], function(x){overlaps_any(x,blocks.mel.int[[3]])}), sum, simplify = T))  / 
                                 (nrow(blocks.mel.int[[1]])))*100 
sim10k.mel.e.overlap.propB <- ((sapply(lapply(sim10k.mel[[3]], function(x){overlaps_any(x,blocks.mel.int[[1]])}), sum, simplify = T))  / 
                                 (nrow(blocks.mel.int[[3]])))*100 

sim10k.mel.w.overlap.propA <- ((sapply(lapply(sim10k.mel[[2]], function(x){overlaps_any(x,blocks.mel.int[[4]])}), sum, simplify = T))  / 
                                 (nrow(blocks.mel.int[[2]])))*100 
sim10k.mel.w.overlap.propB <- ((sapply(lapply(sim10k.mel[[4]], function(x){overlaps_any(x,blocks.mel.int[[2]])}), sum, simplify = T))  / 
                                 (nrow(blocks.mel.int[[4]])))*100 


# do overlapsAB (both directions), counts, prop on the fly
sim10k.mel.e.overlap.prop <- (((sapply(lapply(sim10k.mel[[1]], function(x){overlaps_any(x,blocks.mel.int[[3]])}), sum, simplify = T))+
                                 (sapply(lapply(sim10k.mel[[3]], function(x){overlaps_any(x,blocks.mel.int[[1]])}), sum, simplify = T))) / 
                                (nrow(blocks.mel.int[[1]])+nrow(blocks.mel.int[[3]])))*100 

sim10k.mel.w.overlap.prop <- (((sapply(lapply(sim10k.mel[[2]], function(x){overlaps_any(x,blocks.mel.int[[4]])}), sum, simplify = T))+
                                 (sapply(lapply(sim10k.mel[[4]], function(x){overlaps_any(x,blocks.mel.int[[2]])}), sum, simplify = T))) / 
                                (nrow(blocks.mel.int[[2]])+nrow(blocks.mel.int[[4]])))*100 


### simulations for mel sharing across sides, only two backgrounds
# use overlaps unions for number and size of intervals
sim10k.mel.e.across <- lapply(1:10000, function(x){rand_non_overlapping_intervals(genome_size,size(mel.e.overlap.union))})
sim10k.mel.w.across <- lapply(1:10000, function(x){rand_non_overlapping_intervals(genome_size,size(mel.w.overlap.union))})

# how many overlaps in the simulations
# the union intervals needs to be Z type
sim10k.mel.overlap.across.propA <- ((sapply(lapply(sim10k.mel.e.across, function(x){overlaps_any(x,mel.w.overlap.union)}), sum, simplify = T))  / 
                                      (nrow(mel.e.overlap.union)))*100 

sim10k.mel.overlap.across.propB <- ((sapply(lapply(sim10k.mel.w.across, function(x){overlaps_any(x,mel.e.overlap.union)}), sum, simplify = T))  / 
                                      (nrow(mel.w.overlap.union)))*100 


sim10k.mel.overlap.across.prop <- (((sapply(lapply(sim10k.mel.e.across, function(x){overlaps_any(x,mel.w.overlap.union)}), sum, simplify = T))+
                                      (sapply(lapply(sim10k.mel.w.across, function(x){overlaps_any(x,mel.e.overlap.union)}), sum, simplify = T))) / 
                                     (nrow(mel.e.overlap.union)+nrow(mel.w.overlap.union)))*100 



########### 5. prop overlap block jacknife #########
### era ####
### random jacknife blocks across genome
genome_size <- tail(scafEnds.era,1)
#genome_size <- genome_size-bp.removed.era.mean 
block_size <- 1000000
block_starts <- seq(1,genome_size, block_size)
n_blocks <- length(block_starts)
block_Intervals <- Intervals(as.matrix(cbind(block_starts,block_starts+block_size-1)), closed = c(TRUE,TRUE), type="R")

#jackknifed trio intervals (remove 1 blcok in each case)
era.int.jackknifed <- list()
for (i in 1:4) {
  era.int.jackknifed[[i]] <- lapply(1:n_blocks, function(x){interval_difference(blocks.era.int[[i]], block_Intervals[x,])}) }

#overlap proportions for each jackknife replicate, do both ways as with above (aggregate overlaps)
# era.e
era.e.overlap.counts.jackA <- unlist(lapply(1:n_blocks, function(x){sum(overlaps_any(era.int.jackknifed[[1]][[x]],blocks.era.int[[3]]))}))
era.e.overlaps.int.prop.jackknifedA <-  (era.e.overlap.counts.jackA / (nrow(blocks.era.int[[1]])))*100 
era.e.overlap.counts.jackB <- unlist(lapply(1:n_blocks, function(x){sum(overlaps_any(era.int.jackknifed[[3]][[x]],blocks.era.int[[1]]))}))
era.e.overlaps.int.prop.jackknifedB <-  (era.e.overlap.counts.jackB / (nrow(blocks.era.int[[3]])))*100 
era.e.overlaps.int.prop.jackknifed <-  ((era.e.overlap.counts.jackA +era.e.overlap.counts.jackB) / (nrow(blocks.era.int[[1]])+nrow(blocks.era.int[[3]])))*100 

era.e.overlap.prop.pseudoA <- era.e.overlaps.int.propA*n_blocks - era.e.overlaps.int.prop.jackknifedA*(n_blocks-1)
era.e.overlap.prop.sdA <- sd(era.e.overlap.prop.pseudoA)
era.e.overlap.prop.pseudo.std_errA <- sqrt(var(era.e.overlap.prop.pseudoA)/length(era.e.overlap.prop.pseudoA))
era.e.overlap.prop.pseudoB <- era.e.overlaps.int.propB*n_blocks - era.e.overlaps.int.prop.jackknifedB*(n_blocks-1)
era.e.overlap.prop.sdB <- sd(era.e.overlap.prop.pseudoB)
era.e.overlap.prop.pseudo.std_errB <- sqrt(var(era.e.overlap.prop.pseudoB)/length(era.e.overlap.prop.pseudoB))

era.e.overlap.prop.pseudo <- era.e.overlaps.int.prop*n_blocks - era.e.overlaps.int.prop.jackknifed*(n_blocks-1)
era.e.overlap.prop.sd <- sd(era.e.overlap.prop.pseudo)
era.e.overlap.prop.pseudo.std_err <- sqrt(var(era.e.overlap.prop.pseudo)/length(era.e.overlap.prop.pseudo))

# era.w
era.w.overlap.counts.jackA <- unlist(lapply(1:n_blocks, function(x){sum(overlaps_any(era.int.jackknifed[[2]][[x]],blocks.era.int[[4]]))}))
era.w.overlaps.int.prop.jackknifedA <-  (era.w.overlap.counts.jackA / (nrow(blocks.era.int[[2]])))*100 
era.w.overlap.counts.jackB <- unlist(lapply(1:n_blocks, function(x){sum(overlaps_any(era.int.jackknifed[[4]][[x]],blocks.era.int[[2]]))}))
era.w.overlaps.int.prop.jackknifedB <-  (era.w.overlap.counts.jackB / (nrow(blocks.era.int[[4]])))*100 
era.w.overlaps.int.prop.jackknifed <-  ((era.w.overlap.counts.jackA +era.w.overlap.counts.jackB) / (nrow(blocks.era.int[[2]])+nrow(blocks.era.int[[4]])))*100 

era.w.overlap.prop.pseudoA <- era.w.overlaps.int.propA*n_blocks - era.w.overlaps.int.prop.jackknifedA*(n_blocks-1)
era.w.overlap.prop.sdA <- sd(era.w.overlap.prop.pseudoA)
era.w.overlap.prop.pseudo.std_errA <- sqrt(var(era.w.overlap.prop.pseudoA)/length(era.w.overlap.prop.pseudoA))
era.w.overlap.prop.pseudoB <- era.w.overlaps.int.propB*n_blocks - era.w.overlaps.int.prop.jackknifedB*(n_blocks-1)
era.w.overlap.prop.sdB <- sd(era.w.overlap.prop.pseudoB)
era.w.overlap.prop.pseudo.std_errB <- sqrt(var(era.w.overlap.prop.pseudoB)/length(era.w.overlap.prop.pseudoB))


era.w.overlap.prop.pseudo <- era.w.overlaps.int.prop*n_blocks - era.w.overlaps.int.prop.jackknifed*(n_blocks-1)
era.w.overlap.prop.sd <- sd(era.w.overlap.prop.pseudo)
era.w.overlap.prop.pseudo.std_err <- sqrt(var(era.w.overlap.prop.pseudo)/length(era.w.overlap.prop.pseudo))

## across sides
#jackknifed trio intervals (remove 1 blcok in each case)
era.int.jackknifed.e.across <- lapply(1:n_blocks, function(x){interval_difference(era.e.overlap.union, block_Intervals[x,])})
era.int.jackknifed.w.across <- lapply(1:n_blocks, function(x){interval_difference(era.w.overlap.union, block_Intervals[x,])})

era.overlap.counts.jackA <- unlist(lapply(1:n_blocks, function(x){sum(overlaps_any(era.int.jackknifed.e.across[[x]],era.w.overlap.union))}))
era.overlaps.int.prop.jackknifedA <-  (era.overlap.counts.jackA / (nrow(era.e.overlap.union)))*100 
era.overlap.counts.jackB <- unlist(lapply(1:n_blocks, function(x){sum(overlaps_any(era.int.jackknifed.w.across[[x]],era.e.overlap.union))}))
era.overlaps.int.prop.jackknifedB <-  (era.overlap.counts.jackB / (nrow(era.w.overlap.union)))*100 
era.overlaps.int.prop.jackknifed <-  ((era.overlap.counts.jackA +era.overlap.counts.jackB) / (nrow(era.e.overlap.union)+nrow(era.w.overlap.union)))*100 

era.overlap.prop.pseudoA <- era.overlaps.int.propA*n_blocks - era.overlaps.int.prop.jackknifedA*(n_blocks-1)
era.overlap.prop.sdA <- sd(era.overlap.prop.pseudoA)
era.overlap.prop.pseudo.std_errA <- sqrt(var(era.overlap.prop.pseudoA)/length(era.overlap.prop.pseudoA))
era.overlap.prop.pseudoB <- era.overlaps.int.propB*n_blocks - era.overlaps.int.prop.jackknifedB*(n_blocks-1)
era.overlap.prop.sdB <- sd(era.overlap.prop.pseudoB)
era.overlap.prop.pseudo.std_errB <- sqrt(var(era.overlap.prop.pseudoB)/length(era.overlap.prop.pseudoB))

era.overlap.prop.pseudo <- era.overlaps.int.prop*n_blocks - era.overlaps.int.prop.jackknifed*(n_blocks-1)
era.overlap.prop.sd <- sd(era.overlap.prop.pseudo)
era.overlap.prop.pseudo.std_err <- sqrt(var(era.overlap.prop.pseudo)/length(era.overlap.prop.pseudo))


### mel ######
### random jacknife blocks across genome
genome_size <- tail(scafEnds.mel,1)
block_size <- 1000000
block_starts <- seq(1,genome_size, block_size)
n_blocks <- length(block_starts)
block_Intervals <- Intervals(as.matrix(cbind(block_starts,block_starts+block_size-1)), closed = c(TRUE,TRUE), type="R")

#jackknifed trio intervals (remove 1 blcok in each case)
mel.int.jackknifed <- list()
for (i in 1:4) {
  mel.int.jackknifed[[i]] <- lapply(1:n_blocks, function(x){interval_difference(blocks.mel.int[[i]], block_Intervals[x,])}) }

#overlap proportions for each jackknife replicate, do both ways as with above (aggregate overlaps)
# mel.e
mel.e.overlap.counts.jackA <- unlist(lapply(1:n_blocks, function(x){sum(overlaps_any(mel.int.jackknifed[[1]][[x]],blocks.mel.int[[3]]))}))
mel.e.overlaps.int.prop.jackknifedA <-  (mel.e.overlap.counts.jackA / (nrow(blocks.mel.int[[1]])))*100 
mel.e.overlap.counts.jackB <- unlist(lapply(1:n_blocks, function(x){sum(overlaps_any(mel.int.jackknifed[[3]][[x]],blocks.mel.int[[1]]))}))
mel.e.overlaps.int.prop.jackknifedB <-  (mel.e.overlap.counts.jackB / (nrow(blocks.mel.int[[3]])))*100 
mel.e.overlaps.int.prop.jackknifed <-  ((mel.e.overlap.counts.jackA +mel.e.overlap.counts.jackB) / (nrow(blocks.mel.int[[1]])+nrow(blocks.mel.int[[3]])))*100 

mel.e.overlap.prop.pseudoA <- mel.e.overlaps.int.propA*n_blocks - mel.e.overlaps.int.prop.jackknifedA*(n_blocks-1)
mel.e.overlap.prop.sdA <- sd(mel.e.overlap.prop.pseudoA)
mel.e.overlap.prop.pseudo.std_errA <- sqrt(var(mel.e.overlap.prop.pseudoA)/length(mel.e.overlap.prop.pseudoA))
mel.e.overlap.prop.pseudoB <- mel.e.overlaps.int.propB*n_blocks - mel.e.overlaps.int.prop.jackknifedB*(n_blocks-1)
mel.e.overlap.prop.sdB <- sd(mel.e.overlap.prop.pseudoB)
mel.e.overlap.prop.pseudo.std_errB <- sqrt(var(mel.e.overlap.prop.pseudoB)/length(mel.e.overlap.prop.pseudoB))

mel.e.overlap.prop.pseudo <- mel.e.overlaps.int.prop*n_blocks - mel.e.overlaps.int.prop.jackknifed*(n_blocks-1)
mel.e.overlap.prop.sd <- sd(mel.e.overlap.prop.pseudo)
mel.e.overlap.prop.pseudo.std_err <- sqrt(var(mel.e.overlap.prop.pseudo)/length(mel.e.overlap.prop.pseudo))

# mel.w
mel.w.overlap.counts.jackA <- unlist(lapply(1:n_blocks, function(x){sum(overlaps_any(mel.int.jackknifed[[2]][[x]],blocks.mel.int[[4]]))}))
mel.w.overlaps.int.prop.jackknifedA <-  (mel.w.overlap.counts.jackA / (nrow(blocks.mel.int[[2]])))*100 
mel.w.overlap.counts.jackB <- unlist(lapply(1:n_blocks, function(x){sum(overlaps_any(mel.int.jackknifed[[4]][[x]],blocks.mel.int[[2]]))}))
mel.w.overlaps.int.prop.jackknifedB <-  (mel.w.overlap.counts.jackB / (nrow(blocks.mel.int[[4]])))*100 
mel.w.overlaps.int.prop.jackknifed <-  ((mel.w.overlap.counts.jackA +mel.w.overlap.counts.jackB) / (nrow(blocks.mel.int[[2]])+nrow(blocks.mel.int[[4]])))*100 

mel.w.overlap.prop.pseudoA <- mel.w.overlaps.int.propA*n_blocks - mel.w.overlaps.int.prop.jackknifedA*(n_blocks-1)
mel.w.overlap.prop.sdA <- sd(mel.w.overlap.prop.pseudoA)
mel.w.overlap.prop.pseudo.std_errA <- sqrt(var(mel.w.overlap.prop.pseudoA)/length(mel.w.overlap.prop.pseudoA))
mel.w.overlap.prop.pseudoB <- mel.w.overlaps.int.propB*n_blocks - mel.w.overlaps.int.prop.jackknifedB*(n_blocks-1)
mel.w.overlap.prop.sdB <- sd(mel.w.overlap.prop.pseudoB)
mel.w.overlap.prop.pseudo.std_errB <- sqrt(var(mel.w.overlap.prop.pseudoB)/length(mel.w.overlap.prop.pseudoB))

mel.w.overlap.prop.pseudo <- mel.w.overlaps.int.prop*n_blocks - mel.w.overlaps.int.prop.jackknifed*(n_blocks-1)
mel.w.overlap.prop.sd <- sd(mel.w.overlap.prop.pseudo)
mel.w.overlap.prop.pseudo.std_err <- sqrt(var(mel.w.overlap.prop.pseudo)/length(mel.w.overlap.prop.pseudo))

## across sides
#jackknifed trio intervals (remove 1 blcok in each case)
mel.int.jackknifed.e.across <- lapply(1:n_blocks, function(x){interval_difference(mel.e.overlap.union, block_Intervals[x,])})
mel.int.jackknifed.w.across <- lapply(1:n_blocks, function(x){interval_difference(mel.w.overlap.union, block_Intervals[x,])})

mel.overlap.counts.jackA <- unlist(lapply(1:n_blocks, function(x){sum(overlaps_any(mel.int.jackknifed.e.across[[x]],mel.w.overlap.union))}))
mel.overlaps.int.prop.jackknifedA <-  (mel.overlap.counts.jackA / (nrow(mel.e.overlap.union)))*100 
mel.overlap.counts.jackB <- unlist(lapply(1:n_blocks, function(x){sum(overlaps_any(mel.int.jackknifed.w.across[[x]],mel.e.overlap.union))}))
mel.overlaps.int.prop.jackknifedB <-  (mel.overlap.counts.jackB / (nrow(mel.w.overlap.union)))*100 
mel.overlaps.int.prop.jackknifed <-  ((mel.overlap.counts.jackA +mel.overlap.counts.jackB) / (nrow(mel.e.overlap.union)+nrow(mel.w.overlap.union)))*100 


mel.overlap.prop.pseudoA <- mel.overlaps.int.propA*n_blocks - mel.overlaps.int.prop.jackknifedA*(n_blocks-1)
mel.overlap.prop.sdA <- sd(mel.overlap.prop.pseudoA)
mel.overlap.prop.pseudo.std_errA <- sqrt(var(mel.overlap.prop.pseudoA)/length(mel.overlap.prop.pseudoA))
mel.overlap.prop.pseudoB <- mel.overlaps.int.propB*n_blocks - mel.overlaps.int.prop.jackknifedB*(n_blocks-1)
mel.overlap.prop.sdB <- sd(mel.overlap.prop.pseudoB)
mel.overlap.prop.pseudo.std_errB <- sqrt(var(mel.overlap.prop.pseudoB)/length(mel.overlap.prop.pseudoB))

mel.overlap.prop.pseudo <- mel.overlaps.int.prop*n_blocks - mel.overlaps.int.prop.jackknifed*(n_blocks-1)
mel.overlap.prop.sd <- sd(mel.overlap.prop.pseudo)
mel.overlap.prop.pseudo.std_err <- sqrt(var(mel.overlap.prop.pseudo)/length(mel.overlap.prop.pseudo))



########### 5. PLOT ###########

#### save densities and data for plotting fig2  ####
# only save those that will be plotted, i.e.  for each species: 2 within sides, 1 across sides
# use custom function to store density values
# densities
d1 <- density.df(sim10k.era.e.overlap.prop)
d2 <- density.df(sim10k.era.w.overlap.prop)
d3 <- density.df(sim10k.era.overlap.across.prop)

d4 <- density.df(sim10k.mel.e.overlap.prop)
d5 <- density.df(sim10k.mel.w.overlap.prop)
d6 <- density.df(sim10k.mel.overlap.across.prop)


density.10ksim.50kb <- rbind(d1,d2,d3,d4,d5,d6)
density.10ksim.50kb$species <- c(rep("era", 512*3), rep("mel", 512*3))
density.10ksim.50kb$type.comp <- c(rep("within", 512), rep("within", 512),rep("across", 512),
                                    rep("within", 512), rep("within", 512),rep("across", 512))
density.10ksim.50kb$type.side <- c(rep("e", 512), rep("w", 512),rep("across", 512),
                                    rep("e", 512), rep("w", 512),rep("across", 512))
write.csv(density.10ksim.50kb, "local/data/local/data/data/sharing/pbs.out/density.10ksim.50kb.csv", row.names = FALSE)

# vlines
sharing.percs<-NA
sharing.percs.50kb <- data.frame(overlap.int.prop=c(era.e.overlaps.int.prop,era.w.overlaps.int.prop, era.overlaps.int.prop, 
                                               mel.e.overlaps.int.prop,mel.w.overlaps.int.prop, mel.overlaps.int.prop),
                            overlap.prop.pseudo.std_err=c(era.e.overlap.prop.pseudo.std_err, era.w.overlap.prop.pseudo.std_err, era.overlap.prop.pseudo.std_err, 
                                               mel.e.overlap.prop.pseudo.std_err, mel.w.overlap.prop.pseudo.std_err, mel.overlap.prop.pseudo.std_err),
                            species= c(rep("era", 3), rep("mel", 3)),
                            type.comp=c(rep(c("within", "within", "across"),2)),
                            type.side=c(rep(c("e", "w", "across"),2))); sharing.percs.50kb 

write.csv(sharing.percs.50kb, "local/data/local/data/data/sharing/pbs.out/sharing.percs.50kb.csv", row.names = FALSE)

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

# 6 plots separately 
# era.w.within
era.w.within <- ggplot(data = density.10ksim.50kb.era.w.within, aes(x=perc, y=density)) + 
  geom_line() +
  geom_ribbon(aes(ymin=0, ymax=pmax(density,0)), fill="grey", col="grey", alpha=0.5)+
  scale_x_continuous(limits = c(0,75), expand = c(0, 0)) +
  scale_y_continuous( expand = c(0, 0)) +
  theme_classic()+
  geom_vline(aes(xintercept =sharing.percs.50kb.era.w.within$overlap.int.prop))+
  geom_vline(aes(xintercept =sharing.percs.50kb.era.w.within$overlap.int.prop+sharing.percs.50kb.era.w.within$overlap.prop.pseudo.std_err*1.96), lty = 3)+
  geom_vline(aes(xintercept =sharing.percs.50kb.era.w.within$overlap.int.prop-sharing.percs.50kb.era.w.within$overlap.prop.pseudo.std_err*1.96), lty = 3); era.w.within

# era.w.within
era.e.within <- ggplot(data = density.10ksim.50kb.era.e.within, aes(x=perc, y=density)) + 
  geom_line() +
  geom_ribbon(aes(ymin=0, ymax=pmax(density.10ksim.50kb.era.e.within$density,0)), fill="grey", col="grey", alpha=0.5)+
  scale_x_continuous(limits = c(0,75), expand = c(0, 0)) +
  scale_y_continuous( expand = c(0, 0)) +
  theme_classic()+
  geom_vline(aes(xintercept =sharing.percs.50kb.era.e.within$overlap.int.prop))+
  geom_vline(aes(xintercept =sharing.percs.50kb.era.e.within$overlap.int.prop+sharing.percs.50kb.era.e.within$overlap.prop.pseudo.std_err*1.96), lty = 3)+
  geom_vline(aes(xintercept =sharing.percs.50kb.era.e.within$overlap.int.prop-sharing.percs.50kb.era.e.within$overlap.prop.pseudo.std_err*1.96), lty = 3); era.e.within

# mel.w.within
mel.w.within <- ggplot(data = density.10ksim.50kb.mel.w.within, aes(x=perc, y=density)) + 
  geom_line() +
  geom_ribbon(aes(ymin=0, ymax=pmax(density.10ksim.50kb.mel.w.within$density,0)), fill="grey", col="grey", alpha=0.5)+
  scale_x_continuous(limits = c(0,75), expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme_classic()+
  geom_vline(aes(xintercept =sharing.percs.50kb.mel.w.within$overlap.int.prop))+
  geom_vline(aes(xintercept =sharing.percs.50kb.mel.w.within$overlap.int.prop+sharing.percs.50kb.mel.w.within$overlap.prop.pseudo.std_err*1.96), lty = 3)+
  geom_vline(aes(xintercept =sharing.percs.50kb.mel.w.within$overlap.int.prop-sharing.percs.50kb.mel.w.within$overlap.prop.pseudo.std_err*1.96), lty = 3); mel.w.within

# mel.w.within
mel.e.within <- ggplot(data = density.10ksim.50kb.mel.e.within, aes(x=perc, y=density)) + 
  geom_line() +
  geom_ribbon(aes(ymin=0, ymax=pmax(density.10ksim.50kb.mel.e.within$density,0)), fill="grey", col="grey", alpha=0.5)+
  scale_x_continuous(limits = c(0,75), expand = c(0, 0)) +
  scale_y_continuous( expand = c(0, 0)) +
  theme_classic()+
  geom_vline(aes(xintercept =sharing.percs.50kb.mel.e.within$overlap.int.prop))+
  geom_vline(aes(xintercept =sharing.percs.50kb.mel.e.within$overlap.int.prop+sharing.percs.50kb.mel.e.within$overlap.prop.pseudo.std_err*1.96), lty = 3)+
  geom_vline(aes(xintercept =sharing.percs.50kb.mel.e.within$overlap.int.prop-sharing.percs.50kb.mel.e.within$overlap.prop.pseudo.std_err*1.96), lty = 3); mel.e.within

# era.across

era.across<- ggplot(data = density.10ksim.50kb.era.across.of.total, aes(x=perc, y=density)) + 
  geom_line() +
  geom_ribbon(aes(ymin=0, ymax=pmax(density.10ksim.50kb.era.across.of.total$density,0)), fill="grey", col="grey", alpha=0.5)+
  scale_x_continuous(limits = c(0,75), expand = c(0, 0)) +
  scale_y_continuous( expand = c(0, 0)) +
  theme_classic()+
  geom_vline(aes(xintercept =sharing.percs.50kb.era.across.of.total$overlap.int.prop), color="red", lwd=2)+
  geom_vline(aes(xintercept =sharing.percs.50kb.era.across.of.total$overlap.int.prop + sharing.percs.50kb.era.across.of.total$overlap.prop.pseudo.std_err*1.96), lty = 3, color="red") +
  geom_vline(aes(xintercept =sharing.percs.50kb.era.across.of.total$overlap.int.prop - sharing.percs.50kb.era.across.of.total$overlap.prop.pseudo.std_err*1.96), lty = 3, color="red") ; era.across


# mel.across
mel.across<- ggplot(data = density.10ksim.50kb.mel.across.of.total, aes(x=perc, y=density)) + 
  geom_line() +
  geom_ribbon(aes(ymin=0, ymax=pmax(density.10ksim.50kb.mel.across.of.total$density,0)), fill="grey", col="grey", alpha=0.5)+
  scale_x_continuous(limits = c(0,75), expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme_classic()+
  geom_vline(aes(xintercept =sharing.percs.50kb.mel.across.of.total$overlap.int.prop), color="red", lwd=2)+
  geom_vline(aes(xintercept =sharing.percs.50kb.mel.across.of.total$overlap.int.prop + sharing.percs.50kb.mel.across.of.total$overlap.prop.pseudo.std_err*1.96), lty = 3, color="red") +
  geom_vline(aes(xintercept =sharing.percs.50kb.mel.across.of.total$overlap.int.prop - sharing.percs.50kb.mel.across.of.total$overlap.prop.pseudo.std_err*1.96), lty = 3, color="red") ; mel.across


sharing.all.50kb <- plot_grid(era.w.within, era.across, era.e.within, mel.w.within, mel.across, mel.e.within, nrow = 2,
                               labels = c("era.w.within","era.across", "era.e.within", "mel.w.within", "mel.across", "mel.e.within"), label_y = 1)
save_plot("local/plots/sharing.sims/pbs.out/sharing.all.50kb.pdf", sharing.all.50kb, base_width  = 12, base_height = 6)



########### 6. add sharing data to  tables and save ###########
block.summ.era$overlap.int.prop <-NA; block.summ.era$overlap.prop.pseudo.std_err <- NA
block.summ.era$overlap.int.prop[1] <- era.e.overlaps.int.propA; block.summ.era$overlap.int.prop[3] <- era.e.overlaps.int.propB
block.summ.era$overlap.int.prop[2] <- era.w.overlaps.int.propA; block.summ.era$overlap.int.prop[4] <- era.w.overlaps.int.propB
block.summ.era$overlap.prop.pseudo.std_err[1] <- era.e.overlap.prop.pseudo.std_errA; block.summ.era$overlap.prop.pseudo.std_err[3] <- era.e.overlap.prop.pseudo.std_errB
block.summ.era$overlap.prop.pseudo.std_err[2] <- era.w.overlap.prop.pseudo.std_errA; block.summ.era$overlap.prop.pseudo.std_err[4] <- era.w.overlap.prop.pseudo.std_errB
block.summ.era$sim10k.overlap.prop <-  NA
block.summ.era$sim10k.overlap.prop[1] <- mean(sim10k.era.e.overlap.propA); block.summ.era$sim10k.overlap.prop[3] <- mean(sim10k.era.e.overlap.propB)
block.summ.era$sim10k.overlap.prop[2] <- mean(sim10k.era.w.overlap.propA); block.summ.era$sim10k.overlap.prop[4] <- mean(sim10k.era.w.overlap.propB)

# add across sides data 
block.summ.era <-  rbind(block.summ.era, 
                         c("e", (mean(as.matrix(era.e.overlap.union)[,2] - as.matrix(era.e.overlap.union)[,1] ))/1000, 
                         nrow(era.e.overlap.union),
                         (1-(sum(ref.scaff.era[,2]) - sum(as.matrix(era.e.overlap.union)[,2]- as.matrix(era.e.overlap.union)[,1] ))/sum(ref.scaff.era[,2]))*100,
                         era.overlaps.int.propA,
                         era.overlap.prop.pseudo.std_errA,
                         mean(sim10k.era.overlap.across.propA)))

block.summ.era <-  rbind(block.summ.era, 
                         c("w", (mean(as.matrix(era.w.overlap.union)[,2] - as.matrix(era.w.overlap.union)[,1] ))/1000, 
                           nrow(era.w.overlap.union),
                           (1-(sum(ref.scaff.era[,2]) - sum(as.matrix(era.w.overlap.union)[,2]- as.matrix(era.w.overlap.union)[,1] ))/sum(ref.scaff.era[,2]))*100,
                           era.overlaps.int.propA,
                           era.overlap.prop.pseudo.std_errB,
                           mean(sim10k.era.overlap.across.propB)))

block.summ.era$sim10k.overlap.prop.first.quart <- c(quantile(sim10k.era.e.overlap.propA)[2],quantile(sim10k.era.w.overlap.propA)[2],
                                                    quantile(sim10k.era.e.overlap.propB)[2],quantile(sim10k.era.w.overlap.propB)[2],
                                                    quantile(sim10k.era.overlap.across.propA)[2],quantile(sim10k.era.overlap.across.propB)[2]); block.summ.era
block.summ.era$sim10k.overlap.prop.third.quart <- c(quantile(sim10k.era.e.overlap.propA)[4],quantile(sim10k.era.w.overlap.propA)[4],
                                                    quantile(sim10k.era.e.overlap.propB)[4],quantile(sim10k.era.w.overlap.propB)[4],
                                                    quantile(sim10k.era.overlap.across.propA)[4],quantile(sim10k.era.overlap.across.propB)[4]); block.summ.era
block.summ.era$type.comp <- c(rep("within", 4), rep("across", 2))
block.summ.era$country.comp <- c(rep("co", 2), rep("ec", 2), rep("across", 2)); block.summ.era
write.table(block.summ.era, "local/plots/sharing.sims/pbs.out/block.summ.era.50kb.txt", row.names = FALSE)


block.summ.mel$overlap.int.prop <-NA; block.summ.mel$overlap.prop.pseudo.std_err <- NA
block.summ.mel$overlap.int.prop[1] <- mel.e.overlaps.int.propA; block.summ.mel$overlap.int.prop[3] <- mel.e.overlaps.int.propB
block.summ.mel$overlap.int.prop[2] <- mel.w.overlaps.int.propA; block.summ.mel$overlap.int.prop[4] <- mel.w.overlaps.int.propB
block.summ.mel$overlap.prop.pseudo.std_err[1] <- mel.e.overlap.prop.pseudo.std_errA; block.summ.mel$overlap.prop.pseudo.std_err[3] <- mel.e.overlap.prop.pseudo.std_errB
block.summ.mel$overlap.prop.pseudo.std_err[2] <- mel.w.overlap.prop.pseudo.std_errA; block.summ.mel$overlap.prop.pseudo.std_err[4] <- mel.w.overlap.prop.pseudo.std_errB
block.summ.mel$sim10k.overlap.prop <-  NA
block.summ.mel$sim10k.overlap.prop[1] <- mean(sim10k.mel.e.overlap.propA); block.summ.mel$sim10k.overlap.prop[3] <- mean(sim10k.mel.e.overlap.propB)
block.summ.mel$sim10k.overlap.prop[2] <- mean(sim10k.mel.w.overlap.propA); block.summ.mel$sim10k.overlap.prop[4] <- mean(sim10k.mel.w.overlap.propB)
# add across sides data 
block.summ.mel <-  rbind(block.summ.mel, 
                         c("e", (mean(as.matrix(mel.e.overlap.union)[,2] - as.matrix(mel.e.overlap.union)[,1] ))/1000, 
                           nrow(mel.e.overlap.union),
                           (1-(sum(ref.scaff.mel[,2]) - sum(as.matrix(mel.e.overlap.union)[,2]- as.matrix(mel.e.overlap.union)[,1] ))/sum(ref.scaff.mel[,2]))*100,
                           mel.overlaps.int.propA,
                           mel.overlap.prop.pseudo.std_errA,
                           mean(sim10k.mel.overlap.across.propA)))

block.summ.mel <-  rbind(block.summ.mel, 
                         c("w", (mean(as.matrix(mel.w.overlap.union)[,2] - as.matrix(mel.w.overlap.union)[,1] ))/1000, 
                           nrow(mel.w.overlap.union),
                           (1-(sum(ref.scaff.mel[,2]) - sum(as.matrix(mel.w.overlap.union)[,2]- as.matrix(mel.w.overlap.union)[,1] ))/sum(ref.scaff.mel[,2]))*100,
                           mel.overlaps.int.propA,
                           mel.overlap.prop.pseudo.std_errB,
                           mean(sim10k.mel.overlap.across.propB)))

block.summ.mel$sim10k.overlap.prop.first.quart <- c(quantile(sim10k.mel.e.overlap.propA)[2],quantile(sim10k.mel.w.overlap.propA)[2],
                                                    quantile(sim10k.mel.e.overlap.propB)[2],quantile(sim10k.mel.w.overlap.propB)[2],
                                                    quantile(sim10k.mel.overlap.across.propA)[2],quantile(sim10k.mel.overlap.across.propB)[2]); block.summ.mel
block.summ.mel$sim10k.overlap.prop.third.quart <- c(quantile(sim10k.mel.e.overlap.propA)[4],quantile(sim10k.mel.w.overlap.propA)[4],
                                                    quantile(sim10k.mel.e.overlap.propB)[4],quantile(sim10k.mel.w.overlap.propB)[4],
                                                    quantile(sim10k.mel.overlap.across.propA)[4],quantile(sim10k.mel.overlap.across.propB)[4]); block.summ.mel

write.table(block.summ.mel, "local/plots/sharing.sims/pbs.out/block.summ.mel.50kb.txt", row.names = FALSE)

