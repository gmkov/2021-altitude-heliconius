##### GMK Jan 2021 #####

rm(list=ls())
dev.off()
#### packages ####

library(windowscanr)
library(devtools)
library("windowscanr")
library("ggplot2")
library(cowplot)
library(data.table)
library(dplyr)
library(tidyr)
library(zoo)
library(ape)
library(stringr)
#BiocManager::install("ggtree")
library(ggtree)
library(treeio)
library(phytools)
library(ggpubr)
library(ggbeeswarm)
library(patchwork)
library(data.table)
library(intervals)
library(ggrepel)
library(dplyr)
library(ggtree)
setwd("/Users/gabrielamontejokovacevich/Dropbox (Cambridge University)/PhD/7_Assocation_studies/9_ANGSD/12.era.mel.altitude/local/")
library(ggpubr)
library(EnvStats)
library(ape)
library(ape); library(adephylo)


library(colormap)
scales::show_col(colormap(colormap=colormaps$temperature, nshades=20))
colormap(colormap=colormaps$temperature, nshades=20)[14] #closely related, instead of orange #eb8055ff
colormap(colormap=colormaps$temperature, nshades=20)[19] # distantly related, instead of pink  "#efe350ff"


############################ functions ############################
se <- function(x) sd(x)/sqrt(length(x))
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
  prefix="local/data/joana.fd/output/eratoRelatives.withSRA.chr1-21.max0.5N.minDP3.biallSNPs.mac2."
  
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


getFd.mel<-function(combination){
  prefix="local/data/joana.fd/output/melpomeneRelatives.withSRA.chr1-21.max0.5N.minDP3.biallSNPs.mac2."
  
  # Read in the fd.df
  fd.df<-read.csv(paste0(prefix,combination,"_50kb"))
  
  # Correct negative values
  fd.df<-corrNegfd(fd.df)
  
  # reformat to datatable
  fd.df<-data.table(fd.df)
  
  # Add chromosome name
  fd.df$chr<-substring(fd.df$scaffold,first = 6,last=7)
  
  # add window information
  fd.df<-data.table(fd.df)
  fd.df$start.WG <-fd.df$start+ref.scaff.mel[match(fd.df$scaffold,ref.scaff.mel$scaffold),"add"]
  fd.df$end.WG <-fd.df$end+ref.scaff.mel[match(fd.df$scaffold,ref.scaff.mel$scaffold),"add"]
  
  # add if east (only) PBS outlier (para SHDR)
  shdr.east.ranges<-data.table(mel.outliersEast[,c("scaff","start.bp","end.bp")])
  names(shdr.east.ranges)<-c("scaffold","start","end") # to match fd.df names
  setkey(shdr.east.ranges)
  # by.x concatenates scaffold and start/end
  oveast<-foverlaps(x = fd.df,y = shdr.east.ranges, which=T, by.x = names(shdr.east.ranges), type="any",nomatch=0L)
  fd.df$shdr.para.east.id<- NA
  # use keys to select corresponding shdr names
  fd.df$shdr.para.east.id[oveast$xid]<- paste(mel.outliersEast[ oveast$yid,]$shdr.para.east.id)
  
  # add if west (only) PBS outlier (para SHDR)
  shdr.west.ranges<-data.table(mel.outliersWest[,c("scaff","start.bp","end.bp")])
  names(shdr.west.ranges)<-c("scaffold","start","end") # to match fd.df names
  setkey(shdr.west.ranges)
  # by.x concatenates scaffold and start/end
  ovWest<-foverlaps(x = fd.df,y = shdr.west.ranges, which=T, by.x = names(shdr.west.ranges), type="any",nomatch=0L)
  fd.df$shdr.para.west.id<- NA
  # use keys to select corresponding shdr names
  fd.df$shdr.para.west.id[ovWest$xid]<-  paste(mel.outliersWest[ ovWest$yid,]$shdr.para.west.id)
  
  # add if across (only) PBS outlier (para SHDR)
  shdr.across.ranges<-data.table(mel.outliersAcross[,c("scaff","start.bp","end.bp")])
  names(shdr.across.ranges)<-c("scaffold","start","end") # to match fd.df names
  setkey(shdr.across.ranges)
  # by.x concatenates scaffold and start/end
  ovacross<-foverlaps(x = fd.df,y = shdr.across.ranges, which=T, by.x = names(shdr.across.ranges), type="any",nomatch=0L)
  fd.df$shdr.allo.id<- NA
  # use keys to select corresponding shdr names
  fd.df$shdr.allo.id[ovacross$xid]<-  paste(mel.outliersAcross[ ovacross$yid,]$shdr.allo.id)
  
  
  return(fd.df)
}

# Find fd windows that overlap with East outliers
eastOverlap<-function(dataset){
  # get outlier columns:
  rangesB<-data.table(outliersEast[,c("scaff","start","end")])
  names(rangesB)[1]<-"scaffold"
  
  # Sorts the data.table by the columns provided (here by all columns)
  setkey(rangesB)
  
  dataOverlap<-foverlaps(x = dataset,y = rangesB,  by.x = names(rangesB), type="any",nomatch=0L)
  return(dataOverlap)
}

# Find fd windows that overlap with West outliers
westOverlap<-function(dataset){
  # get outlier columns:
  rangesB<-data.table(outliersWest[,c("scaff","start","end")])
  names(rangesB)[1]<-"scaffold"
  
  # Sorts the data.table by the columns provided (here by all columns)
  setkey(rangesB)
  
  dataOverlap<-foverlaps(x = dataset,y = rangesB,  by.x = names(rangesB), type="any",nomatch=0L)
  return(dataOverlap)
}

# Find fd windows that overlap with outliers shared across sides of the Andes
acrossOverlap<-function(dataset){
  # get outlier columns:
  rangesB<-data.table(outliersAcross[,c("scaff","start","end")])
  names(rangesB)[1]<-"scaffold"
  
  # Sorts the data.table by the columns provided (here by all columns)
  setkey(rangesB)
  
  dataOverlap<-foverlaps(x = dataset,y = rangesB,  by.x = names(rangesB), type="any",nomatch=0L)
  return(dataOverlap)
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

# era twisst has 20k windows (smaller than SHDR ranges)
# import twisst and whether overlaps with SHDRs
getTWISST.era<-function(combination){

# Read in weights twisstset
twisst<-read.table(paste0(prefix,combination,".weights.csv"),header=T)

# add window information
twisst<-cbind(info,twisst)
twisst<-data.table(twisst)
twisst$start.WG <-twisst$start+ref.scaff.era[match(twisst$scaffold,ref.scaff.era$scaffold),"add"]
twisst$end.WG <-twisst$end+ref.scaff.era[match(twisst$scaffold,ref.scaff.era$scaffold),"add"]

# topo weightings (0-1) per window
twisst$topo1.freq <- round(twisst$topo1/rowSums(twisst[,c("topo1","topo2","topo3")]),digits=4)
twisst$topo2.freq <- round(twisst$topo2/rowSums(twisst[,c("topo1","topo2","topo3")]),digits=4)
twisst$topo3.freq <- round(twisst$topo3/rowSums(twisst[,c("topo1","topo2","topo3")]),digits=4)
# freq topo2 -topo3
twisst$topo2.minus.topo3.freq <- twisst$topo2.freq - twisst$topo3.freq


# add if east (only) PBS outlier (para SHDR)
shdr.east.ranges<-data.table(era.outliersEast[,c("scaff","start.bp","end.bp")])
names(shdr.east.ranges)<-c("scaffold","start","end") # to match twisst names
setkey(shdr.east.ranges)
# by.x concatenates scaffold and start/end
oveast<-foverlaps(x = twisst,y = shdr.east.ranges, which=T, by.x = names(shdr.east.ranges), type="any",nomatch=0L)
twisst$shdr.para.east.id<- NA
# use keys to select corresponding shdr names
twisst$shdr.para.east.id[oveast$xid]<- paste(era.outliersEast[ oveast$yid,]$shdr.para.east.id)

# add if west (only) PBS outlier (para SHDR)
shdr.west.ranges<-data.table(era.outliersWest[,c("scaff","start.bp","end.bp")])
names(shdr.west.ranges)<-c("scaffold","start","end") # to match twisst names
setkey(shdr.west.ranges)
# by.x concatenates scaffold and start/end
ovWest<-foverlaps(x = twisst,y = shdr.west.ranges, which=T, by.x = names(shdr.west.ranges), type="any",nomatch=0L)
twisst$shdr.para.west.id<- NA
# use keys to select corresponding shdr names
twisst$shdr.para.west.id[ovWest$xid]<-  paste(era.outliersWest[ ovWest$yid,]$shdr.para.west.id)

# add if across (only) PBS outlier (para SHDR)
shdr.across.ranges<-data.table(era.outliersAcross[,c("scaff","start.bp","end.bp")])
names(shdr.across.ranges)<-c("scaffold","start","end") # to match twisst names
setkey(shdr.across.ranges)
# by.x concatenates scaffold and start/end
ovacross<-foverlaps(x = twisst,y = shdr.across.ranges, which=T, by.x = names(shdr.across.ranges), type="any",nomatch=0L)
twisst$shdr.allo.id<- NA
# use keys to select corresponding shdr names
twisst$shdr.allo.id[ovacross$xid]<-  paste(era.outliersAcross[ ovacross$yid,]$shdr.allo.id)

return(twisst)
}

getTWISST.era<-function(combination){
  
  # Read in weights twisstset
  twisst<-read.table(paste0(prefix,combination,".weights.csv"),header=T)
  
  # add window information
  twisst<-cbind(info,twisst)
  twisst<-data.table(twisst)
  twisst$start.WG <-twisst$start+ref.scaff.era[match(twisst$scaffold,ref.scaff.era$scaffold),"add"]
  twisst$end.WG <-twisst$end+ref.scaff.era[match(twisst$scaffold,ref.scaff.era$scaffold),"add"]
  
  # topo weightings (0-1) per window
  twisst$topo1.freq <- round(twisst$topo1/rowSums(twisst[,c("topo1","topo2","topo3")]),digits=4)
  twisst$topo2.freq <- round(twisst$topo2/rowSums(twisst[,c("topo1","topo2","topo3")]),digits=4)
  twisst$topo3.freq <- round(twisst$topo3/rowSums(twisst[,c("topo1","topo2","topo3")]),digits=4)
  # freq topo2 -topo3
  twisst$topo2.minus.topo3.freq <- twisst$topo2.freq - twisst$topo3.freq
  
  
  # add if east (only) PBS outlier (para SHDR)
  shdr.east.ranges<-data.table(era.outliersEast[,c("scaff","start.bp","end.bp")])
  names(shdr.east.ranges)<-c("scaffold","start","end") # to match twisst names
  setkey(shdr.east.ranges)
  # by.x concatenates scaffold and start/end
  oveast<-foverlaps(x = twisst,y = shdr.east.ranges, which=T, by.x = names(shdr.east.ranges), type="any",nomatch=0L)
  twisst$shdr.para.east.id<- NA
  # use keys to select corresponding shdr names
  twisst$shdr.para.east.id[oveast$xid]<- paste(era.outliersEast[ oveast$yid,]$shdr.para.east.id)
  
  # add if west (only) PBS outlier (para SHDR)
  shdr.west.ranges<-data.table(era.outliersWest[,c("scaff","start.bp","end.bp")])
  names(shdr.west.ranges)<-c("scaffold","start","end") # to match twisst names
  setkey(shdr.west.ranges)
  # by.x concatenates scaffold and start/end
  ovWest<-foverlaps(x = twisst,y = shdr.west.ranges, which=T, by.x = names(shdr.west.ranges), type="any",nomatch=0L)
  twisst$shdr.para.west.id<- NA
  # use keys to select corresponding shdr names
  twisst$shdr.para.west.id[ovWest$xid]<-  paste(era.outliersWest[ ovWest$yid,]$shdr.para.west.id)
  
  # add if across (only) PBS outlier (para SHDR)
  shdr.across.ranges<-data.table(era.outliersAcross[,c("scaff","start.bp","end.bp")])
  names(shdr.across.ranges)<-c("scaffold","start","end") # to match twisst names
  setkey(shdr.across.ranges)
  # by.x concatenates scaffold and start/end
  ovacross<-foverlaps(x = twisst,y = shdr.across.ranges, which=T, by.x = names(shdr.across.ranges), type="any",nomatch=0L)
  twisst$shdr.allo.id<- NA
  # use keys to select corresponding shdr names
  twisst$shdr.allo.id[ovacross$xid]<-  paste(era.outliersAcross[ ovacross$yid,]$shdr.allo.id)
  
  return(twisst)
}

## to be able to add a secondary axes to discrete x....
guide_axis_label_trans <- function(label_trans = identity, ...) {
  axis_guide <- guide_axis(...)
  axis_guide$label_trans <- rlang::as_function(label_trans)
  class(axis_guide) <- c("guide_axis_trans", class(axis_guide))
  axis_guide}

guide_train.guide_axis_trans <- function(x, ...) {
  trained <- NextMethod()
  trained$key$.label <- x$label_trans(trained$key$.label)
  trained}

#######################################################################################################   ERATO fdm ################################################################

######################################################## 0. data prep ############################
#### Read in outlier regions and chromosome information ####

# Get erato outliers shared by either or both East and West:
era.shdr.para.east<-read.csv("local/data/sharing/era.shdr.para.east.csv",header=T); head(era.shdr.para.east)
era.shdr.para.west<-read.csv("local/data/sharing/era.shdr.para.west.csv",header=T)
era.shdr.allo.all <-read.csv("local/data/sharing/era.shdr.allo.all.csv",header=T)

# remove 7 weird shdr erato
era.shdr.para.east <- subset(era.shdr.para.east, is.max.pbs.hig.above4=="yes")
era.shdr.para.west <- subset(era.shdr.para.west, is.max.pbs.hig.above4=="yes")
era.shdr.allo.all <- subset(era.shdr.allo.all, is.max.pbs.hig.above4=="yes")

# add 
era.shdr.allo.all.summ <-read.csv("local/data/sharing/era.shdr.allo.summ.csv",header=T)


# unlike when working with Tajima'sD etc, where we were comparing selection statistics obtained from pops in different sides of the Andes
# twisst includes all populations (both sides), so we dont want to divide the analyses between SHDR west / east
# so best we use SHDR ranges of para.east/west (only, no allo) and allopatric ranges
era.outliersEast<-era.shdr.para.east[era.shdr.para.east$is.allopatric=="no",]
# use era.shdr.allo.all because it'll have max start/ends (all overlap between east/west)
era.outliersAcross<-era.shdr.allo.all
era.outliersWest<-era.shdr.para.west[era.shdr.para.west$is.allopatric=="no",]


# Read in chromosome information
ref.scaff.era <- read.table("local/data/ref/Heliconius_erato_demophoon_v1_-_scaffolds.fa.fai", row.names = NULL)
ref.scaff.era<-ref.scaff.era[,1:2]
names(ref.scaff.era)<-c("scaffold","length")
ref.scaff.era$add<-c(0,cumsum(ref.scaff.era$length)[-length(ref.scaff.era$length)])
ref.scaff.era$chr<-substring(ref.scaff.era$scaffold,first = 7,last=8);head(ref.scaff.era)

scafEnds.era <- cumsum(ref.scaff.era[,2])
era.genome_size <- tail(scafEnds.era,1)


################################################################################################################ fig6  VERSION 2 ########################################################
tree.kozak <- read.tree(file = "local/data/joana.tree/Kozak2015tree.tip.sp.nwk")
tree.kozak$tip.label
tips.to.keep <- c("Heliconius_erato", "Heliconius_himera", "Heliconius_clysonymus", "Heliconius_telesiphe",
                  "Heliconius_melpomene", "Heliconius_timareta", "Heliconius_cydno")
tips.to.drop <- tree.kozak$tip.label[!(tree.kozak$tip.label %in% tips.to.keep)]


tree.kozak.sub <- drop.tip(tree.kozak, tips.to.drop, trim.internal = TRUE, subtree = FALSE, root.edge = 0, rooted = is.rooted(tree.kozak), collapse.singles = TRUE, interactive = FALSE); plot(tree.kozak.sub)

tree.kozak.sub$tip.label <- c("H. cydno W"  ,  "H. timareta E", "H. melpomene W", "H. melpomene FG", "H. melpomene E", "H. himera E" ,"H. erato E" , "H. erato FG",  "H. erato W"  ,  "H. clysonymus W/E","H. telesiphe E"); tree.kozak.sub$tip.label
tree.kozak.sub <- drop.tip(tree.kozak.sub, c("H. melpomene FG", "H. erato FG"), trim.internal = TRUE, subtree = FALSE,
                           root.edge = 0, rooted = is.rooted(tree.kozak.sub), collapse.singles = TRUE,
                           interactive = FALSE); plot(tree.kozak.sub)

tree.kozak.sub$tip.label <- c("cydno W"  ,  "timareta E", "melpomene W", "melpomene E", "himera E" ,"erato E" , "erato W"  ,  "clysonymus W/E","telesiphe E"); tree.kozak.sub$tip.label

samples.tree.info <- data_frame(tip.labels=  c("cydno W"  ,  "timareta E", "melpomene W", "melpomene E", "himera E" ,"erato E" , "erato W"  ,  "clysonymus W/E","telesiphe E"),
                                sample.type=c("close", "close","mel.era", "mel.era", "close", "mel.era", "mel.era", "distant",   "distant"   )); samples.tree.info


tree.kozak.sub$tip.label <- c("cydW"  ,  "timE", "melW", "melE", "himE" ,"eraE" , "eraW"  ,  "clysW/E","telE"); tree.kozak.sub$tip.label

samples.tree.info <- data_frame(tip.labels=  c("cydW"  ,  "timE", "melW", "melE", "himE" ,"eraE" , "eraW"  ,  "clysW/E","telE"),
                                sample.type=c("close", "close","mel.era", "mel.era", "close", "mel.era", "mel.era", "distant",   "distant"   )); samples.tree.info


ggtree(tree.kozak.sub, ladderize =F , aes(color=sample.type), size=2) %<+% samples.tree.info  + 
  #geom_tiplab(aes(label=paste0('italic(', genus, ')~bolditalic(', species, ')~', geo)), parse=T)
  geom_tiplab( align = T, hjust = -0.05, size=7, color="black")+ xlim(0,26)+
  #geom_tippoint(aes(fill=sample.type),  size=5,shape = 24 )+
  scale_color_manual(values = c("#eb8055ff", "#efe350ff",'blue'))+
  theme(legend.position = "none")

ggsave("plots/joana.tree/era.mel.kozak.main.tree.png", width = 5.5, height = 4, bg="white")

## split cly e/w, flip some
str(tree.kozak.sub)
tree.kozak.sub$edge.length
node <- which(tree.kozak.sub$tip.label=="clysW/E")
plot(tree.kozak.sub) ; tree <- bind.tip(tree.kozak.sub, tip.label="clyE",  where=node, position=0.3); plot(tree)
tree$tip.label[8] <- c("clyW")
which(tree.kozak.sub$tip.label=="clysW/E")
flip(tree, 2,4) 

samples.tree.info1 <- data_frame(tip.labels=  c("cydW"  ,  "timE", "melW", "melE", "himE" ,"eraE" , "eraW"  ,  "clyW", "clyE","telE"),
                                sample.type=c("close", "close","mel.era", "mel.era", "close", "mel.era", "mel.era", "distant",   "distant", "distant"   )); samples.tree.info1



p <-ggtree(tree, ladderize =F , aes(color=sample.type), size=2) %<+% samples.tree.info1   + 
  #geom_tiplab(aes(label=paste0('italic(', genus, ')~bolditalic(', species, ')~', geo)), parse=T)
  #geom_text(aes(label=node), size=6)+
  geom_tiplab( align = T, hjust = -0.05, size=7, color="black")+ xlim(0,14)+
  #geom_tippoint(aes(fill=sample.type),  size=5,shape = 24 )+
  scale_color_manual(values = c("#eb8055ff", "#efe350ff",'blue'))+
  theme(legend.position = "none")+
  theme(panel.border = element_rect(colour = "NA", fill=NA, size=.8),legend.position="none", axis.ticks.length=unit(0.15, "cm"),
        axis.title = element_text(size=16), panel.spacing = unit(1, "lines"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),  strip.text = element_blank(), strip.background = element_blank(),
        panel.background = element_rect(fill = "transparent", color = NA), plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
        legend.background = element_rect(fill = "transparent", color = NA), legend.box.background = element_rect(fill = "transparent", color = NA), legend.key = element_rect(fill = "transparent", color = NA),
        strip.background.y = element_blank(), rect = element_rect(fill = "transparent")) ; plot(p); 

era.mel.tree <- flip(p, 18,16) %>% flip( 4,3); era.mel.tree




##### read in shdr/ fdm outlier data #####
shdr.era.east.outlier.df <- read.csv("local/data/shdr.summ/shdr.era.east.outlier.df.csv")
shdr.era.west.outlier.df <- read.csv("local/data/shdr.summ/shdr.era.west.outlier.df.csv")

# first get rid off shitty hdrs
names(shdr.era.east.outlier.df)
shdr.era.east.outlier.df <- subset(shdr.era.east.outlier.df, is.max.pbs.hig.above4=="yes")
shdr.era.west.outlier.df <- subset(shdr.era.west.outlier.df, is.max.pbs.hig.above4=="yes")

##### read in simulated shdr fdm thresholds #####
shdr.sims.fdm.era <- read.csv("local/data/joana.fd/era.within.and.other.species.1ksims.fdm.threhsolds.means.csv"); head(shdr.sims.fdm.era)


######################################################## 2. era fd/fdm from other species (parapatric,  n=8 combos (2 tele, 2 him, 4 clys)) ############################


############################ 2.4 era plot fdm across pops and species ############################
#### change in % fdm outliers with genetic distance ####
era.fd.list.files
names(shdr.era.east.outlier.df); shdr.era.east.outlier.df$is.allo.hdr

plot(tree.era)
tree.era$tip.label

era.shdr.SAV.other.sp.summ.percs <- data_frame(donor.p3=c(rep("era.west.high", 2), rep("era.east.high", 2), rep("him.east.high", 2),  rep("tel.east.high", 2),  rep("cly.east.high", 2),  rep("cly.west.high", 2)),
                                         recipient.p2=c(rep("era.west.high", 2), rep("era.east.high", 2), rep("him.east.high", 2),  rep("tel.east.high", 2),  rep("cly.east.high", 2),  rep("cly.west.high", 2)),
                                         country.p2= c(rep(c("col", "ecu" ), 6)),
                                         side.p2= c(rep("P2 = era. high East", 2), rep("P2 = era. high West", 2), rep("P2 = era. high East", 6), rep("P2 = era. high West", 2)),
                                         pop.p3= c(rep("1.era", 4), rep("2.him", 2), rep("3.tel", 2), rep("4.cly", 4)),
                                         percentage.shdr.with.excess.allele.sharing= c(
                                           nrow(subset(shdr.era.east.outlier.df, fdm.max.more.90th.abs.fdm.min.per.SHDR_fd_eraCvlowE_eraChighE_eraChighW_eleHW=="yes" ))/ nrow(subset(shdr.era.east.outlier.df)) *100,
                                           nrow(subset(shdr.era.east.outlier.df, fdm.max.more.90th.abs.fdm.min.per.SHDR_fd_eraEvlowE_eraEhighE_eraEhighW_eleHW=="yes"))/ nrow(subset(shdr.era.east.outlier.df)) *100,
                                           nrow(subset(shdr.era.west.outlier.df, fdm.max.more.90th.abs.fdm.min.per.SHDR_fd_eraCvlowW_eraChighW_eraChighE_eleHW=="yes"))/ nrow(subset(shdr.era.west.outlier.df)) *100,
                                           nrow(subset(shdr.era.west.outlier.df, fdm.max.more.90th.abs.fdm.min.per.SHDR_fd_eraEvlowW_eraEhighW_eraEhighE_eleHW=="yes"))/ nrow(subset(shdr.era.west.outlier.df)) *100,
                                           nrow(subset(shdr.era.east.outlier.df, fdm.max.more.90th.abs.fdm.min.per.SHDR_fd_eraCvlowE_eraChighE_himHE_eleHW=="yes"))/ nrow(subset(shdr.era.east.outlier.df)) *100,
                                           nrow(subset(shdr.era.east.outlier.df, fdm.max.more.90th.abs.fdm.min.per.SHDR_fd_eraEvlowE_eraEhighE_himHE_eleHW=="yes"))/ nrow(subset(shdr.era.east.outlier.df)) *100,
                                           nrow(subset(shdr.era.east.outlier.df, fdm.max.more.90th.abs.fdm.min.per.SHDR_fd_eraCvlowE_eraChighE_telhighE_eleHW=="yes"))/ nrow(subset(shdr.era.east.outlier.df)) *100,
                                           nrow(subset(shdr.era.east.outlier.df, fdm.max.more.90th.abs.fdm.min.per.SHDR_fd_eraEvlowE_eraEhighE_telhighE_eleHW=="yes"))/ nrow(subset(shdr.era.east.outlier.df)) *100,
                                           nrow(subset(shdr.era.east.outlier.df, fdm.max.more.90th.abs.fdm.min.per.SHDR_fd_eraCvlowE_eraChighE_clyEhighE_eleHW=="yes"))/ nrow(subset(shdr.era.east.outlier.df)) *100,
                                           nrow(subset(shdr.era.east.outlier.df, fdm.max.more.90th.abs.fdm.min.per.SHDR_fd_eraEvlowE_eraEhighE_clyEhighE_eleHW=="yes"))/ nrow(subset(shdr.era.east.outlier.df)) *100,
                                           nrow(subset(shdr.era.west.outlier.df, fdm.max.more.90th.abs.fdm.min.per.SHDR_fd_eraCvlowW_eraChighW_clyEhighW_eleHW=="yes"))/ nrow(subset(shdr.era.west.outlier.df)) *100,
                                           nrow(subset(shdr.era.west.outlier.df, fdm.max.more.90th.abs.fdm.min.per.SHDR_fd_eraEvlowW_eraEhighW_clyEhighW_eleHW=="yes"))/ nrow(subset(shdr.era.west.outlier.df)) *100),
                                         # dist.p2.p3= c(
                                         #   dists.era[, "CAM040106_eratoChighE"][c("15N120_eratoChighW")],
                                         #   dists.era[, "CAM011441_eratoEhighE"][c("CAM040545_eratoEhighW")],
                                         #   dists.era[, "15N120_eratoChighW"][c("CAM040106_eratoChighE")],
                                         #   dists.era[, "CAM040545_eratoEhighW"][c("CAM011441_eratoEhighE")],
                                         #   dists.era[, "CAM040106_eratoChighE"][c("SRR4032038_himeraHE")],
                                         #   dists.era[, "CAM011441_eratoEhighE"][c("SRR4032038_himeraHE")],
                                         #   dists.era[, "CAM040106_eratoChighE"][c("CAM040276_telesiphEhighE")],
                                         #   dists.era[, "CAM011441_eratoEhighE"][c("CAM040276_telesiphEhighE")],
                                         #   dists.era[, "CAM040106_eratoChighE"][c("CAM040271_clysonymusHE")],
                                         #   dists.era[, "CAM011441_eratoEhighE"][c("CAM040271_clysonymusHE")],
                                         #   dists.era[, "15N120_eratoChighW"][c("CAM002846_clysonymusHW")],
                                         #   dists.era[, "CAM040545_eratoEhighW"][c("CAM002846_clysonymusHW")]
                                         # ),
                                         donor.p3.country= paste(donor.p3, country.p2, sep="_"),
                                         side.country.p2= paste(side.p2, country.p2, sep="_"),
                                         type.p3=c(rep("within.species", 4), rep("closely.related", 2), rep("distantly.related",6)),
                                         type.p3.side.country.p2=paste(type.p3, side.country.p2, sep="_" ),
                                         sp.p3=substr(donor.p3, 0,3),
                                         sp.side.P3 = paste(sp.p3,  toupper(substr(donor.p3, 5,5)), sep = "")
); era.shdr.SAV.other.sp.summ.percs 

era.shdr.SAV.other.sp.summ.percs $type.p3 <- factor(era.shdr.SAV.other.sp.summ.percs $type.p3, levels=c("within.species", "closely.related" ,  "distantly.related"))
era.shdr.SAV.other.sp.summ.percs $side.p2 <- factor(era.shdr.SAV.other.sp.summ.percs $side.p2, levels=c("P2 = era. high West", "P2 = era. high East"))

#### plot ####
ggplot(data=subset(era.shdr.SAV.other.sp.summ.percs  ), 
       aes(x=type.p3 , y=percentage.shdr.with.excess.allele.sharing,  group=side.country.p2, fill=type.p3,  shape=country.p2))+
  #geom_line(inherit.aes = F, color= "black", data=subset(era.shdr.SAV.other.sp.summ, donor.p3!="era.west.high" & donor.p3!="era.east.high" ), aes(x=type.p3 , y=percentage.shdr.with.excess.allele.sharing, group=side.country.p2, linetype=country.p2), size=1)+
  geom_text_repel(aes(label=sp.side.P3), fill="transparent", size=6, segment.size  = 0.8, segment.color = "grey50",nudge_x       = c(-0.4, 0.5))+
  geom_point( size=5, stroke=1, width = 0.1)  +
  scale_shape_manual(values=c(24, 25))+
  scale_fill_manual(values=c("blue", "#eb8055ff","#efe350ff"  ))+
  #scale_color_manual(values=c("black",  "black"))+
  #geom_text(aes(label=donor.p3),hjust=-0.2, vjust=0) +
  #geom_point(inherit.aes = F, data=subset(era.shdr.SAV.other.sp.summ, (donor.p3=="era.west.high" | donor.p3=="era.east.high")  ), aes(x=type.p3 , y=percentage.shdr.with.excess.allele.sharing,  group=side.country.p2),
  #           shape=24, fill="blue", size=5, stroke=1)+
  #geom_text(inherit.aes = F, data=subset(era.shdr.SAV.other.sp.summ, (donor.p3=="era.west.high" | donor.p3=="era.east.high" )&country.p2=="col" ), aes(x=dist.p2.p3, y=percentage.shdr.with.excess.allele.sharing,  group=side.country.p2,  label=donor.p3),hjust=0, vjust=0) +
  #scale_linetype_manual(values=c("solid", "11"))+
  facet_wrap(~side.p2, ncol = 1)+
  xlab("Phylogenetic distance: P2 - P3")+ ylab("% SHDR with evidence for excess allele sharing")+
  scale_y_continuous(expand = c(0.1,0.1), breaks=c(0,20,40))+
  scale_x_discrete( labels=c("Within species\nbut opposite sides of Andes","close","distant"),position = "top")+
  #scale_x_continuous(expand = c(0.002,0.002))+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=.8),legend.position="none", axis.ticks.length=unit(0.15, "cm"),
        plot.margin = unit(c(1,1,1,1), "lines"),axis.text.x = element_blank(), axis.text.y = element_text(size=14), 
        axis.title.x =  element_blank(),axis.title.y = element_text(size=16), panel.spacing = unit(1, "lines"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),  strip.text = element_blank(), strip.background = element_blank(),
        panel.background = element_rect(fill = "transparent", color = NA), plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
        legend.background = element_rect(fill = "transparent", color = NA), legend.box.background = element_rect(fill = "transparent", color = NA), legend.key = element_rect(fill = "transparent", color = NA),
        strip.background.y = element_blank(), rect = element_rect(fill = "transparent"))


ggsave("plots/joana.fdm/fdm.perc.shdrs.outliers.era.other.sp.new.png", height = 7, width = 4.5, bg="transparent")



#### ks tests #####
# null hypothesis that the true distribution function of x is  not greater than the distribution function of y
# alternative hypothesis is that dist of x is greater than dist of y
ks.test(shdr.era.west.outlier.df$fdm.max_fd_eraCvlowW_eraChighW_eraChighE_eleHW, abs(shdr.era.west.outlier.df$fdm.min_fd_eraCvlowW_eraChighW_eraChighE_eleHW), alternative = "less" ) # D^+ = 0.27586, p-value = 0.01211
ks.test(shdr.era.west.outlier.df$fdm.max_fd_eraEvlowW_eraEhighW_eraEhighE_eleHW, abs(shdr.era.west.outlier.df$fdm.min_fd_eraEvlowW_eraEhighW_eraEhighE_eleHW), alternative = "less" ) # D^+ = 0.086207, p-value = 0.6498

ks.test(shdr.era.east.outlier.df$fdm.max_fd_eraCvlowE_eraChighE_eraChighW_eleHW, abs(shdr.era.east.outlier.df$fdm.min_fd_eraCvlowE_eraChighE_eraChighW_eleHW), alternative = "less" )
ks.test(shdr.era.east.outlier.df$fdm.max_fd_eraEvlowE_eraEhighE_eraEhighW_eleHW, abs(shdr.era.east.outlier.df$fdm.min_fd_eraEvlowE_eraEhighE_eraEhighW_eleHW), alternative = "less" )

# across sp
ks.test(shdr.era.east.outlier.df$fdm.max_fd_eraCvlowE_eraChighE_himHE_eleHW, abs(shdr.era.east.outlier.df$fdm.min_fd_eraCvlowE_eraChighE_himHE_eleHW), alternative = "less" ) 
ks.test(shdr.era.east.outlier.df$fdm.max_fd_eraEvlowE_eraEhighE_himHE_eleHW, abs(shdr.era.east.outlier.df$fdm.min_fd_eraEvlowE_eraEhighE_himHE_eleHW), alternative = "less" )

ks.test(shdr.era.east.outlier.df$fdm.max_fd_eraCvlowE_eraChighE_telhighE_eleHW, abs(shdr.era.east.outlier.df$fdm.min_fd_eraCvlowE_eraChighE_telhighE_eleHW), alternative = "less" ) 
ks.test(shdr.era.east.outlier.df$fdm.max_fd_eraEvlowE_eraEhighE_telhighE_eleHW, abs(shdr.era.east.outlier.df$fdm.min_fd_eraEvlowE_eraEhighE_telhighE_eleHW), alternative = "less" )

ks.test(shdr.era.east.outlier.df$fdm.max_fd_eraCvlowE_eraChighE_clyEhighE_eleHW, abs(shdr.era.east.outlier.df$fdm.min_fd_eraCvlowE_eraChighE_clyEhighE_eleHW), alternative = "less" ) 
ks.test(shdr.era.east.outlier.df$fdm.max_fd_eraEvlowE_eraEhighE_clyEhighE_eleHW, abs(shdr.era.east.outlier.df$fdm.min_fd_eraEvlowE_eraEhighE_clyEhighE_eleHW), alternative = "less" )

ks.test(shdr.era.west.outlier.df$fdm.max_fd_eraCvlowW_eraChighW_clyEhighW_eleHW, abs(shdr.era.west.outlier.df$fdm.min_fd_eraCvlowW_eraChighW_clyEhighW_eleHW), alternative = "less" ) # D^+ = 0.27586, p-value = 0.01211
ks.test(shdr.era.west.outlier.df$fdm.max_fd_eraEvlowW_eraEhighW_clyEhighW_eleHW, abs(shdr.era.west.outlier.df$fdm.min_fd_eraEvlowW_eraEhighW_clyEhighW_eleHW), alternative = "less" ) # D^+ = 0.086207, p-value = 0.6498

############################ 2.4 mel plot fdm across pops and species ############################
##### read in shdr/ fdm outlier data #####
shdr.mel.east.outlier.df <- read.csv("local/data/shdr.summ/shdr.mel.east.outlier.df.csv")
shdr.mel.west.outlier.df <- read.csv("local/data/shdr.summ/shdr.mel.west.outlier.df.csv")

# first get rid off shitty hdrs
names(shdr.mel.east.outlier.df)
shdr.mel.east.outlier.df <- subset(shdr.mel.east.outlier.df, is.max.pbs.hig.above4=="yes")
shdr.mel.west.outlier.df <- subset(shdr.mel.west.outlier.df, is.max.pbs.hig.above4=="yes")

##### read in simulated shdr fdm thresholds #####
shdr.sims.fdm.mel <- read.csv("local/data/joana.fd/mel.within.and.other.species.1ksims.fdm.threhsolds.means.csv"); head(shdr.sims.fdm.mel)

#### change in % fdm outliers with genetic distance ####
mel.fd.list.files
names(shdr.mel.east.outlier.df); shdr.mel.east.outlier.df$is.allo.hdr
shdr.mel.west.outlier.df$fdm.max.more.90th.abs.fdm.min.per.SHDR_fd_melCvlowW_melChighW_clyEhighW_eleHW
tree.mel$tip.label

mel.shdr.SAV.other.sp.summ.percs <- data_frame(donor.p3=c( rep("mel.west.high", 2) ,rep("mel.east.high", 2), rep("tim.east.high", 2),  rep("cyd.west.high", 2)),
                                         country.p2= c(rep(c("col", "ecu" ), 4)),
                                         side.p2= c(rep("P2 = mel. high East", 2), rep("P2 = mel. high West", 2), rep("P2 = mel. high East", 2), rep("P2 = mel. high West", 2)),
                                         dist= c(rep("2.mel", 4), rep("1.tim.E", 2), rep("3.cyd.W", 2)),
                                         percentage.shdr.with.excess.allele.sharing= c(
                                           nrow(subset(shdr.mel.east.outlier.df, fdm.max.more.90th.abs.fdm.min.per.SHDR_fd_melCvlowE_melChighE_melChighW_hecale=="yes" ))/ nrow(subset(shdr.mel.east.outlier.df)) *100,
                                           nrow(subset(shdr.mel.east.outlier.df, fdm.max.more.90th.abs.fdm.min.per.SHDR_fd_melEvlowE_melEhighE_melEhighW_hecale=="yes"))/ nrow(subset(shdr.mel.east.outlier.df)) *100,
                                           nrow(subset(shdr.mel.west.outlier.df, fdm.max.more.90th.abs.fdm.min.per.SHDR_fd_melCvlowW_melChighW_melChighE_hecale=="yes"))/ nrow(subset(shdr.mel.west.outlier.df)) *100,
                                           nrow(subset(shdr.mel.west.outlier.df, fdm.max.more.90th.abs.fdm.min.per.SHDR_fd_melElowW_melEhighW_melEhighE_hecale=="yes"))/ nrow(subset(shdr.mel.west.outlier.df)) *100,
                                           nrow(subset(shdr.mel.east.outlier.df, fdm.max.more.90th.abs.fdm.min.per.SHDR_fd_melCvlowE_melChighE_timEhighE_hecale=="yes"))/ nrow(subset(shdr.mel.east.outlier.df)) *100,
                                           nrow(subset(shdr.mel.east.outlier.df, fdm.max.more.90th.abs.fdm.min.per.SHDR_fd_melEvlowE_melEhighE_timEhighE_hecale=="yes"))/ nrow(subset(shdr.mel.east.outlier.df)) *100,
                                           nrow(subset(shdr.mel.west.outlier.df, fdm.max.more.90th.abs.fdm.min.per.SHDR_fd_melCvlowW_melChighW_cydEhighW_hecale=="yes"))/ nrow(subset(shdr.mel.west.outlier.df)) *100,
                                           nrow(subset(shdr.mel.west.outlier.df, fdm.max.more.90th.abs.fdm.min.per.SHDR_fd_melElowW_melEhighW_cydEhighW_hecale=="yes"))/ nrow(subset(shdr.mel.west.outlier.df)) *100),
                                         # dist.p2.p3= c(
                                         #   dists.mel[, "CAM040043_melChighE"][c("CS003603_melChighW")],
                                         #   dists.mel[, "CAM041101_melEhighE"][c("14N004_melEhighW")],
                                         #   dists.mel[, "CS003603_melChighW"][c("CAM040043_melChighE")],
                                         #   dists.mel[, "14N004_melEhighW"][c("CAM041101_melEhighE")],
                                         #   dists.mel[, "CAM040043_melChighE"][c("CAM008521_timEhighE")],
                                         #   dists.mel[, "CAM041101_melEhighE"][c("CAM008521_timEhighE")],
                                         #   dists.mel[, "CS003603_melChighW"][c("CAM040419_cydEhighW")],
                                         #   dists.mel[, "14N004_melEhighW"][c("CAM040419_cydEhighW")]
                                         # ),
                                         donor.p3.country= paste(donor.p3, country.p2, sep="_"),
                                         side.country.p2= paste(side.p2, country.p2, sep="_"),
                                         type.p3=c(rep("within.species", 4), rep("closely.related", 4)),
                                         type.p3.side.country.p2=paste(type.p3, side.country.p2, sep="_" ),
                                         sp.p3=substr(donor.p3, 0,3),
                                         sp.side.P3 = paste(sp.p3,  toupper(substr(donor.p3, 5,5)), sep = "")
); mel.shdr.SAV.other.sp.summ.percs

#### ks tests #####
# null hypothesis that the true distribution function of x is  not greater than the distribution function of y
# alternative hypothesis is that dist of x is greater than dist of y
ks.test(shdr.mel.west.outlier.df$fdm.max_fd_melCvlowW_melChighW_melChighE_hecale, abs(shdr.mel.west.outlier.df$fdm.min_fd_melCvlowW_melChighW_melChighE_hecale), alternative = "less" ) # D^+ = 0.27586, p-value = 0.01211
ks.test(shdr.mel.west.outlier.df$fdm.max_fd_melElowW_melEhighW_melEhighE_hecale, abs(shdr.mel.west.outlier.df$fdm.min_fd_melElowW_melEhighW_melEhighE_hecale), alternative = "less" ) # D^+ = 0.086207, p-value = 0.6498


ks.test(shdr.mel.east.outlier.df$fdm.max_fd_melCvlowE_melChighE_melChighW_hecale, abs(shdr.mel.east.outlier.df$fdm.min_fd_melCvlowE_melChighE_melChighW_hecale), alternative = "less" )
ks.test(shdr.mel.east.outlier.df$fdm.max_fd_melEvlowE_melEhighE_melEhighW_hecale, abs(shdr.mel.east.outlier.df$fdm.min_fd_melEvlowE_melEhighE_melEhighW_hecale), alternative = "less" )

# across sp
ks.test(shdr.mel.east.outlier.df$fdm.max_fd_melCvlowW_melChighW_cydEhighW_hecale, abs(shdr.mel.east.outlier.df$fdm.min_fd_melCvlowW_melChighW_cydEhighW_hecale), alternative = "less" ) 
ks.test(shdr.mel.east.outlier.df$fdm.max_fd_melElowW_melEhighW_cydEhighW_hecale, abs(shdr.mel.east.outlier.df$fdm.min_fd_melElowW_melEhighW_cydEhighW_hecale), alternative = "less" )


ks.test(shdr.mel.east.outlier.df$fdm.max_fd_melCvlowE_melChighE_timEhighE_hecale, abs(shdr.mel.east.outlier.df$fdm.min_fd_melCvlowE_melChighE_timEhighE_hecale), alternative = "less" ) 
ks.test(shdr.mel.east.outlier.df$fdm.max_fd_melEvlowE_melEhighE_timEhighE_hecale, abs(shdr.mel.east.outlier.df$fdm.min_fd_melEvlowE_melEhighE_timEhighE_hecale), alternative = "less" )



#### plot ####
mel.shdr.SAV.other.sp.summ.percs$type.p3 <- factor(mel.shdr.SAV.other.sp.summ.percs$type.p3, levels=c("within.species", "closely.related" ,  "distantly.related"))
mel.shdr.SAV.other.sp.summ.percs$side.p2 <- factor(mel.shdr.SAV.other.sp.summ.percs$side.p2, levels=c("P2 = mel. high West", "P2 = mel. high East"))

ggplot(data=subset(mel.shdr.SAV.other.sp.summ.percs ), 
       aes(x=type.p3 , y=percentage.shdr.with.excess.allele.sharing,  group=side.country.p2, fill=type.p3,  shape=country.p2))+
  #geom_line(inherit.aes = F, color= "black", data=subset(mel.shdr.SAV.other.sp.summ, donor.p3!="mel.west.high" & donor.p3!="mel.east.high" ), aes(x=type.p3 , y=percentage.shdr.with.excess.allele.sharing, group=side.country.p2, linetype=country.p2), size=1)+
  geom_text_repel(aes(label=sp.side.P3), fill="transparent", size=6, segment.size  = 0, segment.color = "grey50",nudge_x       = c(0.4, 0.4))+
  geom_point( size=5, stroke=1, width = 0.1)  +
  scale_shape_manual(values=c(24, 25))+
  scale_fill_manual(values=c("blue", "#eb8055ff","#efe350ff"  ))+
  #scale_color_manual(values=c("black",  "black"))+
  #geom_text(aes(label=donor.p3),hjust=-0.2, vjust=0) +
  #geom_point(inherit.aes = F, data=subset(mel.shdr.SAV.other.sp.summ, (donor.p3=="mel.west.high" | donor.p3=="mel.east.high")  ), aes(x=type.p3 , y=percentage.shdr.with.excess.allele.sharing,  group=side.country.p2),
  #           shape=24, fill="blue", size=5, stroke=1)+
  #geom_text(inherit.aes = F, data=subset(mel.shdr.SAV.other.sp.summ, (donor.p3=="mel.west.high" | donor.p3=="mel.east.high" )&country.p2=="col" ), aes(x=dist.p2.p3, y=percentage.shdr.with.excess.allele.sharing,  group=side.country.p2,  label=donor.p3),hjust=0, vjust=0) +
  #scale_linetype_manual(values=c("solid", "11"))+
  facet_wrap(~side.p2, ncol = 1)+
  xlab("Phylogenetic distance: P2 - P3")+ ylab("% SHDR with evidence for excess allele sharing")+
  scale_y_continuous(expand = c(0.1,0.1), breaks=c(0,20,40), limits = c(0,55))+
  scale_x_discrete(expand = c(0.5,0.3), labels=c("Within species\nbut opposite sides of Andes","close","distant"),position = "top")+
  #scale_x_continuous(expand = c(0.002,0.002))+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=.8),legend.position="none", axis.ticks.length=unit(0.15, "cm"),
        plot.margin = unit(c(1,1,1,1), "lines"),axis.text.x = element_blank(), axis.text.y = element_text(size=14), 
        axis.title.x =  element_blank(),axis.title.y = element_text(size=16), panel.spacing = unit(1, "lines"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),  strip.text = element_blank(), strip.background = element_blank(),
        panel.background = element_rect(fill = "transparent", color = NA), plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
        legend.background = element_rect(fill = "transparent", color = NA), legend.box.background = element_rect(fill = "transparent", color = NA), legend.key = element_rect(fill = "transparent", color = NA),
        strip.background.y = element_blank(), rect = element_rect(fill = "transparent"))


ggsave("plots/joana.fdm/fdm.perc.shdrs.outliers.mel.other.sp.new.png", height = 7, width = 4, bg="transparent")





################################################################################################################ fig 6 new version, ######################################
############################################################## fig 6 C- mean fdm max per SHDR ######################################

############################ 2.6 era fdm per SHDR barplot  ############################
names(shdr.era.east.outlier.df)

era.shdr.SAV.other.sp.summ <- data_frame(donor.p3=c(rep("era.west.high", 2), rep("era.east.high", 2), rep("him.east.high", 2),  rep("tel.east.high", 2),  rep("cly.east.high", 2),  rep("cly.west.high", 2)),
                                         recipient.p2=c(rep("era.west.high", 2), rep("era.east.high", 2), rep("him.east.high", 2),  rep("tel.east.high", 2),  rep("cly.east.high", 2),  rep("cly.west.high", 2)),
                                         country.p2= c(rep(c("col", "ecu" ), 6)),
                                         side.p2= c(rep("P2 = era. high East", 2), rep("P2 = era. high West", 2), rep("P2 = era. high East", 6), rep("P2 = era. high West", 2)),
                                         pop.p3= c(rep("1.era", 4), rep("2.him", 2), rep("3.tel", 2), rep("4.cly", 4)),
                                         min.fdm.per.shdr=c(
                                           mean(shdr.era.east.outlier.df$fdm.min_fd_eraCvlowE_eraChighE_eraChighW_eleHW),
                                           mean(shdr.era.east.outlier.df$fdm.min_fd_eraEvlowE_eraEhighE_eraEhighW_eleHW),
                                           mean(shdr.era.west.outlier.df$fdm.min_fd_eraCvlowW_eraChighW_eraChighE_eleHW),
                                           mean(shdr.era.west.outlier.df$fdm.min_fd_eraEvlowW_eraEhighW_eraEhighE_eleHW),
                                           mean(shdr.era.east.outlier.df$fdm.min_fd_eraCvlowE_eraChighE_himHE_eleHW),
                                           mean(shdr.era.east.outlier.df$fdm.min_fd_eraEvlowE_eraEhighE_himHE_eleHW),
                                           mean(shdr.era.east.outlier.df$fdm.min_fd_eraCvlowE_eraChighE_telhighE_eleHW),
                                           mean(shdr.era.east.outlier.df$fdm.min_fd_eraEvlowE_eraEhighE_telhighE_eleHW),
                                           mean(shdr.era.east.outlier.df$fdm.min_fd_eraCvlowE_eraChighE_clyEhighE_eleHW),
                                           mean(shdr.era.east.outlier.df$fdm.min_fd_eraEvlowE_eraEhighE_clyEhighE_eleHW),
                                           mean(shdr.era.west.outlier.df$fdm.min_fd_eraCvlowW_eraChighW_clyEhighW_eleHW),
                                           mean(shdr.era.west.outlier.df$fdm.min_fd_eraEvlowW_eraEhighW_clyEhighW_eleHW)),
                                         max.fdm.per.shdr=c(
                                           mean(shdr.era.east.outlier.df$fdm.max_fd_eraCvlowE_eraChighE_eraChighW_eleHW),
                                           mean(shdr.era.east.outlier.df$fdm.max_fd_eraEvlowE_eraEhighE_eraEhighW_eleHW),
                                           mean(shdr.era.west.outlier.df$fdm.max_fd_eraCvlowW_eraChighW_eraChighE_eleHW),
                                           mean(shdr.era.west.outlier.df$fdm.max_fd_eraEvlowW_eraEhighW_eraEhighE_eleHW),
                                           mean(shdr.era.east.outlier.df$fdm.max_fd_eraCvlowE_eraChighE_himHE_eleHW),
                                           mean(shdr.era.east.outlier.df$fdm.max_fd_eraEvlowE_eraEhighE_himHE_eleHW),
                                           mean(shdr.era.east.outlier.df$fdm.max_fd_eraCvlowE_eraChighE_telhighE_eleHW),
                                           mean(shdr.era.east.outlier.df$fdm.max_fd_eraEvlowE_eraEhighE_telhighE_eleHW),
                                           mean(shdr.era.east.outlier.df$fdm.max_fd_eraCvlowE_eraChighE_clyEhighE_eleHW),
                                           mean(shdr.era.east.outlier.df$fdm.max_fd_eraEvlowE_eraEhighE_clyEhighE_eleHW),
                                           mean(shdr.era.west.outlier.df$fdm.max_fd_eraCvlowW_eraChighW_clyEhighW_eleHW),
                                           mean(shdr.era.west.outlier.df$fdm.max_fd_eraEvlowW_eraEhighW_clyEhighW_eleHW)),
                                         se.min.fdm.per.shdr=c(
                                           se(shdr.era.east.outlier.df$fdm.min_fd_eraCvlowE_eraChighE_eraChighW_eleHW),
                                           se(shdr.era.east.outlier.df$fdm.min_fd_eraEvlowE_eraEhighE_eraEhighW_eleHW),
                                           se(shdr.era.west.outlier.df$fdm.min_fd_eraCvlowW_eraChighW_eraChighE_eleHW),
                                           se(shdr.era.west.outlier.df$fdm.min_fd_eraEvlowW_eraEhighW_eraEhighE_eleHW),
                                           se(shdr.era.east.outlier.df$fdm.min_fd_eraCvlowE_eraChighE_himHE_eleHW),
                                           se(shdr.era.east.outlier.df$fdm.min_fd_eraEvlowE_eraEhighE_himHE_eleHW),
                                           se(shdr.era.east.outlier.df$fdm.min_fd_eraCvlowE_eraChighE_telhighE_eleHW),
                                           se(shdr.era.east.outlier.df$fdm.min_fd_eraEvlowE_eraEhighE_telhighE_eleHW),
                                           se(shdr.era.east.outlier.df$fdm.min_fd_eraCvlowE_eraChighE_clyEhighE_eleHW),
                                           se(shdr.era.east.outlier.df$fdm.min_fd_eraEvlowE_eraEhighE_clyEhighE_eleHW),
                                           se(shdr.era.west.outlier.df$fdm.min_fd_eraCvlowW_eraChighW_clyEhighW_eleHW),
                                           se(shdr.era.west.outlier.df$fdm.min_fd_eraEvlowW_eraEhighW_clyEhighW_eleHW)),
                                         se.max.fdm.per.shdr=c(
                                           se(shdr.era.east.outlier.df$fdm.max_fd_eraCvlowE_eraChighE_eraChighW_eleHW),
                                           se(shdr.era.east.outlier.df$fdm.max_fd_eraEvlowE_eraEhighE_eraEhighW_eleHW),
                                           se(shdr.era.west.outlier.df$fdm.max_fd_eraCvlowW_eraChighW_eraChighE_eleHW),
                                           se(shdr.era.west.outlier.df$fdm.max_fd_eraEvlowW_eraEhighW_eraEhighE_eleHW),
                                           se(shdr.era.east.outlier.df$fdm.max_fd_eraCvlowE_eraChighE_himHE_eleHW),
                                           se(shdr.era.east.outlier.df$fdm.max_fd_eraEvlowE_eraEhighE_himHE_eleHW),
                                           se(shdr.era.east.outlier.df$fdm.max_fd_eraCvlowE_eraChighE_telhighE_eleHW),
                                           se(shdr.era.east.outlier.df$fdm.max_fd_eraEvlowE_eraEhighE_telhighE_eleHW),
                                           se(shdr.era.east.outlier.df$fdm.max_fd_eraCvlowE_eraChighE_clyEhighE_eleHW),
                                           se(shdr.era.east.outlier.df$fdm.max_fd_eraEvlowE_eraEhighE_clyEhighE_eleHW),
                                           se(shdr.era.west.outlier.df$fdm.max_fd_eraCvlowW_eraChighW_clyEhighW_eleHW),
                                           se(shdr.era.west.outlier.df$fdm.max_fd_eraEvlowW_eraEhighW_clyEhighW_eleHW)),
                                         sims.max.fdm.mean=c(
                                           shdr.sims.fdm.era[shdr.sims.fdm.era$combo=="fd_eraCvlowE_eraChighE_eraChighW_eleHW" & shdr.sims.fdm.era$combo.allo.type=="shdr.east",]$max.fdm.mean,
                                           shdr.sims.fdm.era[shdr.sims.fdm.era$combo=="fd_eraEvlowE_eraEhighE_eraEhighW_eleHW" & shdr.sims.fdm.era$combo.allo.type=="shdr.east",]$max.fdm.mean,
                                           shdr.sims.fdm.era[shdr.sims.fdm.era$combo=="fd_eraCvlowW_eraChighW_eraChighE_eleHW" & shdr.sims.fdm.era$combo.allo.type=="shdr.west",]$max.fdm.mean,
                                           shdr.sims.fdm.era[shdr.sims.fdm.era$combo=="fd_eraEvlowW_eraEhighW_eraEhighE_eleHW" & shdr.sims.fdm.era$combo.allo.type=="shdr.west",]$max.fdm.mean,
                                           shdr.sims.fdm.era[shdr.sims.fdm.era$combo=="fd_eraCvlowE_eraChighE_himHE_eleHW" & shdr.sims.fdm.era$combo.allo.type==    "shdr.east",]$max.fdm.mean,
                                           shdr.sims.fdm.era[shdr.sims.fdm.era$combo=="fd_eraEvlowE_eraEhighE_himHE_eleHW" & shdr.sims.fdm.era$combo.allo.type==    "shdr.east",]$max.fdm.mean,
                                           shdr.sims.fdm.era[shdr.sims.fdm.era$combo=="fd_eraCvlowE_eraChighE_telhighE_eleHW" & shdr.sims.fdm.era$combo.allo.type== "shdr.east",]$max.fdm.mean,
                                           shdr.sims.fdm.era[shdr.sims.fdm.era$combo=="fd_eraEvlowE_eraEhighE_telhighE_eleHW" & shdr.sims.fdm.era$combo.allo.type== "shdr.east",]$max.fdm.mean,
                                           shdr.sims.fdm.era[shdr.sims.fdm.era$combo=="fd_eraCvlowE_eraChighE_clyEhighE_eleHW" & shdr.sims.fdm.era$combo.allo.type=="shdr.east",]$max.fdm.mean,
                                           shdr.sims.fdm.era[shdr.sims.fdm.era$combo=="fd_eraEvlowE_eraEhighE_clyEhighE_eleHW" & shdr.sims.fdm.era$combo.allo.type=="shdr.east",]$max.fdm.mean,
                                           shdr.sims.fdm.era[shdr.sims.fdm.era$combo=="fd_eraCvlowW_eraChighW_clyEhighW_eleHW" & shdr.sims.fdm.era$combo.allo.type=="shdr.west",]$max.fdm.mean,
                                           shdr.sims.fdm.era[shdr.sims.fdm.era$combo=="fd_eraEvlowW_eraEhighW_clyEhighW_eleHW" & shdr.sims.fdm.era$combo.allo.type=="shdr.west",]$max.fdm.mean),
                                         sims.max.fdm.se=c(
                                           shdr.sims.fdm.era[shdr.sims.fdm.era$combo=="fd_eraCvlowE_eraChighE_eraChighW_eleHW" & shdr.sims.fdm.era$combo.allo.type=="shdr.east",]$max.fdm.se,
                                           shdr.sims.fdm.era[shdr.sims.fdm.era$combo=="fd_eraEvlowE_eraEhighE_eraEhighW_eleHW" & shdr.sims.fdm.era$combo.allo.type=="shdr.east",]$max.fdm.se,
                                           shdr.sims.fdm.era[shdr.sims.fdm.era$combo=="fd_eraCvlowW_eraChighW_eraChighE_eleHW" & shdr.sims.fdm.era$combo.allo.type=="shdr.west",]$max.fdm.se,
                                           shdr.sims.fdm.era[shdr.sims.fdm.era$combo=="fd_eraEvlowW_eraEhighW_eraEhighE_eleHW" & shdr.sims.fdm.era$combo.allo.type=="shdr.west",]$max.fdm.se,
                                           shdr.sims.fdm.era[shdr.sims.fdm.era$combo=="fd_eraCvlowE_eraChighE_himHE_eleHW" & shdr.sims.fdm.era$combo.allo.type==    "shdr.east",]$max.fdm.se,
                                           shdr.sims.fdm.era[shdr.sims.fdm.era$combo=="fd_eraEvlowE_eraEhighE_himHE_eleHW" & shdr.sims.fdm.era$combo.allo.type==    "shdr.east",]$max.fdm.se,
                                           shdr.sims.fdm.era[shdr.sims.fdm.era$combo=="fd_eraCvlowE_eraChighE_telhighE_eleHW" & shdr.sims.fdm.era$combo.allo.type== "shdr.east",]$max.fdm.se,
                                           shdr.sims.fdm.era[shdr.sims.fdm.era$combo=="fd_eraEvlowE_eraEhighE_telhighE_eleHW" & shdr.sims.fdm.era$combo.allo.type== "shdr.east",]$max.fdm.se,
                                           shdr.sims.fdm.era[shdr.sims.fdm.era$combo=="fd_eraCvlowE_eraChighE_clyEhighE_eleHW" & shdr.sims.fdm.era$combo.allo.type=="shdr.east",]$max.fdm.se,
                                           shdr.sims.fdm.era[shdr.sims.fdm.era$combo=="fd_eraEvlowE_eraEhighE_clyEhighE_eleHW" & shdr.sims.fdm.era$combo.allo.type=="shdr.east",]$max.fdm.se,
                                           shdr.sims.fdm.era[shdr.sims.fdm.era$combo=="fd_eraCvlowW_eraChighW_clyEhighW_eleHW" & shdr.sims.fdm.era$combo.allo.type=="shdr.west",]$max.fdm.se,
                                           shdr.sims.fdm.era[shdr.sims.fdm.era$combo=="fd_eraEvlowW_eraEhighW_clyEhighW_eleHW" & shdr.sims.fdm.era$combo.allo.type=="shdr.west",]$max.fdm.se),
                                         donor.p3.country= paste(donor.p3, country.p2, sep="_"),
                                         side.country.p2= paste(side.p2, country.p2, sep="_"),
                                         type.p3=c(rep("within.species", 4), rep("closely.related", 2), rep("distantly.related",6)),
                                         type.p3.side.country.p2=paste(type.p3, side.country.p2, sep="_" ),
                                         sp.p3=substr(donor.p3, 0,3),
                                         sp.side.P3 = paste(sp.p3,  toupper(substr(donor.p3, 5,5)), sep = "")
); era.shdr.SAV.other.sp.summ

era.shdr.SAV.other.sp.summ$type.p3 <- factor(era.shdr.SAV.other.sp.summ$type.p3, levels=c("within.species", "closely.related" ,  "distantly.related"))
era.shdr.SAV.other.sp.summ$side.p2 <- factor(era.shdr.SAV.other.sp.summ$side.p2, levels=c("P2 = era. high West", "P2 = era. high East"))


################## plot option1- includes minimum fdm ##################
era.shdr.SAV.other.sp.summ.summ.long <- gather(era.shdr.SAV.other.sp.summ[,c("sp.side.P3", "recipient.p2", "type.p3", "min.fdm.per.shdr", "max.fdm.per.shdr", "country.p2", "se.min.fdm.per.shdr", "se.max.fdm.per.shdr", "side.p2")], type.min.max,min.max, -sp.side.P3, -recipient.p2, -type.p3, -country.p2, -side.p2 ); era.shdr.SAV.other.sp.summ.summ.long 
era.shdr.SAV.other.sp.summ.summ.long$side.p2 <- substr(era.shdr.SAV.other.sp.summ.summ.long$side.p2, 16,16)
era.shdr.SAV.other.sp.summ.summ.long$sp.side.P1.P2 <- paste("era", era.shdr.SAV.other.sp.summ.summ.long$side.p2, sep=""); era.shdr.SAV.other.sp.summ.summ.long$sp.side.P1.P2

unique(era.shdr.SAV.other.sp.summ.summ.long$type.min.max)
era.shdr.SAV.other.sp.summ.summ.long$se.mean <- if_else(substr(era.shdr.SAV.other.sp.summ.summ.long$type.min.max, 0,2 )=="se", "se", "mean")  ; era.shdr.SAV.other.sp.summ.summ.long
era.shdr.SAV.other.sp.summ.summ.long$type.min.max<- substr(era.shdr.SAV.other.sp.summ.summ.long$type.min.max, nchar(era.shdr.SAV.other.sp.summ.summ.long$type.min.max)-15,nchar(era.shdr.SAV.other.sp.summ.summ.long$type.min.max)); era.shdr.SAV.other.sp.summ.summ.long
era.shdr.SAV.other.sp.summ.summ.long$p1.low.p2.high <- if_else(era.shdr.SAV.other.sp.summ.summ.long$type.min.max=="min.fdm.per.shdr", "low", "high"); era.shdr.SAV.other.sp.summ.summ.long$p1.low.p2.high
era.shdr.SAV.other.sp.summ.summ.long$p1.low.p2.high.country <- paste(era.shdr.SAV.other.sp.summ.summ.long$p1.low.p2.high, era.shdr.SAV.other.sp.summ.summ.long$country.p2, sep="_"); era.shdr.SAV.other.sp.summ.summ.long$p1.low.p2.high.country
era.shdr.SAV.other.sp.summ.summ.long.wide <- spread(era.shdr.SAV.other.sp.summ.summ.long,  se.mean, min.max); era.shdr.SAV.other.sp.summ.summ.long.wide


the_order <- rev(c("eraW","eraE"  , "himE","telE", "clyE","clyW"))
era.fdm.min.max.dotplot <-era.shdr.SAV.other.sp.summ.summ.long.wide %>% 
  ggplot(aes(x = sp.side.P3, y = mean, group = recipient.p2, shape=p1.low.p2.high.country)) +
  #geom_bar(stat = "identity", width = 0.75) +
  #geom_errorbar(aes(ymin=mean-se, ymax=mean+se), colour="black",  width=0, alpha=0.7) +
  scale_shape_manual(values=c(17,24,19,21))+
  geom_point(size=3, fill="white",alpha=0.99)+
  coord_flip() +
  scale_x_discrete(limits = the_order, position = "top") +
  scale_y_continuous(limits = c(-0.1,0.1), breaks = seq(-1, 1, 0.1),  labels = round(seq(-1, 1, 0.1), digits = 2) ) +
  labs(x = "", y = "") +
  geom_hline(yintercept = 0, color="black")+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=.8),legend.position="none", axis.ticks.length=unit(0.15, "cm"),
        plot.margin = unit(c(1,1,1,1), "lines"),axis.text.x = element_text(size=10), axis.text.y = element_text(size=14, face = "bold", color=c( "#efe350ff","#efe350ff", "#efe350ff",  "#eb8055ff","blue", "blue")),
        axis.title = element_text(size=16), panel.spacing = unit(1, "lines"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),  strip.text = element_blank(), strip.background = element_blank(),
        panel.background = element_rect(fill = "transparent", color = NA), plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
        legend.background = element_rect(fill = "transparent", color = NA), legend.box.background = element_rect(fill = "transparent", color = NA), legend.key = element_rect(fill = "transparent", color = NA),
        strip.background.y = element_blank(), rect = element_rect(fill = "transparent")); era.fdm.min.max.dotplot

ggsave("plots/joana.fdm/dotplot.fdm.era.png", height = 3.5, width = 4)

################## plot option2- mean fdm max + simulated ##################

names(era.shdr.SAV.other.sp.summ)
levels(era.shdr.SAV.other.sp.summ$country.p2)
era.shdr.SAV.other.sp.summ$country.p2 <- factor(era.shdr.SAV.other.sp.summ$country.p2, levels = c("col","ecu"))



era.the_order <- rev(c("eraW","eraE"  , "himE","telE", "clyE","clyW")); era.the_order_p2 <- rev(c("eraE","eraW"  , "eraE","eraE", "eraE","eraW"))
era.fdm.max.sims.dotplot <- era.shdr.SAV.other.sp.summ %>% 
  ggplot(aes(x = sp.side.P3, y = max.fdm.per.shdr, group = recipient.p2, shape=country.p2)) +
  # add simulations, SE are tiny
  geom_point(inherit.aes = F, data= era.shdr.SAV.other.sp.summ, aes(x = sp.side.P3, y = sims.max.fdm.mean, group = recipient.p2, shape=country.p2), 
             size=3, fill="white",alpha=0.99, color="grey", position = position_nudge(c(0.15, -0.15)))+
  
  geom_errorbar(aes(ymin=max.fdm.per.shdr- se.max.fdm.per.shdr, ymax=max.fdm.per.shdr +se.max.fdm.per.shdr), colour="black",  width=0, alpha=0.7, position = position_nudge(c(0.15, -0.15))) +
  geom_point(size=4, fill="white",alpha=0.99, position = position_nudge(c( 0.15, -0.15)))+
  coord_flip() +
  scale_shape_manual(values=c(17,24))+
  scale_x_discrete(limits = era.the_order, position = "bottom") +
  scale_y_continuous(limits = c(0,0.1), breaks = seq(-1, 1, 0.05),  labels = round(seq(-1, 1, 0.05), digits = 2), expand = c(0.001,0) ) +
  guides(y.sec = guide_axis_label_trans(~paste(era.the_order_p2)))+
  labs(x = "", y = "") +
  guides(shape = guide_legend(reverse = TRUE)) +
  #geom_hline(yintercept = 0, color="black")+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=.8),legend.position="none", axis.ticks.length=unit(0.15, "cm"),
        plot.margin =  unit(c(0,0,0,0), "lines"),axis.text.x = element_text(size=10), 
        axis.text.y.left = element_text(size=18, face = "bold", color=c( "#efe350ff","#efe350ff", "#efe350ff",  "#eb8055ff","blue", "blue")),
        axis.text.y.right = element_blank(),
        axis.title = element_blank(),
         panel.spacing = unit(1, "lines"),panel.grid.major.y = element_line(color = "grey90"),
        panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(),  strip.text = element_blank(), strip.background = element_blank(),
        panel.background = element_rect(fill = "transparent", color = NA), plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
        legend.background = element_rect(fill = "transparent", color = NA), legend.box.background = element_rect(fill = "transparent", color = NA), legend.key = element_rect(fill = "transparent", color = NA),
        strip.background.y = element_blank(), rect = element_rect(fill = "transparent")); era.fdm.max.sims.dotplot

#ggsave("plots/joana.fdm/dotplot.fdm.era.png", height = 3.5, width = 4)

################## plot option2- % shdr ##################
head(era.shdr.SAV.other.sp.summ.percs  )
names(era.shdr.SAV.other.sp.summ.percs)
era.the_order <- rev(c("eraW","eraE"  , "himE","telE", "clyE","clyW"))
era.the_order_p2 <- rev(c("eraE","eraW"  , "eraE","eraE", "eraE","eraW"))
era.perc.shdr.outlier.dotplot <- era.shdr.SAV.other.sp.summ.percs   %>% 
  ggplot(aes(x=sp.side.P3 , y=percentage.shdr.with.excess.allele.sharing,   group = recipient.p2, shape=side.country.p2, fill=country.p2)) +
  geom_bar(inherit.aes = F, aes(x=sp.side.P3 , y=percentage.shdr.with.excess.allele.sharing, shape=side.country.p2, fill=country.p2),
           alpha=0.99, stat="identity", position=position_dodge(c(-.7)), width = 0.7, color="black")+
  scale_fill_manual(values = c("black","white"))+
  coord_flip() +
  scale_x_discrete(limits = era.the_order, position = "bottom") +
  scale_y_continuous(limits = c(0,55), breaks = seq(0, 50, 25),  labels = round(seq(0, 50, 25), digits = 2), expand = c(0.001,0) ) +
  guides(y.sec = guide_axis_label_trans(~paste(era.the_order_p2)))+
  labs(x = "", y = "") +
  #geom_hline(yintercept = 0, color="black")+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=.8),legend.position="none", axis.ticks.length=unit(0.15, "cm"),
        plot.margin = unit(c(0,0,0,0), "lines"),axis.text.x = element_text(size=10), 
        axis.text.y.left = element_blank(),
        axis.text.y.right = element_text(size=18, face = "bold", color="black"),
        axis.title = element_blank(),
         panel.spacing = unit(1, "lines"),panel.grid.major.y = element_line(color = "grey90"),
        panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(),  strip.text = element_blank(), strip.background = element_blank(),
        panel.background = element_rect(fill = "transparent", color = NA), plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
        legend.background = element_rect(fill = "transparent", color = NA), legend.box.background = element_rect(fill = "transparent", color = NA), legend.key = element_rect(fill = "transparent", color = NA),
        strip.background.y = element_blank(), rect = element_rect(fill = "transparent")); era.perc.shdr.outlier.dotplot


plot_grid(era.fdm.max.sims.dotplot, era.perc.shdr.outlier.dotplot)
ggsave("plots/joana.fdm/dotplot.fdm.max.sims.percs.era.png", height = 3.7, width = 5)

############################ 2.7 mel fdm per SHDR barplot  ############################
names(shdr.mel.east.outlier.df)
se <- function(x) sd(x)/sqrt(length(x))

mel.shdr.SAV.other.sp.summ <- data_frame(donor.p3=c( rep("mel.west.high", 2) ,rep("mel.east.high", 2), rep("tim.east.high", 2),  rep("cyd.west.high", 2)),
                                         country.p2= c(rep(c("col", "ecu" ), 4)),
                                         side.p2= c(rep("P2 = mel. high East", 2), rep("P2 = mel. high West", 2), rep("P2 = mel. high East", 2), rep("P2 = mel. high West", 2)),
                                         
                                         dist= c(rep("2.mel", 4), rep("1.tim.E", 2), rep("3.cyd.W", 2)),
                                         min.fdm.per.shdr= c(
                                           mean(shdr.mel.east.outlier.df$fdm.min_fd_melCvlowE_melChighE_melChighW_hecale),
                                           mean(shdr.mel.east.outlier.df$fdm.min_fd_melEvlowE_melEhighE_melEhighW_hecale),
                                           mean(shdr.mel.west.outlier.df$fdm.min_fd_melCvlowW_melChighW_melChighE_hecale),
                                           mean(shdr.mel.west.outlier.df$fdm.min_fd_melElowW_melEhighW_melEhighE_hecale),
                                           mean(shdr.mel.east.outlier.df$fdm.min_fd_melCvlowE_melChighE_timEhighE_hecale),
                                           mean(shdr.mel.east.outlier.df$fdm.min_fd_melEvlowE_melEhighE_timEhighE_hecale),
                                           mean(shdr.mel.west.outlier.df$fdm.min_fd_melCvlowW_melChighW_cydEhighW_hecale),
                                           mean(shdr.mel.west.outlier.df$fdm.min_fd_melElowW_melEhighW_cydEhighW_hecale)),
                                         max.fdm.per.shdr= c(
                                           mean(shdr.mel.east.outlier.df$fdm.max_fd_melCvlowE_melChighE_melChighW_hecale),
                                           mean(shdr.mel.east.outlier.df$fdm.max_fd_melEvlowE_melEhighE_melEhighW_hecale),
                                           mean(shdr.mel.west.outlier.df$fdm.max_fd_melCvlowW_melChighW_melChighE_hecale),
                                           mean(shdr.mel.west.outlier.df$fdm.max_fd_melElowW_melEhighW_melEhighE_hecale),
                                           mean(shdr.mel.east.outlier.df$fdm.max_fd_melCvlowE_melChighE_timEhighE_hecale),
                                           mean(shdr.mel.east.outlier.df$fdm.max_fd_melEvlowE_melEhighE_timEhighE_hecale),
                                           mean(shdr.mel.west.outlier.df$fdm.max_fd_melCvlowW_melChighW_cydEhighW_hecale),
                                           mean(shdr.mel.west.outlier.df$fdm.max_fd_melElowW_melEhighW_cydEhighW_hecale)),
                                         se.min.fdm.per.shdr= c(
                                           se(shdr.mel.east.outlier.df$fdm.min_fd_melCvlowE_melChighE_melChighW_hecale),
                                           se(shdr.mel.east.outlier.df$fdm.min_fd_melEvlowE_melEhighE_melEhighW_hecale),
                                           se(shdr.mel.west.outlier.df$fdm.min_fd_melCvlowW_melChighW_melChighE_hecale),
                                           se(shdr.mel.west.outlier.df$fdm.min_fd_melElowW_melEhighW_melEhighE_hecale),
                                           se(shdr.mel.east.outlier.df$fdm.min_fd_melCvlowE_melChighE_timEhighE_hecale),
                                           se(shdr.mel.east.outlier.df$fdm.min_fd_melEvlowE_melEhighE_timEhighE_hecale),
                                           se(shdr.mel.west.outlier.df$fdm.min_fd_melCvlowW_melChighW_cydEhighW_hecale),
                                           se(shdr.mel.west.outlier.df$fdm.min_fd_melElowW_melEhighW_cydEhighW_hecale)),
                                         se.max.fdm.per.shdr= c(
                                           se(shdr.mel.east.outlier.df$fdm.max_fd_melCvlowE_melChighE_melChighW_hecale),
                                           se(shdr.mel.east.outlier.df$fdm.max_fd_melEvlowE_melEhighE_melEhighW_hecale),
                                           se(shdr.mel.west.outlier.df$fdm.max_fd_melCvlowW_melChighW_melChighE_hecale),
                                           se(shdr.mel.west.outlier.df$fdm.max_fd_melElowW_melEhighW_melEhighE_hecale),
                                           se(shdr.mel.east.outlier.df$fdm.max_fd_melCvlowE_melChighE_timEhighE_hecale),
                                           se(shdr.mel.east.outlier.df$fdm.max_fd_melEvlowE_melEhighE_timEhighE_hecale),
                                           se(shdr.mel.west.outlier.df$fdm.max_fd_melCvlowW_melChighW_cydEhighW_hecale),
                                           se(shdr.mel.west.outlier.df$fdm.max_fd_melElowW_melEhighW_cydEhighW_hecale)),
                                         sims.max.fdm.mean=c(
                                           shdr.sims.fdm.mel[shdr.sims.fdm.mel$combo=="fd_melCvlowE_melChighE_melChighW_hecale" & shdr.sims.fdm.mel$combo.para.allo.type=="shdr.para.east",]$max.fdm.mean,
                                           shdr.sims.fdm.mel[shdr.sims.fdm.mel$combo=="fd_melEvlowE_melEhighE_melEhighW_hecale" & shdr.sims.fdm.mel$combo.para.allo.type=="shdr.para.east",]$max.fdm.mean,
                                           shdr.sims.fdm.mel[shdr.sims.fdm.mel$combo=="fd_melCvlowW_melChighW_melChighE_hecale" & shdr.sims.fdm.mel$combo.para.allo.type=="shdr.para.west",]$max.fdm.mean,
                                           shdr.sims.fdm.mel[shdr.sims.fdm.mel$combo=="fd_melElowW_melEhighW_melEhighE_hecale" &  shdr.sims.fdm.mel$combo.para.allo.type== "shdr.para.west",]$max.fdm.mean,
                                           shdr.sims.fdm.mel[shdr.sims.fdm.mel$combo=="fd_melCvlowE_melChighE_timEhighE_hecale" & shdr.sims.fdm.mel$combo.para.allo.type=="shdr.para.east",]$max.fdm.mean,
                                           shdr.sims.fdm.mel[shdr.sims.fdm.mel$combo=="fd_melEvlowE_melEhighE_timEhighE_hecale" & shdr.sims.fdm.mel$combo.para.allo.type=="shdr.para.east",]$max.fdm.mean,
                                           shdr.sims.fdm.mel[shdr.sims.fdm.mel$combo=="fd_melCvlowW_melChighW_cydEhighW_hecale" & shdr.sims.fdm.mel$combo.para.allo.type=="shdr.para.west",]$max.fdm.mean,
                                           shdr.sims.fdm.mel[shdr.sims.fdm.mel$combo=="fd_melElowW_melEhighW_cydEhighW_hecale" &  shdr.sims.fdm.mel$combo.para.allo.type== "shdr.para.west",]$max.fdm.mean),
                                         sims.max.fdm.se=c(
                                           shdr.sims.fdm.mel[shdr.sims.fdm.mel$combo=="fd_melCvlowE_melChighE_melChighW_hecale" & shdr.sims.fdm.mel$combo.para.allo.type=="shdr.para.east",]$max.fdm.se,
                                           shdr.sims.fdm.mel[shdr.sims.fdm.mel$combo=="fd_melEvlowE_melEhighE_melEhighW_hecale" & shdr.sims.fdm.mel$combo.para.allo.type=="shdr.para.east",]$max.fdm.se,
                                           shdr.sims.fdm.mel[shdr.sims.fdm.mel$combo=="fd_melCvlowW_melChighW_melChighE_hecale" & shdr.sims.fdm.mel$combo.para.allo.type=="shdr.para.west",]$max.fdm.se,
                                           shdr.sims.fdm.mel[shdr.sims.fdm.mel$combo=="fd_melElowW_melEhighW_melEhighE_hecale" &  shdr.sims.fdm.mel$combo.para.allo.type== "shdr.para.west",]$max.fdm.se,
                                           shdr.sims.fdm.mel[shdr.sims.fdm.mel$combo=="fd_melCvlowE_melChighE_timEhighE_hecale" & shdr.sims.fdm.mel$combo.para.allo.type=="shdr.para.east",]$max.fdm.se,
                                           shdr.sims.fdm.mel[shdr.sims.fdm.mel$combo=="fd_melEvlowE_melEhighE_timEhighE_hecale" & shdr.sims.fdm.mel$combo.para.allo.type=="shdr.para.east",]$max.fdm.se,
                                           shdr.sims.fdm.mel[shdr.sims.fdm.mel$combo=="fd_melCvlowW_melChighW_cydEhighW_hecale" & shdr.sims.fdm.mel$combo.para.allo.type=="shdr.para.west",]$max.fdm.se,
                                           shdr.sims.fdm.mel[shdr.sims.fdm.mel$combo=="fd_melElowW_melEhighW_cydEhighW_hecale" &  shdr.sims.fdm.mel$combo.para.allo.type== "shdr.para.west",]$max.fdm.se),
                                         donor.p3.country= paste(donor.p3, country.p2, sep="_"),
                                         side.country.p2= paste(side.p2, country.p2, sep="_"),
                                         type.p3=c(rep("within.species", 4), rep("closely.related", 4)),
                                         type.p3.side.country.p2=paste(type.p3, side.country.p2, sep="_" ),
                                         sp.p3=substr(donor.p3, 0,3),
                                         sp.side.P3 = paste(sp.p3,  toupper(substr(donor.p3, 5,5)), sep = "")
); mel.shdr.SAV.other.sp.summ



mel.shdr.SAV.other.sp.summ$type.p3 <- factor(mel.shdr.SAV.other.sp.summ$type.p3, levels=c("within.species", "closely.related" ,  "distantly.related"))
mel.shdr.SAV.other.sp.summ$side.p2 <- factor(mel.shdr.SAV.other.sp.summ$side.p2, levels=c("P2 = mel. high West", "P2 = mel. high East"))


mel.shdr.SAV.other.sp.summ.summ <- summarise(group_by(mel.shdr.SAV.other.sp.summ, sp.side.P3, type.p3),
                                             mean.min.fdm=mean(min.fdm.per.shdr),
                                             se.min.fdm=se(min.fdm.per.shdr),
                                             se.mean.min.fdm=mean(se.min.fdm.per.shdr),
                                             mean.max.fdm=mean(max.fdm.per.shdr),
                                             se.max.fdm=se(max.fdm.per.shdr),
                                             se.mean.max.fdm=mean(se.max.fdm.per.shdr)); mel.shdr.SAV.other.sp.summ.summ


#### with one point for colombia and ecu ####
mel.shdr.SAV.other.sp.summ.summ.long <- gather(mel.shdr.SAV.other.sp.summ[,c("sp.side.P3", "type.p3", "min.fdm.per.shdr", "max.fdm.per.shdr", "country.p2", "se.min.fdm.per.shdr", "se.max.fdm.per.shdr", "side.p2")],
                                               type.min.max,min.max, -sp.side.P3, -type.p3, -country.p2, -side.p2 ); mel.shdr.SAV.other.sp.summ.summ.long 
mel.shdr.SAV.other.sp.summ.summ.long$side.p2 <- substr(mel.shdr.SAV.other.sp.summ.summ.long$side.p2, 16,16)
mel.shdr.SAV.other.sp.summ.summ.long$sp.side.P1.P2 <- paste("mel", mel.shdr.SAV.other.sp.summ.summ.long$side.p2, sep=""); mel.shdr.SAV.other.sp.summ.summ.long$sp.side.P1.P2

unique(mel.shdr.SAV.other.sp.summ.summ.long$type.min.max)
mel.shdr.SAV.other.sp.summ.summ.long$se.mean <- if_else(substr(mel.shdr.SAV.other.sp.summ.summ.long$type.min.max, 0,2 )=="se", "se", "mean")  ; mel.shdr.SAV.other.sp.summ.summ.long
mel.shdr.SAV.other.sp.summ.summ.long$type.min.max<- substr(mel.shdr.SAV.other.sp.summ.summ.long$type.min.max, nchar(mel.shdr.SAV.other.sp.summ.summ.long$type.min.max)-15,nchar(mel.shdr.SAV.other.sp.summ.summ.long$type.min.max)); mel.shdr.SAV.other.sp.summ.summ.long
mel.shdr.SAV.other.sp.summ.summ.long$p1.low.p2.high <- if_else(mel.shdr.SAV.other.sp.summ.summ.long$type.min.max=="min.fdm.per.shdr", "low", "high"); mel.shdr.SAV.other.sp.summ.summ.long$p1.low.p2.high
mel.shdr.SAV.other.sp.summ.summ.long$p1.low.p2.high.country <- paste(mel.shdr.SAV.other.sp.summ.summ.long$p1.low.p2.high, mel.shdr.SAV.other.sp.summ.summ.long$country.p2, sep="_"); mel.shdr.SAV.other.sp.summ.summ.long$p1.low.p2.high.country
mel.shdr.SAV.other.sp.summ.summ.long.wide <- spread(mel.shdr.SAV.other.sp.summ.summ.long,  se.mean, min.max); mel.shdr.SAV.other.sp.summ.summ.long.wide


the_order <- rev(c("melW","melE"  , "timE","cydW"))
mel.fdm.min.max.dotplot <- mel.shdr.SAV.other.sp.summ.summ.long.wide %>% 
  ggplot(aes(x = sp.side.P3, y = mean, group = sp.side.P3, shape=p1.low.p2.high.country)) +
  #geom_bar(stat = "identity", width = 0.75) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), colour="black",  width=0, alpha=0.7) +
  scale_shape_manual(values=c(17,24,19,21))+
  geom_point(size=3, fill="white",alpha=0.99)+
  coord_flip() +
  scale_x_discrete(limits = the_order, position = "top") +
  scale_y_continuous(limits = c(-0.15,0.15), breaks = seq(-1, 1, 0.1),  labels = round(seq(-1, 1, 0.1), digits = 2) ) +
  labs(x = "", y = "") +
  geom_hline(yintercept = 0, color="black")+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=.8),legend.position="none", axis.ticks.length=unit(0.15, "cm"),
        plot.margin = unit(c(1,1,1,1), "lines"),axis.text.x = element_text(size=10), axis.text.y = element_text(size=14, face = "bold", color=c( "#eb8055ff","#eb8055ff","blue", "blue")),
        axis.title = element_text(size=16), panel.spacing = unit(1, "lines"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),  strip.text = element_blank(), strip.background = element_blank(),
        panel.background = element_rect(fill = "transparent", color = NA), plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
        legend.background = element_rect(fill = "transparent", color = NA), legend.box.background = element_rect(fill = "transparent", color = NA), legend.key = element_rect(fill = "transparent", color = NA),
        strip.background.y = element_blank(), rect = element_rect(fill = "transparent")); mel.fdm.min.max.dotplot

ggsave("plots/joana.fdm/dotplot.fdm.mel.png", height = 3, width = 4)


################## plot option2- mean fdm max + simulated ##################
# col first, then ecu
levels(mel.shdr.SAV.other.sp.summ$side.country.p2)
mel.shdr.SAV.other.sp.summ$country.p2 <- factor(mel.shdr.SAV.other.sp.summ$country.p2, levels = c("col","ecu"))

library(forcats)
names(mel.shdr.SAV.other.sp.summ)
mel.shdr.SAV.other.sp.summ$country.p2
mel.the_order <- rev(c("melW","melE"  , "timE","cydW"))
mel.the_order_p2 <- rev(c("melE","melW"  , "melE","melW"))
mel.fdm.max.sims.dotplot <- mel.shdr.SAV.other.sp.summ %>% 
  ggplot(aes(x = sp.side.P3, y = max.fdm.per.shdr, group =  donor.p3, shape=country.p2)) +
  # add simulations, SE are tiny
  geom_point(inherit.aes = F, data= mel.shdr.SAV.other.sp.summ, aes(x = sp.side.P3, y = sims.max.fdm.mean, group =  donor.p3, shape=country.p2), 
             size=3, fill="white",alpha=0.99, color="grey", position = position_nudge(c(0.15,-0.15)))+
  
  geom_errorbar(aes(ymin=max.fdm.per.shdr- se.max.fdm.per.shdr, ymax=max.fdm.per.shdr +se.max.fdm.per.shdr), colour="black",  width=0, alpha=0.7, position = position_nudge(c(0.15,-0.15))) +
  scale_shape_manual(values=c(17,24))+
  geom_point(size=4, fill="white",alpha=0.99, position = position_nudge(c(0.15,-0.15)))+
  coord_flip() +
  scale_x_discrete(limits = mel.the_order, position = "bottom") +
  scale_y_continuous(limits = c(0,0.15), breaks = seq(-1, 1, 0.05),  labels = round(seq(-1, 1, 0.05), digits = 2), expand = c(0.001,0) ) +
  guides(y.sec = guide_axis_label_trans(~paste(mel.the_order_p2)))+
  labs(x = "", y = "") +
  #geom_hline(yintercept = 0, color="black")+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=.8),legend.position="none", axis.ticks.length=unit(0.15, "cm"),
        plot.margin =  unit(c(0,0,0,0), "lines"),axis.text.x = element_text(size=10), 
        axis.text.y.left = element_text(size=18, face = "bold", color=c(  "#eb8055ff","#eb8055ff","blue", "blue")),
        axis.text.y.right = element_blank(),
        axis.title = element_blank(),
        panel.spacing = unit(1, "lines"),panel.grid.major.y = element_line(color = "grey90"),
        panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(),  strip.text = element_blank(), strip.background = element_blank(),
        panel.background = element_rect(fill = "transparent", color = NA), plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
        legend.background = element_rect(fill = "transparent", color = NA), legend.box.background = element_rect(fill = "transparent", color = NA), legend.key = element_rect(fill = "transparent", color = NA),
        strip.background.y = element_blank(), rect = element_rect(fill = "transparent")); mel.fdm.max.sims.dotplot

#ggsave("plots/joana.fdm/dotplot.fdm.mel.png", height = 3.5, width = 4)

################## plot option2- % shdr ##################
head(mel.shdr.SAV.other.sp.summ.percs  )
names(mel.shdr.SAV.other.sp.summ.percs)
mel.the_order <- rev(c("melW","melE"  , "timE","cydW"))
mel.the_order_p2 <- rev(c("melE","melW"  , "melE","melW"))

mel.perc.shdr.outlier.dotplot <- mel.shdr.SAV.other.sp.summ.percs   %>% 
  ggplot(aes(x=sp.side.P3 , y=percentage.shdr.with.excess.allele.sharing,   group = recipient.p2, shape=country.p2, fill=country.p2)) +
  geom_bar(inherit.aes = F, aes(x=sp.side.P3 , y=percentage.shdr.with.excess.allele.sharing, shape=country.p2, fill=country.p2),
           alpha=0.99, stat="identity", position=position_dodge(c(-.7)), width = 0.7, color="black")+
  scale_fill_manual(values = c("black","white"))+
  coord_flip() +
  scale_x_discrete(limits = mel.the_order, position = "bottom") +
  scale_y_continuous(limits = c(0,55), breaks = seq(0, 50, 25),  labels = round(seq(0, 50, 25), digits = 2), expand = c(0.001,0) ) +
  guides(y.sec = guide_axis_label_trans(~paste(mel.the_order_p2)))+
  labs(x = "", y = "") +
  #geom_hline(yintercept = 0, color="black")+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=.8),legend.position="none", axis.ticks.length=unit(0.15, "cm"),
        plot.margin = unit(c(0,0,0,0), "lines"),axis.text.x = element_text(size=10), 
        axis.text.y.left = element_blank(),
        axis.text.y.right = element_text(size=18, face = "bold", color="black"),
        axis.title = element_blank(),
        panel.spacing = unit(1, "lines"),panel.grid.major.y = element_line(color = "grey90"),
        panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(),  strip.text = element_blank(), strip.background = element_blank(),
        panel.background = element_rect(fill = "transparent", color = NA), plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
        legend.background = element_rect(fill = "transparent", color = NA), legend.box.background = element_rect(fill = "transparent", color = NA), legend.key = element_rect(fill = "transparent", color = NA),
        strip.background.y = element_blank(), rect = element_rect(fill = "transparent")); mel.perc.shdr.outlier.dotplot


plot_grid(mel.fdm.max.sims.dotplot, mel.perc.shdr.outlier.dotplot)
ggsave("plots/joana.fdm/dotplot.fdm.max.sims.percs.mel.png", height = 2.5, width = 5)




############################################################## fig 6 C- % SHDR above abs min fdm threshold ######################################
# run data loading above
#### plot ####
ggplot(data=subset(era.shdr.SAV.other.sp.summ ), 
       aes(x=type.p3 , y=percentage.shdr.with.excess.allele.sharing,  group=side.country.p2, fill=type.p3,  shape=country.p2))+
  #geom_line(inherit.aes = F, color= "black", data=subset(era.shdr.SAV.other.sp.summ, donor.p3!="era.west.high" & donor.p3!="era.east.high" ), aes(x=type.p3 , y=percentage.shdr.with.excess.allele.sharing, group=side.country.p2, linetype=country.p2), size=1)+
  geom_text_repel(aes(label=sp.side.P3), fill="transparent", size=6, segment.size  = 0.8, segment.color = "grey50",nudge_x       = c(-0.4, 0.5))+
  geom_point( size=5, stroke=1, width = 0.1)  +
  scale_shape_manual(values=c(24, 25))+
  scale_fill_manual(values=c("blue", "#eb8055ff","#efe350ff"  ))+
  #scale_color_manual(values=c("black",  "black"))+
  #geom_text(aes(label=donor.p3),hjust=-0.2, vjust=0) +
  #geom_point(inherit.aes = F, data=subset(era.shdr.SAV.other.sp.summ, (donor.p3=="era.west.high" | donor.p3=="era.east.high")  ), aes(x=type.p3 , y=percentage.shdr.with.excess.allele.sharing,  group=side.country.p2),
  #           shape=24, fill="blue", size=5, stroke=1)+
  #geom_text(inherit.aes = F, data=subset(era.shdr.SAV.other.sp.summ, (donor.p3=="era.west.high" | donor.p3=="era.east.high" )&country.p2=="col" ), aes(x=dist.p2.p3, y=percentage.shdr.with.excess.allele.sharing,  group=side.country.p2,  label=donor.p3),hjust=0, vjust=0) +
  #scale_linetype_manual(values=c("solid", "11"))+
  facet_wrap(~side.p2, ncol = 1)+
  xlab("Phylogenetic distance: P2 - P3")+ ylab("% SHDR with evidence for excess allele sharing")+
  scale_y_continuous(expand = c(0.1,0.1), breaks=c(0,20,40))+
  scale_x_discrete( labels=c("Within species\nbut opposite sides of Andes","close","distant"),position = "top")+
  #scale_x_continuous(expand = c(0.002,0.002))+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=.8),legend.position="none", axis.ticks.length=unit(0.15, "cm"),
        plot.margin = unit(c(1,1,1,1), "lines"),axis.text.x = element_blank(), axis.text.y = element_text(size=14), 
        axis.title.x =  element_blank(),axis.title.y = element_text(size=16), panel.spacing = unit(1, "lines"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),  strip.text = element_blank(), strip.background = element_blank(),
        panel.background = element_rect(fill = "transparent", color = NA), plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
        legend.background = element_rect(fill = "transparent", color = NA), legend.box.background = element_rect(fill = "transparent", color = NA), legend.key = element_rect(fill = "transparent", color = NA),
        strip.background.y = element_blank(), rect = element_rect(fill = "transparent"))


############################ combine fdm min max dotplot ecu/col era/mel ############################
era.mel.tree | (era.fdm.min.max.dotplot/ mel.fdm.min.max.dotplot)

plot_grid(era.mel.tree, plot_grid(era.fdm.min.max.dotplot, mel.fdm.min.max.dotplot, rel_heights = c(1,0.7), ncol = 1),rel_widths = c(1,1) )
ggsave("plots/joana.fdm/dotplot.era.mel.tree.png", height = 8, width = 10)

plot_grid(era.fdm.min.max.dotplot, mel.fdm.min.max.dotplot, rel_heights = c(1,0.8), ncol = 1)
ggsave("plots/joana.fdm/dotplot.era.mel.png", height = 6, width = 4)


plot_grid(era.mel.tree, plot_grid(era.fdm.max.sims.dotplot, mel.fdm.max.sims.dotplot, rel_heights = c(1,0.7), ncol = 1),rel_widths = c(1,1) )
ggsave("plots/joana.fdm/dotplot.era.mel.tree.png", height = 8, width = 10)

plot_grid(era.fdm.min.max.dotplot, mel.fdm.min.max.dotplot, rel_heights = c(1,0.8), ncol = 1)
ggsave("plots/joana.fdm/dotplot.era.mel.png", height = 6, width = 4)


################################################################################################################ fig6C VERSION 1 ########################################################
######################################################## 1. trees ############################
################ tree ((p1, p2), p3, O) ################ 
# https://4va.github.io/biodatasci/r-ggtree.html#labeling_clades
# high blue "#0000FF" , shape 24
# low green "#5EC55B" , shape 22
# low distant green "#23FF06" , shape 21
ggplot()+
  geom_segment(aes(x=0, y=0, xend=1.05, yend=1.05), size=2)+
  geom_segment(aes(x=1, y=1, xend=2, yend=0), size=2)+
  geom_segment(aes(x=0.75, y=0.75, xend=1.5, yend=0), size=2) +
  geom_segment(aes(x=0.375, y=0.375, xend=.75, yend=0), size=2) +
  geom_point(aes(x=0, y=0), size=10, shape=21, color="black", fill="#23FF06") +
  geom_point(aes(x=0.75, y=0), size=10, shape=24, color="black", fill="#0000FF")+
  geom_point(aes(x=1.5, y=0), size=10, shape=24, color="black", fill="#0000FF")+
  #annotate(geom = "text", x = 0, y=-0.08, label="low distant") + 
  #annotate(geom = "text", x = 0.75, y=-0.08, label="highlands") + 
  #annotate(geom = "text", x = 1.5, y=-0.08, label="highlands\nother species") + 
  #annotate(geom = "text", x = 2, y=-0.08, label="outgroup") + 
  ylim(-0.1, 1.1)+ xlim(-0.1, 2.1) + theme_transparent() 

ggplot()+
  geom_segment(aes(x=0, y=0, xend=1.05, yend=1.05), size=2)+
  geom_segment(aes(x=1, y=1, xend=2, yend=0), size=2)+
  geom_segment(aes(x=0.75, y=0.75, xend=1.5, yend=0), size=2) +
  geom_segment(aes(x=0.375, y=0.375, xend=.75, yend=0), size=2) +
  geom_point(aes(x=0, y=0), size=10, shape=21, color="black", fill="#23FF06") +
  geom_point(aes(x=0.75, y=0), size=10, shape=24, color="black", fill="#0000FF")+
  geom_point(aes(x=1.5, y=0), size=10, shape=24, color="black", fill="white")+
  #annotate(geom = "text", x = 0, y=-0.08, label="low distant") + 
  #annotate(geom = "text", x = 0.75, y=-0.08, label="highlands") + 
  #annotate(geom = "text", x = 1.5, y=-0.08, label="highlands\nother species") + 
  #annotate(geom = "text", x = 2, y=-0.08, label="outgroup") + 
  ylim(-0.1, 1.1)+ xlim(-0.1, 2.1) + theme_transparent() 


################ tree era ################ 
tree.era <- read.tree(file = "local/data/joana.tree/new.tree/RAxML_bipartitions.eratoRelatives.withSRA.chr1-21.max0.5N.minDP3.balSubset.thin100bp.singleSamples.renamed")
#patristic: patristic distance, i.e. sum of branch lengths 
?distTips
dists.era <- as.matrix(distTips(tree.era, tips = "all", method = c("patristic"), useC = TRUE)); dists.era
rownames(dists.era) <- tree.era$tip.label; colnames(dists.era) <- tree.era$tip.label; dists.era

################ tree mel ################ 
tree.mel <- read.tree(file = "local/data/joana.tree/new.tree/RAxML_bipartitions.melpomeneRelatives.withSRA.chr1-21.max0.5N.balTree.100bpThinned.singleSamples.renamed")
#patristic: patristic distance, i.e. sum of branch lengths 
#?distTips
dists.mel <- as.matrix(distTips(tree.mel, tips = "all", method = c("patristic"), useC = TRUE)); dists.mel
rownames(dists.mel) <- tree.mel$tip.label; colnames(dists.mel) <- tree.mel$tip.label; dists.mel


######################################################## 2. era fd/fdm from other species (parapatric,  n=8 combos (2 tele, 2 him, 4 clys)) ############################


############################ 2.4 era plot fdm across pops and species ############################
##### read in shdr/ fdm outlier data #####
shdr.era.east.outlier.df <- read.csv("local/data/shdr.summ/shdr.era.east.outlier.df.csv")
shdr.era.west.outlier.df <- read.csv("local/data/shdr.summ/shdr.era.west.outlier.df.csv")

# first get rid off shitty hdrs
names(shdr.era.east.outlier.df)
shdr.era.east.outlier.df <- subset(shdr.era.east.outlier.df, is.max.pbs.hig.above4=="yes")
shdr.era.west.outlier.df <- subset(shdr.era.west.outlier.df, is.max.pbs.hig.above4=="yes")


#### change in % fdm outliers with genetic distance ####
era.fd.list.files
names(shdr.era.east.outlier.df); shdr.era.east.outlier.df$is.allo.hdr

plot(tree.era)
tree.era$tip.label

era.shdr.SAV.other.sp.summ <- data_frame(donor.p3=c(rep("era.west.high", 2), rep("era.east.high", 2), rep("him.east.high", 2),  rep("tel.east.high", 2),  rep("cly.east.high", 2),  rep("cly.west.high", 2)),
                                         recipient.p2=c(rep("era.west.high", 2), rep("era.east.high", 2), rep("him.east.high", 2),  rep("tel.east.high", 2),  rep("cly.east.high", 2),  rep("cly.west.high", 2)),
                                country.p2= c(rep(c("col", "ecu" ), 6)),
                                side.p2= c(rep("P2 = era. high East", 2), rep("P2 = era. high West", 2), rep("P2 = era. high East", 6), rep("P2 = era. high West", 2)),
                                pop.p3= c(rep("1.era", 4), rep("2.him", 2), rep("3.tel", 2), rep("4.cly", 4)),
                                percentage.shdr.with.excess.allele.sharing= c(
                                  nrow(subset(shdr.era.east.outlier.df, fdm.max.more.90th.abs.fdm.min.per.SHDR_fd_eraCvlowE_eraChighE_eraChighW_eleHW=="yes" ))/ nrow(subset(shdr.era.east.outlier.df)) *100,
                                  nrow(subset(shdr.era.east.outlier.df, fdm.max.more.90th.abs.fdm.min.per.SHDR_fd_eraEvlowE_eraEhighE_eraEhighW_eleHW=="yes"))/ nrow(subset(shdr.era.east.outlier.df)) *100,
                                  nrow(subset(shdr.era.west.outlier.df, fdm.max.more.90th.abs.fdm.min.per.SHDR_fd_eraCvlowW_eraChighW_eraChighE_eleHW=="yes"))/ nrow(subset(shdr.era.west.outlier.df)) *100,
                                  nrow(subset(shdr.era.west.outlier.df, fdm.max.more.90th.abs.fdm.min.per.SHDR_fd_eraEvlowW_eraEhighW_eraEhighE_eleHW=="yes"))/ nrow(subset(shdr.era.west.outlier.df)) *100,
                                  nrow(subset(shdr.era.east.outlier.df, fdm.max.more.90th.abs.fdm.min.per.SHDR_fd_eraCvlowE_eraChighE_himHE_eleHW=="yes"))/ nrow(subset(shdr.era.east.outlier.df)) *100,
                                  nrow(subset(shdr.era.east.outlier.df, fdm.max.more.90th.abs.fdm.min.per.SHDR_fd_eraEvlowE_eraEhighE_himHE_eleHW=="yes"))/ nrow(subset(shdr.era.east.outlier.df)) *100,
                                  nrow(subset(shdr.era.east.outlier.df, fdm.max.more.90th.abs.fdm.min.per.SHDR_fd_eraCvlowE_eraChighE_telhighE_eleHW=="yes"))/ nrow(subset(shdr.era.east.outlier.df)) *100,
                                  nrow(subset(shdr.era.east.outlier.df, fdm.max.more.90th.abs.fdm.min.per.SHDR_fd_eraEvlowE_eraEhighE_telhighE_eleHW=="yes"))/ nrow(subset(shdr.era.east.outlier.df)) *100,
                                  nrow(subset(shdr.era.east.outlier.df, fdm.max.more.90th.abs.fdm.min.per.SHDR_fd_eraCvlowE_eraChighE_clyEhighE_eleHW=="yes"))/ nrow(subset(shdr.era.east.outlier.df)) *100,
                                  nrow(subset(shdr.era.east.outlier.df, fdm.max.more.90th.abs.fdm.min.per.SHDR_fd_eraEvlowE_eraEhighE_clyEhighE_eleHW=="yes"))/ nrow(subset(shdr.era.east.outlier.df)) *100,
                                  nrow(subset(shdr.era.west.outlier.df, fdm.max.more.90th.abs.fdm.min.per.SHDR_fd_eraCvlowW_eraChighW_clyEhighW_eleHW=="yes"))/ nrow(subset(shdr.era.west.outlier.df)) *100,
                                  nrow(subset(shdr.era.west.outlier.df, fdm.max.more.90th.abs.fdm.min.per.SHDR_fd_eraEvlowW_eraEhighW_clyEhighW_eleHW=="yes"))/ nrow(subset(shdr.era.west.outlier.df)) *100),
                                dist.p2.p3= c(
                                  dists.era[, "CAM040106_eratoChighE"][c("15N120_eratoChighW")],
                                  dists.era[, "CAM011441_eratoEhighE"][c("CAM040545_eratoEhighW")],
                                  dists.era[, "15N120_eratoChighW"][c("CAM040106_eratoChighE")],
                                  dists.era[, "CAM040545_eratoEhighW"][c("CAM011441_eratoEhighE")],
                                  dists.era[, "CAM040106_eratoChighE"][c("SRR4032038_himeraHE")],
                                  dists.era[, "CAM011441_eratoEhighE"][c("SRR4032038_himeraHE")],
                                  dists.era[, "CAM040106_eratoChighE"][c("CAM040276_telesiphEhighE")],
                                  dists.era[, "CAM011441_eratoEhighE"][c("CAM040276_telesiphEhighE")],
                                  dists.era[, "CAM040106_eratoChighE"][c("CAM040271_clysonymusHE")],
                                  dists.era[, "CAM011441_eratoEhighE"][c("CAM040271_clysonymusHE")],
                                  dists.era[, "15N120_eratoChighW"][c("CAM002846_clysonymusHW")],
                                  dists.era[, "CAM040545_eratoEhighW"][c("CAM002846_clysonymusHW")]
                                ),
                                donor.p3.country= paste(donor.p3, country.p2, sep="_"),
                                side.country.p2= paste(side.p2, country.p2, sep="_"),
                                type.p3=c(rep("within.species", 4), rep("closely.related", 2), rep("distantly.related",6)),
                                type.p3.side.country.p2=paste(type.p3, side.country.p2, sep="_" ),
); era.shdr.SAV.other.sp.summ

dists.era[, "CAM040106_eratoChighE"][c("MAC004050_eratoCvlowE")]



#### plot ####


ggplot(data=subset(era.shdr.SAV.other.sp.summ, donor.p3!="era.west.high" & donor.p3!="era.east.high" ), aes(x=dist.p2.p3, y=percentage.shdr.with.excess.allele.sharing,  
                                                                                                            group=side.country.p2, fill=type.p3, color=type.p3))+
  geom_line(inherit.aes = F, color= "black", data=subset(era.shdr.SAV.other.sp.summ, donor.p3!="era.west.high" & donor.p3!="era.east.high" ), aes(x=dist.p2.p3, y=percentage.shdr.with.excess.allele.sharing, group=side.country.p2, linetype=country.p2), size=1)+
  geom_point(shape=24, size=5, stroke=1)  +
  scale_fill_manual(values=c("#eb8055ff","#efe350ff"  ))+
  scale_color_manual(values=c("black",  "black"))+
  
  #geom_text(aes(label=donor.p3),hjust=-0.2, vjust=0) +
  geom_point(inherit.aes = F, data=subset(era.shdr.SAV.other.sp.summ, (donor.p3=="era.west.high" | donor.p3=="era.east.high")  ), aes(x=dist.p2.p3, y=percentage.shdr.with.excess.allele.sharing,  group=side.country.p2),
             shape=24, fill="blue", size=5, stroke=1)+
  #geom_text(inherit.aes = F, data=subset(era.shdr.SAV.other.sp.summ, (donor.p3=="era.west.high" | donor.p3=="era.east.high" )&country.p2=="col" ), aes(x=dist.p2.p3, y=percentage.shdr.with.excess.allele.sharing,  group=side.country.p2,  label=donor.p3),hjust=0, vjust=0) +
  scale_linetype_manual(values=c("solid", "11"))+
  facet_wrap(~side.p2, ncol = 1)+
  xlab("Phylogenetic distance: P2 - P3")+ ylab("% SHDR with outlier maximum FdM")+
  scale_y_continuous(expand = c(0.1,0.1))+scale_x_continuous(expand = c(0.002,0.002))+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=.5),legend.position="none", axis.ticks.length=unit(-0.15, "cm"),
        plot.margin = unit(c(0.1,0.1,0,0.1), "lines"),axis.text = element_text(size=12), axis.text.y = element_text(size=10), 
        axis.title.x = element_text(size=16),axis.title.y = element_text(size=16),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),  strip.text = element_text(size=14, face="bold"), strip.background = element_blank(),
        panel.background = element_rect(fill = "transparent", color = NA), plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
        legend.background = element_rect(fill = "transparent", color = NA), legend.box.background = element_rect(fill = "transparent", color = NA), legend.key = element_rect(fill = "transparent", color = NA),
        strip.background.y = element_blank(), rect = element_rect(fill = "transparent"))
ggsave("plots/joana.fdm/fdm.perc.shdrs.outliers.era.other.sp.png", height = 7, width = 3.7, bg="transparent")

ggplot(data=era.shdr.SAV.other.sp.summ, aes(x=dist.p2.p3, y=percentage.shdr.with.excess.allele.sharing,  group=side.country.p2))+
  geom_point() +
  # geom_point(inherit.aes = F, data=subset(era.shdr.SAV.other.sp.summ, donor.p3=="era.west.high" | donor.p3=="era.east.high" ), aes(x=order(dist.p2.p3), y=percentage.shdr.with.excess.allele.sharing,  group=side.country.p2),
  #            shape=24, fill="blue", size=5)+
  geom_text(aes(label=donor.p3),hjust=0, vjust=0) +
  scale_linetype_manual(values=c("solid", "11"))+
  geom_line(aes(linetype=country.p2), size=1)+
  facet_wrap(~side.p2)+
  xlab("Phylogenetic distance between P2 and P3")+ ylab("% SHDR with outlier maximum FdM")+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=.5),legend.position="bottom", axis.ticks.length=unit(-0.15, "cm"),
        plot.margin = unit(c(0.1,0.1,0,0.1), "lines"),axis.text = element_text(size=12), axis.text.y = element_text(size=10), 
        axis.title.x = element_text(size=12),axis.title.y = element_text(size=12),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),  strip.text = element_text(size=12), strip.background = element_blank(),
        panel.background = element_rect(fill = "transparent", color = NA), plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
        legend.background = element_rect(fill = "transparent", color = NA), legend.box.background = element_rect(fill = "transparent", color = NA), legend.key = element_rect(fill = "transparent", color = NA),
        strip.background.y = element_blank(), rect = element_rect(fill = "transparent"))

ggsave("plots/joana.fdm/fdm.boxplot.era.ancestral.variation.png", height = 4, width = 6, bg="transparent")

#### ks tests #####
# null hypothesis that the true distribution function of x is  not greater than the distribution function of y
# alternative hypothesis is that dist of x is greater than dist of y
ks.test(shdr.era.west.outlier.df$fdm.max_fd_eraCvlowW_eraChighW_eraChighE_eleHW, abs(shdr.era.west.outlier.df$fdm.min_fd_eraCvlowW_eraChighW_eraChighE_eleHW), alternative = "less" ) # D^+ = 0.27586, p-value = 0.01211
ks.test(shdr.era.west.outlier.df$fdm.max_fd_eraEvlowW_eraEhighW_eraEhighE_eleHW, abs(shdr.era.west.outlier.df$fdm.min_fd_eraEvlowW_eraEhighW_eraEhighE_eleHW), alternative = "less" ) # D^+ = 0.086207, p-value = 0.6498

ks.test(shdr.era.east.outlier.df$fdm.max_fd_eraCvlowE_eraChighE_eraChighW_eleHW, abs(shdr.era.east.outlier.df$fdm.min_fd_eraCvlowE_eraChighE_eraChighW_eleHW), alternative = "less" )
ks.test(shdr.era.east.outlier.df$fdm.max_fd_eraEvlowE_eraEhighE_eraEhighW_eleHW, abs(shdr.era.east.outlier.df$fdm.min_fd_eraEvlowE_eraEhighE_eraEhighW_eleHW), alternative = "less" )

# across sp
ks.test(shdr.era.east.outlier.df$fdm.max_fd_eraCvlowE_eraChighE_himHE_eleHW, abs(shdr.era.east.outlier.df$fdm.min_fd_eraCvlowE_eraChighE_himHE_eleHW), alternative = "less" ) 
ks.test(shdr.era.east.outlier.df$fdm.max_fd_eraEvlowE_eraEhighE_himHE_eleHW, abs(shdr.era.east.outlier.df$fdm.min_fd_eraEvlowE_eraEhighE_himHE_eleHW), alternative = "less" )

ks.test(shdr.era.east.outlier.df$fdm.max_fd_eraCvlowE_eraChighE_telhighE_eleHW, abs(shdr.era.east.outlier.df$fdm.min_fd_eraCvlowE_eraChighE_telhighE_eleHW), alternative = "less" ) 
ks.test(shdr.era.east.outlier.df$fdm.max_fd_eraEvlowE_eraEhighE_telhighE_eleHW, abs(shdr.era.east.outlier.df$fdm.min_fd_eraEvlowE_eraEhighE_telhighE_eleHW), alternative = "less" )

ks.test(shdr.era.east.outlier.df$fdm.max_fd_eraCvlowE_eraChighE_clyEhighE_eleHW, abs(shdr.era.east.outlier.df$fdm.min_fd_eraCvlowE_eraChighE_clyEhighE_eleHW), alternative = "less" ) 
ks.test(shdr.era.east.outlier.df$fdm.max_fd_eraEvlowE_eraEhighE_clyEhighE_eleHW, abs(shdr.era.east.outlier.df$fdm.min_fd_eraEvlowE_eraEhighE_clyEhighE_eleHW), alternative = "less" )

ks.test(shdr.era.west.outlier.df$fdm.max_fd_eraCvlowW_eraChighW_clyEhighW_eleHW, abs(shdr.era.west.outlier.df$fdm.min_fd_eraCvlowW_eraChighW_clyEhighW_eleHW), alternative = "less" ) # D^+ = 0.27586, p-value = 0.01211
ks.test(shdr.era.west.outlier.df$fdm.max_fd_eraEvlowW_eraEhighW_clyEhighW_eleHW, abs(shdr.era.west.outlier.df$fdm.min_fd_eraEvlowW_eraEhighW_clyEhighW_eleHW), alternative = "less" ) # D^+ = 0.086207, p-value = 0.6498




############################ 2.4 mel plot fdm across pops and species ############################
##### read in shdr/ fdm outlier data #####
shdr.mel.east.outlier.df <- read.csv("local/data/shdr.summ/shdr.mel.east.outlier.df.csv")
shdr.mel.west.outlier.df <- read.csv("local/data/shdr.summ/shdr.mel.west.outlier.df.csv")

# first get rid off shitty hdrs
names(shdr.mel.east.outlier.df)
shdr.mel.east.outlier.df <- subset(shdr.mel.east.outlier.df, is.max.pbs.hig.above4=="yes")
shdr.mel.west.outlier.df <- subset(shdr.mel.west.outlier.df, is.max.pbs.hig.above4=="yes")


#### change in % fdm outliers with genetic distance ####
mel.fd.list.files
names(shdr.mel.east.outlier.df); shdr.mel.east.outlier.df$is.allo.hdr
shdr.mel.west.outlier.df$fdm.max.more.90th.abs.fdm.min.per.SHDR_fd_melCvlowW_melChighW_clyEhighW_eleHW
tree.mel$tip.label

mel.shdr.SAV.other.sp.summ <- data_frame(donor.p3=c( rep("mel.west.high", 2) ,rep("mel.east.high", 2), rep("tim.east.high", 2),  rep("cyd.west.high", 2)),
                                         country.p2= c(rep(c("col", "ecu" ), 4)),
                                         side.p2= c(rep("P2 = mel. high East", 2), rep("P2 = mel. high West", 2), rep("P2 = mel. high East", 2), rep("P2 = mel. high West", 2)),
                                         dist= c(rep("2.mel", 4), rep("1.tim.E", 2), rep("3.cyd.W", 2)),
                                         percentage.shdr.with.excess.allele.sharing= c(
                                           nrow(subset(shdr.mel.east.outlier.df, fdm.max.more.90th.abs.fdm.min.per.SHDR_fd_melCvlowE_melChighE_melChighW_hecale=="yes" ))/ nrow(subset(shdr.mel.east.outlier.df)) *100,
                                           nrow(subset(shdr.mel.east.outlier.df, fdm.max.more.90th.abs.fdm.min.per.SHDR_fd_melEvlowE_melEhighE_melEhighW_hecale=="yes"))/ nrow(subset(shdr.mel.east.outlier.df)) *100,
                                           nrow(subset(shdr.mel.west.outlier.df, fdm.max.more.90th.abs.fdm.min.per.SHDR_fd_melCvlowW_melChighW_melChighE_hecale=="yes"))/ nrow(subset(shdr.mel.west.outlier.df)) *100,
                                           nrow(subset(shdr.mel.west.outlier.df, fdm.max.more.90th.abs.fdm.min.per.SHDR_fd_melElowW_melEhighW_melEhighE_hecale=="yes"))/ nrow(subset(shdr.mel.west.outlier.df)) *100,
                                           nrow(subset(shdr.mel.east.outlier.df, fdm.max.more.90th.abs.fdm.min.per.SHDR_fd_melCvlowE_melChighE_timEhighE_hecale=="yes"))/ nrow(subset(shdr.mel.east.outlier.df)) *100,
                                           nrow(subset(shdr.mel.east.outlier.df, fdm.max.more.90th.abs.fdm.min.per.SHDR_fd_melEvlowE_melEhighE_timEhighE_hecale=="yes"))/ nrow(subset(shdr.mel.east.outlier.df)) *100,
                                           nrow(subset(shdr.mel.west.outlier.df, fdm.max.more.90th.abs.fdm.min.per.SHDR_fd_melCvlowW_melChighW_cydEhighW_hecale=="yes"))/ nrow(subset(shdr.mel.west.outlier.df)) *100,
                                           nrow(subset(shdr.mel.west.outlier.df, fdm.max.more.90th.abs.fdm.min.per.SHDR_fd_melElowW_melEhighW_cydEhighW_hecale=="yes"))/ nrow(subset(shdr.mel.west.outlier.df)) *100),
                                         dist.p2.p3= c(
                                             dists.mel[, "CAM040043_melChighE"][c("CS003603_melChighW")],
                                             dists.mel[, "CAM041101_melEhighE"][c("14N004_melEhighW")],
                                             dists.mel[, "CS003603_melChighW"][c("CAM040043_melChighE")],
                                             dists.mel[, "14N004_melEhighW"][c("CAM041101_melEhighE")],
                                             dists.mel[, "CAM040043_melChighE"][c("CAM008521_timEhighE")],
                                             dists.mel[, "CAM041101_melEhighE"][c("CAM008521_timEhighE")],
                                             dists.mel[, "CS003603_melChighW"][c("CAM040419_cydEhighW")],
                                             dists.mel[, "14N004_melEhighW"][c("CAM040419_cydEhighW")]
                                           ),
                                         donor.p3.country= paste(donor.p3, country.p2, sep="_"),
                                         side.country.p2= paste(side.p2, country.p2, sep="_"),
                                         type.p3=c(rep("within.species", 4), rep("closely.related", 4)),
                                         type.p3.side.country.p2=paste(type.p3, side.country.p2, sep="_" )
); mel.shdr.SAV.other.sp.summ

#### ks tests #####
# null hypothesis that the true distribution function of x is  not greater than the distribution function of y
# alternative hypothesis is that dist of x is greater than dist of y
ks.test(shdr.mel.west.outlier.df$fdm.max_fd_melCvlowW_melChighW_melChighE_hecale, abs(shdr.mel.west.outlier.df$fdm.min_fd_melCvlowW_melChighW_melChighE_hecale), alternative = "less" ) # D^+ = 0.27586, p-value = 0.01211
ks.test(shdr.mel.west.outlier.df$fdm.max_fd_melElowW_melEhighW_melEhighE_hecale, abs(shdr.mel.west.outlier.df$fdm.min_fd_melElowW_melEhighW_melEhighE_hecale), alternative = "less" ) # D^+ = 0.086207, p-value = 0.6498


ks.test(shdr.mel.east.outlier.df$fdm.max_fd_melCvlowE_melChighE_melChighW_hecale, abs(shdr.mel.east.outlier.df$fdm.min_fd_melCvlowE_melChighE_melChighW_hecale), alternative = "less" )
ks.test(shdr.mel.east.outlier.df$fdm.max_fd_melEvlowE_melEhighE_melEhighW_hecale, abs(shdr.mel.east.outlier.df$fdm.min_fd_melEvlowE_melEhighE_melEhighW_hecale), alternative = "less" )

# across sp
ks.test(shdr.mel.east.outlier.df$fdm.max_fd_melCvlowW_melChighW_cydEhighW_hecale, abs(shdr.mel.east.outlier.df$fdm.min_fd_melCvlowW_melChighW_cydEhighW_hecale), alternative = "less" ) 
ks.test(shdr.mel.east.outlier.df$fdm.max_fd_melElowW_melEhighW_cydEhighW_hecale, abs(shdr.mel.east.outlier.df$fdm.min_fd_melElowW_melEhighW_cydEhighW_hecale), alternative = "less" )


ks.test(shdr.mel.east.outlier.df$fdm.max_fd_melCvlowE_melChighE_timEhighE_hecale, abs(shdr.mel.east.outlier.df$fdm.min_fd_melCvlowE_melChighE_timEhighE_hecale), alternative = "less" ) 
ks.test(shdr.mel.east.outlier.df$fdm.max_fd_melEvlowE_melEhighE_timEhighE_hecale, abs(shdr.mel.east.outlier.df$fdm.min_fd_melEvlowE_melEhighE_timEhighE_hecale), alternative = "less" )



#### plot ####

ggplot(data=subset(mel.shdr.SAV.other.sp.summ, donor.p3!="mel.west.high" & donor.p3!="mel.east.high" ), aes(x=dist.p2.p3, y=percentage.shdr.with.excess.allele.sharing,  
                                                                                                            group=side.country.p2, fill=type.p3, color=type.p3))+
  geom_line(inherit.aes = F, color= "black", data=subset(mel.shdr.SAV.other.sp.summ, donor.p3!="mel.west.high" & donor.p3!="mel.east.high" ), aes(x=dist.p2.p3, y=percentage.shdr.with.excess.allele.sharing, group=side.country.p2, linetype=country.p2), size=1)+
  geom_point(shape=24, size=5, stroke=1)  +
  scale_fill_manual(values=c("#eb8055ff","#efe350ff"  ))+
  scale_color_manual(values=c("black",  "black"))+
  
  #geom_text(aes(label=donor.p3),hjust=-0.2, vjust=0) +
  geom_point(inherit.aes = F, data=subset(mel.shdr.SAV.other.sp.summ, (donor.p3=="mel.west.high" | donor.p3=="mel.east.high")  ), aes(x=dist.p2.p3, y=percentage.shdr.with.excess.allele.sharing,  group=side.country.p2),
             shape=24, fill="blue", size=5, stroke=1)+
  #geom_text(inherit.aes = F, data=subset(mel.shdr.SAV.other.sp.summ, (donor.p3=="mel.west.high" | donor.p3=="mel.east.high" )&country.p2=="col" ), aes(x=dist.p2.p3, y=percentage.shdr.with.excess.allele.sharing,  group=side.country.p2,  label=donor.p3),hjust=0, vjust=0) +
  scale_linetype_manual(values=c("solid", "11"))+
  facet_wrap(~side.p2, ncol = 1)+
  xlab("Phylogenetic distance: P2 - P3")+ ylab("% SHDR with outlier maximum FdM")+
  scale_y_continuous(expand = c(0.1,0.1))+scale_x_continuous(expand = c(0.002,0.002))+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=.5),legend.position="none", axis.ticks.length=unit(-0.15, "cm"),
        plot.margin = unit(c(0.1,0.1,0,0.1), "lines"),axis.text = element_text(size=12), axis.text.y = element_text(size=10), 
        axis.title.x = element_text(size=16),axis.title.y = element_text(size=16),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),  strip.text = element_text(size=14, face="bold"), strip.background = element_blank(),
        panel.background = element_rect(fill = "transparent", color = NA), plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
        legend.background = element_rect(fill = "transparent", color = NA), legend.box.background = element_rect(fill = "transparent", color = NA), legend.key = element_rect(fill = "transparent", color = NA),
        strip.background.y = element_blank(), rect = element_rect(fill = "transparent"))
ggsave("plots/joana.fdm/fdm.perc.shdrs.outliers.mel.other.sp.png", height = 7, width = 3.7, bg="transparent")

ggplot(data=subset(mel.shdr.SAV.other.sp.summ, donor.p3!="mel.west.high" & donor.p3!="mel.east.high" ), aes(x=dist.p2.p3, y=percentage.shdr.with.excess.allele.sharing,  
                                                                                                            group=side.country.p2, fill=type.p3, color=type.p3))+
  geom_line(inherit.aes = F, color= "black", data=subset(mel.shdr.SAV.other.sp.summ, donor.p3!="mel.west.high" & donor.p3!="mel.east.high" ), aes(x=dist.p2.p3, y=percentage.shdr.with.excess.allele.sharing, group=side.country.p2, linetype=country.p2), size=1)+
  geom_point(shape=24, size=5, stroke=1)  +
  scale_fill_manual(values=c("#eb8055ff","#efe350ff"  ))+
  scale_color_manual(values=c("black",  "black"))+
  
  geom_text(aes(label=donor.p3.country),hjust=-0.2, vjust=0) +
  geom_point(inherit.aes = F, data=subset(mel.shdr.SAV.other.sp.summ, (donor.p3=="mel.west.high" | donor.p3=="mel.east.high")  ), aes(x=dist.p2.p3, y=percentage.shdr.with.excess.allele.sharing,  group=side.country.p2),
             shape=24, fill="blue", size=5, stroke=1)+
  geom_text(inherit.aes = F, data=subset(mel.shdr.SAV.other.sp.summ, (donor.p3=="mel.west.high" | donor.p3=="mel.east.high" ) ), aes(x=dist.p2.p3, y=percentage.shdr.with.excess.allele.sharing,  group=side.country.p2,  label=donor.p3.country),hjust=0, vjust=0) +
  scale_linetype_manual(values=c("solid", "11"))+
  facet_wrap(~side.p2, ncol = 1)+
  xlab("Phylogenetic distance: P2 - P3")+ ylab("% SHDR with outlier maximum FdM")+
  scale_y_continuous(expand = c(0.1,0.1))+scale_x_continuous(expand = c(0.002,0.002))+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=.5),legend.position="none", axis.ticks.length=unit(-0.15, "cm"),
        plot.margin = unit(c(0.1,0.1,0,0.1), "lines"),axis.text = element_text(size=12), axis.text.y = element_text(size=10), 
        axis.title.x = element_text(size=16),axis.title.y = element_text(size=16),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),  strip.text = element_text(size=14, face="bold"), strip.background = element_blank(),
        panel.background = element_rect(fill = "transparent", color = NA), plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
        legend.background = element_rect(fill = "transparent", color = NA), legend.box.background = element_rect(fill = "transparent", color = NA), legend.key = element_rect(fill = "transparent", color = NA),
        strip.background.y = element_blank(), rect = element_rect(fill = "transparent"))

