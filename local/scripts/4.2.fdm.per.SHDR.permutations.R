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
library(phangorn)
setwd("/Users/gabrielamontejokovacevich/Dropbox (Cambridge University)/PhD/7_Assocation_studies/9_ANGSD/12.era.mel.altitude/local/")


############################ functions ############################

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

getTWISST.mel<-function(combination){
  
  # Read in weights twisstset
  twisst<-read.table(paste0(prefix,combination,".weights.csv"),header=T)
  
  # add window information
  twisst<-cbind(info,twisst)
  twisst<-data.table(twisst)
  twisst$start.WG <-twisst$start+ref.scaff.mel[match(twisst$scaffold,ref.scaff.mel$scaffold),"add"]
  twisst$end.WG <-twisst$end+ref.scaff.mel[match(twisst$scaffold,ref.scaff.mel$scaffold),"add"]
  
  # topo weightings (0-1) per window
  twisst$topo1.freq <- round(twisst$topo1/rowSums(twisst[,c("topo1","topo2","topo3")]),digits=4)
  twisst$topo2.freq <- round(twisst$topo2/rowSums(twisst[,c("topo1","topo2","topo3")]),digits=4)
  twisst$topo3.freq <- round(twisst$topo3/rowSums(twisst[,c("topo1","topo2","topo3")]),digits=4)
  # freq topo2 -topo3
  twisst$topo2.minus.topo3.freq <- twisst$topo2.freq - twisst$topo3.freq
  
  
  # add if east (only) PBS outlier (para SHDR)
  shdr.east.ranges<-data.table(mel.outliersEast[,c("scaff","start.bp","end.bp")])
  names(shdr.east.ranges)<-c("scaffold","start","end") # to match twisst names
  setkey(shdr.east.ranges)
  # by.x concatenates scaffold and start/end
  oveast<-foverlaps(x = twisst,y = shdr.east.ranges, which=T, by.x = names(shdr.east.ranges), type="any",nomatch=0L)
  twisst$shdr.para.east.id<- NA
  # use keys to select corresponding shdr names
  twisst$shdr.para.east.id[oveast$xid]<- paste(mel.outliersEast[ oveast$yid,]$shdr.para.east.id)
  
  # add if west (only) PBS outlier (para SHDR)
  shdr.west.ranges<-data.table(mel.outliersWest[,c("scaff","start.bp","end.bp")])
  names(shdr.west.ranges)<-c("scaffold","start","end") # to match twisst names
  setkey(shdr.west.ranges)
  # by.x concatenates scaffold and start/end
  ovWest<-foverlaps(x = twisst,y = shdr.west.ranges, which=T, by.x = names(shdr.west.ranges), type="any",nomatch=0L)
  twisst$shdr.para.west.id<- NA
  # use keys to select corresponding shdr names
  twisst$shdr.para.west.id[ovWest$xid]<-  paste(mel.outliersWest[ ovWest$yid,]$shdr.para.west.id)
  
  # add if across (only) PBS outlier (para SHDR)
  shdr.across.ranges<-data.table(mel.outliersAcross[,c("scaff","start.bp","end.bp")])
  names(shdr.across.ranges)<-c("scaffold","start","end") # to match twisst names
  setkey(shdr.across.ranges)
  # by.x concatenates scaffold and start/end
  ovacross<-foverlaps(x = twisst,y = shdr.across.ranges, which=T, by.x = names(shdr.across.ranges), type="any",nomatch=0L)
  twisst$shdr.allo.id<- NA
  # use keys to select corresponding shdr names
  twisst$shdr.allo.id[ovacross$xid]<-  paste(mel.outliersAcross[ ovacross$yid,]$shdr.allo.id)
  
  return(twisst)
}
#######################################################################################################   ERATO fd ################################################################

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

# new version Oct21- use combined para/allo shdrs
era.outliersEast<-era.shdr.para.east
era.outliersWest<-era.shdr.para.west



# Read in chromosome information
ref.scaff.era <- read.table("local/data/ref/Heliconius_erato_demophoon_v1_-_scaffolds.fa.fai", row.names = NULL)
ref.scaff.era<-ref.scaff.era[,1:2]
names(ref.scaff.era)<-c("scaffold","length")
ref.scaff.era$add<-c(0,cumsum(ref.scaff.era$length)[-length(ref.scaff.era$length)])
ref.scaff.era$chr<-substring(ref.scaff.era$scaffold,first = 7,last=8);head(ref.scaff.era)

scafEnds.era <- cumsum(ref.scaff.era[,2])
era.genome_size <- tail(scafEnds.era,1)

######################################################## 1. trees ############################
################################ main text single samples ################################

#### read in trees, then save tip labels (indivs used) ####
tree.era <- read.tree(file = "local/data/joana.tree/new.tree/RAxML_bipartitions.eratoRelatives.withSRA.chr1-21.max0.5N.minDP3.balSubset.thin100bp.singleSamples.renamed")
plotBS(tree.era, type = "phylogram")
tip.labels.tree.era <- tree.era$tip.label; tip.labels.tree.era
tip.labels.tree.era <- cbind(as.data.frame(str_split_fixed(tip.labels.tree.era, "_", 2)),tip.labels.tree.era ); names(tip.labels.tree.era) <- c("id", "type.joana", "tip.labels.tree.era"); head(tip.labels.tree.era)
tip.labels.tree.era
# rename just with id
tree.era$tip.label <- as.character(tip.labels.tree.era$id); plot(tree.era)


#write.csv(tip.labels.tree.era, "local/data/joana.tree/era.samples.info.tree.csv", row.names = F)
plotBS(tree.era, type = "phylogram")

tree.mel <- read.tree(file = "local/data/joana.tree/new.tree/RAxML_bipartitions.melpomeneRelatives.withSRA.chr1-21.max0.5N.balTree.100bpThinned.singleSamples.renamed")
plotBS(tree.mel, type = "phylogram")
tip.labels.tree.mel <- tree.mel$tip.label; tip.labels.tree.mel
tip.labels.tree.mel <- cbind(as.data.frame(str_split_fixed(tip.labels.tree.mel, "_", 2)),tip.labels.tree.mel ); names(tip.labels.tree.mel) <- c("id", "type.joana", "tip.labels.tree.mel"); head(tip.labels.tree.mel)
tip.labels.tree.mel
# rename just with id
tree.mel$tip.label <- as.character(tip.labels.tree.mel$id); plot(tree.mel)

#### read in info from all samples used by joana #### 
era.samples.info <- read.csv(file = "local/data/joana.fd/erato.samples.info.csv"); head(era.samples.info)
mel.samples.info <- read.csv(file = "local/data/joana.fd/melpomene.samples.info.csv"); head(mel.samples.info)

# subset info in trees
era.samples.tree.info <- subset(era.samples.info, id.joana %in% tip.labels.tree.era$id )
mel.samples.tree.info <- subset(mel.samples.info, id.joana %in% tip.labels.tree.mel$id )

#### read in fbranch ####

#### plot tree with arrows ####
library(ggrepel)
library(dplyr)
library(ggtree)

era.full.tree <- ggtree(tree.era, ladderize =F ) %<+% era.samples.tree.info + 
  geom_tiplab(aes(label=tip.label), align = T, hjust = -0.2) +
  #geom_nodelab(geom = "text",aes(label = node) )+
  geom_label2(aes(label=label, subset = !is.na(as.numeric(label)) & as.numeric(label) > 20), alpha=0.7, size=2)+
  geom_tippoint(aes(shape = alt.type, fill = type.tree ),   size=3) + 
  scale_shape_manual(values=c(NA,24,21))+
  scale_fill_manual(values = c('blue','green','pink','orange',  NA),label = c("Highland pop.","Lowland pop.","Distantly related\nhigh altitude sp.","Closely related\nhigh altitude sp.", "Outgroup") ) + 
  
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
ggsave("plots/joana.tree/era.main.tree.png", width = 5, height = 4, bg="transparent")

mel.full.tree <- ggtree(tree.mel, ladderize =F ) %<+% mel.samples.tree.info + 
  geom_tiplab(aes(label=sp.country.alt.side), align = T, hjust = -0.2) +
  #geom_nodelab(geom = "text",aes(label = node) )+
  geom_label2(aes(label=label, subset = !is.na(as.numeric(label)) & as.numeric(label) > 20), alpha=0.7, size=2)+
  geom_tippoint(aes(shape = alt.type, fill = type.tree ),   size=3) + 
  scale_shape_manual(values=c(NA,24,21))+
  scale_fill_manual(values = c('pink','orange',  'blue','green',NA) ) + 
  scale_fill_manual(values = c('pink','orange',  'blue','green',NA),label = c("Distantly related\nhigh altitude sp.","Closely related\nhigh altitude sp.","Highland pop.","Lowland pop.", "Outgroup") ) + 
  
  guides(fill = guide_legend(override.aes = list(shape = 21), title = "Type" ),
         shape = guide_legend(override.aes = list(fill = "black"), title = "Altitude"  ) ) +
  theme(legend.position = "right", legend.text = element_text(size = 14), legend.title = element_text(size = 14))  +
  xlim(0,0.07)+
  theme(axis.title.y = element_text(size=14,vjust=-2),
        plot.margin = unit(c(0,0.2,0,0.2), "lines"),axis.text.x = element_blank(), rect = element_rect(fill = "transparent") ,#axis.text.x = element_blank(),
        axis.text.y = element_blank(),axis.title.x = element_blank(),  axis.ticks.length=unit(-0.2, "cm"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_rect(fill = "transparent", color = NA),
        panel.background = element_rect(fill = "transparent", color = NA), plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
        legend.background = element_rect(fill = "transparent", color = NA), legend.box.background = element_rect(fill = "transparent", color = NA), legend.key = element_rect(fill = "transparent", color = NA),strip.background.y = element_blank(),
        strip.text = element_text(size=14, face = "bold"), strip.background = element_blank(), legend.position = "right")   ; mel.full.tree

ggsave("plots/joana.tree/mel.main.tree.png", width = 8, height = 6)




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
  geom_point(aes(x=0, y=0), size=5, shape=21, color="black", fill="#23FF06") +
  geom_point(aes(x=0.75, y=0), size=5, shape=24, color="black", fill="#0000FF")+
  geom_point(aes(x=1.5, y=0), size=5, shape=24, color="black", fill="orange")+
  annotate(geom = "text", x = 0, y=-0.08, label="low distant") + 
  annotate(geom = "text", x = 0.75, y=-0.08, label="highlands") + 
  annotate(geom = "text", x = 1.5, y=-0.08, label="highlands\nother species") + 
  annotate(geom = "text", x = 2, y=-0.08, label="outgroup") + 
  ylim(-0.1, 1.1)+ xlim(-0.1, 2.1) + theme_transparent() 


######################################################## 2. era fd/fdm from other species (parapatric,  n=8 combos (2 tele, 2 him, 4 clys)) ############################

############################ 2.1. era fd prep data ############################

#### Read in fd datasets ####
era.comparisons.info <- read.csv("local/data/joana.fd/era.comparisons.csv")

# within species allopatric from highlands
# include in simulations
fd_eraCvlowE_eraChighE_eraChighW_eleHW <-getFd.era("fd_eraCvlowE_eraChighE_eraChighW_eleHW")
fd_eraEvlowE_eraEhighE_eraEhighW_eleHW <-getFd.era("fd_eraEvlowE_eraEhighE_eraEhighW_eleHW")

fd_eraCvlowW_eraChighW_eraChighE_eleHW <-getFd.era("fd_eraCvlowW_eraChighW_eraChighE_eleHW")
fd_eraEvlowW_eraEhighW_eraEhighE_eleHW <-getFd.era("fd_eraEvlowW_eraEhighW_eraEhighE_eleHW")


# other.species.parapatric
# clysonymus east into highland col/ecu E
fd_eraCvlowE_eraChighE_clyEhighE_eleHW <-getFd.era("fd_eraCvlowE_eraChighE_clyEhighE_eleHW")
fd_eraEvlowE_eraEhighE_clyEhighE_eleHW <-getFd.era("fd_eraEvlowE_eraEhighE_clyEhighE_eleHW")

# clysonymus west into highland col/ecu W
fd_eraCvlowW_eraChighW_clyEhighW_eleHW<-getFd.era("fd_eraCvlowW_eraChighW_clyEhighW_eleHW")
fd_eraEvlowW_eraEhighW_clyEhighW_eleHW<-getFd.era("fd_eraEvlowW_eraEhighW_clyEhighW_eleHW")

# himera into highland col/ecu E
fd_eraCvlowE_eraChighE_himHE_eleHW<-getFd.era("fd_eraCvlowE_eraChighE_himHE_eleHW")
fd_eraEvlowE_eraEhighE_himHE_eleHW<-getFd.era("fd_eraEvlowE_eraEhighE_himHE_eleHW")

# telesiphe into highland col/ecu E
fd_eraCvlowE_eraChighE_telhighE_eleHW<-getFd.era("fd_eraCvlowE_eraChighE_telhighE_eleHW")
fd_eraEvlowE_eraEhighE_telhighE_eleHW<-getFd.era("fd_eraEvlowE_eraEhighE_telhighE_eleHW")

# # other.species.allopatric, these should all be lower
# # these make no sense, clysonymus from high west into erato high east or vicecersa?
# # cant serve as control as this might be true, if the high altitude allele in clysonymus is fixed in both sides (same)
# fd_eraCvlowE_eraChighE_clyEhighW_eleHW <-getFd.era("fd_eraCvlowE_eraChighE_clyEhighW_eleHW")
# fd_eraCvlowW_eraChighW_clyEhighE_eleHW<-getFd.era("fd_eraCvlowW_eraChighW_clyEhighE_eleHW")
# fd_eraEvlowE_eraEhighE_clyEhighW_eleHW <-getFd.era("fd_eraEvlowE_eraEhighE_clyEhighW_eleHW")
# fd_eraEvlowW_eraEhighW_clyEhighE_eleHW<-getFd.era("fd_eraEvlowW_eraEhighW_clyEhighE_eleHW")
# 
# # himera/telesiphe east into erato west (col/ecu)
# fd_eraCvlowW_eraChighW_himHE_eleHW<-getFd.era("fd_eraCvlowW_eraChighW_himHE_eleHW")
# fd_eraCvlowW_eraChighW_telhighE_eleHW<-getFd.era("fd_eraCvlowW_eraChighW_telhighE_eleHW")
# fd_eraEvlowW_eraEhighW_himHE_eleHW<-getFd.era("fd_eraEvlowW_eraEhighW_himHE_eleHW")
# fd_eraEvlowW_eraEhighW_telhighE_eleHW<-getFd.era("fd_eraEvlowW_eraEhighW_telhighE_eleHW")

# list files in environment that match, then creaste list
era.fd.list.files <-grep("fd_era",names(.GlobalEnv),value=TRUE)
era.fd.list<-do.call("list",mget(era.fd.list.files)); str(era.fd.list)
for (i in 1:length(era.fd.list)) {
  era.fd.list[[i]]$file.name.r <- c(era.fd.list.files[i]) }

# remove all individual files
rm(list = ls()[grep("fd_era*", ls())])

# bind and add info
era.fd.df <- rbindlist(era.fd.list)
names(era.comparisons.info)
era.fd.df<-merge(era.fd.df, era.comparisons.info[c(2,19:22,24)], by="file.name.r")

############################ 2.2. summarise each combo per hdr ############################
# start with parapatric sharing only, n=8 combos (2 tele, 2 him, 4 clys)
unique(era.fd.df$type.comparison)
era.other.species.parapatric.list <- unique(subset(era.fd.df, type.comparison=="other.species.parapatric" | type.comparison=="era.allopatric.from.hig")$file.name.r); era.other.species.parapatric.list

# max fd per hdr, max/min fdM per hdr
era.fd.df.shdr.para.east.summ.list<- list(); era.fd.df.shdr.para.west.summ.list<- list(); era.fd.df.shdr.allo.summ.list<- list()

for (i in 1:length(era.other.species.parapatric.list)) {
  combo <- era.other.species.parapatric.list[i]
  era.fd.df.shdr.para.east.summ.list[[i]] <- summarise(group_by(subset(era.fd.df, file.name.r==combo &shdr.para.east.id!=""),file.name.r, p1.newname,p2.newname,p3.newname, shdr.para.east.id),
                                                       no.windows=n(),fd.max=max(fd),fdm.max=max(fdM),fdm.min=min(fdM))
  era.fd.df.shdr.para.west.summ.list[[i]] <- summarise(group_by(subset(era.fd.df, file.name.r==combo &shdr.para.west.id!=""),file.name.r, p1.newname,p2.newname,p3.newname, shdr.para.west.id),
                                                       no.windows=n(),fd.max=max(fd),fdm.max=max(fdM),fdm.min=min(fdM))
  era.fd.df.shdr.allo.summ.list[[i]] <- summarise(group_by(subset(era.fd.df, file.name.r==combo &shdr.allo.id!=""),file.name.r, p1.newname,p2.newname,p3.newname, shdr.allo.id),
                                                       no.windows=n(),fd.max=max(fd),fdm.max=max(fdM),fdm.min=min(fdM))
}

# bind 
era.fd.df.shdr.para.east.summ.df <- rbindlist(era.fd.df.shdr.para.east.summ.list)
era.fd.df.shdr.para.west.summ.df <- rbindlist(era.fd.df.shdr.para.west.summ.list)
era.fd.df.shdr.allo.summ.df <- rbindlist(era.fd.df.shdr.allo.summ.list)


# each shdr.para.east/west.summ will have 8x3 new columns, for each combo max fd, max/min fdm is 10th-90th perc outlier? yes/no. 
# but! some of these combos more relevant than others, if it's an eastern combo it'll be more relevant for shdr.para.east than for shdr.para.west (compare!)

############################ 2.3. sims east/west para and allo COMBINED SHDR ###################
#### prep sims datasets ######
# extract the same number of shdr para east and of the same width, random locations get tajima min, tajima.5thperc.no.windows, tajima.5thperc.mean.taj
# first create intervals with real shdr. with Intervals we must work with BP.wg (not bp per scaff)
era.shdr.east.int <- reduce(Intervals(as.matrix(era.outliersEast[,c("start","end")]), closed=c(TRUE,TRUE), type="R"))
era.shdr.west.int <- reduce(Intervals(as.matrix(era.outliersWest[,c("start","end")]), closed=c(TRUE,TRUE), type="R"))

# random intervals for each data set
sim1k.era.east <- lapply(1:1000, function(x){rand_non_overlapping_intervals(era.genome_size,size(era.shdr.east.int))})
sim1k.era.west <- lapply(1:1000, function(x){rand_non_overlapping_intervals(era.genome_size,size(era.shdr.west.int))})


# run simulations on each ecuador/ colombia dataset separately
# careful because we are checking for overlaps between two sets of ranges (twisst is in 20kb windows)
sim.combo.fd.tmp.shdr.east.summ.list<-list(); sim.combo.fd.tmp.shdr.west.summ.list <-list()
sim.combo.fd.tmp.shdr.east.summ.list.list<-list(); sim.combo.fd.tmp.shdr.west.summ.list.list<-list()

for (i in 1:1000) {
  # prep random intervals for action
  sim1k.era.east.df <- as.data.frame(sim1k.era.east[i]); sim1k.era.east.df$ran.east.hdr.name <- paste("ran", i, ".hdr.", seq(1:nrow(sim1k.era.east.df )), sep = "")
  sim1k.era.west.df <- as.data.frame(sim1k.era.west[i]); sim1k.era.west.df$ran.west.hdr.name <- paste("ran", i, ".hdr.", seq(1:nrow(sim1k.era.west.df )), sep = "")

  # prepare all combo tmps, run for east/west/allo shdrs
  for (j in 1:length(era.other.species.parapatric.list)) {
    combo <- era.other.species.parapatric.list[j]
    combo.fd.tmp <- subset(era.fd.df, file.name.r==combo)[,c( "file.name.r", "start.WG", "end.WG", "fd","fdM")]
    
    #### east shdr ####
    # add ran hdr names to combo.tmp
    sim1k.era.east.dt <- as.data.table(sim1k.era.east.df[,1:2])
    names(sim1k.era.east.dt)[1:2]<-c("start.WG","end.WG") # to match fd tmp names
    
    setkey(sim1k.era.east.dt)
    # find overlaps
    sim.col.oveast<-foverlaps(x = combo.fd.tmp, y = sim1k.era.east.dt, which=T, by.x = names(sim1k.era.east.dt), type="any",nomatch=0L)
    # use keys to select corresponding shdr names
    combo.fd.tmp$ran.east.hdr.name <- NA; combo.fd.tmp$ran.east.hdr.name[sim.col.oveast$xid]<- sim1k.era.east.df[sim.col.oveast$yid,]$ran.east.hdr.name
    
    #### west shdr ####
    # add ran hdr names to combo.tmp
    sim1k.era.west.dt <- as.data.table(sim1k.era.west.df[,1:2])
    names(sim1k.era.west.dt)[1:2]<-c("start.WG","end.WG") # to match fd tmp names
    
    setkey(sim1k.era.west.dt)
    # find overlaps
    sim.col.oveast<-foverlaps(x = combo.fd.tmp, y = sim1k.era.west.dt, which=T, by.x = names(sim1k.era.west.dt), type="any",nomatch=0L)
    # use keys to select corresponding shdr names
    combo.fd.tmp$ran.west.hdr.name <- NA; combo.fd.tmp$ran.west.hdr.name[sim.col.oveast$xid]<- sim1k.era.west.df[sim.col.oveast$yid,]$ran.west.hdr.name
    
    #### summarise combo.fd.tmp by shdr type ####
    # get rid of NA ran hdr (background)
    sim.combo.fd.tmp.shdr.east.summ.list[[j]] <- summarise(group_by(subset(combo.fd.tmp, file.name.r==combo& ran.east.hdr.name!=""),file.name.r, ran.east.hdr.name),
                                                                no.windows=n(),fd.max=max(fd),fdm.max=max(fdM),fdm.min=min(fdM))
    sim.combo.fd.tmp.shdr.west.summ.list[[j]] <- summarise(group_by(subset(combo.fd.tmp, file.name.r==combo& ran.west.hdr.name!=""),file.name.r, ran.west.hdr.name ),
                                                                no.windows=n(),fd.max=max(fd),fdm.max=max(fdM),fdm.min=min(fdM))
    
  }
  
  # store combo summaries
  sim.combo.fd.tmp.shdr.east.summ.list.list[[i]] <- rbindlist(sim.combo.fd.tmp.shdr.east.summ.list)
  sim.combo.fd.tmp.shdr.west.summ.list.list[[i]] <- rbindlist(sim.combo.fd.tmp.shdr.west.summ.list)

}

sim.combo.fd.tmp.shdr.east.summ.list.list.df <- rbindlist(sim.combo.fd.tmp.shdr.east.summ.list.list); head(sim.combo.fd.tmp.shdr.east.summ.list.list.df)
sim.combo.fd.tmp.shdr.west.summ.list.list.df <- rbindlist(sim.combo.fd.tmp.shdr.west.summ.list.list); head(sim.combo.fd.tmp.shdr.west.summ.list.list.df)

# get rid of shdrs that are smaller than 100kb (only have 1 50kb window)
sim.combo.fd.tmp.shdr.east.summ.list.list.df.sub <- subset(sim.combo.fd.tmp.shdr.east.summ.list.list.df, no.windows>1)
sim.combo.fd.tmp.shdr.west.summ.list.list.df.sub <- subset(sim.combo.fd.tmp.shdr.west.summ.list.list.df, no.windows>1)

# store thresholds from simulations
se <- function(x) sd(x)/sqrt(length(x))
options(scipen = 999)
era.other.species.east.sims.threhsolds <- data.frame(combo.allo.type=NA, combo=NA, min.fdm.mean=NA, max.fdm.mean=NA, min.fdm.se=NA, max.fdm.se=NA, min.fdm.5th=NA, min.fdm.10th=NA, max.fdm.90th=NA, max.fdm.95th=NA ); era.other.species.east.sims.threhsolds 
for (i in 1:length(era.other.species.parapatric.list)) {
  combo <- era.other.species.parapatric.list[i]
  era.other.species.east.sims.threhsolds[i,]$combo.allo.type <- "shdr.east"
  era.other.species.east.sims.threhsolds[i,]$combo <- combo
  era.other.species.east.sims.threhsolds[i,]$min.fdm.mean <- mean(subset(sim.combo.fd.tmp.shdr.east.summ.list.list.df.sub, file.name.r==combo)$fdm.min) 
  era.other.species.east.sims.threhsolds[i,]$max.fdm.mean <- mean(subset(sim.combo.fd.tmp.shdr.east.summ.list.list.df.sub, file.name.r==combo)$fdm.max) 
  era.other.species.east.sims.threhsolds[i,]$min.fdm.se <- se(subset(sim.combo.fd.tmp.shdr.east.summ.list.list.df.sub, file.name.r==combo)$fdm.min) 
  era.other.species.east.sims.threhsolds[i,]$max.fdm.se <- se(subset(sim.combo.fd.tmp.shdr.east.summ.list.list.df.sub, file.name.r==combo)$fdm.max) 
  era.other.species.east.sims.threhsolds[i,]$min.fdm.5th <- quantile(subset(sim.combo.fd.tmp.shdr.east.summ.list.list.df.sub, file.name.r==combo)$fdm.min, 0.05) 
  era.other.species.east.sims.threhsolds[i,]$min.fdm.10th <- quantile(subset(sim.combo.fd.tmp.shdr.east.summ.list.list.df.sub, file.name.r==combo)$fdm.min, 0.1) 
  era.other.species.east.sims.threhsolds[i,]$max.fdm.90th <- quantile(subset(sim.combo.fd.tmp.shdr.east.summ.list.list.df.sub, file.name.r==combo)$fdm.max, 0.9) 
  era.other.species.east.sims.threhsolds[i,]$max.fdm.95th <- quantile(subset(sim.combo.fd.tmp.shdr.east.summ.list.list.df.sub, file.name.r==combo)$fdm.max, 0.95) 
}; era.other.species.east.sims.threhsolds

era.other.species.west.sims.threhsolds <- data.frame(combo.allo.type=NA, combo=NA, min.fdm.mean=NA, max.fdm.mean=NA, min.fdm.se=NA, max.fdm.se=NA, min.fdm.5th=NA, min.fdm.10th=NA, max.fdm.90th=NA, max.fdm.95th=NA ); era.other.species.sims.threhsolds
for (i in 1:length(era.other.species.parapatric.list)) {
  combo <- era.other.species.parapatric.list[i]
  era.other.species.west.sims.threhsolds[i,]$combo.allo.type <- "shdr.west"
  era.other.species.west.sims.threhsolds[i,]$combo <- combo
  era.other.species.west.sims.threhsolds[i,]$min.fdm.mean <- mean(subset(sim.combo.fd.tmp.shdr.west.summ.list.list.df.sub, file.name.r==combo)$fdm.min) 
  era.other.species.west.sims.threhsolds[i,]$max.fdm.mean <- mean(subset(sim.combo.fd.tmp.shdr.west.summ.list.list.df.sub, file.name.r==combo)$fdm.max) 
  era.other.species.west.sims.threhsolds[i,]$min.fdm.se <- se(subset(sim.combo.fd.tmp.shdr.west.summ.list.list.df.sub, file.name.r==combo)$fdm.min) 
  era.other.species.west.sims.threhsolds[i,]$max.fdm.se <- se(subset(sim.combo.fd.tmp.shdr.west.summ.list.list.df.sub, file.name.r==combo)$fdm.max) 
  era.other.species.west.sims.threhsolds[i,]$min.fdm.5th <- quantile(subset(sim.combo.fd.tmp.shdr.west.summ.list.list.df.sub, file.name.r==combo)$fdm.min, 0.05) 
  era.other.species.west.sims.threhsolds[i,]$min.fdm.10th <- quantile(subset(sim.combo.fd.tmp.shdr.west.summ.list.list.df.sub, file.name.r==combo)$fdm.min, 0.1) 
  era.other.species.west.sims.threhsolds[i,]$max.fdm.90th <- quantile(subset(sim.combo.fd.tmp.shdr.west.summ.list.list.df.sub, file.name.r==combo)$fdm.max, 0.9) 
  era.other.species.west.sims.threhsolds[i,]$max.fdm.95th <- quantile(subset(sim.combo.fd.tmp.shdr.west.summ.list.list.df.sub, file.name.r==combo)$fdm.max, 0.95) 
}; era.other.species.west.sims.threhsolds


era.other.species.sims.threhsolds <- rbind(era.other.species.east.sims.threhsolds, era.other.species.west.sims.threhsolds )
write.csv(era.other.species.sims.threhsolds, "local/data/joana.fd/era.within.and.other.species.1ksims.fdm.threhsolds.means.csv", row.names = T)


############################ 2.4. plot each combo, obs vs simulation thresholds ############################
###### plot fdm sims  east para SHDR ######
head(era.fd.df.shdr.para.east.summ.df); head(sim.combo.fd.tmp.shdr.para.east.summ.list.list.df)

fd.para.east.sims.plot <- list(); for (i in 1:length(era.other.species.parapatric.list)) {
  combo <- era.other.species.parapatric.list[i]
  fd.para.east.sims.plot[[i]] <-  ggplot(subset(era.fd.df, file.name.r==combo ), aes(x=fdM)) + 
  labs(colour="Altitude")+
  geom_vline(xintercept = 0, colour="lightgrey", lty="dashed")+
  geom_vline(xintercept = subset(era.fd.df.shdr.para.east.summ.df, file.name.r==combo & no.windows>1)$fdm.max, colour="#049E73", lty="solid", alpha=.9, size=0.5)+
    geom_vline(xintercept = subset(era.fd.df.shdr.para.east.summ.df, file.name.r==combo & no.windows>1 & shdr.para.east.id %in% subset(era.outliersEast, overlaps.with.inversion=="yes")$shdr.para.east.id )$fdm.max, 
               colour="yellow", lty="solid", alpha=1, size=0.7)+
  #geom_vline(xintercept = mean(subset(sim.combo.fd.tmp.shdr.para.east.summ.list.list.df.sub, file.name.r==combo)$fdm.max), colour="black", lty="solid", alpha=.9, size=1.5)+
  geom_vline(xintercept = quantile(subset(sim.combo.fd.tmp.shdr.para.east.summ.list.list.df.sub, file.name.r==combo)$fdm.max, 0.9), colour="black", lty="dashed", alpha=.9, size=1.5)+
  geom_vline(xintercept = quantile(subset(sim.combo.fd.tmp.shdr.para.east.summ.list.list.df.sub, file.name.r==combo)$fdm.max, 0.1), colour="black", lty="dashed", alpha=.9, size=1.5)+

  geom_vline(xintercept = subset(era.fd.df.shdr.para.east.summ.df, file.name.r==combo & no.windows>1)$fdm.min, colour="grey", lty="solid", alpha=.9, size=0.5)+
  geom_vline(xintercept = subset(era.fd.df.shdr.para.east.summ.df, file.name.r==combo & no.windows>1 & shdr.para.east.id %in% subset(era.outliersEast, overlaps.with.inversion=="yes")$shdr.para.east.id )$fdm.min, 
             colour="yellow", lty="solid", alpha=1, size=0.7)+
  #geom_vline(xintercept = mean(subset(sim.combo.fd.tmp.shdr.para.east.summ.list.list.df.sub, file.name.r==combo)$fdm.min), colour="black", lty="solid", alpha=.9, size=1.5)+
  geom_vline(xintercept = quantile(subset(sim.combo.fd.tmp.shdr.para.east.summ.list.list.df.sub, file.name.r==combo)$fdm.min, 0.9), colour="darkgrey", lty="dashed", alpha=.9, size=1.5)+
  geom_vline(xintercept = quantile(subset(sim.combo.fd.tmp.shdr.para.east.summ.list.list.df.sub, file.name.r==combo)$fdm.min, 0.1), colour="darkgrey", lty="dashed", alpha=.9, size=1.5)+
  geom_rect(aes(xmin=-.35, xmax=.35, ymin=100, ymax=120 ), fill="white")+
    # geom_density(fill="transparent", alpha=0, size=1)+
  xlab(expression(paste("Topo2-Topo3",  italic(" H. erato"))))+
  theme_classic()+
  theme(axis.line.x = element_line(color="black", size = 0.5),axis.line.y = element_line(color="black", size = 0.5),axis.title.y = element_text(size=14),
        plot.margin = unit(c(0.5,0.2,0.2,0.2), "lines"),axis.text = element_text(size=10), #axis.text.x = element_blank(),
        axis.text.y = element_text(size=10),axis.title.x = element_blank(), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_rect(fill = "transparent"), plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
        legend.background = element_rect(fill = "transparent"), legend.box.background = element_rect(fill = "transparent"), strip.background.y = element_blank(),
        strip.text = element_text(size=14))+
  scale_y_continuous(limits = c(0,120), expand = c(0,0))+
  scale_x_continuous(limits = c(-.35,.35), expand = c(0,0)) }; 

t <- c(fd.para.east.sims.plot[1],fd.para.east.sims.plot[5],fd.para.east.sims.plot[4], fd.para.east.sims.plot[8], fd.para.east.sims.plot[2], fd.para.east.sims.plot[6], fd.para.east.sims.plot[3], fd.para.east.sims.plot[7])
plot_grid(plotlist = t  , ncol = 2,
          labels = list(era.other.species.parapatric.list[1],era.other.species.parapatric.list[5],era.other.species.parapatric.list[4], era.other.species.parapatric.list[8], era.other.species.parapatric.list[2], era.other.species.parapatric.list[6], era.other.species.parapatric.list[3], era.other.species.parapatric.list[7]), label_x = -0.3)
ggsave("plots/joana.fd/era.other.sp.para.east.90th.fdm.max.min.png", width = 10, height = 9)

fd.para.west.sims.plot <- list(); for (i in 1:length(era.other.species.parapatric.list)) {
  combo <- era.other.species.parapatric.list[i]
  fd.para.west.sims.plot[[i]] <-  ggplot(subset(era.fd.df, file.name.r==combo ), aes(x=fdM)) + 
    labs(colour="Altitude")+
    geom_vline(xintercept = 0, colour="lightgrey", lty="dashed")+
    geom_vline(xintercept = subset(era.fd.df.shdr.para.west.summ.df, file.name.r==combo & no.windows>1)$fdm.max, colour="#0372B2", lty="solid", alpha=.9, size=0.5)+
    #geom_vline(xintercept = mean(subset(sim.combo.fd.tmp.shdr.para.west.summ.list.list.df.sub, file.name.r==combo)$fdm.max), colour="black", lty="solid", alpha=.9, size=1.5)+
    geom_vline(xintercept = quantile(subset(sim.combo.fd.tmp.shdr.para.west.summ.list.list.df.sub, file.name.r==combo)$fdm.max, 0.9), colour="black", lty="dashed", alpha=.9, size=1.5)+
    geom_vline(xintercept = quantile(subset(sim.combo.fd.tmp.shdr.para.west.summ.list.list.df.sub, file.name.r==combo)$fdm.max, 0.1), colour="black", lty="dashed", alpha=.9, size=1.5)+
    
    geom_vline(xintercept = subset(era.fd.df.shdr.para.west.summ.df, file.name.r==combo & no.windows>1)$fdm.min, colour="grey", lty="solid", alpha=.9, size=0.5)+
    # geom_vline(xintercept = mean(subset(sim.combo.fd.tmp.shdr.para.west.summ.list.list.df.sub, file.name.r==combo)$fdm.min), colour="black", lty="solid", alpha=.9, size=1.5)+
    geom_vline(xintercept = quantile(subset(sim.combo.fd.tmp.shdr.para.west.summ.list.list.df.sub, file.name.r==combo)$fdm.min, 0.9), colour="black", lty="dashed", alpha=.9, size=1.5)+
    geom_vline(xintercept = quantile(subset(sim.combo.fd.tmp.shdr.para.west.summ.list.list.df.sub, file.name.r==combo)$fdm.min, 0.1), colour="black", lty="dashed", alpha=.9, size=1.5)+
    geom_rect(aes(xmin=-.35, xmax=.35, ymin=100, ymax=120 ), fill="white")+
    
    geom_density(fill="transparent", alpha=0, size=1)+
    xlab(expression(paste("Topo2-Topo3",  italic(" H. erato"))))+
    theme_classic()+
    theme(axis.line.x = element_line(color="black", size = 0.5),axis.line.y = element_line(color="black", size = 0.5),axis.title.y = element_text(size=14),
          plot.margin = unit(c(0.5,0.2,0.2,0.2), "lines"),axis.text = element_text(size=10), #axis.text.x = element_blank(),
          axis.text.y = element_text(size=10),axis.title.x = element_blank(), 
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          panel.background = element_rect(fill = "transparent"), plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
          legend.background = element_rect(fill = "transparent"), legend.box.background = element_rect(fill = "transparent"), strip.background.y = element_blank(),
          strip.text = element_text(size=14))+
    scale_y_continuous(limits = c(0,120), expand = c(0,0))+
    scale_x_continuous(limits = c(-.35,.35), expand = c(0,0)) } ;
t <- c(fd.para.west.sims.plot[1],fd.para.west.sims.plot[5],fd.para.west.sims.plot[4], fd.para.west.sims.plot[8], fd.para.west.sims.plot[2], fd.para.west.sims.plot[6], fd.para.west.sims.plot[3], fd.para.west.sims.plot[7])
plot_grid(plotlist = t  , ncol = 2,
          labels = list(era.other.species.parapatric.list[1],era.other.species.parapatric.list[5],era.other.species.parapatric.list[4], era.other.species.parapatric.list[8], era.other.species.parapatric.list[2], era.other.species.parapatric.list[6], era.other.species.parapatric.list[3], era.other.species.parapatric.list[7]), label_x = -0.3)
#plot_grid(plotlist = fd.para.west.sims.plot, labels = era.other.species.parapatric.list, label_x = -0.3)
ggsave("plots/joana.fd/era.other.sp.para.west.90th.fdm.max.min.png", width = 10, height = 9)


fd.allo.sims.plot <- list()
for (i in 1:length(era.other.species.parapatric.list)) {
  combo <- era.other.species.parapatric.list[i]
  fd.allo.sims.plot[[i]] <-  ggplot(subset(era.fd.df, file.name.r==combo ), aes(x=fdM)) + 
    labs(colour="Altitude")+
    geom_vline(xintercept = 0, colour="lightgrey", lty="dashed")+
    geom_vline(xintercept = subset(era.fd.df.shdr.allo.summ.df, file.name.r==combo & no.windows>1)$fdm.max, colour="#D65D00", lty="solid", alpha=.9, size=0.5)+
    geom_vline(xintercept = subset(era.fd.df.shdr.allo.summ.df, file.name.r==combo & no.windows>1 & shdr.allo.id %in% subset(era.shdr.allo.all.summ, shdr.para.east %in% subset(era.shdr.para.east, overlaps.with.inversion=="yes")$shdr.para.east.id )$shdr.allo.id )$fdm.max, 
               colour="yellow", lty="solid", alpha=.9, size=0.7)+
    #geom_vline(xintercept = mean(subset(sim.combo.fd.tmp.shdr.allo.summ.list.list.df.sub, file.name.r==combo)$fdm.max), colour="black", lty="solid", alpha=.9, size=1.5)+
    geom_vline(xintercept = quantile(subset(sim.combo.fd.tmp.shdr.allo.summ.list.list.df.sub, file.name.r==combo)$fdm.max, 0.9), colour="black", lty="dashed", alpha=.9, size=1.5)+
    geom_vline(xintercept = quantile(subset(sim.combo.fd.tmp.shdr.allo.summ.list.list.df.sub, file.name.r==combo)$fdm.max, 0.1), colour="black", lty="dashed", alpha=.9, size=1.5)+
    geom_vline(xintercept = subset(era.fd.df.shdr.allo.summ.df, file.name.r==combo & no.windows>1)$fdm.min, colour="grey", lty="solid", alpha=.9, size=0.5)+
    geom_vline(xintercept = subset(era.fd.df.shdr.allo.summ.df, file.name.r==combo & no.windows>1 & shdr.allo.id %in% subset(era.shdr.allo.all.summ, shdr.para.east %in% subset(era.shdr.para.east, overlaps.with.inversion=="yes")$shdr.para.east.id )$shdr.allo.id )$fdm.min, 
               colour="yellow", lty="solid", alpha=1, size=0.7)+
    #geom_vline(xintercept = mean(subset(sim.combo.fd.tmp.shdr.allo.summ.list.list.df.sub, file.name.r==combo)$fdm.min), colour="black", lty="solid", alpha=.9, size=1.5)+
    geom_vline(xintercept = quantile(subset(sim.combo.fd.tmp.shdr.allo.summ.list.list.df.sub, file.name.r==combo)$fdm.min, 0.9), colour="darkgrey", lty="dashed", alpha=.9, size=1.5)+
    geom_vline(xintercept = quantile(subset(sim.combo.fd.tmp.shdr.allo.summ.list.list.df.sub, file.name.r==combo)$fdm.min, 0.1), colour="darkgrey", lty="dashed", alpha=.9, size=1.5)+
    
    geom_density(fill="transparent", alpha=0, size=1)+
    xlab(expression(paste("Topo2-Topo3",  italic(" H. erato"))))+
    theme_classic()+
    theme(axis.line.x = element_line(color="black", size = 0.5),axis.line.y = element_line(color="black", size = 0.5),axis.title.y = element_text(size=14),
          plot.margin = unit(c(0.5,0.2,0.2,0.2), "lines"),axis.text = element_text(size=10), #axis.text.x = element_blank(),
          axis.text.y = element_text(size=10),axis.title.x = element_blank(), 
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          panel.background = element_rect(fill = "transparent"), plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
          legend.background = element_rect(fill = "transparent"), legend.box.background = element_rect(fill = "transparent"), strip.background.y = element_blank(),
          strip.text = element_text(size=14))+
    scale_y_continuous(limits = c(0,120), expand = c(0,0))+
    scale_x_continuous(limits = c(-.3,.3), expand = c(0,0)) }
t <- c(fd.allo.sims.plot[1],fd.allo.sims.plot[5],fd.allo.sims.plot[4], fd.allo.sims.plot[8], fd.allo.sims.plot[2], fd.allo.sims.plot[6], fd.allo.sims.plot[3], fd.allo.sims.plot[7])
plot_grid(plotlist = t  , ncol = 2,
          labels = list(era.other.species.parapatric.list[1],era.other.species.parapatric.list[5],era.other.species.parapatric.list[4], era.other.species.parapatric.list[8], era.other.species.parapatric.list[2], era.other.species.parapatric.list[6], era.other.species.parapatric.list[3], era.other.species.parapatric.list[7]), label_x = -0.3)
ggsave("plots/joana.fd/era.other.sp.allo.90th.fdm.max.min.png", width = 10, height = 9)




###### plot fd sims  east para SHDR ######
head(era.fd.df.shdr.para.east.summ.df); head(sim.combo.fd.tmp.shdr.para.east.summ.list.list.df)

fd.para.east.sims.plot <- list(); for (i in 1:length(era.other.species.parapatric.list)) {
  combo <- era.other.species.parapatric.list[i]
  fd.para.east.sims.plot[[i]] <-  ggplot(subset(era.fd.df, file.name.r==combo ), aes(x=fd)) + 
    labs(colour="Altitude")+
    geom_vline(xintercept = 0, colour="lightgrey", lty="dashed")+
    geom_vline(xintercept = subset(era.fd.df.shdr.para.east.summ.df, file.name.r==combo & no.windows>1)$fd.max, colour="#049E73", lty="solid", alpha=.9, size=0.5)+
    geom_vline(xintercept = subset(era.fd.df.shdr.para.east.summ.df, file.name.r==combo & no.windows>1 & shdr.para.east.id %in% subset(era.outliersEast, overlaps.with.inversion=="yes")$shdr.para.east.id )$fd.max, 
               colour="yellow", lty="solid", alpha=1, size=0.7)+
    #geom_vline(xintercept = mean(subset(sim.combo.fd.tmp.shdr.para.east.summ.list.list.df.sub, file.name.r==combo)$fd.max), colour="black", lty="solid", alpha=.9, size=1.5)+
    geom_vline(xintercept = quantile(subset(sim.combo.fd.tmp.shdr.para.east.summ.list.list.df.sub, file.name.r==combo)$fd.max, 0.9), colour="black", lty="dashed", alpha=.9, size=1.5)+
    geom_vline(xintercept = quantile(subset(sim.combo.fd.tmp.shdr.para.east.summ.list.list.df.sub, file.name.r==combo)$fd.max, 0.1), colour="black", lty="dashed", alpha=.9, size=1.5)+
    geom_rect(aes(xmin=-0, xmax=.35, ymin=100, ymax=120 ), fill="white")+
    # geom_density(fill="transparent", alpha=0, size=1)+
    xlab(expression(paste("Topo2-Topo3",  italic(" H. erato"))))+
    theme_classic()+
    theme(axis.line.x = element_line(color="black", size = 0.5),axis.line.y = element_line(color="black", size = 0.5),axis.title.y = element_text(size=14),
          plot.margin = unit(c(0.5,0.2,0.2,0.2), "lines"),axis.text = element_text(size=10), #axis.text.x = element_blank(),
          axis.text.y = element_text(size=10),axis.title.x = element_blank(), 
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          panel.background = element_rect(fill = "transparent"), plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
          legend.background = element_rect(fill = "transparent"), legend.box.background = element_rect(fill = "transparent"), strip.background.y = element_blank(),
          strip.text = element_text(size=14))+
    scale_y_continuous(limits = c(0,120), expand = c(0,0))+
    scale_x_continuous(limits = c(0,.45), expand = c(0,0)) }; 

t <- c(fd.para.east.sims.plot[1],fd.para.east.sims.plot[5],fd.para.east.sims.plot[4], fd.para.east.sims.plot[8], fd.para.east.sims.plot[2], fd.para.east.sims.plot[6], fd.para.east.sims.plot[3], fd.para.east.sims.plot[7])
plot_grid(plotlist = t  , ncol = 2,
          labels = list(era.other.species.parapatric.list[1],era.other.species.parapatric.list[5],era.other.species.parapatric.list[4], era.other.species.parapatric.list[8], era.other.species.parapatric.list[2], era.other.species.parapatric.list[6], 
                        era.other.species.parapatric.list[3], era.other.species.parapatric.list[7]), label_x = -0.2)
ggsave("plots/joana.fd/era.fd.other.sp.para.east.90th.fdm.max.min.png", width = 10, height = 9)


fd.para.west.sims.plot <- list(); for (i in 1:length(era.other.species.parapatric.list)) {
  combo <- era.other.species.parapatric.list[i]
  fd.para.west.sims.plot[[i]] <-  ggplot(subset(era.fd.df, file.name.r==combo ), aes(x=fd)) + 
    labs(colour="Altitude")+
    geom_vline(xintercept = 0, colour="lightgrey", lty="dashed")+
    geom_vline(xintercept = subset(era.fd.df.shdr.para.west.summ.df, file.name.r==combo & no.windows>1)$fd.max, colour="#0372B2", lty="solid", alpha=.9, size=0.5)+
    #geom_vline(xintercept = mean(subset(sim.combo.fd.tmp.shdr.para.west.summ.list.list.df.sub, file.name.r==combo)$fd.max), colour="black", lty="solid", alpha=.9, size=1.5)+
    geom_vline(xintercept = quantile(subset(sim.combo.fd.tmp.shdr.para.west.summ.list.list.df.sub, file.name.r==combo)$fd.max, 0.9), colour="black", lty="dashed", alpha=.9, size=1.5)+
    geom_vline(xintercept = quantile(subset(sim.combo.fd.tmp.shdr.para.west.summ.list.list.df.sub, file.name.r==combo)$fd.max, 0.1), colour="black", lty="dashed", alpha=.9, size=1.5)+
    geom_rect(aes(xmin=-0, xmax=.35, ymin=100, ymax=120 ), fill="white")+
    
    geom_density(fill="transparent", alpha=0, size=1)+
    xlab(expression(paste("Topo2-Topo3",  italic(" H. erato"))))+
    theme_classic()+
    theme(axis.line.x = element_line(color="black", size = 0.5),axis.line.y = element_line(color="black", size = 0.5),axis.title.y = element_text(size=14),
          plot.margin = unit(c(0.5,0.2,0.2,0.2), "lines"),axis.text = element_text(size=10), #axis.text.x = element_blank(),
          axis.text.y = element_text(size=10),axis.title.x = element_blank(), 
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          panel.background = element_rect(fill = "transparent"), plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
          legend.background = element_rect(fill = "transparent"), legend.box.background = element_rect(fill = "transparent"), strip.background.y = element_blank(),
          strip.text = element_text(size=14))+
    scale_y_continuous(limits = c(0,120), expand = c(0,0))+
    scale_x_continuous(limits = c(0,.35), expand = c(0,0)) } ;
t <- c(fd.para.west.sims.plot[1],fd.para.west.sims.plot[5],fd.para.west.sims.plot[4], fd.para.west.sims.plot[8], fd.para.west.sims.plot[2], fd.para.west.sims.plot[6], fd.para.west.sims.plot[3], fd.para.west.sims.plot[7])
plot_grid(plotlist = t  , ncol = 2,
          labels = list(era.other.species.parapatric.list[1],era.other.species.parapatric.list[5],era.other.species.parapatric.list[4], era.other.species.parapatric.list[8], 
                        era.other.species.parapatric.list[2], era.other.species.parapatric.list[6], era.other.species.parapatric.list[3], era.other.species.parapatric.list[7]), label_x = -0.2)
ggsave("plots/joana.fd/era.fd.other.sp.para.west.90th.fdm.max.min.png", width = 10, height = 9)


fd.allo.sims.plot <- list()
for (i in 1:length(era.other.species.parapatric.list)) {
  combo <- era.other.species.parapatric.list[i]
  fd.allo.sims.plot[[i]] <-  ggplot(subset(era.fd.df, file.name.r==combo ), aes(x=fd)) + 
    labs(colour="Altitude")+
    geom_vline(xintercept = 0, colour="lightgrey", lty="dashed")+
    geom_vline(xintercept = subset(era.fd.df.shdr.allo.summ.df, file.name.r==combo & no.windows>1)$fd.max, colour="#D65D00", lty="solid", alpha=.9, size=0.5)+
    geom_vline(xintercept = subset(era.fd.df.shdr.allo.summ.df, file.name.r==combo & no.windows>1 & shdr.allo.id %in% subset(era.shdr.allo.all.summ, shdr.para.east %in% subset(era.shdr.para.east, overlaps.with.inversion=="yes")$shdr.para.east.id )$shdr.allo.id )$fd.max, 
               colour="yellow", lty="solid", alpha=.9, size=0.7)+
    #geom_vline(xintercept = mean(subset(sim.combo.fd.tmp.shdr.allo.summ.list.list.df.sub, file.name.r==combo)$fd.max), colour="black", lty="solid", alpha=.9, size=1.5)+
    geom_vline(xintercept = quantile(subset(sim.combo.fd.tmp.shdr.allo.summ.list.list.df.sub, file.name.r==combo)$fd.max, 0.9), colour="black", lty="dashed", alpha=.9, size=1.5)+
    geom_vline(xintercept = quantile(subset(sim.combo.fd.tmp.shdr.allo.summ.list.list.df.sub, file.name.r==combo)$fd.max, 0.1), colour="black", lty="dashed", alpha=.9, size=1.5)+
    geom_density(fill="transparent", alpha=0, size=1)+
    geom_rect(aes(xmin=-0, xmax=.35, ymin=100, ymax=120 ), fill="white")+
    
    xlab(expression(paste("Topo2-Topo3",  italic(" H. erato"))))+
    theme_classic()+
    theme(axis.line.x = element_line(color="black", size = 0.5),axis.line.y = element_line(color="black", size = 0.5),axis.title.y = element_text(size=14),
          plot.margin = unit(c(0.5,0.2,0.2,0.2), "lines"),axis.text = element_text(size=10), #axis.text.x = element_blank(),
          axis.text.y = element_text(size=10),axis.title.x = element_blank(), 
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          panel.background = element_rect(fill = "transparent"), plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
          legend.background = element_rect(fill = "transparent"), legend.box.background = element_rect(fill = "transparent"), strip.background.y = element_blank(),
          strip.text = element_text(size=14))+
    scale_y_continuous(limits = c(0,120), expand = c(0,0))+
    scale_x_continuous(limits = c(0,.3), expand = c(0,0)) }
t <- c(fd.allo.sims.plot[1],fd.allo.sims.plot[5],fd.allo.sims.plot[4], fd.allo.sims.plot[8], fd.allo.sims.plot[2], fd.allo.sims.plot[6], fd.allo.sims.plot[3], fd.allo.sims.plot[7])
plot_grid(plotlist = t  , ncol = 2,
          labels = list(era.other.species.parapatric.list[1],era.other.species.parapatric.list[5],era.other.species.parapatric.list[4], era.other.species.parapatric.list[8], 
                        era.other.species.parapatric.list[2], era.other.species.parapatric.list[6], era.other.species.parapatric.list[3], era.other.species.parapatric.list[7]), label_x = -0.2)
ggsave("plots/joana.fd/era.fd.other.sp.allo.90th.fdm.max.min.png", width = 10, height = 9)




############################ 2.5. assess which hdr >< simulation thresholds ############################
era.other.species.parapatric.list <- unique(subset(era.fd.df, type.comparison=="other.species.parapatric")$file.name.r); era.other.species.parapatric.list

# max fd per hdr, max/min fdM per hdr
era.fd.df.shdr.para.east.summ.df <- as.data.frame(era.fd.df.shdr.para.east.summ.df)
era.fd.df.shdr.para.east.summ.df$fdm.max.more.90th.perc.sims <- "NA"; era.fd.df.shdr.para.east.summ.df$fdm.max.more.95th.perc.sims <- "NA"; era.fd.df.shdr.para.east.summ.df$fd.max.more.90th.perc.sims <- "NA"; era.fd.df.shdr.para.east.summ.df$fd.max.more.95th.perc.sims <- "NA"
era.fd.df.shdr.para.east.summ.df$fdm.min.less.10th.perc.sims <- "NA"; era.fd.df.shdr.para.east.summ.df$fdm.min.less.5th.perc.sims <- "NA"; era.fd.df.shdr.para.east.summ.df$fd.min.less.10th.perc.sims <- "NA"; era.fd.df.shdr.para.east.summ.df$fd.min.less.5th.perc.sims <- "NA"
era.fd.df.shdr.para.east.summ.df$comb.new.name <- paste(era.fd.df.shdr.para.east.summ.df$p1.newname, era.fd.df.shdr.para.east.summ.df$p2.newname, era.fd.df.shdr.para.east.summ.df$p3.newname, sep = "_")

era.fd.df.shdr.para.west.summ.df <- as.data.frame(era.fd.df.shdr.para.west.summ.df)
era.fd.df.shdr.para.west.summ.df$fdm.max.more.90th.perc.sims <- "NA"; era.fd.df.shdr.para.west.summ.df$fdm.max.more.95th.perc.sims <- "NA"; era.fd.df.shdr.para.west.summ.df$fd.max.more.90th.perc.sims <- "NA"; era.fd.df.shdr.para.west.summ.df$fd.max.more.95th.perc.sims <- "NA"
era.fd.df.shdr.para.west.summ.df$fdm.min.less.10th.perc.sims <- "NA"; era.fd.df.shdr.para.west.summ.df$fdm.min.less.5th.perc.sims <- "NA"; era.fd.df.shdr.para.west.summ.df$fd.min.less.10th.perc.sims <- "NA"; era.fd.df.shdr.para.west.summ.df$fd.min.less.5th.perc.sims <- "NA"
era.fd.df.shdr.para.west.summ.df$comb.new.name <- paste(era.fd.df.shdr.para.west.summ.df$p1.newname, era.fd.df.shdr.para.west.summ.df$p2.newname, era.fd.df.shdr.para.west.summ.df$p3.newname, sep = "_")

era.fd.df.shdr.allo.summ.df <- as.data.frame(era.fd.df.shdr.allo.summ.df)
era.fd.df.shdr.allo.summ.df$fdm.max.more.90th.perc.sims <- "NA"; era.fd.df.shdr.allo.summ.df$fdm.max.more.95th.perc.sims <- "NA"; era.fd.df.shdr.allo.summ.df$fd.max.more.90th.perc.sims <- "NA"; era.fd.df.shdr.allo.summ.df$fd.max.more.95th.perc.sims <- "NA"
era.fd.df.shdr.allo.summ.df$fdm.min.less.10th.perc.sims <- "NA"; era.fd.df.shdr.allo.summ.df$fdm.min.less.5th.perc.sims <- "NA"; era.fd.df.shdr.allo.summ.df$fd.min.less.10th.perc.sims <- "NA"; era.fd.df.shdr.allo.summ.df$fd.min.less.5th.perc.sims <- "NA"
era.fd.df.shdr.allo.summ.df$comb.new.name <- paste(era.fd.df.shdr.allo.summ.df$p1.newname, era.fd.df.shdr.allo.summ.df$p2.newname, era.fd.df.shdr.allo.summ.df$p3.newname, sep = "_")


for (i in 1:length(era.other.species.parapatric.list)) {
  combo <- era.other.species.parapatric.list[i]
  # east fdm.max, fdm.min, fd.max
  era.fd.df.shdr.para.east.summ.df[era.fd.df.shdr.para.east.summ.df$file.name.r==combo,]$fdm.max.more.90th.perc.sims <- if_else(era.fd.df.shdr.para.east.summ.df[era.fd.df.shdr.para.east.summ.df$file.name.r==combo,]$fdm.max > 
                                                                            quantile(subset(sim.combo.fd.tmp.shdr.para.east.summ.list.list.df.sub, file.name.r==combo)$fdm.max, 0.9), "yes", "no")
  era.fd.df.shdr.para.east.summ.df[era.fd.df.shdr.para.east.summ.df$file.name.r==combo,]$fdm.max.more.95th.perc.sims <- if_else(era.fd.df.shdr.para.east.summ.df[era.fd.df.shdr.para.east.summ.df$file.name.r==combo,]$fdm.max > 
                                                                                                                                  quantile(subset(sim.combo.fd.tmp.shdr.para.east.summ.list.list.df.sub, file.name.r==combo)$fdm.max, 0.95), "yes", "no")
  era.fd.df.shdr.para.east.summ.df[era.fd.df.shdr.para.east.summ.df$file.name.r==combo,]$fdm.min.less.10th.perc.sims  <- if_else(era.fd.df.shdr.para.east.summ.df[era.fd.df.shdr.para.east.summ.df$file.name.r==combo,]$fdm.min < 
                                                                                                                                  quantile(subset(sim.combo.fd.tmp.shdr.para.east.summ.list.list.df.sub, file.name.r==combo)$fdm.min, 0.1), "yes", "no")
  era.fd.df.shdr.para.east.summ.df[era.fd.df.shdr.para.east.summ.df$file.name.r==combo,]$fdm.min.less.5th.perc.sims  <- if_else(era.fd.df.shdr.para.east.summ.df[era.fd.df.shdr.para.east.summ.df$file.name.r==combo,]$fdm.min < 
                                                                                                                                   quantile(subset(sim.combo.fd.tmp.shdr.para.east.summ.list.list.df.sub, file.name.r==combo)$fdm.min, 0.05), "yes", "no")
  era.fd.df.shdr.para.east.summ.df[era.fd.df.shdr.para.east.summ.df$file.name.r==combo,]$fd.max.more.90th.perc.sims <- if_else(era.fd.df.shdr.para.east.summ.df[era.fd.df.shdr.para.east.summ.df$file.name.r==combo,]$fd.max > 
                                                                                                                                  quantile(subset(sim.combo.fd.tmp.shdr.para.east.summ.list.list.df.sub, file.name.r==combo)$fd.max, 0.9), "yes", "no")
  era.fd.df.shdr.para.east.summ.df[era.fd.df.shdr.para.east.summ.df$file.name.r==combo,]$fd.max.more.95th.perc.sims <- if_else(era.fd.df.shdr.para.east.summ.df[era.fd.df.shdr.para.east.summ.df$file.name.r==combo,]$fd.max > 
                                                                                                                                  quantile(subset(sim.combo.fd.tmp.shdr.para.east.summ.list.list.df.sub, file.name.r==combo)$fd.max, 0.95), "yes", "no")
  
  # west fdm.max, fdm.min, fd.max
  era.fd.df.shdr.para.west.summ.df[era.fd.df.shdr.para.west.summ.df$file.name.r==combo,]$fdm.max.more.90th.perc.sims <- if_else(era.fd.df.shdr.para.west.summ.df[era.fd.df.shdr.para.west.summ.df$file.name.r==combo,]$fdm.max > 
                                                                                                                                  quantile(subset(sim.combo.fd.tmp.shdr.para.west.summ.list.list.df.sub, file.name.r==combo)$fdm.max, 0.9), "yes", "no")
  era.fd.df.shdr.para.west.summ.df[era.fd.df.shdr.para.west.summ.df$file.name.r==combo,]$fdm.max.more.95th.perc.sims <- if_else(era.fd.df.shdr.para.west.summ.df[era.fd.df.shdr.para.west.summ.df$file.name.r==combo,]$fdm.max > 
                                                                                                                                  quantile(subset(sim.combo.fd.tmp.shdr.para.west.summ.list.list.df.sub, file.name.r==combo)$fdm.max, 0.95), "yes", "no")
  era.fd.df.shdr.para.west.summ.df[era.fd.df.shdr.para.west.summ.df$file.name.r==combo,]$fdm.min.less.10th.perc.sims  <- if_else(era.fd.df.shdr.para.west.summ.df[era.fd.df.shdr.para.west.summ.df$file.name.r==combo,]$fdm.min < 
                                                                                                                                   quantile(subset(sim.combo.fd.tmp.shdr.para.west.summ.list.list.df.sub, file.name.r==combo)$fdm.min, 0.1), "yes", "no")
  era.fd.df.shdr.para.west.summ.df[era.fd.df.shdr.para.west.summ.df$file.name.r==combo,]$fdm.min.less.5th.perc.sims  <- if_else(era.fd.df.shdr.para.west.summ.df[era.fd.df.shdr.para.west.summ.df$file.name.r==combo,]$fdm.min < 
                                                                                                                                  quantile(subset(sim.combo.fd.tmp.shdr.para.west.summ.list.list.df.sub, file.name.r==combo)$fdm.min, 0.05), "yes", "no")
  era.fd.df.shdr.para.west.summ.df[era.fd.df.shdr.para.west.summ.df$file.name.r==combo,]$fd.max.more.90th.perc.sims <- if_else(era.fd.df.shdr.para.west.summ.df[era.fd.df.shdr.para.west.summ.df$file.name.r==combo,]$fd.max > 
                                                                                                                                 quantile(subset(sim.combo.fd.tmp.shdr.para.west.summ.list.list.df.sub, file.name.r==combo)$fd.max, 0.9), "yes", "no")
  era.fd.df.shdr.para.west.summ.df[era.fd.df.shdr.para.west.summ.df$file.name.r==combo,]$fd.max.more.95th.perc.sims <- if_else(era.fd.df.shdr.para.west.summ.df[era.fd.df.shdr.para.west.summ.df$file.name.r==combo,]$fd.max > 
                                                                                                                                 quantile(subset(sim.combo.fd.tmp.shdr.para.west.summ.list.list.df.sub, file.name.r==combo)$fd.max, 0.95), "yes", "no")
  # allo fdm.max, fdm.min, fd.max
  era.fd.df.shdr.allo.summ.df[era.fd.df.shdr.allo.summ.df$file.name.r==combo,]$fdm.max.more.90th.perc.sims <- if_else(era.fd.df.shdr.allo.summ.df[era.fd.df.shdr.allo.summ.df$file.name.r==combo,]$fdm.max > 
                                                                                                                                  quantile(subset(sim.combo.fd.tmp.shdr.allo.summ.list.list.df.sub, file.name.r==combo)$fdm.max, 0.9), "yes", "no")
  era.fd.df.shdr.allo.summ.df[era.fd.df.shdr.allo.summ.df$file.name.r==combo,]$fdm.max.more.95th.perc.sims <- if_else(era.fd.df.shdr.allo.summ.df[era.fd.df.shdr.allo.summ.df$file.name.r==combo,]$fdm.max > 
                                                                                                                                  quantile(subset(sim.combo.fd.tmp.shdr.allo.summ.list.list.df.sub, file.name.r==combo)$fdm.max, 0.95), "yes", "no")
  era.fd.df.shdr.allo.summ.df[era.fd.df.shdr.allo.summ.df$file.name.r==combo,]$fdm.min.less.10th.perc.sims  <- if_else(era.fd.df.shdr.allo.summ.df[era.fd.df.shdr.allo.summ.df$file.name.r==combo,]$fdm.min < 
                                                                                                                                   quantile(subset(sim.combo.fd.tmp.shdr.allo.summ.list.list.df.sub, file.name.r==combo)$fdm.min, 0.1), "yes", "no")
  era.fd.df.shdr.allo.summ.df[era.fd.df.shdr.allo.summ.df$file.name.r==combo,]$fdm.min.less.5th.perc.sims  <- if_else(era.fd.df.shdr.allo.summ.df[era.fd.df.shdr.allo.summ.df$file.name.r==combo,]$fdm.min < 
                                                                                                                                  quantile(subset(sim.combo.fd.tmp.shdr.allo.summ.list.list.df.sub, file.name.r==combo)$fdm.min, 0.05), "yes", "no")
  
  era.fd.df.shdr.allo.summ.df[era.fd.df.shdr.allo.summ.df$file.name.r==combo,]$fd.max.more.90th.perc.sims <- if_else(era.fd.df.shdr.allo.summ.df[era.fd.df.shdr.allo.summ.df$file.name.r==combo,]$fd.max > 
                                                                                                                                 quantile(subset(sim.combo.fd.tmp.shdr.allo.summ.list.list.df.sub, file.name.r==combo)$fd.max, 0.9), "yes", "no")
  era.fd.df.shdr.allo.summ.df[era.fd.df.shdr.allo.summ.df$file.name.r==combo,]$fd.max.more.95th.perc.sims <- if_else(era.fd.df.shdr.allo.summ.df[era.fd.df.shdr.allo.summ.df$file.name.r==combo,]$fd.max > 
                                                                                                                                 quantile(subset(sim.combo.fd.tmp.shdr.allo.summ.list.list.df.sub, file.name.r==combo)$fd.max, 0.95), "yes", "no")
  
  
}

#### add to hdr outlier summaries #### 
# add allo 
era.fd.df.shdr.allo.summ.df$shdr.para.east.id <- era.shdr.allo.all.summ$shdr.para.east[match(era.fd.df.shdr.allo.summ.df$shdr.allo.id,era.shdr.allo.all.summ$shdr.allo.id )]
era.fd.df.shdr.allo.summ.df$shdr.para.west.id <- era.shdr.allo.all.summ$shdr.para.west[match(era.fd.df.shdr.allo.summ.df$shdr.allo.id,era.shdr.allo.all.summ$shdr.allo.id )]

shdr.era.east.outlier.df <- read.csv("local/data/shdr.summ/shdr.era.east.outlier.df.csv")
shdr.era.west.outlier.df <- read.csv("local/data/shdr.summ/shdr.era.west.outlier.df.csv")

for (i in 1:length(era.other.species.parapatric.list)) {
  combo <- era.other.species.parapatric.list[i]
  # 48 new columns, 8 combs x 2 (90,95th) x 3(fdm.max, fdm.min and fdmax), split into east west and allo, but add to para.east and para.west only
  ##### store fdm max min values per combo and cline ######
  shdr.era.east.outlier.df[paste("fdm.max", combo, sep="_")] <- NA
  shdr.era.east.outlier.df[paste("fdm.min", combo, sep="_")] <- NA
  # east
  shdr.era.east.outlier.df[shdr.era.east.outlier.df$is.allo.hdr=="no",][paste("fdm.max", combo, sep="_")] <- era.fd.df.shdr.para.east.summ.df[era.fd.df.shdr.para.east.summ.df$file.name.r==combo,]$fdm.max[
    match(shdr.era.east.outlier.df[shdr.era.east.outlier.df$is.allo.hdr=="no",]$shdr.para.east.id, era.fd.df.shdr.para.east.summ.df[era.fd.df.shdr.para.east.summ.df$file.name.r==combo,]$shdr.para.east.id )]
  shdr.era.east.outlier.df[shdr.era.east.outlier.df$is.allo.hdr=="no",][paste("fdm.min", combo, sep="_")] <- era.fd.df.shdr.para.east.summ.df[era.fd.df.shdr.para.east.summ.df$file.name.r==combo,]$fdm.min[
    match(shdr.era.east.outlier.df[shdr.era.east.outlier.df$is.allo.hdr=="no",]$shdr.para.east.id, era.fd.df.shdr.para.east.summ.df[era.fd.df.shdr.para.east.summ.df$file.name.r==combo,]$shdr.para.east.id )]
  
  # allo east
  shdr.era.east.outlier.df[shdr.era.east.outlier.df$is.allo.hdr=="yes",][paste("fdm.max", combo, sep="_")] <- era.fd.df.shdr.allo.summ.df[era.fd.df.shdr.allo.summ.df$file.name.r==combo,]$fdm.max[
    match(shdr.era.east.outlier.df[shdr.era.east.outlier.df$is.allo.hdr=="yes",]$shdr.para.east.id, era.fd.df.shdr.allo.summ.df[era.fd.df.shdr.allo.summ.df$file.name.r==combo,]$shdr.para.east.id )]
  shdr.era.east.outlier.df[shdr.era.east.outlier.df$is.allo.hdr=="yes",][paste("fdm.min", combo, sep="_")] <- era.fd.df.shdr.allo.summ.df[era.fd.df.shdr.allo.summ.df$file.name.r==combo,]$fdm.min[
    match(shdr.era.east.outlier.df[shdr.era.east.outlier.df$is.allo.hdr=="yes",]$shdr.para.east.id, era.fd.df.shdr.allo.summ.df[era.fd.df.shdr.allo.summ.df$file.name.r==combo,]$shdr.para.east.id )]
  
  # west
  shdr.era.west.outlier.df[paste("fdm.max", combo, sep="_")] <- NA
  shdr.era.west.outlier.df[paste("fdm.min", combo, sep="_")] <- NA
  shdr.era.west.outlier.df[shdr.era.west.outlier.df$is.allo.hdr=="no",][paste("fdm.max", combo, sep="_")] <- era.fd.df.shdr.para.west.summ.df[era.fd.df.shdr.para.west.summ.df$file.name.r==combo,]$fdm.max[
    match(shdr.era.west.outlier.df[shdr.era.west.outlier.df$is.allo.hdr=="no",]$shdr.para.west.id, era.fd.df.shdr.para.west.summ.df[era.fd.df.shdr.para.west.summ.df$file.name.r==combo,]$shdr.para.west.id )]
  shdr.era.west.outlier.df[shdr.era.west.outlier.df$is.allo.hdr=="no",][paste("fdm.min", combo, sep="_")] <- era.fd.df.shdr.para.west.summ.df[era.fd.df.shdr.para.west.summ.df$file.name.r==combo,]$fdm.min[
    match(shdr.era.west.outlier.df[shdr.era.west.outlier.df$is.allo.hdr=="no",]$shdr.para.west.id, era.fd.df.shdr.para.west.summ.df[era.fd.df.shdr.para.west.summ.df$file.name.r==combo,]$shdr.para.west.id )]
  
  # allo west
  shdr.era.west.outlier.df[shdr.era.west.outlier.df$is.allo.hdr=="yes",][paste("fdm.max", combo, sep="_")] <- era.fd.df.shdr.allo.summ.df[era.fd.df.shdr.allo.summ.df$file.name.r==combo,]$fdm.max[
    match(shdr.era.west.outlier.df[shdr.era.west.outlier.df$is.allo.hdr=="yes",]$shdr.para.west.id, era.fd.df.shdr.allo.summ.df[era.fd.df.shdr.allo.summ.df$file.name.r==combo,]$shdr.para.west.id )]
  shdr.era.west.outlier.df[shdr.era.west.outlier.df$is.allo.hdr=="yes",][paste("fdm.min", combo, sep="_")] <- era.fd.df.shdr.allo.summ.df[era.fd.df.shdr.allo.summ.df$file.name.r==combo,]$fdm.min[
    match(shdr.era.west.outlier.df[shdr.era.west.outlier.df$is.allo.hdr=="yes",]$shdr.para.west.id, era.fd.df.shdr.allo.summ.df[era.fd.df.shdr.allo.summ.df$file.name.r==combo,]$shdr.para.west.id )]
  
  ##### store whether they are outliers at diff threholds or not, yes no ######
  # east
  shdr.era.east.outlier.df[paste("fdm.max.more.90th.perc.sims", combo, sep="_")] <- NA; shdr.era.east.outlier.df[paste("fdm.max.more.95th.perc.sims", combo, sep="_")] <- NA
  shdr.era.east.outlier.df[paste("fdm.min.less.10th.perc.sims", combo, sep="_")] <- NA; shdr.era.east.outlier.df[paste("fdm.min.less.5th.perc.sims", combo, sep="_")] <- NA
  shdr.era.east.outlier.df[shdr.era.east.outlier.df$is.allo.hdr=="no",][paste("fdm.max.more.90th.perc.sims", combo, sep="_")] <- era.fd.df.shdr.para.east.summ.df[era.fd.df.shdr.para.east.summ.df$file.name.r==combo,]$fdm.max.more.90th.perc.sims[
    match(shdr.era.east.outlier.df[shdr.era.east.outlier.df$is.allo.hdr=="no",]$shdr.para.east.id, era.fd.df.shdr.para.east.summ.df[era.fd.df.shdr.para.east.summ.df$file.name.r==combo,]$shdr.para.east.id )]
  shdr.era.east.outlier.df[shdr.era.east.outlier.df$is.allo.hdr=="no",][paste("fdm.max.more.95th.perc.sims", combo, sep="_")] <- era.fd.df.shdr.para.east.summ.df[era.fd.df.shdr.para.east.summ.df$file.name.r==combo,]$fdm.max.more.95th.perc.sims[
    match(shdr.era.east.outlier.df[shdr.era.east.outlier.df$is.allo.hdr=="no",]$shdr.para.east.id, era.fd.df.shdr.para.east.summ.df[era.fd.df.shdr.para.east.summ.df$file.name.r==combo,]$shdr.para.east.id )]
  shdr.era.east.outlier.df[shdr.era.east.outlier.df$is.allo.hdr=="no",][paste("fdm.min.less.10th.perc.sims", combo, sep="_")] <- era.fd.df.shdr.para.east.summ.df[era.fd.df.shdr.para.east.summ.df$file.name.r==combo,]$fdm.min.less.10th.perc.sims[
    match(shdr.era.east.outlier.df[shdr.era.east.outlier.df$is.allo.hdr=="no",]$shdr.para.east.id, era.fd.df.shdr.para.east.summ.df[era.fd.df.shdr.para.east.summ.df$file.name.r==combo,]$shdr.para.east.id )]
  shdr.era.east.outlier.df[shdr.era.east.outlier.df$is.allo.hdr=="no",][paste("fdm.min.less.5th.perc.sims", combo, sep="_")] <- era.fd.df.shdr.para.east.summ.df[era.fd.df.shdr.para.east.summ.df$file.name.r==combo,]$fdm.min.less.5th.perc.sims[
    match(shdr.era.east.outlier.df[shdr.era.east.outlier.df$is.allo.hdr=="no",]$shdr.para.east.id, era.fd.df.shdr.para.east.summ.df[era.fd.df.shdr.para.east.summ.df$file.name.r==combo,]$shdr.para.east.id )]
  
  # allo east
  shdr.era.east.outlier.df[shdr.era.east.outlier.df$is.allo.hdr=="yes",][paste("fdm.max.more.90th.perc.sims", combo, sep="_")] <- era.fd.df.shdr.allo.summ.df[era.fd.df.shdr.allo.summ.df$file.name.r==combo,]$fdm.max.more.90th.perc.sims[
    match(shdr.era.east.outlier.df[shdr.era.east.outlier.df$is.allo.hdr=="yes",]$shdr.para.east.id, era.fd.df.shdr.allo.summ.df[era.fd.df.shdr.allo.summ.df$file.name.r==combo,]$shdr.para.east.id )]
  shdr.era.east.outlier.df[shdr.era.east.outlier.df$is.allo.hdr=="yes",][paste("fdm.max.more.95th.perc.sims", combo, sep="_")] <- era.fd.df.shdr.allo.summ.df[era.fd.df.shdr.allo.summ.df$file.name.r==combo,]$fdm.max.more.95th.perc.sims[
    match(shdr.era.east.outlier.df[shdr.era.east.outlier.df$is.allo.hdr=="yes",]$shdr.para.east.id, era.fd.df.shdr.allo.summ.df[era.fd.df.shdr.allo.summ.df$file.name.r==combo,]$shdr.para.east.id )]
  shdr.era.east.outlier.df[shdr.era.east.outlier.df$is.allo.hdr=="yes",][paste("fdm.min.less.10th.perc.sims", combo, sep="_")] <- era.fd.df.shdr.allo.summ.df[era.fd.df.shdr.allo.summ.df$file.name.r==combo,]$fdm.min.less.10th.perc.sims[
    match(shdr.era.east.outlier.df[shdr.era.east.outlier.df$is.allo.hdr=="yes",]$shdr.para.east.id, era.fd.df.shdr.allo.summ.df[era.fd.df.shdr.allo.summ.df$file.name.r==combo,]$shdr.para.east.id )]
  shdr.era.east.outlier.df[shdr.era.east.outlier.df$is.allo.hdr=="yes",][paste("fdm.min.less.5th.perc.sims", combo, sep="_")] <- era.fd.df.shdr.allo.summ.df[era.fd.df.shdr.allo.summ.df$file.name.r==combo,]$fdm.min.less.5th.perc.sims[
    match(shdr.era.east.outlier.df[shdr.era.east.outlier.df$is.allo.hdr=="yes",]$shdr.para.east.id, era.fd.df.shdr.allo.summ.df[era.fd.df.shdr.allo.summ.df$file.name.r==combo,]$shdr.para.east.id )]
  
  # west
  shdr.era.west.outlier.df[paste("fdm.max.more.90th.perc.sims", combo, sep="_")] <- NA; shdr.era.west.outlier.df[paste("fdm.max.more.95th.perc.sims", combo, sep="_")] <- NA
  shdr.era.west.outlier.df[paste("fdm.min.less.10th.perc.sims", combo, sep="_")] <- NA; shdr.era.west.outlier.df[paste("fdm.min.less.5th.perc.sims", combo, sep="_")] <- NA
  shdr.era.west.outlier.df[shdr.era.west.outlier.df$is.allo.hdr=="no",][paste("fdm.max.more.90th.perc.sims", combo, sep="_")] <- era.fd.df.shdr.para.west.summ.df[era.fd.df.shdr.para.west.summ.df$file.name.r==combo,]$fdm.max.more.90th.perc.sims[
    match(shdr.era.west.outlier.df[shdr.era.west.outlier.df$is.allo.hdr=="no",]$shdr.para.west.id, era.fd.df.shdr.para.west.summ.df[era.fd.df.shdr.para.west.summ.df$file.name.r==combo,]$shdr.para.west.id )]
  shdr.era.west.outlier.df[shdr.era.west.outlier.df$is.allo.hdr=="no",][paste("fdm.max.more.95th.perc.sims", combo, sep="_")] <- era.fd.df.shdr.para.west.summ.df[era.fd.df.shdr.para.west.summ.df$file.name.r==combo,]$fdm.max.more.95th.perc.sims[
    match(shdr.era.west.outlier.df[shdr.era.west.outlier.df$is.allo.hdr=="no",]$shdr.para.west.id, era.fd.df.shdr.para.west.summ.df[era.fd.df.shdr.para.west.summ.df$file.name.r==combo,]$shdr.para.west.id )]
  shdr.era.west.outlier.df[shdr.era.west.outlier.df$is.allo.hdr=="no",][paste("fdm.min.less.10th.perc.sims", combo, sep="_")] <- era.fd.df.shdr.para.west.summ.df[era.fd.df.shdr.para.west.summ.df$file.name.r==combo,]$fdm.min.less.10th.perc.sims[
    match(shdr.era.west.outlier.df[shdr.era.west.outlier.df$is.allo.hdr=="no",]$shdr.para.west.id, era.fd.df.shdr.para.west.summ.df[era.fd.df.shdr.para.west.summ.df$file.name.r==combo,]$shdr.para.west.id )]
  shdr.era.west.outlier.df[shdr.era.west.outlier.df$is.allo.hdr=="no",][paste("fdm.min.less.5th.perc.sims", combo, sep="_")] <- era.fd.df.shdr.para.west.summ.df[era.fd.df.shdr.para.west.summ.df$file.name.r==combo,]$fdm.min.less.5th.perc.sims[
    match(shdr.era.west.outlier.df[shdr.era.west.outlier.df$is.allo.hdr=="no",]$shdr.para.west.id, era.fd.df.shdr.para.west.summ.df[era.fd.df.shdr.para.west.summ.df$file.name.r==combo,]$shdr.para.west.id )]
  
  # allo west
  shdr.era.west.outlier.df[shdr.era.west.outlier.df$is.allo.hdr=="yes",][paste("fdm.max.more.90th.perc.sims", combo, sep="_")] <- era.fd.df.shdr.allo.summ.df[era.fd.df.shdr.allo.summ.df$file.name.r==combo,]$fdm.max.more.90th.perc.sims[
    match(shdr.era.west.outlier.df[shdr.era.west.outlier.df$is.allo.hdr=="yes",]$shdr.para.west.id, era.fd.df.shdr.allo.summ.df[era.fd.df.shdr.allo.summ.df$file.name.r==combo,]$shdr.para.west.id )]
  shdr.era.west.outlier.df[shdr.era.west.outlier.df$is.allo.hdr=="yes",][paste("fdm.max.more.95th.perc.sims", combo, sep="_")] <- era.fd.df.shdr.allo.summ.df[era.fd.df.shdr.allo.summ.df$file.name.r==combo,]$fdm.max.more.95th.perc.sims[
    match(shdr.era.west.outlier.df[shdr.era.west.outlier.df$is.allo.hdr=="yes",]$shdr.para.west.id, era.fd.df.shdr.allo.summ.df[era.fd.df.shdr.allo.summ.df$file.name.r==combo,]$shdr.para.west.id )]
  shdr.era.west.outlier.df[shdr.era.west.outlier.df$is.allo.hdr=="yes",][paste("fdm.min.less.10th.perc.sims", combo, sep="_")] <- era.fd.df.shdr.allo.summ.df[era.fd.df.shdr.allo.summ.df$file.name.r==combo,]$fdm.min.less.10th.perc.sims[
    match(shdr.era.west.outlier.df[shdr.era.west.outlier.df$is.allo.hdr=="yes",]$shdr.para.west.id, era.fd.df.shdr.allo.summ.df[era.fd.df.shdr.allo.summ.df$file.name.r==combo,]$shdr.para.west.id )]
  shdr.era.west.outlier.df[shdr.era.west.outlier.df$is.allo.hdr=="yes",][paste("fdm.min.less.5th.perc.sims", combo, sep="_")] <- era.fd.df.shdr.allo.summ.df[era.fd.df.shdr.allo.summ.df$file.name.r==combo,]$fdm.min.less.5th.perc.sims[
    match(shdr.era.west.outlier.df[shdr.era.west.outlier.df$is.allo.hdr=="yes",]$shdr.para.west.id, era.fd.df.shdr.allo.summ.df[era.fd.df.shdr.allo.summ.df$file.name.r==combo,]$shdr.para.west.id )]
  
}

names(shdr.era.west.outlier.df)
test <- shdr.era.west.outlier.df[,c(1:10,35,70:102)]; head(test)

# save
write.csv(shdr.era.east.outlier.df , "local/data/shdr.summ/shdr.era.east.outlier.df.csv", row.names = F)
write.csv(shdr.era.west.outlier.df, "local/data/shdr.summ/shdr.era.west.outlier.df.csv", row.names = F)



#######################################################################################################  3. MELPOMENE fd ############################
ref.scaff.mel <- read.table("local/data/ref/Hmel2.5.scaffolds.fa.fai", row.names = NULL)

ref.scaff.mel<-ref.scaff.mel[,1:2]
names(ref.scaff.mel)<-c("scaffold","length")
ref.scaff.mel$add<-c(0,cumsum(ref.scaff.mel$length)[-length(ref.scaff.mel$length)])
ref.scaff.mel$chr<-substring(ref.scaff.mel$scaffold,first = 6,last=7);head(ref.scaff.mel)

scafEnds.mel <- cumsum(ref.scaff.mel[,2])
mel.genome_size <- tail(scafEnds.mel,1)

# within species allopatric from highlands
# include
fd_melCvlowE_melChighE_melChighW_hecale <-getFd.mel("fd_melCvlowE_melChighE_melChighW_hecale")
fd_melEvlowE_melEhighE_melEhighW_hecale <-getFd.mel("fd_melEvlowE_melEhighE_melEhighW_hecale")

fd_melCvlowW_melChighW_melChighE_hecale <-getFd.mel("fd_melCvlowW_melChighW_melChighE_hecale")
fd_melElowW_melEhighW_melEhighE_hecale <-getFd.mel("fd_melElowW_melEhighW_melEhighE_hecale")


# timareta introgression from highland tim, into era high/low ()
fd_eraCvlowW_eraChighW_timEhighE_eleHW<-getFd.mel("fd_eraCvlowW_eraChighW_timEhighE_eleHW")
fd_eraCvlowE_eraChighE_timEhighE_eleHW<-getFd.mel("fd_eraCvlowE_eraChighE_timEhighE_eleHW")
fd_eraElowW_eraEhighW_timEhighE_eleHW<-getFd.mel("fd_eraElowW_eraEhighW_timEhighE_eleHW")
fd_eraEvlowE_eraEhighE_timEhighE_eleHW<-getFd.mel("fd_eraEvlowE_eraEhighE_timEhighE_eleHW")

# cydno introgression
fd_eraCvlowE_eraChighE_cydEhighW_eleHW<-getFd.mel("fd_eraCvlowE_eraChighE_cydEhighW_eleHW")
fd_eraCvlowW_eraChighW_cydEhighW_eleHW<-getFd.mel("fd_eraCvlowW_eraChighW_cydEhighW_eleHW")
fd_eraElowW_eraEhighW_cydEhighW_eleHW<-getFd.mel("fd_eraElowW_eraEhighW_cydEhighW_eleHW")
fd_eraEvlowE_eraEhighE_cydEhighW_eleHW<-getFd.mel("fd_eraEvlowE_eraEhighE_cydEhighW_eleHW")


# #### joana Plot fd ####
# 
# # Plot melpomene allele sharing across sides of the Andes
# 
# pdf(file="melpomene_acrossAndes.pdf",height=7.5,width=15)
#   par(mfrow=c(8,1),mar=c(0,4,1.5,1),oma=c(2,1,1,1))
#   plotOutl(fd_melCvlowW_melChighW_melChighE_eleHW,title="fd_melCvlowW_melChighW_melChighE_eleHW",maxy = 0.5,horiz=T)
#   legend("topleft",pch=19,col=c("grey","darkorange","cornflowerblue","black"),
#          legend=c("normal","East outlier","West outlier","across outlier"),
#          bty="n",ncol=2)
#   plotOutl(fd_melCvlowW_melChighW_melCvlowE_eleHW,title="fd_melCvlowW_melChighW_melCvlowE_eleHW",maxy = 0.5,horiz=T)
#   plotOutl(fd_melElowW_melEhighW_melEhighE_eleHW,title="fd_melElowW_melEhighW_melEhighE_eleHW",maxy = 0.5,horiz=T)
#   plotOutl(fd_melElowW_melEhighW_melEvlowE_eleHW,title="fd_melElowW_melEhighW_melEvlowE_eleHW",maxy = 0.5,horiz=T)
#   plotOutl(fd_melCvlowE_melChighE_melChighW_eleHW,title="fd_melCvlowE_melChighE_melChighW_eleHW",maxy = 0.5,horiz=T)
#   plotOutl(fd_melCvlowE_melChighE_melCvlowW_eleHW,title="fd_melCvlowE_melChighE_melCvlowW_eleHW",maxy = 0.5,horiz=T)
#   plotOutl(fd_melEvlowE_melEhighE_melEhighW_eleHW,title="fd_melEvlowE_melEhighE_melEhighW_eleHW",maxy = 0.5,horiz=T)
#   plotOutl(fd_melEvlowE_melEhighE_melElowW_eleHW,title="fd_melEvlowE_melEhighE_melElowW_eleHW",maxy = 0.5,horiz=T)
#   axis(1,chrom$add[!duplicated(chrom$chr)][-22]+(diff(chrom$add[!duplicated(chrom$chr)])/2),
#        labels = paste0("chr",chrom$chr[!duplicated(chrom$chr)][-22]),tick = F)
# dev.off()
# 
# 
# # Plot melpomene allele sharing on the same side of the Andes
# pdf(file="melpomene_sameSideAndes.pdf",height=7.5,width=15)
#   par(mfrow=c(8,1),mar=c(0,4,2,1),oma=c(2,1,1,1))
#   plotOutl(fd_melCvlowW_melChighW_melEhighW_eleHW,title="fd_melCvlowW_melChighW_melEhighW_eleHW",maxy = 1,horiz=T)
#   legend("topleft",pch=19,col=c("grey","darkorange","cornflowerblue","black"),legend=c("normal","East outlier","West outlier","across outlier"),horiz = F,ncol=2,bty="n")
#   plotOutl(fd_melCvlowW_melChighW_melElowW_eleHW,title="fd_melCvlowW_melChighW_melElowW_eleHW",maxy = 1,horiz=T)
#   plotOutl(fd_melElowW_melEhighW_melChighW_eleHW,title="fd_melElowW_melEhighW_melChighW_eleHW",maxy = 1,horiz=T)
#   plotOutl(fd_melElowW_melEhighW_melCvlowW_eleHW,title="fd_melElowW_melEhighW_melCvlowW_eleHW",maxy = 1,horiz=T)
#   plotOutl(fd_melCvlowE_melChighE_melEhighE_eleHW,title="fd_melCvlowE_melChighE_melEhighE_eleHW",maxy = 1,horiz=T)
#   plotOutl(fd_melCvlowE_melChighE_melEvlowE_eleHW,title="fd_melCvlowE_melChighE_melEvlowE_eleHW",maxy = 1,horiz=T)
#   plotOutl(fd_melEvlowE_melEhighE_melChighE_eleHW,title="fd_melEvlowE_melEhighE_melChighE_eleHW",maxy = 1,horiz=T)
#   plotOutl(fd_melEvlowE_melEhighE_melCvlowE_eleHW,title="fd_melEvlowE_melEhighE_melCvlowE_eleHW",maxy = 1,horiz=T)
#   axis(1,chrom$add[!duplicated(chrom$chr)][-22]+(diff(chrom$add[!duplicated(chrom$chr)])/2),labels = paste0("chr",chrom$chr[!duplicated(chrom$chr)][-22]),tick = F)
# dev.off()
# 
# # timareta introgression
# pdf(file="melpomene_timareta.pdf",height=7.5,width=15)
# 
# par(mfrow=c(4,1),mar=c(0,4,2,1),oma=c(2,1,1,1))
# 
# plotOutl(fd_melCvlowW_melChighW_timEhighE_eleHW,title="fd_melCvlowW_melChighW_timEhighE_eleHW",maxy = 0.7)
# legend("topleft",pch=19,col=c("grey","darkorange","cornflowerblue","black"),legend=c("normal","East outlier","West outlier","across outlier"),bty="n")
# plotOutl(fd_melCvlowE_melChighE_timEhighE_eleHW,title="fd_melCvlowE_melChighE_timEhighE_eleHW",maxy = 0.7)
# plotOutl(fd_melElowW_melEhighW_timEhighE_eleHW,title="fd_melElowW_melEhighW_timEhighE_eleHW",maxy = 0.7)
# plotOutl(fd_melEvlowE_melEhighE_timEhighE_eleHW,title="fd_melEvlowE_melEhighE_timEhighE_eleHW",maxy = 0.7)
# axis(1,chrom$add[!duplicated(chrom$chr)][-22]+(diff(chrom$add[!duplicated(chrom$chr)])/2),labels = paste0("chr",chrom$chr[!duplicated(chrom$chr)][-22]),tick = F)
# 
# dev.off()
# 
# 
# # cydno introgression
# pdf(file="melpomene_cydno.pdf",height=7.5,width=15)
# 
# par(mfrow=c(4,1),mar=c(0,4,2,1),oma=c(2,1,1,1))
# 
# plotOutl(fd_melCvlowE_melChighE_cydEhighW_eleHW,maxy = 0.7)
# legend("topleft",pch=19,col=c("grey","darkorange","cornflowerblue","black"),legend=c("normal","East outlier","West outlier","across outlier"),bty="n")
# plotOutl(fd_melCvlowW_melChighW_cydEhighW_eleHW,maxy = 0.7)
# plotOutl(fd_melElowW_melEhighW_cydEhighW_eleHW,maxy = 0.7)
# plotOutl(fd_melEvlowE_melEhighE_cydEhighW_eleHW,maxy = 0.7)
# axis(1,chrom$add[!duplicated(chrom$chr)][-22]+(diff(chrom$add[!duplicated(chrom$chr)])/2),labels = paste0("chr",chrom$chr[!duplicated(chrom$chr)][-22]),tick = F)
# 
# dev.off()
# 
# 
# 
# 
# 
# # Difference in fd between P3=high - P3=low
# diffCE<-fd_melCvlowE_melChighE_melEhighE_eleHW
# diffCE$fd<-diffCE$fd-fd_melCvlowE_melChighE_melEvlowE_eleHW$fd
# diffEE<-fd_melEvlowE_melEhighE_melChighE_eleHW
# diffEE$fd<-diffEE$fd-fd_melEvlowE_melEhighE_melCvlowE_eleHW$fd
# diffCW<-fd_melCvlowW_melChighW_melEhighW_eleHW
# diffCW$fd<-diffCW$fd-fd_melCvlowW_melChighW_melElowW_eleHW$fd
# diffEW<-fd_melElowW_melEhighW_melChighW_eleHW[paste0(fd_melElowW_melEhighW_melChighW_eleHW$scaffold,fd_melElowW_melEhighW_melChighW_eleHW$start)%in%paste0(fd_melElowW_melEhighW_melCvlowW_eleHW$scaffold,fd_melElowW_melEhighW_melCvlowW_eleHW$start),]
# diffEW$fd<-diffEW$fd-fd_melElowW_melEhighW_melCvlowW_eleHW$fd
# 
# 
# pdf(file="melpomene_P3high-low.pdf",height=7.5,width=15)
# 
# par(mfrow=c(4,1),mar=c(0,4,2,1),oma=c(2,1,1,1))
# 
# plotOutl(diffCE,title="fd_melCvlowE_melChighE_melEhighE-ELE_eleHW",maxy = 0.6)
# legend("topleft",pch=19,col=c("grey","darkorange","cornflowerblue","black"),legend=c("normal","East outlier","West outlier","across outlier"),horiz = F,ncol=2,bty="n")
# plotOutl(diffEE,title="melEvlowE_melEhighE_melChighE-CLE_eleHW",maxy = 0.6)
# plotOutl(diffCW,title="melCvlowW_melChighW_melEhighW-ELW_eleHW",maxy = 0.6)
# plotOutl(diffEW,title="melElowW_melEhighW_melChighW-CLW_eleHW",maxy = 0.6)
# axis(1,chrom$add[!duplicated(chrom$chr)][-22]+(diff(chrom$add[!duplicated(chrom$chr)])/2),labels = paste0("chr",chrom$chr[!duplicated(chrom$chr)][-22]),tick = F)
# 
# dev.off()
# 
# 
# 

#### Read in outlier regions and chromosome information ####

# Get melpomene outliers shared by either or both East and West:
mel.shdr.para.east<-read.csv("local/data/sharing/mel.shdr.para.east.csv",header=T); head(mel.shdr.para.east)
mel.shdr.para.west<-read.csv("local/data/sharing/mel.shdr.para.west.csv",header=T)
mel.shdr.allo.all <-read.csv("local/data/sharing/mel.shdr.allo.all.csv",header=T)


# add 
mel.shdr.allo.all.summ <-read.csv("local/data/sharing/mel.shdr.allo.summ.csv",header=T)


# unlike when working with Tajima'sD etc, where we were comparing selection statistics obtained from pops in different sides of the Andes
# twisst includes all populations (both sides), so we dont want to divide the analyses between SHDR west / east
# so best we use SHDR ranges of para.east/west (only, no allo) and allopatric ranges
mel.outliersEast<-mel.shdr.para.east[mel.shdr.para.east$is.allopatric=="no",]
# use mel.shdr.allo.all because it'll have max start/ends (all overlap between east/west)
mel.outliersAcross<-mel.shdr.allo.all
mel.outliersWest<-mel.shdr.para.west[mel.shdr.para.west$is.allopatric=="no",]

# new version Oct21- use combined para/allo shdrs
mel.outliersEast<-mel.shdr.para.east
mel.outliersWest<-mel.shdr.para.west


# Read in chromosome information
ref.scaff.mel <- read.table("local/data/ref/Hmel2.5.scaffolds.fa.fai", row.names = NULL)
ref.scaff.mel<-ref.scaff.mel[,1:2]
names(ref.scaff.mel)<-c("scaffold","length")
ref.scaff.mel$add<-c(0,cumsum(ref.scaff.mel$length)[-length(ref.scaff.mel$length)])
ref.scaff.mel$chr<-substring(ref.scaff.mel$scaffold,first = 6,last=7);head(ref.scaff.mel)

scafEnds.mel <- cumsum(ref.scaff.mel[,2])
mel.genome_size <- tail(scafEnds.mel,1); mel.genome_size

######################################################## 3. mel fd/fdm from other species (parapatric,  n=4 combos (2 cydno, 2 timareta)) ############################

############################ 3.1. mel fd prep data ############################

#### Read in fd datasets ####
mel.comparisons.info <- read.csv("local/data/joana.fd/mel.comparisons.csv")
mel.comparisons.info$type.comparison
subset(mel.comparisons.info, type.comparison=="other.species.parapatric" | type.comparison=="mel.allopatric.from.hig")$file.name.r

# within species allopatric from highlands
# include
fd_melCvlowE_melChighE_melChighW_hecale <-getFd.mel("fd_melCvlowE_melChighE_melChighW_hecale")
fd_melEvlowE_melEhighE_melEhighW_hecale <-getFd.mel("fd_melEvlowE_melEhighE_melEhighW_hecale")

fd_melCvlowW_melChighW_melChighE_hecale <-getFd.mel("fd_melCvlowW_melChighW_melChighE_hecale")
fd_melElowW_melEhighW_melEhighE_hecale <-getFd.mel("fd_melElowW_melEhighW_melEhighE_hecale")

# other.species.parapatric
# cydno west into highland col/ecu W
fd_melCvlowW_melChighW_cydEhighW_hecale <-getFd.mel("fd_melCvlowW_melChighW_cydEhighW_hecale")
fd_melElowW_melEhighW_cydEhighW_hecale <-getFd.mel("fd_melElowW_melEhighW_cydEhighW_hecale")

# timareta west into highland col/ecu W
fd_melCvlowE_melChighE_timEhighE_hecale <-getFd.mel("fd_melCvlowE_melChighE_timEhighE_hecale")
fd_melEvlowE_melEhighE_timEhighE_hecale<-getFd.mel("fd_melEvlowE_melEhighE_timEhighE_hecale")

# other species allopatric
# fd_melClowW_melChighW_timEhighE_hecale
# fd_melCvlowE_melChighE_cydEhighW_hecale
# fd_melCvlowW_melChighW_timEhighE_hecale
# fd_melElowW_melEhighW_timEhighE_hecale
# fd_melEvlowE_melEhighE_cydEhighW_hecale

# list files in environment that match, then creaste list
mel.fd.list.files <-grep("fd_mel",names(.GlobalEnv),value=TRUE); mel.fd.list.files
mel.fd.list<-do.call("list",mget(mel.fd.list.files)); str(mel.fd.list)
for (i in 1:length(mel.fd.list)) {
  mel.fd.list[[i]]$file.name.r <- c(mel.fd.list.files[i]) }

# remove all individual files
rm(list = ls()[grep("fd_mel*", ls())])

# bind and add info
mel.fd.df <- rbindlist(mel.fd.list)
names(mel.comparisons.info)
mel.fd.df<-merge(mel.fd.df, mel.comparisons.info[c(2,19:23)], by="file.name.r"); head(mel.fd.df)

############################ 3.3. summarise each combo per hdr ############################
# start with parapatric sharing only, n=8 combos (2 tele, 2 him, 4 clys)

mel.other.species.parapatric.list <- unique(subset(mel.fd.df, type.comparison=="other.species.parapatric"| type.comparison=="mel.allopatric.from.hig")$file.name.r); mel.other.species.parapatric.list

# max fd per hdr, max/min fdM per hdr
mel.fd.df.shdr.para.east.summ.list<- list(); mel.fd.df.shdr.para.west.summ.list<- list(); mel.fd.df.shdr.allo.summ.list<- list()

for (i in 1:length(mel.other.species.parapatric.list)) {
  combo <- mel.other.species.parapatric.list[i]
  mel.fd.df.shdr.para.east.summ.list[[i]] <- summarise(group_by(subset(mel.fd.df, file.name.r==combo &shdr.para.east.id!=""),file.name.r, p1.newname,p2.newname,p3.newname, shdr.para.east.id),
                                                       no.windows=n(),fd.max=max(fd),fdm.max=max(fdM),fdm.min=min(fdM))
  mel.fd.df.shdr.para.west.summ.list[[i]] <- summarise(group_by(subset(mel.fd.df, file.name.r==combo &shdr.para.west.id!=""),file.name.r, p1.newname,p2.newname,p3.newname, shdr.para.west.id),
                                                       no.windows=n(),fd.max=max(fd),fdm.max=max(fdM),fdm.min=min(fdM))
}

# bind 
mel.fd.df.shdr.para.east.summ.df <- rbindlist(mel.fd.df.shdr.para.east.summ.list)
mel.fd.df.shdr.para.west.summ.df <- rbindlist(mel.fd.df.shdr.para.west.summ.list)


# each shdr.para.east/west.summ will have 8x3 new columns, for each combo max fd, max/min fdm is 10th-90th perc outlier? yes/no. 
# but! some of these combos more relevant than others:
# if it's an eastern combo it'll be more relevant for shdr.para.east than for shdr.para.west (compare!)

############################ 3.3. sims east/west para and allo SHDR ###################
#### prep sims datasets ######
# extract the same number of shdr para east and of the same width, random locations get tajima min, tajima.5thperc.no.windows, tajima.5thperc.mean.taj
# first create intervals with real shdr. with Intervals we must work with BP.wg (not bp per scaff)
mel.shdr.para.east.int <- reduce(Intervals(as.matrix(mel.outliersEast[,c("start","end")]), closed=c(TRUE,TRUE), type="R"))
mel.shdr.para.west.int <- reduce(Intervals(as.matrix(mel.outliersWest[,c("start","end")]), closed=c(TRUE,TRUE), type="R"))

# random intervals for each data set
sim1k.mel.para.east <- lapply(1:1000, function(x){rand_non_overlapping_intervals(mel.genome_size,size(mel.shdr.para.east.int))})
sim1k.mel.para.west <- lapply(1:1000, function(x){rand_non_overlapping_intervals(mel.genome_size,size(mel.shdr.para.west.int))})


# run simulations on each ecuador/ colombia dataset separately
# careful because we are checking for overlaps between two sets of ranges (twisst is in 20kb windows)
sim.combo.fd.tmp.shdr.para.east.summ.list<-list(); sim.combo.fd.tmp.shdr.para.west.summ.list <-list()
sim.combo.fd.tmp.shdr.para.east.summ.list.list<-list(); sim.combo.fd.tmp.shdr.para.west.summ.list.list<-list()

for (i in 1:1000) {
  # prep random intervals for action
  sim1k.mel.para.east.df <- as.data.frame(sim1k.mel.para.east[i]); sim1k.mel.para.east.df$ran.para.east.hdr.name <- paste("ran", i, ".hdr.", seq(1:nrow(sim1k.mel.para.east.df )), sep = "")
  sim1k.mel.para.west.df <- as.data.frame(sim1k.mel.para.west[i]); sim1k.mel.para.west.df$ran.para.west.hdr.name <- paste("ran", i, ".hdr.", seq(1:nrow(sim1k.mel.para.west.df )), sep = "")

  # prepare all combo tmps, run for east/west/allo shdrs
  for (j in 1:length(mel.other.species.parapatric.list)) {
    combo <- mel.other.species.parapatric.list[j]
    combo.fd.tmp <- subset(mel.fd.df, file.name.r==combo)[,c( "file.name.r", "start.WG", "end.WG", "fd","fdM")]
    
    #### east shdr ####
    # add ran hdr names to combo.tmp
    sim1k.mel.para.east.dt <- as.data.table(sim1k.mel.para.east.df[,1:2])
    names(sim1k.mel.para.east.dt)[1:2]<-c("start.WG","end.WG") # to match fd tmp names
    
    setkey(sim1k.mel.para.east.dt)
    # find overlaps
    sim.col.oveast<-foverlaps(x = combo.fd.tmp, y = sim1k.mel.para.east.dt, which=T, by.x = names(sim1k.mel.para.east.dt), type="any",nomatch=0L)
    # use keys to select corresponding shdr names
    combo.fd.tmp$ran.para.east.hdr.name <- NA; combo.fd.tmp$ran.para.east.hdr.name[sim.col.oveast$xid]<- sim1k.mel.para.east.df[sim.col.oveast$yid,]$ran.para.east.hdr.name
    
    #### east shdr ####
    # add ran hdr names to combo.tmp
    sim1k.mel.para.west.dt <- as.data.table(sim1k.mel.para.west.df[,1:2])
    names(sim1k.mel.para.west.dt)[1:2]<-c("start.WG","end.WG") # to match fd tmp names
    
    setkey(sim1k.mel.para.west.dt)
    # find overlaps
    sim.col.oveast<-foverlaps(x = combo.fd.tmp, y = sim1k.mel.para.west.dt, which=T, by.x = names(sim1k.mel.para.west.dt), type="any",nomatch=0L)
    # use keys to select corresponding shdr names
    combo.fd.tmp$ran.para.west.hdr.name <- NA; combo.fd.tmp$ran.para.west.hdr.name[sim.col.oveast$xid]<- sim1k.mel.para.west.df[sim.col.oveast$yid,]$ran.para.west.hdr.name
    
    #### summarise combo.fd.tmp by shdr type ####
    # get rid of NA ran hdr (background)
    sim.combo.fd.tmp.shdr.para.east.summ.list[[j]] <- summarise(group_by(subset(combo.fd.tmp, file.name.r==combo& ran.para.east.hdr.name!=""),file.name.r, ran.para.east.hdr.name),
                                                                no.windows=n(),fd.max=max(fd),fdm.max=max(fdM),fdm.min=min(fdM))
    sim.combo.fd.tmp.shdr.para.west.summ.list[[j]] <- summarise(group_by(subset(combo.fd.tmp, file.name.r==combo& ran.para.west.hdr.name!=""),file.name.r, ran.para.west.hdr.name ),
                                                                no.windows=n(),fd.max=max(fd),fdm.max=max(fdM),fdm.min=min(fdM))
    
  }
  # store combo summaries
  sim.combo.fd.tmp.shdr.para.east.summ.list.list[[i]] <- rbindlist(sim.combo.fd.tmp.shdr.para.east.summ.list)
  sim.combo.fd.tmp.shdr.para.west.summ.list.list[[i]] <- rbindlist(sim.combo.fd.tmp.shdr.para.west.summ.list)

}

sim.combo.fd.tmp.shdr.para.east.summ.list.list.df <- rbindlist(sim.combo.fd.tmp.shdr.para.east.summ.list.list); head(sim.combo.fd.tmp.shdr.para.east.summ.list.list.df)
sim.combo.fd.tmp.shdr.para.west.summ.list.list.df <- rbindlist(sim.combo.fd.tmp.shdr.para.west.summ.list.list); head(sim.combo.fd.tmp.shdr.para.west.summ.list.list.df)

# get rid of shdrs that are smaller than 100kb (only have 1 50kb window)
sim.combo.fd.tmp.shdr.para.east.summ.list.list.df.sub <- subset(sim.combo.fd.tmp.shdr.para.east.summ.list.list.df, no.windows>1)
sim.combo.fd.tmp.shdr.para.west.summ.list.list.df.sub <- subset(sim.combo.fd.tmp.shdr.para.west.summ.list.list.df, no.windows>1)

# store thresholds from simulations
mel.other.species.east.para.sims.threhsolds <- data.frame(combo.para.allo.type=NA, combo=NA, min.fdm.mean=NA, max.fdm.mean=NA,  min.fdm.se=NA, max.fdm.se=NA, min.fdm.5th=NA, min.fdm.10th=NA, max.fdm.90th=NA, max.fdm.95th=NA ); mel.other.species.para.sims.threhsolds
for (i in 1:length(mel.other.species.parapatric.list)) {
  combo <- mel.other.species.parapatric.list[i]
  mel.other.species.east.para.sims.threhsolds[i,]$combo.para.allo.type <- "shdr.para.east"
  mel.other.species.east.para.sims.threhsolds[i,]$combo <- combo
  mel.other.species.east.para.sims.threhsolds[i,]$min.fdm.mean <- mean(subset(sim.combo.fd.tmp.shdr.para.east.summ.list.list.df.sub, file.name.r==combo)$fdm.min) 
  mel.other.species.east.para.sims.threhsolds[i,]$max.fdm.mean <- mean(subset(sim.combo.fd.tmp.shdr.para.east.summ.list.list.df.sub, file.name.r==combo)$fdm.max) 
  mel.other.species.east.para.sims.threhsolds[i,]$min.fdm.se <- se(subset(sim.combo.fd.tmp.shdr.para.east.summ.list.list.df.sub, file.name.r==combo)$fdm.min) 
  mel.other.species.east.para.sims.threhsolds[i,]$max.fdm.se <- se(subset(sim.combo.fd.tmp.shdr.para.east.summ.list.list.df.sub, file.name.r==combo)$fdm.max) 
  mel.other.species.east.para.sims.threhsolds[i,]$min.fdm.5th <- quantile(subset(sim.combo.fd.tmp.shdr.para.east.summ.list.list.df.sub, file.name.r==combo)$fdm.min, 0.05) 
  mel.other.species.east.para.sims.threhsolds[i,]$min.fdm.10th <- quantile(subset(sim.combo.fd.tmp.shdr.para.east.summ.list.list.df.sub, file.name.r==combo)$fdm.min, 0.1) 
  mel.other.species.east.para.sims.threhsolds[i,]$max.fdm.90th <- quantile(subset(sim.combo.fd.tmp.shdr.para.east.summ.list.list.df.sub, file.name.r==combo)$fdm.max, 0.9) 
  mel.other.species.east.para.sims.threhsolds[i,]$max.fdm.95th <- quantile(subset(sim.combo.fd.tmp.shdr.para.east.summ.list.list.df.sub, file.name.r==combo)$fdm.max, 0.95) 
}; mel.other.species.east.para.sims.threhsolds

mel.other.species.west.para.sims.threhsolds <- data.frame(combo.para.allo.type=NA, combo=NA, min.fdm.mean=NA, max.fdm.mean=NA,  min.fdm.se=NA, max.fdm.se=NA, min.fdm.5th=NA, min.fdm.10th=NA, max.fdm.90th=NA, max.fdm.95th=NA ); mel.other.species.para.sims.threhsolds
for (i in 1:length(mel.other.species.parapatric.list)) {
  combo <- mel.other.species.parapatric.list[i]
  mel.other.species.west.para.sims.threhsolds[i,]$combo.para.allo.type <- "shdr.para.west"
  mel.other.species.west.para.sims.threhsolds[i,]$combo <- combo
  mel.other.species.west.para.sims.threhsolds[i,]$min.fdm.mean <- mean(subset(sim.combo.fd.tmp.shdr.para.west.summ.list.list.df.sub, file.name.r==combo)$fdm.min) 
  mel.other.species.west.para.sims.threhsolds[i,]$max.fdm.mean <- mean(subset(sim.combo.fd.tmp.shdr.para.west.summ.list.list.df.sub, file.name.r==combo)$fdm.max) 
  mel.other.species.west.para.sims.threhsolds[i,]$min.fdm.se <- se(subset(sim.combo.fd.tmp.shdr.para.west.summ.list.list.df.sub, file.name.r==combo)$fdm.min) 
  mel.other.species.west.para.sims.threhsolds[i,]$max.fdm.se <- se(subset(sim.combo.fd.tmp.shdr.para.west.summ.list.list.df.sub, file.name.r==combo)$fdm.max) 
  mel.other.species.west.para.sims.threhsolds[i,]$min.fdm.5th <- quantile(subset(sim.combo.fd.tmp.shdr.para.west.summ.list.list.df.sub, file.name.r==combo)$fdm.min, 0.05) 
  mel.other.species.west.para.sims.threhsolds[i,]$min.fdm.10th <- quantile(subset(sim.combo.fd.tmp.shdr.para.west.summ.list.list.df.sub, file.name.r==combo)$fdm.min, 0.1) 
  mel.other.species.west.para.sims.threhsolds[i,]$max.fdm.90th <- quantile(subset(sim.combo.fd.tmp.shdr.para.west.summ.list.list.df.sub, file.name.r==combo)$fdm.max, 0.9) 
  mel.other.species.west.para.sims.threhsolds[i,]$max.fdm.95th <- quantile(subset(sim.combo.fd.tmp.shdr.para.west.summ.list.list.df.sub, file.name.r==combo)$fdm.max, 0.95) 
}; mel.other.species.west.para.sims.threhsolds


mel.other.species.sims.threhsolds <- rbind(mel.other.species.east.para.sims.threhsolds, mel.other.species.west.para.sims.threhsolds)
write.csv(mel.other.species.sims.threhsolds, "local/data/joana.fd/mel.within.and.other.species.1ksims.fdm.threhsolds.means.csv", row.names = T)

############################ 3.4. plot each combo, obs vs simulation thresholds ############################

###### plot fdm sims  east para / west /allo SHDR ######
head(mel.fd.df.shdr.para.east.summ.df); head(sim.combo.fd.tmp.shdr.para.east.summ.list.list.df)

fd.para.east.sims.plot <- list(); for (i in 1:length(mel.other.species.parapatric.list)) {
  combo <- mel.other.species.parapatric.list[i]
  fd.para.east.sims.plot[[i]] <-  ggplot(subset(mel.fd.df, file.name.r==combo ), aes(x=fdM)) + 
    labs(colour="Altitude")+
    geom_vline(xintercept = 0, colour="lightgrey", lty="dashed")+
    geom_vline(xintercept = subset(mel.fd.df.shdr.para.east.summ.df, file.name.r==combo & no.windows>1)$fdm.max, colour="#049E73", lty="solid", alpha=.9, size=0.5)+
    geom_vline(xintercept = subset(mel.fd.df.shdr.para.east.summ.df, file.name.r==combo & no.windows>1 & shdr.para.east.id %in% subset(mel.outliersEast, overlaps.with.inversion=="yes")$shdr.para.east.id )$fdm.max, 
               colour="yellow", lty="solid", alpha=1, size=0.7)+
    #geom_vline(xintercept = mean(subset(sim.combo.fd.tmp.shdr.para.east.summ.list.list.df.sub, file.name.r==combo)$fdm.max), colour="black", lty="solid", alpha=.9, size=1.5)+
    geom_vline(xintercept = quantile(subset(sim.combo.fd.tmp.shdr.para.east.summ.list.list.df.sub, file.name.r==combo)$fdm.max, 0.9), colour="black", lty="dashed", alpha=.9, size=1.5)+
    geom_vline(xintercept = quantile(subset(sim.combo.fd.tmp.shdr.para.east.summ.list.list.df.sub, file.name.r==combo)$fdm.max, 0.1), colour="black", lty="dashed", alpha=.9, size=1.5)+
    
    geom_vline(xintercept = subset(mel.fd.df.shdr.para.east.summ.df, file.name.r==combo & no.windows>1)$fdm.min, colour="grey", lty="solid", alpha=.9, size=0.5)+
    geom_vline(xintercept = subset(mel.fd.df.shdr.para.east.summ.df, file.name.r==combo & no.windows>1 & shdr.para.east.id %in% subset(mel.outliersEast, overlaps.with.inversion=="yes")$shdr.para.east.id )$fdm.min, 
               colour="yellow", lty="solid", alpha=1, size=0.7)+
    #geom_vline(xintercept = mean(subset(sim.combo.fd.tmp.shdr.para.east.summ.list.list.df.sub, file.name.r==combo)$fdm.min), colour="black", lty="solid", alpha=.9, size=1.5)+
    geom_vline(xintercept = quantile(subset(sim.combo.fd.tmp.shdr.para.east.summ.list.list.df.sub, file.name.r==combo)$fdm.min, 0.9), colour="darkgrey", lty="dashed", alpha=.9, size=1.5)+
    geom_vline(xintercept = quantile(subset(sim.combo.fd.tmp.shdr.para.east.summ.list.list.df.sub, file.name.r==combo)$fdm.min, 0.1), colour="darkgrey", lty="dashed", alpha=.9, size=1.5)+
    geom_rect(aes(xmin=-.35, xmax=.35, ymin=100, ymax=120 ), fill="white")+
    # geom_density(fill="transparent", alpha=0, size=1)+
    xlab(expression(paste("Topo2-Topo3",  italic(" H. melpomene"))))+
    theme_classic()+
    theme(axis.line.x = element_line(color="black", size = 0.5),axis.line.y = element_line(color="black", size = 0.5),axis.title.y = element_text(size=14),
          plot.margin = unit(c(0.5,0.2,0.2,0.2), "lines"),axis.text = element_text(size=10), #axis.text.x = element_blank(),
          axis.text.y = element_text(size=10),axis.title.x = element_blank(), 
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          panel.background = element_rect(fill = "transparent"), plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
          legend.background = element_rect(fill = "transparent"), legend.box.background = element_rect(fill = "transparent"), strip.background.y = element_blank(),
          strip.text = element_text(size=14))+
    scale_y_continuous(limits = c(0,120), expand = c(0,0))+
    scale_x_continuous(limits = c(-.35,.35), expand = c(0,0)) }; 
t <- c(fd.para.east.sims.plot[2],fd.para.east.sims.plot[3],fd.para.east.sims.plot[1], fd.para.east.sims.plot[4])
plot_grid(plotlist = t  , ncol = 2,
          labels = list(mel.other.species.parapatric.list[2],mel.other.species.parapatric.list[3],mel.other.species.parapatric.list[1],mel.other.species.parapatric.list[4]), label_x = -0.3)
ggsave("plots/joana.fd/mel.fdm.other.sp.para.east.90th.fdm.max.min.png", width = 10, height = 9)

fd.para.west.sims.plot <- list(); for (i in 1:length(mel.other.species.parapatric.list)) {
  combo <- mel.other.species.parapatric.list[i]
  fd.para.west.sims.plot[[i]] <-  ggplot(subset(mel.fd.df, file.name.r==combo ), aes(x=fdM)) + 
    labs(colour="Altitude")+
    geom_vline(xintercept = 0, colour="lightgrey", lty="dashed")+
    geom_vline(xintercept = subset(mel.fd.df.shdr.para.west.summ.df, file.name.r==combo & no.windows>1)$fdm.max, colour="#0372B2", lty="solid", alpha=.9, size=0.5)+
    #geom_vline(xintercept = mean(subset(sim.combo.fd.tmp.shdr.para.west.summ.list.list.df.sub, file.name.r==combo)$fdm.max), colour="black", lty="solid", alpha=.9, size=1.5)+
    geom_vline(xintercept = quantile(subset(sim.combo.fd.tmp.shdr.para.west.summ.list.list.df.sub, file.name.r==combo)$fdm.max, 0.9), colour="black", lty="dashed", alpha=.9, size=1.5)+
    geom_vline(xintercept = quantile(subset(sim.combo.fd.tmp.shdr.para.west.summ.list.list.df.sub, file.name.r==combo)$fdm.max, 0.1), colour="black", lty="dashed", alpha=.9, size=1.5)+
    
    geom_vline(xintercept = subset(mel.fd.df.shdr.para.west.summ.df, file.name.r==combo & no.windows>1)$fdm.min, colour="grey", lty="solid", alpha=.9, size=0.5)+
    # geom_vline(xintercept = mean(subset(sim.combo.fd.tmp.shdr.para.west.summ.list.list.df.sub, file.name.r==combo)$fdm.min), colour="black", lty="solid", alpha=.9, size=1.5)+
    geom_vline(xintercept = quantile(subset(sim.combo.fd.tmp.shdr.para.west.summ.list.list.df.sub, file.name.r==combo)$fdm.min, 0.9), colour="darkgrey", lty="dashed", alpha=.9, size=1.5)+
    geom_vline(xintercept = quantile(subset(sim.combo.fd.tmp.shdr.para.west.summ.list.list.df.sub, file.name.r==combo)$fdm.min, 0.1), colour="darkgrey", lty="dashed", alpha=.9, size=1.5)+
    geom_rect(aes(xmin=-.35, xmax=.35, ymin=100, ymax=120 ), fill="white")+
    
    geom_density(fill="transparent", alpha=0, size=1)+
    xlab(expression(paste("Topo2-Topo3",  italic(" H. melpomene"))))+
    theme_classic()+
    theme(axis.line.x = element_line(color="black", size = 0.5),axis.line.y = element_line(color="black", size = 0.5),axis.title.y = element_text(size=14),
          plot.margin = unit(c(0.5,0.2,0.2,0.2), "lines"),axis.text = element_text(size=10), #axis.text.x = element_blank(),
          axis.text.y = element_text(size=10),axis.title.x = element_blank(), 
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          panel.background = element_rect(fill = "transparent"), plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
          legend.background = element_rect(fill = "transparent"), legend.box.background = element_rect(fill = "transparent"), strip.background.y = element_blank(),
          strip.text = element_text(size=14))+
    scale_y_continuous(limits = c(0,120), expand = c(0,0))+
    scale_x_continuous(limits = c(-.35,.35), expand = c(0,0)) } ;
t <- c(fd.para.west.sims.plot[2],fd.para.west.sims.plot[3],fd.para.west.sims.plot[1], fd.para.west.sims.plot[4])
plot_grid(plotlist = t  , ncol = 2,
          labels = list(mel.other.species.parapatric.list[2],mel.other.species.parapatric.list[3],mel.other.species.parapatric.list[1],mel.other.species.parapatric.list[4]), label_x = -0.3)
ggsave("plots/joana.fd/mel.fdm.other.sp.para.west.90th.fdm.max.min.png", width = 10, height = 9)


fd.allo.sims.plot <- list()
for (i in 1:length(mel.other.species.parapatric.list)) {
  combo <- mel.other.species.parapatric.list[i]
  fd.allo.sims.plot[[i]] <-  ggplot(subset(mel.fd.df, file.name.r==combo ), aes(x=fdM)) + 
    labs(colour="Altitude")+
    geom_vline(xintercept = 0, colour="lightgrey", lty="dashed")+
    geom_vline(xintercept = subset(mel.fd.df.shdr.allo.summ.df, file.name.r==combo & no.windows>1)$fdm.max, colour="#D65D00", lty="solid", alpha=.9, size=0.5)+
    geom_vline(xintercept = subset(mel.fd.df.shdr.allo.summ.df, file.name.r==combo & no.windows>1 & shdr.allo.id %in% subset(mel.shdr.allo.all.summ, shdr.para.east %in% subset(mel.shdr.para.east, overlaps.with.inversion=="yes")$shdr.para.east.id )$shdr.allo.id )$fdm.max, 
               colour="yellow", lty="solid", alpha=.9, size=0.7)+
    #geom_vline(xintercept = mean(subset(sim.combo.fd.tmp.shdr.allo.summ.list.list.df.sub, file.name.r==combo)$fdm.max), colour="black", lty="solid", alpha=.9, size=1.5)+
    geom_vline(xintercept = quantile(subset(sim.combo.fd.tmp.shdr.allo.summ.list.list.df.sub, file.name.r==combo)$fdm.max, 0.9), colour="black", lty="dashed", alpha=.9, size=1.5)+
    geom_vline(xintercept = quantile(subset(sim.combo.fd.tmp.shdr.allo.summ.list.list.df.sub, file.name.r==combo)$fdm.max, 0.1), colour="black", lty="dashed", alpha=.9, size=1.5)+
    geom_vline(xintercept = subset(mel.fd.df.shdr.allo.summ.df, file.name.r==combo & no.windows>1)$fdm.min, colour="grey", lty="solid", alpha=.9, size=0.5)+
    geom_vline(xintercept = subset(mel.fd.df.shdr.allo.summ.df, file.name.r==combo & no.windows>1 & shdr.allo.id %in% subset(mel.shdr.allo.all.summ, shdr.para.east %in% subset(mel.shdr.para.east, overlaps.with.inversion=="yes")$shdr.para.east.id )$shdr.allo.id )$fdm.min, 
               colour="yellow", lty="solid", alpha=1, size=0.7)+
    #geom_vline(xintercept = mean(subset(sim.combo.fd.tmp.shdr.allo.summ.list.list.df.sub, file.name.r==combo)$fdm.min), colour="black", lty="solid", alpha=.9, size=1.5)+
    geom_vline(xintercept = quantile(subset(sim.combo.fd.tmp.shdr.allo.summ.list.list.df.sub, file.name.r==combo)$fdm.min, 0.9), colour="darkgrey", lty="dashed", alpha=.9, size=1.5)+
    geom_vline(xintercept = quantile(subset(sim.combo.fd.tmp.shdr.allo.summ.list.list.df.sub, file.name.r==combo)$fdm.min, 0.1), colour="darkgrey", lty="dashed", alpha=.9, size=1.5)+
    
    geom_density(fill="transparent", alpha=0, size=1)+
    xlab(expression(paste("Topo2-Topo3",  italic(" H. melpomene"))))+
    theme_classic()+
    theme(axis.line.x = element_line(color="black", size = 0.5),axis.line.y = element_line(color="black", size = 0.5),axis.title.y = element_text(size=14),
          plot.margin = unit(c(0.5,0.2,0.2,0.2), "lines"),axis.text = element_text(size=10), #axis.text.x = element_blank(),
          axis.text.y = element_text(size=10),axis.title.x = element_blank(), 
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          panel.background = element_rect(fill = "transparent"), plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
          legend.background = element_rect(fill = "transparent"), legend.box.background = element_rect(fill = "transparent"), strip.background.y = element_blank(),
          strip.text = element_text(size=14))+
    scale_y_continuous(limits = c(0,120), expand = c(0,0))+
    scale_x_continuous(limits = c(-.35,.35), expand = c(0,0)) }
t <- c(fd.allo.sims.plot[2],fd.allo.sims.plot[3],fd.allo.sims.plot[1], fd.allo.sims.plot[4])
plot_grid(plotlist = t  , ncol = 2,
          labels = list(mel.other.species.parapatric.list[2],mel.other.species.parapatric.list[3],mel.other.species.parapatric.list[1],mel.other.species.parapatric.list[4]), label_x = -0.3)
ggsave("plots/joana.fd/mel.fdm.other.sp.allo.90th.fdm.max.min.png", width = 10, height = 9)




###### plot fd sims  east para / west /allo  SHDR ######
head(mel.fd.df.shdr.para.east.summ.df); head(sim.combo.fd.tmp.shdr.para.east.summ.list.list.df)

fd.para.east.sims.plot <- list(); for (i in 1:length(mel.other.species.parapatric.list)) {
  combo <- mel.other.species.parapatric.list[i]
  fd.para.east.sims.plot[[i]] <-  ggplot(subset(mel.fd.df, file.name.r==combo ), aes(x=fd)) + 
    labs(colour="Altitude")+
    geom_vline(xintercept = 0, colour="lightgrey", lty="dashed")+
    geom_vline(xintercept = subset(mel.fd.df.shdr.para.east.summ.df, file.name.r==combo & no.windows>1)$fd.max, colour="#049E73", lty="solid", alpha=.9, size=0.5)+
    geom_vline(xintercept = subset(mel.fd.df.shdr.para.east.summ.df, file.name.r==combo & no.windows>1 & shdr.para.east.id %in% subset(mel.outliersEast, overlaps.with.inversion=="yes")$shdr.para.east.id )$fd.max, 
               colour="yellow", lty="solid", alpha=1, size=0.7)+
    #geom_vline(xintercept = mean(subset(sim.combo.fd.tmp.shdr.para.east.summ.list.list.df.sub, file.name.r==combo)$fd.max), colour="black", lty="solid", alpha=.9, size=1.5)+
    geom_vline(xintercept = quantile(subset(sim.combo.fd.tmp.shdr.para.east.summ.list.list.df.sub, file.name.r==combo)$fd.max, 0.9), colour="black", lty="dashed", alpha=.9, size=1.5)+
    geom_vline(xintercept = quantile(subset(sim.combo.fd.tmp.shdr.para.east.summ.list.list.df.sub, file.name.r==combo)$fd.max, 0.1), colour="black", lty="dashed", alpha=.9, size=1.5)+
    geom_rect(aes(xmin=-0, xmax=.45, ymin=100, ymax=120 ), fill="white")+
    # geom_density(fill="transparent", alpha=0, size=1)+
    xlab(expression(paste("Topo2-Topo3",  italic(" H. melpomene"))))+
    theme_classic()+
    theme(axis.line.x = element_line(color="black", size = 0.5),axis.line.y = element_line(color="black", size = 0.5),axis.title.y = element_text(size=14),
          plot.margin = unit(c(0.5,0.2,0.2,0.2), "lines"),axis.text = element_text(size=10), #axis.text.x = element_blank(),
          axis.text.y = element_text(size=10),axis.title.x = element_blank(), 
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          #panel.background = element_rect(fill = "transparent"), plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
          legend.background = element_rect(fill = "transparent"), legend.box.background = element_rect(fill = "transparent"), strip.background.y = element_blank(),
          strip.text = element_text(size=14))+
    scale_y_continuous(limits = c(0,120), expand = c(0,0))+
    scale_x_continuous(limits = c(0,.45), expand = c(0,0)) }; 

t <- c(fd.para.east.sims.plot[2],fd.para.east.sims.plot[3],fd.para.east.sims.plot[1], fd.para.east.sims.plot[4])
plot_grid(plotlist = t  , ncol = 2,
          labels = list(mel.other.species.parapatric.list[2],mel.other.species.parapatric.list[3],mel.other.species.parapatric.list[1],mel.other.species.parapatric.list[4]), label_x = -0.3)
ggsave("plots/joana.fd/mel.fd.other.sp.para.east.90th.fd.max.min.png", width = 10, height = 9)


fd.para.west.sims.plot <- list(); for (i in 1:length(mel.other.species.parapatric.list)) {
  combo <- mel.other.species.parapatric.list[i]
  fd.para.west.sims.plot[[i]] <-  ggplot(subset(mel.fd.df, file.name.r==combo ), aes(x=fd)) + 
    labs(colour="Altitude")+
    geom_vline(xintercept = 0, colour="lightgrey", lty="dashed")+
    geom_vline(xintercept = subset(mel.fd.df.shdr.para.west.summ.df, file.name.r==combo & no.windows>1)$fd.max, colour="#0372B2", lty="solid", alpha=.9, size=0.5)+
    #geom_vline(xintercept = mean(subset(sim.combo.fd.tmp.shdr.para.west.summ.list.list.df.sub, file.name.r==combo)$fd.max), colour="black", lty="solid", alpha=.9, size=1.5)+
    geom_vline(xintercept = quantile(subset(sim.combo.fd.tmp.shdr.para.west.summ.list.list.df.sub, file.name.r==combo)$fd.max, 0.9), colour="black", lty="dashed", alpha=.9, size=1.5)+
    geom_vline(xintercept = quantile(subset(sim.combo.fd.tmp.shdr.para.west.summ.list.list.df.sub, file.name.r==combo)$fd.max, 0.1), colour="black", lty="dashed", alpha=.9, size=1.5)+
    geom_rect(aes(xmin=-0, xmax=.45, ymin=100, ymax=120 ), fill="white")+
    
    geom_density(fill="transparent", alpha=0, size=1)+
    xlab(expression(paste("Topo2-Topo3",  italic(" H. melpomene"))))+
    theme_classic()+
    theme(axis.line.x = element_line(color="black", size = 0.5),axis.line.y = element_line(color="black", size = 0.5),axis.title.y = element_text(size=14),
          plot.margin = unit(c(0.5,0.2,0.2,0.2), "lines"),axis.text = element_text(size=10), #axis.text.x = element_blank(),
          axis.text.y = element_text(size=10),axis.title.x = element_blank(), 
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          panel.background = element_rect(fill = "transparent"), plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
          legend.background = element_rect(fill = "transparent"), legend.box.background = element_rect(fill = "transparent"), strip.background.y = element_blank(),
          strip.text = element_text(size=14))+
    scale_y_continuous(limits = c(0,120), expand = c(0,0))+
    scale_x_continuous(limits = c(0,.45), expand = c(0,0)) } ;
t <- c(fd.para.west.sims.plot[2],fd.para.west.sims.plot[3],fd.para.west.sims.plot[1], fd.para.west.sims.plot[4])
plot_grid(plotlist = t  , ncol = 2,
          labels = list(mel.other.species.parapatric.list[2],mel.other.species.parapatric.list[3],mel.other.species.parapatric.list[1],mel.other.species.parapatric.list[4]), label_x = -0.3)
ggsave("plots/joana.fd/mel.fd.other.sp.para.west.90th.fd.max.min.png", width = 10, height = 9)


fd.allo.sims.plot <- list()
for (i in 1:length(mel.other.species.parapatric.list)) {
  combo <- mel.other.species.parapatric.list[i]
  fd.allo.sims.plot[[i]] <-  ggplot(subset(mel.fd.df, file.name.r==combo ), aes(x=fd)) + 
    labs(colour="Altitude")+
    geom_vline(xintercept = 0, colour="lightgrey", lty="dashed")+
    geom_vline(xintercept = subset(mel.fd.df.shdr.allo.summ.df, file.name.r==combo & no.windows>1)$fd.max, colour="#D65D00", lty="solid", alpha=.9, size=0.5)+
    geom_vline(xintercept = subset(mel.fd.df.shdr.allo.summ.df, file.name.r==combo & no.windows>1 & shdr.allo.id %in% subset(mel.shdr.allo.all.summ, shdr.para.east %in% subset(mel.shdr.para.east, overlaps.with.inversion=="yes")$shdr.para.east.id )$shdr.allo.id )$fd.max, 
               colour="yellow", lty="solid", alpha=.9, size=0.7)+
    #geom_vline(xintercept = mean(subset(sim.combo.fd.tmp.shdr.allo.summ.list.list.df.sub, file.name.r==combo)$fd.max), colour="black", lty="solid", alpha=.9, size=1.5)+
    geom_vline(xintercept = quantile(subset(sim.combo.fd.tmp.shdr.allo.summ.list.list.df.sub, file.name.r==combo)$fd.max, 0.9), colour="black", lty="dashed", alpha=.9, size=1.5)+
    geom_vline(xintercept = quantile(subset(sim.combo.fd.tmp.shdr.allo.summ.list.list.df.sub, file.name.r==combo)$fd.max, 0.1), colour="black", lty="dashed", alpha=.9, size=1.5)+
    geom_density(fill="transparent", alpha=0, size=1)+
    geom_rect(aes(xmin=-0, xmax=.45, ymin=100, ymax=120 ), fill="white")+
    
    xlab(expression(paste("Topo2-Topo3",  italic(" H. melpomene"))))+
    theme_classic()+
    theme(axis.line.x = element_line(color="black", size = 0.5),axis.line.y = element_line(color="black", size = 0.5),axis.title.y = element_text(size=14),
          plot.margin = unit(c(0.5,0.2,0.2,0.2), "lines"),axis.text = element_text(size=10), #axis.text.x = element_blank(),
          axis.text.y = element_text(size=10),axis.title.x = element_blank(), 
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          panel.background = element_rect(fill = "transparent"), plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
          legend.background = element_rect(fill = "transparent"), legend.box.background = element_rect(fill = "transparent"), strip.background.y = element_blank(),
          strip.text = element_text(size=14))+
    scale_y_continuous(limits = c(0,120), expand = c(0,0))+
    scale_x_continuous(limits = c(0,.45), expand = c(0,0)) }
t <- c(fd.allo.sims.plot[2],fd.allo.sims.plot[3],fd.allo.sims.plot[1], fd.allo.sims.plot[4])
plot_grid(plotlist = t  , ncol = 2, labels = list(mel.other.species.parapatric.list[2],mel.other.species.parapatric.list[3],mel.other.species.parapatric.list[1],mel.other.species.parapatric.list[4]), label_x = -0.3)
ggsave("plots/joana.fd/mel.fd.other.sp.allo.90th.fd.max.min.png", width = 10, height = 9)




############################ 3.5. assess which hdr >< simulation thresholds ############################
mel.other.species.parapatric.list <- unique(subset(mel.fd.df, type.comparison=="other.species.parapatric")$file.name.r); mel.other.species.parapatric.list

# max fd per hdr, max/min fdM per hdr
mel.fd.df.shdr.para.east.summ.df <- as.data.frame(mel.fd.df.shdr.para.east.summ.df)
mel.fd.df.shdr.para.east.summ.df$fdm.max.more.90th.perc.sims <- "NA"; mel.fd.df.shdr.para.east.summ.df$fdm.max.more.95th.perc.sims <- "NA"; mel.fd.df.shdr.para.east.summ.df$fd.max.more.90th.perc.sims <- "NA"; mel.fd.df.shdr.para.east.summ.df$fd.max.more.95th.perc.sims <- "NA"
mel.fd.df.shdr.para.east.summ.df$fdm.min.less.10th.perc.sims <- "NA"; mel.fd.df.shdr.para.east.summ.df$fdm.min.less.5th.perc.sims <- "NA"; mel.fd.df.shdr.para.east.summ.df$fd.min.less.10th.perc.sims <- "NA"; mel.fd.df.shdr.para.east.summ.df$fd.min.less.5th.perc.sims <- "NA"
mel.fd.df.shdr.para.east.summ.df$comb.new.name <- paste(mel.fd.df.shdr.para.east.summ.df$p1.newname, mel.fd.df.shdr.para.east.summ.df$p2.newname, mel.fd.df.shdr.para.east.summ.df$p3.newname, sep = "_")

mel.fd.df.shdr.para.west.summ.df <- as.data.frame(mel.fd.df.shdr.para.west.summ.df)
mel.fd.df.shdr.para.west.summ.df$fdm.max.more.90th.perc.sims <- "NA"; mel.fd.df.shdr.para.west.summ.df$fdm.max.more.95th.perc.sims <- "NA"; mel.fd.df.shdr.para.west.summ.df$fd.max.more.90th.perc.sims <- "NA"; mel.fd.df.shdr.para.west.summ.df$fd.max.more.95th.perc.sims <- "NA"
mel.fd.df.shdr.para.west.summ.df$fdm.min.less.10th.perc.sims <- "NA"; mel.fd.df.shdr.para.west.summ.df$fdm.min.less.5th.perc.sims <- "NA"; mel.fd.df.shdr.para.west.summ.df$fd.min.less.10th.perc.sims <- "NA"; mel.fd.df.shdr.para.west.summ.df$fd.min.less.5th.perc.sims <- "NA"
mel.fd.df.shdr.para.west.summ.df$comb.new.name <- paste(mel.fd.df.shdr.para.west.summ.df$p1.newname, mel.fd.df.shdr.para.west.summ.df$p2.newname, mel.fd.df.shdr.para.west.summ.df$p3.newname, sep = "_")

mel.fd.df.shdr.allo.summ.df <- as.data.frame(mel.fd.df.shdr.allo.summ.df)
mel.fd.df.shdr.allo.summ.df$fdm.max.more.90th.perc.sims <- "NA"; mel.fd.df.shdr.allo.summ.df$fdm.max.more.95th.perc.sims <- "NA"; mel.fd.df.shdr.allo.summ.df$fd.max.more.90th.perc.sims <- "NA"; mel.fd.df.shdr.allo.summ.df$fd.max.more.95th.perc.sims <- "NA"
mel.fd.df.shdr.allo.summ.df$fdm.min.less.10th.perc.sims <- "NA"; mel.fd.df.shdr.allo.summ.df$fdm.min.less.5th.perc.sims <- "NA"; mel.fd.df.shdr.allo.summ.df$fd.min.less.10th.perc.sims <- "NA"; mel.fd.df.shdr.allo.summ.df$fd.min.less.5th.perc.sims <- "NA"
mel.fd.df.shdr.allo.summ.df$comb.new.name <- paste(mel.fd.df.shdr.allo.summ.df$p1.newname, mel.fd.df.shdr.allo.summ.df$p2.newname, mel.fd.df.shdr.allo.summ.df$p3.newname, sep = "_")


for (i in 1:length(mel.other.species.parapatric.list)) {
  combo <- mel.other.species.parapatric.list[i]
  # east fdm.max, fdm.min, fd.max
  mel.fd.df.shdr.para.east.summ.df[mel.fd.df.shdr.para.east.summ.df$file.name.r==combo,]$fdm.max.more.90th.perc.sims <- if_else(mel.fd.df.shdr.para.east.summ.df[mel.fd.df.shdr.para.east.summ.df$file.name.r==combo,]$fdm.max > 
                                                                                                                                  quantile(subset(sim.combo.fd.tmp.shdr.para.east.summ.list.list.df.sub, file.name.r==combo)$fdm.max, 0.9), "yes", "no")
  mel.fd.df.shdr.para.east.summ.df[mel.fd.df.shdr.para.east.summ.df$file.name.r==combo,]$fdm.max.more.95th.perc.sims <- if_else(mel.fd.df.shdr.para.east.summ.df[mel.fd.df.shdr.para.east.summ.df$file.name.r==combo,]$fdm.max > 
                                                                                                                                  quantile(subset(sim.combo.fd.tmp.shdr.para.east.summ.list.list.df.sub, file.name.r==combo)$fdm.max, 0.95), "yes", "no")
  mel.fd.df.shdr.para.east.summ.df[mel.fd.df.shdr.para.east.summ.df$file.name.r==combo,]$fdm.min.less.10th.perc.sims  <- if_else(mel.fd.df.shdr.para.east.summ.df[mel.fd.df.shdr.para.east.summ.df$file.name.r==combo,]$fdm.min < 
                                                                                                                                   quantile(subset(sim.combo.fd.tmp.shdr.para.east.summ.list.list.df.sub, file.name.r==combo)$fdm.min, 0.1), "yes", "no")
  mel.fd.df.shdr.para.east.summ.df[mel.fd.df.shdr.para.east.summ.df$file.name.r==combo,]$fdm.min.less.5th.perc.sims  <- if_else(mel.fd.df.shdr.para.east.summ.df[mel.fd.df.shdr.para.east.summ.df$file.name.r==combo,]$fdm.min < 
                                                                                                                                  quantile(subset(sim.combo.fd.tmp.shdr.para.east.summ.list.list.df.sub, file.name.r==combo)$fdm.min, 0.05), "yes", "no")
  mel.fd.df.shdr.para.east.summ.df[mel.fd.df.shdr.para.east.summ.df$file.name.r==combo,]$fd.max.more.90th.perc.sims <- if_else(mel.fd.df.shdr.para.east.summ.df[mel.fd.df.shdr.para.east.summ.df$file.name.r==combo,]$fd.max > 
                                                                                                                                 quantile(subset(sim.combo.fd.tmp.shdr.para.east.summ.list.list.df.sub, file.name.r==combo)$fd.max, 0.9), "yes", "no")
  mel.fd.df.shdr.para.east.summ.df[mel.fd.df.shdr.para.east.summ.df$file.name.r==combo,]$fd.max.more.95th.perc.sims <- if_else(mel.fd.df.shdr.para.east.summ.df[mel.fd.df.shdr.para.east.summ.df$file.name.r==combo,]$fd.max > 
                                                                                                                                 quantile(subset(sim.combo.fd.tmp.shdr.para.east.summ.list.list.df.sub, file.name.r==combo)$fd.max, 0.95), "yes", "no")
  
  # west fdm.max, fdm.min, fd.max
  mel.fd.df.shdr.para.west.summ.df[mel.fd.df.shdr.para.west.summ.df$file.name.r==combo,]$fdm.max.more.90th.perc.sims <- if_else(mel.fd.df.shdr.para.west.summ.df[mel.fd.df.shdr.para.west.summ.df$file.name.r==combo,]$fdm.max > 
                                                                                                                                  quantile(subset(sim.combo.fd.tmp.shdr.para.west.summ.list.list.df.sub, file.name.r==combo)$fdm.max, 0.9), "yes", "no")
  mel.fd.df.shdr.para.west.summ.df[mel.fd.df.shdr.para.west.summ.df$file.name.r==combo,]$fdm.max.more.95th.perc.sims <- if_else(mel.fd.df.shdr.para.west.summ.df[mel.fd.df.shdr.para.west.summ.df$file.name.r==combo,]$fdm.max > 
                                                                                                                                  quantile(subset(sim.combo.fd.tmp.shdr.para.west.summ.list.list.df.sub, file.name.r==combo)$fdm.max, 0.95), "yes", "no")
  mel.fd.df.shdr.para.west.summ.df[mel.fd.df.shdr.para.west.summ.df$file.name.r==combo,]$fdm.min.less.10th.perc.sims  <- if_else(mel.fd.df.shdr.para.west.summ.df[mel.fd.df.shdr.para.west.summ.df$file.name.r==combo,]$fdm.min < 
                                                                                                                                   quantile(subset(sim.combo.fd.tmp.shdr.para.west.summ.list.list.df.sub, file.name.r==combo)$fdm.min, 0.1), "yes", "no")
  mel.fd.df.shdr.para.west.summ.df[mel.fd.df.shdr.para.west.summ.df$file.name.r==combo,]$fdm.min.less.5th.perc.sims  <- if_else(mel.fd.df.shdr.para.west.summ.df[mel.fd.df.shdr.para.west.summ.df$file.name.r==combo,]$fdm.min < 
                                                                                                                                  quantile(subset(sim.combo.fd.tmp.shdr.para.west.summ.list.list.df.sub, file.name.r==combo)$fdm.min, 0.05), "yes", "no")
  mel.fd.df.shdr.para.west.summ.df[mel.fd.df.shdr.para.west.summ.df$file.name.r==combo,]$fd.max.more.90th.perc.sims <- if_else(mel.fd.df.shdr.para.west.summ.df[mel.fd.df.shdr.para.west.summ.df$file.name.r==combo,]$fd.max > 
                                                                                                                                 quantile(subset(sim.combo.fd.tmp.shdr.para.west.summ.list.list.df.sub, file.name.r==combo)$fd.max, 0.9), "yes", "no")
  mel.fd.df.shdr.para.west.summ.df[mel.fd.df.shdr.para.west.summ.df$file.name.r==combo,]$fd.max.more.95th.perc.sims <- if_else(mel.fd.df.shdr.para.west.summ.df[mel.fd.df.shdr.para.west.summ.df$file.name.r==combo,]$fd.max > 
                                                                                                                                 quantile(subset(sim.combo.fd.tmp.shdr.para.west.summ.list.list.df.sub, file.name.r==combo)$fd.max, 0.95), "yes", "no")
  # allo fdm.max, fdm.min, fd.max
  mel.fd.df.shdr.allo.summ.df[mel.fd.df.shdr.allo.summ.df$file.name.r==combo,]$fdm.max.more.90th.perc.sims <- if_else(mel.fd.df.shdr.allo.summ.df[mel.fd.df.shdr.allo.summ.df$file.name.r==combo,]$fdm.max > 
                                                                                                                        quantile(subset(sim.combo.fd.tmp.shdr.allo.summ.list.list.df.sub, file.name.r==combo)$fdm.max, 0.9), "yes", "no")
  mel.fd.df.shdr.allo.summ.df[mel.fd.df.shdr.allo.summ.df$file.name.r==combo,]$fdm.max.more.95th.perc.sims <- if_else(mel.fd.df.shdr.allo.summ.df[mel.fd.df.shdr.allo.summ.df$file.name.r==combo,]$fdm.max > 
                                                                                                                        quantile(subset(sim.combo.fd.tmp.shdr.allo.summ.list.list.df.sub, file.name.r==combo)$fdm.max, 0.95), "yes", "no")
  mel.fd.df.shdr.allo.summ.df[mel.fd.df.shdr.allo.summ.df$file.name.r==combo,]$fdm.min.less.10th.perc.sims  <- if_else(mel.fd.df.shdr.allo.summ.df[mel.fd.df.shdr.allo.summ.df$file.name.r==combo,]$fdm.min < 
                                                                                                                         quantile(subset(sim.combo.fd.tmp.shdr.allo.summ.list.list.df.sub, file.name.r==combo)$fdm.min, 0.1), "yes", "no")
  mel.fd.df.shdr.allo.summ.df[mel.fd.df.shdr.allo.summ.df$file.name.r==combo,]$fdm.min.less.5th.perc.sims  <- if_else(mel.fd.df.shdr.allo.summ.df[mel.fd.df.shdr.allo.summ.df$file.name.r==combo,]$fdm.min < 
                                                                                                                        quantile(subset(sim.combo.fd.tmp.shdr.allo.summ.list.list.df.sub, file.name.r==combo)$fdm.min, 0.05), "yes", "no")
  
  mel.fd.df.shdr.allo.summ.df[mel.fd.df.shdr.allo.summ.df$file.name.r==combo,]$fd.max.more.90th.perc.sims <- if_else(mel.fd.df.shdr.allo.summ.df[mel.fd.df.shdr.allo.summ.df$file.name.r==combo,]$fd.max > 
                                                                                                                       quantile(subset(sim.combo.fd.tmp.shdr.allo.summ.list.list.df.sub, file.name.r==combo)$fd.max, 0.9), "yes", "no")
  mel.fd.df.shdr.allo.summ.df[mel.fd.df.shdr.allo.summ.df$file.name.r==combo,]$fd.max.more.95th.perc.sims <- if_else(mel.fd.df.shdr.allo.summ.df[mel.fd.df.shdr.allo.summ.df$file.name.r==combo,]$fd.max > 
                                                                                                                       quantile(subset(sim.combo.fd.tmp.shdr.allo.summ.list.list.df.sub, file.name.r==combo)$fd.max, 0.95), "yes", "no")
  
  
}

#### add to hdr outlier summaries #### 
# add allo 
mel.fd.df.shdr.allo.summ.df$shdr.para.east.id <- mel.shdr.allo.all.summ$shdr.para.east[match(mel.fd.df.shdr.allo.summ.df$shdr.allo.id,mel.shdr.allo.all.summ$shdr.allo.id )]
mel.fd.df.shdr.allo.summ.df$shdr.para.west.id <- mel.shdr.allo.all.summ$shdr.para.west[match(mel.fd.df.shdr.allo.summ.df$shdr.allo.id,mel.shdr.allo.all.summ$shdr.allo.id )]

shdr.mel.east.outlier.df <- read.csv("local/data/shdr.summ/shdr.mel.east.outlier.df.csv")
shdr.mel.west.outlier.df <- read.csv("local/data/shdr.summ/shdr.mel.west.outlier.df.csv")

for (i in 1:length(mel.other.species.parapatric.list)) {
  combo <- mel.other.species.parapatric.list[i]
  # 48 new columns, 8 combs x 2 (90,95th) x 3(fdm.max, fdm.min and fdmax), split into east west and allo, but add to para.east and para.west only
  
  ##### store fdm max min values per combo and cline ######
  shdr.mel.east.outlier.df[paste("fdm.max", combo, sep="_")] <- NA
  shdr.mel.east.outlier.df[paste("fdm.min", combo, sep="_")] <- NA
  # east
  shdr.mel.east.outlier.df[shdr.mel.east.outlier.df$is.allo.hdr=="no",][paste("fdm.max", combo, sep="_")] <- mel.fd.df.shdr.para.east.summ.df[mel.fd.df.shdr.para.east.summ.df$file.name.r==combo,]$fdm.max[
    match(shdr.mel.east.outlier.df[shdr.mel.east.outlier.df$is.allo.hdr=="no",]$shdr.para.east.id, mel.fd.df.shdr.para.east.summ.df[mel.fd.df.shdr.para.east.summ.df$file.name.r==combo,]$shdr.para.east.id )]
  shdr.mel.east.outlier.df[shdr.mel.east.outlier.df$is.allo.hdr=="no",][paste("fdm.min", combo, sep="_")] <- mel.fd.df.shdr.para.east.summ.df[mel.fd.df.shdr.para.east.summ.df$file.name.r==combo,]$fdm.min[
    match(shdr.mel.east.outlier.df[shdr.mel.east.outlier.df$is.allo.hdr=="no",]$shdr.para.east.id, mel.fd.df.shdr.para.east.summ.df[mel.fd.df.shdr.para.east.summ.df$file.name.r==combo,]$shdr.para.east.id )]

  # allo east
  shdr.mel.east.outlier.df[shdr.mel.east.outlier.df$is.allo.hdr=="yes",][paste("fdm.max", combo, sep="_")] <- mel.fd.df.shdr.allo.summ.df[mel.fd.df.shdr.allo.summ.df$file.name.r==combo,]$fdm.max[
    match(shdr.mel.east.outlier.df[shdr.mel.east.outlier.df$is.allo.hdr=="yes",]$shdr.para.east.id, mel.fd.df.shdr.allo.summ.df[mel.fd.df.shdr.allo.summ.df$file.name.r==combo,]$shdr.para.east.id )]
  shdr.mel.east.outlier.df[shdr.mel.east.outlier.df$is.allo.hdr=="yes",][paste("fdm.min", combo, sep="_")] <- mel.fd.df.shdr.allo.summ.df[mel.fd.df.shdr.allo.summ.df$file.name.r==combo,]$fdm.min[
    match(shdr.mel.east.outlier.df[shdr.mel.east.outlier.df$is.allo.hdr=="yes",]$shdr.para.east.id, mel.fd.df.shdr.allo.summ.df[mel.fd.df.shdr.allo.summ.df$file.name.r==combo,]$shdr.para.east.id )]

  # west
  shdr.mel.west.outlier.df[paste("fdm.max", combo, sep="_")] <- NA
  shdr.mel.west.outlier.df[paste("fdm.min", combo, sep="_")] <- NA
  shdr.mel.west.outlier.df[shdr.mel.west.outlier.df$is.allo.hdr=="no",][paste("fdm.max", combo, sep="_")] <- mel.fd.df.shdr.para.west.summ.df[mel.fd.df.shdr.para.west.summ.df$file.name.r==combo,]$fdm.max[
    match(shdr.mel.west.outlier.df[shdr.mel.west.outlier.df$is.allo.hdr=="no",]$shdr.para.west.id, mel.fd.df.shdr.para.west.summ.df[mel.fd.df.shdr.para.west.summ.df$file.name.r==combo,]$shdr.para.west.id )]
  shdr.mel.west.outlier.df[shdr.mel.west.outlier.df$is.allo.hdr=="no",][paste("fdm.min", combo, sep="_")] <- mel.fd.df.shdr.para.west.summ.df[mel.fd.df.shdr.para.west.summ.df$file.name.r==combo,]$fdm.min[
    match(shdr.mel.west.outlier.df[shdr.mel.west.outlier.df$is.allo.hdr=="no",]$shdr.para.west.id, mel.fd.df.shdr.para.west.summ.df[mel.fd.df.shdr.para.west.summ.df$file.name.r==combo,]$shdr.para.west.id )]
  
  # allo west
  shdr.mel.west.outlier.df[shdr.mel.west.outlier.df$is.allo.hdr=="yes",][paste("fdm.max", combo, sep="_")] <- mel.fd.df.shdr.allo.summ.df[mel.fd.df.shdr.allo.summ.df$file.name.r==combo,]$fdm.max[
    match(shdr.mel.west.outlier.df[shdr.mel.west.outlier.df$is.allo.hdr=="yes",]$shdr.para.west.id, mel.fd.df.shdr.allo.summ.df[mel.fd.df.shdr.allo.summ.df$file.name.r==combo,]$shdr.para.west.id )]
  shdr.mel.west.outlier.df[shdr.mel.west.outlier.df$is.allo.hdr=="yes",][paste("fdm.min", combo, sep="_")] <- mel.fd.df.shdr.allo.summ.df[mel.fd.df.shdr.allo.summ.df$file.name.r==combo,]$fdm.min[
    match(shdr.mel.west.outlier.df[shdr.mel.west.outlier.df$is.allo.hdr=="yes",]$shdr.para.west.id, mel.fd.df.shdr.allo.summ.df[mel.fd.df.shdr.allo.summ.df$file.name.r==combo,]$shdr.para.west.id )]
  
  ##### store whether they are outliers at diff threholds or not, yes no ######
  
  # east
  shdr.mel.east.outlier.df[paste("fdm.max.more.90th.perc.sims", combo, sep="_")] <- NA; shdr.mel.east.outlier.df[paste("fdm.max.more.95th.perc.sims", combo, sep="_")] <- NA
  shdr.mel.east.outlier.df[paste("fdm.min.less.10th.perc.sims", combo, sep="_")] <- NA; shdr.mel.east.outlier.df[paste("fdm.min.less.5th.perc.sims", combo, sep="_")] <- NA
  shdr.mel.east.outlier.df[shdr.mel.east.outlier.df$is.allo.hdr=="no",][paste("fdm.max.more.90th.perc.sims", combo, sep="_")] <- mel.fd.df.shdr.para.east.summ.df[mel.fd.df.shdr.para.east.summ.df$file.name.r==combo,]$fdm.max.more.90th.perc.sims[
    match(shdr.mel.east.outlier.df[shdr.mel.east.outlier.df$is.allo.hdr=="no",]$shdr.para.east.id, mel.fd.df.shdr.para.east.summ.df[mel.fd.df.shdr.para.east.summ.df$file.name.r==combo,]$shdr.para.east.id )]
  shdr.mel.east.outlier.df[shdr.mel.east.outlier.df$is.allo.hdr=="no",][paste("fdm.max.more.95th.perc.sims", combo, sep="_")] <- mel.fd.df.shdr.para.east.summ.df[mel.fd.df.shdr.para.east.summ.df$file.name.r==combo,]$fdm.max.more.95th.perc.sims[
    match(shdr.mel.east.outlier.df[shdr.mel.east.outlier.df$is.allo.hdr=="no",]$shdr.para.east.id, mel.fd.df.shdr.para.east.summ.df[mel.fd.df.shdr.para.east.summ.df$file.name.r==combo,]$shdr.para.east.id )]
  shdr.mel.east.outlier.df[shdr.mel.east.outlier.df$is.allo.hdr=="no",][paste("fdm.min.less.10th.perc.sims", combo, sep="_")] <- mel.fd.df.shdr.para.east.summ.df[mel.fd.df.shdr.para.east.summ.df$file.name.r==combo,]$fdm.min.less.10th.perc.sims[
    match(shdr.mel.east.outlier.df[shdr.mel.east.outlier.df$is.allo.hdr=="no",]$shdr.para.east.id, mel.fd.df.shdr.para.east.summ.df[mel.fd.df.shdr.para.east.summ.df$file.name.r==combo,]$shdr.para.east.id )]
  shdr.mel.east.outlier.df[shdr.mel.east.outlier.df$is.allo.hdr=="no",][paste("fdm.min.less.5th.perc.sims", combo, sep="_")] <- mel.fd.df.shdr.para.east.summ.df[mel.fd.df.shdr.para.east.summ.df$file.name.r==combo,]$fdm.min.less.5th.perc.sims[
    match(shdr.mel.east.outlier.df[shdr.mel.east.outlier.df$is.allo.hdr=="no",]$shdr.para.east.id, mel.fd.df.shdr.para.east.summ.df[mel.fd.df.shdr.para.east.summ.df$file.name.r==combo,]$shdr.para.east.id )]
  
  # allo east
  shdr.mel.east.outlier.df[shdr.mel.east.outlier.df$is.allo.hdr=="yes",][paste("fdm.max.more.90th.perc.sims", combo, sep="_")] <- mel.fd.df.shdr.allo.summ.df[mel.fd.df.shdr.allo.summ.df$file.name.r==combo,]$fdm.max.more.90th.perc.sims[
    match(shdr.mel.east.outlier.df[shdr.mel.east.outlier.df$is.allo.hdr=="yes",]$shdr.para.east.id, mel.fd.df.shdr.allo.summ.df[mel.fd.df.shdr.allo.summ.df$file.name.r==combo,]$shdr.para.east.id )]
  shdr.mel.east.outlier.df[shdr.mel.east.outlier.df$is.allo.hdr=="yes",][paste("fdm.max.more.95th.perc.sims", combo, sep="_")] <- mel.fd.df.shdr.allo.summ.df[mel.fd.df.shdr.allo.summ.df$file.name.r==combo,]$fdm.max.more.95th.perc.sims[
    match(shdr.mel.east.outlier.df[shdr.mel.east.outlier.df$is.allo.hdr=="yes",]$shdr.para.east.id, mel.fd.df.shdr.allo.summ.df[mel.fd.df.shdr.allo.summ.df$file.name.r==combo,]$shdr.para.east.id )]
  shdr.mel.east.outlier.df[shdr.mel.east.outlier.df$is.allo.hdr=="yes",][paste("fdm.min.less.10th.perc.sims", combo, sep="_")] <- mel.fd.df.shdr.allo.summ.df[mel.fd.df.shdr.allo.summ.df$file.name.r==combo,]$fdm.min.less.10th.perc.sims[
    match(shdr.mel.east.outlier.df[shdr.mel.east.outlier.df$is.allo.hdr=="yes",]$shdr.para.east.id, mel.fd.df.shdr.allo.summ.df[mel.fd.df.shdr.allo.summ.df$file.name.r==combo,]$shdr.para.east.id )]
  shdr.mel.east.outlier.df[shdr.mel.east.outlier.df$is.allo.hdr=="yes",][paste("fdm.min.less.5th.perc.sims", combo, sep="_")] <- mel.fd.df.shdr.allo.summ.df[mel.fd.df.shdr.allo.summ.df$file.name.r==combo,]$fdm.min.less.5th.perc.sims[
    match(shdr.mel.east.outlier.df[shdr.mel.east.outlier.df$is.allo.hdr=="yes",]$shdr.para.east.id, mel.fd.df.shdr.allo.summ.df[mel.fd.df.shdr.allo.summ.df$file.name.r==combo,]$shdr.para.east.id )]
  
  # west
  shdr.mel.west.outlier.df[paste("fdm.max.more.90th.perc.sims", combo, sep="_")] <- NA; shdr.mel.west.outlier.df[paste("fdm.max.more.95th.perc.sims", combo, sep="_")] <- NA
  shdr.mel.west.outlier.df[paste("fdm.min.less.10th.perc.sims", combo, sep="_")] <- NA; shdr.mel.west.outlier.df[paste("fdm.min.less.5th.perc.sims", combo, sep="_")] <- NA
  shdr.mel.west.outlier.df[shdr.mel.west.outlier.df$is.allo.hdr=="no",][paste("fdm.max.more.90th.perc.sims", combo, sep="_")] <- mel.fd.df.shdr.para.west.summ.df[mel.fd.df.shdr.para.west.summ.df$file.name.r==combo,]$fdm.max.more.90th.perc.sims[
    match(shdr.mel.west.outlier.df[shdr.mel.west.outlier.df$is.allo.hdr=="no",]$shdr.para.west.id, mel.fd.df.shdr.para.west.summ.df[mel.fd.df.shdr.para.west.summ.df$file.name.r==combo,]$shdr.para.west.id )]
  shdr.mel.west.outlier.df[shdr.mel.west.outlier.df$is.allo.hdr=="no",][paste("fdm.max.more.95th.perc.sims", combo, sep="_")] <- mel.fd.df.shdr.para.west.summ.df[mel.fd.df.shdr.para.west.summ.df$file.name.r==combo,]$fdm.max.more.95th.perc.sims[
    match(shdr.mel.west.outlier.df[shdr.mel.west.outlier.df$is.allo.hdr=="no",]$shdr.para.west.id, mel.fd.df.shdr.para.west.summ.df[mel.fd.df.shdr.para.west.summ.df$file.name.r==combo,]$shdr.para.west.id )]
  shdr.mel.west.outlier.df[shdr.mel.west.outlier.df$is.allo.hdr=="no",][paste("fdm.min.less.10th.perc.sims", combo, sep="_")] <- mel.fd.df.shdr.para.west.summ.df[mel.fd.df.shdr.para.west.summ.df$file.name.r==combo,]$fdm.min.less.10th.perc.sims[
    match(shdr.mel.west.outlier.df[shdr.mel.west.outlier.df$is.allo.hdr=="no",]$shdr.para.west.id, mel.fd.df.shdr.para.west.summ.df[mel.fd.df.shdr.para.west.summ.df$file.name.r==combo,]$shdr.para.west.id )]
  shdr.mel.west.outlier.df[shdr.mel.west.outlier.df$is.allo.hdr=="no",][paste("fdm.min.less.5th.perc.sims", combo, sep="_")] <- mel.fd.df.shdr.para.west.summ.df[mel.fd.df.shdr.para.west.summ.df$file.name.r==combo,]$fdm.min.less.5th.perc.sims[
    match(shdr.mel.west.outlier.df[shdr.mel.west.outlier.df$is.allo.hdr=="no",]$shdr.para.west.id, mel.fd.df.shdr.para.west.summ.df[mel.fd.df.shdr.para.west.summ.df$file.name.r==combo,]$shdr.para.west.id )]
  
  # allo west
  shdr.mel.west.outlier.df[shdr.mel.west.outlier.df$is.allo.hdr=="yes",][paste("fdm.max.more.90th.perc.sims", combo, sep="_")] <- mel.fd.df.shdr.allo.summ.df[mel.fd.df.shdr.allo.summ.df$file.name.r==combo,]$fdm.max.more.90th.perc.sims[
    match(shdr.mel.west.outlier.df[shdr.mel.west.outlier.df$is.allo.hdr=="yes",]$shdr.para.west.id, mel.fd.df.shdr.allo.summ.df[mel.fd.df.shdr.allo.summ.df$file.name.r==combo,]$shdr.para.west.id )]
  shdr.mel.west.outlier.df[shdr.mel.west.outlier.df$is.allo.hdr=="yes",][paste("fdm.max.more.95th.perc.sims", combo, sep="_")] <- mel.fd.df.shdr.allo.summ.df[mel.fd.df.shdr.allo.summ.df$file.name.r==combo,]$fdm.max.more.95th.perc.sims[
    match(shdr.mel.west.outlier.df[shdr.mel.west.outlier.df$is.allo.hdr=="yes",]$shdr.para.west.id, mel.fd.df.shdr.allo.summ.df[mel.fd.df.shdr.allo.summ.df$file.name.r==combo,]$shdr.para.west.id )]
  shdr.mel.west.outlier.df[shdr.mel.west.outlier.df$is.allo.hdr=="yes",][paste("fdm.min.less.10th.perc.sims", combo, sep="_")] <- mel.fd.df.shdr.allo.summ.df[mel.fd.df.shdr.allo.summ.df$file.name.r==combo,]$fdm.min.less.10th.perc.sims[
    match(shdr.mel.west.outlier.df[shdr.mel.west.outlier.df$is.allo.hdr=="yes",]$shdr.para.west.id, mel.fd.df.shdr.allo.summ.df[mel.fd.df.shdr.allo.summ.df$file.name.r==combo,]$shdr.para.west.id )]
  shdr.mel.west.outlier.df[shdr.mel.west.outlier.df$is.allo.hdr=="yes",][paste("fdm.min.less.5th.perc.sims", combo, sep="_")] <- mel.fd.df.shdr.allo.summ.df[mel.fd.df.shdr.allo.summ.df$file.name.r==combo,]$fdm.min.less.5th.perc.sims[
    match(shdr.mel.west.outlier.df[shdr.mel.west.outlier.df$is.allo.hdr=="yes",]$shdr.para.west.id, mel.fd.df.shdr.allo.summ.df[mel.fd.df.shdr.allo.summ.df$file.name.r==combo,]$shdr.para.west.id )]
  
}

names(shdr.mel.west.outlier.df)
test <- shdr.mel.west.outlier.df[,c(1,35,60:70)]; head(test)

# save
write.csv(shdr.mel.east.outlier.df , "local/data/shdr.summ/shdr.mel.east.outlier.df.csv", row.names = F)
write.csv(shdr.mel.west.outlier.df, "local/data/shdr.summ/shdr.mel.west.outlier.df.csv", row.names = F)


shdr.mel.east.outlier.df[paste("fdm.max", combo, sep="_")]
