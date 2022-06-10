##### obtain SHDRs #####
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



# ### regions ###
# era.e <- read.csv("../sharing/shared.e.regions.era.75.csv")
# head(era.e)
# 
# unique(era.e$scaff)
# ggplot(aes(x=BP.wg, y=zPBS0.1), data=subset(era.e, scaff=="Herato0209"))+
#   geom_point()


########## 0. data prep #####
setwd("git/2021-altitude-heliconius/")
ref.scaff.era <- read.table("local/data/ref/Heliconius_erato_demophoon_v1_-_scaffolds.fa.fai", row.names = NULL)
ref.scaff.mel <- read.table("local/data/ref/Hmel2.5.scaffolds.fa.fai", row.names = NULL)
# with mt
# ref.scaff.mel <- read.table("../../../../0.ref/Hmel2.5.fa.fai", row.names = NULL)

# read in files and modify accordingly 
era.tests <- list()
era.listpaths <- dir(pattern = ".*era.*txt"); era.listpaths
for (i in 1:4) {
  era.tests[[i]] <- read.table(era.listpaths[i])
  # fix headers
  names(era.tests[[i]]) <- lapply(era.tests[[i]][1,], as.character)
  era.tests[[i]] <- era.tests[[i]][-1,]
  era.tests[[i]]$SNP <- as.character(era.tests[[i]]$SNP)
  if(names(era.tests[[i]])[8]=="PBS0"){
    # make variables that will be used numeric
    era.tests[[i]]$PBS0 <- as.numeric(as.character(era.tests[[i]]$PBS0))
    era.tests[[i]]$PBS1 <- as.numeric(as.character(era.tests[[i]]$PBS1))
    era.tests[[i]]$PBS2 <- as.numeric(as.character(era.tests[[i]]$PBS2))
    era.tests[[i]]$PBE0 <- as.numeric(as.character(era.tests[[i]]$PBE0))
    era.tests[[i]]$zPBS0 <- zscore(era.tests[[i]]$PBS0)
    era.tests[[i]]$zPBS1 <- zscore(era.tests[[i]]$PBS1)
    era.tests[[i]]$zPBS2 <- zscore(era.tests[[i]]$PBS2)
    era.tests[[i]]$zPBE0 <- zscore(era.tests[[i]]$PBE0)
  } else if(names(era.tests[[i]])[5]=="fst"){
    era.tests[[i]]$fst <- as.numeric(as.character(era.tests[[i]]$fst))
    era.tests[[i]]$zfst <- zscore(era.tests[[i]]$fst)
  }
}

mel.tests <- list()
mel.listpaths <- dir(pattern = ".*mel.*txt"); mel.listpaths
for (i in 1:4) {
  mel.tests[[i]] <- read.table(mel.listpaths[i])
  names(mel.tests[[i]]) <- lapply(mel.tests[[i]][1,], as.character)
  mel.tests[[i]] <- mel.tests[[i]][-1,]
  mel.tests[[i]]$SNP <- as.character(mel.tests[[i]]$SNP)
  if(names(mel.tests[[i]])[8]=="PBS0"){
    mel.tests[[i]]$PBS0 <- as.numeric(as.character(mel.tests[[i]]$PBS0))
    mel.tests[[i]]$PBS1 <- as.numeric(as.character(mel.tests[[i]]$PBS1))
    mel.tests[[i]]$PBS2 <- as.numeric(as.character(mel.tests[[i]]$PBS2))
    mel.tests[[i]]$PBE0 <- as.numeric(as.character(mel.tests[[i]]$PBE0))
    mel.tests[[i]]$zPBS0 <- zscore(mel.tests[[i]]$PBS0)
    mel.tests[[i]]$zPBS1 <- zscore(mel.tests[[i]]$PBS1)
    mel.tests[[i]]$zPBS2 <- zscore(mel.tests[[i]]$PBS2)
    mel.tests[[i]]$zPBE0 <- zscore(mel.tests[[i]]$PBE0)
  } else if(names(mel.tests[[i]])[5]=="fst"){
    mel.tests[[i]]$fst <- as.numeric(as.character(mel.tests[[i]]$fst))
    mel.tests[[i]]$zfst <- zscore(mel.tests[[i]]$fst)
  }
}

# mel lift
mel.tests.lift  <- list()
mel.listpaths.lift <- dir(pattern = ".*mel.*txt", path = "liftover", full.names = T); mel.listpaths.lift
for (i in 1:4) {
  mel.tests.lift[[i]] <- read.table(mel.listpaths.lift[i])
  names(mel.tests.lift[[i]]) <- lapply(mel.tests.lift[[i]][1,], as.character)
  mel.tests.lift[[i]] <- mel.tests.lift[[i]][-1,]
  mel.tests.lift[[i]]$era.SNP <- as.character(mel.tests.lift[[i]]$era.SNP)
  mel.tests.lift[[i]]$mel.SNP <- as.character(mel.tests.lift[[i]]$mel.SNP)
  if(names(mel.tests.lift[[i]])[10]=="PBS0"){
    mel.tests.lift[[i]]$PBS0 <- as.numeric(as.character(mel.tests.lift[[i]]$PBS0))
    #mel.tests.lift[[i]]$PBE0 <- as.numeric(as.character(mel.tests.lift[[i]]$PBE0))
    mel.tests.lift[[i]]$zPBS0 <- zscore(mel.tests.lift[[i]]$PBS0)
  } else if(names(mel.tests.lift[[i]])[7]=="fst"){
    mel.tests.lift[[i]]$fst <- as.numeric(as.character(mel.tests.lift[[i]]$fst))
    mel.tests.lift[[i]]$zfst <- zscore(mel.tests.lift[[i]]$fst)
  }
}


## add chr.bp variable to tests for later sharing (my method)
# use original datasets for sharing estimations
# BP is continuous within chromosomes, so just add 0 dummy zeroes between chromosome and BP
for (i in 1:4) {
  era.tests[[i]]$chr.bp <- paste(era.tests[[i]]$CHR, "0000", era.tests[[i]]$BP, sep = "")}
for (i in 1:4) {
  mel.tests[[i]]$chr.bp <- paste(mel.tests[[i]]$CHR, "0000", mel.tests[[i]]$BP, sep = "")}

#### BP whole genome offsets ######
# create BP.wg first, use all along- STICHING scaff together
# CHANGE HERE to using chr.bp instead of BP.wg if we dont want chr/scaffold stiched
# get scaffold offsets from fai
scafEnds.era <- cumsum(ref.scaff.era[,2]); offset.era <- scafEnds.era - ref.scaff.era[,2]
plot(scafEnds.era,offset.era); names(offset.era) <- ref.scaff.era[,1]
scafEnds.mel <- cumsum(ref.scaff.mel[,2]); offset.mel <- scafEnds.mel - ref.scaff.mel[,2]
plot(scafEnds.mel,offset.mel); names(offset.mel) <- ref.scaff.mel[,1]
# check levels present in data
levels(era.tests[[1]]$scaff); levels(mel.tests[[1]]$scaff)
#  sub offset by scaffolds that are present
offset.era.sub <- offset.era[names(offset.era) %in% unique(era.tests[[1]]$scaff)]
offset.mel.sub <- offset.mel[names(offset.mel) %in% unique(mel.tests[[1]]$scaff)]
offset.mel.sub.lift <- offset.era[names(offset.era) %in% unique(mel.tests.lift[[1]]$era.scaff)]


# adjust midpos of window by adding whole genome offset, to make whole genome BP - to all tests
# before we had added scaff offset within chromosome (here no need for BP)
for (i in 1:4) {
  era.tests[[i]]$BP.wg <- as.numeric(as.character(era.tests[[i]]$midPos)) + offset.era.sub[era.tests[[i]]$scaff]}
for (i in 1:4) {
  mel.tests[[i]]$BP.wg <- as.numeric(as.character(mel.tests[[i]]$midPos)) + offset.mel.sub[mel.tests[[i]]$scaff]}

for (i in 1:4) {
  # drop weird "era.scaff" level
  mel.tests.lift[[i]]$era.scaff <- droplevels(subset(mel.tests.lift[[i]], era.scaff!="era.scaff", drop = TRUE)$era.scaff)
  mel.tests.lift[[i]]$era.BP.wg <- as.numeric(as.character(mel.tests.lift[[i]]$midPos)) + offset.mel.sub.lift[mel.tests.lift[[i]]$era.scaff]}


###### above 4std from mean ######

### era zPBS0
era.tests.top1.pbs0 <- list()
for (i in seq_along(era.tests)) {
  n<-1
  if(names(era.tests[[i]])[8]=="PBS0"){
    # get sites above 4std from mean
    era.tests.top1.pbs0[[i]] <- subset(era.tests[[i]], zPBS0 > 4)
  } else if(names(era.tests[[i]])[5]=="fst"){
    # get sites above 4std from mean
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


### BUFFERS EXPAND #######
# BP is continuous within chromosome, not within scaffold- so best use this!! already numeric, just need to add the CHR value
# use BP+- alone, until doing the sharing across datasets, then do share BPNUM+CHR (concatenated)
## list of top snps +-20 sites, have 1 list per test

#generic function that expands dataset +-x windows
bp.list.expanded <- function(dat, x){
  # use chr.bp which has 0000 added between scaffolds
  chr.bp <- dat[["chr.bp"]]
  # x= windows (1kb), convert to bp
  y <- x*1000
  seq(from=as.numeric(as.character(chr.bp))-y, to=as.numeric(as.character(chr.bp))+y, by=1000)}

list.pos.to.include <- list()
era.tests.top1.pbs0.expanded50 <- list()
for (i in 1:4) {
  list.pos.to.include[[i]] <- unique((apply(era.tests.top1.pbs0[[i]], 1 , bp.list.expanded, 50))[TRUE])
  era.tests.top1.pbs0.expanded50[[i]] <-subset(era.tests[[i]], chr.bp %in% list.pos.to.include[[i]] )
}

# mel- - buffer50 top1
list.pos.to.include <- list()
mel.tests.top1.pbs0.expanded50 <- list()
for (i in 1:4) {
  list.pos.to.include[[i]] <- unique((apply(mel.tests.top1.pbs0[[i]], 1 , bp.list.expanded, 50))[TRUE])
  mel.tests.top1.pbs0.expanded50[[i]] <-subset(mel.tests[[i]], chr.bp %in% list.pos.to.include[[i]] )
}



#########  all windows pbs/fst  #############################################
## prep shared regions data ####
## all the same positions, but PBS0 different in each
# grab positions that in average, for all, have PBS0>4std
era.all <- era.tests[[1]][c("scaff", "midPos","BP.wg", "zPBS0")]; names(era.all)[4] <- "zPBS0.1"
era.all$zfst.2 <-era.tests[[2]]$zfst[match(era.all$BP.wg, era.tests[[2]]$BP.wg)]; era.all$zPBS0.3 <-era.tests[[3]]$zPBS0[match(era.all$BP.wg, era.tests[[3]]$BP.wg)]
era.all$zPBS0.4 <-era.tests[[4]]$zPBS0[match(era.all$BP.wg, era.tests[[4]]$BP.wg)]; era.all$fst.2 <-era.tests[[2]]$fst[match(era.all$BP.wg, era.tests[[2]]$BP.wg)]
era.all$Fst01.1 <-era.tests[[1]]$Fst01[match(era.all$BP.wg, era.tests[[1]]$BP.wg)]; era.all$Fst01.3 <-era.tests[[3]]$Fst01[match(era.all$BP.wg, era.tests[[3]]$BP.wg)]; era.all$Fst01.4 <-era.tests[[4]]$Fst01[match(era.all$BP.wg, era.tests[[4]]$BP.wg)]
era.all$Fst02.1 <-era.tests[[1]]$Fst02[match(era.all$BP.wg, era.tests[[1]]$BP.wg)]; era.all$Fst02.3 <-era.tests[[3]]$Fst02[match(era.all$BP.wg, era.tests[[3]]$BP.wg)]; era.all$Fst02.4 <-era.tests[[4]]$Fst02[match(era.all$BP.wg, era.tests[[4]]$BP.wg)]
era.all$Fst12.1 <-era.tests[[1]]$Fst12[match(era.all$BP.wg, era.tests[[1]]$BP.wg)]; era.all$Fst12.3 <-era.tests[[3]]$Fst12[match(era.all$BP.wg, era.tests[[3]]$BP.wg)]; era.all$Fst12.4 <-era.tests[[4]]$Fst12[match(era.all$BP.wg, era.tests[[4]]$BP.wg)]
era.all$zPBE0.1 <-era.tests[[1]]$zPBE0[match(era.all$BP.wg, era.tests[[1]]$BP.wg)]; era.all$zPBE0.3 <-era.tests[[3]]$zPBE0[match(era.all$BP.wg, era.tests[[3]]$BP.wg)]; era.all$zPBE0.4 <-era.tests[[4]]$zPBE0[match(era.all$BP.wg, era.tests[[4]]$BP.wg)]
era.all$zPBS1.1 <-era.tests[[1]]$zPBS1[match(era.all$BP.wg, era.tests[[1]]$BP.wg)]; era.all$zPBS1.3 <-era.tests[[3]]$zPBS1[match(era.all$BP.wg, era.tests[[3]]$BP.wg)]; era.all$zPBS1.4 <-era.tests[[4]]$zPBS1[match(era.all$BP.wg, era.tests[[4]]$BP.wg)]
era.all$zPBS2.1 <-era.tests[[1]]$zPBS2[match(era.all$BP.wg, era.tests[[1]]$BP.wg)]; era.all$zPBS2.3 <-era.tests[[3]]$zPBS2[match(era.all$BP.wg, era.tests[[3]]$BP.wg)]; era.all$zPBS2.4 <-era.tests[[4]]$zPBS2[match(era.all$BP.wg, era.tests[[4]]$BP.wg)]; era.all

# use sharing sims intervals approach (overlaps)
# ## add region ID to each block
# unique(sites.east.era$shared.type)
# era.all$east.sharing.type <- sites.east.era$shared.type[match(era.all$BP.wg, sites.east.era$BP.wg)]
# era.all$east.sharing.id <-sites.east.era$shared.id[match(era.all$BP.wg, sites.east.era$BP.wg)]
# era.all$west.sharing.type <- sites.west.era$shared.type[match(era.all$BP.wg, sites.west.era$BP.wg)]
# era.all$west.sharing.id <-sites.west.era$shared.id[match(era.all$BP.wg, sites.west.era$BP.wg)]; unique(era.all$west.sharing.id)
# 
# head(era.all); unique(era.all$east.sharing.type)

# era.all.long <- gather(era.all, pbs, value, -scaff, -BP.wg, -BP.wg); era.all.long
# era.all.long$pbs.type <- substr(era.all.long$pbs,0,5); unique(era.all.long$pbs.type )
# era.all.long$pop <-sapply(strsplit(era.all.long$pbs,".", fixed=TRUE), "[[", 2); unique(era.all.long$pop )
# era.all.long$value.type <- if_else(substr(era.all.long$pbs.type,0,1)=="z", "z", "raw"); unique(era.all.long$value.type )
write.csv(era.all, "local/data/pbs.out/era.all.pbs.csv", row.names = F)
str(era.all)

mel.all <- mel.tests[[1]][c("scaff", "midPos","BP.wg", "zPBS0")]; names(mel.all)[4] <- "zPBS0.1"
mel.all$zfst.4 <-mel.tests[[4]]$zfst[match(mel.all$BP.wg, mel.tests[[4]]$BP.wg)]; mel.all$zPBS0.3 <-mel.tests[[3]]$zPBS0[match(mel.all$BP.wg, mel.tests[[3]]$BP.wg)]
mel.all$zPBS0.2 <-mel.tests[[2]]$zPBS0[match(mel.all$BP.wg, mel.tests[[2]]$BP.wg)]; mel.all$fst.4 <-mel.tests[[4]]$fst[match(mel.all$BP.wg, mel.tests[[4]]$BP.wg)]
mel.all$Fst01.1 <-mel.tests[[1]]$Fst01[match(mel.all$BP.wg, mel.tests[[1]]$BP.wg)]; mel.all$Fst01.3 <-mel.tests[[3]]$Fst01[match(mel.all$BP.wg, mel.tests[[3]]$BP.wg)]; mel.all$Fst01.2 <-mel.tests[[2]]$Fst01[match(mel.all$BP.wg, mel.tests[[2]]$BP.wg)]
mel.all$Fst02.1 <-mel.tests[[1]]$Fst02[match(mel.all$BP.wg, mel.tests[[1]]$BP.wg)]; mel.all$Fst02.3 <-mel.tests[[3]]$Fst02[match(mel.all$BP.wg, mel.tests[[3]]$BP.wg)]; mel.all$Fst02.2 <-mel.tests[[2]]$Fst02[match(mel.all$BP.wg, mel.tests[[2]]$BP.wg)]
mel.all$Fst12.1 <-mel.tests[[1]]$Fst12[match(mel.all$BP.wg, mel.tests[[1]]$BP.wg)]; mel.all$Fst12.3 <-mel.tests[[3]]$Fst12[match(mel.all$BP.wg, mel.tests[[3]]$BP.wg)]; mel.all$Fst12.2 <-mel.tests[[2]]$Fst12[match(mel.all$BP.wg, mel.tests[[2]]$BP.wg)]
mel.all$zPBE0.1 <-mel.tests[[1]]$zPBE0[match(mel.all$BP.wg, mel.tests[[1]]$BP.wg)]; mel.all$zPBE0.3 <-mel.tests[[3]]$zPBE0[match(mel.all$BP.wg, mel.tests[[3]]$BP.wg)]; mel.all$zPBE0.2 <-mel.tests[[2]]$zPBE0[match(mel.all$BP.wg, mel.tests[[2]]$BP.wg)]
mel.all$zPBS1.1 <-mel.tests[[1]]$zPBS1[match(mel.all$BP.wg, mel.tests[[1]]$BP.wg)]; mel.all$zPBS1.3 <-mel.tests[[3]]$zPBS1[match(mel.all$BP.wg, mel.tests[[3]]$BP.wg)]; mel.all$zPBS1.2 <-mel.tests[[2]]$zPBS1[match(mel.all$BP.wg, mel.tests[[2]]$BP.wg)]
mel.all$zPBS2.1 <-mel.tests[[1]]$zPBS2[match(mel.all$BP.wg, mel.tests[[1]]$BP.wg)]; mel.all$zPBS2.3 <-mel.tests[[3]]$zPBS2[match(mel.all$BP.wg, mel.tests[[3]]$BP.wg)]; mel.all$zPBS2.2 <-mel.tests[[2]]$zPBS2[match(mel.all$BP.wg, mel.tests[[2]]$BP.wg)]; mel.all

# add erato scaff, start and end positions to mel gwas file
mel.all$midPos.era <- mel.tests.lift[[1]]$midPos[match(paste(mel.all$scaff, mel.all$midPos, sep="_"), paste(mel.tests.lift[[1]]$mel.scaff, mel.tests.lift[[1]]$HmelmidPos, sep="_"))]
mel.all$scaff.era <- mel.tests.lift[[1]]$era.scaff[match(paste(mel.all$scaff, mel.all$midPos, sep="_"), paste(mel.tests.lift[[1]]$mel.scaff, mel.tests.lift[[1]]$HmelmidPos, sep="_"))]
mel.all$era.BP <- mel.tests.lift[[1]]$era.BP[match(paste(mel.all$scaff, mel.all$midPos, sep="_"), paste(mel.tests.lift[[1]]$mel.scaff, mel.tests.lift[[1]]$HmelmidPos, sep="_"))]
mel.all$era.CHR <- mel.tests.lift[[1]]$era.CHR[match(paste(mel.all$scaff, mel.all$midPos, sep="_"), paste(mel.tests.lift[[1]]$mel.scaff, mel.tests.lift[[1]]$HmelmidPos, sep="_"))]
mel.all$era.BP.wg <- mel.tests.lift[[1]]$era.BP.wg[match(paste(mel.all$scaff, mel.all$midPos, sep="_"), paste(mel.tests.lift[[1]]$mel.scaff, mel.tests.lift[[1]]$HmelmidPos, sep="_"))]

head(mel.all)
# mel.all.long <- gather(mel.all, pbs, value, -scaff, -BP.wg, -BP.wg); mel.all.long
# mel.all.long$pbs.type <- substr(mel.all.long$pbs,0,5); unique(mel.all.long$pbs.type )
# mel.all.long$pop <-sapply(strsplit(mel.all.long$pbs,".", fixed=TRUE), "[[", 2); unique(mel.all.long$pop )
# mel.all.long$value.type <- if_else(substr(mel.all.long$pbs.type,0,1)=="z", "z", "raw"); unique(mel.all.long$value.type )
write.csv(mel.all, "local/data/pbs.out/mel.all.pbs.csv", row.names = F)


