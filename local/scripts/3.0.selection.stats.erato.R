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
setwd("/git/2021-altitude-heliconius/")
pop.info <- read.csv("02.info/pop.short.info.csv")

# ## all windows
# era.all.pbs <- read.csv("local/data/pbs.out/era.all.pbs.recomb.csv")
# mel.all.pbs <- read.csv("local/data/pbs.out/mel.all.pbs.recomb.csv")
# 
####### sort shdr para, allo #######
# era.shdr.para.east <- read.csv("local/data/sharing/era.shdr.para.east.csv")
# era.shdr.para.west <- read.csv("local/data/sharing/era.shdr.para.west.csv")
# era.shdr.allo.all  <- read.csv("local/data/sharing/era.shdr.allo.all.csv")
# 
# era.shdr.para.east$scaff <- era.all.pbs$scaff[match(era.shdr.para.east$start, era.all.pbs$BP.wg)]
# era.shdr.para.west$scaff <- era.all.pbs$scaff[match(era.shdr.para.west$start, era.all.pbs$BP.wg)]
# era.shdr.allo.all$scaff <- era.all.pbs$scaff[match(era.shdr.allo.all$start, era.all.pbs$BP.wg)]
# 
# era.shdr.para.east$chr <- as.numeric(as.character(substr(era.shdr.para.east$scaff, 7,8)))
# era.shdr.para.west$chr <- as.numeric(as.character(substr(era.shdr.para.west$scaff, 7,8)))
# era.shdr.allo.all$chr <- as.numeric(as.character(substr(era.shdr.allo.all$scaff, 7,8)))
# 
# era.shdr.para.east$start.bp <- era.all.pbs$midPos[match(era.shdr.para.east$start, era.all.pbs$BP.wg)]
# era.shdr.para.west$start.bp <- era.all.pbs$midPos[match(era.shdr.para.west$start, era.all.pbs$BP.wg)]
# era.shdr.allo.all$start.bp <-era.all.pbs$midPos[match(era.shdr.allo.all$start, era.all.pbs$BP.wg)]
# era.shdr.para.east$end.bp <- era.all.pbs$midPos[match(era.shdr.para.east$end, era.all.pbs$BP.wg)]
# era.shdr.para.west$end.bp <- era.all.pbs$midPos[match(era.shdr.para.west$end, era.all.pbs$BP.wg)]
# era.shdr.allo.all$end.bp <-era.all.pbs$midPos[match(era.shdr.allo.all$end, era.all.pbs$BP.wg)]
# 
# # obtain position of max PBS0
# names(era.all.pbs)
# era.shdr.para.east.max <- summarise(group_by(subset(era.all.pbs, shdr.para.east.id!=""), shdr.para.east.id ),
#           zPBS.hig.co.e.max.value=max(zPBS0.1),
#           zPBS.hig.co.e.max.pos.BP.wg=BP.wg[zPBS0.1==zPBS.hig.co.e.max.value],
#           zPBS.hig.co.e.max.pos.midPos=midPos[zPBS0.1==zPBS.hig.co.e.max.value],
#           zPBS.hig.ec.e.max.value=max(zPBS0.3),
#           zPBS.hig.ec.e.max.pos.BP.wg=BP.wg[zPBS0.3==zPBS.hig.ec.e.max.value],
#           zPBS.hig.ec.e.max.pos.midPos=midPos[zPBS0.3==zPBS.hig.ec.e.max.value],); era.shdr.para.east.max
# era.shdr.para.west.max <- summarise(group_by(subset(era.all.pbs, shdr.para.west.id!="" & zfst.2!="" & zPBS0.4!="" ), shdr.para.west.id ),
#                                 zFst.co.w.max.value=max(zfst.2, na.rm = T),
#                                 zFst.co.w.max.pos.BP.wg=BP.wg[zfst.2==zFst.co.w.max.value],
#                                 zFst.co.w.max.pos.midPos=midPos[zfst.2==zFst.co.w.max.value],
#                                 zPBS.hig.ec.w.max.value=max(zPBS0.4, na.rm = T),
#                                 zPBS.hig.ec.w.max.pos.BP.wg=BP.wg[zPBS0.4==zPBS.hig.ec.w.max.value],
#                                 zPBS.hig.ec.w.max.pos.midPos=midPos[zPBS0.4==zPBS.hig.ec.w.max.value],); era.shdr.para.west.max
# 
# era.shdr.para.east <- merge(era.shdr.para.east, era.shdr.para.east.max, by ="shdr.para.east.id" )
# era.shdr.para.west <- merge(era.shdr.para.west, era.shdr.para.west.max, by ="shdr.para.west.id" )
# 
# #careful with tiny hdrs, remove from future
# era.shdr.para.east$no.windows <- (era.shdr.para.east$end - era.shdr.para.east$start)/1000
# era.shdr.para.west$no.windows <- (era.shdr.para.west$end - era.shdr.para.west$start)/1000
# 
# mel.shdr.para.east <- read.csv("local/data/sharing/mel.shdr.para.east.csv")
# mel.shdr.para.west <- read.csv("local/data/sharing/mel.shdr.para.west.csv")
# mel.shdr.allo.all  <- read.csv("local/data/sharing/mel.shdr.allo.all.csv")
# 
# mel.shdr.para.east$scaff <- mel.all.pbs$scaff[match(mel.shdr.para.east$start, mel.all.pbs$BP.wg)]
# mel.shdr.para.west$scaff <- mel.all.pbs$scaff[match(mel.shdr.para.west$start, mel.all.pbs$BP.wg)]
# mel.shdr.allo.all$scaff <- mel.all.pbs$scaff[match(mel.shdr.allo.all$start, mel.all.pbs$BP.wg)]
# 
# mel.shdr.para.east$chr <- as.numeric(as.character(substr(mel.shdr.para.east$scaff, 6,7)))
# mel.shdr.para.west$chr <- as.numeric(as.character(substr(mel.shdr.para.west$scaff,  6,7)))
# mel.shdr.allo.all$chr <- as.numeric(as.character(substr(mel.shdr.allo.all$scaff, 6,7)))
# 
# mel.shdr.para.east$start.bp <- mel.all.pbs$midPos[match(mel.shdr.para.east$start, mel.all.pbs$BP.wg)]
# mel.shdr.para.west$start.bp <- mel.all.pbs$midPos[match(mel.shdr.para.west$start, mel.all.pbs$BP.wg)]
# mel.shdr.allo.all$start.bp <-mel.all.pbs$midPos[match(mel.shdr.allo.all$start, mel.all.pbs$BP.wg)]
# mel.shdr.para.east$end.bp <- mel.all.pbs$midPos[match(mel.shdr.para.east$end, mel.all.pbs$BP.wg)]
# mel.shdr.para.west$end.bp <- mel.all.pbs$midPos[match(mel.shdr.para.west$end, mel.all.pbs$BP.wg)]
# mel.shdr.allo.all$end.bp <-mel.all.pbs$midPos[match(mel.shdr.allo.all$end, mel.all.pbs$BP.wg)]
# 
# # obtain position of max PBS0
# mel.shdr.para.east.max <- summarise(group_by(subset(mel.all.pbs, shdr.para.east.id!=""), shdr.para.east.id ),
#                                     zPBS.hig.co.e.max.value=max(zPBS0.1),
#                                     zPBS.hig.co.e.max.pos.BP.wg=BP.wg[zPBS0.1==zPBS.hig.co.e.max.value],
#                                     zPBS.hig.co.e.max.pos.midPos=midPos[zPBS0.1==zPBS.hig.co.e.max.value],
#                                     zPBS.hig.ec.e.max.value=max(zPBS0.3),
#                                     zPBS.hig.ec.e.max.pos.BP.wg=BP.wg[zPBS0.3==zPBS.hig.ec.e.max.value],
#                                     zPBS.hig.ec.e.max.pos.midPos=midPos[zPBS0.3==zPBS.hig.ec.e.max.value],); mel.shdr.para.east.max
# mel.shdr.para.west.max <- summarise(group_by(subset(mel.all.pbs, shdr.para.west.id!=""), shdr.para.west.id ),
#                                     zFst.ec.w.max.value=max(zfst.4),
#                                     zFst.ec.w.max.pos.BP.wg=BP.wg[zfst.4==zFst.ec.w.max.value],
#                                     zFst.ec.w.max.pos.midPos=midPos[zfst.4==zFst.ec.w.max.value],
#                                     zPBS.hig.co.w.max.value=max(zPBS0.2),
#                                     zPBS.hig.co.w.max.pos.BP.wg=BP.wg[zPBS0.2==zPBS.hig.co.w.max.value],
#                                     zPBS.hig.co.w.max.pos.midPos=midPos[zPBS0.2==zPBS.hig.co.w.max.value],); mel.shdr.para.west.max
# 
# mel.shdr.para.east <- merge(mel.shdr.para.east, mel.shdr.para.east.max, by ="shdr.para.east.id" )
# mel.shdr.para.west <- merge(mel.shdr.para.west, mel.shdr.para.west.max, by ="shdr.para.west.id" )
# mel.shdr.para.east$no.windows <- (mel.shdr.para.east$end - mel.shdr.para.east$start)/1000
# mel.shdr.para.west$no.windows <- (mel.shdr.para.west$end - mel.shdr.para.west$start)/1000
# 
# 
# write.csv(era.shdr.para.east , "local/data/sharing/era.shdr.para.east.csv", row.names = F)
# write.csv(era.shdr.para.west , "local/data/sharing/era.shdr.para.west.csv", row.names = F)
# write.csv(era.shdr.allo.all , "local/data/sharing/era.shdr.allo.all.csv", row.names = F)
# 
# write.csv(mel.shdr.para.east , "local/data/sharing/mel.shdr.para.east.csv", row.names = F)
# write.csv(mel.shdr.para.west , "local/data/sharing/mel.shdr.para.west.csv", row.names = F)
# write.csv(mel.shdr.allo.all , "local/data/sharing/mel.shdr.allo.all.csv", row.names = F)


####### read in formatted #######
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
ref.scaff.era <- read.table("local/data/ref/Heliconius_erato_demophoon_v1_-_scaffolds.fa.fai", row.names = NULL)
ref.scaff.mel <- read.table("local/data/ref/Hmel2.5.scaffolds.fa.fai", row.names = NULL)
scafEnds.era <- cumsum(ref.scaff.era[,2]); offset.era <- scafEnds.era - ref.scaff.era[,2]
scafEnds.mel <- cumsum(ref.scaff.mel[,2]); offset.mel <- scafEnds.mel - ref.scaff.mel[,2]
era.genome_size <- tail(scafEnds.era,1)
mel.genome_size <- tail(scafEnds.mel,1)

#######  sort thetas ######
# ### load all subpop thetas #####
# listpaths <- dir(pattern = ".*PG", path = "local/data/local/data/output/thetas/by.subpop.folded", full.names = T); listpaths
# thetas <- list()
# for (i in 1:length(listpaths)) {
#   thetas[[i]] <- read.table(listpaths[i], header = TRUE)}
# for (i in 1:length(listpaths)) {
#   names(thetas[[i]]) <- c("(indexStart,indexStop)(firstPos_withData,lastPos_withData)(WinStart,WinStop)",	"scaff",	"WinCenter",	"tW",	"tP",	"tF",	"tH",	"tL",	"Tajima",	"fuf",	"fud",	"fayh",	"zeng",	"nSites")}
# 
# # divide by species
# thetas.era <- list(); thetas.era <- thetas[c(1:18)]
# thetas.mel<- list(); thetas.mel <- thetas[c(19:34)]
# 
# ## add BP.wg, remove lines from unassigned chr in mel, add sharing ids
# # add cline info
# listpaths.era <- listpaths[1:18]
# listpaths.mel <- listpaths[19:34]
# 
# for (i in 1:length(thetas.era)) {
#   thetas.era[[i]]$chr <- as.numeric(as.character(substr(thetas.era[[i]]$scaff, 7,8)))
#   thetas.era[[i]]$BP.wg <- era.all.pbs$BP.wg[match(paste(thetas.era[[i]]$scaff, thetas.era[[i]]$WinCenter), paste(era.all.pbs$scaff, era.all.pbs$midPos))]
#   thetas.era[[i]]$shdr.para.east.id <- era.all.pbs$shdr.para.east.id[match(thetas.era[[i]]$BP.wg, era.all.pbs$BP.wg)]
#   thetas.era[[i]]$shdr.para.west.id <- era.all.pbs$shdr.para.west.id[match(thetas.era[[i]]$BP.wg, era.all.pbs$BP.wg)]
#   thetas.era[[i]]$shdr.allo.id<- era.all.pbs$shdr.allo.id[match(thetas.era[[i]]$BP.wg, era.all.pbs$BP.wg)]
#   thetas.era[[i]]$sharing.type<- era.all.pbs$sharing.type[match(thetas.era[[i]]$BP.wg, era.all.pbs$BP.wg)]
#   thetas.era[[i]]$cline <- substr(listpaths.era[i], 38,45)
#   thetas.era[[i]]$pop.short <- substr(listpaths.era[i], 47,49)
#   thetas.era[[i]]$pop.full <- substr(listpaths.era[i], 38,49)}
# 
# for (i in 1:length(thetas.mel)) {
#   thetas.mel[[i]]$chr <- as.numeric(as.character(substr(thetas.mel[[i]]$scaff, 6,7)))
#   thetas.mel[[i]] <- subset( thetas.mel[[i]],  thetas.mel[[i]]$chr!=0)
#   thetas.mel[[i]]$BP.wg <- mel.all.pbs$BP.wg[match(paste(thetas.mel[[i]]$scaff, thetas.mel[[i]]$WinCenter), paste(mel.all.pbs$scaff, mel.all.pbs$midPos))]
#   thetas.mel[[i]]$shdr.para.east.id <- mel.all.pbs$shdr.para.east.id[match(thetas.mel[[i]]$BP.wg, mel.all.pbs$BP.wg)]
#   thetas.mel[[i]]$shdr.para.west.id <- mel.all.pbs$shdr.para.west.id[match(thetas.mel[[i]]$BP.wg, mel.all.pbs$BP.wg)]
#   thetas.mel[[i]]$shdr.allo.id<- mel.all.pbs$shdr.allo.id[match(thetas.mel[[i]]$BP.wg, mel.all.pbs$BP.wg)]
#   thetas.mel[[i]]$sharing.type<- mel.all.pbs$sharing.type[match(thetas.mel[[i]]$BP.wg, mel.all.pbs$BP.wg)]
#   thetas.mel[[i]]$cline <- substr(listpaths.mel[i], 38,45)
#   thetas.mel[[i]]$pop.short <- substr(listpaths.mel[i], 47,49)
#   thetas.mel[[i]]$pop.full <- substr(listpaths.mel[i], 38,49)}
# 
# # cocatenate
# thetas.era.all <- rbindlist(thetas.era)
# thetas.era.all$alt.type <- pop.info$alt.type[match(thetas.era.all$pop.full, pop.info$type.an.pop)]
# 
# thetas.mel.all <- rbindlist(thetas.mel)
# thetas.mel.all$alt.type <- pop.info$alt.type[match(thetas.mel.all$pop.full, pop.info$type.an.pop)]
# 
# thetas.era.all$background <- if_else(thetas.era.all$sharing.type=="", "background", "shared.region")
# thetas.mel.all$background <- if_else(thetas.mel.all$sharing.type=="", "background", "shared.region")
# 
# # save
# write.csv(thetas.era.all, "local/data/thetas/thetas.era.all.csv", row.names = F)
# write.csv(thetas.mel.all, "local/data/thetas/thetas.mel.all.csv", row.names = F)
# 

### load thetas sorted #####
thetas.era.all <- read.csv("local/data/thetas/thetas.era.all.csv")

thetas.era.all$background <- if_else(thetas.era.all$sharing.type=="", "background", "shared.region")

## test plots
head(thetas.era.all)
# &alt.type=="hig" &cline=="era.ec.e"
ggplot(subset(thetas.era.all, chr==1 &alt.type=="hig"&cline=="era.ec.e"), aes(y=Tajima, x=BP.wg))+
  geom_rect(inherit.aes = F, data=subset(era.shdr.para.west, !(shdr.para.west.id %in% era.shdr.allo.summ$shdr.para.west)) , aes(xmin=start, xmax=end, ymin=-2.5,ymax=1), colour="transparent", fill="#0372B2", alpha=0.4 )+
  geom_rect(inherit.aes = F, data=subset(era.shdr.para.east, !(shdr.para.east.id %in% era.shdr.allo.summ$shdr.para.east)), aes(xmin=start, xmax=end, ymin=-2.5,ymax=1), colour="transparent", fill="#049E73", alpha=0.4 )+
  geom_rect(inherit.aes = F, data=era.shdr.allo.all, aes(xmin=start, xmax=end, ymin=-2.5,ymax=1), colour="transparent", fill="#D65D00", alpha=0.6 )+
  geom_point(inherit.aes = F, data=subset(thetas.era.all, chr==1 &alt.type=="hig"&cline=="era.ec.e" & Tajima < quantile(Tajima , 0.05)), colour="orange",aes(y=Tajima, x=BP.wg))+
  geom_point(size=.1)+
  #geom_rect(inherit.aes = F, data=subset(era.shared.regions, sharing.type=="within.east"), aes(xmin=start.BP.wg-100000, xmax=end.BP.wg+100000, ymin=0,ymax=7), colour="transparent", fill="#049E73", alpha=0.7 )+
  #geom_rect(inherit.aes = F, data=subset(era.shared.regions, sharing.type=="across"), aes(xmin=start.BP.wg-100000, xmax=end.BP.wg+100000, ymin=0,ymax=7), colour="transparent", fill="#049E73", alpha=1 )+
  facet_wrap(~cline+alt.type+pop.short, ncol = 1)+
  xlim(10000000,20000000)+theme_bw()



ggplot(aes(x=WinCenter, y=Tajima, colour=alt.type, group=pop.full, shape=sharing.type), 
       data=subset(thetas.era.all, !(is.na(shdr.para.east.id!="")) & 
                     (shdr.para.east.id=="shdr.east.001"|shdr.para.east.id=="shdr.east.002"|shdr.para.east.id=="shdr.east.003"|
                        shdr.allo.id=="shdr.all.001"|shdr.allo.id=="shdr.all.004"|shdr.allo.id=="shdr.all.005")))+
  #geom_ma(lty=1,ma_fun = ZLEMA, n=10)+
  geom_point(alpha=0.9)+
  geom_point(inherit.aes = F, data= d[d$Point < quantile(d$Point, 0.95), ])
theme_bw()+
  facet_wrap(~shdr.para.east.id+sharing.type, scales = "free_x")

###  population means: tajima, pi #####
head(thetas[[1]])
taj <- tibble(tajima.mean=rep(0,length(listpaths) ),tp.mean=rep(0,length(listpaths) ), cline=rep(0,length(listpaths) ), pop=rep(0,length(listpaths) ))
for (i in 1:length(listpaths)) {
  taj$tajima.mean[i] <- mean(thetas[[i]]$Tajima)
  taj$tp.mean[i] <- mean(thetas[[i]]$tP/thetas[[i]]$nSites, na.rm=TRUE)
  taj$cline[i] <- substr(listpaths[i],38,45)
  taj$type.an.pop[i] <- substr(listpaths[i],38,49)
  taj$pop[i] <- substr(listpaths[i],47,49)}; taj

taj$sp <- substr(taj$cline,0, 3) ; taj
taj$alt.type <- pop.info$alt.type[match(taj$type.an.pop, pop.info$type.an.pop)]; taj



# #######  sort dxy ######
# 
# listpaths <- dir(pattern = ".txt", path = "local/data/local/data/output/dxy", full.names = T); listpaths
# dxy <- list()
# for (i in 1:length(listpaths)) {
#   dxy[[i]] <- read.table(listpaths[i], header = TRUE)
#   dxy[[i]]$cline <- substr(listpaths[i], 31,38)
#   dxy[[i]]$pop.short1 <- substr(listpaths[i], 27,29)
#   dxy[[i]]$pop.full1 <- substr(listpaths[i], 18,29)
#   dxy[[i]]$pop.short2 <- substr(listpaths[i], 40,42)
#   dxy[[i]]$pop.full2 <- substr(listpaths[i], 31,42)}
# 
# # divide by species
# dxy.era <- list(); dxy.era <- dxy[c(1:10)]
# dxy.mel <- list(); dxy.mel <- dxy[c(10:20)]
# 
# 
# ## add BP.wg, remove lines from unassigned chr in mel, add sharing ids
# head(dxy.era[[i]])
# 
# for (i in 1:length(dxy.era)) {
#   dxy.era[[i]]$chr <- as.numeric(as.character(substr(dxy.era[[i]]$chromo, 7,8)))
#   dxy.era[[i]]$BP.wg <- era.all.pbs$BP.wg[match(paste(dxy.era[[i]]$chromo, dxy.era[[i]]$win_mid), paste(era.all.pbs$scaff, era.all.pbs$midPos))]
#   dxy.era[[i]]$shdr.para.east.id <- era.all.pbs$shdr.para.east.id[match(dxy.era[[i]]$BP.wg, era.all.pbs$BP.wg)]
#   dxy.era[[i]]$shdr.para.west.id <- era.all.pbs$shdr.para.west.id[match(dxy.era[[i]]$BP.wg, era.all.pbs$BP.wg)]
#   dxy.era[[i]]$shdr.allo.id<- era.all.pbs$shdr.allo.id[match(dxy.era[[i]]$BP.wg, era.all.pbs$BP.wg)]
#   dxy.era[[i]]$sharing.type<- era.all.pbs$sharing.type[match(dxy.era[[i]]$BP.wg, era.all.pbs$BP.wg)]}
# 
# for (i in 1:length(dxy.mel)) {
#   dxy.mel[[i]]$chr <- as.numeric(as.character(substr(dxy.mel[[i]]$chromo, 6,7)))
#   dxy.mel[[i]] <- subset( dxy.mel[[i]],  dxy.mel[[i]]$chr!=0)
#   dxy.mel[[i]]$BP.wg <- mel.all.pbs$BP.wg[match(paste(dxy.mel[[i]]$chromo, dxy.mel[[i]]$win_mid), paste(mel.all.pbs$scaff, mel.all.pbs$midPos))]
#   dxy.mel[[i]]$shdr.para.east.id <- mel.all.pbs$shdr.para.east.id[match(dxy.mel[[i]]$BP.wg, mel.all.pbs$BP.wg)]
#   dxy.mel[[i]]$shdr.para.west.id <- mel.all.pbs$shdr.para.west.id[match(dxy.mel[[i]]$BP.wg, mel.all.pbs$BP.wg)]
#   dxy.mel[[i]]$shdr.allo.id<- mel.all.pbs$shdr.allo.id[match(dxy.mel[[i]]$BP.wg, mel.all.pbs$BP.wg)]
#   dxy.mel[[i]]$sharing.type<- mel.all.pbs$sharing.type[match(dxy.mel[[i]]$BP.wg, mel.all.pbs$BP.wg)]}
# 
# # concatenate
# dxy.era.all <- rbindlist(dxy.era); head(dxy.era.all)
# dxy.era.all$pop1.alt.type <- pop.info$alt.type[match(dxy.era.all$pop.full1, pop.info$type.an.pop)]
# dxy.era.all$pop2.alt.type <- pop.info$alt.type[match(dxy.era.all$pop.full2, pop.info$type.an.pop)]
# dxy.era.all$pop1.pop2.alt.type <- paste(dxy.era.all$pop1.alt.type, dxy.era.all$pop2.alt.type, sep = "_"); unique(dxy.era.all$pop1.pop2.alt.type)
# 
# dxy.mel.all <- rbindlist(dxy.mel); head(dxy.mel.all)
# dxy.mel.all$pop1.alt.type <- pop.info$alt.type[match(dxy.mel.all$pop.full1, pop.info$type.an.pop)]
# dxy.mel.all$pop2.alt.type <- pop.info$alt.type[match(dxy.mel.all$pop.full2, pop.info$type.an.pop)]
# dxy.mel.all$pop1.pop2.alt.type <- paste(dxy.mel.all$pop1.alt.type, dxy.mel.all$pop2.alt.type, sep = "_"); unique(dxy.mel.all$pop1.pop2.alt.type)
# 
# write.csv(dxy.era.all, "local/data/dxy/dxy.era.all.csv", row.names = F)
# write.csv(dxy.mel.all, "local/data/dxy/dxy.mel.all.csv", row.names = F)

### load dxy ####
dxy.era.all <- read.csv("local/data/dxy/dxy.era.all.csv")

######################################## 1. PLOT DENSITIES GENOME WIDE ##########
############## tajima #######
## era
thetas.era.all$country <- as.factor(substr(thetas.era.all$cline, 5,6))
thetas.era.all$side <- as.factor(substr(thetas.era.all$cline, 8,9))

# CAREFUL THIS FUCKS UP SIDES
levels(thetas.era.all$country) <- c("Colombia", "Ecuador")
levels(thetas.era.all$side) <- c("East", "West"); levels(thetas.era.all$side) <- c("West", "East" )

taji.dens.era <-  ggplot(subset(thetas.era.all), aes(x=Tajima,  colour=alt.type, shape=pop.full)) + 
  geom_density(fill="transparent", alpha=0, size=.5)+
  scale_color_manual(values=c("#0000ff", "#5ec55b","#23ff06"), labels=c("High", "Low", "Low distant"))+
  labs(colour="Altitude")+
  geom_vline(xintercept = 0, colour="lightgrey", lty="dashed")+
  xlab(expression(paste("Tajima's D",  italic(" H. erato"))))+
  theme_classic()+
  facet_wrap(~country+side)+
  theme(axis.line.x = element_line(color="black", size = 0.5),axis.line.y = element_line(color="black", size = 0.5),axis.title.y = element_text(size=14),
        plot.margin = unit(c(0.5,0.2,0.2,0.2), "lines"),axis.text = element_text(size=10), #axis.text.x = element_blank(),
        axis.text.y = element_text(size=10),axis.title.x = element_blank(), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_rect(fill = "transparent"), plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
        legend.background = element_rect(fill = "transparent"), legend.box.background = element_rect(fill = "transparent"), strip.background.y = element_blank(),
        strip.text = element_text(size=14))+
  scale_y_continuous(limits = c(0,3), expand = c(0,0))+
  scale_x_continuous(limits = c(-3,2.5), expand = c(0,0)); taji.dens.era

## mel
thetas.mel.all$country <- as.factor(substr(thetas.mel.all$cline, 5,6))
thetas.mel.all$side <- as.factor(substr(thetas.mel.all$cline, 8,9))

levels(thetas.mel.all$country) <- c("Colombia", "Ecuador")
levels(thetas.mel.all$side) <- c("East", "West"); levels(thetas.mel.all$side) <- c("West", "East" )

taji.dens.mel <-  ggplot(subset(thetas.mel.all), aes(x=Tajima,  colour=alt.type, shape=pop.full)) + 
  geom_density(fill="transparent", alpha=0, size=.5)+
  scale_color_manual(values=c("#0000ff", "#5ec55b","#23ff06"), labels=c("High", "Low", "Low distant"))+
  labs(colour="Altitude")+
  geom_vline(xintercept = 0, colour="lightgrey", lty="dashed")+
  xlab(expression(paste("Tajima's D",  italic(" H. melto"))))+
  theme_classic()+
  facet_wrap(~country+side)+
  theme(axis.line.x = element_line(color="black", size = 0.5),axis.line.y = element_line(color="black", size = 0.5),axis.title.y = element_text(size=14),
        plot.margin = unit(c(0.5,0.2,0.2,0.2), "lines"),axis.text = element_text(size=10), #axis.text.x = element_blank(),
        axis.text.y = element_text(size=10),axis.title.x = element_blank(), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_rect(fill = "transparent"), plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
        legend.background = element_rect(fill = "transparent"), legend.box.background = element_rect(fill = "transparent"), strip.background.y = element_blank(),
        strip.text = element_text(size=14))+
  scale_y_continuous(limits = c(0,2), expand = c(0,0))+
  scale_x_continuous(limits = c(-3,3), expand = c(0,0)); taji.dens.mel

## plot
plot_grid(taji.dens.era , taji.dens.mel, ncol = 1, labels = c("A", "B"), align = "v")
ggsave("~/Dropbox (Cambridge University)/PhD/22_pop.gen.paper/figures/supplementary/tajima.density.era.mel.png", height = 10, width=9)

############## pi #######
## era
pi.dens.era <-  ggplot(subset(thetas.era.all), aes(x=tP/nSites,  colour=alt.type, shape=pop.full)) + 
  geom_density(fill="transparent", alpha=0, size=.5)+
  scale_color_manual(values=c("#0000ff", "#5ec55b","#23ff06"), labels=c("High", "Low", "Low distant"))+
  labs(colour="Altitude")+
  #geom_vline(xintercept = 0, colour="lightgrey", lty="dashed")+
  xlab(expression(paste("Nucleotide diversity",  italic(" H. erato"))))+
  theme_classic()+
  facet_wrap(~country+side)+
  theme(axis.line.x = element_line(color="black", size = 0.5),axis.line.y = element_line(color="black", size = 0.5),axis.title.y = element_text(size=14),
        plot.margin = unit(c(0.5,0.2,0.2,0.2), "lines"),axis.text = element_text(size=10), #axis.text.x = element_blank(),
        axis.text.y = element_text(size=10),axis.title.x = element_blank(), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_rect(fill = "transparent"), plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
        legend.background = element_rect(fill = "transparent"), legend.box.background = element_rect(fill = "transparent"), strip.background.y = element_blank(),
        strip.text = element_text(size=14))+
  scale_y_continuous(limits = c(0,60), expand = c(0,0))+
  scale_x_continuous(limits = c(0,0.1), expand = c(0,0)); pi.dens.era

## mel
pi.dens.mel <-  ggplot(subset(thetas.mel.all), aes(x=tP/nSites,  colour=alt.type, shape=pop.full)) + 
  geom_density(fill="transparent", alpha=0, size=.5)+
  scale_color_manual(values=c("#0000ff", "#5ec55b","#23ff06"), labels=c("High", "Low", "Low distant"))+
  labs(colour="Altitude")+
  #geom_vline(xintercept = 0, colour="lightgrey", lty="dashed")+
  xlab(expression(paste("Nucleotide diversity",  italic(" H. melto"))))+
  theme_classic()+
  facet_wrap(~country+side)+
  theme(axis.line.x = element_line(color="black", size = 0.5),axis.line.y = element_line(color="black", size = 0.5),axis.title.y = element_text(size=14),
        plot.margin = unit(c(0.5,0.2,0.2,0.2), "lines"),axis.text = element_text(size=10), #axis.text.x = element_blank(),
        axis.text.y = element_text(size=10),axis.title.x = element_blank(), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_rect(fill = "transparent"), plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
        legend.background = element_rect(fill = "transparent"), legend.box.background = element_rect(fill = "transparent"), strip.background.y = element_blank(),
        strip.text = element_text(size=14))+
  scale_y_continuous(limits = c(0,60), expand = c(0,0))+
  scale_x_continuous(limits = c(0,0.1), expand = c(0,0)); pi.dens.mel


## plot
plot_grid(pi.dens.era , pi.dens.mel, ncol = 1, labels = c("A", "B"), align = "v")
ggsave("~/Dropbox (Cambridge University)/PhD/22_pop.gen.paper/figures/supplementary/pi.density.era.mel.png", height = 10, width=9)


############## dxy #######
ggplot(data=dxy.era.all, aes(colour=pop1.pop2.alt.type, x=dxy_mean))+
  geom_density()+facet_wrap(~substr(pop.full1, 0,6))

ggplot(data=dxy.mel.all, aes(colour=pop1.pop2.alt.type, x=dxy_mean))+
  geom_density()+facet_wrap(~substr(pop.full1, 0,6))
## era
dxy.dens.era <-  ggplot(subset(dxy.era.all), aes(x=dxy_mean,  colour=pop1.pop2.alt.type)) + 
  geom_density(fill="transparent", alpha=0, size=1.5)+
  scale_color_manual(values=c("#0000ff", "#5ec55b","#23ff06"), labels=c("High-low", "Hig-vLow", "low-vlow"))+
  labs(colour="Altitude")+
  geom_vline(xintercept = 0, colour="lightgrey", lty="dashed")+
  xlab(expression(paste("dxy's D",  italic(" H. erato"))))+
  theme_classic()+
  facet_wrap(~substr(pop.full1, 0,9))+
  theme(axis.line.x = element_line(color="black", size = 0.5),axis.line.y = element_line(color="black", size = 0.5),axis.title.y = element_text(size=14),
        plot.margin = unit(c(0.5,0.2,0.2,0.2), "lines"),axis.text = element_text(size=10), #axis.text.x = element_blank(),
        axis.text.y = element_text(size=10),axis.title.x = element_blank(), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_rect(fill = "transparent"), plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
        legend.background = element_rect(fill = "transparent"), legend.box.background = element_rect(fill = "transparent"), strip.background.y = element_blank(),
        strip.text = element_text(size=14))+
  scale_y_continuous(expand = c(0,0))+
  scale_x_continuous(expand = c(0,0)); dxy.dens.era

## mel
dxy.dens.mel <-  ggplot(subset(dxy.mel.all), aes(x=dxy_mean,  colour=pop1.pop2.alt.type)) + 
  geom_density(fill="transparent", alpha=0, size=1)+
  scale_color_manual(values=c("#0000ff", "#5ec55b","#23ff06"), labels=c("High-low", "Hig-vLow", "low-vlow"))+
  labs(colour="Altitude")+
  geom_vline(xintercept = 0, colour="lightgrey", lty="dashed")+
  xlab(expression(paste("dxy's D",  italic(" H. melto"))))+
  theme_classic()+
  facet_wrap(~substr(pop.full1, 0,9))+
  theme(axis.line.x = element_line(color="black", size = 0.5),axis.line.y = element_line(color="black", size = 0.5),axis.title.y = element_text(size=14),
        plot.margin = unit(c(0.5,0.2,0.2,0.2), "lines"),axis.text = element_text(size=10), #axis.text.x = element_blank(),
        axis.text.y = element_text(size=10),axis.title.x = element_blank(), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_rect(fill = "transparent"), plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
        legend.background = element_rect(fill = "transparent"), legend.box.background = element_rect(fill = "transparent"), strip.background.y = element_blank(),
        strip.text = element_text(size=14))+
  scale_y_continuous(expand = c(0,0))+
  scale_x_continuous(expand = c(0,0)); dxy.dens.mel

## plot
plot_grid(dxy.dens.era , dxy.dens.mel, ncol = 1, labels = c("A", "B"), align = "v")
ggsave("~/Dropbox (Cambridge University)/PhD/22_pop.gen.paper/figures/supplementary/dxy.density.era.mel.png", height = 10, width=9)


######################################## 2. TAJIMA PER SHDR ##########
################# 2.1 prep data ##########
# # only use thetas from 1 highland pop per , subset to use less memory
# # are there 5 windows in each hdr that are in the 5th percentile? yes no
# list.high.pops <- subset(pop.info, alt.type=="hig" )$type.an.pop; list.high.pops
# pops.to.use <- read.csv("local/data/local/data/output/dxy/pop.to.compare.csv")
# 
# thetas.era <- subset(thetas.era.all,pop.full %in% pops.to.use$pop1 & pop.full %in% list.high.pops )
# thetas.mel <- subset(thetas.mel.all,pop.full %in% hig.pops.to.use$pop1 & pop.full %in% list.high.pops)
# 
# # 4 high alt pops in each species - check
# unique(subset(thetas.era.all,pop.full %in% pops.to.use$pop1 & pop.full %in% list.high.pops )$pop.full)
# unique(subset(thetas.mel.all,pop.full %in% pops.to.use$pop1 & pop.full %in% list.high.pops )$pop.full)
# 
# # save
# write.csv(thetas.era, "local/data/thetas/thetas.era.hig.pop.csv", row.names = F)
# write.csv(thetas.mel, "local/data/thetas/thetas.mel.hig.pop.csv", row.names = F)

## load hig alt thetas only 
thetas.era.hig <- read.csv("local/data/thetas/thetas.era.hig.pop.csv")
thetas.era.hig$side <- as.factor(substr(thetas.era.hig$cline, 8,9))

################# 2.1 SHDR east ######
str(thetas.era.hig)
# lowest in ecuador or colombia
## get summaries

# calculate genomewide quantile threshold for each pop (one high alt per country)
theta.5th.quant.east.col <- quantile(subset(thetas.era.hig, side=="e"&country=="Colombia")$Tajima, 0.05)
theta.5th.quant.east.ecu <- quantile(subset(thetas.era.hig, side=="e"&country=="Ecuador")$Tajima, 0.05)

# assess in each cline tajima values per hdr
thetas.era.hig.para.east.col.summ <- summarise(group_by(subset(thetas.era.hig, side=="e"&country=="Colombia"&shdr.para.east.id!="" ),shdr.para.east.id, pop.full , country ),
                                               start=min(BP.wg),
                                               end=max(BP.wg),
                                               size.hdr.bp=max(BP.wg)- min(BP.wg),
                                               no.windows=n(), scaff=unique(scaff),
                                               tajima.mean=mean(Tajima), tajima.min=min(Tajima),
                                               tajima.5thperc.no.windows= length(Tajima[Tajima < theta.5th.quant.east.col ]),
                                               tajima.5thperc.mean.taj= mean(Tajima[Tajima < theta.5th.quant.east.col]) ); thetas.era.hig.para.east.col.summ 

thetas.era.hig.para.east.ecu.summ <- summarise(group_by(subset(thetas.era.hig, side=="e"&country=="Ecuador"&shdr.para.east.id!="" ),shdr.para.east.id,pop.full , country ),
                                               start=min(BP.wg),
                                               end=max(BP.wg),
                                               size.hdr.bp=max(BP.wg)- min(BP.wg),
                                               no.windows=n(), scaff=unique(scaff),tajima.mean=mean(Tajima), tajima.min=min(Tajima),
                                               tajima.5thperc.no.windows= length(Tajima[Tajima < theta.5th.quant.east.ecu ]),
                                               tajima.5thperc.mean.taj= mean(Tajima[Tajima < theta.5th.quant.east.ecu]) ); thetas.era.hig.para.east.ecu.summ 

thetas.era.hig.para.east.summ <- rbind(thetas.era.hig.para.east.col.summ, thetas.era.hig.para.east.ecu.summ)
thetas.era.hig.para.east.summ$tajima.5thperc.perc.windows <- (thetas.era.hig.para.east.summ$tajima.5thperc.no.windows / thetas.era.hig.para.east.summ$no.windows)*100

# obtain mean tajima across clines
thetas.era.hig.para.east.summ.summ <- summarise(group_by(thetas.era.hig.para.east.summ ,shdr.para.east.id,start, end, size.hdr.bp,no.windows ),
                                                scaff=unique(scaff),
                                                tajima.mean.mean=mean(tajima.mean),
                                                tajima.min.mean=mean(tajima.min),
                                                tajima.5thperc.no.windows.mean=mean(tajima.5thperc.no.windows),
                                                tajima.5thperc.mean.taj.mean=mean(tajima.5thperc.mean.taj),
                                                tajima.5thperc.perc.windows.mean=mean(tajima.5thperc.perc.windows )); thetas.era.hig.para.east.summ.summ

# which of these para SHDR also allo SHDR?
thetas.era.hig.para.east.summ.summ$is.allo.shdr  <- if_else(thetas.era.hig.para.east.summ.summ$shdr.para.east.id %in% era.shdr.allo.summ$shdr.para.east, "yes", "no")

###### sims ######
# # extract the same number of shdr para east and of the same width, random locations get tajima min, tajima.5thperc.no.windows, tajima.5thperc.mean.taj
# # first create intervals with real shdr
# era.shdr.para.east.int <- reduce(Intervals(as.matrix(era.shdr.para.east[,c("start","end")]), closed=c(TRUE,TRUE), type="R"))
# 
# # random intervals
# sim1k.era.east <- lapply(1:1000, function(x){rand_non_overlapping_intervals(era.genome_size,size(era.shdr.para.east.int))})
# 
# sim.thetas.era.hig.para.east.summ.summ <-list()
# for (i in 1:1000) {
#   # prep random intervals for action
#   sim1k.era.df <- as.data.frame(sim1k.era.east[i])
#   sim1k.era.df$ran.hdr.name <- paste("ran", i, ".hdr.", seq(1:nrow(sim1k.era.df )), sep = "")
#   
#   # all thetas tmp , use only east pops here
#   theta.era.tmp <- thetas.era.hig[,c("BP.wg", "Tajima", "pop.full","side", "country")]
#   theta.era.tmp <- subset(theta.era.tmp, side=="e")
#   
#   # add ran hdr names to theta.tmp
#   theta.era.tmp$ran.hdr.name <- sim1k.era.df$ran.hdr.name[mapply(match.by.range.1k, theta.era.tmp$BP.wg)]
#   
#   # use observed genome-wide quantile thresholds
#   # get min td in col and ecu
#   # subset by hdr being present to avoid getting an NA category
#   sim.thetas.era.hig.para.east.col.summ <- summarise(group_by(subset(theta.era.tmp , country=="Colombia"&ran.hdr.name!=""), ran.hdr.name, pop.full  ),
#                                                  no.windows=n(),
#                                                  tajima.mean=mean(Tajima),tajima.min=min(Tajima),
#                                                  tajima.5thperc.no.windows= length(Tajima[Tajima < theta.5th.quant.east.col ]),
#                                                  tajima.5thperc.mean.taj= mean(Tajima[Tajima < theta.5th.quant.east.col ]),
#                                                  tajima.5thperc.perc.windows=(tajima.5thperc.no.windows/ no.windows)*100 )
#   
#   sim.thetas.era.hig.para.east.ecu.summ <- summarise(group_by(subset(theta.era.tmp , country=="Ecuador"&ran.hdr.name!=""), ran.hdr.name, pop.full  ),
#                                                      no.windows=n(),
#                                                      tajima.mean=mean(Tajima),tajima.min=min(Tajima),
#                                                      tajima.5thperc.no.windows= length(Tajima[Tajima < theta.5th.quant.east.ecu ]),
#                                                      tajima.5thperc.mean.taj= mean(Tajima[Tajima < theta.5th.quant.east.ecu ]),
#                                                      tajima.5thperc.perc.windows=(tajima.5thperc.no.windows/ no.windows)*100 )
#   sim.thetas.era.hig.para.east.summ <- rbind( sim.thetas.era.hig.para.east.col.summ,  sim.thetas.era.hig.para.east.ecu.summ)
#   
#   # get summ summ
#   sim.thetas.era.hig.para.east.summ.summ[[i]] <- summarise(group_by(sim.thetas.era.hig.para.east.summ  , ran.hdr.name),
#                                                   tajima.mean.mean=mean(tajima.mean),
#                                                   tajima.min.mean=mean(tajima.min),
#                                                   tajima.5thperc.no.windows.mean=mean(tajima.5thperc.no.windows),
#                                                   tajima.5thperc.mean.taj.mean=mean(tajima.5thperc.mean.taj),
#                                                   tajima.5thperc.perc.windows.mean= mean(tajima.5thperc.perc.windows) )
# }
# 
# sim.thetas.era.hig.para.east.summ.summ.df <- bind_rows(sim.thetas.era.hig.para.east.summ.summ, .id = "column_label")
# write.csv(sim.thetas.era.hig.para.east.summ.summ.df, "local/data/thetas/sim.1k.thetas.era.hig.para.east.summ.summ.df.csv")

### load them 
sim.thetas.era.hig.para.east.summ.summ.df <- read.csv("local/data/thetas/sim.1k.thetas.era.hig.para.east.summ.summ.df.csv")


###### plot #####
library(ggbeeswarm)
ggplot(data=thetas.era.hig.para.east.summ.summ, aes(y=tajima.5thperc.mean.taj.mean,x= is.allo.shdr))+
  geom_boxplot()+ theme_classic() + geom_beeswarm(color="grey")+
  geom_boxplot(inherit.aes = F, data=sim.thetas.era.hig.para.east.summ.summ.df, aes(y=tajima.5thperc.mean.taj.mean, x=3))+
  geom_hline(yintercept = mean(subset(thetas.era.hig, side=="e")$Tajima), lty="dashed", color="darkgrey")

ggplot(data=thetas.era.hig.para.east.summ.summ, aes(y=tajima.5thperc.perc.windows.mean ,x= is.allo.shdr, fill=is.allo.shdr))+
  geom_boxplot()+ theme_classic() + geom_beeswarm(color="grey")+ ylim(0,60)+
  scale_fill_manual(values=c("#049E73", "#D65D00"))+
  geom_violin(inherit.aes = F, data=sim.thetas.era.hig.para.east.summ.summ.df, aes(y=tajima.5thperc.perc.windows.mean, x=3))+
  geom_boxplot(inherit.aes = F, data=sim.thetas.era.hig.para.east.summ.summ.df, aes(y=tajima.5thperc.perc.windows.mean, x=3), fill="grey", alpha=.5)+
  geom_hline(yintercept = quantile(sim.thetas.era.hig.para.east.summ.summ.df$tajima.5thperc.perc.windows.mean, 0.95))
#geom_hline(yintercept = mean(subset(thetas.era.hig, side=="e")$Tajima), lty="dashed", color="darkgrey")



# plot min tajima's d per HDR vs simlated min per HDR
taji.dens.era.hig.e <-  ggplot(subset(thetas.era.hig, side=="e"), aes(x=Tajima,  colour=alt.type, shape=pop.full)) + 
  scale_color_manual(values=c("#0000ff", "#5ec55b","#23ff06"), labels=c("High", "Low", "Low distant"))+
  labs(colour="Altitude")+
  geom_vline(xintercept = 0, colour="lightgrey", lty="dashed")+
  geom_vline(xintercept = subset(thetas.era.hig.para.east.summ.summ, is.allo.shdr=="no")$tajima.min.mean, colour="#049E73", lty="solid", alpha=.9, size=0.2)+
  geom_vline(xintercept = subset(thetas.era.hig.para.east.summ.summ, is.allo.shdr=="yes")$tajima.min.mean, colour="#D65D00", lty="solid", alpha=.9, size=0.2)+
  geom_vline(xintercept = mean(subset(sim.thetas.era.hig.para.east.summ.summ.df)$tajima.min.mean), colour="black", lty="solid", alpha=.9, size=1.5)+
  geom_vline(xintercept = quantile(subset(sim.thetas.era.hig.para.east.summ.summ.df)$tajima.min.mean, 0.9), colour="black", lty="dashed", alpha=.9, size=1.5)+
  geom_vline(xintercept = quantile(subset(sim.thetas.era.hig.para.east.summ.summ.df)$tajima.min.mean, 0.1), colour="black", lty="dashed", alpha=.9, size=1.5)+
  geom_density(fill="transparent", alpha=0, size=1)+
    xlab(expression(paste("Tajima's D",  italic(" H. erato"))))+
  theme_classic()+
  theme(axis.line.x = element_line(color="black", size = 0.5),axis.line.y = element_line(color="black", size = 0.5),axis.title.y = element_text(size=14),
        plot.margin = unit(c(0.5,0.2,0.2,0.2), "lines"),axis.text = element_text(size=10), #axis.text.x = element_blank(),
        axis.text.y = element_text(size=10),axis.title.x = element_blank(), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_rect(fill = "transparent"), plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
        legend.background = element_rect(fill = "transparent"), legend.box.background = element_rect(fill = "transparent"), strip.background.y = element_blank(),
        strip.text = element_text(size=14))+
  scale_y_continuous(limits = c(0,3), expand = c(0,0))+
  scale_x_continuous(limits = c(-3,2.5), expand = c(0,0)); taji.dens.era.hig.e

ggplot(data=thetas.era.hig.para.east.summ.summ, aes(y=tajima.5thperc.perc.windows.mean ,x= is.allo.shdr,  colour=is.allo.shdr))+
  geom_beeswarm(cex = 2, size=2)+ ylim(0,60)+
  scale_colour_manual(values=c("#049E73", "#D65D00", labels=c("Within east\n(parapatric)", "Across sides\n(allopatric)")))+
  theme_classic() + 
  geom_rect(fill="grey80",color=NA, alpha=0.02, aes(xmin=.7,xmax=1.3, ymin=quantile(sim.thetas.era.hig.para.east.summ.summ.df$tajima.5thperc.perc.windows.mean,0.05),
                                                    ymax=quantile(sim.thetas.era.hig.para.east.summ.summ.df$tajima.5thperc.perc.windows.mean, 0.9)))+
  geom_rect(fill="grey80",color=NA, alpha=0.02, aes(xmin=1.7,xmax=2.3, ymin=quantile(sim.thetas.era.hig.para.east.summ.summ.df$tajima.5thperc.perc.windows.mean,0.05),
                                                    ymax=quantile(sim.thetas.era.hig.para.east.summ.summ.df$tajima.5thperc.perc.windows.mean, 0.9)))+
  xlab("SHDR type") + scale_x_discrete(labels=c("no" = "Within east\n(parapatric)", "yes"="Across sides\n(allopatric)"))+
  ylab("% Tajima outlier windows within SHDR")+
  theme(legend.position = "none", axis.text = element_text(size=14), axis.title = element_text(size=14))


ggplot(data=thetas.era.hig.para.east.summ.summ, aes(x=tajima.5thperc.perc.windows.mean, y=tajima.5thperc.mean.taj.mean))+ geom_point()+
  stat_cor()+stat_smooth(method = "lm")+theme_bw()+ xlab("% of HDR windows <5th quatile")+ ylab("Mean Tajima of windows <5th quantile")
# very few have % above 0
ggplot(data=sim.thetas.era.hig.para.east.summ.summ.df, aes(x=tajima.5thperc.perc.windows.mean, y=tajima.5thperc.mean.taj.mean))+ geom_point()+
  stat_cor()+stat_smooth(method = "lm")+theme_bw()+ xlab("% of HDR windows <5th quatile")+ ylab("Mean Tajima of windows <5th quantile")



quantile(sim.thetas.era.hig.para.east.summ.summ.df$tajima.5thperc.perc.windows.mean, 0.95)


################# 2.2 SHDR west ######
str(thetas.era.hig)
# lowest in ecuador or colombia
## get summaries

# calculate genomewide quantile threshold for each pop (one high alt per country)
theta.5th.quant.west.col <- quantile(subset(thetas.era.hig, side=="w"&country=="Colombia")$Tajima, 0.05)
theta.5th.quant.west.ecu <- quantile(subset(thetas.era.hig, side=="w"&country=="Ecuador")$Tajima, 0.05)

thetas.era.hig.para.west.col.summ <- summarise(group_by(subset(thetas.era.hig, side=="w"&country=="Colombia"&shdr.para.west.id!="" ),shdr.para.west.id, pop.full , country ),
                                               start=min(BP.wg),
                                               end=max(BP.wg),
                                               size.hdr.bp=max(BP.wg)- min(BP.wg),
                                               no.windows=n(), scaff=unique(scaff),
                                               tajima.mean=mean(Tajima), tajima.min=min(Tajima),
                                               tajima.5thperc.no.windows= length(Tajima[Tajima < theta.5th.quant.west.col ]),
                                               tajima.5thperc.mean.taj= mean(Tajima[Tajima < theta.5th.quant.west.col]) ); thetas.era.hig.para.west.col.summ 

thetas.era.hig.para.west.ecu.summ <- summarise(group_by(subset(thetas.era.hig, side=="w"&country=="Ecuador"&shdr.para.west.id!="" ),shdr.para.west.id,pop.full , country ),
                                               start=min(BP.wg),
                                               end=max(BP.wg),
                                               size.hdr.bp=max(BP.wg)- min(BP.wg),
                                               no.windows=n(), scaff=unique(scaff),tajima.mean=mean(Tajima), tajima.min=min(Tajima),
                                               tajima.5thperc.no.windows= length(Tajima[Tajima < theta.5th.quant.west.ecu ]),
                                               tajima.5thperc.mean.taj= mean(Tajima[Tajima < theta.5th.quant.west.ecu]) ); thetas.era.hig.para.west.ecu.summ 

thetas.era.hig.para.west.summ <- rbind(thetas.era.hig.para.west.col.summ, thetas.era.hig.para.west.ecu.summ)
thetas.era.hig.para.west.summ$tajima.5thperc.perc.windows <- (thetas.era.hig.para.west.summ$tajima.5thperc.no.windows / thetas.era.hig.para.west.summ$no.windows)*100

# summarise to get mean lowest
thetas.era.hig.para.west.summ.summ <- summarise(group_by(thetas.era.hig.para.west.summ ,shdr.para.west.id,start, end, size.hdr.bp,no.windows ),
                                                scaff=unique(scaff),
                                                tajima.mean.mean=mean(tajima.mean),
                                                tajima.min.mean=mean(tajima.min),
                                                tajima.5thperc.no.windows.mean=mean(tajima.5thperc.no.windows),
                                                tajima.5thperc.mean.taj.mean=mean(tajima.5thperc.mean.taj),
                                                tajima.5thperc.perc.windows.mean=mean(tajima.5thperc.perc.windows )); thetas.era.hig.para.west.summ.summ

# which of these para SHDR also allo SHDR?
thetas.era.hig.para.west.summ.summ$is.allo.shdr  <- if_else(thetas.era.hig.para.west.summ.summ$shdr.para.west.id %in% era.shdr.allo.summ$shdr.para.west, "yes", "no")

###### sims ######
# # extract the same number of shdr para west and of the same width, random locations get tajima min, tajima.5thperc.no.windows, tajima.5thperc.mean.taj
# # first create intervals with real shdr
# era.shdr.para.west.int <- reduce(Intervals(as.matrix(era.shdr.para.west[,c("start","end")]), closed=c(TRUE,TRUE), type="R"))
# 
# # random intervals
# sim1k.era.west <- lapply(1:1000, function(x){rand_non_overlapping_intervals(era.genome_size,size(era.shdr.para.west.int))})
# 
# sim.thetas.era.hig.para.west.summ.summ <-list()
# for (i in 1:1000) {
#   # prep random intervals for action
#   sim1k.era.df <- as.data.frame(sim1k.era.west[i])
#   sim1k.era.df$ran.hdr.name <- paste("ran", i, ".hdr.", seq(1:nrow(sim1k.era.df )), sep = "")
#   
#   # all thetas tmp , use only west pops here
#   theta.era.tmp <- thetas.era.hig[,c("BP.wg", "Tajima", "pop.full","side", "country")]
#   theta.era.tmp <- subset(theta.era.tmp, side=="w")
#   
#   # add ran hdr names to theta.tmp
#   theta.era.tmp$ran.hdr.name <- sim1k.era.df$ran.hdr.name[mapply(match.by.range.1k, theta.era.tmp$BP.wg)]
#   
#   # use observed genome-wide quantile thresholds
#   # get min td in col and ecu
#   # subset by hdr being present to avoid getting an NA category
#   sim.thetas.era.hig.para.west.col.summ <- summarise(group_by(subset(theta.era.tmp , country=="Colombia"&ran.hdr.name!=""), ran.hdr.name, pop.full  ),
#                                                      no.windows=n(),
#                                                      tajima.mean=mean(Tajima),tajima.min=min(Tajima),
#                                                      tajima.5thperc.no.windows= length(Tajima[Tajima < theta.5th.quant.west.col ]),
#                                                      tajima.5thperc.mean.taj= mean(Tajima[Tajima < theta.5th.quant.west.col ]),
#                                                      tajima.5thperc.perc.windows=(tajima.5thperc.no.windows/ no.windows)*100 )
#   
#   sim.thetas.era.hig.para.west.ecu.summ <- summarise(group_by(subset(theta.era.tmp , country=="Ecuador"&ran.hdr.name!=""), ran.hdr.name, pop.full  ),
#                                                      no.windows=n(),
#                                                      tajima.mean=mean(Tajima),tajima.min=min(Tajima),
#                                                      tajima.5thperc.no.windows= length(Tajima[Tajima < theta.5th.quant.west.ecu ]),
#                                                      tajima.5thperc.mean.taj= mean(Tajima[Tajima < theta.5th.quant.west.ecu ]),
#                                                      tajima.5thperc.perc.windows=(tajima.5thperc.no.windows/ no.windows)*100 )
#   sim.thetas.era.hig.para.west.summ <- rbind( sim.thetas.era.hig.para.west.col.summ,  sim.thetas.era.hig.para.west.ecu.summ)
#   
#   # get summ summ
#   sim.thetas.era.hig.para.west.summ.summ[[i]] <- summarise(group_by(sim.thetas.era.hig.para.west.summ  , ran.hdr.name),
#                                                            tajima.mean.mean=mean(tajima.mean),
#                                                            tajima.min.mean=mean(tajima.min),
#                                                            tajima.5thperc.no.windows.mean=mean(tajima.5thperc.no.windows),
#                                                            tajima.5thperc.mean.taj.mean=mean(tajima.5thperc.mean.taj),
#                                                            tajima.5thperc.perc.windows.mean= mean(tajima.5thperc.perc.windows) )
# }
# 
# sim.thetas.era.hig.para.west.summ.summ.df <- bind_rows(sim.thetas.era.hig.para.west.summ.summ, .id = "column_label")
# write.csv(sim.thetas.era.hig.para.west.summ.summ.df, "local/data/thetas/sim.1k.thetas.era.hig.para.west.summ.summ.df.csv")

### load them 
sim.thetas.era.hig.para.west.summ.summ.df <- read.csv("local/data/thetas/sim.1k.thetas.era.hig.para.west.summ.summ.df.csv")


###### plot #####
library(ggbeeswarm)
ggplot(data=thetas.era.hig.para.west.summ.summ, aes(y=tajima.5thperc.mean.taj.mean,x= is.allo.shdr))+
  geom_boxplot()+ theme_classic() + geom_beeswarm(color="grey")+
  geom_boxplot(inherit.aes = F, data=sim.thetas.era.hig.para.west.summ.summ.df, aes(y=tajima.5thperc.mean.taj.mean, x=3))+
  geom_hline(yintercept = mean(subset(thetas.era.hig, side=="w")$Tajima), lty="dashed", color="darkgrey")

ggplot(data=thetas.era.hig.para.west.summ.summ, aes(y=tajima.5thperc.perc.windows.mean ,x= is.allo.shdr,  colour=is.allo.shdr))+
  geom_beeswarm(cex = 2, size=2)+ ylim(0,60)+
  scale_colour_manual(values=c("#0372B2", "#D65D00", labels=c("Within west\n(parapatric)", "Across sides\n(allopatric)")))+
  theme_classic() + 
  geom_rect(fill="grey80",color=NA, alpha=0.02, aes(xmin=.7,xmax=1.3, ymin=quantile(sim.thetas.era.hig.para.west.summ.summ.df$tajima.5thperc.perc.windows.mean,0.05),
                        ymax=quantile(sim.thetas.era.hig.para.west.summ.summ.df$tajima.5thperc.perc.windows.mean, 0.9)))+
  geom_rect(fill="grey80",color=NA, alpha=0.02, aes(xmin=1.7,xmax=2.3, ymin=quantile(sim.thetas.era.hig.para.west.summ.summ.df$tajima.5thperc.perc.windows.mean,0.05),
                                                    ymax=quantile(sim.thetas.era.hig.para.west.summ.summ.df$tajima.5thperc.perc.windows.mean, 0.9)))+
  xlab("SHDR type") + scale_x_discrete(labels=c("no" = "Within west\n(parapatric)", "yes"="Across sides\n(allopatric)"))+
  theme(legend.position = "none")



# plot min tajima's d per HDR vs simlated min per HDR
taji.dens.era.hig.w <-  ggplot(subset(thetas.era.hig, side=="w"), aes(x=Tajima,  colour=alt.type, shape=pop.full)) + 
  scale_color_manual(values=c("#0000ff", "#5ec55b","#23ff06"), labels=c("High", "Low", "Low distant"))+
  labs(colour="Altitude")+
  geom_vline(xintercept = 0, colour="lightgrey", lty="dashed")+
  geom_vline(xintercept = subset(thetas.era.hig.para.west.summ.summ, is.allo.shdr=="no")$tajima.min.mean, colour="#0372B2", lty="solid", alpha=.9, size=0.2)+
  geom_vline(xintercept = subset(thetas.era.hig.para.west.summ.summ, is.allo.shdr=="yes")$tajima.min.mean, colour="#D65D00", lty="solid", alpha=.9, size=0.2)+
  geom_vline(xintercept = mean(subset(sim.thetas.era.hig.para.west.summ.summ.df)$tajima.min.mean), colour="black", lty="solid", alpha=.9, size=1.5)+
  geom_vline(xintercept = quantile(subset(sim.thetas.era.hig.para.west.summ.summ.df)$tajima.min.mean, 0.9), colour="black", lty="dashed", alpha=.9, size=1.5)+
  geom_vline(xintercept = quantile(subset(sim.thetas.era.hig.para.west.summ.summ.df)$tajima.min.mean, 0.1), colour="black", lty="dashed", alpha=.9, size=1.5)+
  geom_density(fill="transparent", alpha=0, size=1)+
  xlab(expression(paste("Tajima's D",  italic(" H. erato"))))+
  theme_classic()+
  theme(axis.line.x = element_line(color="black", size = 0.5),axis.line.y = element_line(color="black", size = 0.5),axis.title.y = element_text(size=14),
        plot.margin = unit(c(0.5,0.2,0.2,0.2), "lines"),axis.text = element_text(size=10), #axis.text.x = element_blank(),
        axis.text.y = element_text(size=10),axis.title.x = element_blank(), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_rect(fill = "transparent"), plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
        legend.background = element_rect(fill = "transparent"), legend.box.background = element_rect(fill = "transparent"), strip.background.y = element_blank(),
        strip.text = element_text(size=14))+
  scale_y_continuous(limits = c(0,3), expand = c(0,0))+
  scale_x_continuous(limits = c(-3,2.5), expand = c(0,0)); taji.dens.era.hig.w


ggplot(data=thetas.era.hig.para.west.summ.summ, aes(x=tajima.5thperc.perc.windows.mean, y=tajima.5thperc.mean.taj.mean))+ geom_point()+
  stat_cor()+stat_smooth(method = "lm")+theme_bw()+ xlab("% of HDR windows <5th quatile")+ ylab("Mean Tajima of windows <5th quantile")
# very few have % above 0
ggplot(data=sim.thetas.era.hig.para.west.summ.summ.df, aes(x=tajima.5thperc.perc.windows.mean, y=tajima.5thperc.mean.taj.mean))+ geom_point()+
  stat_cor()+stat_smooth(method = "lm")+theme_bw()+ xlab("% of HDR windows <5th quatile")+ ylab("Mean Tajima of windows <5th quantile")




################# 2.3 which SHDR are tajima outliers? ######
library(patchwork)
taji.dens.era.hig.e  / taji.dens.era.hig.w 
ggsave("local/plots/thetas/tajima.era.min.per.HDR.90-10thpercentile.dens.east.west.png", width = 5, height = 5)


# check thresholds for min tajima per HDR and % windows with tajima <5th perc
nrow(subset(thetas.era.hig.para.east.summ.summ, tajima.min.mean < quantile(subset(sim.thetas.era.hig.para.east.summ.summ.df)$tajima.min.mean, 0.1)))
nrow(subset(thetas.era.hig.para.east.summ.summ, tajima.5thperc.perc.windows.mean > quantile(subset(sim.thetas.era.hig.para.east.summ.summ.df)$tajima.5thperc.perc.windows.mean, 0.9)))

nrow(subset(thetas.era.hig.para.east.summ.summ, tajima.min.mean < quantile(subset(sim.thetas.era.hig.para.east.summ.summ.df)$tajima.min.mean, 0.05)))
nrow(subset(thetas.era.hig.para.east.summ.summ, tajima.5thperc.perc.windows.mean > quantile(subset(sim.thetas.era.hig.para.east.summ.summ.df)$tajima.5thperc.perc.windows.mean, 0.95)))

nrow(subset(thetas.era.hig.para.west.summ.summ, tajima.min.mean < quantile(subset(sim.thetas.era.hig.para.west.summ.summ.df)$tajima.min.mean, 0.1)))
nrow(subset(thetas.era.hig.para.west.summ.summ, tajima.5thperc.perc.windows.mean > quantile(subset(sim.thetas.era.hig.para.west.summ.summ.df)$tajima.5thperc.perc.windows.mean, 0.9)))

nrow(subset(thetas.era.hig.para.west.summ.summ, tajima.min.mean < quantile(subset(sim.thetas.era.hig.para.west.summ.summ.df)$tajima.min.mean, 0.05)))
nrow(subset(thetas.era.hig.para.west.summ.summ, tajima.5thperc.perc.windows.mean > quantile(subset(sim.thetas.era.hig.para.west.summ.summ.df)$tajima.5thperc.perc.windows.mean, 0.95)))

## add whether SHDRs are outliers or not
thetas.era.hig.para.east.summ.summ$tajima.min.less.5th.perc.sims <- if_else(thetas.era.hig.para.east.summ.summ$tajima.min.mean < quantile(sim.thetas.era.hig.para.east.summ.summ.df$tajima.min.mean, 0.05), "yes", "no")
thetas.era.hig.para.east.summ.summ$tajima.min.less.10th.perc.sims <- if_else(thetas.era.hig.para.east.summ.summ$tajima.min.mean < quantile(sim.thetas.era.hig.para.east.summ.summ.df$tajima.min.mean, 0.1), "yes", "no")
thetas.era.hig.para.east.summ.summ$tajima.5thperc.perc.windows.mean.less.95th.perc.sims <- if_else(thetas.era.hig.para.east.summ.summ$tajima.5thperc.perc.windows.mean > quantile(sim.thetas.era.hig.para.east.summ.summ.df$tajima.5thperc.perc.windows.mean, 0.95), "yes", "no")
thetas.era.hig.para.east.summ.summ$tajima.5thperc.perc.windows.mean.less.90th.perc.sims <- if_else(thetas.era.hig.para.east.summ.summ$tajima.5thperc.perc.windows.mean > quantile(sim.thetas.era.hig.para.east.summ.summ.df$tajima.5thperc.perc.windows.mean, 0.9), "yes", "no")

thetas.era.hig.para.west.summ.summ$tajima.min.less.5th.perc.sims <- if_else(thetas.era.hig.para.west.summ.summ$tajima.min.mean < quantile(sim.thetas.era.hig.para.west.summ.summ.df$tajima.min.mean, 0.05), "yes", "no")
thetas.era.hig.para.west.summ.summ$tajima.min.less.10th.perc.sims <- if_else(thetas.era.hig.para.west.summ.summ$tajima.min.mean < quantile(sim.thetas.era.hig.para.west.summ.summ.df$tajima.min.mean, 0.1), "yes", "no")
thetas.era.hig.para.west.summ.summ$tajima.5thperc.perc.windows.mean.less.95th.perc.sims <- if_else(thetas.era.hig.para.west.summ.summ$tajima.5thperc.perc.windows.mean > quantile(sim.thetas.era.hig.para.west.summ.summ.df$tajima.5thperc.perc.windows.mean, 0.95), "yes", "no")
thetas.era.hig.para.west.summ.summ$tajima.5thperc.perc.windows.mean.less.90th.perc.sims <- if_else(thetas.era.hig.para.west.summ.summ$tajima.5thperc.perc.windows.mean > quantile(sim.thetas.era.hig.para.west.summ.summ.df$tajima.5thperc.perc.windows.mean, 0.9), "yes", "no"); thetas.era.hig.para.west.summ.summ

# save
thetas.era.hig.para.east.summ.summ; write.csv(thetas.era.hig.para.east.summ.summ, "local/data/thetas/shdr.thetas.era.hig.para.east.summ.summ.csv")
thetas.era.hig.para.west.summ.summ; write.csv(thetas.era.hig.para.west.summ.summ, "local/data/thetas/shdr.thetas.era.hig.para.west.summ.summ.csv")


######################################## 3. DELTAPI PER SHDR ##########
################# 3.1 prep data ##########
# # only use thetas from 1 highland pop per , subset to use less memory
# # are there 5 windows in each hdr that are in the 5th percentile? yes no
# list.high.pops <- subset(pop.info, alt.type=="hig" )$type.an.pop; list.high.pops
# list.vlo.pops <- subset(pop.info, alt.type=="vlo" )$type.an.pop; list.vlo.pops
# 
# # same pops as dxy for delta pi
# pops.to.use <- read.csv("local/data/local/data/output/dxy/pop.to.compare.csv")
# 
# thetas.era <- subset(thetas.era.all, (pop.full %in% pops.to.use$pop1 & pop.full %in% list.high.pops) | (pop.full %in% pops.to.use$pop2 & pop.full %in% list.vlo.pops) )
# thetas.mel <- subset(thetas.mel.all,(pop.full %in% pops.to.use$pop1 & pop.full %in% list.high.pops) | (pop.full %in% pops.to.use$pop2 & pop.full %in% list.vlo.pops)  )
# 
# # 4 high alt and 4 vlo pops in each species - check
# unique(subset(thetas.era.all, (pop.full %in% pops.to.use$pop1 & pop.full %in% list.high.pops) | (pop.full %in% pops.to.use$pop2 & pop.full %in% list.vlo.pops)  )$pop.full)
# unique(subset(thetas.mel.all,(pop.full %in% pops.to.use$pop1 & pop.full %in% list.high.pops) | (pop.full %in% pops.to.use$pop2 & pop.full %in% list.vlo.pops) )$pop.full)
# 
# # save
# write.csv(thetas.era, "local/data/thetas/thetas.era.hig.vlo.pop.csv", row.names = F)
# write.csv(thetas.mel, "local/data/thetas/thetas.mel.hig.vlo.pop.csv", row.names = F)

#### read prepped data ####
thetas.era <- read.csv( "local/data/thetas/thetas.era.hig.vlo.pop.csv")

### summarise (so mean pi of pops) by alt.type per pos, then spread by altitude, then plot deltapi
pi.era.all.alt <- summarise(group_by(thetas.era, chr, WinCenter,BP.wg, shdr.para.east.id, shdr.para.west.id, alt.type, cline,  background),
          #Tajima=mean(Tajima),
          pi=mean(tP/nSites)); head(pi.era.all.alt)

# spread
pi.era.all.alt.wide <-spread(pi.era.all.alt, alt.type, pi)
pi.era.all.alt.wide$side <- substr(pi.era.all.alt.wide$cline, 8,9)

# calculate delta pi (if not subsetting hig low, can do the others)
#pi.era.all.alt.wide$delta.pi.hig.low <- pi.era.all.alt.wide$hig-pi.era.all.alt.wide$low
#pi.era.all.alt.wide$delta.pi.low.vlo <- pi.era.all.alt.wide$low-pi.era.all.alt.wide$vlo
pi.era.all.alt.wide$delta.pi.hig.vlo <- pi.era.all.alt.wide$hig-pi.era.all.alt.wide$vlo
head(pi.era.all.alt.wide)

###### plots test ####
ggplot(subset(pi.era.all.alt.wide, chr==1 ), aes(y=delta.pi.hig.vlo, x=BP.wg))+
  geom_rect(inherit.aes = F, data=subset(era.shdr.para.west, !(shdr.para.west.id %in% era.shdr.allo.summ$shdr.para.west)) , aes(xmin=start, xmax=end, ymin=-.01,ymax=.01), colour="transparent", fill="#0372B2", alpha=0.4 )+
  geom_rect(inherit.aes = F, data=subset(era.shdr.para.east, !(shdr.para.east.id %in% era.shdr.allo.summ$shdr.para.east)), aes(xmin=start, xmax=end,  ymin=-.01,ymax=.01), colour="transparent", fill="#049E73", alpha=0.4 )+
  geom_rect(inherit.aes = F, data=era.shdr.allo.all, aes(xmin=start, xmax=end, ymin=-.01,ymax=.01), colour="transparent", fill="#D65D00", alpha=0.6 )+
  geom_point(inherit.aes = F, data=subset(pi.era.all.alt.wide, chr==1 & cline=="era.ec.e" & delta.pi.hig.vlo < quantile(delta.pi.hig.vlo , 0.05, na.rm = T)), colour="orange",aes(y=delta.pi.hig.vlo, x=BP.wg))+
  geom_point(size=.1)+
  scale_y_continuous( limits=c(-.025,.015) , expand = c(0,0))+
  #geom_rect(inherit.aes = F, data=subset(era.shared.regions, sharing.type=="within.east"), aes(xmin=start.BP.wg-100000, xmax=end.BP.wg+100000, ymin=0,ymax=7), colour="transparent", fill="#049E73", alpha=0.7 )+
  #geom_rect(inherit.aes = F, data=subset(era.shared.regions, sharing.type=="across"), aes(xmin=start.BP.wg-100000, xmax=end.BP.wg+100000, ymin=0,ymax=7), colour="transparent", fill="#049E73", alpha=1 )+
  facet_wrap(~cline, ncol = 1)+
  xlim(10000000,20000000)+theme_bw()

ggplot(subset(pi.era.all.alt.wide), aes(y=delta.pi.hig.vlo, x=BP.wg))+
  geom_rect(inherit.aes = F, data=subset(era.shdr.para.west, !(shdr.para.west.id %in% era.shdr.allo.summ$shdr.para.west)) , aes(xmin=start, xmax=end, ymin=-.01,ymax=.01), colour="transparent", fill="#0372B2", alpha=0.4 )+
  geom_rect(inherit.aes = F, data=subset(era.shdr.para.east, !(shdr.para.east.id %in% era.shdr.allo.summ$shdr.para.east)), aes(xmin=start, xmax=end,  ymin=-.01,ymax=.01), colour="transparent", fill="#049E73", alpha=0.4 )+
  geom_rect(inherit.aes = F, data=era.shdr.allo.all, aes(xmin=start, xmax=end, ymin=-.01,ymax=.01), colour="transparent", fill="#D65D00", alpha=0.6 )+
  geom_point(inherit.aes = F, data=subset(pi.era.all.alt.wide,  cline=="era.ec.e" & delta.pi.hig.vlo < quantile(delta.pi.hig.vlo , 0.05, na.rm = T)), colour="orange",aes(y=delta.pi.hig.vlo, x=BP.wg))+
  geom_point(size=.1)+
  scale_y_continuous( limits=c(-.025,.015) , expand = c(0,0))+
  #geom_rect(inherit.aes = F, data=subset(era.shared.regions, sharing.type=="within.east"), aes(xmin=start.BP.wg-100000, xmax=end.BP.wg+100000, ymin=0,ymax=7), colour="transparent", fill="#049E73", alpha=0.7 )+
  #geom_rect(inherit.aes = F, data=subset(era.shared.regions, sharing.type=="across"), aes(xmin=start.BP.wg-100000, xmax=end.BP.wg+100000, ymin=0,ymax=7), colour="transparent", fill="#049E73", alpha=1 )+
  facet_wrap(~cline, ncol = 1)+
  xlim(25000000,27000000)+theme_bw()

ggplot(subset(pi.era.all.alt.wide, chr==2), aes(y=delta.pi.hig.vlo, x=BP.wg))+
  geom_point(inherit.aes = F, data=subset(pi.era.all.alt.wide, chr==2 &cline=="era.ec.e" & delta.pi.hig.vlo < quantile(delta.pi.hig.vlo , 0.05, na.rm = T)), colour="orange",aes(y=delta.pi.hig.vlo, x=BP.wg))+
  geom_point(size=.1)+scale_y_continuous( limits=c(-.025,.015) , expand = c(0,0))+facet_wrap(~cline, ncol = 1)+theme_bw()

ggplot(subset(pi.era.all.alt.wide, chr==2), aes(y=delta.pi.hig.low, x=BP.wg))+
  geom_point(inherit.aes = F, data=subset(pi.era.all.alt.wide, chr==2 &cline=="era.ec.e" & delta.pi.hig.low < quantile(delta.pi.hig.low , 0.05, na.rm = T)), colour="orange",aes(y=delta.pi.hig.low, x=BP.wg))+
  geom_point(size=.1)+
  scale_y_continuous( limits=c(-.025,.015) , expand = c(0,0))+
  facet_wrap(~cline, ncol = 1)+theme_bw()

ggplot(subset(pi.era.all.alt.wide), aes(y=delta.pi.low.vlo, x=BP.wg, colour=cline))+
  geom_rect(inherit.aes = F, data=subset(era.shdr.para.west, !(shdr.para.west.id %in% era.shdr.allo.summ$shdr.para.west)) , aes(xmin=start, xmax=end, ymin=-.01,ymax=.01), colour="transparent", fill="#0372B2", alpha=0.4 )+
  geom_rect(inherit.aes = F, data=subset(era.shdr.para.east, !(shdr.para.east.id %in% era.shdr.allo.summ$shdr.para.east)), aes(xmin=start, xmax=end,  ymin=-.01,ymax=.01), colour="transparent", fill="#049E73", alpha=0.4 )+
  geom_rect(inherit.aes = F, data=era.shdr.allo.all, aes(xmin=start, xmax=end, ymin=-.01,ymax=.01), colour="transparent", fill="#D65D00", alpha=0.6 )+
  geom_point(inherit.aes = F, data=subset(pi.era.all.alt.wide,  cline=="era.ec.e" & delta.pi.low.vlo < quantile(delta.pi.low.vlo , 0.05, na.rm = T)), colour="orange",aes(y=delta.pi.low.vlo, x=BP.wg))+
  geom_point(size=.1)+
  scale_y_continuous( limits=c(-.015,.015) , expand = c(0,0))+
  #geom_rect(inherit.aes = F, data=subset(era.shared.regions, sharing.type=="within.east"), aes(xmin=start.BP.wg-100000, xmax=end.BP.wg+100000, ymin=0,ymax=7), colour="transparent", fill="#049E73", alpha=0.7 )+
  #geom_rect(inherit.aes = F, data=subset(era.shared.regions, sharing.type=="across"), aes(xmin=start.BP.wg-100000, xmax=end.BP.wg+100000, ymin=0,ymax=7), colour="transparent", fill="#049E73", alpha=1 )+
  facet_wrap(~cline, ncol = 1)+
  xlim(10000000,20000000)+theme_bw()

ggplot(subset(pi.era.all.alt.wide, chr==1 &cline=="era.ec.e"), aes(y=delta.pi.low.vlo, x=WinCenter))+
  geom_rect(inherit.aes = F, data=subset(era.shdr.para.west, !(shdr.para.west.id %in% era.shdr.allo.summ$shdr.para.west)) , aes(xmin=start, xmax=end, ymin=-.01,ymax=.01), colour="transparent", fill="#0372B2", alpha=0.4 )+
  geom_rect(inherit.aes = F, data=subset(era.shdr.para.east, !(shdr.para.east.id %in% era.shdr.allo.summ$shdr.para.east)), aes(xmin=start, xmax=end,  ymin=-.01,ymax=.01), colour="transparent", fill="#049E73", alpha=0.4 )+
  geom_rect(inherit.aes = F, data=era.shdr.allo.all, aes(xmin=start, xmax=end, ymin=-.01,ymax=.01), colour="transparent", fill="#D65D00", alpha=0.6 )+
  geom_point(inherit.aes = F, data=subset(pi.era.all.alt.wide, chr==1 & cline=="era.ec.e" & delta.pi.low.vlo < quantile(delta.pi.low.vlo , 0.05, na.rm = T)), colour="orange",aes(y=delta.pi.low.vlo, x=WinCenter))+
  geom_point(size=.1)+
  #geom_rect(inherit.aes = F, data=subset(era.shared.regions, sharing.type=="within.east"), aes(xmin=start.BP.wg-100000, xmax=end.BP.wg+100000, ymin=0,ymax=7), colour="transparent", fill="#049E73", alpha=0.7 )+
  #geom_rect(inherit.aes = F, data=subset(era.shared.regions, sharing.type=="across"), aes(xmin=start.BP.wg-100000, xmax=end.BP.wg+100000, ymin=0,ymax=7), colour="transparent", fill="#049E73", alpha=1 )+
  facet_wrap(~cline, ncol = 1)+
  xlim(10000000,20000000)+theme_bw()


################# 3.2 SHDR east ######
head(pi.era.all.alt.wide)

# calculate genomewide quantile threshold for each pop (one high alt per country)
delta.pi.hig.vlo.5th.quant.east.col <- quantile(subset(pi.era.all.alt.wide, cline=="era.co.e")$delta.pi.hig.vlo, 0.05, na.rm = T )
delta.pi.hig.vlo.5th.quant.east.ecu <- quantile(subset(pi.era.all.alt.wide, cline=="era.ec.e")$delta.pi.hig.vlo, 0.05, na.rm = T)

head(pi.era.all.alt.wide)
pi.era.all.alt.wide.para.east.col.summ <- summarise(group_by(subset(pi.era.all.alt.wide, cline=="era.co.e" &shdr.para.east.id!="" ), shdr.para.east.id, cline ),
                                               start=min(BP.wg),
                                               end=max(BP.wg),
                                               size.hdr.bp=max(BP.wg)- min(BP.wg),
                                               no.windows=n(), chr=unique(chr),
                                               delta.pi.hig.vlo.mean=mean(delta.pi.hig.vlo), delta.pi.hig.vlo.min=min(delta.pi.hig.vlo),
                                               delta.pi.hig.vlo.5thperc.no.windows= length(delta.pi.hig.vlo[delta.pi.hig.vlo  < delta.pi.hig.vlo.5th.quant.east.col ]),
                                               delta.pi.hig.vlo.5thperc.mean.delta.pi= mean(delta.pi.hig.vlo[delta.pi.hig.vlo < delta.pi.hig.vlo.5th.quant.east.col]) ); pi.era.all.alt.wide.para.east.col.summ 

pi.era.all.alt.wide.para.east.ecu.summ <- summarise(group_by(subset(pi.era.all.alt.wide, cline=="era.ec.e" &shdr.para.east.id!="" ),shdr.para.east.id, cline ),
                                                    start=min(BP.wg),
                                                    end=max(BP.wg),
                                                    size.hdr.bp=max(BP.wg)- min(BP.wg),
                                                    no.windows=n(), chr=unique(chr),
                                                    delta.pi.hig.vlo.mean=mean(delta.pi.hig.vlo), delta.pi.hig.vlo.min=min(delta.pi.hig.vlo),
                                                    delta.pi.hig.vlo.5thperc.no.windows= length(delta.pi.hig.vlo[delta.pi.hig.vlo  < delta.pi.hig.vlo.5th.quant.east.ecu ]),
                                                    delta.pi.hig.vlo.5thperc.mean.delta.pi= mean(delta.pi.hig.vlo[delta.pi.hig.vlo < delta.pi.hig.vlo.5th.quant.east.ecu ]) ); pi.era.all.alt.wide.para.east.ecu.summ 

pi.era.all.alt.wide.para.east.summ <- rbind(pi.era.all.alt.wide.para.east.col.summ, pi.era.all.alt.wide.para.east.ecu.summ)
pi.era.all.alt.wide.para.east.summ$delta.pi.hig.vlo.5thperc.perc.windows <- (pi.era.all.alt.wide.para.east.summ$delta.pi.hig.vlo.5thperc.no.windows / pi.era.all.alt.wide.para.east.summ$no.windows)*100

# summarise to get mean loeast
pi.era.all.alt.wide.para.east.summ.summ <- summarise(group_by(pi.era.all.alt.wide.para.east.summ ,shdr.para.east.id,start, end, size.hdr.bp, no.windows ),
                                                     chr=unique(chr),
                                                delta.pi.hig.vlo.mean.mean=mean(delta.pi.hig.vlo.mean),
                                                delta.pi.hig.vlo.min.mean=mean(delta.pi.hig.vlo.min),
                                                delta.pi.hig.vlo.5thperc.no.windows.mean=mean(delta.pi.hig.vlo.5thperc.no.windows),
                                                delta.pi.hig.vlo.5thperc.mean.delta.pi.mean=mean(delta.pi.hig.vlo.5thperc.mean.delta.pi),
                                                delta.pi.hig.vlo.5thperc.perc.windows.mean=mean(delta.pi.hig.vlo.5thperc.perc.windows )); pi.era.all.alt.wide.para.east.summ.summ

# which of these para SHDR also allo SHDR?
pi.era.all.alt.wide.para.east.summ.summ$is.allo.shdr  <- if_else(pi.era.all.alt.wide.para.east.summ.summ$shdr.para.east.id %in% era.shdr.allo.summ$shdr.para.east, "yes", "no")

###### sims ######
# # extract the same number of shdr para east and of the same width, random locations get delta.pi min, delta.pi.5thperc.no.windows, delta.pi.5thperc.mean.taj
# # first create intervals with real shdr
# era.shdr.para.east.int <- reduce(Intervals(as.matrix(era.shdr.para.east[,c("start","end")]), closed=c(TRUE,TRUE), type="R"))
# 
# # random intervals
# sim1k.era.east <- lapply(1:1000, function(x){rand_non_overlapping_intervals(era.genome_size,size(era.shdr.para.east.int))})
# 
# sim.pi.era.all.alt.wide.para.east.summ.summ <-list()
# for (i in 1:1000) {
#   # prep random intervals for action
#   sim1k.era.df <- as.data.frame(sim1k.era.east[i])
#   sim1k.era.df$ran.hdr.name <- paste("ran", i, ".hdr.", seq(1:nrow(sim1k.era.df )), sep = "")
#   
#   # all pi tmp , use only east pops here
#   delta.pi.era.tmp <- pi.era.all.alt.wide[,c("BP.wg", "delta.pi.hig.vlo", "cline")]
#   delta.pi.era.tmp$side <- substr(delta.pi.era.tmp$cline, 8,9)
#   delta.pi.era.tmp <- subset(delta.pi.era.tmp, side=="e")
#   
#   # add ran hdr names to delta.pitmp
#   delta.pi.era.tmp$ran.hdr.name <- sim1k.era.df$ran.hdr.name[mapply(match.by.range.1k, delta.pi.era.tmp$BP.wg)]
#   
#   # use observed genome-wide quantile thresholds
#   # get min td in col and ecu
#   # subset by hdr being present to avoid getting an NA category
#   sim.pi.era.all.alt.wide.para.east.col.summ <- summarise(group_by(subset(delta.pi.era.tmp, cline=="era.co.e" &ran.hdr.name!=""),ran.hdr.name, cline ),
#                                                                     no.windows=n(), 
#                                                                     delta.pi.hig.vlo.mean=mean(delta.pi.hig.vlo), delta.pi.hig.vlo.min=min(delta.pi.hig.vlo),
#                                                                     delta.pi.hig.vlo.5thperc.no.windows= length(delta.pi.hig.vlo[delta.pi.hig.vlo  < delta.pi.hig.vlo.5th.quant.east.col ]),
#                                                                     delta.pi.hig.vlo.5thperc.mean.delta.pi= mean(delta.pi.hig.vlo[delta.pi.hig.vlo < delta.pi.hig.vlo.5th.quant.east.col]),
#                                                           delta.pi.hig.vlo.5thperc.perc.windows=(delta.pi.hig.vlo.5thperc.no.windows/ no.windows)*100 )
#   sim.pi.era.all.alt.wide.para.east.ecu.summ <- summarise(group_by(subset(delta.pi.era.tmp, cline=="era.ec.e" &ran.hdr.name!=""),ran.hdr.name, cline ),
#                                                           no.windows=n(), 
#                                                           delta.pi.hig.vlo.mean=mean(delta.pi.hig.vlo), delta.pi.hig.vlo.min=min(delta.pi.hig.vlo),
#                                                           delta.pi.hig.vlo.5thperc.no.windows= length(delta.pi.hig.vlo[delta.pi.hig.vlo  < delta.pi.hig.vlo.5th.quant.east.ecu ]),
#                                                           delta.pi.hig.vlo.5thperc.mean.delta.pi= mean(delta.pi.hig.vlo[delta.pi.hig.vlo < delta.pi.hig.vlo.5th.quant.east.ecu]),
#                                                           delta.pi.hig.vlo.5thperc.perc.windows=(delta.pi.hig.vlo.5thperc.no.windows/ no.windows)*100 )  
# 
#   sim.pi.era.all.alt.wide.para.east.summ <- rbind( sim.pi.era.all.alt.wide.para.east.col.summ,  sim.pi.era.all.alt.wide.para.east.ecu.summ)
#   
#   # get summ summ
#   sim.pi.era.all.alt.wide.para.east.summ.summ[[i]] <- summarise(group_by(sim.pi.era.all.alt.wide.para.east.summ  , ran.hdr.name),
#                                                            delta.pi.hig.vlo.mean.mean=mean(delta.pi.hig.vlo.mean),
#                                                            delta.pi.hig.vlo.min.mean=mean(delta.pi.hig.vlo.min),
#                                                            delta.pi.hig.vlo.5thperc.no.windows.mean=mean(delta.pi.hig.vlo.5thperc.no.windows),
#                                                            delta.pi.hig.vlo.5thperc.mean.delta.pi.mean=mean(delta.pi.hig.vlo.5thperc.mean.delta.pi),
#                                                            delta.pi.hig.vlo.5thperc.perc.windows.mean= mean(delta.pi.hig.vlo.5thperc.perc.windows) )
# }
# 
# sim.delta.pi.hig.vlo.era.para.east.summ.summ.df <- bind_rows(sim.pi.era.all.alt.wide.para.east.summ.summ, .id = "column_label")
# write.csv(sim.delta.pi.hig.vlo.era.para.east.summ.summ.df, "local/data/thetas/sim.delta.pi.hig.vlo.era.para.east.summ.summ.df.csv")

#### read prepped data ####
sim.delta.pi.hig.vlo.era.para.east.summ.summ.df <- read.csv( "local/data/thetas/sim.delta.pi.hig.vlo.era.para.east.summ.summ.df.csv")


###### plot #####
library(ggbeeswarm)
ggplot(data=pi.era.all.alt.wide.para.east.summ.summ, aes(y=delta.pi.hig.vlo.min.mean,x= is.allo.shdr))+
  geom_boxplot()+ theme_classic() + geom_beeswarm(color="grey")+
  geom_boxplot(inherit.aes = F, data=sim.delta.pi.hig.vlo.era.para.east.summ.summ.df, aes(y=delta.pi.hig.vlo.min.mean, x=3))+
  geom_hline(yintercept = mean(subset(pi.era.all.alt.wide, side=="e")$delta.pi.hig.vlo), lty="dashed", color="darkgrey")

ggplot(data=pi.era.all.alt.wide.para.east.summ.summ, aes(y=delta.pi.hig.vlo.5thperc.perc.windows.mean ,x= is.allo.shdr, fill=is.allo.shdr))+
  geom_boxplot()+ theme_classic() + geom_beeswarm(color="grey")+ ylim(0,60)+
  scale_fill_manual(values=c("#049E73", "#D65D00"))+
  geom_violin(inherit.aes = F, data=sim.delta.pi.hig.vlo.era.para.east.summ.summ.df, aes(y=delta.pi.hig.vlo.5thperc.perc.windows.mean, x=3))+
  geom_boxplot(inherit.aes = F, data=sim.delta.pi.hig.vlo.era.para.east.summ.summ.df, aes(y=delta.pi.hig.vlo.5thperc.perc.windows.mean, x=3), fill="grey", alpha=.5)+
  geom_hline(yintercept = quantile(sim.delta.pi.hig.vlo.era.para.east.summ.summ.df$delta.pi.hig.vlo.5thperc.perc.windows.mean, 0.95))
#geom_hline(yintercept = mean(subset(pi.era.all.alt.wide, side=="e")$Tajima), lty="dashed", color="darkgrey")



# plot min delta.pi's d per HDR vs simlated min per HDR
delta.pi.hig.vlo.dens.era.hig.e <-  ggplot(subset(pi.era.all.alt.wide, side=="e"), aes(x=delta.pi.hig.vlo, shape=cline)) + 
  scale_color_manual(values=c("#0000ff", "#5ec55b","#23ff06"), labels=c("High", "Low", "Low distant"))+
  labs(colour="Altitude")+
  geom_vline(xintercept = 0, colour="lightgrey", lty="dashed")+
  geom_vline(xintercept = subset(pi.era.all.alt.wide.para.east.summ.summ, is.allo.shdr=="no")$delta.pi.hig.vlo.min.mean, colour="#049E73", lty="solid", alpha=.9, size=0.2)+
  geom_vline(xintercept = subset(pi.era.all.alt.wide.para.east.summ.summ, is.allo.shdr=="yes")$delta.pi.hig.vlo.min.mean, colour="#D65D00", lty="solid", alpha=.9, size=0.2)+
  geom_vline(xintercept = mean(subset(sim.delta.pi.hig.vlo.era.para.east.summ.summ.df)$delta.pi.hig.vlo.min.mean,na.rm=T), colour="black", lty="solid", alpha=.9, size=1.5)+
  geom_vline(xintercept = quantile(subset(sim.delta.pi.hig.vlo.era.para.east.summ.summ.df)$delta.pi.hig.vlo.min.mean, 0.9,na.rm=T), colour="black", lty="dashed", alpha=.9, size=1.5)+
  geom_vline(xintercept = quantile(subset(sim.delta.pi.hig.vlo.era.para.east.summ.summ.df)$delta.pi.hig.vlo.min.mean, 0.1,na.rm=T), colour="black", lty="dashed", alpha=.9, size=1.5)+
  geom_density(fill="transparent", alpha=0, size=1, colour="blue")+
  xlab(expression(paste("Tajima's D",  italic(" H. erato"))))+
  theme_classic()+
  theme(axis.line.x = element_line(color="black", size = 0.5),axis.line.y = element_line(color="black", size = 0.5),axis.title.y = element_text(size=14),
        plot.margin = unit(c(0.5,0.2,0.2,0.2), "lines"),axis.text = element_text(size=10), #axis.text.x = element_blank(),
        axis.text.y = element_text(size=10),axis.title.x = element_blank(), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_rect(fill = "transparent"), plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
        legend.background = element_rect(fill = "transparent"), legend.box.background = element_rect(fill = "transparent"), strip.background.y = element_blank(),
        strip.text = element_text(size=14))+
  scale_y_continuous(expand = c(0,0))+
  scale_x_continuous(expand = c(0,0)); delta.pi.hig.vlo.dens.era.hig.e

ggplot(data=pi.era.all.alt.wide.para.east.summ.summ, aes(y=delta.pi.hig.vlo.5thperc.perc.windows.mean ,x= is.allo.shdr,  colour=is.allo.shdr))+
  geom_beeswarm(cex = 2, size=2)+ ylim(0,60)+
  scale_colour_manual(values=c("#049E73", "#D65D00", labels=c("Within east\n(parapatric)", "Across sides\n(allopatric)")))+
  theme_classic() + 
  geom_rect(fill="grey80",color=NA, alpha=0.02, aes(xmin=.7,xmax=1.3, ymin=quantile(sim.delta.pi.hig.vlo.era.para.east.summ.summ.df$delta.pi.hig.vlo.5thperc.perc.windows.mean,0.05),
                                                    ymax=quantile(sim.delta.pi.hig.vlo.era.para.east.summ.summ.df$delta.pi.hig.vlo.5thperc.perc.windows.mean, 0.9)))+
  geom_rect(fill="grey80",color=NA, alpha=0.02, aes(xmin=1.7,xmax=2.3, ymin=quantile(sim.delta.pi.hig.vlo.era.para.east.summ.summ.df$delta.pi.hig.vlo.5thperc.perc.windows.mean,0.05),
                                                    ymax=quantile(sim.delta.pi.hig.vlo.era.para.east.summ.summ.df$delta.pi.hig.vlo.5thperc.perc.windows.mean, 0.9)))+
  xlab("SHDR type") + scale_x_discrete(labels=c("no" = "Within east\n(parapatric)", "yes"="Across sides\n(allopatric)"))+
  ylab("% Tajima outlier windows within SHDR")+
  theme(legend.position = "none", axis.text = element_text(size=14), axis.title = element_text(size=14))



################# 3.3 SHDR west ######
head(pi.era.all.alt.wide)

# calculate genomewide quantile threshold for each pop (one high alt per country)
delta.pi.hig.vlo.5th.quant.west.col <- quantile(subset(pi.era.all.alt.wide, cline=="era.co.w")$delta.pi.hig.vlo, 0.05, na.rm = T )
delta.pi.hig.vlo.5th.quant.west.ecu <- quantile(subset(pi.era.all.alt.wide, cline=="era.ec.w")$delta.pi.hig.vlo, 0.05, na.rm = T)

head(pi.era.all.alt.wide)
pi.era.all.alt.wide.para.west.col.summ <- summarise(group_by(subset(pi.era.all.alt.wide, cline=="era.co.w" &shdr.para.west.id!="" ), shdr.para.west.id, cline ),
                                                    start=min(BP.wg),
                                                    end=max(BP.wg),
                                                    size.hdr.bp=max(BP.wg)- min(BP.wg),
                                                    no.windows=n(), chr=unique(chr),
                                                    delta.pi.hig.vlo.mean=mean(delta.pi.hig.vlo), delta.pi.hig.vlo.min=min(delta.pi.hig.vlo),
                                                    delta.pi.hig.vlo.5thperc.no.windows= length(delta.pi.hig.vlo[delta.pi.hig.vlo  < delta.pi.hig.vlo.5th.quant.west.col ]),
                                                    delta.pi.hig.vlo.5thperc.mean.delta.pi= mean(delta.pi.hig.vlo[delta.pi.hig.vlo < delta.pi.hig.vlo.5th.quant.west.col]) ); pi.era.all.alt.wide.para.west.col.summ 

pi.era.all.alt.wide.para.west.ecu.summ <- summarise(group_by(subset(pi.era.all.alt.wide, cline=="era.ec.w" &shdr.para.west.id!="" ),shdr.para.west.id, cline ),
                                                    start=min(BP.wg),
                                                    end=max(BP.wg),
                                                    size.hdr.bp=max(BP.wg)- min(BP.wg),
                                                    no.windows=n(), chr=unique(chr),
                                                    delta.pi.hig.vlo.mean=mean(delta.pi.hig.vlo), delta.pi.hig.vlo.min=min(delta.pi.hig.vlo),
                                                    delta.pi.hig.vlo.5thperc.no.windows= length(delta.pi.hig.vlo[delta.pi.hig.vlo  < delta.pi.hig.vlo.5th.quant.west.ecu ]),
                                                    delta.pi.hig.vlo.5thperc.mean.delta.pi= mean(delta.pi.hig.vlo[delta.pi.hig.vlo < delta.pi.hig.vlo.5th.quant.west.ecu ]) ); pi.era.all.alt.wide.para.west.ecu.summ 

pi.era.all.alt.wide.para.west.summ <- rbind(pi.era.all.alt.wide.para.west.col.summ, pi.era.all.alt.wide.para.west.ecu.summ)
pi.era.all.alt.wide.para.west.summ$delta.pi.hig.vlo.5thperc.perc.windows <- (pi.era.all.alt.wide.para.west.summ$delta.pi.hig.vlo.5thperc.no.windows / pi.era.all.alt.wide.para.west.summ$no.windows)*100

# summarise to get mean lowest
pi.era.all.alt.wide.para.west.summ.summ <- summarise(group_by(pi.era.all.alt.wide.para.west.summ ,shdr.para.west.id,start, end, size.hdr.bp, no.windows ),
                                                     chr=unique(chr),
                                                     delta.pi.hig.vlo.mean.mean=mean(delta.pi.hig.vlo.mean),
                                                     delta.pi.hig.vlo.min.mean=mean(delta.pi.hig.vlo.min),
                                                     delta.pi.hig.vlo.5thperc.no.windows.mean=mean(delta.pi.hig.vlo.5thperc.no.windows),
                                                     delta.pi.hig.vlo.5thperc.mean.delta.pi.mean=mean(delta.pi.hig.vlo.5thperc.mean.delta.pi),
                                                     delta.pi.hig.vlo.5thperc.perc.windows.mean=mean(delta.pi.hig.vlo.5thperc.perc.windows )); pi.era.all.alt.wide.para.west.summ.summ

# which of these para SHDR also allo SHDR?
pi.era.all.alt.wide.para.west.summ.summ$is.allo.shdr  <- if_else(pi.era.all.alt.wide.para.west.summ.summ$shdr.para.west.id %in% era.shdr.allo.summ$shdr.para.west, "yes", "no")

###### sims ######
# # extract the same number of shdr para west and of the same width, random locations get delta.pi min, delta.pi.5thperc.no.windows, delta.pi.5thperc.mean.taj
# # first create intervals with real shdr
# era.shdr.para.west.int <- reduce(Intervals(as.matrix(era.shdr.para.west[,c("start","end")]), closed=c(TRUE,TRUE), type="R"))
# 
# # random intervals
# sim1k.era.west <- lapply(1:1000, function(x){rand_non_overlapping_intervals(era.genome_size,size(era.shdr.para.west.int))})
# 
# sim.pi.era.all.alt.wide.para.west.summ.summ <-list()
# for (i in 1:1000) {
#   # prep random intervals for action
#   sim1k.era.df <- as.data.frame(sim1k.era.west[i])
#   sim1k.era.df$ran.hdr.name <- paste("ran", i, ".hdr.", seq(1:nrow(sim1k.era.df )), sep = "")
#   
#   # all pi tmp , use only west pops here
#   delta.pi.era.tmp <- pi.era.all.alt.wide[,c("BP.wg", "delta.pi.hig.vlo", "cline")]
#   delta.pi.era.tmp$side <- substr(delta.pi.era.tmp$cline, 8,9)
#   delta.pi.era.tmp <- subset(delta.pi.era.tmp, side=="w")
#   
#   # add ran hdr names to delta.pitmp
#   delta.pi.era.tmp$ran.hdr.name <- sim1k.era.df$ran.hdr.name[mapply(match.by.range.1k, delta.pi.era.tmp$BP.wg)]
#   
#   # use observed genome-wide quantile thresholds
#   # get min td in col and ecu
#   # subset by hdr being present to avoid getting an NA category
#   sim.pi.era.all.alt.wide.para.west.col.summ <- summarise(group_by(subset(delta.pi.era.tmp, cline=="era.co.w" &ran.hdr.name!=""),ran.hdr.name, cline ),
#                                                           no.windows=n(), 
#                                                           delta.pi.hig.vlo.mean=mean(delta.pi.hig.vlo), delta.pi.hig.vlo.min=min(delta.pi.hig.vlo),
#                                                           delta.pi.hig.vlo.5thperc.no.windows= length(delta.pi.hig.vlo[delta.pi.hig.vlo  < delta.pi.hig.vlo.5th.quant.west.col ]),
#                                                           delta.pi.hig.vlo.5thperc.mean.delta.pi= mean(delta.pi.hig.vlo[delta.pi.hig.vlo < delta.pi.hig.vlo.5th.quant.west.col]),
#                                                           delta.pi.hig.vlo.5thperc.perc.windows=(delta.pi.hig.vlo.5thperc.no.windows/ no.windows)*100 )
#   sim.pi.era.all.alt.wide.para.west.ecu.summ <- summarise(group_by(subset(delta.pi.era.tmp, cline=="era.ec.w" &ran.hdr.name!=""),ran.hdr.name, cline ),
#                                                           no.windows=n(), 
#                                                           delta.pi.hig.vlo.mean=mean(delta.pi.hig.vlo), delta.pi.hig.vlo.min=min(delta.pi.hig.vlo),
#                                                           delta.pi.hig.vlo.5thperc.no.windows= length(delta.pi.hig.vlo[delta.pi.hig.vlo  < delta.pi.hig.vlo.5th.quant.west.ecu ]),
#                                                           delta.pi.hig.vlo.5thperc.mean.delta.pi= mean(delta.pi.hig.vlo[delta.pi.hig.vlo < delta.pi.hig.vlo.5th.quant.west.ecu]),
#                                                           delta.pi.hig.vlo.5thperc.perc.windows=(delta.pi.hig.vlo.5thperc.no.windows/ no.windows)*100 )  
#   
#   sim.pi.era.all.alt.wide.para.west.summ <- rbind( sim.pi.era.all.alt.wide.para.west.col.summ,  sim.pi.era.all.alt.wide.para.west.ecu.summ)
#   
#   # get summ summ
#   sim.pi.era.all.alt.wide.para.west.summ.summ[[i]] <- summarise(group_by(sim.pi.era.all.alt.wide.para.west.summ  , ran.hdr.name),
#                                                                 delta.pi.hig.vlo.mean.mean=mean(delta.pi.hig.vlo.mean),
#                                                                 delta.pi.hig.vlo.min.mean=mean(delta.pi.hig.vlo.min),
#                                                                 delta.pi.hig.vlo.5thperc.no.windows.mean=mean(delta.pi.hig.vlo.5thperc.no.windows),
#                                                                 delta.pi.hig.vlo.5thperc.mean.delta.pi.mean=mean(delta.pi.hig.vlo.5thperc.mean.delta.pi),
#                                                                 delta.pi.hig.vlo.5thperc.perc.windows.mean= mean(delta.pi.hig.vlo.5thperc.perc.windows) )
# }
# 
# sim.delta.pi.hig.vlo.era.para.west.summ.summ.df <- bind_rows(sim.pi.era.all.alt.wide.para.west.summ.summ, .id = "column_label")
# write.csv(sim.delta.pi.hig.vlo.era.para.west.summ.summ.df, "local/data/thetas/sim.delta.pi.hig.vlo.era.para.west.summ.summ.df.csv")

#### read prepped data ####
sim.delta.pi.hig.vlo.era.para.west.summ.summ.df <- read.csv( "local/data/thetas/sim.delta.pi.hig.vlo.era.para.west.summ.summ.df.csv")



###### plot #####
library(ggbeeswarm)
ggplot(data=pi.era.all.alt.wide.para.west.summ.summ, aes(y=delta.pi.hig.vlo.min.mean,x= is.allo.shdr))+
  geom_boxplot()+ theme_classic() + geom_beeswarm(color="grey")+
  geom_boxplot(inherit.aes = F, data=sim.delta.pi.hig.vlo.era.para.west.summ.summ.df, aes(y=delta.pi.hig.vlo.min.mean, x=3))+
  geom_hline(yintercept = mean(subset(pi.era.all.alt.wide, side=="e")$delta.pi.hig.vlo), lty="dashed", color="darkgrey")

ggplot(data=pi.era.all.alt.wide.para.west.summ.summ, aes(y=delta.pi.hig.vlo.5thperc.perc.windows.mean ,x= is.allo.shdr, fill=is.allo.shdr))+
  geom_boxplot()+ theme_classic() + geom_beeswarm(color="grey")+ ylim(0,60)+
  scale_fill_manual(values=c("#049E73", "#D65D00"))+
  geom_violin(inherit.aes = F, data=sim.delta.pi.hig.vlo.era.para.west.summ.summ.df, aes(y=delta.pi.hig.vlo.5thperc.perc.windows.mean, x=3))+
  geom_boxplot(inherit.aes = F, data=sim.delta.pi.hig.vlo.era.para.west.summ.summ.df, aes(y=delta.pi.hig.vlo.5thperc.perc.windows.mean, x=3), fill="grey", alpha=.5)+
  geom_hline(yintercept = quantile(sim.delta.pi.hig.vlo.era.para.west.summ.summ.df$delta.pi.hig.vlo.5thperc.perc.windows.mean, 0.95))
#geom_hline(yintercept = mean(subset(pi.era.all.alt.wide, side=="e")$Tajima), lty="dashed", color="darkgrey")



# plot min delta.pi's d per HDR vs simlated min per HDR
delta.pi.hig.vlo.dens.era.hig.w <-  ggplot(subset(pi.era.all.alt.wide, side=="e"), aes(x=delta.pi.hig.vlo, shape=cline)) + 
  scale_color_manual(values=c("#0000ff", "#5ec55b","#23ff06"), labels=c("High", "Low", "Low distant"))+
  labs(colour="Altitude")+
  geom_vline(xintercept = 0, colour="lightgrey", lty="dashed")+
  geom_vline(xintercept = subset(pi.era.all.alt.wide.para.west.summ.summ, is.allo.shdr=="no")$delta.pi.hig.vlo.min.mean, colour="#0372B2", lty="solid", alpha=.9, size=0.2)+
  geom_vline(xintercept = subset(pi.era.all.alt.wide.para.west.summ.summ, is.allo.shdr=="yes")$delta.pi.hig.vlo.min.mean, colour="#D65D00", lty="solid", alpha=.9, size=0.2)+
  geom_vline(xintercept = mean(subset(sim.delta.pi.hig.vlo.era.para.west.summ.summ.df)$delta.pi.hig.vlo.min.mean,na.rm=T), colour="black", lty="solid", alpha=.9, size=1.5)+
  geom_vline(xintercept = quantile(subset(sim.delta.pi.hig.vlo.era.para.west.summ.summ.df)$delta.pi.hig.vlo.min.mean, 0.9,na.rm=T), colour="black", lty="dashed", alpha=.9, size=1.5)+
  geom_vline(xintercept = quantile(subset(sim.delta.pi.hig.vlo.era.para.west.summ.summ.df)$delta.pi.hig.vlo.min.mean, 0.1,na.rm=T), colour="black", lty="dashed", alpha=.9, size=1.5)+
  geom_density(fill="transparent", alpha=0, size=1, colour="blue")+
  xlab(expression(paste("Delta pi high - low distant,",  italic(" H. erato"))))+
  theme_classic()+
  theme(axis.line.x = element_line(color="black", size = 0.5),axis.line.y = element_line(color="black", size = 0.5),axis.title.y = element_text(size=14),
        plot.margin = unit(c(0.5,0.2,0.2,0.2), "lines"),axis.text = element_text(size=10), #axis.text.x = element_blank(),
        axis.text.y = element_text(size=10),#axis.title.x = element_blank(), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_rect(fill = "transparent"), plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
        legend.background = element_rect(fill = "transparent"), legend.box.background = element_rect(fill = "transparent"), strip.background.y = element_blank(),
        strip.text = element_text(size=14))+
  scale_y_continuous(expand = c(0,0))+
  scale_x_continuous(expand = c(0,0)); delta.pi.hig.vlo.dens.era.hig.w

ggplot(data=pi.era.all.alt.wide.para.west.summ.summ, aes(y=delta.pi.hig.vlo.5thperc.perc.windows.mean ,x= is.allo.shdr,  colour=is.allo.shdr))+
  geom_beeswarm(cex = 2, size=2)+ ylim(0,60)+
  scale_colour_manual(values=c("#049E73", "#D65D00", labels=c("Within west\n(parapatric)", "Across sides\n(allopatric)")))+
  theme_classic() + 
  geom_rect(fill="grey80",color=NA, alpha=0.02, aes(xmin=.7,xmax=1.3, ymin=quantile(sim.delta.pi.hig.vlo.era.para.west.summ.summ.df$delta.pi.hig.vlo.5thperc.perc.windows.mean,0.05),
                                                    ymax=quantile(sim.delta.pi.hig.vlo.era.para.west.summ.summ.df$delta.pi.hig.vlo.5thperc.perc.windows.mean, 0.9)))+
  geom_rect(fill="grey80",color=NA, alpha=0.02, aes(xmin=1.7,xmax=2.3, ymin=quantile(sim.delta.pi.hig.vlo.era.para.west.summ.summ.df$delta.pi.hig.vlo.5thperc.perc.windows.mean,0.05),
                                                    ymax=quantile(sim.delta.pi.hig.vlo.era.para.west.summ.summ.df$delta.pi.hig.vlo.5thperc.perc.windows.mean, 0.9)))+
  xlab("SHDR type") + scale_x_discrete(labels=c("no" = "Within west\n(parapatric)", "yes"="Across sides\n(allopatric)"))+
  ylab("% Tajima outlier windows within SHDR")+
  theme(legend.position = "none", axis.text = element_text(size=14), axis.title = element_text(size=14))




################# 3.4 which SHDR are delta pi outliers? ######
library(patchwork)
delta.pi.hig.vlo.dens.era.hig.e / delta.pi.hig.vlo.dens.era.hig.w 
ggsave("local/plots/thetas/delta.pi.era.min.per.HDR.90-10thpercentile.hig.vlo.dens.east.west.png", width = 5, height = 5)

# check thresholds for min delta.pi.hig.vlo per HDR and % windows with delta.pi.hig.vlo <5th perc
nrow(subset(pi.era.all.alt.wide.para.east.summ.summ, delta.pi.hig.vlo.min.mean < quantile(subset(sim.delta.pi.hig.vlo.era.para.east.summ.summ.df)$delta.pi.hig.vlo.min.mean, 0.1, na.rm = T)))
nrow(subset(pi.era.all.alt.wide.para.east.summ.summ, delta.pi.hig.vlo.5thperc.perc.windows.mean > quantile(subset(sim.delta.pi.hig.vlo.era.para.east.summ.summ.df)$delta.pi.hig.vlo.5thperc.perc.windows.mean, 0.9, na.rm = T)))

nrow(subset(pi.era.all.alt.wide.para.east.summ.summ, delta.pi.hig.vlo.min.mean < quantile(subset(sim.delta.pi.hig.vlo.era.para.east.summ.summ.df)$delta.pi.hig.vlo.min.mean, 0.05, na.rm = T)))
nrow(subset(pi.era.all.alt.wide.para.east.summ.summ, delta.pi.hig.vlo.5thperc.perc.windows.mean > quantile(subset(sim.delta.pi.hig.vlo.era.para.east.summ.summ.df)$delta.pi.hig.vlo.5thperc.perc.windows.mean, 0.95, na.rm = T)))

nrow(subset(pi.era.all.alt.wide.para.west.summ.summ, delta.pi.hig.vlo.min.mean < quantile(subset(sim.delta.pi.hig.vlo.era.para.west.summ.summ.df)$delta.pi.hig.vlo.min.mean, 0.1, na.rm = T)))
nrow(subset(pi.era.all.alt.wide.para.west.summ.summ, delta.pi.hig.vlo.5thperc.perc.windows.mean > quantile(subset(sim.delta.pi.hig.vlo.era.para.west.summ.summ.df)$delta.pi.hig.vlo.5thperc.perc.windows.mean, 0.9, na.rm = T)))

nrow(subset(pi.era.all.alt.wide.para.west.summ.summ, delta.pi.hig.vlo.min.mean < quantile(subset(sim.delta.pi.hig.vlo.era.para.west.summ.summ.df)$delta.pi.hig.vlo.min.mean, 0.05, na.rm = T)))
nrow(subset(pi.era.all.alt.wide.para.west.summ.summ, delta.pi.hig.vlo.5thperc.perc.windows.mean > quantile(subset(sim.delta.pi.hig.vlo.era.para.west.summ.summ.df)$delta.pi.hig.vlo.5thperc.perc.windows.mean, 0.95, na.rm = T)))

## add whether SHDRs are outliers or not
pi.era.all.alt.wide.para.east.summ.summ$delta.pi.hig.vlo.min.less.5th.perc.sims <- if_else(pi.era.all.alt.wide.para.east.summ.summ$delta.pi.hig.vlo.min.mean < quantile(sim.delta.pi.hig.vlo.era.para.east.summ.summ.df$delta.pi.hig.vlo.min.mean, 0.05,na.rm=T), "yes", "no")
pi.era.all.alt.wide.para.east.summ.summ$delta.pi.hig.vlo.min.less.10th.perc.sims <- if_else(pi.era.all.alt.wide.para.east.summ.summ$delta.pi.hig.vlo.min.mean < quantile(sim.delta.pi.hig.vlo.era.para.east.summ.summ.df$delta.pi.hig.vlo.min.mean, 0.1,na.rm=T), "yes", "no")
pi.era.all.alt.wide.para.east.summ.summ$delta.pi.hig.vlo.5thperc.perc.windows.mean.less.95th.perc.sims <- if_else(pi.era.all.alt.wide.para.east.summ.summ$delta.pi.hig.vlo.5thperc.perc.windows.mean > quantile(sim.delta.pi.hig.vlo.era.para.east.summ.summ.df$delta.pi.hig.vlo.5thperc.perc.windows.mean, 0.95,na.rm=T), "yes", "no")
pi.era.all.alt.wide.para.east.summ.summ$delta.pi.hig.vlo.5thperc.perc.windows.mean.less.90th.perc.sims <- if_else(pi.era.all.alt.wide.para.east.summ.summ$delta.pi.hig.vlo.5thperc.perc.windows.mean > quantile(sim.delta.pi.hig.vlo.era.para.east.summ.summ.df$delta.pi.hig.vlo.5thperc.perc.windows.mean, 0.9,na.rm=T), "yes", "no")

pi.era.all.alt.wide.para.west.summ.summ$delta.pi.hig.vlo.min.less.5th.perc.sims <- if_else(pi.era.all.alt.wide.para.west.summ.summ$delta.pi.hig.vlo.min.mean < quantile(sim.delta.pi.hig.vlo.era.para.west.summ.summ.df$delta.pi.hig.vlo.min.mean, 0.05,na.rm=T), "yes", "no")
pi.era.all.alt.wide.para.west.summ.summ$delta.pi.hig.vlo.min.less.10th.perc.sims <- if_else(pi.era.all.alt.wide.para.west.summ.summ$delta.pi.hig.vlo.min.mean < quantile(sim.delta.pi.hig.vlo.era.para.west.summ.summ.df$delta.pi.hig.vlo.min.mean, 0.1,na.rm=T), "yes", "no")
pi.era.all.alt.wide.para.west.summ.summ$delta.pi.hig.vlo.5thperc.perc.windows.mean.less.95th.perc.sims <- if_else(pi.era.all.alt.wide.para.west.summ.summ$delta.pi.hig.vlo.5thperc.perc.windows.mean > quantile(sim.delta.pi.hig.vlo.era.para.west.summ.summ.df$delta.pi.hig.vlo.5thperc.perc.windows.mean, 0.95,na.rm=T), "yes", "no")
pi.era.all.alt.wide.para.west.summ.summ$delta.pi.hig.vlo.5thperc.perc.windows.mean.less.90th.perc.sims <- if_else(pi.era.all.alt.wide.para.west.summ.summ$delta.pi.hig.vlo.5thperc.perc.windows.mean > quantile(sim.delta.pi.hig.vlo.era.para.west.summ.summ.df$delta.pi.hig.vlo.5thperc.perc.windows.mean, 0.9,na.rm=T), "yes", "no"); pi.era.all.alt.wide.para.west.summ.summ

# save
pi.era.all.alt.wide.para.east.summ.summ; write.csv(pi.era.all.alt.wide.para.east.summ.summ, "local/data/thetas/shdr.pi.era.all.alt.wide.para.east.summ.summ.csv", row.names = F)
pi.era.all.alt.wide.para.west.summ.summ; write.csv(pi.era.all.alt.wide.para.west.summ.summ, "local/data/thetas/shdr.pi.era.all.alt.wide.para.west.summ.summ.csv", row.names = F)







######################################## 4. DXY PER SHDR ##########
################# 4.1 prep data ##########
## keep only hig_vlow, except in mel.ec.w where vlo not avail, use hig_low
dxy.era.sub <- subset(dxy.era.all, pop1.pop2.alt.type=="hig_vlo"); head(dxy.era.sub)
dxy.mel.sub <- subset(dxy.mel.all, (pop1.pop2.alt.type=="hig_vlo")| (cline=="mel.ec.w"&pop1.pop2.alt.type=="hig_low")); head(dxy.mel.sub)

# load, add within pop pi 
thetas.era.hig.vlo.pop <- read.csv("local/data/thetas/thetas.era.hig.vlo.pop.csv")
thetas.era.hig.vlo.pop$pi <- (thetas.era.hig.vlo.pop$tP/thetas.era.hig.vlo.pop$nSites)


# best use scaled dxy by pi
# add pi for each pop, then mean
# for melpomene will have to find hig low pi (maybe already there)
dxy.era.sub$pi.pop1 <- thetas.era.hig.vlo.pop$pi[match(paste(dxy.era.sub$pop.full1,dxy.era.sub$BP.wg), paste(thetas.era.hig.vlo.pop$pop.full, thetas.era.hig.vlo.pop$BP.wg))]
dxy.era.sub$pi.pop2 <- thetas.era.hig.vlo.pop$pi[match(paste(dxy.era.sub$pop.full2,dxy.era.sub$BP.wg), paste(thetas.era.hig.vlo.pop$pop.full, thetas.era.hig.vlo.pop$BP.wg))]
dxy.era.sub$pi.mean.pop1.pop2 <- rowMeans(dxy.era.sub[,c("pi.pop1", "pi.pop2")], na.rm = T)

thetas.era.hig.vlo.pop <- NA

#
###### plots test ####
ggplot(subset(dxy.era.sub, chr==1 ), aes(y=dxy_mean, x=BP.wg))+
  geom_rect(inherit.aes = F, data=subset(era.shdr.para.west, !(shdr.para.west.id %in% era.shdr.allo.summ$shdr.para.west)) , aes(xmin=start, xmax=end, ymin=0,ymax=.5), colour="transparent", fill="#0372B2", alpha=0.4 )+
  geom_rect(inherit.aes = F, data=subset(era.shdr.para.east, !(shdr.para.east.id %in% era.shdr.allo.summ$shdr.para.east)), aes(xmin=start, xmax=end,  ymin=0,ymax=.5), colour="transparent", fill="#049E73", alpha=0.4 )+
  geom_rect(inherit.aes = F, data=era.shdr.allo.all, aes(xmin=start, xmax=end, ymin=0,ymax=.5), colour="transparent", fill="#D65D00", alpha=0.6 )+
  #geom_point(inherit.aes = F, data=subset(dxy.era.sub, chr==1 & cline=="era.ec.e" & dxy_mean < quantile(dxy_mean , 0.05, na.rm = T)), colour="orange",aes(y=dxy_mean, x=BP.wg))+
  geom_point(size=.1)+
  scale_y_continuous( expand = c(0,0))+
  #geom_rect(inherit.aes = F, data=subset(era.shared.regions, sharing.type=="within.east"), aes(xmin=start.BP.wg-100000, xmax=end.BP.wg+100000, ymin=0,ymax=7), colour="transparent", fill="#049E73", alpha=0.7 )+
  #geom_rect(inherit.aes = F, data=subset(era.shared.regions, sharing.type=="across"), aes(xmin=start.BP.wg-100000, xmax=end.BP.wg+100000, ymin=0,ymax=7), colour="transparent", fill="#049E73", alpha=1 )+
  facet_wrap(~cline, ncol = 1)+
  xlim(10000000,20000000)+theme_bw()

ggplot(subset(dxy.era.sub), aes(y=dxy_mean, x=BP.wg))+
  geom_rect(inherit.aes = F, data=subset(era.shdr.para.west, !(shdr.para.west.id %in% era.shdr.allo.summ$shdr.para.west)) , aes(xmin=start, xmax=end, ymin=0,ymax=.5), colour="transparent", fill="#0372B2", alpha=0.4 )+
  geom_rect(inherit.aes = F, data=subset(era.shdr.para.east, !(shdr.para.east.id %in% era.shdr.allo.summ$shdr.para.east)), aes(xmin=start, xmax=end, ymin=0,ymax=.5), colour="transparent", fill="#049E73", alpha=0.4 )+
  geom_rect(inherit.aes = F, data=era.shdr.allo.all, aes(xmin=start, xmax=end, ymin=0,ymax=.5), colour="transparent", fill="#D65D00", alpha=0.6 )+
  geom_point(size=.1)+
  scale_y_continuous(  expand = c(0,0))+
  #geom_rect(inherit.aes = F, data=subset(era.shared.regions, sharing.type=="within.east"), aes(xmin=start.BP.wg-100000, xmax=end.BP.wg+100000, ymin=0,ymax=7), colour="transparent", fill="#049E73", alpha=0.7 )+
  #geom_rect(inherit.aes = F, data=subset(era.shared.regions, sharing.type=="across"), aes(xmin=start.BP.wg-100000, xmax=end.BP.wg+100000, ymin=0,ymax=7), colour="transparent", fill="#049E73", alpha=1 )+
  facet_wrap(~cline, ncol = 1)+
  xlim(25000000,27000000)+theme_bw()



ggplot(subset(dxy.era.sub, chr==2), aes(y=dxy_mean, x=BP.wg))+
  geom_rect(inherit.aes = F, data=subset(era.shdr.para.west, !(shdr.para.west.id %in% era.shdr.allo.summ$shdr.para.west) &chr==2) , aes(xmin=start, xmax=end, ymin=0,ymax=.5), colour="transparent", fill="#0372B2", alpha=0.4 )+
  geom_rect(inherit.aes = F, data=subset(era.shdr.para.east, !(shdr.para.east.id %in% era.shdr.allo.summ$shdr.para.east) & chr==2), aes(xmin=start, xmax=end, ymin=0,ymax=.5), colour="transparent", fill="#049E73", alpha=0.4 )+
  geom_rect(inherit.aes = F, data=subset(era.shdr.allo.all,chr==2), aes(xmin=start, xmax=end, ymin=0,ymax=.5), colour="transparent", fill="#D65D00", alpha=0.6 )+
  geom_point(inherit.aes = F, data=subset(dxy.era.sub, chr==2 &cline=="era.ec.e" & dxy_mean < quantile(dxy_mean , 0.05, na.rm = T)), colour="orange",aes(y=dxy_mean, x=BP.wg))+
  geom_point(size=.1)+
  scale_y_continuous(  expand = c(0,0))+
  facet_wrap(~cline, ncol = 1)+theme_bw()

ggplot(subset(dxy.era.all, chr==1  & cline=="era.co.e"), aes(y=dxy_mean, x=BP.wg))+
  geom_rect(inherit.aes = F, data=subset(era.shdr.para.west, !(shdr.para.west.id %in% era.shdr.allo.summ$shdr.para.west) &chr==1) , aes(xmin=start, xmax=end, ymin=0,ymax=.5), colour="transparent", fill="#0371B1", alpha=0.4 )+
  geom_rect(inherit.aes = F, data=subset(era.shdr.para.east, !(shdr.para.east.id %in% era.shdr.allo.summ$shdr.para.east) & chr==1), aes(xmin=start, xmax=end, ymin=0,ymax=.5), colour="transparent", fill="#049E73", alpha=0.4 )+
  geom_rect(inherit.aes = F, data=subset(era.shdr.allo.all,chr==1), aes(xmin=start, xmax=end, ymin=0,ymax=.5), colour="transparent", fill="#D65D00", alpha=0.6 )+
  geom_point(inherit.aes = F, data=subset(dxy.era.sub, chr==1 &cline=="era.ec.e" & dxy_mean < quantile(dxy_mean , 0.05, na.rm = T)), colour="orange",aes(y=dxy_mean, x=BP.wg))+
  geom_hline(yintercept = dxy_mean.95th.quant.east.col)+
  geom_point(size=.1)+
  scale_y_continuous(  expand = c(0,0))+
  facet_wrap(~cline+pop1.pop2.alt.type, ncol = 1)+theme_bw()

# dxy pi
ggplot(subset(dxy.era.sub), aes(y=dxy_mean, x=pi.mean.pop1.pop2))+
  geom_abline(intercept = 0, slope = 1, colour="red")+
  geom_point(size=.1)+scale_y_continuous( expand = c(0.05,0.05))+
  geom_point(inherit.aes = F, data=subset(dxy.era.sub, cline=="era.ec.e" & dxy_mean < quantile(dxy_mean , 0.05, na.rm = T)), colour="orange",alpha=.1, aes(y=dxy_mean, x=pi.mean.pop1.pop2))+
  geom_point(inherit.aes = F, data=subset(dxy.era.sub, cline=="era.ec.e" & dxy_mean > quantile(dxy_mean , 0.95, na.rm = T)), colour="orange",alpha=.1, aes(y=dxy_mean, x=pi.mean.pop1.pop2))+
  geom_smooth(method="lm")+
  stat_cor(label.y.npc = 0.6)+ stat_regline_equation()+
  facet_wrap(~cline, ncol = 2)+theme_bw()

ggplot(subset(dxy.era.sub, chr==2), aes(y=dxy_mean/pi.mean.pop1.pop2, x=cline))+
  geom_boxplot()+
  geom_point(size=.1)+scale_y_continuous( expand = c(0.01,0))+
  geom_point(inherit.aes = F, data=subset(dxy.era.sub, cline=="era.ec.e" & dxy_mean < quantile(dxy_mean , 0.05, na.rm = T)), colour="orange",alpha=.1,aes(y=dxy_mean/pi.mean.pop1.pop2, x=cline))+
  geom_point(inherit.aes = F, data=subset(dxy.era.sub,cline=="era.ec.e" & dxy_mean > quantile(dxy_mean , 0.95, na.rm = T)), colour="orange",alpha=.1,aes(y=dxy_mean/pi.mean.pop1.pop2, x=cline))+
  #geom_smooth(method="lm")+stat_cor(label.y.npc = 0.6)+ stat_regline_equation()+
  theme_bw()


################# 4.2 SHDR east ######
head(dxy.era.sub)

# calculate genomewide quantile threshold for each pop (one high alt per country)
dxy_mean.5th.quant.east.col <- quantile(subset(dxy.era.sub, cline=="era.co.e")$dxy_mean, 0.05, na.rm = T ); dxy_mean.5th.quant.east.col
dxy_mean.5th.quant.east.ecu <- quantile(subset(dxy.era.sub, cline=="era.ec.e")$dxy_mean, 0.05, na.rm = T); dxy_mean.5th.quant.east.ecu

dxy_mean.95th.quant.east.col <- quantile(subset(dxy.era.sub, cline=="era.co.e")$dxy_mean, 0.95, na.rm = T ); dxy_mean.95th.quant.east.col
dxy_mean.95th.quant.east.ecu <- quantile(subset(dxy.era.sub, cline=="era.ec.e")$dxy_mean, 0.95, na.rm = T); dxy_mean.95th.quant.east.ecu

# for now only obtain dxy max min outliers
head(dxy.era.sub)
dxy.era.sub.para.east.col.summ <- summarise(group_by(subset(dxy.era.sub, cline=="era.co.e" &shdr.para.east.id!="" ), shdr.para.east.id, cline ),
                                                    start=min(BP.wg),
                                                    end=max(BP.wg),
                                                    size.hdr.bp=max(BP.wg)- min(BP.wg),
                                                    no.windows=n(), chr=unique(chr),
                                                    dxy.mean=mean(dxy_mean, na.rm = T), dxy.min=min(dxy_mean, na.rm = T), dxy.max=max(dxy_mean, na.rm = T),
                                                    dxy.5thperc.perc.windows= (length(dxy_mean[dxy_mean  < dxy_mean.5th.quant.east.col ])/no.windows)*100,
                                            dxy.95thperc.perc.windows= (length(dxy_mean[dxy_mean  > dxy_mean.95th.quant.east.col ])/no.windows)*100,
                                                    dxy.5thperc.mean.dxy= mean(dxy_mean[dxy_mean < dxy_mean.5th.quant.east.col], na.rm = T),
                                            dxy.95thperc.mean.dxy= mean(dxy_mean[dxy_mean > dxy_mean.95th.quant.east.col], na.rm = T)); dxy.era.sub.para.east.col.summ 

dxy.era.sub.para.east.ecu.summ <- summarise(group_by(subset(dxy.era.sub, cline=="era.ec.e" &shdr.para.east.id!="" ), shdr.para.east.id, cline ),
                                            start=min(BP.wg),
                                            end=max(BP.wg),
                                            size.hdr.bp=max(BP.wg)- min(BP.wg),
                                            no.windows=n(), chr=unique(chr),
                                            dxy.mean=mean(dxy_mean, na.rm = T), dxy.min=min(dxy_mean, na.rm = T), dxy.max=max(dxy_mean, na.rm = T),
                                            dxy.5thperc.perc.windows= (length(dxy_mean[dxy_mean  < dxy_mean.5th.quant.east.ecu ])/no.windows)*100,
                                            dxy.95thperc.perc.windows= (length(dxy_mean[dxy_mean  > dxy_mean.95th.quant.east.ecu ])/no.windows)*100,
                                            dxy.5thperc.mean.dxy= mean(dxy_mean[dxy_mean < dxy_mean.5th.quant.east.ecu], na.rm = T),
                                            dxy.95thperc.mean.dxy= mean(dxy_mean[dxy_mean > dxy_mean.95th.quant.east.ecu], na.rm = T)); dxy.era.sub.para.east.ecu.summ 

dxy.era.sub.para.east.summ <- rbind(dxy.era.sub.para.east.col.summ, dxy.era.sub.para.east.ecu.summ)

# summarise to get mean lowest
dxy.era.sub.para.east.summ.summ <- summarise(group_by(dxy.era.sub.para.east.summ ,shdr.para.east.id,start, end, size.hdr.bp, no.windows ),
                                                     chr=unique(chr),
                                                     dxy.mean.mean=mean(dxy.mean),
                                                     dxy.min.mean=mean(dxy.min),
                                             dxy.max.mean=mean(dxy.max),
                                                     dxy.5thperc.mean.dxy.mean=mean(dxy.5thperc.mean.dxy),
                                                     dxy.5thperc.perc.windows.mean=mean(dxy.5thperc.perc.windows ),
                                             dxy.95thperc.mean.dxy.mean=mean(dxy.95thperc.mean.dxy),
                                             dxy.95thperc.perc.windows.mean=mean(dxy.95thperc.perc.windows )); dxy.era.sub.para.east.summ.summ

# which of these para SHDR also allo SHDR?
dxy.era.sub.para.east.summ.summ$is.allo.shdr  <- if_else(dxy.era.sub.para.east.summ.summ$shdr.para.east.id %in% era.shdr.allo.summ$shdr.para.east, "yes", "no")

###### sims ######
# # extract the same number of shdr para east and of the same width, random locations get dxy min, dxy.5thperc.no.windows, dxy.5thperc.mean.taj
# # first create intervals with real shdr
# era.shdr.para.east.int <- reduce(Intervals(as.matrix(era.shdr.para.east[,c("start","end")]), closed=c(TRUE,TRUE), type="R"))
# 
# # random intervals
# sim1k.era.east <- lapply(1:1000, function(x){rand_non_overlapping_intervals(era.genome_size,size(era.shdr.para.east.int))})
# 
# sim.dxy.era.sub.para.east.summ.summ <-list()
# for (i in 1:1000) {
#   # prep random intervals for action
#   sim1k.era.df <- as.data.frame(sim1k.era.east[i])
#   sim1k.era.df$ran.hdr.name <- paste("ran", i, ".hdr.", seq(1:nrow(sim1k.era.df )), sep = "")
#   
#   # all dxy tmp , use only east pops here
#   dxy.era.tmp <- dxy.era.sub[,c("BP.wg", "dxy_mean", "cline")]
#   dxy.era.tmp$side <- substr(dxy.era.tmp$cline, 8,9)
#   dxy.era.tmp <- subset(dxy.era.tmp, side=="e")
#   
#   # add ran hdr names to dxytmp
#   dxy.era.tmp$ran.hdr.name <- sim1k.era.df$ran.hdr.name[mapply(match.by.range.1k, dxy.era.tmp$BP.wg)]
#   
#   # use observed genome-wide quantile thresholds
#   # get min td in col and ecu
#   # subset by hdr being present to avoid getting an NA category
#   sim.dxy.era.sub.para.east.col.summ <- summarise(group_by(subset(dxy.era.tmp, cline=="era.co.e" & ran.hdr.name!=""),ran.hdr.name, cline ),
#                                                   no.windows=n(),
#                                                   dxy.mean=mean(dxy_mean, na.rm = T), dxy.min=min(dxy_mean, na.rm = T), dxy.max=max(dxy_mean, na.rm = T),
#                                                   dxy.5thperc.perc.windows= (length(dxy_mean[dxy_mean  < dxy_mean.5th.quant.east.col ])/no.windows)*100,
#                                                   dxy.95thperc.perc.windows= (length(dxy_mean[dxy_mean  > dxy_mean.95th.quant.east.col ])/no.windows)*100,
#                                                   dxy.5thperc.mean.dxy= mean(dxy_mean[dxy_mean < dxy_mean.5th.quant.east.col], na.rm = T),
#                                                   dxy.95thperc.mean.dxy= mean(dxy_mean[dxy_mean > dxy_mean.95th.quant.east.col], na.rm = T))
#     
#   sim.dxy.era.sub.para.east.ecu.summ <- summarise(group_by(subset(dxy.era.tmp, cline=="era.ec.e" & ran.hdr.name!=""),ran.hdr.name, cline ),
#                                                   no.windows=n(),
#                                                   dxy.mean=mean(dxy_mean, na.rm = T), dxy.min=min(dxy_mean, na.rm = T), dxy.max=max(dxy_mean, na.rm = T),
#                                                   dxy.5thperc.perc.windows= (length(dxy_mean[dxy_mean  < dxy_mean.5th.quant.east.ecu ])/no.windows)*100,
#                                                   dxy.95thperc.perc.windows= (length(dxy_mean[dxy_mean  > dxy_mean.95th.quant.east.ecu ])/no.windows)*100,
#                                                   dxy.5thperc.mean.dxy= mean(dxy_mean[dxy_mean < dxy_mean.5th.quant.east.ecu], na.rm = T),
#                                                   dxy.95thperc.mean.dxy= mean(dxy_mean[dxy_mean > dxy_mean.95th.quant.east.ecu], na.rm = T))   
#   
#   sim.dxy.era.sub.para.east.summ <- rbind( sim.dxy.era.sub.para.east.col.summ,  sim.dxy.era.sub.para.east.ecu.summ)
#   
#   # get summ summ, if one of the cline SHDR summs had 0% of windows as outliers the mean for the mean dxy in top/bottom 5% will be NA (as one is NA,not zero)
#   sim.dxy.era.sub.para.east.summ.summ[[i]] <- summarise(group_by(sim.dxy.era.sub.para.east.summ  , ran.hdr.name),
#                                                         dxy.mean.mean=mean(dxy.mean),
#                                                         dxy.min.mean=mean(dxy.min),
#                                                         dxy.max.mean=mean(dxy.max),
#                                                         dxy.5thperc.mean.dxy.mean=mean(dxy.5thperc.mean.dxy),
#                                                         dxy.5thperc.perc.windows.mean=mean(dxy.5thperc.perc.windows ),
#                                                         dxy.95thperc.mean.dxy.mean=mean(dxy.95thperc.mean.dxy),
#                                                         dxy.95thperc.perc.windows.mean=mean(dxy.95thperc.perc.windows ))
# }
# 
# sim.dxy.era.para.east.summ.summ.df <- bind_rows(sim.dxy.era.sub.para.east.summ.summ, .id = "column_label")
# write.csv(sim.dxy.era.para.east.summ.summ.df, "local/data/dxy/sim.1k.dxy.era.para.east.summ.summ.df.csv")

sim.dxy.era.para.east.summ.summ.df <- read.csv( "local/data/dxy/sim.1k.dxy.era.para.east.summ.summ.df.csv")


###### plot #####
library(ggbeeswarm)
dxy.era.sub.para.east.col.summ 

ggplot(data=dxy.era.sub.para.east.summ.summ, aes(y=dxy.min.mean,x= is.allo.shdr))+
  geom_boxplot()+ theme_classic() + geom_beeswarm(color="grey")+
  geom_boxplot(inherit.aes = F, data=sim.dxy.era.para.east.summ.summ.df, aes(y=dxy.min.mean, x=3))+
  geom_hline(yintercept = mean(subset(dxy.era.sub, substr(cline, 8,8)=="e")$dxy_mean), lty="dashed", color="darkgrey")

ggplot(data=dxy.era.sub.para.east.summ.summ, aes(y=dxy.max.mean, x= is.allo.shdr))+
  geom_boxplot()+ theme_classic() + geom_beeswarm(color="grey")+
  geom_boxplot(inherit.aes = F, data=sim.dxy.era.para.east.summ.summ.df, aes(y=dxy.max.mean, x=3))+
  geom_hline(yintercept = mean(subset(dxy.era.sub, substr(cline, 8,8)=="e")$dxy_mean), lty="dashed", color="darkgrey")

ggplot(data=dxy.era.sub.para.east.summ.summ, aes(y=dxy.5thperc.perc.windows.mean ,x= is.allo.shdr, fill=is.allo.shdr))+
  geom_boxplot()+ theme_classic() + geom_beeswarm(color="grey")+ ylim(0,60)+
  scale_fill_manual(values=c("#049E73", "#D65D00"))+
  geom_violin(inherit.aes = F, data=sim.dxy.era.para.east.summ.summ.df, aes(y=dxy.5thperc.perc.windows.mean, x=3))+
  geom_boxplot(inherit.aes = F, data=sim.dxy.era.para.east.summ.summ.df, aes(y=dxy.5thperc.perc.windows.mean, x=3), fill="grey", alpha=.5)+
  geom_hline(yintercept = quantile(sim.dxy.era.para.east.summ.summ.df$dxy.5thperc.perc.windows.mean, 0.95))
#geom_hline(yintercept = mean(subset(dxy.era.sub, side=="e")$Tajima), lty="dashed", color="darkgrey")

ggplot(data=dxy.era.sub.para.east.summ.summ, aes(y=dxy.95thperc.perc.windows.mean ,x= is.allo.shdr, fill=is.allo.shdr))+
  geom_boxplot()+ theme_classic() + geom_beeswarm(color="grey")+ ylim(0,90)+
  scale_fill_manual(values=c("#049E73", "#D65D00"))+
  geom_violin(inherit.aes = F, data=sim.dxy.era.para.east.summ.summ.df, aes(y=dxy.95thperc.perc.windows.mean, x=3))+
  geom_boxplot(inherit.aes = F, data=sim.dxy.era.para.east.summ.summ.df, aes(y=dxy.95thperc.perc.windows.mean, x=3), fill="grey", alpha=.5)+
  geom_hline(yintercept = quantile(sim.dxy.era.para.east.summ.summ.df$dxy.95thperc.perc.windows.mean, 0.95))
#geom_hline(yintercept = mean(subset(dxy.era.sub, side=="e")$Tajima), lty="dashed", color="darkgrey")



# plot min dxy's d per HDR vs simlated min per HDR
dxy.min.dens.era.e <-  ggplot(subset(dxy.era.sub, substr(cline, 8,8)=="e"), aes(x=dxy_mean, shape=cline)) + 
  scale_color_manual(values=c("#0000ff", "#5ec55b","#23ff06"), labels=c("High", "Low", "Low distant"))+
  labs(colour="Altitude")+
  geom_vline(xintercept = 0, colour="lightgrey", lty="dashed")+
  geom_vline(xintercept = subset(dxy.era.sub.para.east.summ.summ, is.allo.shdr=="no")$dxy.min.mean, colour="#049E73", lty="solid", alpha=.9, size=0.2)+
  geom_vline(xintercept = subset(dxy.era.sub.para.east.summ.summ, is.allo.shdr=="yes")$dxy.min.mean, colour="#D65D00", lty="solid", alpha=.9, size=0.2)+
  geom_vline(xintercept = mean(subset(sim.dxy.era.para.east.summ.summ.df)$dxy.min.mean,na.rm=T), colour="black", lty="solid", alpha=.9, size=1.5)+
  geom_vline(xintercept = quantile(subset(sim.dxy.era.para.east.summ.summ.df)$dxy.min.mean, 0.9,na.rm=T), colour="black", lty="dashed", alpha=.9, size=1.5)+
  geom_vline(xintercept = quantile(subset(sim.dxy.era.para.east.summ.summ.df)$dxy.min.mean, 0.1,na.rm=T), colour="black", lty="dashed", alpha=.9, size=1.5)+
  geom_density(fill="transparent", alpha=0, size=1, colour="blue")+
  xlab(expression(paste("Tajima's D",  italic(" H. erato"))))+
  theme_classic()+
  theme(axis.line.x = element_line(color="black", size = 0.5),axis.line.y = element_line(color="black", size = 0.5),axis.title.y = element_text(size=14),
        plot.margin = unit(c(0.5,0.2,0.2,0.2), "lines"),axis.text = element_text(size=10), #axis.text.x = element_blank(),
        axis.text.y = element_text(size=10),axis.title.x = element_blank(), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_rect(fill = "transparent"), plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
        legend.background = element_rect(fill = "transparent"), legend.box.background = element_rect(fill = "transparent"), strip.background.y = element_blank(),
        strip.text = element_text(size=14))+
  scale_y_continuous(expand = c(0,0))+
  scale_x_continuous(expand = c(0,0)); dxy.min.dens.era.e

dxy.max.dens.era.e <-  ggplot(subset(dxy.era.sub, substr(cline, 8,8)=="e"), aes(x=dxy_mean, shape=cline)) + 
  scale_color_manual(values=c("#0000ff", "#5ec55b","#23ff06"), labels=c("High", "Low", "Low distant"))+
  labs(colour="Altitude")+
  geom_vline(xintercept = 0, colour="lightgrey", lty="dashed")+
  geom_vline(xintercept = subset(dxy.era.sub.para.east.summ.summ, is.allo.shdr=="no")$dxy.max.mean, colour="#049E73", lty="solid", alpha=.9, size=0.2)+
  geom_vline(xintercept = subset(dxy.era.sub.para.east.summ.summ, is.allo.shdr=="yes")$dxy.max.mean, colour="#D65D00", lty="solid", alpha=.9, size=0.2)+
  geom_vline(xintercept = mean(subset(sim.dxy.era.para.east.summ.summ.df)$dxy.max.mean,na.rm=T), colour="black", lty="solid", alpha=.9, size=1.5)+
  geom_vline(xintercept = quantile(subset(sim.dxy.era.para.east.summ.summ.df)$dxy.max.mean, 0.9,na.rm=T), colour="black", lty="dashed", alpha=.9, size=1.5)+
  geom_vline(xintercept = quantile(subset(sim.dxy.era.para.east.summ.summ.df)$dxy.max.mean, 0.1,na.rm=T), colour="black", lty="dashed", alpha=.9, size=1.5)+
  geom_density(fill="transparent", alpha=0, size=1, colour="blue")+
  xlab(expression(paste("Tajima's D",  italic(" H. erato"))))+
  theme_classic()+
  theme(axis.line.x = element_line(color="black", size = 0.5),axis.line.y = element_line(color="black", size = 0.5),axis.title.y = element_text(size=14),
        plot.margin = unit(c(0.5,0.2,0.2,0.2), "lines"),axis.text = element_text(size=10), #axis.text.x = element_blank(),
        axis.text.y = element_text(size=10),axis.title.x = element_blank(), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_rect(fill = "transparent"), plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
        legend.background = element_rect(fill = "transparent"), legend.box.background = element_rect(fill = "transparent"), strip.background.y = element_blank(),
        strip.text = element_text(size=14))+
  scale_y_continuous(expand = c(0,0))+
  scale_x_continuous(expand = c(0,0)); dxy.max.dens.era.e

################# 4.2 SHDR west ######
head(dxy.era.sub)

# calculate genomewide quantile threshold for each pop (one high alt per country)
dxy_mean.5th.quant.west.col <- quantile(subset(dxy.era.sub, cline=="era.co.w")$dxy_mean, 0.05, na.rm = T ); dxy_mean.5th.quant.west.col
dxy_mean.5th.quant.west.ecu <- quantile(subset(dxy.era.sub, cline=="era.ec.w")$dxy_mean, 0.05, na.rm = T); dxy_mean.5th.quant.west.ecu

dxy_mean.95th.quant.west.col <- quantile(subset(dxy.era.sub, cline=="era.co.w")$dxy_mean, 0.95, na.rm = T ); dxy_mean.95th.quant.west.col
dxy_mean.95th.quant.west.ecu <- quantile(subset(dxy.era.sub, cline=="era.ec.w")$dxy_mean, 0.95, na.rm = T); dxy_mean.95th.quant.west.ecu

# for now only obtain dxy max min outliers
head(dxy.era.sub)
dxy.era.sub.para.west.col.summ <- summarise(group_by(subset(dxy.era.sub, cline=="era.co.w" &shdr.para.west.id!="" ), shdr.para.west.id, cline ),
                                            start=min(BP.wg),
                                            end=max(BP.wg),
                                            size.hdr.bp=max(BP.wg)- min(BP.wg),
                                            no.windows=n(), chr=unique(chr),
                                            dxy.mean=mean(dxy_mean, na.rm = T), dxy.min=min(dxy_mean, na.rm = T), dxy.max=max(dxy_mean, na.rm = T),
                                            dxy.5thperc.perc.windows= (length(dxy_mean[dxy_mean  < dxy_mean.5th.quant.west.col ])/no.windows)*100,
                                            dxy.95thperc.perc.windows= (length(dxy_mean[dxy_mean  > dxy_mean.95th.quant.west.col ])/no.windows)*100,
                                            dxy.5thperc.mean.dxy= mean(dxy_mean[dxy_mean < dxy_mean.5th.quant.west.col], na.rm = T),
                                            dxy.95thperc.mean.dxy= mean(dxy_mean[dxy_mean > dxy_mean.95th.quant.west.col], na.rm = T)); dxy.era.sub.para.west.col.summ 

dxy.era.sub.para.west.ecu.summ <- summarise(group_by(subset(dxy.era.sub, cline=="era.ec.w" &shdr.para.west.id!="" ), shdr.para.west.id, cline ),
                                            start=min(BP.wg),
                                            end=max(BP.wg),
                                            size.hdr.bp=max(BP.wg)- min(BP.wg),
                                            no.windows=n(), chr=unique(chr),
                                            dxy.mean=mean(dxy_mean, na.rm = T), dxy.min=min(dxy_mean, na.rm = T), dxy.max=max(dxy_mean, na.rm = T),
                                            dxy.5thperc.perc.windows= (length(dxy_mean[dxy_mean  < dxy_mean.5th.quant.west.ecu ])/no.windows)*100,
                                            dxy.95thperc.perc.windows= (length(dxy_mean[dxy_mean  > dxy_mean.95th.quant.west.ecu ])/no.windows)*100,
                                            dxy.5thperc.mean.dxy= mean(dxy_mean[dxy_mean < dxy_mean.5th.quant.west.ecu], na.rm = T),
                                            dxy.95thperc.mean.dxy= mean(dxy_mean[dxy_mean > dxy_mean.95th.quant.west.ecu], na.rm = T)); dxy.era.sub.para.west.ecu.summ 

dxy.era.sub.para.west.summ <- rbind(dxy.era.sub.para.west.col.summ, dxy.era.sub.para.west.ecu.summ)

# summarise to get mean lowest
dxy.era.sub.para.west.summ.summ <- summarise(group_by(dxy.era.sub.para.west.summ ,shdr.para.west.id,start, end, size.hdr.bp, no.windows ),
                                             chr=unique(chr),
                                             dxy.mean.mean=mean(dxy.mean),
                                             dxy.min.mean=mean(dxy.min),
                                             dxy.max.mean=mean(dxy.max),
                                             dxy.5thperc.mean.dxy.mean=mean(dxy.5thperc.mean.dxy),
                                             dxy.5thperc.perc.windows.mean=mean(dxy.5thperc.perc.windows ),
                                             dxy.95thperc.mean.dxy.mean=mean(dxy.95thperc.mean.dxy),
                                             dxy.95thperc.perc.windows.mean=mean(dxy.95thperc.perc.windows )); dxy.era.sub.para.west.summ.summ

# which of these para SHDR also allo SHDR?
dxy.era.sub.para.west.summ.summ$is.allo.shdr  <- if_else(dxy.era.sub.para.west.summ.summ$shdr.para.west.id %in% era.shdr.allo.summ$shdr.para.west, "yes", "no")

###### sims ######
# # extract the same number of shdr para west and of the same width, random locations get dxy min, dxy.5thperc.no.windows, dxy.5thperc.mean.taj
# # first create intervals with real shdr
# era.shdr.para.west.int <- reduce(Intervals(as.matrix(era.shdr.para.west[,c("start","end")]), closed=c(TRUE,TRUE), type="R"))
# 
# # random intervals
# sim1k.era.west <- lapply(1:1000, function(x){rand_non_overlapping_intervals(era.genome_size,size(era.shdr.para.west.int))})
# 
# sim.dxy.era.sub.para.west.summ.summ <-list()
# for (i in 1:1000) {
#   # prep random intervals for action
#   sim1k.era.df <- as.data.frame(sim1k.era.west[i])
#   sim1k.era.df$ran.hdr.name <- paste("ran", i, ".hdr.", seq(1:nrow(sim1k.era.df )), sep = "")
#   
#   # all dxy tmp , use only west pops here
#   dxy.era.tmp <- dxy.era.sub[,c("BP.wg", "dxy_mean", "cline")]
#   dxy.era.tmp$side <- substr(dxy.era.tmp$cline, 8,9)
#   dxy.era.tmp <- subset(dxy.era.tmp, side=="w")
#   
#   # add ran hdr names to dxytmp
#   dxy.era.tmp$ran.hdr.name <- sim1k.era.df$ran.hdr.name[mapply(match.by.range.1k, dxy.era.tmp$BP.wg)]
#   
#   # use observed genome-wide quantile thresholds
#   # get min td in col and ecu
#   # subset by hdr being present to avoid getting an NA category
#   sim.dxy.era.sub.para.west.col.summ <- summarise(group_by(subset(dxy.era.tmp, cline=="era.co.w" & ran.hdr.name!=""),ran.hdr.name, cline ),
#                                                   no.windows=n(),
#                                                   dxy.mean=mean(dxy_mean, na.rm = T), dxy.min=min(dxy_mean, na.rm = T), dxy.max=max(dxy_mean, na.rm = T),
#                                                   dxy.5thperc.perc.windows= (length(dxy_mean[dxy_mean  < dxy_mean.5th.quant.west.col ])/no.windows)*100,
#                                                   dxy.95thperc.perc.windows= (length(dxy_mean[dxy_mean  > dxy_mean.95th.quant.west.col ])/no.windows)*100,
#                                                   dxy.5thperc.mean.dxy= mean(dxy_mean[dxy_mean < dxy_mean.5th.quant.west.col], na.rm = T),
#                                                   dxy.95thperc.mean.dxy= mean(dxy_mean[dxy_mean > dxy_mean.95th.quant.west.col], na.rm = T))
#   
#   sim.dxy.era.sub.para.west.ecu.summ <- summarise(group_by(subset(dxy.era.tmp, cline=="era.ec.w" & ran.hdr.name!=""),ran.hdr.name, cline ),
#                                                   no.windows=n(),
#                                                   dxy.mean=mean(dxy_mean, na.rm = T), dxy.min=min(dxy_mean, na.rm = T), dxy.max=max(dxy_mean, na.rm = T),
#                                                   dxy.5thperc.perc.windows= (length(dxy_mean[dxy_mean  < dxy_mean.5th.quant.west.ecu ])/no.windows)*100,
#                                                   dxy.95thperc.perc.windows= (length(dxy_mean[dxy_mean  > dxy_mean.95th.quant.west.ecu ])/no.windows)*100,
#                                                   dxy.5thperc.mean.dxy= mean(dxy_mean[dxy_mean < dxy_mean.5th.quant.west.ecu], na.rm = T),
#                                                   dxy.95thperc.mean.dxy= mean(dxy_mean[dxy_mean > dxy_mean.95th.quant.west.ecu], na.rm = T))   
#   
#   sim.dxy.era.sub.para.west.summ <- rbind( sim.dxy.era.sub.para.west.col.summ,  sim.dxy.era.sub.para.west.ecu.summ)
#   
#   # get summ summ, if one of the cline SHDR summs had 0% of windows as outliers the mean for the mean dxy in top/bottom 5% will be NA (as one is NA,not zero)
#   sim.dxy.era.sub.para.west.summ.summ[[i]] <- summarise(group_by(sim.dxy.era.sub.para.west.summ  , ran.hdr.name),
#                                                         dxy.mean.mean=mean(dxy.mean),
#                                                         dxy.min.mean=mean(dxy.min),
#                                                         dxy.max.mean=mean(dxy.max),
#                                                         dxy.5thperc.mean.dxy.mean=mean(dxy.5thperc.mean.dxy),
#                                                         dxy.5thperc.perc.windows.mean=mean(dxy.5thperc.perc.windows ),
#                                                         dxy.95thperc.mean.dxy.mean=mean(dxy.95thperc.mean.dxy),
#                                                         dxy.95thperc.perc.windows.mean=mean(dxy.95thperc.perc.windows ))
# }
# 
# sim.dxy.era.para.west.summ.summ.df <- bind_rows(sim.dxy.era.sub.para.west.summ.summ, .id = "column_label")
# write.csv(sim.dxy.era.para.west.summ.summ.df, "local/data/dxy/sim.1k.dxy.era.para.west.summ.summ.df.csv")

sim.dxy.era.para.west.summ.summ.df <- read.csv( "local/data/dxy/sim.1k.dxy.era.para.west.summ.summ.df.csv")


###### plot #####
library(ggbeeswarm)
dxy.era.sub.para.west.col.summ 

ggplot(data=dxy.era.sub.para.west.summ.summ, aes(y=dxy.min.mean,x= is.allo.shdr))+
  geom_boxplot()+ theme_classic() + geom_beeswarm(color="grey")+
  geom_boxplot(inherit.aes = F, data=sim.dxy.era.para.west.summ.summ.df, aes(y=dxy.min.mean, x=3))+
  geom_hline(yintercept = mean(subset(dxy.era.sub, substr(cline, 8,8)=="w")$dxy_mean), lty="dashed", color="darkgrey")

ggplot(data=dxy.era.sub.para.west.summ.summ, aes(y=dxy.max.mean, x= is.allo.shdr))+
  geom_boxplot()+ theme_classic() + geom_beeswarm(color="grey")+
  geom_boxplot(inherit.aes = F, data=sim.dxy.era.para.west.summ.summ.df, aes(y=dxy.max.mean, x=3))+
  geom_hline(yintercept = mean(subset(dxy.era.sub, substr(cline, 8,8)=="w")$dxy_mean), lty="dashed", color="darkgrey")

ggplot(data=dxy.era.sub.para.west.summ.summ, aes(y=dxy.5thperc.perc.windows.mean ,x= is.allo.shdr, fill=is.allo.shdr))+
  geom_boxplot()+ theme_classic() + geom_beeswarm(color="grey")+ ylim(0,60)+
  scale_fill_manual(values=c("#0372B2", "#D65D00"))+
  geom_violin(inherit.aes = F, data=sim.dxy.era.para.west.summ.summ.df, aes(y=dxy.5thperc.perc.windows.mean, x=3))+
  geom_boxplot(inherit.aes = F, data=sim.dxy.era.para.west.summ.summ.df, aes(y=dxy.5thperc.perc.windows.mean, x=3), fill="grey", alpha=.5)+
  geom_hline(yintercept = quantile(sim.dxy.era.para.west.summ.summ.df$dxy.5thperc.perc.windows.mean, 0.95))
#geom_hline(yintercept = mean(subset(dxy.era.sub, side=="w")$Tajima), lty="dashed", color="darkgrey")

ggplot(data=dxy.era.sub.para.west.summ.summ, aes(y=dxy.95thperc.perc.windows.mean ,x= is.allo.shdr, fill=is.allo.shdr))+
  geom_boxplot()+ theme_classic() + geom_beeswarm(color="grey")+ ylim(0,60)+
  scale_fill_manual(values=c("#0372B2", "#D65D00"))+
  geom_violin(inherit.aes = F, data=sim.dxy.era.para.west.summ.summ.df, aes(y=dxy.95thperc.perc.windows.mean, x=3))+
  geom_boxplot(inherit.aes = F, data=sim.dxy.era.para.west.summ.summ.df, aes(y=dxy.95thperc.perc.windows.mean, x=3), fill="grey", alpha=.5)+
  geom_hline(yintercept = quantile(sim.dxy.era.para.west.summ.summ.df$dxy.95thperc.perc.windows.mean, 0.95))
#geom_hline(yintercept = mean(subset(dxy.era.sub, side=="w")$Tajima), lty="dashed", color="darkgrey")



# plot min dxy's d per HDR vs simlated min per HDR
dxy.min.dens.era.w <-  ggplot(subset(dxy.era.sub, substr(cline, 8,8)=="w"), aes(x=dxy_mean, shape=cline)) + 
  scale_color_manual(values=c("#0000ff", "#5ec55b","#23ff06"), labels=c("High", "Low", "Low distant"))+
  labs(colour="Altitude")+
  geom_vline(xintercept = 0, colour="lightgrey", lty="dashed")+
  geom_vline(xintercept = subset(dxy.era.sub.para.west.summ.summ, is.allo.shdr=="no")$dxy.min.mean, colour="#0372B2", lty="solid", alpha=.9, size=0.2)+
  geom_vline(xintercept = subset(dxy.era.sub.para.west.summ.summ, is.allo.shdr=="yes")$dxy.min.mean, colour="#D65D00", lty="solid", alpha=.9, size=0.2)+
  geom_vline(xintercept = mean(subset(sim.dxy.era.para.west.summ.summ.df)$dxy.min.mean,na.rm=T), colour="black", lty="solid", alpha=.9, size=1.5)+
  geom_vline(xintercept = quantile(subset(sim.dxy.era.para.west.summ.summ.df)$dxy.min.mean, 0.9,na.rm=T), colour="black", lty="dashed", alpha=.9, size=1.5)+
  geom_vline(xintercept = quantile(subset(sim.dxy.era.para.west.summ.summ.df)$dxy.min.mean, 0.1,na.rm=T), colour="black", lty="dashed", alpha=.9, size=1.5)+
  geom_density(fill="transparent", alpha=0, size=1, colour="blue")+
  xlab(expression(paste("Tajima's D",  italic(" H. erato"))))+
  theme_classic()+
  theme(axis.line.x = element_line(color="black", size = 0.5),axis.line.y = element_line(color="black", size = 0.5),axis.title.y = element_text(size=14),
        plot.margin = unit(c(0.5,0.2,0.2,0.2), "lines"),axis.text = element_text(size=10), #axis.text.x = element_blank(),
        axis.text.y = element_text(size=10),axis.title.x = element_blank(), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_rect(fill = "transparent"), plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
        legend.background = element_rect(fill = "transparent"), legend.box.background = element_rect(fill = "transparent"), strip.background.y = element_blank(),
        strip.text = element_text(size=14))+
  scale_y_continuous(expand = c(0,0))+
  scale_x_continuous(expand = c(0,0)); dxy.min.dens.era.w

dxy.max.dens.era.w <-  ggplot(subset(dxy.era.sub, substr(cline, 8,8)=="w"), aes(x=dxy_mean, shape=cline)) + 
  scale_color_manual(values=c("#0000ff", "#5ec55b","#23ff06"), labels=c("High", "Low", "Low distant"))+
  labs(colour="Altitude")+
  geom_vline(xintercept = 0, colour="lightgrey", lty="dashed")+
  geom_vline(xintercept = subset(dxy.era.sub.para.west.summ.summ, is.allo.shdr=="no")$dxy.max.mean, colour="#0372B2", lty="solid", alpha=.9, size=0.2)+
  geom_vline(xintercept = subset(dxy.era.sub.para.west.summ.summ, is.allo.shdr=="yes")$dxy.max.mean, colour="#D65D00", lty="solid", alpha=.9, size=0.2)+
  geom_vline(xintercept = mean(subset(sim.dxy.era.para.west.summ.summ.df)$dxy.max.mean,na.rm=T), colour="black", lty="solid", alpha=.9, size=1.5)+
  geom_vline(xintercept = quantile(subset(sim.dxy.era.para.west.summ.summ.df)$dxy.max.mean, 0.9,na.rm=T), colour="black", lty="dashed", alpha=.9, size=1.5)+
  geom_vline(xintercept = quantile(subset(sim.dxy.era.para.west.summ.summ.df)$dxy.max.mean, 0.1,na.rm=T), colour="black", lty="dashed", alpha=.9, size=1.5)+
  geom_density(fill="transparent", alpha=0, size=1, colour="blue")+
  xlab(expression(paste("Tajima's D",  italic(" H. erato"))))+
  theme_classic()+
  theme(axis.line.x = element_line(color="black", size = 0.5),axis.line.y = element_line(color="black", size = 0.5),axis.title.y = element_text(size=14),
        plot.margin = unit(c(0.5,0.2,0.2,0.2), "lines"),axis.text = element_text(size=10), #axis.text.x = element_blank(),
        axis.text.y = element_text(size=10),axis.title.x = element_blank(), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_rect(fill = "transparent"), plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
        legend.background = element_rect(fill = "transparent"), legend.box.background = element_rect(fill = "transparent"), strip.background.y = element_blank(),
        strip.text = element_text(size=14))+
  scale_y_continuous(expand = c(0,0))+
  scale_x_continuous(expand = c(0,0)); dxy.max.dens.era.w




################# 4.4 which SHDR are dxy outliers? ######
library(patchwork)
dxy.min.dens.era.w/dxy.min.dens.era.e
ggsave("local/plots/dxy/dxy.era.min.per.HDR.90-10thpercentile.hig.vlo.dens.east.west.png", width = 5, height = 5)

dxy.max.dens.era.w/dxy.max.dens.era.e
ggsave("local/plots/dxy/dxy.era.max.per.HDR.90-10thpercentile.hig.vlo.dens.east.west.png", width = 5, height = 5)

# check thresholds for min dxy per HDR and % windows with dxy <5th perc
nrow(subset(dxy.era.sub.para.east.summ.summ, dxy.min.mean < quantile(subset(sim.dxy.era.para.east.summ.summ.df)$dxy.min.mean, 0.1, na.rm = T)))
nrow(subset(dxy.era.sub.para.east.summ.summ, dxy.5thperc.perc.windows.mean > quantile(subset(sim.dxy.era.para.east.summ.summ.df)$dxy.5thperc.perc.windows.mean, 0.9, na.rm = T)))

nrow(subset(dxy.era.sub.para.west.summ.summ, dxy.min.mean < quantile(subset(sim.dxy.era.para.west.summ.summ.df)$dxy.min.mean, 0.1, na.rm = T)))
nrow(subset(dxy.era.sub.para.west.summ.summ, dxy.5thperc.perc.windows.mean > quantile(subset(sim.dxy.era.para.west.summ.summ.df)$dxy.5thperc.perc.windows.mean, 0.9, na.rm = T)))

# max dxy
nrow(subset(dxy.era.sub.para.east.summ.summ, dxy.max.mean > quantile(subset(sim.dxy.era.para.east.summ.summ.df)$dxy.max.mean, 0.9, na.rm = T)))
nrow(subset(dxy.era.sub.para.east.summ.summ, dxy.5thperc.perc.windows.mean > quantile(subset(sim.dxy.era.para.east.summ.summ.df)$dxy.5thperc.perc.windows.mean, 0.9, na.rm = T)))

nrow(subset(dxy.era.sub.para.west.summ.summ, dxy.max.mean > quantile(subset(sim.dxy.era.para.west.summ.summ.df)$dxy.max.mean, 0.9, na.rm = T)))
nrow(subset(dxy.era.sub.para.west.summ.summ, dxy.5thperc.perc.windows.mean > quantile(subset(sim.dxy.era.para.west.summ.summ.df)$dxy.5thperc.perc.windows.mean, 0.9, na.rm = T)))


## add whether SHDRs are outliers or not
dxy.era.sub.para.east.summ.summ$dxy.min.less.10th.perc.sims <- if_else(dxy.era.sub.para.east.summ.summ$dxy.min.mean < quantile(sim.dxy.era.para.east.summ.summ.df$dxy.min.mean, 0.1,na.rm=T), "yes", "no")
dxy.era.sub.para.east.summ.summ$dxy.max.more.90th.perc.sims <- if_else(dxy.era.sub.para.east.summ.summ$dxy.max.mean > quantile(sim.dxy.era.para.east.summ.summ.df$dxy.max.mean, 0.9,na.rm=T), "yes", "no")
dxy.era.sub.para.east.summ.summ$dxy.5thperc.perc.windows.mean.more.90th.perc.sims <- if_else(dxy.era.sub.para.east.summ.summ$dxy.5thperc.perc.windows.mean > quantile(sim.dxy.era.para.east.summ.summ.df$dxy.5thperc.perc.windows.mean, 0.9,na.rm=T), "yes", "no")
dxy.era.sub.para.east.summ.summ$dxy.95thperc.perc.windows.mean.more.90th.perc.sims <- if_else(dxy.era.sub.para.east.summ.summ$dxy.95thperc.perc.windows.mean > quantile(sim.dxy.era.para.east.summ.summ.df$dxy.95thperc.perc.windows.mean, 0.9,na.rm=T), "yes", "no")


# west
dxy.era.sub.para.west.summ.summ$dxy.min.less.10th.perc.sims <- if_else(dxy.era.sub.para.west.summ.summ$dxy.min.mean < quantile(sim.dxy.era.para.west.summ.summ.df$dxy.min.mean, 0.1,na.rm=T), "yes", "no")
dxy.era.sub.para.west.summ.summ$dxy.max.more.90th.perc.sims <- if_else(dxy.era.sub.para.west.summ.summ$dxy.max.mean > quantile(sim.dxy.era.para.west.summ.summ.df$dxy.max.mean, 0.9,na.rm=T), "yes", "no")

dxy.era.sub.para.west.summ.summ$dxy.5thperc.perc.windows.mean.more.95th.perc.sims <- if_else(dxy.era.sub.para.west.summ.summ$dxy.5thperc.perc.windows.mean > quantile(sim.dxy.era.para.west.summ.summ.df$dxy.5thperc.perc.windows.mean, 0.95,na.rm=T), "yes", "no")
dxy.era.sub.para.west.summ.summ$dxy.5thperc.perc.windows.mean.more.90th.perc.sims <- if_else(dxy.era.sub.para.west.summ.summ$dxy.5thperc.perc.windows.mean > quantile(sim.dxy.era.para.west.summ.summ.df$dxy.5thperc.perc.windows.mean, 0.9,na.rm=T), "yes", "no")

dxy.era.sub.para.west.summ.summ$dxy.95thperc.perc.windows.mean.more.95th.perc.sims <- if_else(dxy.era.sub.para.west.summ.summ$dxy.95thperc.perc.windows.mean > quantile(sim.dxy.era.para.west.summ.summ.df$dxy.95thperc.perc.windows.mean, 0.95,na.rm=T), "yes", "no")
dxy.era.sub.para.west.summ.summ$dxy.95thperc.perc.windows.mean.more.90th.perc.sims <- if_else(dxy.era.sub.para.west.summ.summ$dxy.95thperc.perc.windows.mean > quantile(sim.dxy.era.para.west.summ.summ.df$dxy.95thperc.perc.windows.mean, 0.9,na.rm=T), "yes", "no")


# being a dxy max outlier is not independent from being a dxy min outlier
table(dxy.era.sub.para.west.summ.summ$dxy.min.less.10th.perc.sims,dxy.era.sub.para.west.summ.summ$dxy.max.more.90th.perc.sims, dnn = c("dxy min outliers", "dxy max outliers") )
chisq.test(dxy.era.sub.para.west.summ.summ$dxy.min.less.10th.perc.sims,dxy.era.sub.para.west.summ.summ$dxy.max.more.90th.perc.sims )

dxy.era.sub.para.west.summ.summ$dxy.95thperc.perc.windows.mean.more.95th.perc.sims <- if_else(dxy.era.sub.para.west.summ.summ$dxy.95thperc.perc.windows.mean > quantile(sim.dxy.era.para.west.summ.summ.df$dxy.95thperc.perc.windows.mean, 0.95,na.rm=T), "yes", "no")
dxy.era.sub.para.west.summ.summ$dxy.95thperc.perc.windows.mean.more.90th.perc.sims <- if_else(dxy.era.sub.para.west.summ.summ$dxy.95thperc.perc.windows.mean > quantile(sim.dxy.era.para.west.summ.summ.df$dxy.95thperc.perc.windows.mean, 0.9,na.rm=T), "yes", "no")

# save
dxy.era.sub.para.east.summ.summ; write.csv(dxy.era.sub.para.east.summ.summ, "local/data/dxy/shdr.dxy.era.sub.para.east.summ.summ.csv", row.names = F)
dxy.era.sub.para.west.summ.summ; write.csv(dxy.era.sub.para.west.summ.summ, "local/data/dxy/shdr.dxy.era.sub.para.west.summ.summ.csv", row.names = F)









######################################## 5. SUMMARIES PER HDR ####################
# load summaries 
thetas.era.hig.para.east.summ.summ<- read.csv("local/data/thetas/shdr.thetas.era.hig.para.east.summ.summ.csv")
thetas.era.hig.para.west.summ.summ<- read.csv("local/data/thetas/shdr.thetas.era.hig.para.west.summ.summ.csv")

pi.era.all.alt.wide.para.east.summ.summ<- read.csv("local/data/thetas/shdr.pi.era.all.alt.wide.para.east.summ.summ.csv")
pi.era.all.alt.wide.para.west.summ.summ<- read.csv("local/data/thetas/shdr.pi.era.all.alt.wide.para.west.summ.summ.csv")

dxy.era.sub.para.east.summ.summ<- read.csv("local/data/dxy/shdr.dxy.era.sub.para.east.summ.summ.csv")
dxy.era.sub.para.west.summ.summ<- read.csv("local/data/dxy/shdr.dxy.era.sub.para.west.summ.summ.csv")


################# 5.1 create SHDR EAST df with tajima, deltapi, dxy ######


names(dxy.era.sub.para.east.summ.summ)
shdr.era.east.outlier.df <- merge(thetas.era.hig.para.east.summ.summ[,c("shdr.para.east.id","scaff"  , "no.windows", "start", "end", "size.hdr.bp" ,
                                                                        "tajima.min.mean", "tajima.5thperc.mean.taj.mean", "tajima.5thperc.perc.windows.mean",
                                                                        "tajima.min.less.5th.perc.sims","tajima.min.less.10th.perc.sims","tajima.5thperc.perc.windows.mean.less.95th.perc.sims",
                                                                        "tajima.5thperc.perc.windows.mean.less.90th.perc.sims")], 
      pi.era.all.alt.wide.para.east.summ.summ[,c("shdr.para.east.id", "delta.pi.hig.vlo.min.mean","delta.pi.hig.vlo.5thperc.mean.delta.pi.mean","delta.pi.hig.vlo.5thperc.perc.windows.mean",
                                                 "delta.pi.hig.vlo.min.less.5th.perc.sims","delta.pi.hig.vlo.min.less.10th.perc.sims","delta.pi.hig.vlo.5thperc.perc.windows.mean.less.95th.perc.sims", 
                                                 "delta.pi.hig.vlo.5thperc.perc.windows.mean.less.90th.perc.sims")], by = "shdr.para.east.id" ) %>% merge(
                                                   dxy.era.sub.para.east.summ.summ[,c("shdr.para.east.id","dxy.min.mean", "dxy.max.mean","dxy.5thperc.mean.dxy.mean", "dxy.5thperc.perc.windows.mean"                    
                                         ,"dxy.95thperc.mean.dxy.mean", "dxy.95thperc.perc.windows.mean" , "dxy.min.less.5th.perc.sims","dxy.min.less.10th.perc.sims" ,"dxy.max.more.95th.perc.sims",                      
                                         "dxy.max.more.90th.perc.sims","dxy.5thperc.perc.windows.mean.more.95th.perc.sims","dxy.5thperc.perc.windows.mean.more.90th.perc.sims",
                                         "dxy.95thperc.perc.windows.mean.more.95th.perc.sims", "dxy.95thperc.perc.windows.mean.more.90th.perc.sims")], by = "shdr.para.east.id"); head(shdr.era.east.outlier.df )

# add allopatric id if shared allo
shdr.era.east.outlier.df$is.allo.hdr <- if_else(shdr.era.east.outlier.df$shdr.para.east.id %in% era.shdr.allo.summ$shdr.para.east, "yes" , "no")
shdr.era.east.outlier.df$shdr.allo.id <- era.shdr.allo.summ$shdr.allo.id[match(shdr.era.east.outlier.df$shdr.para.east.id, era.shdr.allo.summ$shdr.para.east)]
shdr.era.east.outlier.df[is.na(shdr.era.east.outlier.df$shdr.allo.id),]$shdr.allo.id <- NA
library(forcats); shdr.era.east.outlier.df$shdr.allo.id <- fct_explicit_na(shdr.era.east.outlier.df$shdr.allo.id, "no"); head(shdr.era.east.outlier.df)

# add bp pos
shdr.era.east.outlier.df$start.bp <-era.shdr.para.east$start.bp[match(shdr.era.east.outlier.df$shdr.para.east.id, era.shdr.para.east$shdr.para.east.id)]
shdr.era.east.outlier.df$end.bp <-era.shdr.para.east$end.bp[match(shdr.era.east.outlier.df$shdr.para.east.id, era.shdr.para.east$shdr.para.east.id)]

################# 5.2 create SHDR west df with tajima, deltapi, dxy ######
shdr.era.west.outlier.df <- merge(thetas.era.hig.para.west.summ.summ[,c("shdr.para.west.id","scaff"  , "no.windows", "start", "end", "size.hdr.bp" ,
                                                                        "tajima.min.mean", "tajima.5thperc.mean.taj.mean", "tajima.5thperc.perc.windows.mean",
                                                                        "tajima.min.less.5th.perc.sims","tajima.min.less.10th.perc.sims","tajima.5thperc.perc.windows.mean.less.95th.perc.sims",
                                                                        "tajima.5thperc.perc.windows.mean.less.90th.perc.sims")], 
                                  pi.era.all.alt.wide.para.west.summ.summ[,c("shdr.para.west.id", "delta.pi.hig.vlo.min.mean","delta.pi.hig.vlo.5thperc.mean.delta.pi.mean","delta.pi.hig.vlo.5thperc.perc.windows.mean",
                                                                             "delta.pi.hig.vlo.min.less.5th.perc.sims","delta.pi.hig.vlo.min.less.10th.perc.sims","delta.pi.hig.vlo.5thperc.perc.windows.mean.less.95th.perc.sims", 
                                                                             "delta.pi.hig.vlo.5thperc.perc.windows.mean.less.90th.perc.sims")], by = "shdr.para.west.id" ) %>% merge(
                                                                               dxy.era.sub.para.west.summ.summ[,c("shdr.para.west.id","dxy.min.mean", "dxy.max.mean","dxy.5thperc.mean.dxy.mean", "dxy.5thperc.perc.windows.mean"                    
                                                                                                                  ,"dxy.95thperc.mean.dxy.mean", "dxy.95thperc.perc.windows.mean" , "dxy.min.less.5th.perc.sims","dxy.min.less.10th.perc.sims" ,"dxy.max.more.95th.perc.sims",                      
                                                                                                                  "dxy.max.more.90th.perc.sims","dxy.5thperc.perc.windows.mean.more.95th.perc.sims","dxy.5thperc.perc.windows.mean.more.90th.perc.sims",
                                                                                                                  "dxy.95thperc.perc.windows.mean.more.95th.perc.sims", "dxy.95thperc.perc.windows.mean.more.90th.perc.sims")], by = "shdr.para.west.id"); head(shdr.era.west.outlier.df )

# add allopatric id if shared allo
shdr.era.west.outlier.df$is.allo.hdr <- if_else(shdr.era.west.outlier.df$shdr.para.west.id %in% era.shdr.allo.summ$shdr.para.west, "yes" , "no")
shdr.era.west.outlier.df$shdr.allo.id <- era.shdr.allo.summ$shdr.allo.id[match(shdr.era.west.outlier.df$shdr.para.west.id, era.shdr.allo.summ$shdr.para.west)]
shdr.era.west.outlier.df[is.na(shdr.era.west.outlier.df$shdr.allo.id),]$shdr.allo.id <- NA
library(forcats); shdr.era.west.outlier.df$shdr.allo.id <- fct_explicit_na(shdr.era.west.outlier.df$shdr.allo.id, "no"); head(shdr.era.west.outlier.df)

# add bp pos
shdr.era.west.outlier.df$start.bp <-era.shdr.para.west$start.bp[match(shdr.era.west.outlier.df$shdr.para.west.id, era.shdr.para.west$shdr.para.west.id)]
shdr.era.west.outlier.df$end.bp <-era.shdr.para.west$end.bp[match(shdr.era.west.outlier.df$shdr.para.west.id, era.shdr.para.west$shdr.para.west.id)]

################# 5.3 venn diagrams max/min ##########
######  prep ######
# this worked https://stackoverflow.com/questions/50997636/problems-installing-rgeos-and-rgdal-on-mac-os-x-high-sierra
# install.packages('sp', type='source')
# install.packages("rgeos", repos="http://R-Forge.R-project.org", type="source")
# require(rgeos)
# install.packages("rgdal", repos="http://R-Forge.R-project.org", type="source")
# require(rgdal)
# libraries <- c("rgdal","rgeos","raster","spdep","spatstat")
# install.packages(libraries)
# install.packages("ggVennDiagram")
library(ggVennDiagram)

########### 10th-90th percentiles of simulations #######

###### east all shdr: min tajima, min deltapi (i.e. lower pi at hig), max dxy- <10th percentile or >90th ######
tajima.min.less.10th.perc.sims.east.shdr <- as.character(subset(shdr.era.east.outlier.df, tajima.min.less.10th.perc.sims=="yes")$shdr.para.east.id)
delta.pi.hig.vlo.min.less.10th.perc.sims.east.shdr <- as.character(subset(shdr.era.east.outlier.df, delta.pi.hig.vlo.min.less.10th.perc.sims=="yes")$shdr.para.east.id)
dxy.max.more.90th.perc.sims.east.shdr <- as.character(subset(shdr.era.east.outlier.df, dxy.max.more.90th.perc.sims=="yes")$shdr.para.east.id)
outliers.east <- list(tajima=tajima.min.less.10th.perc.sims.east.shdr, delta.pi=delta.pi.hig.vlo.min.less.10th.perc.sims.east.shdr, dxy=dxy.max.more.90th.perc.sims.east.shdr)
str(outliers.east)

east.all.stats <- ggVennDiagram(outliers.east, color="black")+
  ggplot2::scale_fill_gradient(low="white",high = "black")+
  ggplot2::scale_color_manual(values=c("black","black","black")); east.all.stats

dxy.min.less.10th.perc.sims.east.shdr <- as.character(subset(shdr.era.east.outlier.df, dxy.min.less.10th.perc.sims=="yes")$shdr.para.east.id)
outliers.east <- list(tajima=tajima.min.less.10th.perc.sims.east.shdr, delta.pi=delta.pi.hig.vlo.min.less.10th.perc.sims.east.shdr, 
                      dxy.max=dxy.max.more.90th.perc.sims.east.shdr, dxy.min=dxy.min.less.10th.perc.sims.east.shdr)

east.all.stats <- ggVennDiagram(outliers.east, color="black")+
  ggplot2::scale_fill_gradient(low="white",high = "black")+
  ggplot2::scale_color_manual(values=c("black","black","black","black")); east.all.stats


###### east para shdr: min tajima, min deltapi (i.e. lower pi at hig), max dxy- <10th percentile or >90th ######
tajima.min.less.10th.perc.sims.east.para.shdr <- as.character(subset(shdr.era.east.outlier.df, tajima.min.less.10th.perc.sims=="yes"&is.allo.hdr=="no")$shdr.para.east.id)
delta.pi.hig.vlo.min.less.10th.perc.sims.east.para.shdr <- as.character(subset(shdr.era.east.outlier.df, delta.pi.hig.vlo.min.less.10th.perc.sims=="yes"&is.allo.hdr=="no")$shdr.para.east.id)
dxy.max.more.90th.perc.sims.east.para.shdr <- as.character(subset(shdr.era.east.outlier.df, dxy.max.more.90th.perc.sims=="yes"&is.allo.hdr=="no")$shdr.para.east.id)
outliers.east.para <- list(tajima=tajima.min.less.10th.perc.sims.east.para.shdr, delta.pi=delta.pi.hig.vlo.min.less.10th.perc.sims.east.para.shdr, dxy=dxy.max.more.90th.perc.sims.east.para.shdr)

east.para.stats <- ggVennDiagram(outliers.east.para, color="black")+
  ggplot2::scale_fill_gradient(low="white",high = "#049E73")+
  ggplot2::scale_color_manual(values=c("black","black","black")); east.para.stats

###### east allo shdr: min tajima, min deltapi (i.e. lower pi at hig), max dxy- <10th percentile or >90th ######
tajima.min.less.10th.perc.sims.east.allo.shdr <- as.character(subset(shdr.era.east.outlier.df, tajima.min.less.10th.perc.sims=="yes"&is.allo.hdr=="yes")$shdr.para.east.id); tajima.min.less.10th.perc.sims.east.allo.shdr
delta.pi.hig.vlo.min.less.10th.perc.sims.east.allo.shdr <- as.character(subset(shdr.era.east.outlier.df, delta.pi.hig.vlo.min.less.10th.perc.sims=="yes"&is.allo.hdr=="yes")$shdr.para.east.id)
dxy.max.more.90th.perc.sims.east.allo.shdr <- as.character(subset(shdr.era.east.outlier.df, dxy.max.more.90th.perc.sims=="yes"&is.allo.hdr=="yes")$shdr.para.east.id)
outliers.east.allo <- list(tajima=tajima.min.less.10th.perc.sims.east.allo.shdr, delta.pi=delta.pi.hig.vlo.min.less.10th.perc.sims.east.allo.shdr, dxy=dxy.max.more.90th.perc.sims.east.allo.shdr)

east.allo.stats <- ggVennDiagram(outliers.east.allo, color="black")+
  ggplot2::scale_fill_gradient(low="white",high = "#D65D00")+
  ggplot2::scale_color_manual(values=c("black","black","black")); east.allo.stats

library(patchwork)
east.all.stats | east.para.stats |  east.allo.stats

###### west all shdr: min tajima, min deltapi (i.e. lower pi at hig), max dxy- <10th percentile or >90th ######
tajima.min.less.10th.perc.sims.west.shdr <- as.character(subset(shdr.era.west.outlier.df, tajima.min.less.10th.perc.sims=="yes")$shdr.para.west.id)
delta.pi.hig.vlo.min.less.10th.perc.sims.west.shdr <- as.character(subset(shdr.era.west.outlier.df, delta.pi.hig.vlo.min.less.10th.perc.sims=="yes")$shdr.para.west.id)
dxy.max.more.90th.perc.sims.west.shdr <- as.character(subset(shdr.era.west.outlier.df, dxy.max.more.90th.perc.sims=="yes")$shdr.para.west.id)
outliers.west <- list(tajima=tajima.min.less.10th.perc.sims.west.shdr, delta.pi=delta.pi.hig.vlo.min.less.10th.perc.sims.west.shdr, dxy=dxy.max.more.90th.perc.sims.west.shdr)
str(outliers.west)

west.all.stats <- ggVennDiagram(outliers.west, color="black")+
  ggplot2::scale_fill_gradient(low="white",high = "black")+
  ggplot2::scale_color_manual(values=c("black","black","black")); west.all.stats

###### west para shdr: min tajima, min deltapi (i.e. lower pi at hig), max dxy- <10th percentile or >90th ######
tajima.min.less.10th.perc.sims.west.para.shdr <- as.character(subset(shdr.era.west.outlier.df, tajima.min.less.10th.perc.sims=="yes"&is.allo.hdr=="no")$shdr.para.west.id)
delta.pi.hig.vlo.min.less.10th.perc.sims.west.para.shdr <- as.character(subset(shdr.era.west.outlier.df, delta.pi.hig.vlo.min.less.10th.perc.sims=="yes"&is.allo.hdr=="no")$shdr.para.west.id)
dxy.max.more.90th.perc.sims.west.para.shdr <- as.character(subset(shdr.era.west.outlier.df, dxy.max.more.90th.perc.sims=="yes"&is.allo.hdr=="no")$shdr.para.west.id)
outliers.west.para <- list(tajima=tajima.min.less.10th.perc.sims.west.para.shdr, delta.pi=delta.pi.hig.vlo.min.less.10th.perc.sims.west.para.shdr, dxy=dxy.max.more.90th.perc.sims.west.para.shdr)

west.para.stats <- ggVennDiagram(outliers.west.para, color="black")+
  ggplot2::scale_fill_gradient(low="white",high = "#0372B2")+
  ggplot2::scale_color_manual(values=c("black","black","black")); west.para.stats

###### west allo shdr: min tajima, min deltapi (i.e. lower pi at hig), max dxy- <10th percentile or >90th ######
tajima.min.less.10th.perc.sims.west.allo.shdr <- as.character(subset(shdr.era.west.outlier.df, tajima.min.less.10th.perc.sims=="yes"&is.allo.hdr=="yes")$shdr.para.west.id); tajima.min.less.10th.perc.sims.west.allo.shdr
delta.pi.hig.vlo.min.less.10th.perc.sims.west.allo.shdr <- as.character(subset(shdr.era.west.outlier.df, delta.pi.hig.vlo.min.less.10th.perc.sims=="yes"&is.allo.hdr=="yes")$shdr.para.west.id)
dxy.max.more.90th.perc.sims.west.allo.shdr <- as.character(subset(shdr.era.west.outlier.df, dxy.max.more.90th.perc.sims=="yes"&is.allo.hdr=="yes")$shdr.para.west.id)
outliers.west.allo <- list(tajima=tajima.min.less.10th.perc.sims.west.allo.shdr, delta.pi=delta.pi.hig.vlo.min.less.10th.perc.sims.west.allo.shdr, dxy=dxy.max.more.90th.perc.sims.west.allo.shdr)

west.allo.stats <- ggVennDiagram(outliers.west.allo, color="black")+
  ggplot2::scale_fill_gradient(low="white",high = "#D65D00")+
  ggplot2::scale_color_manual(values=c("black","black","black")); west.allo.stats

library(patchwork)
(west.all.stats | west.para.stats |  west.allo.stats) / (east.all.stats | east.para.stats |  east.allo.stats)

ggsave("local/plots/thetas/venn.era.east.west.tajima.dxy.deltapi.10th.90th.ouliers.png", width = 12, height = 6)

################# 5.4 venn diagrams % of hdr with 5th outliers ##########

########### 10th-90th percentiles of simulations #######

###### east all shdr: min tajima, min deltapi (i.e. lower pi at hig), max dxy- <10th percentile or >90th ######
names(shdr.era.east.outlier.df)
tajima.5thperc.perc.windows.mean.less.90th.perc.sims.east.shdr <- as.character(subset(shdr.era.east.outlier.df, tajima.5thperc.perc.windows.mean.less.90th.perc.sims=="yes")$shdr.para.east.id)
delta.pi.hig.vlo.5thperc.perc.windows.mean.less.90th.perc.sims.east.shdr <- as.character(subset(shdr.era.east.outlier.df, delta.pi.hig.vlo.5thperc.perc.windows.mean.less.90th.perc.sims=="yes")$shdr.para.east.id)
dxy.95thperc.perc.windows.mean.more.90th.perc.sims.east.shdr <- as.character(subset(shdr.era.east.outlier.df, dxy.95thperc.perc.windows.mean.more.90th.perc.sims=="yes")$shdr.para.east.id)
outliers.east <- list(tajima=tajima.5thperc.perc.windows.mean.less.90th.perc.sims.east.shdr, delta.pi=delta.pi.hig.vlo.5thperc.perc.windows.mean.less.90th.perc.sims.east.shdr, dxy=dxy.95thperc.perc.windows.mean.more.90th.perc.sims.east.shdr)
str(outliers.east)

east.all.stats <- ggVennDiagram(outliers.east, color="black")+
  ggplot2::scale_fill_gradient(low="white",high = "black")+
  ggplot2::scale_color_manual(values=c("black","black","black")); east.all.stats

###### east para shdr: min tajima, min deltapi (i.e. lower pi at hig), max dxy- <10th percentile or >90th ######
tajima.5thperc.perc.windows.mean.less.90th.perc.sims.east.para.shdr <- as.character(subset(shdr.era.east.outlier.df, tajima.5thperc.perc.windows.mean.less.90th.perc.sims=="yes"&is.allo.hdr=="no")$shdr.para.east.id)
delta.pi.hig.vlo.5thperc.perc.windows.mean.less.90th.perc.sims.east.para.shdr <- as.character(subset(shdr.era.east.outlier.df, delta.pi.hig.vlo.5thperc.perc.windows.mean.less.90th.perc.sims=="yes"&is.allo.hdr=="no")$shdr.para.east.id)
dxy.95thperc.perc.windows.mean.more.90th.perc.sims.east.para.shdr <- as.character(subset(shdr.era.east.outlier.df, dxy.95thperc.perc.windows.mean.more.90th.perc.sims=="yes"&is.allo.hdr=="no")$shdr.para.east.id)
outliers.east.para <- list(tajima=tajima.5thperc.perc.windows.mean.less.90th.perc.sims.east.para.shdr, delta.pi=delta.pi.hig.vlo.5thperc.perc.windows.mean.less.90th.perc.sims.east.para.shdr, dxy=dxy.95thperc.perc.windows.mean.more.90th.perc.sims.east.para.shdr)

east.para.stats <- ggVennDiagram(outliers.east.para, color="black")+
  ggplot2::scale_fill_gradient(low="white",high = "#049E73")+
  ggplot2::scale_color_manual(values=c("black","black","black")); east.para.stats

###### east allo shdr: min tajima, min deltapi (i.e. lower pi at hig), max dxy- <10th percentile or >90th ######
tajima.5thperc.perc.windows.mean.less.90th.perc.sims.east.allo.shdr <- as.character(subset(shdr.era.east.outlier.df, tajima.5thperc.perc.windows.mean.less.90th.perc.sims=="yes"&is.allo.hdr=="yes")$shdr.para.east.id); tajima.5thperc.perc.windows.mean.less.90th.perc.sims.east.allo.shdr
delta.pi.hig.vlo.5thperc.perc.windows.mean.less.90th.perc.sims.east.allo.shdr <- as.character(subset(shdr.era.east.outlier.df, delta.pi.hig.vlo.5thperc.perc.windows.mean.less.90th.perc.sims=="yes"&is.allo.hdr=="yes")$shdr.para.east.id)
dxy.95thperc.perc.windows.mean.more.90th.perc.sims.east.allo.shdr <- as.character(subset(shdr.era.east.outlier.df, dxy.95thperc.perc.windows.mean.more.90th.perc.sims=="yes"&is.allo.hdr=="yes")$shdr.para.east.id)
outliers.east.allo <- list(tajima=tajima.5thperc.perc.windows.mean.less.90th.perc.sims.east.allo.shdr, delta.pi=delta.pi.hig.vlo.5thperc.perc.windows.mean.less.90th.perc.sims.east.allo.shdr, dxy=dxy.95thperc.perc.windows.mean.more.90th.perc.sims.east.allo.shdr)

east.allo.stats <- ggVennDiagram(outliers.east.allo, color="black")+
  ggplot2::scale_fill_gradient(low="white",high = "#D65D00")+
  ggplot2::scale_color_manual(values=c("black","black","black")); east.allo.stats

library(patchwork)
east.all.stats | east.para.stats |  east.allo.stats

###### west all shdr: min tajima, min deltapi (i.e. lower pi at hig), max dxy- <10th percentile or >90th ######
tajima.5thperc.perc.windows.mean.less.90th.perc.sims.west.shdr <- as.character(subset(shdr.era.west.outlier.df, tajima.5thperc.perc.windows.mean.less.90th.perc.sims=="yes")$shdr.para.west.id)
delta.pi.hig.vlo.5thperc.perc.windows.mean.less.90th.perc.sims.west.shdr <- as.character(subset(shdr.era.west.outlier.df, delta.pi.hig.vlo.5thperc.perc.windows.mean.less.90th.perc.sims=="yes")$shdr.para.west.id)
dxy.95thperc.perc.windows.mean.more.90th.perc.sims.west.shdr <- as.character(subset(shdr.era.west.outlier.df, dxy.95thperc.perc.windows.mean.more.90th.perc.sims=="yes")$shdr.para.west.id)
outliers.west <- list(tajima=tajima.5thperc.perc.windows.mean.less.90th.perc.sims.west.shdr, delta.pi=delta.pi.hig.vlo.5thperc.perc.windows.mean.less.90th.perc.sims.west.shdr, dxy=dxy.95thperc.perc.windows.mean.more.90th.perc.sims.west.shdr)
str(outliers.west)

west.all.stats <- ggVennDiagram(outliers.west, color="black")+
  ggplot2::scale_fill_gradient(low="white",high = "black")+
  ggplot2::scale_color_manual(values=c("black","black","black")); west.all.stats

###### west para shdr: min tajima, min deltapi (i.e. lower pi at hig), max dxy- <10th percentile or >90th ######
tajima.5thperc.perc.windows.mean.less.90th.perc.sims.west.para.shdr <- as.character(subset(shdr.era.west.outlier.df, tajima.5thperc.perc.windows.mean.less.90th.perc.sims=="yes"&is.allo.hdr=="no")$shdr.para.west.id)
delta.pi.hig.vlo.5thperc.perc.windows.mean.less.90th.perc.sims.west.para.shdr <- as.character(subset(shdr.era.west.outlier.df, delta.pi.hig.vlo.5thperc.perc.windows.mean.less.90th.perc.sims=="yes"&is.allo.hdr=="no")$shdr.para.west.id)
dxy.95thperc.perc.windows.mean.more.90th.perc.sims.west.para.shdr <- as.character(subset(shdr.era.west.outlier.df, dxy.95thperc.perc.windows.mean.more.90th.perc.sims=="yes"&is.allo.hdr=="no")$shdr.para.west.id)
outliers.west.para <- list(tajima=tajima.5thperc.perc.windows.mean.less.90th.perc.sims.west.para.shdr, delta.pi=delta.pi.hig.vlo.5thperc.perc.windows.mean.less.90th.perc.sims.west.para.shdr, dxy=dxy.95thperc.perc.windows.mean.more.90th.perc.sims.west.para.shdr)

west.para.stats <- ggVennDiagram(outliers.west.para, color="black")+
  ggplot2::scale_fill_gradient(low="white",high = "#0372B2")+
  ggplot2::scale_color_manual(values=c("black","black","black")); west.para.stats

###### west allo shdr: min tajima, min deltapi (i.e. lower pi at hig), max dxy- <10th percentile or >90th ######
tajima.5thperc.perc.windows.mean.less.90th.perc.sims.west.allo.shdr <- as.character(subset(shdr.era.west.outlier.df, tajima.5thperc.perc.windows.mean.less.90th.perc.sims=="yes"&is.allo.hdr=="yes")$shdr.para.west.id); tajima.5thperc.perc.windows.mean.less.90th.perc.sims.west.allo.shdr
delta.pi.hig.vlo.5thperc.perc.windows.mean.less.90th.perc.sims.west.allo.shdr <- as.character(subset(shdr.era.west.outlier.df, delta.pi.hig.vlo.5thperc.perc.windows.mean.less.90th.perc.sims=="yes"&is.allo.hdr=="yes")$shdr.para.west.id)
dxy.95thperc.perc.windows.mean.more.90th.perc.sims.west.allo.shdr <- as.character(subset(shdr.era.west.outlier.df, dxy.95thperc.perc.windows.mean.more.90th.perc.sims=="yes"&is.allo.hdr=="yes")$shdr.para.west.id)
outliers.west.allo <- list(tajima=tajima.5thperc.perc.windows.mean.less.90th.perc.sims.west.allo.shdr, delta.pi=delta.pi.hig.vlo.5thperc.perc.windows.mean.less.90th.perc.sims.west.allo.shdr, dxy=dxy.95thperc.perc.windows.mean.more.90th.perc.sims.west.allo.shdr)

west.allo.stats <- ggVennDiagram(outliers.west.allo, color="black")+
  ggplot2::scale_fill_gradient(low="white",high = "#D65D00")+
  ggplot2::scale_color_manual(values=c("black","black","black")); west.allo.stats

library(patchwork)
(west.all.stats | west.para.stats |  west.allo.stats) / (east.all.stats | east.para.stats |  east.allo.stats)

ggsave("local/plots/thetas/venn.era.east.west.tajima.dxy.deltapi.10th.90th.perc.windows.10th.90th.outliers.png", width = 12, height = 6)

################# 5.5 number and % of SHDR with outliers for 1,2,3 stats ################# 
#### prep ####
# summarise by whether is para/allo, and by the number of SHDR with 0/1/2/3 outlier stats
shdr.era.east.outlier.df$no.stats.outlier.min.max.10th.90th.sims <- if_else(shdr.era.east.outlier.df$tajima.min.less.10th.perc.sims=="yes", 1, 0) +
  if_else(shdr.era.east.outlier.df$delta.pi.hig.vlo.min.less.10th.perc.sims=="yes", 1, 0) +
  if_else(shdr.era.east.outlier.df$dxy.max.more.90th.perc.sims=="yes", 1, 0); shdr.era.east.outlier.df$no.stats.outlier.min.max.10th.90th.sims

shdr.era.east.outlier.df$no.stats.outlier.min.max.5th.95th.sims <- if_else(shdr.era.east.outlier.df$tajima.min.less.5th.perc.sims=="yes", 1, 0) +
  if_else(shdr.era.east.outlier.df$delta.pi.hig.vlo.min.less.5th.perc.sims=="yes", 1, 0) +
  if_else(shdr.era.east.outlier.df$dxy.max.more.95th.perc.sims=="yes", 1, 0); shdr.era.east.outlier.df$no.stats.outlier.min.max.5th.95th.sims

shdr.era.east.outlier.df.summ <- summarise(group_by(shdr.era.east.outlier.df, is.allo.hdr, no.stats.outlier.min.max.10th.90th.sims),
          n=n()); shdr.era.east.outlier.df.summ
shdr.era.east.no <- summarise(group_by(shdr.era.east.outlier.df, is.allo.hdr),
                                           n=n()); shdr.era.east.no
shdr.era.east.outlier.df.summ$total.shdr <- shdr.era.east.no$n[match(shdr.era.east.outlier.df.summ$is.allo.hdr, shdr.era.east.no$is.allo.hdr)]; shdr.era.east.outlier.df.summ

# summarise by whether is para/allo, and by the number of SHDR with 0/1/2/3 outlier stats
shdr.era.west.outlier.df$no.stats.outlier.min.max.10th.90th.sims <- if_else(shdr.era.west.outlier.df$tajima.min.less.10th.perc.sims=="yes", 1, 0) +
  if_else(shdr.era.west.outlier.df$delta.pi.hig.vlo.min.less.10th.perc.sims=="yes", 1, 0) +
  if_else(shdr.era.west.outlier.df$dxy.max.more.90th.perc.sims=="yes", 1, 0); shdr.era.west.outlier.df$no.stats.outlier.min.max.10th.90th.sims

shdr.era.west.outlier.df$no.stats.outlier.min.max.5th.95th.sims <- if_else(shdr.era.west.outlier.df$tajima.min.less.5th.perc.sims=="yes", 1, 0) +
  if_else(shdr.era.west.outlier.df$delta.pi.hig.vlo.min.less.5th.perc.sims=="yes", 1, 0) +
  if_else(shdr.era.west.outlier.df$dxy.max.more.95th.perc.sims=="yes", 1, 0); shdr.era.west.outlier.df$no.stats.outlier.min.max.5th.95th.sims

shdr.era.west.outlier.df.summ <- summarise(group_by(shdr.era.west.outlier.df, is.allo.hdr, no.stats.outlier.min.max.10th.90th.sims),
                                           n=n()); shdr.era.west.outlier.df.summ
shdr.era.west.no <- summarise(group_by(shdr.era.west.outlier.df, is.allo.hdr),
                              n=n()); shdr.era.west.no
shdr.era.west.outlier.df.summ$total.shdr <- shdr.era.west.no$n[match(shdr.era.west.outlier.df.summ$is.allo.hdr, shdr.era.west.no$is.allo.hdr)]; shdr.era.west.outlier.df.summ

shdr.era.east.outlier.df.summ$side <- "East"
shdr.era.west.outlier.df.summ$side <- "West"

shdr.era.outlier.df.summ <- rbind(shdr.era.east.outlier.df.summ, shdr.era.west.outlier.df.summ); shdr.era.outlier.df.summ
shdr.era.outlier.df.summ$is.allo.hdr.side <- paste(shdr.era.outlier.df.summ$is.allo.hdr, shdr.era.outlier.df.summ$side)

shdr.era.outlier.df.summ$side <- factor(shdr.era.outlier.df.summ$side, levels=c("West", "East"))
shdr.era.outlier.df.summ$is.allo.hdr <- factor(shdr.era.outlier.df.summ$is.allo.hdr, levels=c("no", "yes"), labels = c("Parapatric", "Allopatric"))
str(shdr.era.outlier.df.summ )

# add max pbs per region and its position
shdr.era.east.outlier.df <-  merge(shdr.era.east.outlier.df,era.shdr.para.east.max,by ="shdr.para.east.id" )
shdr.era.west.outlier.df <-  merge(shdr.era.west.outlier.df,era.shdr.para.west.max,by ="shdr.para.west.id" )
# to-do
shdr.mel.east.outlier.df <-  merge(shdr.mel.east.outlier.df,mel.shdr.para.east.max,by ="shdr.para.east.id" )
shdr.mel.west.outlier.df <-  merge(shdr.mel.west.outlier.df,mel.shdr.para.west.max,by ="shdr.para.west.id" )

# is max pbs or fst >4? SHOULD be!
shdr.era.east.outlier.df$is.max.pbs.hig.above4 <-if_else((shdr.era.east.outlier.df$zPBS.hig.co.e.max.value>=4 & shdr.era.east.outlier.df$zPBS.hig.ec.e.max.value>=4), "yes", "no" )
shdr.era.west.outlier.df$is.max.pbs.hig.above4 <-if_else((shdr.era.west.outlier.df$zFst.co.w.max.value>=4 & shdr.era.west.outlier.df$zPBS.hig.ec.w.max.value >=4), "yes", "no" )
nrow(subset(shdr.era.east.outlier.df, is.max.pbs.hig.above4=="no")); nrow(subset(shdr.era.west.outlier.df, is.max.pbs.hig.above4=="no"))

# save summaries
head(shdr.era.east.outlier.df); head(shdr.era.west.outlier.df)
write.csv(shdr.era.east.outlier.df, "local/data/shdr.summ/shdr.era.east.outlier.df.csv", row.names = F )
write.csv(shdr.era.west.outlier.df, "local/data/shdr.summ/shdr.era.west.outlier.df.csv", row.names = F )

write.csv(shdr.era.east.outlier.df, "local/data/shdr.summ/shdr.era.east.outlier.df.csv", row.names = F )
write.csv(shdr.era.west.outlier.df, "local/data/shdr.summ/shdr.era.west.outlier.df.csv", row.names = F )
write.csv(shdr.era.outlier.df.summ, "local/data/shdr.summ/shdr.era.outlier.df.summ.csv", row.names = F )



#### load prepped summs  ####
shdr.era.east.outlier.df <- read.csv("local/data/shdr.summ/shdr.era.east.outlier.df.csv"); head(shdr.era.east.outlier.df)
shdr.era.west.outlier.df <- read.csv("local/data/shdr.summ/shdr.era.west.outlier.df.csv"); head(shdr.era.west.outlier.df)
shdr.era.outlier.df.summ <- read.csv( "local/data/shdr.summ/shdr.era.outlier.df.summ.csv")
shdr.era.outlier.df.summ$side <- factor(shdr.era.outlier.df.summ$side, levels=c("West", "East"))
shdr.era.outlier.df.summ$is.allo.hdr <- factor(shdr.era.outlier.df.summ$is.allo.hdr, levels=c("Parapatric", "Allopatric"))


# what percentage of hdrs have 1 2 or 3
library(ggtext)
ggplot(data=shdr.era.outlier.df.summ, aes(fill=no.stats.outlier.min.max.10th.90th.sims, x=is.allo.hdr, y=n, label=n)) + 
  geom_col(position = "fill", aes(color=is.allo.hdr.side), size=1)+ 
  ylab("% of SHDR with selection statistic outliers") + 
  geom_textbox(width = unit(0.1, "npc"), halign=.5,valign=.5, position = position_fill(vjust = 0.5), fill = c("cornsilk")) +
  scale_color_manual(values=c( "#049E73","#0372B2","#D65D00","#D65D00" )) +
  scale_y_continuous( labels = function(x) paste0(x*100, "%"), expand = c(0.05,0))+
  scale_fill_continuous(low="white", high="black") +theme_bw()+facet_wrap(~side)+
  theme(axis.line.x = element_line(color="black", size = 0.5),axis.line.y = element_line(color="black", size = 0.5),axis.title.y = element_text(size=14),
        plot.margin = unit(c(0.5,0.2,0.2,0.2), "lines"),axis.text = element_text(size=14), #axis.text.x = element_blank(),
        axis.text.y = element_text(size=10),axis.title.x = element_blank(), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_rect(fill = "transparent"), plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
        legend.background = element_rect(fill = "transparent"), legend.box.background = element_rect(fill = "transparent"), strip.background.y = element_blank(),
        strip.text = element_text(size=14, face = "bold"), strip.background = element_blank(), legend.position = "none")
ggsave("local/plots/shdr.summ/era.e.w.outlier.min.max.10th.90th.sims.tajima.deltapi.dxy.barplot.png", height = 6, width = 6)

ggplot(data=shdr.era.outlier.df.summ, aes(fill=no.stats.outlier.min.max.10th.90th.sims, x=is.allo.hdr, y=n, label=total.shdr)) + 
  geom_col(aes(color=is.allo.hdr.side), size=1)+ 
  ylab("Number of SHDR with selection statistic outliers") + 
  scale_color_manual(values=c( "#049E73","#0372B2","#D65D00","#D65D00" )) +
  scale_y_continuous( expand = c(0,0), limits = c(0,100))+
  scale_fill_continuous(low="white", high="black") +theme_bw()+facet_wrap(~side)+
  theme(axis.line.x = element_line(color="black", size = 0.5),axis.line.y = element_line(color="black", size = 0.5),axis.title.y = element_text(size=14),
        plot.margin = unit(c(0.5,0.2,0.2,0.2), "lines"),axis.text = element_text(size=14), #axis.text.x = element_blank(),
        axis.text.y = element_text(size=10),axis.title.x = element_blank(), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_rect(fill = "transparent"), plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
        legend.background = element_rect(fill = "transparent"), legend.box.background = element_rect(fill = "transparent"), strip.background.y = element_blank(),
        strip.text = element_text(size=14, face = "bold"), strip.background = element_blank(), legend.position = "none")
ggsave("local/plots/shdr.summ/era.e.w.outlier.number.min.max.10th.90th.sims.tajima.deltapi.dxy.barplot.png", height = 6, width = 6)


