##### get sites per SHDR for local PCA #####
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



############################## prep gff ######
# load gff, get genes, add offsets, get intervals
ref.scaff.era <- read.table("local/data/ref/Heliconius_erato_demophoon_v1_-_scaffolds.fa.fai", row.names = NULL)
era.scafEnds <- cumsum(ref.scaff.era[,2]); era.offset <- era.scafEnds - ref.scaff.era[,2]; names(era.offset) <- ref.scaff.era[,1]
era.gff <- read.table(pipe("cut -f 1,3,4,5 local/data/ref/heliconius_erato_demophoon_v1_core_32_85_1.gff"),
                      sep = "\t", as.is = T, col.names = c("scaffold", "feature", "start", "end"))
era.gff_gene <- subset(era.gff, feature == "gene")
era.gff_gene$start.offset <- era.gff_gene$start + era.offset[era.gff_gene$scaffold]; era.gff_gene$end.offset <- era.gff_gene$end + era.offset[era.gff_gene$scaffold]
era.gff_gene$chr <- as.integer(substr(era.gff_gene$scaffold, 7,8))

############################## not intersects, any site >=4stdv ##########
########## ERATO data #####
setwd("git/2021-altitude-heliconius/")
ref.scaff.era <- read.table("local/data/ref/Heliconius_erato_demophoon_v1_-_scaffolds.fa.fai", row.names = NULL)
all.bams.list.info <- read.csv("02.info/all.bams.list.info.csv"); names(all.bams.list.info)
era.e.list.info <- subset(all.bams.list.info,to.use.pbs=="yes" &species=="erato"&side.short=="e" )
era.w.list.info <- subset(all.bams.list.info,to.use.pbs=="yes" &species=="erato"&side.short=="w" )

era.all.pbs <- read.csv("local/data/pbs.out/era.all.pbs.recomb.csv")
era.all.pbs$CHR <- as.numeric(as.character(substr(era.all.pbs$scaff,7,8))); unique(era.all.pbs$CHR)
names(era.all.pbs)
unique(era.all.pbs$shdr.para.east.id)

summarise(group_by(era.all.pbs, shdr.para.east.id),
          bp.min=min(midPos),
          bp.max=max(midPos))

##### era.e subset so that either clines have >=4stdv #####
era.all.pbs.within.shdr.and.above.4stdv <- subset(era.all.pbs, !(is.na(shdr.para.east.id) ) &( zPBS0.1>=4 | zPBS0.3>=4 ))
summarise(group_by(era.all.pbs.within.shdr.and.above.4stdv, shdr.para.east.id, scaff),
          bp.min=min(midPos),
          bp.max=max(midPos),
          BP.wg.min=min(BP.wg),
          BP.wg.max=max(BP.wg))

## expand 1000 each window or 5000, subset to uniques, export for pcangsd
names(era.all.pbs.within.shdr.and.above.4stdv)
era.all.pbs.within.shdr.and.above.4stdv.angsd <- era.all.pbs.within.shdr.and.above.4stdv[c("scaff","midPos","BP.wg", "shdr.para.east.id")]; head(era.all.pbs.within.shdr.and.above.4stdv)

# windows are 5kb
era.all.pbs.within.shdr.and.above.4stdv.angsd$start <-  as.numeric(as.character(era.all.pbs.within.shdr.and.above.4stdv.angsd$midPos))-2500; head(era.all.pbs.within.shdr.and.above.4stdv.angsd)
era.all.pbs.within.shdr.and.above.4stdv.angsd$end <-  as.numeric(as.character(era.all.pbs.within.shdr.and.above.4stdv.angsd$midPos))+2500; head(era.all.pbs.within.shdr.and.above.4stdv.angsd)

era.all.pbs.within.shdr.and.above.4stdv.angsd$start.BP.wg <-  as.numeric(as.character(era.all.pbs.within.shdr.and.above.4stdv.angsd$BP.wg))-2500; head(era.all.pbs.within.shdr.and.above.4stdv.angsd)
era.all.pbs.within.shdr.and.above.4stdv.angsd$end.BP.wg <-  as.numeric(as.character(era.all.pbs.within.shdr.and.above.4stdv.angsd$BP.wg))+2500; head(era.all.pbs.within.shdr.and.above.4stdv.angsd)

era.all.pbs.within.shdr.and.above.4stdv.angsd.summ <- era.all.pbs.within.shdr.and.above.4stdv.angsd %>% 
  rowwise() %>% 
  do(data.frame(scaff=.$scaff, midPos= .$midPos, pos = .$start:.$end, shdr.para.east.id=.$shdr.para.east.id, BP.wg=.$start.BP.wg:.$end.BP.wg )); era.all.pbs.within.shdr.and.above.4stdv.angsd.summ

era.e.summ <- summarise(group_by(era.all.pbs.within.shdr.and.above.4stdv.angsd.summ, shdr.para.east.id, scaff ),
                        bp.min=min(pos ),
                        bp.max=max(pos ),
                        region.size=bp.max-bp.min,
                        BP.wg.min=min(BP.wg),
                        BP.wg.max=max(BP.wg)); era.e.summ

write.csv(era.e.summ, "02.info/sites/shared.exp50kb.era.e.std4.new.era.e.summ.csv", row.names = F)

sites.shared.ids <- list(); regions.shared.ids <- list(); for (i in unique(era.all.pbs.within.shdr.and.above.4stdv.angsd.summ$shdr.para.east.id )) {
  sites.shared.ids[[i]] <- subset(era.all.pbs.within.shdr.and.above.4stdv.angsd.summ, shdr.para.east.id ==i)
  sites.shared.ids[[i]] <-  sites.shared.ids[[i]][c("scaff", "pos")]
  regions.shared.ids[[i]] <- unique(sites.shared.ids[[i]]$scaff)
  write.table(sites.shared.ids[[i]], paste0("02.info/sites/shared.exp50kb.era.e.std4/sites.shared.region.", i,".txt"), 
              col.names = F, row.names = F, quote =F )
  write.table(regions.shared.ids[[i]], paste0("02.info/sites/shared.exp50kb.era.e.std4/regions.shared.region.", i,".txt"), 
              col.names = F, row.names = F, quote =F )
}





##### era.w subset so that either clines have >=4stdv #####
era.all.pbs.within.shdr.w.and.above.4stdv <- subset(era.all.pbs, !(is.na(shdr.para.west.id) ) & (zfst.2>=4 | zPBS0.4>=4 ))
summarise(group_by(era.all.pbs.within.shdr.w.and.above.4stdv, shdr.para.west.id, scaff),
          bp.min=min(midPos),
          bp.max=max(midPos),
          BP.wg.min=min(BP.wg),
          BP.wg.max=max(BP.wg))

## expand 1000 each window or 5000, subset to uniques, export for pcangsd
names(era.all.pbs.within.shdr.w.and.above.4stdv)
era.all.pbs.within.shdr.w.and.above.4stdv.angsd <- era.all.pbs.within.shdr.w.and.above.4stdv[c("scaff","midPos","BP.wg", "shdr.para.west.id")]; head(era.all.pbs.within.shdr.w.and.above.4stdv)

# windows are 5kb
era.all.pbs.within.shdr.w.and.above.4stdv.angsd$start <-  as.numeric(as.character(era.all.pbs.within.shdr.w.and.above.4stdv.angsd$midPos))-2500; head(era.all.pbs.within.shdr.w.and.above.4stdv.angsd)
era.all.pbs.within.shdr.w.and.above.4stdv.angsd$end <-  as.numeric(as.character(era.all.pbs.within.shdr.w.and.above.4stdv.angsd$midPos))+2500; head(era.all.pbs.within.shdr.w.and.above.4stdv.angsd)

era.all.pbs.within.shdr.w.and.above.4stdv.angsd$start.BP.wg <-  as.numeric(as.character(era.all.pbs.within.shdr.w.and.above.4stdv.angsd$BP.wg))-2500; head(era.all.pbs.within.shdr.w.and.above.4stdv.angsd)
era.all.pbs.within.shdr.w.and.above.4stdv.angsd$end.BP.wg <-  as.numeric(as.character(era.all.pbs.within.shdr.w.and.above.4stdv.angsd$BP.wg))+2500; head(era.all.pbs.within.shdr.w.and.above.4stdv.angsd)

era.all.pbs.within.shdr.w.and.above.4stdv.angsd.summ <- era.all.pbs.within.shdr.w.and.above.4stdv.angsd %>% 
  rowwise() %>% 
  do(data.frame(scaff=.$scaff, midPos= .$midPos, pos = .$start:.$end, shdr.para.west.id=.$shdr.para.west.id, BP.wg=.$start.BP.wg:.$end.BP.wg )); era.all.pbs.within.shdr.w.and.above.4stdv.angsd.summ

era.w.summ <- summarise(group_by(era.all.pbs.within.shdr.w.and.above.4stdv.angsd.summ, shdr.para.west.id, scaff ),
                        bp.min=min(pos ),
                        bp.max=max(pos ),
                        region.size=bp.max-bp.min,
                        BP.wg.min=min(BP.wg),
                        BP.wg.max=max(BP.wg)); era.w.summ

write.csv(era.w.summ, "02.info/sites/shared.exp50kb.era.w.std4.new.era.w.summ.csv", row.names = F)

sites.shared.ids <- list(); regions.shared.ids <- list(); for (i in unique(era.all.pbs.within.shdr.w.and.above.4stdv.angsd.summ$shdr.para.west.id )) {
  sites.shared.ids[[i]] <- subset(era.all.pbs.within.shdr.w.and.above.4stdv.angsd.summ, shdr.para.west.id ==i)
  sites.shared.ids[[i]] <-  sites.shared.ids[[i]][c("scaff", "pos")]
  regions.shared.ids[[i]] <- unique(sites.shared.ids[[i]]$scaff)
  write.table(sites.shared.ids[[i]], paste0("02.info/sites/shared.exp50kb.era.w.std4/sites.shared.region.", i,".txt"), 
              col.names = F, row.names = F, quote =F )
  write.table(regions.shared.ids[[i]], paste0("02.info/sites/shared.exp50kb.era.w.std4/regions.shared.region.", i,".txt"), 
              col.names = F, row.names = F, quote =F )
}


########## MELPOMENE data #####
setwd("/Users/gabrielamontejokovacevich/Dropbox (Cambridge University)/PhD/7_Assocation_studies/9_ANGSD/12.mel.mel.altitude/local/")
ref.scaff.mel <- read.table("local/data/ref/Hmel2.5.fa.fai", row.names = NULL)
all.bams.list.info <- read.csv("02.info/all.bams.list.info.csv"); names(all.bams.list.info)
mel.e.list.info <- subset(all.bams.list.info,to.use.pbs=="yes" &species=="melpomene"&side.short=="e" )
mel.w.list.info <- subset(all.bams.list.info,to.use.pbs=="yes" &species=="melpomene"&side.short=="w" )

mel.all.pbs <- read.csv("local/data/pbs.out/mel.all.pbs.recomb.csv")
mel.all.pbs$CHR <- as.numeric(as.character(substr(mel.all.pbs$scaff,6,7))); unique(mel.all.pbs$CHR)
names(mel.all.pbs)
unique(mel.all.pbs$shdr.para.east.id)

summarise(group_by(mel.all.pbs, shdr.para.east.id),
          bp.min=min(midPos),
          bp.max=max(midPos))



##### mel.e subset so that either clines have >=4stdv #####
mel.all.pbs.within.shdr.and.above.4stdv <- subset(mel.all.pbs, !(is.na(shdr.para.east.id) ) &( zPBS0.1>=4 | zPBS0.3>=4 ))
summarise(group_by(mel.all.pbs.within.shdr.and.above.4stdv, shdr.para.east.id, scaff),
          bp.min=min(midPos),
          bp.max=max(midPos),
          BP.wg.min=min(BP.wg),
          BP.wg.max=max(BP.wg))

## expand 1000 each window or 5000, subset to uniques, export for pcangsd
names(mel.all.pbs.within.shdr.and.above.4stdv)
mel.all.pbs.within.shdr.and.above.4stdv.angsd <- mel.all.pbs.within.shdr.and.above.4stdv[c("scaff","midPos","BP.wg", "shdr.para.east.id")]; head(mel.all.pbs.within.shdr.and.above.4stdv)

# windows are 5kb
mel.all.pbs.within.shdr.and.above.4stdv.angsd$start <-  as.numeric(as.character(mel.all.pbs.within.shdr.and.above.4stdv.angsd$midPos))-2500; head(mel.all.pbs.within.shdr.and.above.4stdv.angsd)
mel.all.pbs.within.shdr.and.above.4stdv.angsd$end <-  as.numeric(as.character(mel.all.pbs.within.shdr.and.above.4stdv.angsd$midPos))+2500; head(mel.all.pbs.within.shdr.and.above.4stdv.angsd)

mel.all.pbs.within.shdr.and.above.4stdv.angsd$start.BP.wg <-  as.numeric(as.character(mel.all.pbs.within.shdr.and.above.4stdv.angsd$BP.wg))-2500; head(mel.all.pbs.within.shdr.and.above.4stdv.angsd)
mel.all.pbs.within.shdr.and.above.4stdv.angsd$end.BP.wg <-  as.numeric(as.character(mel.all.pbs.within.shdr.and.above.4stdv.angsd$BP.wg))+2500; head(mel.all.pbs.within.shdr.and.above.4stdv.angsd)

mel.all.pbs.within.shdr.and.above.4stdv.angsd.summ <- mel.all.pbs.within.shdr.and.above.4stdv.angsd %>% 
  rowwise() %>% 
  do(data.frame(scaff=.$scaff, midPos= .$midPos, pos = .$start:.$end, shdr.para.east.id=.$shdr.para.east.id, BP.wg=.$start.BP.wg:.$end.BP.wg )); mel.all.pbs.within.shdr.and.above.4stdv.angsd.summ

mel.e.summ <- summarise(group_by(mel.all.pbs.within.shdr.and.above.4stdv.angsd.summ, shdr.para.east.id, scaff ),
                        bp.min=min(pos ),
                        bp.max=max(pos ),
                        region.size=bp.max-bp.min,
                        BP.wg.min=min(BP.wg),
                        BP.wg.max=max(BP.wg)); mel.e.summ

write.csv(mel.e.summ, "02.info/sites/shared.exp50kb.mel.e.std4.new.mel.e.summ.csv", row.names = F)

sites.shared.ids <- list(); regions.shared.ids <- list(); for (i in unique(mel.all.pbs.within.shdr.and.above.4stdv.angsd.summ$shdr.para.east.id )) {
  sites.shared.ids[[i]] <- subset(mel.all.pbs.within.shdr.and.above.4stdv.angsd.summ, shdr.para.east.id ==i)
  sites.shared.ids[[i]] <-  sites.shared.ids[[i]][c("scaff", "pos")]
  regions.shared.ids[[i]] <- unique(sites.shared.ids[[i]]$scaff)
  write.table(sites.shared.ids[[i]], paste0("02.info/sites/shared.exp50kb.mel.e.std4/sites.shared.region.", i,".txt"), 
              col.names = F, row.names = F, quote =F )
  write.table(regions.shared.ids[[i]], paste0("02.info/sites/shared.exp50kb.mel.e.std4/regions.shared.region.", i,".txt"), 
              col.names = F, row.names = F, quote =F )
}


##### mel.w subset so that either clines have >=4stdv #####
mel.all.pbs.within.shdr.w.and.above.4stdv <- subset(mel.all.pbs, !(is.na(shdr.para.west.id) ) & (zPBS0.2>=4 | zfst.4>=4 ))
summarise(group_by(mel.all.pbs.within.shdr.w.and.above.4stdv, shdr.para.west.id, scaff),
          bp.min=min(midPos),
          bp.max=max(midPos),
          BP.wg.min=min(BP.wg),
          BP.wg.max=max(BP.wg))

## expand 1000 each window or 5000, subset to uniques, export for pcangsd
names(mel.all.pbs.within.shdr.w.and.above.4stdv)
mel.all.pbs.within.shdr.w.and.above.4stdv.angsd <- mel.all.pbs.within.shdr.w.and.above.4stdv[c("scaff","midPos","BP.wg", "shdr.para.west.id")]; head(mel.all.pbs.within.shdr.w.and.above.4stdv)

# windows are 5kb
mel.all.pbs.within.shdr.w.and.above.4stdv.angsd$start <-  as.numeric(as.character(mel.all.pbs.within.shdr.w.and.above.4stdv.angsd$midPos))-2500; head(mel.all.pbs.within.shdr.w.and.above.4stdv.angsd)
mel.all.pbs.within.shdr.w.and.above.4stdv.angsd$end <-  as.numeric(as.character(mel.all.pbs.within.shdr.w.and.above.4stdv.angsd$midPos))+2500; head(mel.all.pbs.within.shdr.w.and.above.4stdv.angsd)

mel.all.pbs.within.shdr.w.and.above.4stdv.angsd$start.BP.wg <-  as.numeric(as.character(mel.all.pbs.within.shdr.w.and.above.4stdv.angsd$BP.wg))-2500; head(mel.all.pbs.within.shdr.w.and.above.4stdv.angsd)
mel.all.pbs.within.shdr.w.and.above.4stdv.angsd$end.BP.wg <-  as.numeric(as.character(mel.all.pbs.within.shdr.w.and.above.4stdv.angsd$BP.wg))+2500; head(mel.all.pbs.within.shdr.w.and.above.4stdv.angsd)

mel.all.pbs.within.shdr.w.and.above.4stdv.angsd.summ <- mel.all.pbs.within.shdr.w.and.above.4stdv.angsd %>% 
  rowwise() %>% 
  do(data.frame(scaff=.$scaff, midPos= .$midPos, pos = .$start:.$end, shdr.para.west.id=.$shdr.para.west.id, BP.wg=.$start.BP.wg:.$end.BP.wg )); mel.all.pbs.within.shdr.w.and.above.4stdv.angsd.summ

mel.w.summ <- summarise(group_by(mel.all.pbs.within.shdr.w.and.above.4stdv.angsd.summ, shdr.para.west.id, scaff ),
                        bp.min=min(pos ),
                        bp.max=max(pos ),
                        region.size=bp.max-bp.min,
                        BP.wg.min=min(BP.wg),
                        BP.wg.max=max(BP.wg)); mel.w.summ

write.csv(mel.w.summ, "02.info/sites/shared.exp50kb.mel.w.std4.new.mel.w.summ.csv", row.names = F)

sites.shared.ids <- list(); regions.shared.ids <- list(); for (i in unique(mel.all.pbs.within.shdr.w.and.above.4stdv.angsd.summ$shdr.para.west.id )) {
  sites.shared.ids[[i]] <- subset(mel.all.pbs.within.shdr.w.and.above.4stdv.angsd.summ, shdr.para.west.id ==i)
  sites.shared.ids[[i]] <-  sites.shared.ids[[i]][c("scaff", "pos")]
  regions.shared.ids[[i]] <- unique(sites.shared.ids[[i]]$scaff)
  write.table(sites.shared.ids[[i]], paste0("02.info/sites/shared.exp50kb.mel.w.std4/sites.shared.region.", i,".txt"), 
              col.names = F, row.names = F, quote =F )
  write.table(regions.shared.ids[[i]], paste0("02.info/sites/shared.exp50kb.mel.w.std4/regions.shared.region.", i,".txt"), 
              col.names = F, row.names = F, quote =F )
}




############################## genome-wide neutral sites ###############
########## ERATO data #####
setwd("/Users/gabrielamontejokovacevich/Dropbox (Cambridge University)/PhD/7_Assocation_studies/9_ANGSD/12.era.mel.altitude/local/")
ref.scaff.era <- read.table("local/data/ref/Heliconius_erato_demophoon_v1_-_scaffolds.fa.fai", row.names = NULL)
all.bams.list.info <- read.csv("02.info/all.bams.list.info.csv"); names(all.bams.list.info)
era.e.list.info <- subset(all.bams.list.info,to.use.pbs=="yes" &species=="erato"&side.short=="e" )
era.w.list.info <- subset(all.bams.list.info,to.use.pbs=="yes" &species=="erato"&side.short=="w" )

era.all.pbs <- read.csv("local/data/pbs.out/era.all.pbs.recomb.csv")
era.all.pbs$CHR <- as.numeric(as.character(substr(era.all.pbs$scaff,7,8))); unique(era.all.pbs$CHR)
names(era.all.pbs)
unique(era.all.pbs$shdr.para.east.id)

##### era.e subset so that neither cline has PBS or any Fst > 4 stdv #####
era.all.pbs.neutral <- subset(era.all.pbs, is.hdr=="background")

## expand 1000 each window or 5000, subset to uniques, export for pcangsd
era.all.pbs.neutral.angsd <- era.all.pbs.neutral[c("scaff","midPos")]; head(era.all.pbs.neutral)

# windows are 5kb
era.all.pbs.neutral.angsd$start <-  as.numeric(as.character(era.all.pbs.neutral.angsd$midPos))-2500; head(era.all.pbs.neutral.angsd)
era.all.pbs.neutral.angsd$end <-  as.numeric(as.character(era.all.pbs.neutral.angsd$midPos))+2500; head(era.all.pbs.neutral.angsd)

# randomly take 10% of windows to reduce memory
# exclude scaffold with inversion: 
exclude <- c(paste(rep("Herato02",8), str_pad(seq(7, 15), 2, pad = "0") , sep="" ))
era.all.pbs.neutral.angsd <- subset(era.all.pbs.neutral.angsd, !(scaff %in% exclude))
era.all.pbs.neutral.angsd.sub <- era.all.pbs.neutral.angsd[sample(nrow(era.all.pbs.neutral.angsd), 0.1*nrow(era.all.pbs.neutral.angsd)), ]

# or just chr1!!!!!
era.all.pbs.neutral.angsd.sub <- subset(era.all.pbs.neutral.angsd, substr(scaff, 7,8 )=="01")

# get all sites within those windows- careful MEMORY
era.all.pbs.neutral.angsd.sub.summ <- NA
era.all.pbs.neutral.angsd.sub.summ <- era.all.pbs.neutral.angsd.sub %>% 
  rowwise() %>% 
  do(data.frame(scaff=.$scaff, midPos= .$midPos, pos = .$start:.$end)); era.all.pbs.neutral.angsd.sub.summ


# sample every 1000th row to avoid LD, these will be the sites analysed by pcangsd
nrow(era.all.pbs.neutral.angsd.sub.summ)
era.all.pbs.neutral.angsd.sub.summ <- era.all.pbs.neutral.angsd.sub.summ[seq(1, nrow(era.all.pbs.neutral.angsd.sub.summ), 10000), ]; head(era.all.pbs.neutral.angsd.sub.summ)
nrow(era.all.pbs.neutral.angsd.sub.summ)
#era.all.pbs.neutral.angsd.sub.summ <- era.all.pbs.neutral.angsd.sub.summ[seq(1, nrow(era.all.pbs.neutral.angsd.sub.summ), 100000), ]; head(era.all.pbs.neutral.angsd.sub.summ)

# get sites
era.all.pbs.neutral.angsd.sub.summ  <- era.all.pbs.neutral.angsd.sub.summ[c("scaff", "pos")]
# order
era.all.pbs.neutral.angsd.sub.summ <- era.all.pbs.neutral.angsd.sub.summ[order( era.all.pbs.neutral.angsd.sub.summ[,1], era.all.pbs.neutral.angsd.sub.summ[,2] ), ]

regions.shared.ids <- unique(era.all.pbs.neutral.angsd.sub.summ$scaff)

# save same sites/ regions for 3 sets of indivs
nrow(era.all.pbs.neutral.angsd.sub.summ)
substr(era.all.pbs.neutral.angsd.sub.summ$scaff, 7,8)

# # for chr1 only
# write.table(subset(era.all.pbs.neutral.angsd.sub.summ, substr(scaff, 7,8 )=="01" ), "02.info/sites/era.e/sites.era.e.txt", col.names = F, row.names = F, quote =F )
# write.table(unique(subset(era.all.pbs.neutral.angsd.sub.summ, substr(scaff, 7,8 )=="01" )$scaff), "02.info/sites/era.e/regions.era.e.txt", col.names = F, row.names = F, quote =F )


write.table(era.all.pbs.neutral.angsd.sub.summ, "02.info/sites/era.a/sites.era.a.txt", col.names = F, row.names = F, quote =F )
write.table(era.all.pbs.neutral.angsd.sub.summ, "02.info/sites/era.e/sites.era.e.txt", col.names = F, row.names = F, quote =F )
write.table(era.all.pbs.neutral.angsd.sub.summ, "02.info/sites/era.w/sites.era.w.txt", col.names = F, row.names = F, quote =F )

write.table(regions.shared.ids, "02.info/sites/era.a/regions.era.a.txt", col.names = F, row.names = F, quote =F )
write.table(regions.shared.ids, "02.info/sites/era.e/regions.era.e.txt", col.names = F, row.names = F, quote =F )
write.table(regions.shared.ids, "02.info/sites/era.w/regions.era.w.txt", col.names = F, row.names = F, quote =F )

# get region file for new angsd method, not just -region -sites, but -rf directly
write.table(paste(era.all.pbs.neutral.angsd.sub.summ$scaff, era.all.pbs.neutral.angsd.sub.summ$pos, sep=":"), "02.info/sites/era.a/rf.era.a.txt", col.names = F, row.names = F, quote =F )
write.table(paste(era.all.pbs.neutral.angsd.sub.summ$scaff, era.all.pbs.neutral.angsd.sub.summ$pos, sep=":"), "02.info/sites/era.e/rf.era.e.txt", col.names = F, row.names = F, quote =F )
write.table(paste(era.all.pbs.neutral.angsd.sub.summ$scaff, era.all.pbs.neutral.angsd.sub.summ$pos, sep=":"), "02.info/sites/era.w/rf.era.w.txt", col.names = F, row.names = F, quote =F )


########## melpomene data #####
setwd("/Users/gabrielamontejokovacevich/Dropbox (Cambridge University)/PhD/7_Assocation_studies/9_ANGSD/12.era.mel.altitude/local/")
ref.scaff.mel <- read.table("local/data/ref/Heliconius_melto_demophoon_v1_-_scaffolds.fa.fai", row.names = NULL)
all.bams.list.info <- read.csv("02.info/all.bams.list.info.csv"); names(all.bams.list.info)
mel.e.list.info <- subset(all.bams.list.info,to.use.pbs=="yes" &species=="melpomene"&side.short=="e" )
mel.w.list.info <- subset(all.bams.list.info,to.use.pbs=="yes" &species=="melpomene"&side.short=="w" )

mel.all.pbs <- read.csv("local/data/pbs.out/mel.all.pbs.recomb.csv")
mel.all.pbs$CHR <- as.numeric(as.character(substr(mel.all.pbs$scaff,6,7))); unique(mel.all.pbs$CHR)
names(mel.all.pbs)
unique(mel.all.pbs$shdr.para.east.id)

##### mel.e subset so that neither cline has PBS or any Fst > 4 stdv #####
mel.all.pbs.neutral <- subset(mel.all.pbs, is.hdr=="background")

## expand 1000 each window or 5000, subset to uniques, export for pcangsd
mel.all.pbs.neutral.angsd <- mel.all.pbs.neutral[c("scaff","midPos")]; head(mel.all.pbs.neutral)

# windows are 5kb
mel.all.pbs.neutral.angsd$start <-  as.numeric(as.character(mel.all.pbs.neutral.angsd$midPos))-2500; head(mel.all.pbs.neutral.angsd)
mel.all.pbs.neutral.angsd$end <-  as.numeric(as.character(mel.all.pbs.neutral.angsd$midPos))+2500; head(mel.all.pbs.neutral.angsd)

# randomly take 10% of windows to reduce memory
mel.all.pbs.neutral.angsd.sub <- mel.all.pbs.neutral.angsd[sample(nrow(mel.all.pbs.neutral.angsd), 0.1*nrow(mel.all.pbs.neutral.angsd)), ]

# # or just chr1!!!!!
# unique(substr(mel.all.pbs.neutral.angsd$scaff, 6,7 ))
# mel.all.pbs.neutral.angsd.sub <- subset(mel.all.pbs.neutral.angsd, substr(scaff, 6,7 )=="01")


# get all sites within those windows- careful MEMORY
mel.all.pbs.neutral.angsd.sub.summ <- NA
mel.all.pbs.neutral.angsd.sub.summ <- mel.all.pbs.neutral.angsd.sub %>% 
  rowwise() %>% 
  do(data.frame(scaff=.$scaff, midPos= .$midPos, pos = .$start:.$end)); mel.all.pbs.neutral.angsd.sub.summ


# sample every 1000th row to avoid LD, these will be the sites analysed by pcangsd
nrow(mel.all.pbs.neutral.angsd.sub.summ)
mel.all.pbs.neutral.angsd.sub.summ <- mel.all.pbs.neutral.angsd.sub.summ[seq(1, nrow(mel.all.pbs.neutral.angsd.sub.summ), 10000), ]; head(mel.all.pbs.neutral.angsd.sub.summ)
nrow(mel.all.pbs.neutral.angsd.sub.summ)

# get sites
mel.all.pbs.neutral.angsd.sub.summ  <- mel.all.pbs.neutral.angsd.sub.summ[c("scaff", "pos")]
# order
mel.all.pbs.neutral.angsd.sub.summ <- mel.all.pbs.neutral.angsd.sub.summ[order( mel.all.pbs.neutral.angsd.sub.summ[,1], mel.all.pbs.neutral.angsd.sub.summ[,2] ), ]
regions.shared.ids <- unique(mel.all.pbs.neutral.angsd.sub.summ$scaff)

# save same sites/ regions for 3 sets of indivs
nrow(mel.all.pbs.neutral.angsd.sub.summ)
write.table(mel.all.pbs.neutral.angsd.sub.summ, "02.info/sites/mel.a/sites.mel.a.txt", col.names = F, row.names = F, quote =F )
write.table(mel.all.pbs.neutral.angsd.sub.summ, "02.info/sites/mel.e/sites.mel.e.txt", col.names = F, row.names = F, quote =F )
write.table(mel.all.pbs.neutral.angsd.sub.summ, "02.info/sites/mel.w/sites.mel.w.txt", col.names = F, row.names = F, quote =F )

write.table(regions.shared.ids, "02.info/sites/mel.a/regions.mel.a.txt", col.names = F, row.names = F, quote =F )
write.table(regions.shared.ids, "02.info/sites/mel.e/regions.mel.e.txt", col.names = F, row.names = F, quote =F )
write.table(regions.shared.ids, "02.info/sites/mel.w/regions.mel.w.txt", col.names = F, row.names = F, quote =F )

# get region file for new angsd method, not just -region -sites, but -rf directly
write.table(paste(mel.all.pbs.neutral.angsd.sub.summ$scaff, mel.all.pbs.neutral.angsd.sub.summ$pos, sep=":"), "02.info/sites/mel.a/rf.mel.a.txt", col.names = F, row.names = F, quote =F )
write.table(paste(mel.all.pbs.neutral.angsd.sub.summ$scaff, mel.all.pbs.neutral.angsd.sub.summ$pos, sep=":"), "02.info/sites/mel.e/rf.mel.e.txt", col.names = F, row.names = F, quote =F )
write.table(paste(mel.all.pbs.neutral.angsd.sub.summ$scaff, mel.all.pbs.neutral.angsd.sub.summ$pos, sep=":"), "02.info/sites/mel.w/rf.mel.w.txt", col.names = F, row.names = F, quote =F )




############ 


