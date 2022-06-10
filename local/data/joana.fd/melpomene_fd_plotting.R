
# fd melpomene relatives
setwd("D:/Dropbox/Heliconius/HighAltitude/fd/")

require(data.table)

#### Read in outlier regions and chromosome information ####

# Get melpomene outliers shared by either or both East and West:
outliersRegions<-read.csv("D:/Dropbox/Heliconius/HighAltitude/GabbyRegionsOfInterest/mel.summ.regions.shared.csv",header=T)
outliersEast<-outliersRegions[!is.na(outliersRegions$east.sharing.id),]
outliersWest<-outliersRegions[!is.na(outliersRegions$west.sharing.id),]
outliersAcross<-outliersRegions[!is.na(outliersRegions$west.sharing.id)&!is.na(outliersRegions$east.sharing.id),]

# Read in chromosome information
chrom<-read.table("Hmel2.5.fa.fai")
chrom<-chrom[,1:2]
names(chrom)<-c("scaffold","length")
chrom$add<-c(0,cumsum(chrom$length)[-length(chrom$length)])
chrom$chr<-substring(chrom$scaffold,first = 6,last=7)



#### Functions ####

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
getFd<-function(combination){
  prefix="melpomeneRelatives.withSRA.chr1-21.max0.5N.minDP3.biallSNPs.mac2."
  
  # Read in the dataset
  dataset<-read.csv(paste0(prefix,combination,"_50kb"))
  
  # Correct negative values
  dataset<-corrNegfd(dataset)
  
  # reformat to datatable
  dataset<-data.table(dataset)
  
  # Add chromosome name
  dataset$chr<-substring(dataset$scaffold,first = 7,last=8)
  
  return(dataset)
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

# Plot fd along the genome highlighting outlier regions
plotOutl<-function(data,title="",maxy=0.5,horiz=F){
  
  # get outlier overlap datasets
  dataEastOv<-eastOverlap(data)
  dataWestOv<-westOverlap(data)
  dataAcrossOv<-acrossOverlap(data)
  
  # If no title given, use name of the fd dataset
  if(title=="") title=deparse(substitute(data))
  
  # Plot the fd values along the genome
  plot(x=data$start+chrom[match(data$scaffold,chrom$scaffold),"add"],y=data$fd,
       col="grey",cex=0.3,xaxt="n",ylab="fd",ylim=c(0,maxy),xaxs='i',yaxs='i')
  points(x=dataEastOv$start+chrom[match(dataEastOv$scaffold,chrom$scaffold),"add"],y=dataEastOv$fd,
         col="darkorange",cex=ifelse(dataEastOv$fd>quantile(data$fd,probs = 0.95),1,0.7),pch=19)
  points(x=dataWestOv$start+chrom[match(dataWestOv$scaffold,chrom$scaffold),"add"],y=dataWestOv$fd,
         col="cornflowerblue",cex=ifelse(dataWestOv$fd>quantile(data$fd,probs = 0.95),1,0.7),pch=19)
  points(x=dataAcrossOv$start+chrom[match(dataAcrossOv$scaffold,chrom$scaffold),"add"],y=dataAcrossOv$fd,
         col="black",cex=ifelse(dataAcrossOv$fd>quantile(data$fd,probs = 0.95),1,0.7),pch=19)
  abline(v=chrom$add[!duplicated(chrom$chr)],col="darkgrey")
  title(main=paste0(title," (mean all windows fd: ",round(mean(data$fd),digits=4),
                    ", WestOutl: ",round(mean(dataWestOv$fd),digits=4),
                    ", EastOutl: ",round(mean(dataEastOv$fd),digits=4),
                    ")"))
  legend("topright",legend=c(paste0("highfdM WestOutl prop:",round(length(which(dataWestOv$fdM>quantile(data$fdM,probs = 0.95)))/length(dataWestOv$fdM),digits=2)),
                             paste0("lowfdM WestOutl prop:",round(length(which(dataWestOv$fdM<quantile(data$fdM,probs = 0.05)))/length(dataWestOv$fdM),digits=2)),
                             paste0("highfdM EastOutl prop:",round(length(which(dataEastOv$fdM>quantile(data$fdM,probs = 0.95)))/length(dataEastOv$fdM),digits=2)),
                             paste0("lowfdM EastOutl prop:",round(length(which(dataEastOv$fdM<quantile(data$fdM,probs = 0.05)))/length(dataEastOv$fdM),digits=2))),bty="n",horiz=horiz)

}



#### Read in fd datasets ####

prefix="melpomeneRelatives.withSRA.chr1-21.max0.5N.minDP3.biallSNPs.mac2."

# melpomene across sides of the Andes
fd_melCvlowW_melChighW_melChighE_hecale<-getFd("fd_melCvlowW_melChighW_melChighE_hecale")
fd_melElowW_melEhighW_melEhighE_hecale<-getFd("fd_melElowW_melEhighW_melEhighE_hecale")
fd_melCvlowE_melChighE_melChighW_hecale<-getFd("fd_melCvlowE_melChighE_melChighW_hecale")
fd_melEvlowE_melEhighE_melEhighW_hecale<-getFd("fd_melEvlowE_melEhighE_melEhighW_hecale")
fd_melCvlowW_melChighW_melCvlowE_hecale<-getFd("fd_melCvlowW_melChighW_melCvlowE_hecale")
fd_melElowW_melEhighW_melEvlowE_hecale<-getFd("fd_melElowW_melEhighW_melEvlowE_hecale")
fd_melCvlowE_melChighE_melCvlowW_hecale<-getFd("fd_melCvlowE_melChighE_melCvlowW_hecale")
fd_melEvlowE_melEhighE_melElowW_hecale<-getFd("fd_melEvlowE_melEhighE_melElowW_hecale")


# melpomene on the same Andes side
fd_melCvlowW_melChighW_melEhighW_hecale<-getFd("fd_melCvlowW_melChighW_melEhighW_hecale")
fd_melElowW_melEhighW_melChighW_hecale<-getFd("fd_melElowW_melEhighW_melChighW_hecale")
fd_melCvlowW_melChighW_melElowW_hecale<-getFd("fd_melCvlowW_melChighW_melElowW_hecale")
fd_melElowW_melEhighW_melCvlowW_hecale<-getFd("fd_melElowW_melEhighW_melCvlowW_hecale")
fd_melCvlowE_melChighE_melEhighE_hecale<-getFd("fd_melCvlowE_melChighE_melEhighE_hecale")
fd_melEvlowE_melEhighE_melChighE_hecale<-getFd("fd_melEvlowE_melEhighE_melChighE_hecale")
fd_melCvlowE_melChighE_melEvlowE_hecale<-getFd("fd_melCvlowE_melChighE_melEvlowE_hecale")
fd_melEvlowE_melEhighE_melCvlowE_hecale<-getFd("fd_melEvlowE_melEhighE_melCvlowE_hecale")

# timareta introgression
fd_melCvlowW_melChighW_timEhighE_hecale<-getFd("fd_melCvlowW_melChighW_timEhighE_hecale")
fd_melCvlowE_melChighE_timEhighE_hecale<-getFd("fd_melCvlowE_melChighE_timEhighE_hecale")
fd_melElowW_melEhighW_timEhighE_hecale<-getFd("fd_melElowW_melEhighW_timEhighE_hecale")
fd_melEvlowE_melEhighE_timEhighE_hecale<-getFd("fd_melEvlowE_melEhighE_timEhighE_hecale")

# cydno introgression
fd_melCvlowE_melChighE_cydEhighW_hecale<-getFd("fd_melCvlowE_melChighE_cydEhighW_hecale")
fd_melCvlowW_melChighW_cydEhighW_hecale<-getFd("fd_melCvlowW_melChighW_cydEhighW_hecale")
fd_melElowW_melEhighW_cydEhighW_hecale<-getFd("fd_melElowW_melEhighW_cydEhighW_hecale")
fd_melEvlowE_melEhighE_cydEhighW_hecale<-getFd("fd_melEvlowE_melEhighE_cydEhighW_hecale")


#### Plot fd ####

# Plot erato allele sharing across sides of the Andes

pdf(file="melpomene_acrossAndes.pdf",height=7.5,width=15)
  par(mfrow=c(8,1),mar=c(0,4,1.5,1),oma=c(2,1,1,1))
  plotOutl(fd_melCvlowW_melChighW_melChighE_hecale,title="fd_melCvlowW_melChighW_melChighE_hecale",maxy = 0.5,horiz=T)
  legend("topleft",pch=19,col=c("grey","darkorange","cornflowerblue","black"),
         legend=c("normal","East outlier","West outlier","across outlier"),
         bty="n",ncol=2)
  plotOutl(fd_melCvlowW_melChighW_melCvlowE_hecale,title="fd_melCvlowW_melChighW_melCvlowE_hecale",maxy = 0.5,horiz=T)
  plotOutl(fd_melElowW_melEhighW_melEhighE_hecale,title="fd_melElowW_melEhighW_melEhighE_hecale",maxy = 0.5,horiz=T)
  plotOutl(fd_melElowW_melEhighW_melEvlowE_hecale,title="fd_melElowW_melEhighW_melEvlowE_hecale",maxy = 0.5,horiz=T)
  plotOutl(fd_melCvlowE_melChighE_melChighW_hecale,title="fd_melCvlowE_melChighE_melChighW_hecale",maxy = 0.5,horiz=T)
  plotOutl(fd_melCvlowE_melChighE_melCvlowW_hecale,title="fd_melCvlowE_melChighE_melCvlowW_hecale",maxy = 0.5,horiz=T)
  plotOutl(fd_melEvlowE_melEhighE_melEhighW_hecale,title="fd_melEvlowE_melEhighE_melEhighW_hecale",maxy = 0.5,horiz=T)
  plotOutl(fd_melEvlowE_melEhighE_melElowW_hecale,title="fd_melEvlowE_melEhighE_melElowW_hecale",maxy = 0.5,horiz=T)
  axis(1,chrom$add[!duplicated(chrom$chr)][-22]+(diff(chrom$add[!duplicated(chrom$chr)])/2),
       labels = paste0("chr",chrom$chr[!duplicated(chrom$chr)][-22]),tick = F)
dev.off()


# Plot melpomene allele sharing on the same side of the Andes
pdf(file="melpomene_sameSideAndes.pdf",height=7.5,width=15)
  par(mfrow=c(8,1),mar=c(0,4,2,1),oma=c(2,1,1,1))
  plotOutl(fd_melCvlowW_melChighW_melEhighW_hecale,title="fd_melCvlowW_melChighW_melEhighW_hecale",maxy = 1,horiz=T)
  legend("topleft",pch=19,col=c("grey","darkorange","cornflowerblue","black"),legend=c("normal","East outlier","West outlier","across outlier"),horiz = F,ncol=2,bty="n")
  plotOutl(fd_melCvlowW_melChighW_melElowW_hecale,title="fd_melCvlowW_melChighW_melElowW_hecale",maxy = 1,horiz=T)
  plotOutl(fd_melElowW_melEhighW_melChighW_hecale,title="fd_melElowW_melEhighW_melChighW_hecale",maxy = 1,horiz=T)
  plotOutl(fd_melElowW_melEhighW_melCvlowW_hecale,title="fd_melElowW_melEhighW_melCvlowW_hecale",maxy = 1,horiz=T)
  plotOutl(fd_melCvlowE_melChighE_melEhighE_hecale,title="fd_melCvlowE_melChighE_melEhighE_hecale",maxy = 1,horiz=T)
  plotOutl(fd_melCvlowE_melChighE_melEvlowE_hecale,title="fd_melCvlowE_melChighE_melEvlowE_hecale",maxy = 1,horiz=T)
  plotOutl(fd_melEvlowE_melEhighE_melChighE_hecale,title="fd_melEvlowE_melEhighE_melChighE_hecale",maxy = 1,horiz=T)
  plotOutl(fd_melEvlowE_melEhighE_melCvlowE_hecale,title="fd_melEvlowE_melEhighE_melCvlowE_hecale",maxy = 1,horiz=T)
  axis(1,chrom$add[!duplicated(chrom$chr)][-22]+(diff(chrom$add[!duplicated(chrom$chr)])/2),labels = paste0("chr",chrom$chr[!duplicated(chrom$chr)][-22]),tick = F)
dev.off()

# timareta introgression
pdf(file="melpomene_timareta.pdf",height=7.5,width=15)

par(mfrow=c(4,1),mar=c(0,4,2,1),oma=c(2,1,1,1))

plotOutl(fd_melCvlowW_melChighW_timEhighE_hecale,title="fd_melCvlowW_melChighW_timEhighE_hecale",maxy = 0.7)
legend("topleft",pch=19,col=c("grey","darkorange","cornflowerblue","black"),legend=c("normal","East outlier","West outlier","across outlier"),bty="n")
plotOutl(fd_melCvlowE_melChighE_timEhighE_hecale,title="fd_melCvlowE_melChighE_timEhighE_hecale",maxy = 0.7)
plotOutl(fd_melElowW_melEhighW_timEhighE_hecale,title="fd_melElowW_melEhighW_timEhighE_hecale",maxy = 0.7)
plotOutl(fd_melEvlowE_melEhighE_timEhighE_hecale,title="fd_melEvlowE_melEhighE_timEhighE_hecale",maxy = 0.7)
axis(1,chrom$add[!duplicated(chrom$chr)][-22]+(diff(chrom$add[!duplicated(chrom$chr)])/2),labels = paste0("chr",chrom$chr[!duplicated(chrom$chr)][-22]),tick = F)

dev.off()


# cydno introgression
pdf(file="melpomene_cydno.pdf",height=7.5,width=15)

par(mfrow=c(4,1),mar=c(0,4,2,1),oma=c(2,1,1,1))

plotOutl(fd_melCvlowE_melChighE_cydEhighW_hecale,maxy = 0.7)
legend("topleft",pch=19,col=c("grey","darkorange","cornflowerblue","black"),legend=c("normal","East outlier","West outlier","across outlier"),bty="n")
plotOutl(fd_melCvlowW_melChighW_cydEhighW_hecale,maxy = 0.7)
plotOutl(fd_melElowW_melEhighW_cydEhighW_hecale,maxy = 0.7)
plotOutl(fd_melEvlowE_melEhighE_cydEhighW_hecale,maxy = 0.7)
axis(1,chrom$add[!duplicated(chrom$chr)][-22]+(diff(chrom$add[!duplicated(chrom$chr)])/2),labels = paste0("chr",chrom$chr[!duplicated(chrom$chr)][-22]),tick = F)

dev.off()





# Difference in fd between P3=high - P3=low
diffCE<-fd_melCvlowE_melChighE_melEhighE_hecale
diffCE$fd<-diffCE$fd-fd_melCvlowE_melChighE_melEvlowE_hecale$fd
diffEE<-fd_melEvlowE_melEhighE_melChighE_hecale
diffEE$fd<-diffEE$fd-fd_melEvlowE_melEhighE_melCvlowE_hecale$fd
diffCW<-fd_melCvlowW_melChighW_melEhighW_hecale
diffCW$fd<-diffCW$fd-fd_melCvlowW_melChighW_melElowW_hecale$fd
diffEW<-fd_melElowW_melEhighW_melChighW_hecale[paste0(fd_melElowW_melEhighW_melChighW_hecale$scaffold,fd_melElowW_melEhighW_melChighW_hecale$start)%in%paste0(fd_melElowW_melEhighW_melCvlowW_hecale$scaffold,fd_melElowW_melEhighW_melCvlowW_hecale$start),]
diffEW$fd<-diffEW$fd-fd_melElowW_melEhighW_melCvlowW_hecale$fd


pdf(file="melpomene_P3high-low.pdf",height=7.5,width=15)

par(mfrow=c(4,1),mar=c(0,4,2,1),oma=c(2,1,1,1))

plotOutl(diffCE,title="fd_melCvlowE_melChighE_melEhighE-ELE_hecale",maxy = 0.6)
legend("topleft",pch=19,col=c("grey","darkorange","cornflowerblue","black"),legend=c("normal","East outlier","West outlier","across outlier"),horiz = F,ncol=2,bty="n")
plotOutl(diffEE,title="melEvlowE_melEhighE_melChighE-CLE_hecale",maxy = 0.6)
plotOutl(diffCW,title="melCvlowW_melChighW_melEhighW-ELW_hecale",maxy = 0.6)
plotOutl(diffEW,title="melElowW_melEhighW_melChighW-CLW_hecale",maxy = 0.6)
axis(1,chrom$add[!duplicated(chrom$chr)][-22]+(diff(chrom$add[!duplicated(chrom$chr)])/2),labels = paste0("chr",chrom$chr[!duplicated(chrom$chr)][-22]),tick = F)

dev.off()



