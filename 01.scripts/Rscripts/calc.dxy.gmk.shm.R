####
library(dplyr)
library(optparse)
library(windowscanr) 
library(data.table)
option_list = list(
  make_option(c("-a","--popA"), type="character",default=NULL,help="path to unzipped mafs file for pop 1",metavar="character"),
  make_option(c("-b","--popB"), type="character",default=NULL,help="path to unzipped mafs file for pop 2",metavar="character"),
  make_option(c("-c","--popAB"), type="character",default=NULL,help="path to unzipped mafs file for pop 2",metavar="character"),
  make_option(c("-d","--nIndA"), type="numeric",default=NULL,help="number of individuals in pop",metavar="numeric"),
  make_option(c("-e","--nIndB"), type="numeric",default=NULL,help="number of individuals in pop",metavar="numeric"),
  make_option(c("-f","--nIndAB"), type="numeric",default=NULL,help="number of individuals in joint pop",metavar="numeric"),
  make_option(c("-w","--winSize"), type="numeric",default=NULL,help="window size (in kb) to split dxy per site [optional]",metavar="numeric"),
  make_option(c("-s","--winStep"), type="numeric",default=NULL,help="window step (in kb) to split dxy per site [optional]",metavar="numeric"),
  make_option(c("-t","--cpu"), type="numeric", default=NULL,help="how many cpu for sliding windows scan [optional]",metavar="numeric")
)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

# input sample sizes for testing
# nAB=35; nA=14; nB=21
nA=as.numeric(as.character(opt$nIndA)); nB=as.numeric(as.character(opt$nIndB)); nAB=as.numeric(as.character(opt$nIndAB))
str(nA); str(nB); str(nAB)
str(opt)

# read in thetas
thetaAB <- fread(paste0(opt$popAB, ".thetas.txt", sep="") )
thetaA <- fread(paste0(opt$popA, ".thetas.txt", sep="") )
thetaB <- fread(paste0(opt$popB, ".thetas.txt", sep=""))

# name columns
names(thetaA) <- c("Chromo",	"Pos",	"Watterson",	"Pairwise",	"thetaSingleton",	"thetaH",	"thetaL"); thetaA <- thetaA[,c(1,2,4)]
names(thetaB) <- c("Chromo",	"Pos",	"Watterson",	"Pairwise",	"thetaSingleton",	"thetaH",	"thetaL"); thetaB <- thetaB[,c(1,2,4)]
names(thetaAB) <- c("Chromo",	"Pos",	"Watterson",	"Pairwise",	"thetaSingleton",	"thetaH",	"thetaL"); thetaAB <- thetaAB[,c(1,2,4)]

# join datasets
theta <- left_join(thetaA, thetaB, by = c("Chromo", "Pos")) %>%
  left_join(thetaAB, by = c("Chromo", "Pos")) %>%  rename(thetaA.p = Pairwise.x, thetaB.p = Pairwise.y, thetaAB.p = Pairwise )
str(theta)

## dxy calculation ##
# normalise diversity by pop size
print("going to calculate dxy")
theta$dxy <-( exp(theta$thetaAB.p)*(factorial(nAB)/(factorial(2)*factorial(nAB-2))) -
              exp(theta$thetaA.p)*(factorial(nA)/(factorial(2)*factorial(nA-2))) -
              exp(theta$thetaB.p)*(factorial(nB)/(factorial(2)*factorial(nB-2)))  )/ (nA * nB)

theta.sub <- subset(theta, dxy!="")

print("going to start making sliding windows")

#### windows #####
# 5kb 1kb windows, position fixed
#win_size = opt$winSize*1000, , cores=opt$cpu
dxy.windows <- winScan(x = theta.sub,
                       position = "Pos",
                       values = "dxy",
                       groups = "Chromo",
                       win_size = opt$winSize*1000,
                       win_step = opt$winStep*1000,
                       funs = c("mean"), cores=opt$cpu)

write.table(dxy.windows, file=paste0(opt$popAB, ".dxy.", opt$winSize,"kb.", opt$winStep, "kb.txt"),quote=FALSE, row.names=FALSE, sep='\t')

print("created sliding window dxy file")

