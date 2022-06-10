#this R script uses the mafs file provided by the analyses on all individuals, with filter for quality and maf 0.05
#it simply extract the first columns with chr, position of each SNP and Major/minor alleles as determined in step 03
#output is a "sites" files that will allwo restraining subsequent analyses 
#(for instance maf by pop, FSt by pop pairs, etc to a limited number of loci

argv <- commandArgs(T)
POP_LIST <- argv[1]
MIN_MAF <- argv[2]
PERCENT_IND <- argv[3]
MAX_DEPTH <- argv [4]
DOWNSAMPLING_PERC <- argv [5]

#maf<-read.table( paste0(POP_LIST,MIN_MAF,"_pctind",PERCENT_IND,"_maxdepth",MAX_DEPTH,"_downsample",DOWNSAMPLING_PERC,".unlinked.header.mafs"), header=T)
maf<-read.table( paste0(POP_LIST,".filtered.mafs"), header=T)

head (maf)

### minMAC (optional, if want to subset by minMAC)
# obtain allele counts
maf$mac <- maf$knownEM * maf$nInd


### no minMAC
# grab site names to index
sites<-maf[,1:4]
# remove mt sites
sites <- subset(sites, chromo!="Herato_mt")
sites_order <- sites[order(sites$chromo),] 
regions_order<-unique(sites_order$chromo)
write.table(sites_order, paste0("../../../02.info/sites/sites.bedassle.maf.", POP_LIST,".txt"), row.names=F, col.names=F, sep="\t", quote=F)
write.table(regions_order, paste0("../../../02.info/sites/regions.bedassle.maf.", POP_LIST,".txt"), row.names=F, col.names=F, sep="\t", quote=F)


