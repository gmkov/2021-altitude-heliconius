#this R script uses the mafs file provided by the analyses on all individuals, with filter for quality and maf 0.05
#it simply extract the first columns with chr, position of each SNP and Major/minor alleles as determined in step 03
#output is a "sites" files that will allwo restraining subsequent analyses 
#(for instance maf by pop, FSt by pop pairs, etc to a limited number of loci

argv <- commandArgs(T)
MIN_MAF <- argv[1]
PERCENT_IND <- argv[2]
MAX_DEPTH_FACTOR <- argv [3]
DOWNSAMPLING_PERC <- argv [4]

maf<-read.table( paste0("../../03.saf.maf.gl.all/era/all_maf",MIN_MAF,"_pctind",PERCENT_IND,"_maxdepth",MAX_DEPTH_FACTOR,"_downsample",DOWNSAMPLING_PERC,".mafs"), header=T)
head (maf)

### minMAC
# obtain allele counts
maf$mac <- maf$knownEM * maf$nInd

# apply minMAC =>2
maf.minmac2 <- subset(maf, mac>=2)

# grab site names to index
sites<-maf.minmac2[,1:4]

sites_order <- sites[order(sites$chromo),] 
regions_order<-levels(sites_order$chromo)
write.table(sites_order, paste0("../../02.info/sites_era_maf", MIN_MAF, "_pctind",PERCENT_IND,"_maxdepth",MAX_DEPTH_FACTOR,"_downsample",DOWNSAMPLING_PERC,"_minMAC2"), row.names=F, col.names=F, sep="\t", quote=F)
write.table(regions_order, paste0("../../02.info/regions_era_maf", MIN_MAF, "_pctind",PERCENT_IND,"_maxdepth",MAX_DEPTH_FACTOR,"_downsample",DOWNSAMPLING_PERC,"_minMAC2"), row.names=F, col.names=F, sep="\t", quote=F)

### no minMAC
# grab site names to index
sites<-maf[,1:4]
sites_order <- sites[order(sites$chromo),] 
regions_order<-levels(sites_order$chromo)
write.table(sites_order, paste0("../../02.info/sites_era_maf", MIN_MAF, "_pctind",PERCENT_IND,"_maxdepth",MAX_DEPTH_FACTOR,"_downsample",DOWNSAMPLING_PERC), row.names=F, col.names=F, sep="\t", quote=F)
write.table(regions_order, paste0("../../02.info/regions_era_maf", MIN_MAF, "_pctind",PERCENT_IND,"_maxdepth",MAX_DEPTH_FACTOR,"_downsample",DOWNSAMPLING_PERC), row.names=F, col.names=F, sep="\t", quote=F)
