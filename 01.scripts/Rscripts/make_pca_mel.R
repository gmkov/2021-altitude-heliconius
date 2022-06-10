#this R script performs a PCa on the covariance matrix and export it labelling rows with bamfile name and columns as PC
#it also makes a first visualisation in a pdf file with PC1-2 and PC3-4

# install.packages("Rcpp")
# install.packages("RcppCNPy")
library(RcppCNPy)

argv <- commandArgs(T)
MIN_MAF_MEL <- argv[1]
PERCENT_IND <- argv[2]
BAM <- argv[3]
MAX_DEPTH <-argv[4]
DOWNSAMPLING_PERC <- argv [5]

#perform a pca on covariance matrix, import as npy
cov_mat<-npyLoad(file=paste0("../../04.pca/mel/all_maf",MIN_MAF_MEL,"_maxdepth",MAX_DEPTH,"_downsample",DOWNSAMPLING_PERC,".cov.npy"))

#perform a pca on covariance matrix
#cov_mat<-as.matrix(read.table(paste0("../../04.pca/mel/all_maf",MIN_MAF_MEL,"_maxdepth",MAX_DEPTH,"_downsample",DOWNSAMPLING_PERC,".cov"), header=F))
pca<-prcomp(cov_mat)
print(summary(pca))

#add column names
nPC<-dim(pca$x)[2]
col_PC<-vector(length=nPC)
for (i in 1 : nPC) {col_PC[i]<-paste0("PC",i)}
colnames(pca$x)<-c(col_PC)

#add rownames
bam_names<-read.table(BAM,header=F)
rownames(pca$x)<-bam_names$V1

#plot pca
pdf(file=paste0("../../04.pca/mel/all_maf",MIN_MAF_MEL,"_maxdepth",MAX_DEPTH,"_downsample",DOWNSAMPLING_PERC,".pca.pdf"))
par(mfrow=c(1,1))
plot(pca$x[,1], pca$x[,2], pch=20, ylab="PC2", xlab="PC1")
#plot(pca$x[,3], pca$x[,4], pch=20, ylab="PC4", xlab="PC3")
dev.off()

write.table(pca$x, paste0("../../04.pca/mel/all_maf",MIN_MAF_MEL,"_maxdepth",MAX_DEPTH,"_downsample",DOWNSAMPLING_PERC,".pca"), quote=F)
