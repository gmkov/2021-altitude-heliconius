library("BEDASSLE") 

#setwd("../../17_bedassle/output/")

#### find files with maf.ac ####
listpaths <- dir(pattern = ".*mafs.ac$"); listpaths
ac <- list()
for (i in 1:length(listpaths)) {
  ac[[i]] <- read.table(listpaths[i])
  colnames(ac[[i]]) <- c("scaff", "pos","major", "minor", "knownEM", "nInd", "ac")
  ac[[i]]$chr.pos <- paste(ac[[i]]$scaff,ac[[i]]$pos, sep = "" )
  ac[[i]]$pop <- substr(listpaths[i],1,12)}

AREA  <- paste(substr(listpaths[1],1,3),substr(listpaths[1],7,8), sep=""); AREA

#### INTERSECT- so that only sites with data for all pops are kept #####
# grabs eigth column (chr.pos) of any df within list
chr.pos.all <- Reduce(intersect, list(c(lapply(ac, "[[", 8)))[[1]] ); 
print(paste("sites that intersect across all pops = ", length(chr.pos.all)))

#### subset all to find intersects ####
ac <- lapply(ac, subset, chr.pos %in% chr.pos.all)
# check that all pos the same
pos.names <-  lapply(ac, "[[", 8); 
if(all(pos.names[[1]]==pos.names[[2]])) print(" sites pop1 and pop2 equal") else print(" sites pop1 and pop2 NOT equal");
if(all(pos.names[[3]]==pos.names[[4]])) print(" sites pop3 and pop4 equal") else print(" sites pop3 and pop4 NOT equal")

#### subset for unlinked sites across ac, at least 2kb apart #####
library(dplyr)
pos.apart.subset <-  function(df, pos.diff) {
  # create new dfs to store output
  nrows <- 1
  previous_match <- 1
  # iterate through each row of df
  for (i in 2:nrow(df)) {
    if(df$pos[i] - df$pos[previous_match] > pos.diff | df$pos[i] - df$pos[previous_match] < 0) {
      # add rows that meet conditions to list
      nrows <- c(nrows, i)
      previous_match <- i
    }
  } 
  df1 <- df[nrows,]
  return(df1)
}

ac.sub.2k <-list()
ac.sub.2k  <- lapply(ac, pos.apart.subset, 2000)
print("check that the number of sites matches in all pops")
for (i in 1:length(ac.sub.2k)){ print(dim(ac.sub.2k[[i]])) }

#### matrix #### 
ac.mat <- do.call(rbind, lapply(ac.sub.2k, "[[", 7)); dim(ac.mat)
rownames(ac.mat) <- unique(unlist(lapply(ac.sub.2k, "[[", 9))); ac.mat[,1:10]
colnames(ac.mat) <- unique(unlist(lapply(ac.sub.2k, "[[", 8))); ac.mat[,1:10]
n.mat <- do.call(rbind, lapply(ac.sub.2k, "[[", 6)); 
rownames(n.mat) <- unique(unlist(lapply(ac.sub.2k, "[[", 9))); 
colnames(n.mat) <- unique(unlist(lapply(ac.sub.2k, "[[", 8))); n.mat[,1:3]

# double allele counts and n
ac.mat <- ac.mat*2
n.mat <- n.mat*2

# ## make matrix if need to fill with 0 #####
# print("obtain ac mat and n mat")
# # prep loci names (columns)
# chr.pos.all <- unique(unlist(c(lapply(ac, "[[", 8)))); chr.pos.all[1:10]
# # prep matrix to populate
# ac.mat <- matrix(NA, ncol = length(chr.pos.all), nrow = length( unique(unlist(lapply(ac, "[[", 9)))),
#                   dimnames = list(c(unique(unlist(lapply(ac, "[[", 9)))),c(chr.pos.all))); ac.mat[,1:10]
# # populate with allele counts
# for (i in 1:length(listpaths)) {
#   ac.mat[i,] <- ac[[i]]$ac[match(colnames(ac.mat), ac[[i]]$chr.pos)]*2
# }
# # fill with 0s for those without data (potential deletions)
# ac.mat[is.na(ac.mat)] <- 0; ac.mat[,1:10]
# 
# # prep sample size matrix
# n.mat <- ac.mat
# for (i in 1:length(listpaths)) {
#   n.mat[i,] <- ac[[i]]$nInd[match(colnames(n.mat), ac[[i]]$chr.pos)]*2
#   #n.mat[i,][is.na(n.mat[i,])] <- unique(sort((n.mat[i,])))[1]*2
# }
# n.mat [is.na(n.mat)] <- 0
# n.mat[,1:10]
# 
# ## remove sites with missing data in any population ###
# ## fst wont use sites with missing data (bedassle not sure)
# n.mat <- n.mat[,colSums(n.mat==0)==0]
# ac.mat <- ac.mat[,colSums(ac.mat==0)==0]

# ## REMOVE COL PATTERN AND INV ###
# ac.mat <-ac.mat[,c(substr(colnames(ac.mat), 1,8)!="Herato02"&substr(colnames(ac.mat), 1,8)!="Herato10"&
#                      substr(colnames(ac.mat), 1,8)!="Herato13"&substr(colnames(ac.mat), 1,8)!="Herato15"&substr(colnames(ac.mat), 1,8)!="Herato18")]
# n.mat <-n.mat[,c(substr(colnames(n.mat), 1,8)!="Herato02"&substr(colnames(n.mat), 1,8)!="Herato10"&
#                      substr(colnames(n.mat), 1,8)!="Herato13"&substr(colnames(n.mat), 1,8)!="Herato15"&substr(colnames(n.mat), 1,8)!="Herato18")]

write.table(n.mat, paste0("./n.mat.unlinked.", AREA, ".txt"), row.names=T, col.names=NA, sep="\t", quote=F)
write.table(ac.mat, paste0("./ac.mat.unlinked.", AREA, ".txt"), row.names=T, col.names=NA, sep="\t", quote=F)


## FST #####
print("run fst")
fst.pop.mat <- calculate.all.pairwise.Fst(ac.mat, n.mat) ; fst.pop.mat
rownames(fst.pop.mat) <- unique(unlist(lapply(ac.sub.2k, "[[", 9))); 
colnames(fst.pop.mat) <- unique(unlist(lapply(ac.sub.2k, "[[", 9))); fst.pop.mat
write.table(fst.pop.mat , paste0("./fst.pop.mat.unlinked.", AREA, ".txt"), row.names=T, col.names=T, sep="\t", quote=F)

# ## FST test, just 100 sites #####
# print("run fst with 1000 sites")
# fst.pop.mat.1000 <- calculate.all.pairwise.Fst(ac.mat[sample(colnames(ac.mat), 1000)], n.mat[sample(colnames(n.mat), 1000)]) ; fst.pop.mat.1000
# write.table(fst.pop.mat.1000 , paste0("./fst.pop.mat.1000.", AREA, ".txt"), row.names=T, col.names=NA, sep="\t", quote=F)

# ## FST with deletions #####
# print("run fst including deletions")
# #change n.mat so that all have nInd
# for (i in 1:length(listpaths)) {
#   n.mat[i,][n.mat[i,]==0] <- unique(n.mat[i,][n.mat[i,]!=0])[1]
# }
# n.mat[,1:10]
# fst.pop.mat.dels <- calculate.all.pairwise.Fst(ac.mat, n.mat) ; fst.pop.mat.dels
# write.table(fst.pop.mat , paste0("./fst.pop.mat.dels.", AREA, ".txt"), row.names=T, col.names=NA, sep="\t", quote=F)

