# want to grab positions at least 10
# allow shift len to vary until finds a row that meets the criteria
# input name of sites list to subset, and the unlinked position spacing (in bp)
# prep data
argv <- commandArgs(T)
df.input.name <- argv[1]
pos.diff.input <- argv[2]; pos.diff.input <-  as.numeric(as.character(pos.diff.input))

df <- read.table(paste0("./", df.input.name), header = F); head(df)
names(df)[1] <- "scaff"
names(df)[2] <- "pos"
df$pos <- as.numeric(as.character(df$pos)); str(df)

# subset by 
library(dplyr)
nrows <- 1
previous_match <- 1
for(i in 2:nrow(df)) {
  # condition 1: that the difference between row i and the previous match is more than pos.diff
  # OR condition 2: that the difference between row i and the previous match is <0 , non-ascending (for scaffold changes)
  if(df$pos[i] - df$pos[previous_match] > pos.diff.input | df$pos[i] - df$pos[previous_match] < 0) {
    nrows <- c(nrows, i)
    previous_match <- i
  }
}
df1 <- df[nrows, ]; head(df1)

sites.df <- df1[order(df1$scaff),] 
regions.df <-unique(sites.df$scaff)


# save
write.table(sites.df , paste0("./", substr(df.input.name, 0, 24), ".sub.", pos.diff.input, ".txt"), col.names = F, row.names=F, sep="\t", quote=F)
write.table(regions.df, paste0("./regions.", substr(df.input.name, 7, 24), ".sub.", pos.diff.input, ".txt"), col.names = F, row.names=F, sep="\t", quote=F)

#testing
 # df1 <- data.frame(scaff=rep(1,1000), x=c(1:1000), pos=sort(sample(1:10000, 1000))); df1
 # df2 <- data.frame(scaff=rep(2,1000), x=c(1:1000), pos=sort(sample(1:10000, 1000))); df2
 # df.input <- rbind(df1, df2)
 # head(df.input, n = 20); pos.diff.input <-10
 # head(subset(df.input, c(TRUE, diff(pos) > pos.diff.input)), n = 20)
 # head(df.input %>% filter(pos - lag(pos, default = -Inf) > 10), n = 20)
 
