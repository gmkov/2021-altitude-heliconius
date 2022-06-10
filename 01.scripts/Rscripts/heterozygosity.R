# from angsd sfs (ml) calculate global hete
# allow shift len to vary until finds a row that meets the criteria
# input name of sites list to subset, and the unlinked position spacing (in bp)
# prep data
argv <- commandArgs(T)
pop.name <- argv[1]


sfs <- scan(paste0("./", pop.name, ".sfs"))
sfs[2]/sum(sfs) # should print result
