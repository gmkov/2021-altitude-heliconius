#filter : will keep SNP above this allele frequency (over all individuals)
# best to use alelle frequency counts instead
# 1/(2*nSamples)
MIN_MAF_era=0.0015
MIN_MAF_mel=0.0025

# filter: will keep SNP with at least a coverage of this factor multiplied by the number of ind - across all ind. usually set 2-4 (Chao et al 19
# times the coverage to remove repeated regions
# we will use pre-calculated mean cov. per pop
MAX_DEPTH_FACTOR=2

# downsampling factor, percentage of genome to keep for some analyses
DOWNSAMPLING_PERC=0.01

# and mean all era/mel depth for first maf filtering
ERA_MEAN_DEPTH=8.96
MEL_MEAN_DEPTH=12.4

# window size for sliding window FST & Thetas
WINDOW=5000

#window step
WINDOW_STEP=1000

# recombination rates
RECOMB_mel=3.7
RECOMB_era=3