#!/bin/bash
module load python-3.5.2-gcc-5.4.0-rdp6q5l gcc

# this script will work on all individuals using the beagle genotype likelihood and calculate a covariance matrix with angsd & a pca with R
# this requires pcangsd to be cloned and a version of Python v2 with alias python2

#prepare variables - avoid to modify
source ../../01.scripts/01_config.sh
PERCENT_IND=0.5
#prepare variables - avoid to modify
N_IND=$(wc -l ../../02.info/era.bam.list | cut -d " " -f 1)
MIN_IND_FLOAT=$(echo "($N_IND * $PERCENT_IND)"| bc -l)
MIN_IND=${MIN_IND_FLOAT%.*} 
MAX_DEPTH=$(echo "($ERA_MEAN_DEPTH * $MAX_DEPTH_FACTOR * $N_IND)" |bc -l)

#this is the list of file we are working on
BAM_LIST=../../02.info/era.bam.list

#this is the input file for the pca
INPUT=../../03.saf.maf.gl.all/era/all_maf"$MIN_MAF_ERA"_pctind"$PERCENT_IND"_maxdepth"$MAX_DEPTH"_downsample"$DOWNSAMPLING_PERC".beagle.gz

echo "analyse covariance matrix on all individuals"
python3 ~/programmes/pcangsd/pcangsd.py \
	-beagle $INPUT -o ../../04.pca/era/all_maf"$MIN_MAF_ERA"_pctind"$PERCENT_IND"_maxdepth"$MAX_DEPTH"_downsample"$DOWNSAMPLING_PERC"

echo "transform covariance matrix into PCA"
Rscript ../../01.scripts/Rscripts/make_pca_era_maxdepth.R "$MIN_MAF_ERA" "$PERCENT_IND" "$BAM_LIST" "$MAX_DEPTH" "$DOWNSAMPLING_PERC"

