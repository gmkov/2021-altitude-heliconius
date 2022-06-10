#!/bin/bash
#SBATCH -p skylake-himem
#SBATCH --time=08:00:00
#SBATCH --ntasks=1 
#SBATCH --cpus-per-task=6
#SBATCH --mail-type=ALL

###this script will work on all bamfiles and calculate saf, maf & genotype likelihood
# input will be shared.region.id (e.g. 1)
POP=$1

###this script will work on all bamfiles and calculate saf, maf & genotype likelihood

#prepare variables - avoid to modify
source 01.scripts/01_config.sh

# based on the directory and parents name, look for correct pop bam list
# for running in east/west
SIDE=$(basename "`pwd`")
SP=$(basename $(dirname "$PWD"))
BAM_LIST=/home/mgm49/rds/rds-cj107-jiggins-rds/projects/mgm49/era.mel.altitude/02.info/pop.list/"$SP"."$SIDE".txt
DEPTH_POP=$(cat /home/mgm49/rds/rds-cj107-jiggins-rds/projects/mgm49/era.mel.altitude/02.info/pop.depth/"$SP"."$SIDE".txt)

#prepare variables - avoid to modify
N_IND=$(wc -l $BAM_LIST | cut -d " " -f 1)
MAX_DEPTH_FLOAT=$(echo "($DEPTH_POP * $MAX_DEPTH_FACTOR * $N_IND)" |bc -l)
MAX_DEPTH=${MAX_DEPTH_FLOAT%.*} 
MIN_DEPTH_FLOAT=$(echo "(2 * $N_IND)" |bc -l) 
MIN_DEPTH=${MIN_DEPTH_FLOAT%.*} 

echo "working on pop $POP_LIST, which has $N_IND individuals, and mean depth $DEPTH_POP"
echo "will filter for sites with at least one read in $MIN_IND individuals, which is $PERCENT_IND of the total"

# Do saf/maf for all population listed, loop if
REF_era=/home/mgm49/Hera1/Heliconius_erato_demophoon_v1_-_scaffolds.fa
REF_mel=REF=/home/mgm49/rds/rds-cj107-jiggins-rds/genomes/Hmel/Hmel2.5.fa
REF="REF_$SP"

# take specific list of sites (and corresponding regions) of SHDR
SITES=/home/mgm49/rds/rds-cj107-jiggins-rds/projects/mgm49/era.mel.altitude/02.info/sites/shared.exp50kb."$SP"."$SIDE".std4/sites.shared.region."$POP".txt
REGIONS=/home/mgm49/rds/rds-cj107-jiggins-rds/projects/mgm49/era.mel.altitude/02.info/sites/shared.exp50kb."$SP"."$SIDE".std4/regions.shared.region."$POP".txt

echo " Calculate the beagle GL for all individuals od $SP"
echo "filter out sites with more than $MEL_MEAN_DEPTH * $MAX_DEPTH_FACTOR * $N_IND = $MAX_DEPTH"

#### beagle
angsd -P 6 -anc ${!REF}  -ref ${!REF} \
-rf $REGIONS -sites $SITES  \
-GL 1 -doGlf 2 -doMajorMinor 4 -doMaf 1 -doCounts 1 \
-remove_bads 1 -baq 2 -minQ 20 -setMinDepth $MIN_DEPTH -setMaxDepth $MAX_DEPTH  \
-b $BAM_LIST  \
-out "$SP"."$SIDE".region."$POP"

### pcangsd
module load python-3.5.2-gcc-5.4.0-rdp6q5l gcc

echo "analyse covariance matrix on all individuals"
python3 ~/programmes/pcangsd/pcangsd.py -threads 6 \
	-beagle "$SP"."$SIDE".region."$POP".beagle.gz -o ./"$SP"."$SIDE".region."$POP"

## plot pca or  do locally
#echo "transform covariance matrix into PCA"
#Rscript 01.scripts/Rscripts/make_pca_region.R "$SP" "$SIDE" "$POP" "$BAM_LIST"

