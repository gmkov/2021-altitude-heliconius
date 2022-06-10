#!/bin/bash
#SBATCH -p skylake-himem
#SBATCH --time=12:00:00
#SBATCH --ntasks=1 
#SBATCH --cpus-per-task=32
#SBATCH --mail-type=ALL

module load python-3.5.2-gcc-5.4.0-rdp6q5l gcc

# this script will work on all individuals using the beagle genotype likelihood and calculate a covariance matrix with angsd & a pca with R
# this requires pcangsd to be cloned and a version of Python v2 with alias python2

#prepare variables - avoid to modify
source /home/mgm49/rds/rds-cj107-jiggins-rds/projects/mgm49/era.mel.altitude/01.scripts/01_config.sh
PERCENT_IND=0.5
SP=$1 #era
SIDE=$2 #e

#prepare variables - avoid to modify
N_IND=$(wc -l /home/mgm49/rds/rds-cj107-jiggins-rds/projects/mgm49/era.mel.altitude/02.info/pop.list/"$SP"."$SIDE".txt | cut -d " " -f 1)
MIN_IND_FLOAT=$(echo "($N_IND * $PERCENT_IND)"| bc -l)
MIN_IND=${MIN_IND_FLOAT%.*} 
MAX_DEPTH=$(echo "($ERA_MEAN_DEPTH * $MAX_DEPTH_FACTOR * $N_IND)" |bc -l)

#use beagle input from bedassle, already filtered with snp val
INPUT=/home/mgm49/rds/rds-cj107-jiggins-rds/projects/mgm49/era.mel.altitude/11.bedassle/"$SP"/"$SIDE"/"$SP"."$SIDE".filtered.beagle.gz

echo "analyse covariance matrix on all individuals"
python3 ~/programmes/pcangsd/pcangsd.py -threads 2 \
	-beagle $INPUT -o /home/mgm49/rds/rds-cj107-jiggins-rds/projects/mgm49/era.mel.altitude/04.pca/sides/"$SP"."$SIDE".bedassle.filtered
