#!/bin/bash
#SBATCH -p skylake-himem
#SBATCH --time=08:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mail-type=ALL

#prepare variables - avoid to modify
source ../../../01.scripts/01_config.sh

# based on the directory and parents name, look for correct pop bam list
POP_LIST=$1
POP=$(echo $1 | cut -c10-12)  # subpop, e.g. "nar" 
AREA=$(basename $(dirname "$PWD"))$(echo $POP_LIST | cut -c7-8) # era.e
BAM_LIST=../../../02.info/pop.list/"$POP_LIST".txt
DEPTH_POP=$(cat ../../../02.info/pop.depth/"$POP_LIST".txt) # mean_depth pop
SP=$(basename $(dirname "$PWD"))

#prepare variables - avoid to modify
N_IND=$(wc -l $BAM_LIST | cut -d " " -f 1)
MIN_IND_FLOAT=$(echo "($N_IND * $PERCENT_IND)"| bc -l) #percent ind 0.2
MIN_IND=${MIN_IND_FLOAT%.*} 
MAX_DEPTH_FLOAT=$(echo "($DEPTH_POP * $MAX_DEPTH_FACTOR * $N_IND)" |bc -l) #max depth factor 3
MAX_DEPTH=${MAX_DEPTH_FLOAT%.*} 
MIN_DEPTH_FLOAT=$(echo "(2 * $N_IND)" |bc -l) 
MIN_DEPTH=${MIN_DEPTH_FLOAT%.*} 
MIN_MAF=$(echo "(1/(2*$N_IND))" |bc -l)

echo "obtain MAFs from $POP_LIST with the bam list found in: $BAM_LIST and selected sites"

SITES=../../../02.info/sites/sites.bedassle.maf."$AREA".txt
REGIONS=../../../02.info/sites/regions.bedassle.maf."$AREA".txt

# Do saf/maf for all population listed, loop if
REF_era=/home/mgm49/Hera1/Heliconius_erato_demophoon_v1_-_scaffolds.fa
REF_mel=/home/mgm49/Hmel2.5/Hmel2.5.scaffolds.fa
REF="REF_$SP"

echo "choose what ref to use accoding to species name ($SP), in this case $REF"
# ${!} needed to to expand the variable

angsd -P 1  \
-doMaf 1 -doMajorMinor 3 -GL 1 -doCounts 1  \
-rf $REGIONS -sites $SITES -bam $BAM_LIST -out "$POP_LIST"  \

## create allele count
# gunzip force to overwrite
gunzip -f "$POP_LIST".mafs.gz 
# obtain allele counts from allele frequencies 
#(by multiplying the knownEM(known freq) and Nind(number of individuals with data) at each loci
awk '$(NF+1)=$5*$6' "$POP_LIST".mafs > "$POP_LIST".mafs.ac 




