#!/bin/bash
#SBATCH -p skylake
#SBATCH --time=18:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mail-type=ALL


###this script will work on all bamfiles and calculate saf, maf & genotype likelihood

#prepare variables - avoid to modify
source ../../../01.scripts/01_config.sh
PERCENT_IND=0.5
# input
#POP_LIST=$1 # pop, e.g. "era.ec.e.narupa"
POP_LIST=$(basename $(dirname "$PWD")).$(basename "`pwd`")

# based on the directory and parents name, look for correct pop bam list
BAM_LIST=../../../02.info/pop.list/"$POP_LIST".txt
DEPTH_POP=$(cat ../../../02.info/pop.depth/"$POP_LIST".txt) # mean_depth pop
SP=$(basename $(dirname "$PWD"))

#prepare variables - avoid to modify
N_IND=$(wc -l $BAM_LIST | cut -d " " -f 1)
MIN_IND_FLOAT=$(echo "($N_IND * $PERCENT_IND)"| bc -l) #percent ind 0.2
MIN_IND=${MIN_IND_FLOAT%.*} 
MAX_DEPTH_FLOAT=$(echo "($DEPTH_POP * $MAX_DEPTH_FACTOR * $N_IND)" |bc -l) #max depth factor 2
MAX_DEPTH=${MAX_DEPTH_FLOAT%.*} 
MIN_DEPTH_FLOAT=$(echo "(2 * $N_IND)" |bc -l) 
MIN_DEPTH=${MIN_DEPTH_FLOAT%.*} 
MIN_MAF=0.05
 
echo "working on pop $POP_LIST, which has $N_IND individuals, and mean depth $DEPTH_POP"
echo "will filter for sites with at least one read in $MIN_IND individuals, which is $PERCENT_IND of the total"

# Do saf/maf for all population listed, loop if
REF_era=/home/mgm49/Hera1/Heliconius_erato_demophoon_v1_-_scaffolds.fa
REF_mel=/home/mgm49/Hmel2.5/Hmel2.5.scaffolds.fa
REF="REF_$SP"

echo "choose what ref to use accoding to species name ($SP), in this case $REF"
# ${!} needed to to expand the variable

angsd -P 32 -bam $BAM_LIST -ref ${!REF} \
-doMaf 1 -doMajorMinor 4 -GL 1 -doCounts 1 \
-skipTriallelic 1 -rmTriallelic 1 -SNP_pval 0.000001 \
-minMaf $MIN_MAF \
-minQ 20 -setMinDepth $MIN_DEPTH -setMaxDepth $MAX_DEPTH \
-uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 \
-doGlf 2 \
-out "$POP_LIST".filtered  \

#extract SNP which passed the PERCENT_IND filters & their Major-minor alleles
#order the sites file by chromosome names 
#makes a region file matching the sites files and with same order
#index sites file

gunzip -f "$POP_LIST".filtered.mafs.gz 

echo "from the maf file, calculate MAC, extract a list of SNP chr, position, major all, minor all"
Rscript ../../../01.scripts/Rscripts/make_sites_list_bedassle_mafs.R "$POP_LIST"
angsd sites index ../../../02.info/sites/sites.bedassle.maf."$POP_LIST".txt

