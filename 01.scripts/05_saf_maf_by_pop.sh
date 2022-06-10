#!/bin/bash
#SBATCH -p skylake
#SBATCH --time=24:00:00
#SBATCH --ntasks=1 
#SBATCH --cpus-per-task=6
#SBATCH --mail-type=ALL

###this script will work on bamfiles by population and calculate saf  & maf 
# input
POP=$1 # first input after calling script e.g."hig"

#prepare variables - avoid to modify
source ../../../01.scripts/01_config.sh

# based on the directory and parents name, look for correct pop bam list
BAM_LIST=../../../02.info/pop.list/$(basename $(dirname "$PWD")).$(basename "`pwd`").$POP.txt
POP_LIST=$(basename $(dirname "$PWD")).$(basename "`pwd`").$POP
DEPTH_POP=$(cat ../../../02.info/pop.depth/"$POP_LIST".txt)
SP=$(basename $(dirname "$PWD"))

#prepare variables - avoid to modify
N_IND=$(wc -l $BAM_LIST | cut -d " " -f 1)
MAX_DEPTH_FLOAT=$(echo "($DEPTH_POP * $MAX_DEPTH_FACTOR * $N_IND)" |bc -l)
MAX_DEPTH=${MAX_DEPTH_FLOAT%.*} 

echo "working on pop $POP_LIST, which has $N_IND individuals, and mean depth $DEPTH_POP"

# Do saf/maf for all population listed, loop if
REF_era=/home/mgm49/Hera1/Heliconius_erato_demophoon_v1_-_scaffolds.fa
REF_mel=/home/mgm49/rds/rds-cj107-jiggins-rds/genomes/Hmel/Hmel2.5.fa
REF="REF_$SP"

echo "choose what ref to use accoding to species name ($SP), in this case $REF"
# ${!} needed to to expand the variable

angsd -P 6 \
-anc ${!REF}  -ref ${!REF} \
-dosaf 1 -gl 1 -doCounts 1 \
-remove_bads 1  -minQ 20 -baq 2 -minMapQ 10 \
-setMaxDepth $MAX_DEPTH \
-b $BAM_LIST -out "$POP_LIST"_maxdepth"$MAX_DEPTH"
