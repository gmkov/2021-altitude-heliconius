#!/bin/bash
#SBATCH -p skylake-himem
#SBATCH --time=12:00:00
#SBATCH --ntasks=1 
#SBATCH --cpus-per-task=32
#SBATCH --mail-type=ALL

###this script will use the saf by population in step 02 and calculate 2dsfs
# input
POP1=$1 # first input after calling script e.g."hig"
POP2=$2

#prepare variables - avoid to modify
source ../../../01.scripts/01_config.sh

# based on the directory and parents name, look for correct pop bam list
POP_LIST1=$(basename $(dirname "$PWD")).$(basename "`pwd`").$POP1
DEPTH_POP1=$(cat ../../../02.info/pop.depth/"$POP_LIST1".txt)
BAM_LIST_POP1=../../../02.info/pop.list/$(basename $(dirname "$PWD")).$(basename "`pwd`").$POP1.txt
#prepare variables - avoid to modify
N_IND_POP1=$(wc -l $BAM_LIST_POP1 | cut -d " " -f 1)
MAX_DEPTH_FLOAT=$(echo "($DEPTH_POP1 * $MAX_DEPTH_FACTOR * $N_IND_POP1)" |bc -l)
MAX_DEPTH_POP1=${MAX_DEPTH_FLOAT%.*} 

# based on the directory and parents name, look for correct pop bam list
POP_LIST2=$(basename $(dirname "$PWD")).$(basename "`pwd`").$POP2
DEPTH_POP2=$(cat ../../../02.info/pop.depth/"$POP_LIST2".txt)
BAM_LIST_POP2=../../../02.info/pop.list/$(basename $(dirname "$PWD")).$(basename "`pwd`").$POP2.txt
#prepare variables - avoid to modify
N_IND_POP2=$(wc -l $BAM_LIST_POP2 | cut -d " " -f 1)
MAX_DEPTH_FLOAT=$(echo "($DEPTH_POP2 * $MAX_DEPTH_FACTOR * $N_IND_POP2)" |bc -l)
MAX_DEPTH_POP2=${MAX_DEPTH_FLOAT%.*} 

# based on the directory and parents name, look for correct pop area (e.g. ec.e)
POP_AREA=$(basename $(dirname "$PWD")).$(basename "`pwd`")

echo "doing 2dsfs of $POP1 and $POP2 from $POP_AREA, which has $N_IND_POP1 and $N_IND_POP2 individuals"

realSFS  -P 32 "$POP_AREA"."$POP1"_maxdepth"$MAX_DEPTH_POP1".saf.idx \
"$POP_AREA"."$POP2"_maxdepth"$MAX_DEPTH_POP2".saf.idx > "$POP1"."$POP2".ml


