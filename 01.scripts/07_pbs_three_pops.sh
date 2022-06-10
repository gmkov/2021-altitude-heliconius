#!/bin/bash
#SBATCH -o log_%j
#SBATCH -c 1
#SBATCH -p medium
#SBATCH --mail-type=ALL
#SBATCH --time=3-00:00

## this script will run pbs with three populations given as input
## then run sliding window analyses
# input
POP1=$1 # first input pop 
POP2=$2
POP3=$3

#prepare variables - avoid to modify
source ../../../01.scripts/01_config.sh

# based on the directory and parents name, look for correct pop bam list
POP_LIST1=$(basename $(dirname "$PWD")).$(basename "`pwd`").$POP1
DEPTH_POP1=$(cat ../../../02.info/pop.depth/"$POP_LIST1".txt)
BAM_LIST_POP1=../../../02.info/pop.list/$(basename $(dirname "$PWD")).$(basename "`pwd`").$POP1.txt
SP=$(basename $(dirname "$PWD"))

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

# based on the directory and parents name, look for correct pop bam list
POP_LIST3=$(basename $(dirname "$PWD")).$(basename "`pwd`").$POP3
DEPTH_POP3=$(cat ../../../02.info/pop.depth/"$POP_LIST3".txt)
BAM_LIST_POP3=../../../02.info/pop.list/$(basename $(dirname "$PWD")).$(basename "`pwd`").$POP3.txt
#prepare variables - avoid to modify
N_IND_POP3=$(wc -l $BAM_LIST_POP3 | cut -d " " -f 1)
MAX_DEPTH_FLOAT=$(echo "($DEPTH_POP3 * $MAX_DEPTH_FACTOR * $N_IND_POP3)" |bc -l)
MAX_DEPTH_POP3=${MAX_DEPTH_FLOAT%.*} 

# based on the directory and parents name, look for correct pop area (e.g. ec.e)
POP_AREA=$(basename $(dirname "$PWD")).$(basename "`pwd`")

# calculate PBS
echo "calculate PBS for $POP_LIST1 $POP_LIST2 $POP_LIST3"
realSFS fst index "$POP_AREA"."$POP1"_maxdepth"$MAX_DEPTH_POP1".saf.idx "$POP_AREA"."$POP2"_maxdepth"$MAX_DEPTH_POP2".saf.idx "$POP_AREA"."$POP3"_minind"$MIN_IND_POP3"_maxdepth"$MAX_DEPTH_POP3".saf.idx -sfs "$POP1"."$POP2".ml -sfs "$POP1"."$POP3".ml -sfs "$POP2"."$POP3".ml -fstout "$POP_AREA"."$POP1"."$POP2"."$POP3".pbs -whichFst 1

create windows
echo "calculate sliding window means for $POP_AREA.pbs"
realSFS fst stats2 "$POP_AREA"."$POP1"."$POP2"."$POP3".pbs.fst.idx -win $WINDOW -step $WINDOW_STEP -type 0 -whichFst 1 > "$POP_AREA"."$POP1"."$POP2"."$POP3".pbs.5kb.1kb.idx

echo "done with windows"

# if melpomene, liftover positions from erato ref
if [ $SP=='mel' ]
then
	file="$POP_AREA"."$POP1"."$POP2"."$POP3".pbs.5kb.1kb.idx
	chain_file=/home/mgm49/rds/rds-cj107-jiggins-rds/genomes/Hmel/hmel2.5-helera1_demo.largeGap.14k.filtered.chain.gz
	
	liftOver <(awk 'NR>1{print $2,$3,$3+1000,$2"_"$3"_"$4"_"$5"_"$6"_"$7"_"$8"_"$9"_"$10}' "$file") \
	"$chain_file" "$file".eratoLiftOver.tmp "$file".eratoLiftOver.unmapped -multiple 
	
	sed -e '1ichr\tmidPos\tmidPos.plus1000\tHmelchr\tHmelmidPos\tNsites\tFst01\tFst02\tFst12\tPBS0\tPBS1\tPBS2\tindex' \
	-e 's/_/\t/g' "$file".eratoLiftOver.tmp | awk '{if($13<2 || $13=="index") print}' > "$file".eratoLiftOver
	
	rm $file.eratoLiftOver.tmp 
fi

echo "done with liftover"






