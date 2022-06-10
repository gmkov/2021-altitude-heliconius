#!/bin/bash
#SBATCH -p skylake
#SBATCH --time=24:00:00
#SBATCH --ntasks=1 
#SBATCH --cpus-per-task=32
#SBATCH --mail-type=ALL

###this script will work on bamfiles by population and calculate saf, then thetas, then windows
#maybe edit
POP=$1 #input file name list


#prepare variables - avoid to modify
source ../../01.scripts/01_config.sh

#Do thetas on all samples
echo "estimate thetas for all samples"

REF_era=/home/mgm49/Hera1/Heliconius_erato_demophoon_v1_-_scaffolds.fa
REF_mel=/home/mgm49/rds/rds-cj107-jiggins-rds/genomes/Hmel/Hmel2.5.fa

N_IND=$(wc -l ../../02.info/pop.list/"$POP".txt | cut -d " " -f 1)
SP=$(echo $POP | cut -d "." -f 1)
REF="REF_$SP"

#unsure whether we should specify the sites & filter for the thetas too??
echo "working on pop $POP, $N_IND individuals of the species $SP and reference ${!REF}"
#we need to re-do doSaf because we don't want to filter on maf for thetas calculation

angsd -P 32 \
-anc ${!REF}  -ref ${!REF} \
-dosaf 1 -GL 1 -doMajorMinor 1 \
-remove_bads 1 -baq 2 -minQ 20 \
-b ../../02.info/pop.list/"$POP".txt -out ./"$POP"

echo "estimate real sfs for pop $i"
realSFS -fold 1 ./"$POP"*.saf.idx -P 32 > ./"$POP".sfs

echo "estimate thetas for pop $POP"
realSFS saf2theta -P 32 ./"$POP"*.saf.idx -sfs ./"$POP".sfs -outname ./"$POP"

#Estimate for every Chromosome/scaffold
thetaStat do_stat ./"$POP".thetas.idx

#Do a sliding window analysis
echo "estimate thetas for pop $POP"
thetaStat do_stat ./"$POP".thetas.idx -win $WINDOW -step $WINDOW_STEP \
-outnames ./"$POP".thetas.5kb.1kb

