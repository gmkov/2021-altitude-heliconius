#!/bin/bash
#SBATCH -p skylake
#SBATCH --time=23:00:00
#SBATCH --ntasks=1 
#SBATCH --cpus-per-task=32
#SBATCH --mail-type=ALL
#SBATCH --output=R-%x.%j.out

### this script will get dxy from thetas

JOINT=$1 #input file name list
POPA=$(echo $JOINT | cut -d "." -f 1,2,3,4)
POPB=$(echo $JOINT | cut -d "." -f 5,6,7,8)
POPAB=$(echo $JOINT | cut -d "." -f 1,2,3,4,5,6,7,8)


#prepare variables 
source /home/mgm49/rds/rds-cj107-jiggins-rds/projects/mgm49/era.mel.altitude/01.scripts/01_config.sh

REF_era=/home/mgm49/Hera1/Heliconius_erato_demophoon_v1_-_scaffolds.fa
REF_mel=/home/mgm49/rds/rds-cj107-jiggins-rds/genomes/Hmel/Hmel2.5.fa

N_INDA=$(wc -l /home/mgm49/rds/rds-cj107-jiggins-rds/projects/mgm49/era.mel.altitude/02.info/pop.list/"$POPA".txt | cut -d " " -f 1)
N_INDB=$(wc -l /home/mgm49/rds/rds-cj107-jiggins-rds/projects/mgm49/era.mel.altitude/02.info/pop.list/"$POPB".txt | cut -d " " -f 1)
N_INDAB=$(wc -l /home/mgm49/rds/rds-cj107-jiggins-rds/projects/mgm49/era.mel.altitude/02.info/pop.list/"$POPAB".txt | cut -d " " -f 1)

echo "use thetas from $POPA, $POPB, and joint $POPAB, with $N_INDA and $N_INDB individuals, totaling $N_INDAB"
module load R
Rscript /home/mgm49/rds/rds-cj107-jiggins-rds/projects/mgm49/era.mel.altitude/01.scripts/Rscripts/calc.dxy.gmk.shm.R \
--popA $POPA --popB $POPB --popAB $POPAB  \
--nIndA $N_INDA --nIndB $N_INDB --nIndAB $N_INDAB \
-w 5 -s 1 -t 1



