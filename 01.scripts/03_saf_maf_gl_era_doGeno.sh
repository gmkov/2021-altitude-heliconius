#!/bin/bash

###this script will work on all bamfiles and calculate saf, maf & genotype likelihood

#prepare variables - avoid to modify
source /home/mgm49/rds/rds-cj107-jiggins-rds/projects/mgm49/era.mel.altitude/01.scripts/01_config.sh


# Important: Move to directory where job was submitted
# cd $SLURM_SUBMIT_DIR
PERCENT_IND=0.5

#prepare variables - avoid to modify
N_IND=$(wc -l /home/mgm49/rds/rds-cj107-jiggins-rds/projects/mgm49/era.mel.altitude/02.info/era.bam.list | cut -d " " -f 1)
MIN_IND_FLOAT=$(echo "($N_IND * $PERCENT_IND)"| bc -l)
MIN_IND=${MIN_IND_FLOAT%.*} 
MAX_DEPTH=$(echo "($ERA_MEAN_DEPTH * $MAX_DEPTH_FACTOR * $N_IND)" |bc -l)
SIDE=$1

REF=/home/mgm49/Hera1/Heliconius_erato_demophoon_v1_-_scaffolds.fa

echo " Calculate the SAF, MAF and GL for all individuals listed in 02_info/era.bam.list"
echo "filter on allele frequency = $MIN_MAF_ERA"
echo "filter out sites with more than $ERA_MEAN_DEPTH * $MAX_DEPTH_FACTOR * $N_IND = $MAX_DEPTH"
echo "downsample genome to $DOWNSAMPLING_PERC"

##Calculate the SAF, MAF and GL
angsd -P 15 \
-doMaf 1 -dosaf 1 -GL 2 -doGlf 2 -doMajorMinor 1 -doCounts 1 \
-doGeno 2 -doPost 1 -doDepth 1 -maxDepth $MAX_DEPTH -dumpCounts 2 \
-anc $REF -remove_bads 1 -minMapQ 30 -minQ 20 \
-minInd $MIN_IND -minMaf $MIN_MAF_ERA -setMaxDepth $MAX_DEPTH -downSample $DOWNSAMPLING_PERC \
-b ../../02.info/era."$SIDE".bam.list \
-out ../../03.saf.maf.gl.all/era/"$SIDE"_maf"$MIN_MAF_ERA"_pctind"$PERCENT_IND"_maxdepth"$MAX_DEPTH"_downsample"$DOWNSAMPLING_PERC"


#extract SNP which passed the PERCENT_IND filters & their Major-minor alleles
#order the sites file by chromosome names 
#makes a region file matching the sites files and with same order
#index sites file
echo "from the maf file, calculate MAC, extract a list of SNP chr, position, major all, minor all"

gunzip ../../03.saf.maf.gl.all/era/all_maf"$MIN_MAF_ERA"_pctind"$PERCENT_IND"_maxdepth"$MAX_DEPTH"_downsample"$DOWNSAMPLING_PERC".mafs.gz 
Rscript ../../01.scripts/Rscripts/make_sites_list_maxdepth_era.R "$MIN_MAF_ERA" "$PERCENT_IND" "$MAX_DEPTH" "$DOWNSAMPLING_PERC"
angsd sites index ../../02.info/sites_era_maf"$MIN_MAF_ERA"_pctind"$PERCENT_IND"_maxdepth"$MAX_DEPTH"_downsample"$DOWNSAMPLING_PERC"


