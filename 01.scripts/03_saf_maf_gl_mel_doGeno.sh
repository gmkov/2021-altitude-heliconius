#!/bin/bash

###this script will work on all bamfiles and calculate saf, maf & genotype likelihood

#prepare variables - avoid to modify
source ../../01.scripts/01_config.sh
PERCENT_IND=0.5

#prepare variables - avoid to modify
N_IND=$(wc -l ../../02.info/mel.bam.list | cut -d " " -f 1)
MIN_IND_FLOAT=$(echo "($N_IND * $PERCENT_IND)"| bc -l)
MIN_IND=${MIN_IND_FLOAT%.*} 
MAX_DEPTH=$(echo "($MEL_MEAN_DEPTH * $MAX_DEPTH_FACTOR * $N_IND)" |bc -l)

REF=/home/mgm49/Hmel2.5/Hmel2.5.scaffolds.fa

echo " Calculate the SAF, MAF and GL for all individuals listed in 02_info/mel.bam.list"
echo "keep loci with at leat one read for n individuals = $MIN_IND, which is $PERCENT_IND % of total $N_IND individuals"
echo "filter on allele frequency = $MIN_MAF_MEL"
echo "filter out sites with more than $MEL_MEAN_DEPTH * $MAX_DEPTH_FACTOR * $N_IND = $MAX_DEPTH"
echo "downsample genome to $DOWNSAMPLING_PERC"

####Calculate the SAF, MAF and GL
angsd -P 15 \
-doMaf 1 -dosaf 1 -GL 2 -doGlf 2 -doMajorMinor 1 -doCounts 1 -checkBamHeaders 0 \
-doGeno 2 -doPost 1 -doDepth 1 -maxDepth $MAX_DEPTH -dumpCounts 2 \
-anc $REF -remove_bads 1 -minMapQ 30 -minQ 20 \
-minInd $MIN_IND -minMaf $MIN_MAF_MEL -setMaxDepth $MAX_DEPTH -downSample $DOWNSAMPLING_PERC \
-b ../../02.info/mel.bam.list  \
-out ../../03.saf.maf.gl.all/mel/all_maf"$MIN_MAF_MEL"_pctind"$PERCENT_IND"_maxdepth"$MAX_DEPTH"_downsample"$DOWNSAMPLING_PERC"

#extract SNP which passed the MIN_MAF and PERCENT_IND filters & their Major-minor alleles
#order the sites file by chromosome names 
#makes a region file matching the sites files and with same order
#index sites file
echo "from the maf file, extract a list of SNP chr, position, major all, minor all"
gunzip ../../03.saf.maf.gl.all/mel/all_maf"$MIN_MAF_MEL"_pctind"$PERCENT_IND"_maxdepth"$MAX_DEPTH"_downsample"$DOWNSAMPLING_PERC".mafs.gz 
Rscript ../../01.scripts/Rscripts/make_sites_list_maxdepth_mel.R "$MIN_MAF_MEL" "$PERCENT_IND" "$MAX_DEPTH" "$DOWNSAMPLING_PERC"
angsd sites index ../../02.info/sites_mel_maf"$MIN_MAF_MEL"_pctind"$PERCENT_IND"_maxdepth"$MAX_DEPTH"_downsample"$DOWNSAMPLING_PERC"

