# Overall pipeline

This is an overview of the analyses, indicating which scripts are relevant for each section

Most of these bash scripts work with the directory names they are being run within, i.e. will take information from the path to know what population list to use (erato or melpomene, east or west etc.) as well as to know what reference to use, what filters etc. When the analyses are run at the replicate/cline level the directory structure looks like this, e.g. ‘era/co.e/’ is *H. erato* Colombia East,

```
├── era
│   ├── co.e
│   ├── co.w
│   ├── ec.e
│   └── ec.w
└── mel
    ├── co.e
    ├── co.w
    ├── ec.e
    └── ec.w
```

## 01_config.sh

Configuration file for cluster. Filters, downsampling percentage for PCA, mean depth of species

## 03_saf_maf_gl_XXX_doGeno.sh

Run ANGSD of all individuals of each species or per sides of the Andes, with conservative filters for retaining conserved sites for PCA.
Downsample randomly to 1% of the genome for PCA
Obtain sites above minMAC 2, export sites list


```
nohup nice bash 01.scripts/03_saf_maf_gl_era_doGeno.sh &> 03_saf_maf_gl_era_doGeno.out &
nohup nice bash 01.scripts/03_saf_maf_gl_mel_doGeno.sh &> 03_saf_maf_gl_mel_doGeno.out &


```

## 04_pca_era.sh

Input: beagle file fromstep 03 
After PCAngsd, run Rscript that gives row names to the covariance matrix, runs pca (saves as .pca > export for better plotting), initial plotting (pca.pdf)

```
# example running erato east
nohup nice bash 01.scripts/04_pca_sides.sh era e &> 04_pca_era.e.out &

```

Can also plot locally with local/scripts/5.2.fig5_shdr.chr2_pca_local.pca.R

## 05_saf_maf_by_pop.sh

SFS per population. Obtain SFS per populations with **01.scripts/05_saf_maf_by_pop.sh** 

```
# inside each species
mkdir ec.e; mkdir ec.w; mkdir co.e; mkdir co.w

# running example
nohup nice bash ../../../01.scripts/05_saf_maf_by_pop.sh hig &> 05_saf_maf_by_pop.hig.out & 
nohup nice bash ../../../01.scripts/05_saf_maf_by_pop.sh low &> 05_saf_maf_by_pop.low.out & 
nohup nice bash ../../../01.scripts/05_saf_maf_by_pop.sh vlo &> 05_saf_maf_by_pop.vlo.out &
 
```

## 06_2dsfs_by_pop_pair.sh

Calculates 2dsfs per population pair. needs two inputs (pop1 pop2), with realSFS. When three pops available do all combinations. Input population shortcut (‘hig’ will be all high individuals from that area, or ‘ama’ for all individuals from a population called ‘ama’)- this will use individuals in the corresponding list under [02.info/pop.list/](http://02.info/pop.list/) 


```
# running example
nohup nice bash ../../../01.scripts/06_2dsfs_by_pop_pair.sh hig low &> 06_2dsfs_by_pop_pair.sh.hig.low.out & 
nohup nice bash ../../../01.scripts/06_2dsfs_by_pop_pair.sh hig vlo &> 06_2dsfs_by_pop_pair.sh.hig.vlo.out & 
nohup nice bash ../../../01.scripts/06_2dsfs_by_pop_pair.sh low vlo &> 06_2dsfs_by_pop_pair.sh.low.vlo.out &
 
```

## 07_pbs_three_pops.sh or 07_fst_by_pop_pair.sh

Gets PBS when three populations are available or Fst if only 2 per cline.

```
# example for PBS hig low low distant populations, working within 07.output/era/ec.e/
nohup nice bash ../../../01.scripts/07_pbs_three_pops_maxdepth.sh hig low vlo &> 07_pbs_three_pops_maxdepth.out & 
# example for hig low fst
nohup nice bash ../../../01.scripts/07_fst_by_pop_pair.sh hig low &> 07_fst_by_pop_pair.out & 

```

## 08_thetas.pop.sh

Calculate sfs, then thetas, then windows

1. Obtain folded global site-frequency spectra for each population. 
2. Calculate pairwise nucleotide diversity per site (thetaD, realSFS saf2theta).
3. Perform sliding window analysis of 5kb window size and 1kb steps (thetaStat do_stat) to obtain sum of pairwise differences, Tajima’s D, and total effective number of sites per window. 
4. Nucleotide diversity (pi) was obtained by dividing the sum of pairwise differences by the total number of sites per window.

```
# example to run on list of pops
for i in $(less pops.list.txt)
do
nohup nice bash  ../../01.scripts/08_thetas.pop.sh  $i  &> 08.thetas.pop.$i.out & 
done



```

From this we can get global heterozygosity for each pop, listed in pop.txt (matching to [02.info/pop.list/](http://02.info/pop.list/) )

```
# our ml is saved as .sfs , check R
for i in $(less pops.txt)
do
nohup nice Rscript  ../../01.scripts/Rscripts/heterozygosity.R  $i  &> het.pop.$i.out & 
done

# extract 3rd line, append to file with all global hets
for i in $(ls het.sub*)
do
sed '3q;d' $i >> het.pop.txt
done


```

## 09_thetas.pop.dxy.sh

Absolute divergence (Dxy, Nei, 1978) between high (population A) and low/low distant (population B) populations was estimated by additionally obtaining pairwise nucleotide diversity per site (thetaD) for all individuals pooled from populations A and B (population AB), and then per-site Dxy obtained (see formula in main manuscript). 

Bash script 09_thetas.pop.dxy.sh will find 3 relevant thetas files from step 8, estimate number of individuals per pop (and sum of the two). R script 01/scripts/Rscripts/calc.dxy.gmk.shm.R will read in the 3 files and a load of options (see script for details), calculate per site dxy, then do windows with winscanr based on POSITION (0-5000), but counting the number of sites included (e.g. 4988)

To run you must input the ‘joint population’ pairs that are to be compared, this will match a list in 02.info/pop.list/; for example ‘era.co.e.cam.era.co.e.esc.txt’ corresponds to the individuals from two populations from erato colombia east ‘cam’ and ‘esc’.

```
for i in $(less joint.pops.txt)
do
nohup nice bash  ../../../01.scripts/08_thetas.pop.dxy.sh  $i  &> dxy.pop.$i.out & 
done

```

## 11. pairwise Fst with BEDASSLE (genetic dist matrix), geographic and environmental matrices

We will use ANGSD + custom scripts to obtain an allele count matrix per replicate (with populations as rows and alleles as columns). First we will obtain species-wide SNP list, then individual minor allele frequency at each SNP, then aggregate by population to obtain the allele count matrix.


1. run angsd with 10 random bam files per area to filter fairly heavily, obtain MAF file → sites/regions kept file per species/side : **01.scripts/11_doMaf.all.bedassle.sh** (which **** makes list of sites to guide angsd in next step **01.scripts/Rscripts/make_sites_list_bedassle_mafs.R**)
2. run population specific angsd to get maf (no filters) in sites/regions determined in step 1, use **01.scripts/11_doMaf.pop.bedassle.sh**
3. In R (use same sites for fst, unlinked and present in all pops),  use **01.scripts/Rscripts/mafs.ac.to.ac.pop.mat.R**:
    1. concatenate population mafs, and subset/intersect so that there is no missing info in any population for a given site. convert to counts (allele count matrix (dim: nrows=populations, ncols=loci) 
    2. run pairwise fst, to obtain matrix (within 01.scripts/Rscripts/mafs.ac.to.ac.pop.mat.R)
4. To obtain geographic/distance matrix and environmental/altitude matrix; run **local/scripts/dist.matrix.map.**R and follow least cost path tutorial **extras/tutorial-Realistic topographical distance between populations of Heliconius.pdf**

## 12. local PCA per SHDR

Local PCAs with all outlier windows (zPBS/zFst > 4, i.e. >4 standard deviations from the mean) of each SHDR (total = 370 local PCAs). Use sites and regions specified in [02.info/sites/](http://02.info/sites/) (obtained with local/scripts/12.local.pca.sites). All individuals from replicate transects were included, leading to local PCAs for all Western and Eastern SHDRs of each species that included Colombia and Ecuador samples. 

Obtain genotype likelihoods in beagle format (-doGlf 2) as input for PCAngsd, similarly to the population structure analysis. 

Use **01.scripts/12_local_PCA_per_SHDR.sh**

```
# list of shdr, e.g. 'shdr.east.001'
for i in $(less SHDR.id.list.txt)
do
 sbatch -J $i.pca -A JIGGINS-SL2-CPU ../../../01.scripts/12_local_PCA_per_SHDR.sh $i 
done


```

### Figures

Most intermediate files (outputs from steps above) are within local/data/, and the scripts to make figures, permutations, and analyses are within local/data/scripts.


