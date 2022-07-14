## Description

Data and scripts associated with the article "Repeated genetic adaptation to altitude in two tropical butterflies", accepted in Nature Communications (2022), by Gabriela Montejo-Kovacevich, Joana I. Meier, Caroline N. Bacquet, Ian A. Warren, Yingguang Frank Chan, Marek Kucka, Camilo Salazar, Nicol Rueda, Stephen H. Montgomery, W. Owen McMillan, Krzysztof M. Kozak, Nicola J. Nadeau, Simon Martin, Chris D. Jiggins.

Sequence data is deposited at ENA (primary accession: PRJEB35570	ERP118642) or elsewhere on ENA, as specified in Supplementary Table 1 (see manuscript). 

Extra metadata on each individual butterfly, including photos for most (file names of photos are indicated on 02.info/all.bams.list.info.csv), can be found on the [Earthcape database](https://heliconius.ecdb.io/) and FAQ on how to use it can be found [here](https://heliconius.zoo.cam.ac.uk/databases/earthcape-specimen-database/). 02.info/all.bams.list.info.csv has all information for the 518 individuals of *H. erato* and *H. melpomene* included in this study. Supplementary Table 1 also contains information for the other species included in this study for introgression analyses.

The scripts within 01.scripts/ are thought for a computing cluster. Local scripts and intermediate files can be found in local/. These scripts will require adapting to local requirements (cpus, memory, working directories) and a fair amount of work to understand them. The scripts to produce the figures are in local/scripts/.

Any requests please contact me (email on my github profile)

## Contents

Only the first two directory levels are shown, many more files within some of the directories.


```
├── 01.scripts
│ ├── 01_config.sh
│ ├── 03_saf_maf_gl_era_doGeno.sh
│ ├── 03_saf_maf_gl_mel_doGeno.sh
│ ├── 03_saf_maf_gl_sides_doGeno.sh
│ ├── 04_pca_era.sh
│ ├── 04_pca_mel.sh
│ ├── 04_pca_sides.sh
│ ├── 05_saf_maf_by_pop.sh
│ ├── 06_2dsfs_by_pop_pair.sh
│ ├── 07_pbs_three_pops.sh
│ ├── 08_thetas.subpop.sh
│ ├── 09_thetas.subpop.dxy.sh
│ ├── 11_doMaf.all.bedassle.sh
│ ├── 11_doMaf.pop.bedassle.sh
│ ├── 12_local_PCA_per_SHDR.sh
│ ├── README.md
│ └── Rscripts
├── 02.info
│ ├── all.bams.list.info.csv
│ ├── dist
│ ├── pop.depth
│ ├── pop.info.csv
│ ├── pop.list
│ ├── pop.short.info.csv
│ ├── relatives.info.csv
│ ├── side.short.info.csv
│ ├── sites
│ └── stats
├── local
│ ├── data
│ ├── plots
│ └── scripts
└── pipelines
  └── tutorial-Realistic\ topographical\ distance\ between\ populations\ of\ Heliconius.pdf

```
