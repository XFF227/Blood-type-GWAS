# Blood-type-GWAS
## Objective
### 1. identify tag SNPs in linkage disequilibrium with causative variants that result in the different blood groups (1k-10k per blood group) 
  1. GWAS blood group as phenotype, find non-overlapping SNPs
  2. check GWAS catalog for tag SNPs
  3. create database for tag SNPs 
    1. decide on where to store database and VCF files (aliyun temporarily)
### 2. predict blood group based on 
  1. alternate methods : multiplex Taqman PCR, sequencing
  2. potential challenges with homologous genes that are almost identical like RHD and RHCE, require long read sequencing to distinguish reads between the two genes
  3. for individuals without SNP data for specific blood group, treat as having the default blood group as the reference genome hg19 

#### use `renv` to manage R packages for this repository
Activate specific version of R on NUS HPC. Make sure current directory is Blood-type-GWAS. If in Xg directory, `renv` will fail. Open R then use `renv::restore()` to install packages.

```sh
source /app1/ebenv
# do not use R-bundle-Bioconductor which is R 4.3.2 and cannot install MASS package 
module load R/4.2.2-foss-2022b
```

For new projects
```R
library(renv) 
renv::init() # initialize new renv project in directoy, only do once for new directory

# packages to be installed
packages <- c("dplyr", "ggplot2", "tidyr", "manhattanly",  
		"qqman", "plotly") 

renv::install(packages)
renv::snapshot() # creat renv.lock file to save package versions
renv::restore() # install all packages in renv.lock
```

