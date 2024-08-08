#!/bin/bash

# script to run mutliple scripts 

# concatenate all gwas results into one file for each super population
/hpctmp/xgao32/Blood-type-GWAS/Xg/05_concatenate_gwas_result.sh

# does not work? 
Rscript /hpctmp/xgao32/Blood-type-GWAS/Xg/02b_xg_plot_qq_manhattan.R