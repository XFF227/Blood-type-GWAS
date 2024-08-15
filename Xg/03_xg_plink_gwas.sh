#!/bin/bash
# GWAS on 1KG data from PLINK website, must use plink2 for handling pfam pgen files
# no white space in the variable assignment else error 

# phenotype file with FID IID xg sex population super population etc.
phenotypeFile="all_xg.tsv"
#"../xg_genotype_info_modified.txt"  # this file missing sex and population columns

# ---- Section: GWAS for autosomes 
colName="xg"
out="/hpctmp/xgao32/Blood-type-GWAS/Xg/all1kg/plink_results"
:'
for chr in {1..22}; do
    echo -e "\nChromosome $chr.\n"

    genotypeFile="/hpctmp/xgao32/1KG/phase3_grch37/filtered_vcf/ALL.chr$chr.filtered" 
    echo "Genotype file: $genotypeFile"
    plink2 \
        --bfile "${genotypeFile}" \
        --pheno "${phenotypeFile}" \
        --pheno-name "${colName}" \
        --glm dominant hide-covar firth \
        --out "${out}/xg_gwas_chr${chr}"
done
'

# X chromsome
# --covar-sex not needed if --chr flag specifies X
genotypeFile="/hpctmp/xgao32/1KG/phase3_grch37/filtered_vcf/ALL.chrX.filtered" 
echo "Genotype file: $genotypeFile"
plink2 \
    --bfile "${genotypeFile}" \
    --pheno "${phenotypeFile}" \
    --pheno-name "${colName}" \
    --glm hide-covar firth \
    --chr X \
    --out "${out}/xg_gwas_chrX"