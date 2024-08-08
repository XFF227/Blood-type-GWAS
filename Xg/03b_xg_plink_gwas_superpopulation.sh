#!/bin/bash

# GWAS for each superpopulation in 1KG data from PLINK website
# no white space in the variable assignment else error 

phenotypeFile="all_xg.tsv"
colName="xg"

# Define the super populations
super_pops=("AFR" "AMR" "EAS" "EUR" "SAS")

for super_pop in "${super_pops[@]}"; do
    echo -e "\nSuper Population: $super_pop\n"
    
    out="/hpctmp/xgao32/Blood-type-GWAS/Xg/${super_pop}/plink_results"
    mkdir -p "${out}" # make directory if not existing

    # Initialize the combined output file
    combined_output="${out}/xg_gwas_combined_${super_pop}.txt"
    echo -e "CHR\tSNP\tBP\tA1\tTEST\tOBS_CT\tBETA\tSTAT\tP" > "$combined_output"

    for chr in {1..22}; do
        echo -e "\nChromosome $chr.\n"

        genotypeFile="/hpctmp/xgao32/1KG/phase3_grch37/filtered_vcf/ALL.chr$chr.filtered" 
        echo "Genotype file: $genotypeFile"
        
        # Filter the genotype data for the current super population
        plink2 \
            --bfile "${genotypeFile}" \
            --keep "${super_pop}_individuals.txt" \
            --pheno "${phenotypeFile}" \
            --pheno-name "${colName}" \
            --glm dominant hide-covar firth \
            --out "${out}/xg_gwas_chr${chr}_${super_pop}"

        # Append the results to the combined output file
        cat "${out}/xg_gwas_chr${chr}_${super_pop}.xg.glm.dom" | tail -n +2 >> "$combined_output"
    done

    # X chromosome
    genotypeFile="/hpctmp/xgao32/1KG/phase3_grch37/filtered_vcf/ALL.chrX.filtered" 
    echo "Genotype file: $genotypeFile"
    plink2 \
        --bfile "${genotypeFile}" \
        --keep "${super_pop}_individuals.txt" \
        --pheno "${phenotypeFile}" \
        --pheno-name "${colName}" \
        --glm hide-covar firth \
        --chr X \
        --out "${out}/xg_gwas_chrX_${super_pop}"

    # Append the X chromosome results to the combined output file
    cat "${out}/xg_gwas_chrX_${super_pop}.xg.glm.firth" | tail -n +2 >> "$combined_output"
done