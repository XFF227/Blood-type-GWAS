#!/bin/bash

mkdir 'mytest'
root_dir="mytest"
# Local VCF file path and output prefix
local_vcf_path="ALL.chr9.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz"
output_prefix="mytest/test_chr9"

# Set the path to the VCF file, output phenotype file, and variant information
vcf_file="ALL.chr9.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz"
output_file="${root_dir}/test_chr9.txt"
chrom="9"
pos=136132908
ref="T"
alt="TC"

genotypeFile="${root_dir}/test_chr9" # the clean dataset we generated in previous section
phenotypeFile="${root_dir}/test_chr9.txt" # the phenotype file

colName="PHENO"
threadnum=20
maf_threshold=0.01
# Run the Python script with the specified function and arguments
python -c "
from gwas import convert_vcf_to_plink, create_phenotype_file, run_gwas_hail
mt = convert_vcf_to_plink('$local_vcf_path', '$output_prefix', '$threadnum')
create_phenotype_file(mt, '$chrom', '$pos', '$ref', '$alt', '$output_file')
run_gwas_hail('$genotypeFile', '$phenotypeFile', '$colName', '$maf_threshold', '$output_prefix')
"