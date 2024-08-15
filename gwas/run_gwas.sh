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
pos=132136935
ref="G"
alt="A"

# Run the Python script with the specified function and arguments
python -c "
from gwas import convert_vcf_to_plink, create_phenotype_file
convert_vcf_to_plink('$local_vcf_path', '$output_prefix')
create_phenotype_file('$vcf_file', '$chrom', '$pos', '$ref', '$alt', '$output_file')

mt = hl.read_matrix_table(genotype_file)
phenotypes = hl.import_table(phenotype_file, key='sample_id')
mt = mt.annotate_cols(pheno=phenotypes[mt.s])
pheno_name = col_name
mt = hl.variant_qc(mt)
mt = mt.filter_rows(mt.variant_qc.AF[1] > 0.01)
gwas = hl.linear_regression_rows(
    y=mt.pheno[pheno_name],
    x=mt.GT.n_alt_alleles(),
    covariates=[1.0],
    pass_through=['rsid']
)
gwas.write('1kgeas.gwas.ht', overwrite=True)
gwas_result = gwas.rows()
gwas_result.export('1kgeas.gwas.txt')
"

genotypeFile="${root_dir}/test_chr9" # the clean dataset we generated in previous section
phenotypeFile="${root_dir}/test_chr9.txt" # the phenotype file

colName="PHENO"
threadnum=20

plink2 \
    --bfile ${genotypeFile} \
    --pheno ${phenotypeFile} \
    --pheno-name ${colName} \
    --maf 0.01 \
        --glm allow-no-covars \
    --threads ${threadnum} \
    --out 1kgeas
