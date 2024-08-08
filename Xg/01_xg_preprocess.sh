#!/bin/bash

# 01_xg_process.sh
# convert vcf of chromosome X into other formats

# NUS HPC load modules
module load plink
module load bcftools

# convert to PED/MAP format. MAP file contain all markers, PED file has sample information and genotype columns like 
# FID     IID     PID MID Sex (1 male, 2 female, 0 unknown) Phenotype (-9 unknown) Genotype
# FAM001  IND001  0  0  1  -9  A A  C C  T T
# FAM001  IND002  0  0  2  -9  A G  C C  T T
:'plink \
    --vcf ALL.chrX.phase3_shapeit2_mvncall_integrated_v1c.20130502.genotypes.vcf.gz \
    --recode --out fid_iid_sex
'

# convert to BED/BIM/FAM format after filerting 
echo -e "\n convert vcf file to bed/bim/fam, filter out variants with missing rate > 0.01, not in Hardy-Weinberg equilibrium, minor allele frequency < 0.01 \n"
: 'plink \
    --vcf All_1KG_combined.vcf.gz \
    --geno 0.01 \
    --hwe 1e-6 \
    --maf 0.01 \
    --mind 0.01 \
    --indep-pairwise 50 5 0.2 \
    --make-bed \
    --out \plink_results\All_plink_results \ 
'


# Extract the genotype information for the specific variant
# echo -e "\n  Extract the genotype information for the specific variant \n"
# bcftools query -r X:2666384 -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t[%SAMPLE\t%GT\n]' ALL.chrX.phase3_shapeit2_mvncall_integrated_v1c.20130502.genotypes.vcf.gz > xg_genotype_info.txt


# split PED file into chunks
split -l 1000 output.ped ped_chunk_


# encode Xg phenotype in phenotype column of PED file 
# Recode the phenotype column based on the extracted genotype information. G is reference, C is alternate
echo -e "\n parallelized encode Xg phenotype in phenotype column of PED file  \n"
ls ped_chunk_* | xargs -P 4 -I {} sh -c '
    awk -v geno="xg_genotype_info.txt" '
    BEGIN {FS=OFS="\t"}
    NR==FNR {
        if ($6 == "0|0" || $6 == "0/0") genotype[$1] = "G G";
        else if ($6 == "0|1" || $6 == "1|0" || $6 == "0/1" || $6 == "1/0") genotype[$1] = "G C";
        else if ($6 == "1|1" || $6 == "1/1") genotype[$1] = "C C";
        next
    }
    {
        if (FNR == 1) print $0;
        else {
            if (genotype[$2] == "C C") $6 = 1;
            else if (genotype[$2] == "G C" || genotype[$2] == "G G") $6 = 2;
            print
        }
    }
    ' xg_genotype_info.txt {} > {}.recoded
'
# combine recoded files
cat ped_chunk_*.recoded > parallelized_recoded.ped
