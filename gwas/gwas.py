import hail as hl
import pandas as pd


def convert_vcf_to_plink(local_vcf_path, output_prefix, threadnum):
    """
    Converts a VCF file to PLINK format using Hail.

    :param local_vcf_path: Path to the local VCF file
    :param output_prefix: Prefix for the output PLINK files
    """
    hl.init(spark_conf={'spark.executor.cores': str(threadnum),
                        'spark.driver.memory': '8g'})

    # Import VCF from the local file system
    mt = hl.import_vcf(local_vcf_path, reference_genome='GRCh37', force=True)
    mt = hl.variant_qc(mt)

    # Filter variants with call rate > 0.95
    mt = mt.filter_rows(mt.variant_qc.call_rate > 0.95)

    # Export to PLINK format
    hl.export_plink(mt, output_prefix)
    print(f"Plink files generated with prefix '{output_prefix}'")

def create_phenotype_file(vcf_file, chrom, pos, ref, alt, output_file):

    # read vcf
    mt = hl.import_vcf(vcf_file, reference_genome='GRCh37', force=True)
    print("sucess import")
    # set variant
    variant = mt.filter_rows(
        (mt.locus.contig == chrom) &
        (mt.locus.position == pos) &
        (mt.alleles[0] == ref) &
        (mt.alleles[1] == alt)
    )
    print("success create variant profile")
    # find variant
    if variant.count_rows() == 0:
        print("Variant not found in the VCF file.")
        return

    geno = variant.GT.collect()
    sample_ids = variant.s.collect()
    print("success collect variants")
    # 1|1: 2 (variant), 0|1, 1|0, 0|0: 1 (control)
    phenotype = []
    for gt in geno:
        if gt.is_hom_var():
            phenotype.append(2)
        else:
            phenotype.append(1)
    print("success create phenotypes")
    # create data frame
    pheno_data = pd.DataFrame({
        'FID': 0,
        'IID': sample_ids,
        'PHENO': phenotype
    })

    pheno_data.to_csv(output_file, sep='\t', index=False, header=['FID', 'IID', 'PHENO'])
    print(f"Phenotype file saved to {output_file}")


def run_gwas_hail(genotype_file, phenotype_file, col_name, maf_threshold, output_prefix):

    mt = hl.read_matrix_table(genotype_file)
    phenotypes = hl.import_table(phenotype_file, key='sample_id')
    mt = mt.annotate_cols(pheno=phenotypes[mt.s])
    mt = hl.variant_qc(mt)
    mt = mt.filter_rows(mt.variant_qc.AF[1] > maf_threshold)
    gwas = hl.linear_regression_rows(
        y=mt.pheno[col_name],
        x=mt.GT.n_alt_alleles(),
        covariates=[1.0],
        pass_through=['rsid']
    )
    gwas_result_path = f'{output_prefix}.gwas.ht'
    gwas.write(gwas_result_path, overwrite=True)
    gwas_result = gwas.rows()
    gwas_result.export(f'{output_prefix}.gwas.txt')
    hl.stop()
    print("success create gwas result")