# awk script to make text file of phenotype compatible with PLINK
# first column family ID (FID)
# second column individual ID (IID)
# third column phenotype (for dominant trait, 1 for case, 2 for control)

awk 'BEGIN {
    OFS="\t"
    print "FID", "IID", "xg"
}
{
    # Duplicate the first column
    $5 = $1
    
    # Create the xg column (1 if $2 is "1|1", 2 otherwise)
    $6 = ($2 == "1|1") ? 1 : 2
    
    # Print the columns in the desired order
    print $1, $5, $6, $2
}' xg_genotype_info.txt > xg_genotype_info_modified.txt
