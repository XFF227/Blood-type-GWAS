# script to make Manhattan and QQ plots for Xg GWAS result for 1KG population data

# NUS HPC requires source /app1/ebenv
# do not use R-bundle-Bioconductor which is R 4.3.2 and cannot install MASS package 
# module load R/4.2.2-foss-2022b
# make sure current directory is Blood-type-GWAS. If in Xg directory, renv will fail.
# open R then use renv::restore() to create environment matching the same as in github codespace

# Load Libraries
library(qqman)
library(data.table)

# Define the base directory and super populations
base_directory <- "Xg"
super_pops <- c("AFR", "AMR", "EAS", "EUR", "SAS")

# Function to create plots
create_plots <- function(super_pop) {
    # Define file paths
    DIR_FILE <- file.path(base_directory, super_pop, "plink_results", paste0(super_pop, "_xg_gwas_combined_results.glm.firth.gz"))
    PATH <- file.path(getwd(), DIR_FILE)
    
    # Read GWAS data
    gwas_data <- fread(PATH)
    
    # Check Data Structure
    str(gwas_data)
    
    # Change X to 23, change data type to numeric
    gwas_data$`#CHROM`[gwas_data$`#CHROM` == "X"] <- 23
    gwas_data$`#CHROM` <- as.numeric(gwas_data$`#CHROM`)
    
    # Create Manhattan Plot
    png(file.path(base_directory, super_pop, paste0(super_pop,"_all_xg_chr_manhattan_plot.png")), width = 1600, height = 1000)
    par(cex.main=4, cex.lab=3, cex.axis=2.6, mar=c(6, 6, 6, 10))
    manhattan(gwas_data, chr="#CHROM", bp="POS", snp="ID", p="P", 
                main=paste("Manhattan Plot 1KG All Chromosomes -", super_pop),
                ylim=c(0, ceiling(max(-log10(gwas_data$P), na.rm = TRUE))),
                cex.main=4, cex.lab=3, cex.axis=2.6, cex.names=2.4, cex.axis.names=2.4)
    dev.off()
    
    # Calculate lambda for QQ plot
    calculate_lambda <- function(p_values) {
        chi_squared <- qchisq(1 - p_values, df = 1)
        lambda <- median(chi_squared, na.rm = TRUE) / qchisq(0.5, df = 1)
        return(lambda)
    }

    lambda <- calculate_lambda(gwas_data$P)
    max_logp <- max(-log10(gwas_data$P), na.rm = TRUE)

    # Create QQ Plot
    png(file.path(base_directory, super_pop, paste0(super_pop, "_xg_all_chr_1KG_qqplot_lambda.png")), width = 650, height = 600)
    par(cex.main=2, cex.lab=1.5, cex.axis=1.3, mar=c(6, 6, 6, 10))
    qq(gwas_data$P, main=paste("QQ Plot \n 1KG Xg All Chr -", super_pop))
    text(1, max_logp - 0.3 * max_logp, paste("\u03BB =", round(lambda, 3)), adj=c(0,1), cex=2)
    dev.off()

    # Filter out rows where #CHROM is 23 for autosomes only plots
    filtered_gwas_data <- gwas_data[gwas_data$`#CHROM` != 23, ]

    # Create Manhattan Plot for autosomes
    png(file.path(base_directory, super_pop, paste0(super_pop,"_xg_autosomes_all_1KG_manhattan_plot.png")), width = 1600, height = 1000)
    par(cex.main=4, cex.lab=3, cex.axis=2.6, mar=c(6, 6, 6, 10))
    manhattan(filtered_gwas_data, chr="#CHROM", bp="POS", snp="ID", p="P", 
            main=paste("Manhattan Plot 1KG - Autosomes -", super_pop),
            ylim=c(0, ceiling(max(-log10(filtered_gwas_data$P), na.rm = TRUE))),
            cex.main=4, cex.lab=3, cex.axis=2.6, cex.names=2.4, cex.axis.names=2.4)
    dev.off()

    # Calculate lambda for autosomes QQ plot
    lambda <- calculate_lambda(filtered_gwas_data$P)
    max_logp <- max(-log10(filtered_gwas_data$P), na.rm = TRUE)

    # Create QQ Plot for autosomes
    png(file.path(base_directory, super_pop, paste0(super_pop,"_xg_autosomes_chr_1KG_qqplot_lambda_filtered.png")), width = 650, height = 600)
    par(cex.main=2, cex.lab=1.5, cex.axis=1.3, mar=c(6, 6, 6, 10))
    qq(filtered_gwas_data$P, main=paste("QQ Plot \n 1KG Xg Autosomes -", super_pop))
    text(1, max_logp - 0.3 * max_logp, paste("\u03BB =", round(lambda, 3)), adj=c(0,1), cex=2)
    dev.off()

    # List of significant SNPs
    threshold <- 5e-8
    subset_data <- gwas_data[gwas_data$P < threshold & gwas_data$`#CHROM` != 23, ]
    PATH_AUT_SNPS <- file.path(base_directory, super_pop, paste0(super_pop,"_xg_autosomes_no_chrX_all_1KG_sig_snps.txt"))
    write.table(subset_data, file = PATH_AUT_SNPS, sep = "\t", quote = FALSE, row.names = FALSE)

    subset_data_chrX <- gwas_data[gwas_data$P < threshold & gwas_data$`#CHROM` == 23, ]
    PATH_X_SNPS <- file.path(base_directory, super_pop, paste0(super_pop,"_xg_chrX_all_1KG_sig_snps.txt"))
    write.table(subset_data_chrX, file = PATH_X_SNPS, sep = "\t", quote = FALSE, row.names = FALSE)
    }

# Loop through each super population and create plots
for (super_pop in super_pops) {
    create_plots(super_pop)
    }