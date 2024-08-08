# Manhattan plot and QQ plot for GWAS result for chromosome X SNPs associated with Xg 
# https://sahirbhatnagar.com/manhattanly/articles/web_only/manhattanly_full.html

# NUS HPC requires source /app1/ebenv
# do not use R-bundle-Bioconductor which is R 4.3.2 and cannot install MASS package 
# module load R/4.2.2-foss-2022b
# make sure current directory is Blood-type-GWAS. If in Xg directory, renv will fail.
# open R then use renv::restore() to create environment matching the same as in github codespace

# install.packages(c("plotly", "qqman", "data.table") if renv::restore() doesn't work

# ---- Load Libraries ----
library(qqman)
library(data.table)

# ---- Read GWAS Results ----
CUR_DIR = getwd() # should be in Blood-type-GWAS/ directory
DIR_FILE = "Xg/all1kg/plink_results/all_xg_gwas_combined_results.glm.firth.gz"
PATH = file.path(CUR_DIR,DIR_FILE)
gwas_data <- fread(PATH)

# ---- Check Data Structure ----
str(gwas_data)

# ---- Create Manhattan Plot with qqman ----

# change X to 23, change data type to numeric
gwas_data$`#CHROM`[gwas_data$`#CHROM` == "X"] <- 23
gwas_data$`#CHROM` <- as.numeric(gwas_data$`#CHROM`)

# Open a PNG device
png("Xg/xg_all_chr_manhattan_plot.png", width = 1600, height = 1000)
par(cex.main=4, cex.lab=3, cex.axis=2.6, mar=c(6, 6, 6, 10))
xg_manplot <- manhattan(gwas_data, chr="#CHROM", bp="POS", snp="ID", p="P", 
                        main="Manhattan Plot 1KG",
                        ylim=c(0, ceiling(max(-log10(gwas_data$P), na.rm = TRUE))),
                        cex.main=4, cex.lab=3, cex.axis=2.6, cex.names=2.4, cex.axis.names=2.4)
# Close the device
dev.off()

# ---- QQ plot with genomic inflation factor
# function to compute genomic inflation factor lambda
calculate_lambda <- function(p_values) {
    # Convert p-values to chi-squared statistics
    chi_squared <- qchisq(1 - p_values, df = 1)
    # Calculate lambda
    lambda <- median(chi_squared, na.rm = TRUE) / qchisq(0.5, df = 1)
    return(lambda)
}
lambda <- calculate_lambda(gwas_data$P) # need to handle NA values

# Calculate the maximum -log10(p-value)
max_logp <- max(-log10(gwas_data$P), na.rm = TRUE)

png("Xg/xg_all_chr_1KG_qqplot_lambda.png", width = 650, height = 600)
# Set graphical parameters to increase font sizes and margins
par(cex.main=2, cex.lab=1.5, cex.axis=1.3, mar=c(6, 6, 6, 10))
xg_qqplot <- qq(gwas_data$P, main="QQ Plot \n 1KG Xg All Chr")
text_x <- 1
text_y <- max_logp - 0.3 * max_logp
# Add the text to the plot
text(text_x, text_y, paste("\u03BB =", round(lambda, 3)), adj=c(0,1), cex=2)
dev.off()


# ---- plots for autosomes only ----
# Filter out rows where #CHROM is 23
gwas_data$`#CHROM` <- as.numeric(gwas_data$`#CHROM`)
filtered_gwas_data <- gwas_data[gwas_data$`#CHROM` != 23, ]

# Open a PNG device for the Manhattan plot
png("Xg/xg_autosomes_all_1KG_manhattan_plot.png", width = 1600, height = 1000)
    par(cex.main=4, cex.lab=3, cex.axis=2.6, mar=c(6, 6, 6, 10))
xg_manplot <- manhattan(filtered_gwas_data, chr="#CHROM", bp="POS", snp="ID", p="P", 
                        main="Manhattan Plot 1KG",
                        ylim=c(0, ceiling(max(-log10(filtered_gwas_data$P), na.rm = TRUE))),
                        cex.main=4, cex.lab=3, cex.axis=2.6, cex.names=2.4, cex.axis.names=2.4)
dev.off()

# ---- QQ plot with genomic inflation factor
# Calculate the maximum -log10(p-value) for the filtered data
max_logp <- max(-log10(filtered_gwas_data$P), na.rm = TRUE)

# Calculate lambda for the filtered data
lambda <- calculate_lambda(filtered_gwas_data$P) # need to handle NA values

# Open a PNG device for the QQ plot
png("Xg/xg_autosomes_chr_1KG_qqplot_lambda_filtered.png", width = 650, height = 600)
# Set graphical parameters to increase font sizes and margins
par(cex.main=2, cex.lab=1.5, cex.axis=1.3, mar=c(6, 6, 6, 10))
xg_qqplot <- qq(filtered_gwas_data$P, main="QQ Plot \n 1KG Xg Autosomes")
# Calculate the coordinates for placing the text
text_x <- 1
text_y <- max_logp - 0.3 * max_logp
# Add the text to the plot
text(text_x, text_y, paste("\u03BB =", round(lambda, 3)), adj=c(0,1), cex=2)
dev.off()


# ---- list of significant SNPs ----
# not on chromosome X, Subset the gwas_data table
threshold <- 5e-8;
subset_data <- gwas_data[gwas_data$P < threshold & gwas_data$`#CHROM` != 23, ]

# specific paths for significant X chromosome and autosomes SNPs
PATH_X_SNPS = "Xg/all1kg/xg_chrX_all_1KG_sig_snps.txt"
PATH_AUT_SNPS = "Xg/all1kg/xg_autosomes_no_chrX_all_1KG_sig_snps.txt"

# save subset_data to file
write.table(subset_data, file = file.path(CUR_DIR,PATH_AUT_SNPS), sep = "\t", quote = FALSE, row.names = FALSE)

# subset gwas_data table to contain only rows with #CHROM == 23 and below threshold for column P
subset_data_chrX <- gwas_data[gwas_data$P < threshold & gwas_data$`#CHROM` == 23, ]

# save subset_data_chrX to file
write.table(subset_data_chrX, file = file.path(CUR_DIR,PATH_X_SNPS), sep = "\t", quote = FALSE, row.names = FALSE)
