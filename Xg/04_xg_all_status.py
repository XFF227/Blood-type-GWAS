# python script to create text file with FID IID PAT MAT Xg status

import pandas as pd

# Define the file paths
psam_file = "all_phase3_xg.psam"
txt_file = "xg_genotype_info_modified_trimmed.txt"
output_file = "all_phase3_xg.txt"

## Read the PSAM file
psam_df = pd.read_csv(psam_file, sep="\t")

# Read the TXT file
txt_df = pd.read_csv(txt_file, sep="\t")

# Merge the dataframes on the IID column and select the first 3 columns
merged_df = pd.merge(psam_df, txt_df[["FID", "IID", "xg"]], left_on="#IID", right_on="IID")

# Save the merged dataframe to a file
merged_df.to_csv(output_file, sep="\t", index=False)