#!/bin/bash

# unable to make the script work with the xg_genotype_info_modified.txt file, had to create the desired file
# manually by copying the first three columns of the xg_genotype_info_modified.txt file and the all_phase3_xg_modified.psam file
# in google docs

# Define the file paths
PSAM_FILE="all_phase3_xg_modified.psam"
GENOTYPE_FILE="xg_genotype_info_modified.txt"
TRIMMED_GENOTYPE_FILE="xg_genotype_info_modified_trimmed.txt"
OUTPUT_FILE="merged_output.tsv"

# Step 1: Trim the xg_genotype_info_modified.txt file to keep only the first three columns
awk '{print $1, $2, $3}' "$GENOTYPE_FILE" > "$TRIMMED_GENOTYPE_FILE"

# Step 2: Merge the trimmed file with the PSAM file on the IID column
awk '
BEGIN {
    FS = OFS = "\t"
}

# Read the trimmed genotype file and store the values in an associative array
FNR == NR {
    if (NR > 1) {  # Skip header
        genotype[$2] = $3
    }
    next
}

# Process the PSAM file
{
    if (FNR == 1) {
        # Print the header with additional xg column
        print $0, "xg"
    } else {
        iid = $2
        if (iid in genotype) {
            print $0, genotype[iid]
        } else {
            print $0, ""
        }
    }
}
' "$TRIMMED_GENOTYPE_FILE" "$PSAM_FILE" > "$OUTPUT_FILE"

echo "Merged file created: $OUTPUT_FILE"