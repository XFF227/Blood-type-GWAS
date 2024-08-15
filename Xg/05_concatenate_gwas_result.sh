#!/bin/bash

# Define the base directory and pattern
base_directory="/hpctmp/xgao32/Blood-type-GWAS/Xg"
pattern="*.glm.firth"
exclude_pattern="*.glm.firth.id"

# Define the super populations
super_pops=("AFR" "AMR" "EAS" "EUR" "SAS")

# Loop through each super population directory
for super_pop in "${super_pops[@]}"; do
  directory="${base_directory}/${super_pop}/plink_results"
  output_file="${directory}/${super_pop}_xg_gwas_combined_results.glm.firth.gz"

  # List all files matching the pattern and exclude those ending with .glm.firth.id
  file_list=$(find "$directory" -type f -name "$pattern" ! -name "$exclude_pattern")

  # Initialize a flag to track the first file
  first_file=true

  # Loop through each file and concatenate them
  for file in $file_list; do
    if $first_file; then
      # Print the header and content for the first file
      cat "$file" | gzip > "$output_file"
      first_file=false
    else
      # Skip the header and append the content for subsequent files
      tail -n +2 "$file" | gzip >> "$output_file"
    fi
  done

  echo "Concatenation complete for ${super_pop}. Output saved to $output_file."
done



# Define the directory and pattern
directory="/hpctmp/xgao32/Blood-type-GWAS/Xg/all1kg/plink_results"
pattern="*.glm.firth"
exclude_pattern="*.glm.firth.id"
output_file="${directory}/all_xg_gwas_combined_results.glm.firth.gz"

# List all files matching the pattern and exclude those ending with .glm.firth.id
file_list=$(find "$directory" -type f -name "$pattern" ! -name "$exclude_pattern")

# Initialize a flag to track the first file
first_file=true

# Loop through each file and concatenate them
for file in $file_list; do
  if $first_file; then
    # Print the header and content for the first file
    cat "$file" | gzip > "$output_file"
    first_file=false
  else
    # Skip the header and append the content for subsequent files
    tail -n +2 "$file" | gzip >> "$output_file"
  fi
done

echo "Concatenation complete. Output saved to $output_file."
