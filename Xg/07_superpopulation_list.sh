#!/bin/bash

# create list of individuals for each superpopulation in the 1KG data

# Define the input file
input_file="all_xg.tsv"

# Get the unique super populations
super_pops=$(awk -F'\t' 'NR > 1 {print $6}' "$input_file" | sort | uniq)

# Loop through each super population and create a corresponding file
for super_pop in $super_pops; do
    output_file="${super_pop}_individuals.txt"
    
    # Extract the rows corresponding to the current super population
    awk -F'\t' -v sp="$super_pop" 'NR > 1 && $6 == sp {print $1, $2}' "$input_file" > "$output_file"
    
    echo "Created file: $output_file"
done