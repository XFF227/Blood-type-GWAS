#!/bin/bash

# Define the directory and patterns
directory="$(pwd)/all1kg/plink_results"
start_pattern="xg_gwas"
end_pattern="*.glm.firth"
exclude_pattern="*.glm.firth.id"

# Find all files matching the start and end patterns, excluding the exclude pattern
file_list=$(find "$directory" -type f -name "${start_pattern}*${end_pattern}" ! -name "$exclude_pattern")

# Initialize a counter for the total number of lines
total_lines=0

# Loop through each file and count the number of lines
for file in $file_list; do
  line_count=$(wc -l < "$file")
  total_lines=$((total_lines + line_count))
done

# Print the total number of lines
echo "Total number of lines: $total_lines"