#!/bin/bash

# Output file
output_file="pdb_dataset"

# Clear the output file before starting
> "$output_file"

# Track whether the header has been written
header_written=false

for folder in refined-set/*; do
    if [[ -d "$folder" ]]; then
        for file in "$folder"/*pocket*; do  # Directly match 'pocket' in filename
            if [[ -f "$file" ]]; then
                if ! $header_written; then
                    # Extract headers from the first output and write them
                    python3 project.py <<< "$file" | head -n 1 >> "$output_file"
                    header_written=true
                fi
                # Append only data (excluding header)
                python3 project.py <<< "$file" | tail -n +2 >> "$output_file"
            fi
        done
    fi
done
