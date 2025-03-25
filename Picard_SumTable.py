#!/usr/bin/env python

import os
import glob

# Define input directory and output file
input_dir = input("Enter the directory containing _metrics.txt files: ").strip()
output_file = input("Enter the name of the output file: ").strip()

# Prepare the output file
with open(output_file, "w") as out:
    out.write("Sample\tPercent Duplication\tEstimated Library Size\n")
    
    # Find all metric files
    metric_files = glob.glob(os.path.join(input_dir, "*_metrics.txt"))
    
    for file in metric_files:
        sample_name = os.path.basename(file).replace("_metrics.txt", "")
        percent_duplication = "NA"
        estimated_library_size = "NA"
        
        with open(file, "r") as f:
            for line in f:
                if line.startswith("Unknown Library"):
                    parts = line.strip().split("\t")
                    if len(parts) >= 9:
                        percent_duplication = parts[8]  # PERCENT_DUPLICATION
                    if len(parts) >= 10:
                        estimated_library_size = parts[9]  # ESTIMATED_LIBRARY_SIZE
                    break  # Stop reading after finding the relevant line
        
        # Write extracted values to output file
        out.write(f"{sample_name}\t{percent_duplication}\t{estimated_library_size}\n")

print(f"Summary saved to {output_file}")
