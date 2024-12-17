#!/bin/bash

# Ensure all dependencies are installed
# Usage: bash summarize_bams.sh <folder_with_bams> <output_summary.tsv>

# Input folder and output file
input_folder=$1
output_file=$2

# Check if required arguments are provided
if [ -z "$input_folder" ] || [ -z "$output_file" ]; then
    echo "Usage: bash $0 <folder_with_bams> <output_summary.tsv>"
    exit 1
fi

# Create the output file with headers
echo -e "File\tTotal_Reads\tMapped_Reads\t%_Mapped\t%_Duplicates" > "$output_file"

# Iterate over BAM files in the folder
for bam_file in "$input_folder"/*.bam; do
    echo "Processing $bam_file..."

    # Total reads
    total_reads=$(samtools view -c "$bam_file")

    # Mapped reads
    mapped_reads=$(samtools view -c -F 4 "$bam_file")

    # Percentage of mapping
    if [ "$total_reads" -gt 0 ]; then
        percentage_mapped=$(echo "scale=2; $mapped_reads / $total_reads * 100" | bc)
    else
        percentage_mapped=0
    fi

    # Duplicates using Picard
    metrics_file=$(mktemp)
    picard MarkDuplicates I="$bam_file" O=/dev/null M="$metrics_file" >/dev/null 2>&1
    duplicate_percentage=$(grep "PERCENT_DUPLICATION" -A 1 "$metrics_file" | tail -n 1 | awk '{print $9}')
    rm "$metrics_file"

    # Append results to the summary file
    echo -e "$(basename "$bam_file")\t$total_reads\t$mapped_reads\t$percentage_mapped\t$duplicate_percentage" >> "$output_file"
done

echo "Summary written to $output_file"
