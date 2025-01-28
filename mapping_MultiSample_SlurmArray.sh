#!/bin/bash

### Slurm array script to run mapping and parsing

#SBATCH --job-name=5sp                       # Job name
#SBATCH --output=/local/projects-t3/SerreDLab-3/fdumetz/Pv_ethiopia/MultiSpecies/5sp_aligned/5spHisat2.out   # Standard output log
#SBATCH --error=/local/projects-t3/SerreDLab-3/fdumetz/Pv_ethiopia/MultiSpecies/5sp_aligned/5spHisat2.err    # Standard error log
#SBATCH --mail-type=BEGIN,END --mail-user=fdumetz@som.umaryland.edu
#SBATCH --cpus-per-task=12				# Number of CPUs per task
#SBATCH --mem=40G                                        # Memory per node
#SBATCH --array=1-46                     # Number of array to run according to line 27


# Paths and settings
HISAT2_INDEX="path/hisat/index/folder"  # Replace with your HISAT2 index base name
PARENT_DIR="/path/to/generic/sequencing/folder"  # Parent directory containing subfolders with FASTQ files
ALIGN_OUTPUT_DIR="path/aligned"  # Directory for aligned BAM files
SUM_DIR="5path/sum_dir"
THREADS=16  # Number of threads for HISAT2

# Create necessary directories
mkdir -p "$ALIGN_OUTPUT_DIR"
mkdir -p "$SUM_DIR"

# Get the subfolder name from the corresponding line in "slurm_dirList3.txt"
SUBFOLDER=$(sed "${SLURM_ARRAY_TASK_ID}q;d" /local/projects-t3/SerreDLab-3/fdumetz/Pv_ethiopia/MultiSpecies/slurm_dirList4.txt)

if [ -d "$SUBFOLDER" ]; then
    echo "Processing folder: $SUBFOLDER"

    # Find paired FASTQ files
    FILE1=$(find "$SUBFOLDER" -type f -name "*_R1_trimmed.fastq.gz")
    FILE2=$(find "$SUBFOLDER" -type f -name "*_R2_trimmed.fastq.gz")

    # Check if both files exist
    if [[ -z "$FILE1" || -z "$FILE2" ]]; then
        echo "Skipping $SUBFOLDER: Missing paired files"
        exit 1
    fi

    # Extract the base folder name to use for output naming
    BASENAME=$(basename "$SUBFOLDER")

    echo "Processing $FILE1 and $FILE2..."

    # Step 1: Map the paired-end sequences using HISAT2
    /usr/local/packages/hisat2-2.2.1/hisat2 -x "$HISAT2_INDEX" \
                                            -1 "$FILE1" \
                                            -2 "$FILE2" \
                                            -S "$ALIGN_OUTPUT_DIR/${BASENAME}_5sp.sam" \
                                            --max-intronlen 5000 \
                                            -p "$THREADS" \
                                            --summary-file "$SUM_DIR/${BASENAME}_5sp_summary.txt" \
                                            --new-summary 

    echo "Mapping completed for $BASENAME. Converting SAM to BAM..."

    # Convert SAM to BAM and sort
    samtools view -bhF 2308 "$ALIGN_OUTPUT_DIR/${BASENAME}_5sp.sam" > "$ALIGN_OUTPUT_DIR/${BASENAME}_5sp.bam"
    samtools sort "$ALIGN_OUTPUT_DIR/${BASENAME}_5sp.bam" -o "$ALIGN_OUTPUT_DIR/${BASENAME}_5spsorted.bam"
    samtools index "$ALIGN_OUTPUT_DIR/${BASENAME}_5spsorted.bam"

    # Remove intermediate SAM and unsorted BAM files
    rm "$ALIGN_OUTPUT_DIR/${BASENAME}_5sp.sam" "$ALIGN_OUTPUT_DIR/${BASENAME}_5sp.bam"

    echo "Processing completed for $BASENAME."
else
    echo "Folder $SUBFOLDER does not exist."
    exit 1
fi

echo "All processes completed. Aligned files saved in $ALIGN_OUTPUT_DIR."
