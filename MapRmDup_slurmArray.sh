#!/bin/bash

### Slurm array script to run mapping and parsing

#SBATCH --job-name=human                       # Job name
#SBATCH --output=/local/scratch/fdumetz/human/HsHisat.out   # Standard output log
#SBATCH --error=/local/scratch/fdumetz/human/HsHisat.err    # Standard error log
#SBATCH --mail-type=BEGIN,END --mail-user=fdumetz@som.umaryland.edu
#SBATCH --cpus-per-task=12				# Number of CPUs per task
#SBATCH --mem=60G                                        # Memory per node
#SBATCH --array=1-49                     # Number of array jobs to run according to list size

# Paths and settings
HISAT2_INDEX="/local/projects-t3/SerreDLab-3/fdumetz/Pv_ethiopia/human/Hs_hisat_index/HomoS_hisat_ref"  # HISAT2 index
PARENT_DIR="/local/projects-t4/aberdeen2ro/SerreDLab-4/raw_reads/2024-12-10_Eugenia_Lo"  # Parent directory containing FASTQ files
ALIGN_OUTPUT_DIR="/local/scratch/fdumetz/human"  # Directory for aligned BAM files
SUM_DIR="$ALIGN_OUTPUT_DIR/sum_dir"
THREADS=32  # Number of threads for HISAT2

# Picard tool path (adjust if needed)
PICARD="/usr/local/packages/picard-tools/picard.jar"

# Create necessary directories
mkdir -p "$ALIGN_OUTPUT_DIR"
mkdir -p "$SUM_DIR"

# Get the subfolder name from the corresponding line in "slurm_dirList4.txt"
SUBFOLDER=$(sed "${SLURM_ARRAY_TASK_ID}q;d" /local/projects-t3/SerreDLab-3/fdumetz/Pv_ethiopia/human/slurmList.txt)

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
                                            -S "$ALIGN_OUTPUT_DIR/${BASENAME}_Hs.sam" \
                                            --max-intronlen 5000 \
                                            -p "$THREADS" \
                                            --summary-file "$SUM_DIR/${BASENAME}_Hs_summary.txt" \
                                            --new-summary 

    echo "Mapping completed for $BASENAME. Converting SAM to BAM..."

    # Convert SAM to BAM and sort
    samtools view -bhF 2308 "$ALIGN_OUTPUT_DIR/${BASENAME}_Hs.sam" > "$ALIGN_OUTPUT_DIR/${BASENAME}_Hs.bam"
    samtools sort "$ALIGN_OUTPUT_DIR/${BASENAME}_Hs.bam" -o "$ALIGN_OUTPUT_DIR/${BASENAME}_Hssorted.bam"

    # Step 2: Remove duplicates using Picard
    echo "Running Picard to remove duplicates..."
    java -jar "$PICARD" MarkDuplicates \
        INPUT="$ALIGN_OUTPUT_DIR/${BASENAME}_Hssorted.bam" \
        OUTPUT="$ALIGN_OUTPUT_DIR/${BASENAME}_Hs_dedup.bam" \
        METRICS_FILE="$ALIGN_OUTPUT_DIR/${BASENAME}_dedup_metrics.txt" \
        REMOVE_DUPLICATES=true \
        ASSUME_SORTED=true \
        VALIDATION_STRINGENCY=SILENT \
        TMP_DIR=/tmp

    samtools index "$ALIGN_OUTPUT_DIR/${BASENAME}_Hs_dedup.bam"

    # Cleanup: Remove intermediate files
    rm "$ALIGN_OUTPUT_DIR/${BASENAME}_Hs.sam"
    rm "$ALIGN_OUTPUT_DIR/${BASENAME}_Hs.bam"
    rm "$ALIGN_OUTPUT_DIR/${BASENAME}_Hssorted.bam"

    echo "Processing completed for $BASENAME. Deduplicated BAM saved as ${BASENAME}_Hs_dedup.bam."
else
    echo "Folder $SUBFOLDER does not exist."
    exit 1
fi

echo "All processes completed. Aligned and deduplicated files saved in $ALIGN_OUTPUT_DIR."
