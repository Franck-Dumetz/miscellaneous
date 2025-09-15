#!/bin/bash

### Slurm array script to run mapping and parsing

#SBATCH --job-name=Pv_aligned                       # Job name
#SBATCH --output=PATH/PvHisat.out   # Standard output log
#SBATCH --error=PATH/PvHisat.err    # Standard error log
#SBATCH --mail-type=BEGIN,END --mail-user=
#SBATCH --cpus-per-task=12				# Number of CPUs per task
#SBATCH --mem=10G                                        # Memory per node
#SBATCH --array=1-100                     # Number of array jobs to run according to list size

# Paths and settings
HISAT2_INDEX="PATH/P01_hisat_index/PvivaxP01_hisat_index"  # HISAT2 index
PARENT_DIR="PATH/to/reads.fastq"  # Parent directory containing FASTQ files
ALIGN_OUTPUT_DIR="PATH"  # Directory for aligned BAM files
SUM_DIR="$ALIGN_OUTPUT_DIR/sum_dir"
THREADS=32  # Number of threads for HISAT2

# Picard tool path (adjust if needed)
PICARD="/usr/local/packages/picard-2.9.4/picard.jar"

# Create necessary directories
mkdir -p "$ALIGN_OUTPUT_DIR"
mkdir -p "$SUM_DIR"

# Get the subfolder name from the corresponding line in "slurm_dirList4.txt"
SUBFOLDER=$(sed "${SLURM_ARRAY_TASK_ID}q;d" PATH/slurmlist_Pv.txt)

# Ensure the subfolder exists
if [ ! -d "$SUBFOLDER" ]; then
    echo "ERROR: Folder $SUBFOLDER does not exist. Job $SLURM_ARRAY_TASK_ID cannot proceed."
    exit 1
fi

echo "Processing folder: $SUBFOLDER"

# Find paired FASTQ files
FILE1=$(find "$SUBFOLDER" -type f -name "*_R1_trimmed.fastq.gz")
FILE2=$(find "$SUBFOLDER" -type f -name "*_R2_trimmed.fastq.gz")

# Check if both files exist
if [[ -z "$FILE1" ]]; then
    echo "ERROR: Missing file: *_R1_trimmed.fastq.gz in $SUBFOLDER. Job $SLURM_ARRAY_TASK_ID aborted."
    exit 1
fi

if [[ -z "$FILE2" ]]; then
    echo "ERROR: Missing file: *_R2_trimmed.fastq.gz in $SUBFOLDER. Job $SLURM_ARRAY_TASK_ID aborted."
    exit 1
fi

# Extract the base folder name to use for output naming
BASENAME=$(basename "$SUBFOLDER")

echo "Processing files: $FILE1 and $FILE2..."

# Step 1: Map the paired-end sequences using HISAT2
/usr/local/packages/hisat2-2.2.1/hisat2 -x "$HISAT2_INDEX" \
                                        -1 "$FILE1" \
                                        -2 "$FILE2" \
                                        -S "$ALIGN_OUTPUT_DIR/${BASENAME}_Pv.sam" \
                                        --max-intronlen 5000 \
                                        -p "$THREADS" \
                                        --summary-file "$SUM_DIR/${BASENAME}_Pv_summary.txt" \
                                        --new-summary 

# Verify if SAM file was created successfully
if [ ! -s "$ALIGN_OUTPUT_DIR/${BASENAME}_Pv.sam" ]; then
    echo "ERROR: HISAT2 did not produce an output SAM file for $BASENAME. Job aborted."
    exit 1
fi

echo "Mapping completed for $BASENAME. Converting SAM to BAM..."

# Convert SAM to BAM and sort
samtools view -bhF 2308 "$ALIGN_OUTPUT_DIR/${BASENAME}_Pv.sam" > "$ALIGN_OUTPUT_DIR/${BASENAME}_Pv.bam"

# Verify if BAM file was created
if [ ! -s "$ALIGN_OUTPUT_DIR/${BASENAME}_Pv.bam" ]; then
    echo "ERROR: Failed to create BAM file from SAM for $BASENAME. Job aborted."
    exit 1
fi

samtools sort "$ALIGN_OUTPUT_DIR/${BASENAME}_Pv.bam" -o "$ALIGN_OUTPUT_DIR/${BASENAME}_Pvsorted.bam"

# Verify if sorted BAM file exists
if [ ! -s "$ALIGN_OUTPUT_DIR/${BASENAME}_Pvsorted.bam" ]; then
    echo "ERROR: Sorting of BAM file failed for $BASENAME. Job aborted."
    exit 1
fi

# Step 2: Remove duplicates using Picard
echo "Running Picard to remove duplicates..."
java -jar "$PICARD" MarkDuplicates \
    INPUT="$ALIGN_OUTPUT_DIR/${BASENAME}_Pvsorted.bam" \
    OUTPUT="$ALIGN_OUTPUT_DIR/${BASENAME}_Pv_dedup.bam" \
    METRICS_FILE="$ALIGN_OUTPUT_DIR/${BASENAME}_dedup_metrics.txt" \
    REMOVE_DUPLICATES=true \
    ASSUME_SORTED=true \
    VALIDATION_STRINGENCY=SILENT \
    TMP_DIR=/tmp

# Verify if deduplicated BAM file exists
if [ ! -s "$ALIGN_OUTPUT_DIR/${BASENAME}_Pv_dedup.bam" ]; then
    echo "ERROR: Picard failed to create deduplicated BAM file for $BASENAME. Job aborted."
    exit 1
fi

samtools index "$ALIGN_OUTPUT_DIR/${BASENAME}_Pv_dedup.bam"

# Cleanup: Remove intermediate files
rm "$ALIGN_OUTPUT_DIR/${BASENAME}_Pv.sam"
rm "$ALIGN_OUTPUT_DIR/${BASENAME}_Pv.bam"
rm "$ALIGN_OUTPUT_DIR/${BASENAME}_Pvsorted.bam"

echo "Processing completed for $BASENAME. Deduplicated BAM saved as ${BASENAME}_Pv_dedup.bam."

echo "All processes completed. Aligned and deduplicated files saved in $ALIGN_OUTPUT_DIR."
