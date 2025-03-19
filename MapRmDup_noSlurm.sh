#!/bin/bash

# Paths and settings
HISAT2_INDEX="/usr/local/scratch/fdumetz/Vandana/Pf3d7_v68_index/Pf3d7_v68"  # HISAT2 index
PARENT_DIR="/usr/local/scratch/fdumetz/Vandana"  # Parent directory containing FASTQ files
ALIGN_OUTPUT_DIR="$PARENT_DIR/Pf_bam"  # Directory for aligned BAM files
SUM_DIR="$ALIGN_OUTPUT_DIR/sum_dir"
THREADS=12  # Number of threads for HISAT2

# Picard tool path
PICARD="/usr/local/packages/picard-2.9.4/picard.jar"

# Create necessary directories
mkdir -p "$ALIGN_OUTPUT_DIR"
mkdir -p "$SUM_DIR"

# Get the sample list
SAMPLE_LIST=( $(ls $PARENT_DIR/*_1.fq.gz | xargs -n 1 basename | sed 's/_1.fq.gz//') )

for SAMPLE in "${SAMPLE_LIST[@]}"; do
    FILE1="$PARENT_DIR/${SAMPLE}_1.fq.gz"
    FILE2="$PARENT_DIR/${SAMPLE}_2.fq.gz"

    # Check if both files exist
    if [[ ! -f "$FILE1" || ! -f "$FILE2" ]]; then
        echo "ERROR: Missing FASTQ pair for $SAMPLE. Skipping."
        continue
    fi

    echo "Processing files: $FILE1 and $FILE2..."

    # Step 1: Map the paired-end sequences using HISAT2
    /usr/local/packages/hisat2-2.2.1/hisat2 -x "$HISAT2_INDEX" \
                                            -1 "$FILE1" \
                                            -2 "$FILE2" \
                                            -S "$ALIGN_OUTPUT_DIR/${SAMPLE}.sam" \
                                            --max-intronlen 5000 \
                                            -p "$THREADS" \
                                            --summary-file "$SUM_DIR/${SAMPLE}_summary.txt" \
                                            --new-summary 

    # Verify if SAM file was created successfully
    if [[ ! -s "$ALIGN_OUTPUT_DIR/${SAMPLE}.sam" ]]; then
        echo "ERROR: HISAT2 did not produce an output SAM file for $SAMPLE. Skipping."
        continue
    fi

    # Convert SAM to BAM and sort
    samtools view -bhF 2308 "$ALIGN_OUTPUT_DIR/${SAMPLE}.sam" > "$ALIGN_OUTPUT_DIR/${SAMPLE}.bam"
    samtools sort "$ALIGN_OUTPUT_DIR/${SAMPLE}.bam" -o "$ALIGN_OUTPUT_DIR/${SAMPLE}_sorted.bam"

    # Step 2: Remove duplicates using Picard
    java -jar "$PICARD" MarkDuplicates \
        INPUT="$ALIGN_OUTPUT_DIR/${SAMPLE}_sorted.bam" \
        OUTPUT="$ALIGN_OUTPUT_DIR/${SAMPLE}_dedup.bam" \
        METRICS_FILE="$ALIGN_OUTPUT_DIR/${SAMPLE}_dedup_metrics.txt" \
        REMOVE_DUPLICATES=true \
        ASSUME_SORTED=true \
        VALIDATION_STRINGENCY=SILENT \
        TMP_DIR=/tmp

    # Verify if deduplicated BAM file exists
    if [[ ! -s "$ALIGN_OUTPUT_DIR/${SAMPLE}_dedup.bam" ]]; then
        echo "ERROR: Picard failed to create deduplicated BAM file for $SAMPLE. Skipping."
        continue
    fi

    samtools index "$ALIGN_OUTPUT_DIR/${SAMPLE}_dedup.bam"

    # Cleanup: Remove intermediate files
    rm "$ALIGN_OUTPUT_DIR/${SAMPLE}.sam"
    rm "$ALIGN_OUTPUT_DIR/${SAMPLE}.bam"
    rm "$ALIGN_OUTPUT_DIR/${SAMPLE}_sorted.bam"

    echo "Processing completed for $SAMPLE. Deduplicated BAM saved as ${SAMPLE}_dedup.bam."
    echo "All processes completed for $SAMPLE."

done

echo "All samples processed successfully."
