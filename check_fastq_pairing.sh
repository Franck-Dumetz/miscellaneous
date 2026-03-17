#!/bin/bash
#SBATCH --job-name=check_fastq_pair
#SBATCH --output=check_fastq_pair_%A.out
#SBATCH --error=check_fastq_pair_%A.err
#SBATCH --account=dserre-lab
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G

set -euo pipefail

if [[ $# -ne 3 ]]; then
    echo "Usage: sbatch $0 <R1.fastq.gz> <R2.fastq.gz> <result.txt>"
    exit 1
fi

R1="$1"
R2="$2"
RESULT="$3"

{
    echo "===================================="
    echo "FASTQ pairing check"
    echo "Date: $(date)"
    echo "Host: $(hostname)"
    echo "R1: $R1"
    echo "R2: $R2"
    echo "===================================="
    echo

    echo "Checking files exist"
    for f in "$R1" "$R2"; do
        if [[ ! -f "$f" ]]; then
            echo "ERROR: File not found: $f"
            exit 1
        fi
    done
    echo "Both files found."
    echo

    echo "Testing gzip integrity"
    gzip -t "$R1"
    echo "R1 gzip OK"
    gzip -t "$R2"
    echo "R2 gzip OK"
    echo

    echo "Counting reads"
    R1_LINES=$(zcat "$R1" | wc -l)
    R2_LINES=$(zcat "$R2" | wc -l)

    if (( R1_LINES % 4 != 0 )); then
        echo "ERROR: R1 line count is not divisible by 4: $R1_LINES"
        exit 1
    fi

    if (( R2_LINES % 4 != 0 )); then
        echo "ERROR: R2 line count is not divisible by 4: $R2_LINES"
        exit 1
    fi

    R1_READS=$((R1_LINES / 4))
    R2_READS=$((R2_LINES / 4))

    echo "R1 reads: $R1_READS"
    echo "R2 reads: $R2_READS"

    if [[ "$R1_READS" -eq "$R2_READS" ]]; then
        echo "Read counts match"
    else
        echo "WARNING: Read counts do not match"
    fi
    echo

    echo "Checking that read IDs match in order"
    paste \
      <(zcat "$R1" | awk 'NR%4==1{sub(/^@/,""); sub(/\/1$/,""); sub(/ 1:.*/,""); print}') \
      <(zcat "$R2" | awk 'NR%4==1{sub(/^@/,""); sub(/\/2$/,""); sub(/ 2:.*/,""); print}') \
    | awk '
      $1 != $2 {
        print "Mismatch at pair " NR
        print "R1: " $1
        print "R2: " $2
        bad=1
        exit 1
      }
      END {
        if (!bad) print "All read IDs match and are in the same order"
      }
    '
    echo

    echo "Average read lengths"
    zcat "$R1" | awk 'NR%4==2{sum+=length($0); n++} END{if(n>0) print "R1 avg length:", sum/n; else print "R1 avg length: NA"}'
    zcat "$R2" | awk 'NR%4==2{sum+=length($0); n++} END{if(n>0) print "R2 avg length:", sum/n; else print "R2 avg length: NA"}'
    echo

    echo "DONE"
} > "$RESULT" 2>&1