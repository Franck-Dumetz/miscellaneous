#!/usr/bin/env bash
#SBATCH --job-name=fast5_bc_ama
#SBATCH --output=/path//ONT_rawsig/%x_%j.out
#SBATCH --error=/path//ONT_rawsig/%x_%j.err
#SBATCH --mail-type=BEGIN,END --mail-user=
#SBATCH --cpus-per-task=16
#SBATCH --array=0-3
#SBATCH --mem=100G

set -euo pipefail

BASE_DIR=/path/to/raw/fast5/ONT_rawsig

# If glob matches nothing, expand to empty list rather than literal
shopt -s nullglob

# Collect all *_fast5 directories (not *_fast5_guppy)
fast5_dirs=( "$BASE_DIR"/*_fast5 )

if (( ${#fast5_dirs[@]} == 0 )); then
  echo "No *_fast5 directories found in $BASE_DIR" >&2
  exit 1
fi

echo "Found ${#fast5_dirs[@]} fast5 dirs:"
for idx in "${!fast5_dirs[@]}"; do
  echo "  [$idx] ${fast5_dirs[$idx]}"
done

task_id=${SLURM_ARRAY_TASK_ID:-0}

# Safety check
if (( task_id < 0 || task_id >= ${#fast5_dirs[@]} )); then
  echo "SLURM_ARRAY_TASK_ID=$task_id is out of range (0..$(( ${#fast5_dirs[@]}-1 )))" >&2
  exit 1
fi

fast5_dir="${fast5_dirs[$task_id]}"
sample_name="$(basename "$fast5_dir")"
output_dir="$BASE_DIR/${sample_name}_guppy"

mkdir -p "$output_dir"

echo "[$(date)] Starting guppy for sample: $sample_name"
echo "Input directory : $fast5_dir"
echo "Output directory: $output_dir"

/usr/local/packages/guppy-4.2.2/bin/guppy_basecaller \
  --input_path "$fast5_dir" \
  --save_path "$output_dir" \
  --config /path/to/models/rna_r9.4.1_70bps_hac.cfg \
  --num_callers "$SLURM_CPUS_PER_TASK" \
  --fast5_out

cd "$output_dir"
tar -czvf "${sample_name}_guppy_workspace_fast5.tar.gz" workspace/

echo "[$(date)] Done with sample: $sample_name"
