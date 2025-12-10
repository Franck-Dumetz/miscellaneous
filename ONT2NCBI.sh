#!/usr/bin/env bash
#SBATCH --job-name=fast5_bc_ama
#SBATCH --output=/path//ONT_rawsig/%x_%j.out
#SBATCH --error=/path//ONT_rawsig/%x_%j.err
#SBATCH --mail-type=BEGIN,END --mail-user=
#SBATCH --cpus-per-task=16
#SBATCH --array=0-3
#SBATCH --mem=100G

set -euo pipefail

BASE_DIR=/path/to/ONT_rawsig

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

# Use all CPUs without oversubscribing threads
CPU_THREADS_PER_CALLER=2
NUM_CALLERS=$(( SLURM_CPUS_PER_TASK / CPU_THREADS_PER_CALLER ))
if (( NUM_CALLERS < 1 )); then
  NUM_CALLERS=1
fi

echo "[$(date)] Starting guppy for sample: $sample_name"
echo "Input directory : $fast5_dir"
echo "Output directory: $output_dir"
echo "num_callers=${NUM_CALLERS}, cpu_threads_per_caller=${CPU_THREADS_PER_CALLER} (total threads=$((NUM_CALLERS * CPU_THREADS_PER_CALLER)))"

/usr/local/packages/guppy-4.2.2_cpu/bin/guppy_basecaller \
  --input_path "$fast5_dir" \
  --save_path "$output_dir" \
  --config /usr/local/packages/guppy-4.2.2_gpu/data/rna_r9.4.1_70bps_hac.cfg \
  --num_callers "$NUM_CALLERS" \
  --cpu_threads_per_caller "$CPU_THREADS_PER_CALLER" \
  --fast5_out

cd "$output_dir"

# Only tar if workspace exists and is non-empty
if [ -d workspace ] && [ "$(ls -A workspace | wc -l)" -gt 0 ]; then
  tar -czvf "${sample_name}_guppy_workspace_fast5.tar.gz" workspace/
else
  echo "Warning: workspace/ directory missing or empty for $sample_name" >&2
fi

echo "[$(date)] Done with sample: $sample_name"
