#!/bin/bash

set -euo pipefail
echo "Starting steering pipeline"

# Entry point. Positional args: win_len, left_margin, right_margin; hard-code layer=14 here
win_len="${1:-200}"
left_margin="${2:-2000}"
right_margin="${3:-1000}"
layer="14"
echo "Starting steering pipeline"
echo "left_margin: $left_margin"
echo "right_margin: $right_margin"
echo "layer: $layer"
echo "win_len: $win_len"

# Global FASTA to collect all sequences (original + mutant per aid)
global_fasta="output/queries_layer${layer}.fasta"
query_table="output/queries_layer${layer}.query.tsv"

# Validate required input files used throughout pipeline
if [ ! -f "input/human_contigs_src.fasta" ]; then
  echo "Error: input/human_contigs_src.fasta not found"; exit 1
fi
if [ ! -f "input/updated_data.tsv" ]; then
  echo "Error: input/updated_data.tsv not found"; exit 1
fi

# Create top-level directories once
mkdir -p input output jobs figures

# Generate per-aid FASTA and .tab for original and target codons across all rows
python3 scripts/gen_steering_vars.py \
  --fasta input/human_contigs_src.fasta \
  --updated-data input/updated_data.tsv \
  --left-margin "$left_margin" \
  --right-margin "$right_margin" \
  --output-fasta "$global_fasta" \
  --query-table "$query_table"

## Submit embeddings jobs with separate IDs and downloads
# pre_norm (static ids)
job_id_pre="steering-embeddings-prenorm"
job_version_pre="layer${layer}"
job_dir_pre="jobs/${job_id_pre}-${job_version_pre}"

# # evo_gcp submit --job "$job_id_pre" \
# #   --output_type logits_and_embedding \
# #   --input_fasta "$(pwd)/$global_fasta" \
# #   --query_table "$(pwd)/$query_table" \
# #   --job_version "$job_version_pre" \
# #   --embedding_layers "blocks.${layer}.pre_norm" --wait

mkdir -p "$job_dir_pre"
evo_gcp download --job "$job_id_pre" \
  --job_version "$job_version_pre" \
  --jobs_dir "$(pwd)/jobs"

# mlp.l3 (static ids)
job_id_mlp="get-steering-embeddings-mlpl3"
job_version_mlp="all-aids-layer${layer}"
job_dir_mlp="jobs/${job_id_mlp}-${job_version_mlp}"

# evo_gcp submit --job "$job_id_mlp" \
#   --output_type logits_and_embedding \
#   --input_fasta "$(pwd)/$global_fasta" \
#   --query_table "$(pwd)/$query_table" \
#   --job_version "$job_version_mlp" \
#   --embedding_layers "blocks.${layer}.mlp.l3" --wait

mkdir -p "$job_dir_mlp"
evo_gcp download --job "$job_id_mlp" \
  --job_version "$job_version_mlp" \
  --jobs_dir "$(pwd)/jobs"
  
# Build subject vectors and average steering vector (pre_norm)
echo "Building subject vectors and steering vector from embeddings..."
Rscript scripts/make_steering_vectors.R \
  --jobs_dir "$job_dir_pre" \
  --updated_data input/updated_data.tsv \
  --layer "$layer" \
  --embed pre_norm \
  --vectors_dir vectors/non_norm_diff_vec \
  --win_len "$win_len"

# Also build for mlp_l3 (will skip if embeddings not present)
Rscript scripts/make_steering_vectors.R \
  --jobs_dir "$job_dir_mlp" \
  --updated_data input/updated_data.tsv \
  --layer "$layer" \
  --embed mlp_l3 \
  --vectors_dir vectors/non_norm_diff_vec \
  --win_len "$win_len"

