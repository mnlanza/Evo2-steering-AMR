#!/bin/bash
set -euo pipefail

# Notes:
# - Layer is hard-coded to 14
# - updated_data.tsv is hard-coded to input/updated_data.tsv
# - codon table is hard-coded to input/codon_table

LAYER=14
STEER_VEC=""
SCALE="20"
AID_FILTER=""
LAYER_KIND="pre_norm"   # or mlp_l3

while [[ $# -gt 0 ]]; do
  case "$1" in
    --aid) AID_FILTER="$2"; shift 2;;
    --steering_vec) STEER_VEC="$2"; shift 2;;
    --scale) SCALE="$2"; shift 2;;
    --layer_kind) LAYER_KIND="$2"; shift 2;;
    *) echo "Unknown arg: $1"; exit 1;;
  esac
done

if [[ -z "$AID_FILTER" ]]; then
  echo "Usage: $0 --aid AID [--embed pre_norm|mlp_l3] [--steering_vec file] [--scale val]"; exit 1
fi

mkdir -p output jobs figures results

# Constants
UPDATED="input/updated_data.tsv"
CODON_TABLE="input/codon_table"
FASTA="input/human_contigs_src.fasta"

# Determine steering kind for figure organization based on vector path
STEER_KIND="unsteered"
if [[ -n "$STEER_VEC" ]]; then
  # Classify by path keywords. Prioritize explicit non_norm.
  if [[ "$STEER_VEC" == *"non_norm"* ]]; then
    STEER_KIND="non_normalized_steering"
  elif [[ "$STEER_VEC" == *"/normalized/"* || "$STEER_VEC" == *"vectors/normalized"* || "$STEER_VEC" == *"normalized_diff_vec"* || "$STEER_VEC" == *"/norm_steering"* ]]; then
    STEER_KIND="normalized_steering"
  else
    STEER_KIND="non_normalized_steering"
  fi
fi
# Vector kind slug for job naming
case "$STEER_KIND" in
  normalized_steering) VEC_KIND="norm";;
  non_normalized_steering) VEC_KIND="non-norm";;
  *) VEC_KIND="unsteered";;
esac

# Human-readable label for figure titles
case "$STEER_KIND" in
  normalized_steering) STEER_LABEL="normalized";;
  non_normalized_steering) STEER_LABEL="non-normalized";;
  *) STEER_LABEL="unsteered";;
esac

# Parse vector length from filename if present (expects ..._len<NUM>.tsv)
if [[ -n "$STEER_VEC" && "$STEER_VEC" =~ _len([0-9]+)\.tsv$ ]]; then
  VEC_LEN="${BASH_REMATCH[1]}"
else
  VEC_LEN="na"
fi

# Read updated_data.tsv and process matching AID rows (no subshell)
while IFS=$'\t' read -r AID GENE CONTIG START END STRAND FLIPPED SRC_CODON TGT_CODON MUT_POS; do
  [[ "${AID}" != "${AID_FILTER}" ]] && continue
  aid_lc=$(echo "$AID" | tr '[:upper:]' '[:lower:]')
  seq_id="$CONTIG"
  pos="$MUT_POS"

  # Place outputs under vector kind, then scale, then length, then aid to avoid overwrites across kinds/lengths
  # output/<normalized|non_normalized|unsteered>/steer_scale_<scale>/len_<len>/<aid>
  out_dir="output/${STEER_KIND}/steer_scale_${SCALE}/len_${VEC_LEN}/${aid_lc}"
  mkdir -p "$out_dir"
  fasta_out="${out_dir}/all64_${seq_id}_${pos}.fasta"
  codon_tab_out="${out_dir}/all64_${seq_id}_${pos}.tab"

  echo "[Gen] ${AID} ${seq_id} pos ${pos} → ${fasta_out}"
  python3 scripts/generate_codons_variants.py \
    --fasta "$FASTA" \
    --codon-table "$CODON_TABLE" \
    --aa-coord "$pos" \
    --seq-id "$seq_id" \
    --output-fasta "$fasta_out" \
    --output-codon-table "$codon_tab_out" \
    --left-margin 0 \
    --right-margin 0 \
    --gene-start "$START" \
    --gene-end "$END"

  job_id=$(echo "${aid_lc}" | tr '_' '-')
  # Sanitize scale for job naming: replace '.' with 'p' to satisfy GCP job id rules
  SCALE_SLUG="${SCALE//./p}"
  job_version=$(echo "${LAYER_KIND//_/-}-scale${SCALE_SLUG}-vec-${VEC_KIND}-len${VEC_LEN}" | tr '[:upper:]' '[:lower:]')
  job_dir="jobs/${job_id}-${job_version}"
  mkdir -p "$job_dir"

  echo "[Submit] un/steered summaries for ${AID} scale ${SCALE}"
  # Map layer kind to evo_gcp layer name
  case "$LAYER_KIND" in
    pre_norm) layer_name_cli="pre_norm";;
    mlp_l3) layer_name_cli="mlp.l3";;
    *) echo "Unknown --layer_kind: $LAYER_KIND (expected pre_norm or mlp_l3)"; exit 1;;
  esac
  evo_gcp submit --job "$job_id" \
    --output_type summary_only \
    --input_fasta "$(pwd)/$fasta_out" \
    --job_version "$job_version" \
    --steering_layer "blocks.${LAYER}.${layer_name_cli}" \
    ${STEER_VEC:+--steering_vector_file "$(pwd)/$STEER_VEC"} \
    ${STEER_VEC:+--steering_scales "$SCALE"} \
    --wait

  echo "[Download] ${job_dir}"
  evo_gcp download --job "$job_id" \
    --job_version "$job_version" \
    --jobs_dir "$(pwd)/jobs"

  # Write comparison table beside the FASTA for this aid
  tbl_out="${out_dir}/steered_vs_unsteered_${seq_id}_${pos}.tsv"
  Rscript scripts/create_steering_table.r \
    --jobs_dir "$job_dir" \
    --scale "$SCALE" \
    --out "$tbl_out"

  # Save plots under figures/<normalized|non_normalized|unsteered>/steer_scale_<scale>/<aid>
  plot_dir="figures/${STEER_KIND}/steer_scale_${SCALE}/len_${VEC_LEN}/${aid_lc}"
  plot_args=(
    --table "$tbl_out"
    --title "${AID} ${seq_id}:${pos} (scale ${SCALE}, ${STEER_LABEL}, len ${VEC_LEN})"
    --outdir "$plot_dir"
    --codon_table "$CODON_TABLE"
  )
  if [ -n "${SRC_CODON:-}" ]; then plot_args+=( --src_codon "$SRC_CODON" ); fi
  if [ -n "${TGT_CODON:-}" ]; then plot_args+=( --tgt_codon "$TGT_CODON" ); fi
  Rscript scripts/plot_steered_vs_unsteered.r "${plot_args[@]}"

done < <(tail -n +2 "$UPDATED")



