## Steering mutations pipeline (layer 14)

This repo builds steering vectors from `evo_gcp` embeddings and evaluates how steering changes codon likelihoods across all 64 variants for a given AID and scale.

### Key scripts

- `build_steering_vec.sh`: Orchestrates building steering vectors from embeddings.
  - Inputs: `input/updated_data.tsv`, embeddings under `jobs/<job>/output/*.npy`.
  - Args: positional `win_len` (default 200), optional margins (left/right; defaults 2000/1000).
  - Calls `scripts/make_steering_vectors.R` for `pre_norm` and `mlp_l3`.

- `scripts/make_steering_vectors.R`: Reads `.npy` embeddings, computes per‑AID difference vectors and the average steering vector.
  - Windowing: start index depends on layer (`pre_norm`: 1; `mlp_l3`: 2), length = `--win_len`.
  - Vector formula: per‑position differences (target − source) → column mean → L2‑normalize per‑AID → average across AIDs → L2‑normalize.
  - Outputs (non‑normalized per‑position diffs before averaging):
    - Subject vectors: `vectors/non_norm_diff_vec/non_norm_<aid>_<embed>_len<win_len>.tsv`
    - Steering vector: `vectors/non_norm_diff_vec/non_norm_steering_<embed>_len<win_len>.tsv`

- `scripts/run_64_codons_steering.sh`: Generates all 64 codon variants for an AID, submits `evo_gcp` summaries with a steering vector, downloads, merges, and plots.
  - Detects vector kind from path: normalized vs non‑normalized.
  - Job naming: `<aid>-<layer_kind>-scale<S>-vec-<norm|non-norm>-len<LEN>`
  - Outputs per run:
    - FASTA/TSV: `output/<normalized_steering|non_normalized_steering>/steer_scale_<S>/len_<LEN>/<aid>/`
    - Figure: `figures/<normalized_steering|non_normalized_steering>/steer_scale_<S>/len_<LEN>/<aid>/steered_vs_unsteered.pdf`

- `scripts/create_steering_table.r`: Merges `input_summary_unsteered.txt` with `input_summary_scale_<scale>.txt` (accepts `10` or `10.0`).
  - Output: `steered_vs_unsteered_<seq>_<pos>.tsv` with columns `seq_id, unsteered, steered`.

- `scripts/plot_steered_vs_unsteered.r`: Plots steered vs unsteered totals for 64 codons with labels `AA_CODON`; optionally highlights source/target codons.

### Vectors directory

- Non‑normalized outputs (default): `vectors/non_norm_diff_vec/`
  - Subject: `non_norm_<aid>_<embed>_len<win_len>.tsv`
  - Steering: `non_norm_steering_<embed>_len<win_len>.tsv`
- Normalized aliases (optional): `vectors/normalized_diff_vec/`
  - Steering: `norm_steering_<embed>.tsv` (if created)

### Figures directory

Plots are segregated to avoid overwrite across scale/vector kind/length:

- `figures/non_normalized_steering/steer_scale_<S>/len_<LEN>/<aid>/steered_vs_unsteered.pdf`
- `figures/normalized_steering/steer_scale_<S>/len_<LEN>/<aid>/steered_vs_unsteered.pdf`

### Typical runs

1) Build vectors from existing embeddings with a specific window length (e.g., 100):
```bash
bash build_steering_vec.sh 100
```

2) Run 64‑codon steering for an AID using a non‑normalized pre_norm vector at scale 20:
```bash
LEN=100
bash scripts/run_64_codons_steering.sh \
  --aid BAA \
  --steering_vec vectors/non_norm_diff_vec/non_norm_steering_pre_norm_len${LEN}.tsv \
  --scale 20 \
  --layer_kind pre_norm
```

3) Negative scales
```bash
# Either quote the minus or use the documented n-prefix (n100 == -100)
bash scripts/run_64_codons_steering.sh --aid BAA --steering_vec vectors/non_norm_diff_vec/non_norm_steering_pre_norm_len100.tsv --scale "-100" --layer_kind pre_norm
bash scripts/run_64_codons_steering.sh --aid BAA --steering_vec vectors/non_norm_diff_vec/non_norm_steering_pre_norm_len100.tsv --scale n100 --layer_kind pre_norm
```

### Notes

- Embeddings must exist under the job dirs referenced by `build_steering_vec.sh` (or adjust `--jobs_dir` when calling the R script directly).
- The `evo_gcp` CLI and its steering options are documented here: [evo2_gcp docs](https://github.com/eitanyaffe/evo2_gcp/tree/main).
- Ensure your `input/updated_data.tsv` has: `aid gene contig start end strand flipped src_codon tgt_codon mut_pos`.

