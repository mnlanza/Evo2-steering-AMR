## Mutation‑aware steering with EVO2 (layer 14)

Evaluate whether mutation‑aware steering vectors can increase the model’s attention to biologically meaningful codon changes. This repo builds steering vectors from EVO2 embeddings and tests their effect on total sequence log‑likelihood across all 64 codon variants at a target amino‑acid position.

---

## At a glance

- **Goal**: Build per‑subject/gene(AID identifier) difference vectors around a mutation, average them into a steering vector, and test if steering shifts likelihoods in favor of target mutations.
- **Scope**: Uses EVO2 embeddings at layer 14; supports `pre_norm` and `mlp_l3` model embeddings; runs 64‑codon sweeps per subject across scales.
- **Outputs**: Subject vectors, aggregate steering vectors, subject 64‑codon summaries, and plots comparing steered vs unsteered totals.

---

## Results at a glance

- Working case: `figures/non_normalized_steering/steer_scale_-8/len_100/baa/steered_vs_unsteered.pdf`
- Not yet generalizable for other subject/genes: `figures/non_normalized_steering/steer_scale_-10/len_100/baf/steered_vs_unsteered.pdf`

These examples show steering can shift likelihoods in some subjects but is not yet robust across genes.

---

## Quickstart

### Prereqs

- Python 3, R
- EVO2 GCP CLI configured (see `evo2_gcp` docs)
- Python deps: `pip install -r requirements.txt`
- R deps (in R):
  - `install.packages(c("argparse", "ggplot2", "ggrepel"))`
  - `install.packages("reticulate")`

### Minimal demo

1) Build vectors from existing embeddings (choose window length, e.g. 100):
```bash
bash build_steering_vec.sh 100
```

2) Run 64‑codon steering for one AID using a non‑normalized `pre_norm` vector at scale -8:
```bash
LEN=100
bash scripts/run_64_codons_steering.sh \
  --aid BAA \
  --steering_vec vectors/non_norm_diff_vec/non_norm_steering_pre_norm_len${LEN}.tsv \
  --scale -8 \
  --layer_kind pre_norm
```

3) Try negative scales (inverts steering direction):
```bash
bash scripts/run_64_codons_steering.sh \
  --aid BAA \
  --steering_vec vectors/non_norm_diff_vec/non_norm_steering_pre_norm_len100.tsv \
  --scale -100 \
  --layer_kind pre_norm
```

Runtime varies with inputs and GCP throughput.

---

## Recommended defaults

- Layer kind: `pre_norm`
- Layer index: 14 (middle layer)
- Steering window length (`win_len`): 100
- Steering scale: ~ -8 to -10 (tune per dataset)

---

## Method choices and rationale

- Layer choice: 14 as a middle layer to balance locality with broader context.
- Embedding stream: `pre_norm` reflects activations just before logits, better capturing output‑proximal changes.
- Normalization: compute per‑position diffs (target − source) without per‑row normalization so mutation‑proximal positions retain higher weight; average across positions; L2‑normalize the per‑subject vector; average across subjects; L2‑normalize to form the final steering vector so each subject contributes equally.

---

## How it works

1. Variant batching (`scripts/gen_steering_vars.py`)
   - Reads contigs (`input/human_contigs_src.fasta`) and rows from `input/updated_data.tsv`.
   - Builds per‑row original/target codon sequences with user margins (defaults 2000/1000 nt).
   - Writes a global FASTA and a query table targeting a 201‑nt window centered at the mutation nucleotide for EVO2.

2. Embeddings (via `evo_gcp`)
   - `build_steering_vec.sh` sets `layer=14`, creates job scaffolding, and downloads precomputed outputs for `pre_norm` and `mlp_l3` job IDs.
   - You can switch to submitting jobs by uncommenting the provided `evo_gcp submit` commands.

3. Steering vectors (`scripts/make_steering_vectors.R`)
   - Loads `.npy` arrays with `reticulate` and selects a position window: start at 1 for `pre_norm`, at 2 for `mlp_l3`; length = `--win_len`.
   - Computes per‑position diffs (target − source) → column mean → L2‑normalize per‑AID → average across AIDs → L2‑normalize.
   - Writes subject vectors and an aggregate non‑normalized steering vector under `vectors/non_norm_diff_vec`.

4. 64‑codon steering runs (`scripts/run_64_codons_steering.sh`)
   - Generates all 64 codon variants for one AID with margins.
   - Submits `evo_gcp` summary jobs with optional `--steering_vector_file` and `--steering_scales`.
   - Downloads results and merges un/steered totals with `scripts/create_steering_table.r`.
   - Plots steered vs unsteered totals per codon with `scripts/plot_steered_vs_unsteered.r`.

---

## Inputs

1. `input/human_contigs_src.fasta`: source contigs (wild‑type)
2. `input/codon_table`: codon→aa map (table 11)
3. `input/updated_data.tsv`: columns include `aid, gene, contig, start, end, strand, flipped, src_codon, tgt_codon, mut_pos`

## Outputs

- Vectors (`vectors/non_norm_diff_vec/`):
  - Subject: `non_norm_<aid>_<embed>_len<win_len>.tsv`
  - Steering: `non_norm_steering_<embed>_len<win_len>.tsv`
- Per‑run artifacts:
  - FASTA/TSV: `output/<normalized_steering|non_normalized_steering|unsteered>/steer_scale_<S>/len_<LEN>/<aid>/`
  - Table: `steered_vs_unsteered_<seq>_<pos>.tsv` with `seq_id, unsteered, steered`
  - Figure: `figures/<normalized_steering|non_normalized_steering|unsteered>/steer_scale_<S>/len_<LEN>/<aid>/steered_vs_unsteered.pdf`

## Repo structure

```
build_steering_vec.sh           # end‑to‑end: batch FASTA → evo_gcp download → vectors
scripts/
  gen_steering_vars.py          # builds global FASTA + query table from updated_data.tsv
  make_steering_vectors.R       # per‑AID diffs → subject vectors → averaged steering vector
  run_64_codons_steering.sh     # 64‑codon sweep per AID; submit/download/merge/plot
  create_steering_table.r       # merge unsteered vs steered summaries into one TSV
  generate_codons_variants.py   # generate 64 codon sequences with margins
  plot_steered_vs_unsteered.r   # ggplot2 + ggrepel plot (labels AA_CODON)
```

---

## What I contributed

- Designed and implemented the steering pipeline across Bash, Python, and R
- Built the vector construction method and CLI orchestration for 64‑codon sweeps
- Integrated `evo_gcp` end‑to‑end and produced clear, reproducible figures

---

## Requirements

- Python: see `requirements.txt`
- R: `argparse`, `ggplot2`, `ggrepel`, `reticulate`
- EVO2 GCP CLI configured with access to submit/download jobs

## Public usage notice

Running full jobs requires:

- A Google Cloud Platform (GCP) project with billing enabled
- Appropriate IAM permissions/service account credentials
- Enabling required Google APIs; compute/storage usage incurs costs

---

## Next directions

- Try different nucleotide length steering vectors with pre-norm layer to see if a better steering vector is possible
- Add statistical tests for log‑likelihood shifts
- Use more subjects/genes in the calculation of the steering vector 
- Test if this current steering vector is applicable to other genes that the vector wasn't created from

## Contact

marcolanza@berkeley.edu • LinkedIn: marconlanza
