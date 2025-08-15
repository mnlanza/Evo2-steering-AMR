## Steering workflow (layer14)

This workflow builds per-subject average vectors, creates a single steering vector, runs steered jobs across scales, evaluates Δlog-likelihood for codon position 83 (positive strand), and plots decay curves for two layers: `mlp_l3` (blocks.14.mlp.l3) and `pre_norm` (blocks.14.pre_norm).

### Layout

```
steering/
  scripts/
    lib_common.sh
    lib_plot.R
    build_subject_vectors.py
    average_vectors.py
    run_steering.sh
    eval_pos83.R
    plot_decay_curves.R
  vectors/
  results/
  plots/
  .gitignore
```

### Conventions

- **vectors**: `vectors/subject_<aid>_layer14_<embed>.tsv.gz`, `vectors/steering_avg_layer14_<embed>.tsv.gz`
- **results**: `results/pos83_ll_scale_<scale>_layer14_<embed>.tsv`
- **layers**: keep separate; never mix `mlp_l3` and `pre_norm`.

### Steps

1) Build per-subject average vectors (normalize each input vector to unit length before averaging):

```bash
python steering/scripts/build_subject_vectors.py \
  --raw-dir steering/vectors/raw \
  --embed mlp_l3 \
  --layer 14 \
  --out-dir steering/vectors

python steering/scripts/build_subject_vectors.py \
  --raw-dir steering/vectors/raw \
  --embed pre_norm \
  --layer 14 \
  --out-dir steering/vectors
```

2) Make the 8-subject steering average vector (and optionally emit a raw numeric TSV for evo_gcp):

```bash
python steering/scripts/average_vectors.py \
  --vectors-dir steering/vectors \
  --embed mlp_l3 \
  --layer 14 \
  --emit-raw

python steering/scripts/average_vectors.py \
  --vectors-dir steering/vectors \
  --embed pre_norm \
  --layer 14 \
  --emit-raw
```

3) Run steered jobs (scales: -10, -1, -0.1, 0, 0.1, 1, 10) via evo_gcp. The script will submit once per embed with all scales and parse downloaded summaries into per-scale TSVs expected by downstream steps:

```bash
bash steering/scripts/run_steering.sh
```

4) Evaluate Δlog-likelihood at position 83 (positive strand):

```bash
Rscript steering/scripts/eval_pos83.R --embed mlp_l3 --layer 14 --results-dir steering/results
Rscript steering/scripts/eval_pos83.R --embed pre_norm --layer 14 --results-dir steering/results
```

5) Plot decay curves and ΔlogL vs scale:

```bash
Rscript steering/scripts/plot_decay_curves.R --embed mlp_l3 --layer 14 \
  --results-dir steering/results --out-dir steering/plots

Rscript steering/scripts/plot_decay_curves.R --embed pre_norm --layer 14 \
  --results-dir steering/results --out-dir steering/plots
```

### Notes

- If you already have plotting helpers in `scripts/plot_layer.R`, `steering/scripts/lib_plot.R` can source from there or be used standalone.
- Set `JOBS_DIR` in your shell to point to your evo_gcp job outputs if the default is not correct.
- The evo_gcp steering interface is documented here: [evo2_gcp README](https://github.com/eitanyaffe/evo2_gcp/tree/main). Ensure your local install is up to date (git pull) so `--steering_layer`, `--steering_vector_file`, and `--steering_scales` are supported.

