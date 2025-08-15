#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(argparse)
  library(reticulate)
})

# Compute per-position distances between two embedding sequences
compute_euclidean_distances <- function(ref_mat, mut_mat) {
  sapply(1:nrow(ref_mat), function(i) {
    src_i <- ref_mat[i, ]
    tgt_i <- mut_mat[i, ]
    sqrt(sum((src_i - tgt_i)^2))
  })
}

compute_log_invert_cosine_sim <- function(ref_mat, mut_mat) {
  sapply(1:nrow(ref_mat), function(i) {
    src_i <- ref_mat[i, ]
    tgt_i <- mut_mat[i, ]
    x <- sum(src_i * tgt_i) / (sqrt(sum(src_i^2)) * sqrt(sum(tgt_i^2)))
    x <- min(x, 0.999999)  # Prevent log2(0)
    log2(1 - x)
  })
}

plot_series <- function(values, title, ylab, out_path) {
  df <- data.frame(pos = seq_along(values), val = values)
  p <- ggplot(df, aes(x = pos, y = val)) +
    geom_point(color = "#3366cc", alpha = 0.7, size = 1.1) +
    theme_minimal() +
    labs(title = title, x = "Nucleotide Position", y = ylab)
  dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
  ggsave(out_path, p, width = 8, height = 5)
}

main <- function() {
  parser <- ArgumentParser()
  parser$add_argument("--updated_data", required = TRUE, help = "Path to updated_data.tsv")
  parser$add_argument("--embed_dir", required = TRUE, help = "Path to jobs/<job>-<version> (downloaded)")
  parser$add_argument("--layer", type = "integer", default = 14)
  parser$add_argument("--embeds", default = "pre_norm,mlp_l3", help = "Comma-separated embed keys")
  parser$add_argument("--output_dir", default = "figures/pairwise_checks")
  args <- parser$parse_args()

  np <- import("numpy")

  df <- read.delim(args$updated_data, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  embeds <- strsplit(args$embeds, ",", fixed = TRUE)[[1]]

  for (row_i in seq_len(nrow(df))) {
    aid <- tolower(df$aid[row_i])
    src_codon <- df$src_codon[row_i]
    tgt_codon <- df$tgt_codon[row_i]
    if (is.na(aid) || is.na(src_codon) || is.na(tgt_codon)) next

    src_id <- paste0(aid, "_", src_codon)
    tgt_id <- paste0(aid, "_", tgt_codon)

    for (embed in embeds) {
      src_path <- file.path(args$embed_dir, "output",
                            sprintf("input_%s_embeddings_blocks_%d_%s.npy", src_id, args$layer, embed))
      tgt_path <- file.path(args$embed_dir, "output",
                            sprintf("input_%s_embeddings_blocks_%d_%s.npy", tgt_id, args$layer, embed))
      if (!file.exists(src_path) || !file.exists(tgt_path)) {
        message(sprintf("Skipping %s (%s): missing %s or %s", aid, embed, basename(src_path), basename(tgt_path)))
        next
      }
      src_arr <- np$load(src_path)
      tgt_arr <- np$load(tgt_path)
      src_mat <- py_to_r(src_arr)[1, , ]
      tgt_mat <- py_to_r(tgt_arr)[1, , ]

      eu <- compute_euclidean_distances(src_mat, tgt_mat)
      cs <- compute_log_invert_cosine_sim(src_mat, tgt_mat)

      out_base <- file.path(args$output_dir, aid, paste0("layer", args$layer), embed)
      plot_series(eu, sprintf("Euclidean distance: %s vs %s (%s)", src_codon, tgt_codon, aid),
                  "Euclidean distance", file.path(out_base, sprintf("euclid_%s_%s.pdf", aid, embed)))
      plot_series(cs, sprintf("Log inverse cosine: %s vs %s (%s)", src_codon, tgt_codon, aid),
                  "log2(1 - cos sim)", file.path(out_base, sprintf("cosine_%s_%s.pdf", aid, embed)))
    }
  }
}

if (identical(environment(), globalenv())) main()

