#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(argparse)
  library(reticulate)
})

# --- Helpers ---------------------------------------------------------------
l2_norm <- function(x) sqrt(sum(x * x))

normalize_vec <- function(x, eps = 1e-12) {
  n <- l2_norm(x)
  if (n < eps) return(x)
  x / n
}

## Filenames now assumed to always include suffixes (e.g., _unsteered.npy)
path_with_suffix <- function(base_no_ext) {
  paste0(base_no_ext, "_unsteered.npy")
}

write_subject_vector <- function(out_dir, aid, layer, embed, win_len, vec) {
  # Write one value per line (single column), like the final steering vector
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  out_raw <- file.path(out_dir, sprintf("non_norm_%s_%s_len%d.tsv", aid, embed, win_len))
  write.table(matrix(vec, ncol = 1), out_raw, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
  message(sprintf("Wrote subject vector: %s", out_raw))
}

write_steering_avg <- function(out_dir, layer, embed, win_len, vec) {
  # Only write the single-column TSV used by evo_gcp
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  out_raw <- file.path(out_dir, sprintf("non_norm_steering_%s_len%d.tsv", embed, win_len))
  write.table(matrix(vec, ncol = 1), out_raw, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
  message(sprintf("Wrote steering average: %s", out_raw))
}

get_window_rows <- function(n_rows, embed, win_len) {
  # Preserve prior behavior: start at 1 for pre_norm, 2 for mlp_l3
  start_idx <- if (identical(embed, "mlp_l3")) 2L else 1L
  end_idx <- min(n_rows, start_idx + as.integer(win_len) - 1L)
  if (end_idx < start_idx) return(integer(0))
  seq.int(start_idx, end_idx)
}

# --- Main -----------------------------------------------------------------
main <- function() {
  parser <- ArgumentParser()
  parser$add_argument("--jobs_dir", required = TRUE, help = "Path to jobs/<job>-<version>")
  parser$add_argument("--updated_data", required = TRUE, help = "Path to updated_data.tsv")
  parser$add_argument("--layer", type = "integer", default = 14)
  parser$add_argument("--embed", default = "pre_norm", help = "Embedding key: pre_norm or mlp_l3")
  parser$add_argument("--vectors_dir", default = "vectors/non_norm_diff_vec", help = "Directory to write subject and steering vectors")
  parser$add_argument("--win_len", type = "integer", default = 200, help = "Window length (number of positions from start=1)")
  args <- parser$parse_args()

  message("[Info] Loading numpy via reticulate…")
  np <- import("numpy")

  # Read TSV (expects headers)
  message(sprintf("[Info] Reading updated data: %s", args$updated_data))
  df <- read.delim(args$updated_data, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

  # Collect normalized diffs per AID
  aid_to_diffs <- list()
  seen_shape <- FALSE

  for (i in seq_len(nrow(df))) {
    aid <- tolower(df$aid[i])
    src_codon <- toupper(df$src_codon[i])
    tgt_codon <- toupper(df$tgt_codon[i])
    if (is.na(aid) || is.na(src_codon) || is.na(tgt_codon)) next

    src_id <- sprintf("%s_%s", aid, src_codon)
    tgt_id <- sprintf("%s_%s", aid, tgt_codon)
    src_base <- file.path(args$jobs_dir, "output",
                          sprintf("input_%s_embeddings_blocks_%d_%s", src_id, args$layer, args$embed))
    tgt_base <- file.path(args$jobs_dir, "output",
                          sprintf("input_%s_embeddings_blocks_%d_%s", tgt_id, args$layer, args$embed))
    src_path <- path_with_suffix(src_base)
    tgt_path <- path_with_suffix(tgt_base)

    if (!file.exists(src_path) || !file.exists(tgt_path)) {
      message(sprintf("[Warn] Missing pair for %s: %s or %s", aid, basename(src_path), basename(tgt_path)))
      next
    }

    # Load arrays and enforce 2D shape [positions, dims] per updated evo outputs
    src_arr <- np$load(src_path)
    tgt_arr <- np$load(tgt_path)
    src_r <- py_to_r(src_arr)
    tgt_r <- py_to_r(tgt_arr)
    dsrc <- dim(src_r)
    dtgt <- dim(tgt_r)
    if (!seen_shape) {
      message(sprintf("[Info] Observed shapes src %s, tgt %s (expect c(201,D))",
                      paste(dsrc, collapse = "x"), paste(dtgt, collapse = "x")))
      seen_shape <- TRUE
    }
    if (is.null(dsrc) || length(dsrc) != 2) {
      message(sprintf("[Warn] Unsupported src dims for %s: %s; skipping", src_id, paste(dsrc, collapse = "x")))
      next
    }
    if (is.null(dtgt) || length(dtgt) != 2) {
      message(sprintf("[Warn] Unsupported tgt dims for %s: %s; skipping", tgt_id, paste(dtgt, collapse = "x")))
      next
    }
    src_mat <- src_r
    tgt_mat <- tgt_r

    # Minimal length check: ensure we can take win_len from the chosen start
    min_needed <- if (identical(args$embed, "mlp_l3")) 1L + args$win_len else args$win_len
    if (nrow(src_mat) < min_needed || nrow(tgt_mat) < min_needed) {
      message(sprintf("[Warn] Too few positions; got src %d, tgt %d; skipping", nrow(src_mat), nrow(tgt_mat)))
      next
    }

    # Select window rows and compute per-position diffs
    rows <- get_window_rows(nrow(src_mat), args$embed, args$win_len)
    src_win <- src_mat[rows, , drop = FALSE]
    tgt_win <- tgt_mat[rows, , drop = FALSE]

    # Compute per-position diffs without row-wise normalization, then
    # average across positions and normalize the per-AID vector
    diff_win <- tgt_win - src_win                        # matrix [positions, dims]
    diff_vec <- normalize_vec(colMeans(diff_win))

    # Store a single vector per aid (simplified: assumes one row per aid)
    aid_to_diffs[[aid]] <- diff_vec
  }

  # Write subject vectors (already one vector per aid)
  message("[Info] Writing subject vectors…")
  subject_vectors <- aid_to_diffs
  for (aid in names(subject_vectors)) {
    write_subject_vector(args$vectors_dir, aid, args$layer, args$embed, args$win_len, subject_vectors[[aid]])
  }
  message(sprintf("[Info] Wrote %d subject vectors", length(subject_vectors)))

  # Average steering vector across AIDs
  if (length(subject_vectors) > 0) {
    message("[Info] Computing average steering vector across subjects…")
    m <- do.call(rbind, subject_vectors)
    steer <- normalize_vec(colMeans(m))
    write_steering_avg(args$vectors_dir, args$layer, args$embed, args$win_len, steer)
  } else {
    message("[Warn] No subject vectors built; skipping steering average.")
  }

  message("[Done] Steering vectors created.")
}

if (identical(environment(), globalenv())) main()

