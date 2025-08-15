#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(argparse)
})

main <- function() {
  p <- ArgumentParser()
  p$add_argument("--jobs_dir", required = TRUE, help = "jobs/<job>-<version> directory")
  p$add_argument("--scale", required = TRUE, help = "Steering scale (e.g., 1, -1, 10)")
  p$add_argument("--out", required = TRUE, help = "Output TSV path for steered vs unsteered table")
  args <- p$parse_args()

  jobs_dir_abs <- normalizePath(args$jobs_dir, mustWork = TRUE)
  output_dir_abs <- file.path(jobs_dir_abs, "output")
  if (!dir.exists(output_dir_abs)) {
    # Some downloads place summaries directly under jobs_dir
    output_dir_abs <- jobs_dir_abs
  }

  # Simple, fixed filenames with minimal fallbacks
  unsteered_path <- file.path(output_dir_abs, "input_summary_unsteered.txt")
  if (!file.exists(unsteered_path)) stop(sprintf("Missing unsteered summary: %s", unsteered_path))

  steered_path <- file.path(output_dir_abs, sprintf("input_summary_scale_%s.txt", args$scale))
  if (!file.exists(steered_path)) {
    # Accept numeric-equivalent names like 20 vs 20.0
    all_scale_files <- list.files(output_dir_abs, pattern = "^input_summary_scale_[-0-9eE+.]+\\.txt$", full.names = TRUE)
    if (length(all_scale_files) > 0) {
      want <- suppressWarnings(as.numeric(args$scale))
      if (!is.na(want)) {
        # Extract numeric part and pick closest
        nums <- suppressWarnings(as.numeric(sub("^input_summary_scale_([-0-9eE+.]+)\\.txt$", "\\1", basename(all_scale_files))))
        diffs <- abs(nums - want)
        pick <- which.min(diffs)
        if (length(pick) == 1 && is.finite(diffs[pick])) {
          steered_path <- all_scale_files[pick]
        }
      }
    }
    if (!file.exists(steered_path)) {
      steered_path <- file.path(output_dir_abs, "input_summary.txt")
    }
  }
  if (!file.exists(steered_path)) stop(sprintf("Missing steered summary (searched scale %s): %s", args$scale, steered_path))

  un_df <- read.delim(unsteered_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  st_df <- read.delim(steered_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

  required_cols <- c("seq_id", "total_log_likelihood")
  if (!all(required_cols %in% names(un_df))) stop("Unsteered file missing required columns: seq_id, total_log_likelihood")
  if (!all(required_cols %in% names(st_df))) stop("Steered file missing required columns: seq_id, total_log_likelihood")

  un_df <- un_df[, required_cols]
  st_df <- st_df[, required_cols]
  names(un_df)[2] <- "unsteered"
  names(st_df)[2] <- "steered"

  merged_df <- merge(un_df, st_df, by = "seq_id", all = FALSE)
  merged_df <- merged_df[order(merged_df$seq_id), ]

  write.table(merged_df, args$out, sep = "\t", row.names = FALSE, quote = FALSE)
  cat(sprintf("Wrote %s with %d rows\n", args$out, nrow(merged_df)))
}

if (identical(environment(), globalenv())) main()


