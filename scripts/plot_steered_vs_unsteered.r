#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(argparse)
})

main <- function() {
  p <- ArgumentParser()
  p$add_argument("--table", required = TRUE, help = "TSV with columns: seq_id, unsteered, steered")
  p$add_argument("--title", default = "Steered vs Unsteered")
  p$add_argument("--outdir", default = "figures/steering_checks")
  p$add_argument("--codon_table", default = "input/codon_table", help = "TSV with columns: codon, aa (for labels)")
  p$add_argument("--src_codon", default = "", help = "Source/reference codon to highlight")
  p$add_argument("--tgt_codon", default = "", help = "Target/mutation codon to highlight")
  args <- p$parse_args()

  df <- read.delim(args$table, stringsAsFactors = FALSE)
  # Derive codon as the last token of seq_id
  split_list <- strsplit(df$seq_id, "_")
  codons <- vapply(split_list, function(x) x[length(x)], character(1))
  parts <- data.frame(codon = toupper(codons), stringsAsFactors = FALSE)
  # Map codon -> amino acid using codon_table
  aa_label <- parts$codon
  ct <- tryCatch(read.delim(args$codon_table, header = TRUE, sep = "\t", stringsAsFactors = FALSE), error = function(e) NULL)
  if (!is.null(ct)) {
    ln <- tolower(names(ct))
    cod_col <- which(ln == "codon")
    aa_col <- which(ln == "aa")
    if (length(cod_col) == 1 && length(aa_col) == 1) {
      # Ensure uppercase for matching
      ct_codons <- toupper(ct[[cod_col]])
      ct_aa <- ct[[aa_col]]
      aa_label <- ct_aa[match(parts$codon, ct_codons)]
      # Fallback to codon if no match
      aa_label[is.na(aa_label)] <- parts$codon[is.na(aa_label)]
    }
  }
  df$label <- paste0(aa_label, "_", parts$codon)

  # Highlight source/target codons if provided
  which_vec <- rep("Other", nrow(df))
  if (nzchar(args$src_codon)) which_vec[parts$codon == toupper(args$src_codon)] <- "Source"
  if (nzchar(args$tgt_codon)) which_vec[parts$codon == toupper(args$tgt_codon)] <- "Target"
  df$which <- factor(which_vec, levels = c("Other", "Source", "Target"))
  dir.create(args$outdir, recursive = TRUE, showWarnings = FALSE)

  p <- ggplot(df, aes(x = unsteered, y = steered, label = label)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray60") +
    geom_point(aes(color = which), alpha = 0.85, size = 2.2) +
    scale_color_manual(values = c(Other = "#2c7bb6", Source = "#d95f02", Target = "#1b9e77"), name = "Highlighted") +
    ggrepel::geom_text_repel(size = 3, max.overlaps = 100) +
    theme_minimal() +
    labs(title = args$title, x = "Unsteered total log-likelihood", y = "Steered total log-likelihood")

  ofn <- file.path(args$outdir, "steered_vs_unsteered.pdf")
  ggsave(ofn, p, width = 8, height = 6)
  cat("Saved:", ofn, "\n")
}

if (identical(environment(), globalenv())) main()


