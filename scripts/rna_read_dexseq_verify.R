#!/usr/bin/env Rscript

# Verify that GENtle's flattened GFF and two-column count export form one
# DEXSeqDataSetFromHTSeq-compatible input pair. This script never installs or
# downloads R/Bioconductor packages.

required_packages <- c("DEXSeq")

dexseq_fail <- function(message) {
  message <- gsub("[\r\n\t]+", " ", as.character(message))
  cat("DEXSEQ_FAIL ", trimws(message), "\n", sep = "")
  quit(status = 1)
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  dexseq_fail(
    "usage=Rscript scripts/rna_read_dexseq_verify.R <flattened_gff_path> <count_file_path> [r_library_path ...]"
  )
}

flattened_gff_path <- args[[1]]
count_file_path <- args[[2]]
r_library_paths <- if (length(args) > 2) args[3:length(args)] else character()
existing_r_library_paths <- unique(r_library_paths[dir.exists(r_library_paths)])
if (length(existing_r_library_paths) > 0) {
  .libPaths(unique(c(
    normalizePath(existing_r_library_paths, winslash = "/", mustWork = TRUE),
    .libPaths()
  )))
}

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0) {
  dexseq_fail(paste0(
    "missing_required_package=",
    paste(missing_packages, collapse = ",")
  ))
}

tryCatch(
  {
    sample_data <- data.frame(
      condition = factor(c("pseudo_a", "pseudo_b")),
      row.names = c("pseudo_a", "pseudo_b")
    )
    suppressMessages(
      DEXSeq::DEXSeqDataSetFromHTSeq(
        countfiles = c(count_file_path, count_file_path),
        sampleData = sample_data,
        design = ~ sample + exon + condition:exon,
        flattenedfile = flattened_gff_path
      )
    )

    gff_lines <- readLines(flattened_gff_path, warn = FALSE)
    gff_columns <- strsplit(gff_lines, "\t", fixed = TRUE)
    feature_kinds <- vapply(
      gff_columns,
      function(columns) if (length(columns) >= 3) columns[[3]] else "",
      character(1)
    )
    aggregate_gene_count <- sum(feature_kinds == "aggregate_gene")
    exonic_part_count <- sum(feature_kinds == "exonic_part")

    count_rows <- utils::read.delim(
      count_file_path,
      header = FALSE,
      sep = "\t",
      quote = "",
      comment.char = "",
      stringsAsFactors = FALSE,
      col.names = c("feature_id", "count")
    )
    ordinary_rows <- !startsWith(count_rows$feature_id, "_")
    total_counts <- sum(as.numeric(count_rows$count[ordinary_rows]))

    cat(
      "DEXSEQ_OK genes=", aggregate_gene_count,
      " exonic_parts=", exonic_part_count,
      " total_counts=", format(total_counts, scientific = FALSE, trim = TRUE),
      "\n",
      sep = ""
    )
  },
  error = function(error) {
    dexseq_fail(conditionMessage(error))
  }
)
