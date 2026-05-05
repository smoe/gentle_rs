#!/usr/bin/env Rscript

# Prepare transcript-cluster and probeset-level Clariom D expression tables
# for the Rostock p73 publication resource.
#
# Required Bioconductor packages:
#   BiocManager::install(c("oligo", "limma", "pd.clariom.d.human"))
#
# Usage:
#   Rscript scripts/analyze_p73_clariomd_probe_level.R \
#     data/publication_resources/rostock_p73_clariomd_e_mtab_14704

args <- commandArgs(trailingOnly = TRUE)
dataset_dir <- if (length(args) >= 1) args[[1]] else "data/publication_resources/rostock_p73_clariomd_e_mtab_14704"
output_dir <- if (length(args) >= 2) args[[2]] else file.path(dataset_dir, "analysis", "clariomd_probe_level")
workspace_r_lib <- file.path(getwd(), ".r-lib")
if (dir.exists(workspace_r_lib)) {
  .libPaths(c(normalizePath(workspace_r_lib), .libPaths()))
}

required_packages <- c("oligo", "limma", "Biobase", "DBI", "RSQLite")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
pd_package <- Sys.getenv("GENTLE_CLARIOMD_PD_PACKAGE", "pd.clariom.d.human")
pd_missing <- !requireNamespace(pd_package, quietly = TRUE)
if (length(missing_packages) > 0 || pd_missing) {
  missing_all <- c(missing_packages, if (pd_missing) pd_package else character())
  stop(
    "Missing R/Bioconductor package(s): ", paste(missing_all, collapse = ", "), "\n",
    "Install with: if (!requireNamespace('BiocManager', quietly=TRUE)) install.packages('BiocManager'); ",
    "BiocManager::install(c('oligo', 'limma', 'pd.clariom.d.human'))",
    call. = FALSE
  )
}

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

message("Dataset directory: ", dataset_dir)
message("Output directory: ", output_dir)
message("Platform design package: ", pd_package)

write_tsv <- function(x, path) {
  utils::write.table(
    x,
    file = path,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = TRUE,
    na = ""
  )
}

read_sdrf_sample_table <- function(dataset_dir) {
  sdrf_path <- file.path(dataset_dir, "E-MTAB-14704.sdrf.txt")
  if (!file.exists(sdrf_path)) {
    warning("SDRF file not found at ", sdrf_path, "; deriving sample table from CEL filenames only.")
    sdrf <- data.frame()
    cel_names <- basename(list.files(dataset_dir, pattern = "\\.CEL$", ignore.case = TRUE, full.names = TRUE))
  } else {
    sdrf <- utils::read.delim(sdrf_path, check.names = FALSE, stringsAsFactors = FALSE)
    cel_columns <- names(sdrf)[vapply(sdrf, function(col) {
      any(grepl("\\.CEL$", basename(as.character(col)), ignore.case = TRUE))
    }, logical(1))]
    if (length(cel_columns) == 0) {
      stop("Could not find a CEL filename column in ", sdrf_path, call. = FALSE)
    }
    cel_names <- basename(as.character(sdrf[[cel_columns[[1]]]]))
  }

  condition <- sub("^P_SKMel29_", "", cel_names)
  condition <- sub("_[0-9]+\\.CEL$", "", condition, ignore.case = TRUE)
  replicate <- suppressWarnings(as.integer(sub("^.*_([0-9]+)\\.CEL$", "\\1", cel_names, ignore.case = TRUE)))
  sample_id <- sub("\\.CEL$", "", cel_names, ignore.case = TRUE)
  isoform <- ifelse(
    condition == "AdGFP",
    "control",
    ifelse(condition == "AdTAp73alpha", "TAp73alpha", ifelse(condition == "AdDNp73beta", "DNp73beta", condition))
  )
  cel_path <- file.path(dataset_dir, cel_names)
  table <- data.frame(
    sample_id = sample_id,
    file_name = cel_names,
    condition = condition,
    isoform = isoform,
    replicate = replicate,
    cel_path = cel_path,
    cel_exists = file.exists(cel_path),
    stringsAsFactors = FALSE
  )
  if (nrow(sdrf) == nrow(table)) {
    table$source_name <- if ("Source Name" %in% names(sdrf)) sdrf[["Source Name"]] else ""
    table$assay_name <- if ("Assay Name" %in% names(sdrf)) sdrf[["Assay Name"]] else ""
  }
  table[order(table$condition, table$replicate, table$file_name), , drop = FALSE]
}

sample_table <- read_sdrf_sample_table(dataset_dir)
rownames(sample_table) <- sample_table$sample_id
write_tsv(sample_table, file.path(output_dir, "sample_table.tsv"))

if (!all(sample_table$cel_exists)) {
  missing <- sample_table$file_name[!sample_table$cel_exists]
  stop(
    "Missing CEL file(s): ", paste(missing, collapse = ", "), "\n",
    "Fetch them with: cargo run --quiet --bin gentle_cli -- shell ",
    "'resources prepare-publication-dataset E-MTAB-14704 --category raw_microarray --download-files'",
    call. = FALSE
  )
}

cel_files <- sample_table$cel_path
names(cel_files) <- sample_table$sample_id

raw <- oligo::read.celfiles(cel_files, pkgname = pd_package)
Biobase::sampleNames(raw) <- sample_table$sample_id
Biobase::pData(raw) <- sample_table[Biobase::sampleNames(raw), , drop = FALSE]

load_platform_annotations <- function(pd_package) {
  annotations <- list(transcript = NULL, probeset = NULL)

  transcript_path <- system.file("extdata", "netaffxTranscript.rda", package = pd_package)
  if (nzchar(transcript_path) && file.exists(transcript_path)) {
    env <- new.env(parent = emptyenv())
    load(transcript_path, envir = env)
    if (exists("netaffxTranscript", envir = env, inherits = FALSE)) {
      transcript <- Biobase::pData(get("netaffxTranscript", envir = env))
      keep <- intersect(
        c(
          "transcriptclusterid", "seqname", "strand", "start", "stop", "totalprobes",
          "geneassignment", "mrnaassignment", "category", "locustype", "notes"
        ),
        names(transcript)
      )
      transcript <- transcript[, keep, drop = FALSE]
      transcript$feature_id <- as.character(transcript$transcriptclusterid)
      annotations$transcript <- transcript[, c("feature_id", setdiff(names(transcript), "feature_id")), drop = FALSE]
    }
  }

  sqlite_path <- system.file("extdata", "pd.clariom.d.human.sqlite", package = pd_package)
  if (nzchar(sqlite_path) && file.exists(sqlite_path)) {
    con <- DBI::dbConnect(RSQLite::SQLite(), sqlite_path)
    on.exit(DBI::dbDisconnect(con), add = TRUE)
    probeset <- DBI::dbGetQuery(
      con,
      paste(
        "select",
        "man_fsetid as feature_id,",
        "cast(fsetid as text) as platform_fsetid,",
        "transcript_cluster_id, exon_id, chrom, strand, start, stop,",
        "crosshyb_type, level, type, junction_start_edge, junction_stop_edge,",
        "junction_sequence, has_cds",
        "from featureSet"
      )
    )
    annotations$probeset <- probeset
  }

  annotations
}

platform_annotations <- load_platform_annotations(pd_package)

make_expression_table <- function(eset) {
  expr <- Biobase::exprs(eset)
  out <- data.frame(feature_id = rownames(expr), expr, check.names = FALSE)
  rownames(out) <- NULL
  out
}

make_feature_table <- function(eset, label, platform_annotations) {
  feature <- Biobase::fData(eset)
  if (nrow(feature) == 0) {
    feature <- data.frame(feature_id = rownames(Biobase::exprs(eset)), stringsAsFactors = FALSE)
  } else {
    feature <- data.frame(feature_id = rownames(feature), feature, check.names = FALSE)
    rownames(feature) <- NULL
  }
  if (label == "transcript_cluster" && !is.null(platform_annotations$transcript)) {
    return(merge(feature, platform_annotations$transcript, by = "feature_id", all.x = TRUE, sort = FALSE))
  }
  if (label == "probeset" && !is.null(platform_annotations$probeset)) {
    return(merge(feature, platform_annotations$probeset, by = "feature_id", all.x = TRUE, sort = FALSE))
  }
  feature
}

run_limma_contrasts <- function(expr_table, sample_table, label, output_dir) {
  expr <- as.matrix(expr_table[, setdiff(names(expr_table), "feature_id"), drop = FALSE])
  rownames(expr) <- expr_table$feature_id
  conditions <- factor(sample_table$condition, levels = unique(sample_table$condition))
  design <- stats::model.matrix(~ 0 + conditions)
  colnames(design) <- make.names(levels(conditions))

  available <- colnames(design)
  contrast_defs <- c(
    TAp73alpha_vs_AdGFP = "AdTAp73alpha-AdGFP",
    DNp73beta_vs_AdGFP = "AdDNp73beta-AdGFP",
    TAp73alpha_vs_DNp73beta = "AdTAp73alpha-AdDNp73beta"
  )
  contrast_defs <- contrast_defs[vapply(strsplit(contrast_defs, "-", fixed = TRUE), function(parts) {
    all(parts %in% available)
  }, logical(1))]
  if (length(contrast_defs) == 0) {
    warning("No standard p73 contrasts could be made for ", label)
    return(list())
  }

  fit <- limma::lmFit(expr, design)
  contrast_matrix <- limma::makeContrasts(contrasts = contrast_defs, levels = design)
  fit2 <- limma::eBayes(limma::contrasts.fit(fit, contrast_matrix))

  results <- list()
  for (contrast in colnames(contrast_matrix)) {
    tab <- limma::topTable(fit2, coef = contrast, number = Inf, sort.by = "none")
    tab <- data.frame(feature_id = rownames(tab), tab, check.names = FALSE)
    rownames(tab) <- NULL
    path <- file.path(output_dir, paste0(label, "_limma_", contrast, ".tsv"))
    write_tsv(tab, path)
    add_artifact(path, paste(label, contrast, "limma differential-expression table"))
    results[[contrast]] <- tab
  }
  results
}

guess_group_column <- function(feature_table) {
  candidates <- names(feature_table)[grepl("transcript.*cluster|transcriptcluster|gene|mrna", names(feature_table), ignore.case = TRUE)]
  candidates <- setdiff(candidates, "feature_id")
  if (length(candidates) == 0) "" else candidates[[1]]
}

first_gene_symbol <- function(geneassignment) {
  vapply(as.character(geneassignment), function(value) {
    if (is.na(value) || !nzchar(value)) {
      return("")
    }
    first_assignment <- strsplit(value, " /// ", fixed = TRUE)[[1]][[1]]
    parts <- strsplit(first_assignment, " // ", fixed = TRUE)[[1]]
    if (length(parts) >= 2) parts[[2]] else first_assignment
  }, character(1))
}

write_splice_hint_tables <- function(limma_results, feature_table, output_dir) {
  group_column <- guess_group_column(feature_table)
  if (!nzchar(group_column)) {
    warning("No obvious transcript/gene grouping column found for probeset splice-hint summaries.")
    return(invisible(FALSE))
  }
  for (contrast in names(limma_results)) {
    joined <- merge(limma_results[[contrast]], feature_table, by = "feature_id", all.x = TRUE)
    group <- as.character(joined[[group_column]])
    keep <- !is.na(group) & nzchar(group)
    joined <- joined[keep, , drop = FALSE]
    group <- group[keep]
    if (nrow(joined) == 0) {
      next
    }
    split_logfc <- split(joined$logFC, group)
    summary <- data.frame(
      group_id = names(split_logfc),
      grouping_column = group_column,
      probeset_count = vapply(split_logfc, length, integer(1)),
      logFC_min = vapply(split_logfc, min, numeric(1), na.rm = TRUE),
      logFC_max = vapply(split_logfc, max, numeric(1), na.rm = TRUE),
      logFC_range = vapply(split_logfc, function(x) diff(range(x, na.rm = TRUE)), numeric(1)),
      logFC_sd = vapply(split_logfc, stats::sd, numeric(1), na.rm = TRUE),
      stringsAsFactors = FALSE
    )
    summary <- summary[summary$probeset_count > 1, , drop = FALSE]
    summary <- summary[order(-summary$logFC_range, -summary$probeset_count), , drop = FALSE]
    if (!is.null(platform_annotations$transcript)) {
      tx <- platform_annotations$transcript
      tx$group_id <- tx$feature_id
      tx$gene_symbol <- first_gene_symbol(tx$geneassignment)
      keep_columns <- intersect(
        c("group_id", "gene_symbol", "seqname", "start", "stop", "locustype", "category"),
        names(tx)
      )
      summary <- merge(summary, tx[, keep_columns, drop = FALSE], by = "group_id", all.x = TRUE, sort = FALSE)
    }
    path <- file.path(output_dir, paste0("probeset_splice_hint_", contrast, ".tsv"))
    write_tsv(summary, path)
    add_artifact(path, paste("probeset-level within-group logFC heterogeneity for", contrast))
  }
  invisible(TRUE)
}

safe_file_token <- function(value) {
  token <- gsub("[^A-Za-z0-9_.-]+", "_", value)
  token <- gsub("^_+|_+$", "", token)
  if (!nzchar(token)) "contrast" else token
}

json_escape <- function(value) {
  value <- as.character(value)
  value <- gsub("\\\\", "\\\\\\\\", value)
  value <- gsub("\"", "\\\\\"", value)
  value <- gsub("\n", "\\\\n", value, fixed = TRUE)
  value <- gsub("\r", "\\\\r", value, fixed = TRUE)
  value <- gsub("\t", "\\\\t", value, fixed = TRUE)
  value
}

json_string <- function(value) {
  paste0("\"", json_escape(value), "\"")
}

json_array <- function(values) {
  values <- values[!is.na(values) & nzchar(values)]
  paste0("[", paste(vapply(values, json_string, character(1)), collapse = ", "), "]")
}

manifest_track_entry_json <- function(entry) {
  paste0(
    "    {\n",
    "      \"contrast\": ", json_string(entry$contrast), ",\n",
    "      \"level\": ", json_string(entry$level), ",\n",
    "      \"path\": ", json_string(entry$path), ",\n",
    "      \"row_count\": ", as.integer(entry$row_count), "\n",
    "    }"
  )
}

write_microarray_track_manifest <- function(limma_results, feature_table, output_dir) {
  if (length(limma_results) == 0) {
    warning("No probeset limma results available; microarray track manifest not written.")
    return(invisible(NULL))
  }
  required_coordinate_columns <- c("chrom", "start", "stop")
  if (!all(required_coordinate_columns %in% names(feature_table))) {
    warning(
      "Probeset feature annotation lacks coordinate columns (needed: ",
      paste(required_coordinate_columns, collapse = ", "),
      "); microarray track manifest not written."
    )
    return(invisible(NULL))
  }

  feature_augmented <- feature_table
  if (!"gene_symbol" %in% names(feature_augmented)) {
    feature_augmented$gene_symbol <- ""
  }
  if (!is.null(platform_annotations$transcript) && "transcript_cluster_id" %in% names(feature_augmented)) {
    tx <- platform_annotations$transcript
    if ("geneassignment" %in% names(tx)) {
      tx$transcript_cluster_id <- as.character(tx$feature_id)
      tx$gene_symbol_from_transcript <- first_gene_symbol(tx$geneassignment)
      feature_augmented <- merge(
        feature_augmented,
        tx[, c("transcript_cluster_id", "gene_symbol_from_transcript"), drop = FALSE],
        by = "transcript_cluster_id",
        all.x = TRUE,
        sort = FALSE
      )
      missing_gene <- !nzchar(feature_augmented$gene_symbol) &
        !is.na(feature_augmented$gene_symbol_from_transcript) &
        nzchar(feature_augmented$gene_symbol_from_transcript)
      feature_augmented$gene_symbol[missing_gene] <- feature_augmented$gene_symbol_from_transcript[missing_gene]
      feature_augmented$gene_symbol_from_transcript <- NULL
    }
  }

  track_dir <- file.path(output_dir, "microarray_tracks")
  dir.create(track_dir, recursive = TRUE, showWarnings = FALSE)
  track_entries <- data.frame(
    contrast = character(),
    level = character(),
    path = character(),
    row_count = integer(),
    stringsAsFactors = FALSE
  )

  for (contrast in names(limma_results)) {
    joined <- merge(limma_results[[contrast]], feature_augmented, by = "feature_id", all.x = TRUE, sort = FALSE)
    row <- data.frame(
      chrom = as.character(joined$chrom),
      start_1based = suppressWarnings(as.integer(joined$start)),
      end_1based = suppressWarnings(as.integer(joined$stop)),
      strand = as.character(joined$strand),
      feature_id = as.character(joined$feature_id),
      transcript_cluster_id = if ("transcript_cluster_id" %in% names(joined)) as.character(joined$transcript_cluster_id) else "",
      exon_id = if ("exon_id" %in% names(joined)) as.character(joined$exon_id) else "",
      probe_type = if ("type" %in% names(joined)) as.character(joined$type) else if ("crosshyb_type" %in% names(joined)) as.character(joined$crosshyb_type) else "",
      logFC = joined$logFC,
      AveExpr = joined$AveExpr,
      P.Value = joined$P.Value,
      adj.P.Val = joined$adj.P.Val,
      gene_symbol = as.character(joined$gene_symbol),
      junction_start_edge = if ("junction_start_edge" %in% names(joined)) as.character(joined$junction_start_edge) else "",
      junction_stop_edge = if ("junction_stop_edge" %in% names(joined)) as.character(joined$junction_stop_edge) else "",
      junction_sequence = if ("junction_sequence" %in% names(joined)) as.character(joined$junction_sequence) else "",
      has_cds = if ("has_cds" %in% names(joined)) as.character(joined$has_cds) else "",
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
    keep <- !is.na(row$chrom) & nzchar(row$chrom) &
      !is.na(row$start_1based) & !is.na(row$end_1based) &
      row$start_1based > 0 & row$end_1based >= row$start_1based
    row <- row[keep, , drop = FALSE]
    rel_path <- file.path("microarray_tracks", paste0("clariomd_", safe_file_token(contrast), "_probeset_track.tsv"))
    path <- file.path(output_dir, rel_path)
    write_tsv(row, path)
    add_artifact(path, paste("GENtle microarray probeset track TSV for", contrast))
    track_entries <- rbind(
      track_entries,
      data.frame(
        contrast = contrast,
        level = "probeset",
        path = rel_path,
        row_count = nrow(row),
        stringsAsFactors = FALSE
      )
    )
  }

  coordinate_system <- Sys.getenv("GENTLE_CLARIOMD_COORDINATE_SYSTEM", "unverified")
  supported_genome_ids <- unlist(strsplit(Sys.getenv("GENTLE_CLARIOMD_SUPPORTED_GENOME_IDS", ""), ",", fixed = TRUE))
  supported_genome_ids <- trimws(supported_genome_ids)
  supported_genome_ids <- supported_genome_ids[nzchar(supported_genome_ids)]
  warnings <- character()
  if (!nzchar(coordinate_system) || identical(coordinate_system, "unverified") || length(supported_genome_ids) == 0) {
    warnings <- c(
      warnings,
      "Clariom D platform coordinate build was not verified; GENtle will reject projection until coordinate_system and supported_genome_ids are curated."
    )
  }

  manifest_path <- file.path(output_dir, "clariomd_microarray_track_manifest.json")
  track_json <- if (nrow(track_entries) == 0) {
    ""
  } else {
    paste(vapply(seq_len(nrow(track_entries)), function(i) manifest_track_entry_json(track_entries[i, ]), character(1)), collapse = ",\n")
  }
  manifest_lines <- c(
    "{",
    paste0("  \"schema\": ", json_string("gentle.microarray_track_manifest.v1"), ","),
    paste0("  \"dataset\": ", json_string("E-MTAB-14704"), ","),
    paste0("  \"platform\": ", json_string("Clariom D human"), ","),
    paste0("  \"normalization\": ", json_string("RMA"), ","),
    paste0("  \"coordinate_system\": ", json_string(coordinate_system), ","),
    paste0("  \"supported_genome_ids\": ", json_array(supported_genome_ids), ","),
    paste0("  \"contrast_order\": ", json_array(track_entries$contrast), ","),
    "  \"contrasts\": [",
    track_json,
    "  ],",
    paste0("  \"warnings\": ", json_array(warnings)),
    "}"
  )
  writeLines(manifest_lines, manifest_path, useBytes = TRUE)
  add_artifact(manifest_path, "GENtle microarray track manifest for probeset-level Clariom D contrasts")
  invisible(manifest_path)
}

analysis_manifest <- data.frame(
  artifact = character(),
  path = character(),
  description = character(),
  stringsAsFactors = FALSE
)
add_artifact <- function(path, description) {
  analysis_manifest <<- rbind(
    analysis_manifest,
    data.frame(artifact = basename(path), path = path, description = description, stringsAsFactors = FALSE)
  )
}

run_rma_target <- function(raw, target, label) {
  message("Running oligo::rma(target = '", target, "')")
  eset <- oligo::rma(raw, target = target)
  expr_path <- file.path(output_dir, paste0(label, "_expression_rma.tsv"))
  feature_path <- file.path(output_dir, paste0(label, "_feature_annotation.tsv"))
  expr_table <- make_expression_table(eset)
  feature_table <- make_feature_table(eset, label, platform_annotations)
  write_tsv(expr_table, expr_path)
  write_tsv(feature_table, feature_path)
  add_artifact(expr_path, paste(label, "RMA expression matrix"))
  add_artifact(feature_path, paste(label, "feature annotation exported from ExpressionSet featureData"))
  contrasts <- run_limma_contrasts(expr_table, sample_table, label, output_dir)
  list(eset = eset, expr = expr_table, feature = feature_table, contrasts = contrasts)
}

transcript_result <- run_rma_target(raw, "core", "transcript_cluster")

probeset_result <- tryCatch(
  run_rma_target(raw, "probeset", "probeset"),
  error = function(e) {
    warning("Probeset-level RMA failed: ", conditionMessage(e))
    NULL
  }
)
if (!is.null(probeset_result)) {
  write_splice_hint_tables(probeset_result$contrasts, probeset_result$feature, output_dir)
  write_microarray_track_manifest(probeset_result$contrasts, probeset_result$feature, output_dir)
}

session_path <- file.path(output_dir, "sessionInfo.txt")
capture.output(utils::sessionInfo(), file = session_path)
add_artifact(file.path(output_dir, "sample_table.tsv"), "sample metadata derived from SDRF and CEL filenames")
add_artifact(session_path, "R session information for reproducibility")
write_tsv(analysis_manifest, file.path(output_dir, "analysis_manifest.tsv"))

message("Finished Clariom D probe/probeset-level analysis.")
message("Primary outputs: ", output_dir)
