#!/usr/bin/env Rscript

# Generic Affymetrix CEL probe/probeset-region helper for GENtle.
#
# This script is intentionally external to the Rust engine. It consumes explicit
# CEL files and optional sample metadata, runs an R/oligo RMA backend, and writes
# chromosome-ordered intensity tables that GENtle can ingest or inspect later.
# It never downloads or installs R/Bioconductor packages.

usage <- function(status = 0) {
  cat(
    "Usage:\n",
    "  Rscript scripts/probe_regions_oligo.R --cel sample1.CEL [--cel sample2.CEL ...]\n",
    "    (--gene SYMBOL | --locus CHR:START-END | --transcript-cluster-id ID | --probeset-id ID ...)\n",
    "    [--metadata samples.tsv] [--sample-column NAME] [--condition-column NAME] [--block-column NAME]\n",
    "    [--platform-package pd.clariom.d.human] [--normalization rma]\n",
    "    [--coordinate-system hg38] [--genome-build GRCh38]\n",
    "    [--target probeset|transcript_cluster] [--contrast A-B] --output DIR [--cache-dir DIR]\n",
    "\n",
    "Outputs include region_intensity_chrom_order.csv, feature/expression TSVs, limma contrast TSVs,\n",
    "normalized_feature_matrix_manifest.json, provenance.json, and sessionInfo.txt.\n",
    sep = ""
  )
  quit(status = status)
}

empty_args <- function() {
  list(
    cel = character(),
    metadata = "",
    sample_column = "",
    condition_column = "",
    block_column = "",
    platform_package = Sys.getenv("GENTLE_PROBE_REGION_PD_PACKAGE", "pd.clariom.d.human"),
    platform_name = Sys.getenv("GENTLE_PROBE_REGION_PLATFORM", "Clariom_D_Human"),
    coordinate_system = Sys.getenv("GENTLE_PROBE_REGION_COORDINATE_SYSTEM", ""),
    genome_build = Sys.getenv("GENTLE_PROBE_REGION_GENOME_BUILD", ""),
    normalization = "rma",
    targets = "probeset",
    contrasts = character(),
    output = "analysis/probe_regions",
    cache_dir = "",
    genes = character(),
    loci = character(),
    transcript_cluster_ids = character(),
    probeset_ids = character(),
    allow_all_features = FALSE
  )
}

take_value <- function(args, i, flag) {
  if (i + 1 > length(args)) {
    stop("Missing value after ", flag, call. = FALSE)
  }
  args[[i + 1]]
}

split_csv <- function(value) {
  value <- trimws(value)
  if (!nzchar(value)) {
    character()
  } else {
    trimws(strsplit(value, ",", fixed = TRUE)[[1]])
  }
}

parse_args <- function(argv) {
  if (length(argv) == 0 || any(argv %in% c("--help", "-h"))) {
    usage(0)
  }
  out <- empty_args()
  i <- 1
  while (i <= length(argv)) {
    flag <- argv[[i]]
    if (flag == "--cel") {
      out$cel <- c(out$cel, take_value(argv, i, flag))
      i <- i + 2
    } else if (flag == "--metadata" || flag == "--sdrf") {
      out$metadata <- take_value(argv, i, flag)
      i <- i + 2
    } else if (flag == "--sample-column") {
      out$sample_column <- take_value(argv, i, flag)
      i <- i + 2
    } else if (flag == "--condition-column") {
      out$condition_column <- take_value(argv, i, flag)
      i <- i + 2
    } else if (flag == "--block-column" || flag == "--batch-column") {
      out$block_column <- take_value(argv, i, flag)
      i <- i + 2
    } else if (flag == "--platform-package") {
      out$platform_package <- take_value(argv, i, flag)
      i <- i + 2
    } else if (flag == "--platform-name") {
      out$platform_name <- take_value(argv, i, flag)
      i <- i + 2
    } else if (flag == "--coordinate-system") {
      out$coordinate_system <- take_value(argv, i, flag)
      i <- i + 2
    } else if (flag == "--genome-build" || flag == "--reference-genome-id") {
      out$genome_build <- take_value(argv, i, flag)
      i <- i + 2
    } else if (flag == "--normalization") {
      out$normalization <- tolower(take_value(argv, i, flag))
      i <- i + 2
    } else if (flag == "--target") {
      out$targets <- unique(c(out$targets, split_csv(take_value(argv, i, flag))))
      i <- i + 2
    } else if (flag == "--targets") {
      out$targets <- unique(c(out$targets, split_csv(take_value(argv, i, flag))))
      i <- i + 2
    } else if (flag == "--contrast") {
      out$contrasts <- c(out$contrasts, take_value(argv, i, flag))
      i <- i + 2
    } else if (flag == "--output" || flag == "--output-dir") {
      out$output <- take_value(argv, i, flag)
      i <- i + 2
    } else if (flag == "--cache-dir") {
      out$cache_dir <- take_value(argv, i, flag)
      i <- i + 2
    } else if (flag == "--gene") {
      out$genes <- c(out$genes, take_value(argv, i, flag))
      i <- i + 2
    } else if (flag == "--genes") {
      out$genes <- c(out$genes, split_csv(take_value(argv, i, flag)))
      i <- i + 2
    } else if (flag == "--locus") {
      out$loci <- c(out$loci, take_value(argv, i, flag))
      i <- i + 2
    } else if (flag == "--loci") {
      out$loci <- c(out$loci, split_csv(take_value(argv, i, flag)))
      i <- i + 2
    } else if (flag == "--transcript-cluster-id") {
      out$transcript_cluster_ids <- c(out$transcript_cluster_ids, take_value(argv, i, flag))
      i <- i + 2
    } else if (flag == "--transcript-cluster-ids") {
      out$transcript_cluster_ids <- c(out$transcript_cluster_ids, split_csv(take_value(argv, i, flag)))
      i <- i + 2
    } else if (flag == "--probeset-id") {
      out$probeset_ids <- c(out$probeset_ids, take_value(argv, i, flag))
      i <- i + 2
    } else if (flag == "--probeset-ids") {
      out$probeset_ids <- c(out$probeset_ids, split_csv(take_value(argv, i, flag)))
      i <- i + 2
    } else if (flag == "--allow-all-features") {
      out$allow_all_features <- TRUE
      i <- i + 1
    } else {
      stop("Unknown option: ", flag, call. = FALSE)
    }
  }
  out$targets <- unique(out$targets[nzchar(out$targets)])
  if (length(out$targets) == 0) {
    out$targets <- "probeset"
  }
  out$genes <- unique(out$genes[nzchar(out$genes)])
  out$loci <- unique(out$loci[nzchar(out$loci)])
  out$transcript_cluster_ids <- unique(out$transcript_cluster_ids[nzchar(out$transcript_cluster_ids)])
  out$probeset_ids <- unique(out$probeset_ids[nzchar(out$probeset_ids)])
  out
}

args <- parse_args(commandArgs(trailingOnly = TRUE))

workspace_r_lib <- file.path(getwd(), ".r-lib")
if (dir.exists(workspace_r_lib)) {
  .libPaths(c(normalizePath(workspace_r_lib), .libPaths()))
}

required_packages <- c("oligo", "limma", "Biobase", "DBI", "RSQLite")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (!requireNamespace(args$platform_package, quietly = TRUE)) {
  missing_packages <- c(missing_packages, args$platform_package)
}
if (length(missing_packages) > 0) {
  stop(
    "Missing R/Bioconductor package(s): ", paste(unique(missing_packages), collapse = ", "), "\n",
    "Install outside GENtle, for example: if (!requireNamespace('BiocManager', quietly=TRUE)) ",
    "install.packages('BiocManager'); BiocManager::install(c('oligo','limma','",
    args$platform_package, "'))",
    call. = FALSE
  )
}
if (!identical(args$normalization, "rma")) {
  stop("This helper currently supports --normalization rma only", call. = FALSE)
}
if (length(args$cel) == 0) {
  stop("At least one --cel path is required", call. = FALSE)
}
has_selector <- length(args$genes) > 0 ||
  length(args$loci) > 0 ||
  length(args$transcript_cluster_ids) > 0 ||
  length(args$probeset_ids) > 0
if (!has_selector && !args$allow_all_features) {
  stop("At least one selector is required; use --allow-all-features intentionally for full exports", call. = FALSE)
}

dir.create(args$output, recursive = TRUE, showWarnings = FALSE)
if (nzchar(args$cache_dir)) {
  dir.create(args$cache_dir, recursive = TRUE, showWarnings = FALSE)
}

safe_token <- function(value) {
  token <- gsub("[^A-Za-z0-9_.-]+", "_", value)
  token <- gsub("^_+|_+$", "", token)
  if (!nzchar(token)) "value" else token
}

write_tsv <- function(x, path) {
  utils::write.table(x, file = path, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE, na = "")
}

write_csv <- function(x, path) {
  utils::write.csv(x, file = path, quote = TRUE, row.names = FALSE, na = "")
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

json_bool <- function(value) {
  if (isTRUE(value)) "true" else "false"
}

file_info_json <- function(path) {
  exists <- file.exists(path)
  info <- if (exists) file.info(path) else NULL
  size <- if (exists) as.character(info$size) else "null"
  mtime <- if (exists) as.character(as.integer(as.POSIXct(info$mtime))) else "null"
  paste0(
    "{",
    "\"path\":", json_string(path), ",",
    "\"exists\":", json_bool(exists), ",",
    "\"size_bytes\":", size, ",",
    "\"modified_unix_seconds\":", mtime,
    "}"
  )
}

metadata_delimiter <- function(path) {
  if (grepl("\\.csv$", path, ignore.case = TRUE)) "," else "\t"
}

header_key <- function(value) {
  tolower(gsub("[^A-Za-z0-9]+", "", value))
}

find_column <- function(columns, requested, aliases) {
  if (nzchar(requested)) {
    idx <- which(header_key(columns) == header_key(requested))
    if (length(idx) == 0) stop("Requested metadata column not found: ", requested, call. = FALSE)
    return(columns[[idx[[1]]]])
  }
  alias_keys <- header_key(aliases)
  idx <- which(header_key(columns) %in% alias_keys)
  if (length(idx) == 0) "" else columns[[idx[[1]]]]
}

read_metadata_table <- function(path) {
  if (!nzchar(path)) {
    return(NULL)
  }
  if (!file.exists(path)) {
    stop("Metadata file does not exist: ", path, call. = FALSE)
  }
  utils::read.table(
    path,
    sep = metadata_delimiter(path),
    header = TRUE,
    quote = "\"",
    comment.char = "",
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
}

make_sample_table <- function(args) {
  cel_paths <- normalizePath(args$cel, mustWork = FALSE)
  missing <- cel_paths[!file.exists(cel_paths)]
  if (length(missing) > 0) {
    stop("Missing CEL file(s): ", paste(missing, collapse = ", "), call. = FALSE)
  }
  metadata <- read_metadata_table(args$metadata)
  if (is.null(metadata)) {
    sample_id <- sub("\\.CEL$", "", basename(cel_paths), ignore.case = TRUE)
    return(data.frame(
      sample_id = make.unique(sample_id),
      file_name = basename(cel_paths),
      condition = sample_id,
      block = "",
      cel_path = cel_paths,
      stringsAsFactors = FALSE
    ))
  }
  columns <- names(metadata)
  sample_column <- find_column(columns, args$sample_column, c("file", "cel", "cel_file", "cel file", "array data file", "sample", "sample_name", "source name"))
  condition_column <- find_column(columns, args$condition_column, c("condition", "group", "treatment", "sample group", "factor value condition", "factor value treatment"))
  block_column <- find_column(columns, args$block_column, c("block", "batch", "replicate", "biological replicate"))
  if (!nzchar(sample_column)) {
    stop("Could not infer a sample/CEL column from metadata; pass --sample-column", call. = FALSE)
  }
  sample_values <- as.character(metadata[[sample_column]])
  base_map <- stats::setNames(cel_paths, basename(cel_paths))
  exact_map <- stats::setNames(cel_paths, cel_paths)
  cel_for_sample <- ifelse(sample_values %in% names(exact_map), exact_map[sample_values], "")
  missing_exact <- !nzchar(cel_for_sample)
  cel_for_sample[missing_exact & sample_values %in% names(base_map)] <- base_map[sample_values[missing_exact & sample_values %in% names(base_map)]]
  missing_match <- !nzchar(cel_for_sample)
  if (any(missing_match)) {
    metadata_dir <- dirname(normalizePath(args$metadata, mustWork = FALSE))
    candidate <- file.path(metadata_dir, sample_values[missing_match])
    candidate_exists <- file.exists(candidate)
    missing_idx <- which(missing_match)
    if (any(candidate_exists)) {
      cel_for_sample[missing_idx[candidate_exists]] <- normalizePath(candidate[candidate_exists])
    }
  }
  keep <- nzchar(cel_for_sample)
  if (!any(keep)) {
    stop("Metadata rows did not match any supplied CEL files", call. = FALSE)
  }
  condition <- if (nzchar(condition_column)) as.character(metadata[[condition_column]]) else sub("\\.CEL$", "", basename(cel_for_sample), ignore.case = TRUE)
  block <- if (nzchar(block_column)) as.character(metadata[[block_column]]) else rep("", nrow(metadata))
  sample_id <- make.unique(sub("\\.CEL$", "", basename(cel_for_sample), ignore.case = TRUE))
  table <- data.frame(
    sample_id = sample_id[keep],
    file_name = basename(cel_for_sample[keep]),
    condition = condition[keep],
    block = block[keep],
    cel_path = cel_for_sample[keep],
    stringsAsFactors = FALSE
  )
  table[order(table$condition, table$file_name), , drop = FALSE]
}

first_gene_symbol <- function(geneassignment) {
  vapply(as.character(geneassignment), function(value) {
    if (is.na(value) || !nzchar(value)) return("")
    first_assignment <- strsplit(value, " /// ", fixed = TRUE)[[1]][[1]]
    parts <- strsplit(first_assignment, " // ", fixed = TRUE)[[1]]
    if (length(parts) >= 2) parts[[2]] else first_assignment
  }, character(1))
}

load_platform_annotations <- function(pd_package) {
  annotations <- list(transcript = NULL, probeset = NULL)
  transcript_path <- system.file("extdata", "netaffxTranscript.rda", package = pd_package)
  if (nzchar(transcript_path) && file.exists(transcript_path)) {
    env <- new.env(parent = emptyenv())
    load(transcript_path, envir = env)
    if (exists("netaffxTranscript", envir = env, inherits = FALSE)) {
      transcript <- Biobase::pData(get("netaffxTranscript", envir = env))
      keep <- intersect(c("transcriptclusterid", "seqname", "strand", "start", "stop", "totalprobes", "geneassignment", "mrnaassignment", "category", "locustype", "notes"), names(transcript))
      transcript <- transcript[, keep, drop = FALSE]
      transcript$feature_id <- as.character(transcript$transcriptclusterid)
      transcript$gene_symbol <- if ("geneassignment" %in% names(transcript)) first_gene_symbol(transcript$geneassignment) else ""
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
        "select man_fsetid as feature_id, cast(fsetid as text) as platform_fsetid,",
        "transcript_cluster_id, exon_id, chrom, strand, start, stop,",
        "crosshyb_type, level, type, junction_start_edge, junction_stop_edge,",
        "junction_sequence, has_cds from featureSet"
      )
    )
    annotations$probeset <- probeset
  }
  annotations
}

make_expression_table <- function(eset) {
  expr <- Biobase::exprs(eset)
  out <- data.frame(feature_id = rownames(expr), expr, check.names = FALSE)
  rownames(out) <- NULL
  out
}

make_feature_table <- function(eset, label, annotations) {
  feature <- Biobase::fData(eset)
  if (nrow(feature) == 0) {
    feature <- data.frame(feature_id = rownames(Biobase::exprs(eset)), stringsAsFactors = FALSE)
  } else {
    feature <- data.frame(feature_id = rownames(feature), feature, check.names = FALSE)
    rownames(feature) <- NULL
  }
  if (label == "transcript_cluster" && !is.null(annotations$transcript)) {
    return(merge(feature, annotations$transcript, by = "feature_id", all.x = TRUE, sort = FALSE))
  }
  if (label == "probeset" && !is.null(annotations$probeset)) {
    feature <- merge(feature, annotations$probeset, by = "feature_id", all.x = TRUE, sort = FALSE)
    if (!is.null(annotations$transcript) && "transcript_cluster_id" %in% names(feature)) {
      tx <- annotations$transcript
      tx$transcript_cluster_id <- as.character(tx$feature_id)
      keep <- intersect(c("transcript_cluster_id", "gene_symbol", "geneassignment", "seqname"), names(tx))
      feature <- merge(feature, tx[, keep, drop = FALSE], by = "transcript_cluster_id", all.x = TRUE, sort = FALSE)
    }
    return(feature)
  }
  feature
}

chrom_key <- function(value) {
  tolower(sub("^chr", "", as.character(value), ignore.case = TRUE))
}

chrom_order <- function(value) {
  key <- chrom_key(value)
  number <- suppressWarnings(as.integer(key))
  ifelse(!is.na(number), number, ifelse(key == "x", 23, ifelse(key == "y", 24, ifelse(key %in% c("m", "mt"), 25, 1000))))
}

parse_locus <- function(value) {
  m <- regexec("^([^:]+):(\\d+)-(\\d+)$", gsub(",", "", value))
  parts <- regmatches(value, m)[[1]]
  if (length(parts) != 4) {
    stop("Invalid locus selector '", value, "'; expected CHR:START-END", call. = FALSE)
  }
  data.frame(chrom = parts[[2]], start = as.integer(parts[[3]]), end = as.integer(parts[[4]]), stringsAsFactors = FALSE)
}

selector_mask <- function(table, args) {
  mask <- rep(FALSE, nrow(table))
  if (length(args$probeset_ids) > 0) {
    mask <- mask | as.character(table$feature_id) %in% args$probeset_ids
  }
  if (length(args$transcript_cluster_ids) > 0 && "transcript_cluster_id" %in% names(table)) {
    mask <- mask | as.character(table$transcript_cluster_id) %in% args$transcript_cluster_ids
  }
  if (length(args$genes) > 0) {
    gene_query <- toupper(args$genes)
    gene_columns <- intersect(c("gene_symbol", "gene_symbol.x", "gene_symbol.y", "geneassignment"), names(table))
    for (column in gene_columns) {
      values <- toupper(as.character(table[[column]]))
      mask <- mask | values %in% gene_query
      for (gene in gene_query) {
        mask <- mask | grepl(paste0("(^|[^A-Z0-9])", gene, "([^A-Z0-9]|$)"), values)
      }
    }
  }
  if (length(args$loci) > 0 && all(c("chrom", "start", "stop") %in% names(table))) {
    for (locus_value in args$loci) {
      locus <- parse_locus(locus_value)
      mask <- mask | (chrom_key(table$chrom) == chrom_key(locus$chrom) & table$stop >= locus$start & table$start <= locus$end)
    }
  }
  mask
}

add_condition_summaries <- function(region, sample_table) {
  sample_columns <- intersect(sample_table$sample_id, names(region))
  for (condition in unique(sample_table$condition)) {
    ids <- intersect(sample_table$sample_id[sample_table$condition == condition], sample_columns)
    if (length(ids) == 0) next
    token <- safe_token(condition)
    values <- as.matrix(region[, ids, drop = FALSE])
    region[[paste0("mean_log2_", token)]] <- rowMeans(values, na.rm = TRUE)
    region[[paste0("sd_log2_", token)]] <- apply(values, 1, stats::sd, na.rm = TRUE)
  }
  region
}

default_contrasts <- function(sample_table) {
  conditions <- unique(sample_table$condition)
  if (length(conditions) < 2) return(character())
  baseline <- conditions[[1]]
  paste(conditions[-1], baseline, sep = "-")
}

add_logfc_summaries <- function(region, sample_table, contrasts) {
  for (contrast in contrasts) {
    parts <- strsplit(contrast, "-", fixed = TRUE)[[1]]
    if (length(parts) != 2) next
    numerator <- parts[[1]]
    denominator <- parts[[2]]
    n_ids <- intersect(sample_table$sample_id[sample_table$condition == numerator], names(region))
    d_ids <- intersect(sample_table$sample_id[sample_table$condition == denominator], names(region))
    if (length(n_ids) == 0 || length(d_ids) == 0) next
    n_mean <- rowMeans(as.matrix(region[, n_ids, drop = FALSE]), na.rm = TRUE)
    d_mean <- rowMeans(as.matrix(region[, d_ids, drop = FALSE]), na.rm = TRUE)
    region[[paste0("log2FC_", safe_token(contrast))]] <- n_mean - d_mean
  }
  region
}

run_limma <- function(expr_table, sample_table, contrasts, label, output_dir) {
  if (length(contrasts) == 0) return(character())
  expr <- as.matrix(expr_table[, setdiff(names(expr_table), "feature_id"), drop = FALSE])
  rownames(expr) <- expr_table$feature_id
  conditions <- factor(sample_table$condition, levels = unique(sample_table$condition))
  design <- stats::model.matrix(~ 0 + conditions)
  colnames(design) <- make.names(levels(conditions))
  parsed <- strsplit(contrasts, "-", fixed = TRUE)
  valid <- vapply(parsed, function(parts) length(parts) == 2 && all(make.names(parts) %in% colnames(design)), logical(1))
  contrasts <- contrasts[valid]
  if (length(contrasts) == 0) return(character())
  contrast_exprs <- vapply(parsed[valid], function(parts) paste(make.names(parts), collapse = "-"), character(1))
  names(contrast_exprs) <- safe_token(contrasts)
  fit <- limma::lmFit(expr, design)
  contrast_matrix <- limma::makeContrasts(contrasts = contrast_exprs, levels = design)
  fit2 <- limma::eBayes(limma::contrasts.fit(fit, contrast_matrix))
  paths <- character()
  for (contrast in colnames(contrast_matrix)) {
    tab <- limma::topTable(fit2, coef = contrast, number = Inf, sort.by = "none")
    tab <- data.frame(feature_id = rownames(tab), tab, check.names = FALSE)
    rownames(tab) <- NULL
    path <- file.path(output_dir, paste0(label, "_limma_", safe_token(contrast), ".tsv"))
    write_tsv(tab, path)
    paths <- c(paths, path)
  }
  paths
}

make_region_table <- function(expr_table, feature_table, sample_table, args) {
  joined <- merge(feature_table, expr_table, by = "feature_id", all.x = TRUE, sort = FALSE)
  if (!"chrom" %in% names(joined) && "seqname" %in% names(joined)) joined$chrom <- joined$seqname
  if (!"start" %in% names(joined)) joined$start <- NA_integer_
  if (!"stop" %in% names(joined)) joined$stop <- NA_integer_
  if (!"strand" %in% names(joined)) joined$strand <- ""
  if (!"transcript_cluster_id" %in% names(joined)) joined$transcript_cluster_id <- ""
  if (!"exon_id" %in% names(joined)) joined$exon_id <- ""
  if (!"totalprobes" %in% names(joined)) joined$totalprobes <- ""
  if (!"gene_symbol" %in% names(joined)) joined$gene_symbol <- ""
  mask <- if (args$allow_all_features) rep(TRUE, nrow(joined)) else selector_mask(joined, args)
  joined <- joined[mask, , drop = FALSE]
  sample_ids <- intersect(sample_table$sample_id, names(joined))
  region <- data.frame(
    chromosome = as.character(joined$chrom),
    start = suppressWarnings(as.integer(joined$start)),
    stop = suppressWarnings(as.integer(joined$stop)),
    strand = as.character(joined$strand),
    probeset_or_region_id = as.character(joined$feature_id),
    transcript_cluster_id = as.character(joined$transcript_cluster_id),
    exon_id = as.character(joined$exon_id),
    number_of_probes = as.character(joined$totalprobes),
    gene_symbol = as.character(joined$gene_symbol),
    joined[, sample_ids, drop = FALSE],
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  region <- add_condition_summaries(region, sample_table)
  region <- add_logfc_summaries(region, sample_table, args$contrasts)
  region[order(chrom_order(region$chromosome), region$start, region$stop, region$probeset_or_region_id), , drop = FALSE]
}

sample_table <- make_sample_table(args)
rownames(sample_table) <- sample_table$sample_id
write_tsv(sample_table, file.path(args$output, "sample_table.tsv"))
if (length(args$contrasts) == 0) {
  args$contrasts <- default_contrasts(sample_table)
}

cel_files <- sample_table$cel_path
names(cel_files) <- sample_table$sample_id
message("Reading ", length(cel_files), " CEL file(s) with platform package ", args$platform_package)
raw <- oligo::read.celfiles(cel_files, pkgname = args$platform_package)
Biobase::sampleNames(raw) <- sample_table$sample_id
Biobase::pData(raw) <- sample_table[Biobase::sampleNames(raw), , drop = FALSE]
annotations <- load_platform_annotations(args$platform_package)

artifact_paths <- c(file.path(args$output, "sample_table.tsv"))
region_written <- FALSE
for (target in args$targets) {
  label <- if (target == "transcript_cluster") "transcript_cluster" else "probeset"
  rma_target <- if (label == "transcript_cluster") "core" else "probeset"
  message("Running oligo::rma(target = '", rma_target, "')")
  eset <- oligo::rma(raw, target = rma_target)
  expr_table <- make_expression_table(eset)
  feature_table <- make_feature_table(eset, label, annotations)
  expr_path <- file.path(args$output, paste0(label, "_expression_rma.tsv"))
  feature_path <- file.path(args$output, paste0(label, "_feature_annotation.tsv"))
  write_tsv(expr_table, expr_path)
  write_tsv(feature_table, feature_path)
  artifact_paths <- c(artifact_paths, expr_path, feature_path)
  artifact_paths <- c(artifact_paths, run_limma(expr_table, sample_table, args$contrasts, label, args$output))
  region <- make_region_table(expr_table, feature_table, sample_table, args)
  region_path <- if (!region_written) {
    file.path(args$output, "region_intensity_chrom_order.csv")
  } else {
    file.path(args$output, paste0(label, "_region_intensity_chrom_order.csv"))
  }
  write_csv(region, region_path)
  artifact_paths <- c(artifact_paths, region_path)
  region_written <- TRUE
}

session_path <- file.path(args$output, "sessionInfo.txt")
capture.output(utils::sessionInfo(), file = session_path)
artifact_paths <- c(artifact_paths, session_path)

manifest_path <- file.path(args$output, "normalized_feature_matrix_manifest.json")
manifest_lines <- c(
  "{",
  paste0("  \"schema\": ", json_string("gentle.probe_region_normalized_matrix_manifest.v1"), ","),
  paste0("  \"platform\": ", json_string(args$platform_name), ","),
  paste0("  \"platform_package\": ", json_string(args$platform_package), ","),
  paste0("  \"coordinate_system\": ", json_string(args$coordinate_system), ","),
  paste0("  \"genome_build\": ", json_string(args$genome_build), ","),
  paste0("  \"normalization\": ", json_string(args$normalization), ","),
  paste0("  \"targets\": ", json_array(args$targets), ","),
  paste0("  \"contrasts\": ", json_array(args$contrasts), ","),
  paste0("  \"cel_files\": [", paste(vapply(cel_files, file_info_json, character(1)), collapse = ", "), "],"),
  paste0("  \"metadata\": ", if (nzchar(args$metadata)) file_info_json(args$metadata) else "null", ","),
  paste0("  \"artifacts\": ", json_array(artifact_paths)),
  "}"
)
writeLines(manifest_lines, manifest_path, useBytes = TRUE)
artifact_paths <- c(artifact_paths, manifest_path)

provenance_path <- file.path(args$output, "provenance.json")
provenance_lines <- c(
  "{",
  paste0("  \"schema\": ", json_string("gentle.probe_region_backend_provenance.v1"), ","),
  paste0("  \"backend\": ", json_string("r_oligo"), ","),
  paste0("  \"command\": ", json_string(paste(commandArgs(FALSE), collapse = " ")), ","),
  paste0("  \"platform_package\": ", json_string(args$platform_package), ","),
  paste0("  \"coordinate_system\": ", json_string(args$coordinate_system), ","),
  paste0("  \"genome_build\": ", json_string(args$genome_build), ","),
  paste0("  \"normalization\": ", json_string(args$normalization), ","),
  paste0("  \"output_dir\": ", json_string(args$output), ","),
  paste0("  \"artifacts\": ", json_array(artifact_paths)),
  "}"
)
writeLines(provenance_lines, provenance_path, useBytes = TRUE)

message("Finished GENtle probe-region oligo backend.")
message("Primary output: ", file.path(args$output, "region_intensity_chrom_order.csv"))
