#!/usr/bin/env Rscript

# Generic Affymetrix 3' IVT/CDF CEL helper for GENtle.
#
# This script is intentionally external to the Rust engine. It consumes explicit
# CEL files plus a local CDF package/file and optional local annotation table or
# Bioconductor annotation package, runs an R/affy RMA backend, and writes the
# same helper-output contract used by the Clariom/ST probe-region bridge.
# It never downloads or installs R/Bioconductor packages.

usage <- function(status = 0) {
  cat(
    "Usage:\n",
    "  Rscript scripts/probe_regions_affy.R --cel sample1.CEL [--cel sample2.CEL ...]\n",
    "    (--gene SYMBOL | --locus CHR:START-END | --probeset-id ID ...)\n",
    "    [--metadata samples.tsv] [--sample-column NAME] [--condition-column NAME] [--block-column NAME]\n",
    "    [--cdf-package hgu133plus2cdf] [--cdf-name hgu133plus2]\n",
    "    [--annotation-package hgu133plus2.db] [--annotation-library probeset_coordinates.csv]\n",
    "    [--platform-name HG_U133_Plus_2] [--normalization rma]\n",
    "    [--coordinate-system hg19] [--genome-build GRCh37]\n",
    "    [--contrast A-B] --output DIR [--cache-dir DIR] [--allow-all-features]\n",
    "\n",
    "Outputs include region_intensity_chrom_order.csv, expression/annotation TSVs,\n",
    "normalized_feature_matrix_manifest.json, provenance.json, and sessionInfo.txt.\n",
    "A CDF supplies probe grouping; genome coordinates require an annotation package\n",
    "with CHRLOC/CHRLOCEND-style columns or a user-supplied coordinate table.\n",
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
    cdf_package = Sys.getenv("GENTLE_PROBE_REGION_CDF_PACKAGE", ""),
    cdf_name = Sys.getenv("GENTLE_PROBE_REGION_CDF_NAME", ""),
    annotation_package = Sys.getenv("GENTLE_PROBE_REGION_ANNOTATION_PACKAGE", ""),
    annotation_library = "",
    platform_name = Sys.getenv("GENTLE_PROBE_REGION_PLATFORM", ""),
    coordinate_system = Sys.getenv("GENTLE_PROBE_REGION_COORDINATE_SYSTEM", ""),
    genome_build = Sys.getenv("GENTLE_PROBE_REGION_GENOME_BUILD", ""),
    normalization = "rma",
    contrasts = character(),
    output = "analysis/probe_regions/affy",
    cache_dir = "",
    genes = character(),
    loci = character(),
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
  if (!nzchar(value)) character() else trimws(strsplit(value, ",", fixed = TRUE)[[1]])
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
    } else if (flag == "--cdf-package") {
      out$cdf_package <- take_value(argv, i, flag)
      i <- i + 2
    } else if (flag == "--cdf-name") {
      out$cdf_name <- take_value(argv, i, flag)
      i <- i + 2
    } else if (flag == "--annotation-package") {
      out$annotation_package <- take_value(argv, i, flag)
      i <- i + 2
    } else if (flag == "--annotation-library") {
      out$annotation_library <- take_value(argv, i, flag)
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
  out$genes <- unique(out$genes[nzchar(out$genes)])
  out$loci <- unique(out$loci[nzchar(out$loci)])
  out$probeset_ids <- unique(out$probeset_ids[nzchar(out$probeset_ids)])
  out
}

args <- parse_args(commandArgs(trailingOnly = TRUE))

workspace_r_lib <- file.path(getwd(), ".r-lib")
if (dir.exists(workspace_r_lib)) {
  .libPaths(c(normalizePath(workspace_r_lib), .libPaths()))
}

required_packages <- c("affy", "limma", "Biobase")
if (nzchar(args$annotation_package)) {
  required_packages <- c(required_packages, "AnnotationDbi")
}
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (nzchar(args$cdf_package) && !requireNamespace(args$cdf_package, quietly = TRUE)) {
  missing_packages <- c(missing_packages, args$cdf_package)
}
if (nzchar(args$annotation_package) && !requireNamespace(args$annotation_package, quietly = TRUE)) {
  missing_packages <- c(missing_packages, args$annotation_package)
}
if (length(missing_packages) > 0) {
  stop(
    "Missing R/Bioconductor package(s): ", paste(unique(missing_packages), collapse = ", "), "\n",
    "Install outside GENtle, for example: if (!requireNamespace('BiocManager', quietly=TRUE)) ",
    "install.packages('BiocManager'); BiocManager::install(c('affy','limma','",
    paste(unique(missing_packages), collapse = "','"), "'))",
    call. = FALSE
  )
}
if (!identical(args$normalization, "rma")) {
  stop("This helper currently supports --normalization rma only", call. = FALSE)
}
if (length(args$cel) == 0) {
  stop("At least one --cel path is required", call. = FALSE)
}
has_selector <- length(args$genes) > 0 || length(args$loci) > 0 || length(args$probeset_ids) > 0
if (!has_selector && !args$allow_all_features) {
  stop("At least one selector is required; use --allow-all-features intentionally for full exports", call. = FALSE)
}
if (!nzchar(args$cdf_package) && !nzchar(args$cdf_name)) {
  stop("A legacy 3' IVT run requires --cdf-package or --cdf-name", call. = FALSE)
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

json_string <- function(value) paste0("\"", json_escape(value), "\"")

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

read_delimited_table <- function(path) {
  if (!file.exists(path)) {
    stop("Table does not exist: ", path, call. = FALSE)
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

read_metadata_table <- function(path) {
  if (!nzchar(path)) NULL else read_delimited_table(path)
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
  data.frame(chromosome = parts[[2]], start = as.integer(parts[[3]]), stop = as.integer(parts[[4]]), stringsAsFactors = FALSE)
}

infer_cdf_name <- function(cdf_package, cdf_name) {
  if (nzchar(cdf_name)) return(cdf_name)
  name <- sub("cdf$", "", cdf_package, ignore.case = TRUE)
  if (!nzchar(name)) cdf_package else name
}

available_annotation_columns <- function(package) {
  if (!nzchar(package) || !requireNamespace("AnnotationDbi", quietly = TRUE)) {
    return(character())
  }
  suppressWarnings(tryCatch(AnnotationDbi::columns(package), error = function(e) character()))
}

select_annotation_package <- function(package, probeset_ids) {
  if (!nzchar(package)) {
    return(data.frame(probeset_or_region_id = probeset_ids, stringsAsFactors = FALSE))
  }
  columns <- available_annotation_columns(package)
  wanted <- intersect(c("SYMBOL", "GENENAME", "CHR", "CHRLOC", "CHRLOCEND"), columns)
  if (length(wanted) == 0) {
    warning("Annotation package '", package, "' exposes no SYMBOL/CHR/CHRLOC-style columns")
    return(data.frame(probeset_or_region_id = probeset_ids, stringsAsFactors = FALSE))
  }
  keytype <- if ("PROBEID" %in% AnnotationDbi::keytypes(package)) "PROBEID" else AnnotationDbi::keytypes(package)[[1]]
  selected <- suppressWarnings(AnnotationDbi::select(package, keys = probeset_ids, keytype = keytype, columns = wanted))
  if (nrow(selected) == 0) {
    return(data.frame(probeset_or_region_id = probeset_ids, stringsAsFactors = FALSE))
  }
  id_column <- names(selected)[[1]]
  selected$probeset_or_region_id <- as.character(selected[[id_column]])
  if ("SYMBOL" %in% names(selected)) selected$gene_symbol <- as.character(selected$SYMBOL)
  if ("CHR" %in% names(selected)) selected$chromosome <- as.character(selected$CHR)
  if ("CHRLOC" %in% names(selected)) selected$start <- abs(suppressWarnings(as.integer(selected$CHRLOC)))
  if ("CHRLOCEND" %in% names(selected)) selected$stop <- abs(suppressWarnings(as.integer(selected$CHRLOCEND)))
  if (!"stop" %in% names(selected) && "start" %in% names(selected)) selected$stop <- selected$start
  if ("CHRLOC" %in% names(selected)) selected$strand <- ifelse(suppressWarnings(as.integer(selected$CHRLOC)) < 0, "-", "+")
  keep <- intersect(c("probeset_or_region_id", "chromosome", "start", "stop", "strand", "gene_symbol", "GENENAME"), names(selected))
  selected <- selected[, keep, drop = FALSE]
  selected <- selected[!duplicated(selected$probeset_or_region_id), , drop = FALSE]
  selected
}

normalize_annotation_table <- function(path) {
  if (!nzchar(path)) {
    return(data.frame())
  }
  ext <- tolower(tools::file_ext(path))
  if (!ext %in% c("csv", "tsv", "txt")) {
    warning("Annotation library '", path, "' is not a CSV/TSV coordinate table; treating it as backend-only CDF support")
    return(data.frame())
  }
  table <- read_delimited_table(path)
  cols <- names(table)
  pick <- function(requested, aliases) find_column(cols, "", c(requested, aliases))
  id_col <- pick("probeset_or_region_id", c("probeset_id", "probe_set_id", "probe set id", "id", "feature_id"))
  chr_col <- pick("chromosome", c("chrom", "seqname", "chr"))
  start_col <- pick("start", c("start_1based", "genomic_start", "chrom_start"))
  stop_col <- pick("stop", c("end", "end_1based", "genomic_end", "chrom_stop"))
  strand_col <- pick("strand", c("direction"))
  gene_col <- pick("gene_symbol", c("symbol", "gene", "gene_assignment"))
  if (!nzchar(id_col)) {
    stop("Annotation coordinate table requires a probeset id column", call. = FALSE)
  }
  out <- data.frame(probeset_or_region_id = as.character(table[[id_col]]), stringsAsFactors = FALSE)
  if (nzchar(chr_col)) out$chromosome <- as.character(table[[chr_col]])
  if (nzchar(start_col)) out$start <- suppressWarnings(as.integer(table[[start_col]]))
  if (nzchar(stop_col)) out$stop <- suppressWarnings(as.integer(table[[stop_col]]))
  if (nzchar(strand_col)) out$strand <- as.character(table[[strand_col]])
  if (nzchar(gene_col)) out$gene_symbol <- as.character(table[[gene_col]])
  out
}

merge_annotation <- function(probeset_ids, args) {
  base <- data.frame(probeset_or_region_id = probeset_ids, stringsAsFactors = FALSE)
  package_annotation <- select_annotation_package(args$annotation_package, probeset_ids)
  table_annotation <- normalize_annotation_table(args$annotation_library)
  out <- merge(base, package_annotation, by = "probeset_or_region_id", all.x = TRUE, sort = FALSE)
  if (nrow(table_annotation) > 0) {
    out <- merge(out, table_annotation, by = "probeset_or_region_id", all.x = TRUE, sort = FALSE, suffixes = c("", ".table"))
    for (name in c("chromosome", "start", "stop", "strand", "gene_symbol")) {
      table_name <- paste0(name, ".table")
      if (table_name %in% names(out)) {
        if (!name %in% names(out)) out[[name]] <- NA
        missing <- is.na(out[[name]]) | !nzchar(as.character(out[[name]]))
        out[[name]][missing] <- out[[table_name]][missing]
        out[[table_name]] <- NULL
      }
    }
  }
  for (name in c("chromosome", "start", "stop", "strand", "gene_symbol")) {
    if (!name %in% names(out)) out[[name]] <- if (name %in% c("start", "stop")) NA_integer_ else ""
  }
  out
}

selector_mask <- function(table, args) {
  mask <- rep(FALSE, nrow(table))
  if (length(args$probeset_ids) > 0) {
    mask <- mask | as.character(table$probeset_or_region_id) %in% args$probeset_ids
  }
  if (length(args$genes) > 0 && "gene_symbol" %in% names(table)) {
    gene_query <- toupper(args$genes)
    values <- toupper(as.character(table$gene_symbol))
    mask <- mask | values %in% gene_query
    for (gene in gene_query) {
      mask <- mask | grepl(paste0("(^|[^A-Z0-9])", gene, "([^A-Z0-9]|$)"), values)
    }
  }
  if (length(args$loci) > 0 && all(c("chromosome", "start", "stop") %in% names(table))) {
    for (locus_value in args$loci) {
      locus <- parse_locus(locus_value)
      mask <- mask | (chrom_key(table$chromosome) == chrom_key(locus$chromosome) & table$stop >= locus$start & table$start <= locus$stop)
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

run_limma <- function(expr_table, sample_table, contrasts, output_dir) {
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
    path <- file.path(output_dir, paste0("probeset_limma_", safe_token(contrast), ".tsv"))
    write_tsv(tab, path)
    paths <- c(paths, path)
  }
  paths
}

sample_table <- make_sample_table(args)
rownames(sample_table) <- sample_table$sample_id
write_tsv(sample_table, file.path(args$output, "sample_table.tsv"))
if (length(args$contrasts) == 0) {
  args$contrasts <- default_contrasts(sample_table)
}

cdf_name <- infer_cdf_name(args$cdf_package, args$cdf_name)
cel_files <- sample_table$cel_path
names(cel_files) <- sample_table$sample_id
message("Reading ", length(cel_files), " CEL file(s) with CDF name ", cdf_name)
raw <- affy::ReadAffy(filenames = cel_files, cdfname = cdf_name)
Biobase::sampleNames(raw) <- sample_table$sample_id
Biobase::pData(raw) <- sample_table[Biobase::sampleNames(raw), , drop = FALSE]
message("Running affy::rma()")
eset <- affy::rma(raw)

expr <- Biobase::exprs(eset)
expr_table <- data.frame(feature_id = rownames(expr), expr, check.names = FALSE)
rownames(expr_table) <- NULL
annotation <- merge_annotation(expr_table$feature_id, args)
feature_table <- merge(
  data.frame(feature_id = expr_table$feature_id, probeset_or_region_id = expr_table$feature_id, stringsAsFactors = FALSE),
  annotation,
  by = "probeset_or_region_id",
  all.x = TRUE,
  sort = FALSE
)

mask <- if (args$allow_all_features) rep(TRUE, nrow(feature_table)) else selector_mask(feature_table, args)
selected <- feature_table[mask, , drop = FALSE]
selected_expr <- expr_table[expr_table$feature_id %in% selected$feature_id, , drop = FALSE]
joined <- merge(selected, selected_expr, by = "feature_id", all.x = TRUE, sort = FALSE)
sample_ids <- intersect(sample_table$sample_id, names(joined))
region <- data.frame(
  chromosome = as.character(joined$chromosome),
  start = suppressWarnings(as.integer(joined$start)),
  stop = suppressWarnings(as.integer(joined$stop)),
  strand = as.character(joined$strand),
  probeset_or_region_id = as.character(joined$probeset_or_region_id),
  transcript_cluster_id = "",
  number_of_probes = "",
  gene_symbol = as.character(joined$gene_symbol),
  joined[, sample_ids, drop = FALSE],
  check.names = FALSE,
  stringsAsFactors = FALSE
)
region <- add_condition_summaries(region, sample_table)
region <- add_logfc_summaries(region, sample_table, args$contrasts)
region <- region[order(chrom_order(region$chromosome), region$start, region$stop, region$probeset_or_region_id), , drop = FALSE]

artifact_paths <- c(file.path(args$output, "sample_table.tsv"))
expr_path <- file.path(args$output, "probeset_expression_rma.tsv")
feature_path <- file.path(args$output, "probeset_feature_annotation.tsv")
region_path <- file.path(args$output, "region_intensity_chrom_order.csv")
write_tsv(expr_table, expr_path)
write_tsv(feature_table, feature_path)
write_csv(region, region_path)
artifact_paths <- c(artifact_paths, expr_path, feature_path, region_path)
artifact_paths <- c(artifact_paths, run_limma(expr_table, sample_table, args$contrasts, args$output))

session_path <- file.path(args$output, "sessionInfo.txt")
capture.output(utils::sessionInfo(), file = session_path)
artifact_paths <- c(artifact_paths, session_path)

manifest_path <- file.path(args$output, "normalized_feature_matrix_manifest.json")
manifest_lines <- c(
  "{",
  paste0("  \"schema\": ", json_string("gentle.probe_region_normalized_matrix_manifest.v1"), ","),
  paste0("  \"platform\": ", json_string(args$platform_name), ","),
  paste0("  \"cdf_package\": ", json_string(args$cdf_package), ","),
  paste0("  \"cdf_name\": ", json_string(cdf_name), ","),
  paste0("  \"annotation_package\": ", json_string(args$annotation_package), ","),
  paste0("  \"annotation_library\": ", json_string(args$annotation_library), ","),
  paste0("  \"coordinate_system\": ", json_string(args$coordinate_system), ","),
  paste0("  \"genome_build\": ", json_string(args$genome_build), ","),
  paste0("  \"normalization\": ", json_string(args$normalization), ","),
  paste0("  \"targets\": ", json_array(c("probeset")), ","),
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
  paste0("  \"backend\": ", json_string("r_affy_cdf"), ","),
  paste0("  \"command\": ", json_string(paste(commandArgs(FALSE), collapse = " ")), ","),
  paste0("  \"platform\": ", json_string(args$platform_name), ","),
  paste0("  \"cdf_package\": ", json_string(args$cdf_package), ","),
  paste0("  \"cdf_name\": ", json_string(cdf_name), ","),
  paste0("  \"annotation_package\": ", json_string(args$annotation_package), ","),
  paste0("  \"annotation_library\": ", json_string(args$annotation_library), ","),
  paste0("  \"coordinate_system\": ", json_string(args$coordinate_system), ","),
  paste0("  \"genome_build\": ", json_string(args$genome_build), ","),
  paste0("  \"normalization\": ", json_string(args$normalization), ","),
  paste0("  \"output_dir\": ", json_string(args$output), ","),
  paste0("  \"source_note\": ", json_string("Legacy 3' IVT/CDF helper output; CDF supplies probe grouping, while genomic coordinates require a local annotation package or coordinate table."), ","),
  paste0("  \"artifacts\": ", json_array(artifact_paths)),
  "}"
)
writeLines(provenance_lines, provenance_path, useBytes = TRUE)

message("Finished GENtle probe-region affy/CDF backend.")
message("Primary output: ", region_path)
