# PhyloResist AI: reproduce the study's multi-locus GWAS with mrMLM.
#
# Run from the project root in an R/RStudio terminal:
#   Rscript notebooks/05_mrmlm_gwas.R
#
# This script reads the original CSVs without modifying them. The author
# confirmed trait1 is the isolate-level least-squares mean EC50 across three
# replicates, measured in mM. No log transformation is applied. Some replicate
# curve estimates were excluded by the paper's authors before calculating the
# least-squares means; this script uses the supplied means as-is.
#
# The analysis settings below mirror analysis_2025.R and the preprint:
# character genotype calls, trait 1, the five specified mrMLM methods,
# SearchRadius = 20, and CriLOD = 3. This is a reproduction attempt, not an
# independent validation or a multi-pathogen-ready result.

options(stringsAsFactors = FALSE)

project_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
raw_dir <- file.path(project_root, "data", "raw")
output_root <- file.path(project_root, "results", "gwas", "mrmlm")
genotype_path <- file.path(raw_dir, "genotype.csv")
phenotype_paths <- c(
  propiconazole = file.path(raw_dir, "phenotype_propiconazole.csv"),
  tebuconazole = file.path(raw_dir, "phenotype_tebuconazole.csv")
)

required_paths <- c(genotype_path, unname(phenotype_paths))
missing_paths <- required_paths[!file.exists(required_paths)]
if (length(missing_paths) > 0) {
  stop("Missing required input file(s):\n", paste(missing_paths, collapse = "\n"))
}

if (!requireNamespace("mrMLM", quietly = TRUE)) {
  stop("The R package 'mrMLM' is required. Install it in your R environment, then rerun.")
}
suppressPackageStartupMessages(library(mrMLM))

normalize_id <- function(x) {
  x <- trimws(as.character(x))
  # Match numeric-looking IDs read by different tools as 11157 vs 11157.0.
  x <- sub("\\.0+$", "", x)
  x[is.na(x)] <- ""
  x
}

# Read genotype calls as text so A/C/G/T are kept categorical and isolate
# headers such as 00123 are not silently converted to numbers.
genotype <- read.csv(
  genotype_path,
  check.names = FALSE,
  colClasses = "character",
  na.strings = c("", "NA"),
  strip.white = TRUE
)
if (ncol(genotype) <= 3) {
  stop("Expected three variant metadata columns followed by isolate ID columns.")
}
metadata_columns <- names(genotype)[1:3]
isolate_ids <- normalize_id(names(genotype)[-(1:3)])
if (any(isolate_ids == "") || anyDuplicated(isolate_ids)) {
  stop("Genotype isolate columns contain blank or duplicated IDs.")
}
genotype_calls <- as.matrix(genotype[, -(1:3), drop = FALSE])
genotype_calls <- toupper(trimws(genotype_calls))
observed_calls <- unique(as.vector(genotype_calls))
observed_calls <- observed_calls[!is.na(observed_calls)]
unexpected_calls <- setdiff(observed_calls, c("A", "C", "G", "T"))
if (length(unexpected_calls) > 0) {
  stop("Unexpected genotype codes: ", paste(unexpected_calls, collapse = ", "))
}

# Confirm the phenotype ID and trait columns before passing the original file
# to mrMLM. The phenotype headers differ in capitalization between files.
for (fungicide in names(phenotype_paths)) {
  phenotype_path <- phenotype_paths[[fungicide]]
  phenotype <- read.csv(
    phenotype_path,
    check.names = FALSE,
    colClasses = "character",
    na.strings = c("", "NA"),
    strip.white = TRUE
  )
  header_lower <- tolower(trimws(names(phenotype)))
  id_candidates <- which(grepl("phenotype|isolate|sample|strain|(^|[^a-z])id([^a-z]|$)", header_lower))
  if (length(id_candidates) == 0) id_candidates <- 1L
  id_column <- id_candidates[[1]]
  trait_candidates <- which(header_lower == "trait1")
  if (length(trait_candidates) != 1) {
    stop("Expected exactly one trait1 column in ", basename(phenotype_path))
  }
  phenotype_ids <- normalize_id(phenotype[[id_column]])
  trait_values <- suppressWarnings(as.numeric(phenotype[[trait_candidates[[1]]]]))
  if (any(phenotype_ids == "") || anyDuplicated(phenotype_ids)) {
    stop("Blank or duplicate phenotype isolate IDs in ", basename(phenotype_path))
  }
  if (any(!is.finite(trait_values))) {
    stop("Missing or non-numeric trait1 values in ", basename(phenotype_path))
  }
  if (!setequal(isolate_ids, phenotype_ids)) {
    genotype_only <- setdiff(isolate_ids, phenotype_ids)
    phenotype_only <- setdiff(phenotype_ids, isolate_ids)
    stop(
      "Genotype/phenotype IDs do not match for ", fungicide,
      ". Genotype-only: ", paste(genotype_only, collapse = ", "),
      "; phenotype-only: ", paste(phenotype_only, collapse = ", ")
    )
  }
}

methods_to_run <- c("mrMLM", "ISIS EM-BLASSO", "FASTmrMLM", "pKWmEB", "FASTmrEMMA")

# Keep a small audit manifest with the verified phenotype interpretation and
# the exact settings used for this run.
dir.create(output_root, recursive = TRUE, showWarnings = FALSE)
manifest <- data.frame(
  fungicide = names(phenotype_paths),
  phenotype_file = unname(phenotype_paths),
  trait_column = "trait1",
  phenotype_definition = "Isolate least-squares mean EC50 across three replicates",
  unit = "mM (confirmed by study author)",
  log_transformed = FALSE,
  n_isolates = length(isolate_ids),
  n_variants = nrow(genotype),
  metadata_columns = paste(metadata_columns, collapse = ";"),
  genotype_format = "Cha",
  methods = paste(methods_to_run, collapse = ";"),
  SearchRadius = 20,
  CriLOD = 3,
  mrMLM_version = as.character(utils::packageVersion("mrMLM")),
  stringsAsFactors = FALSE
)
write.csv(manifest, file.path(output_root, "input_manifest.csv"), row.names = FALSE)

for (fungicide in names(phenotype_paths)) {
  phenotype_path <- phenotype_paths[[fungicide]]
  result_dir <- file.path(output_root, fungicide)
  dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)

  message("Starting mrMLM GWAS for ", fungicide, "...")
  message("Using ", length(isolate_ids), " isolates, ", nrow(genotype), " variants, trait1 in mM.")

  # Pass the source CSV unchanged. In this file format, trait=1 selects the
  # first phenotype column after the isolate ID, named trait1.
  mrMLM::mrMLM(
    fileGen = normalizePath(genotype_path, winslash = "/", mustWork = TRUE),
    filePhe = normalizePath(phenotype_path, winslash = "/", mustWork = TRUE),
    fileKin = NULL,
    filePS = NULL,
    dir = normalizePath(result_dir, winslash = "/", mustWork = TRUE),
    Genformat = "Cha",
    method = methods_to_run,
    trait = 1,
    SearchRadius = 20,
    CriLOD = 3,
    DrawPlot = TRUE,
    Plotformat = "pdf"
  )

  message("Finished ", fungicide, ". Results written under: ", result_dir)
}

message("Both fungicide GWAS runs completed. Input manifest: ", file.path(output_root, "input_manifest.csv"))
