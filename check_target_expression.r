# ============================================================
# Script: check_target_expression_uniprot.R
# Purpose: Check if target genes/proteins (partial matches allowed)
#          are expressed in a Seurat object or 10X dataset.
#          Reports expression counts, percentages, and thresholds.
# ============================================================

suppressPackageStartupMessages({
  library(optparse)
  library(Seurat)
  library(dplyr)
  library(httr)
  library(jsonlite)
  library(readr)
  library(Matrix)
})

# ------------------------------------------------------------
# ---- Utility Functions ----
# ------------------------------------------------------------

# Load Seurat object or 10X dataset
load_10x_data <- function(input_path) {
  if (!dir.exists(input_path)) {
    stop("❌ Input path not found or not a directory. Only 10X directories are allowed: ", input_path)
  }

  message("📂 Input detected as directory. Reading 10X data...")
  tryCatch(
    {
      counts <- Read10X(data.dir = input_path)
      seurat_obj <- CreateSeuratObject(counts = counts, project = "Project")
    },
    error = function(e) stop("❌ Failed to read 10X data: ", e$message)
  )
  seurat_obj
}

# Map gene symbol to EnsemblPlants GeneID via UniProt (robust)
map_gene_to_ensemblplants <- function(gene_symbol) {
  query <- URLencode(gene_symbol)
  url <- paste0(
    "https://rest.uniprot.org/uniprotkb/search?query=", query,
    "&fields=xref_ensemblplants&format=json&size=1"
  )

  message("🔍 Querying: ", url)
  resp <- tryCatch(GET(url, add_headers(Accept = "application/json")),
                   error = function(e) { warning("HTTP request failed: ", e$message); return(NULL) })
  if (is.null(resp) || status_code(resp) != 200) {
    warning("❌ UniProt query failed for ", gene_symbol)
    return(NA)
  }

  data <- content(resp, as = "text", encoding = "UTF-8")
  json_data <- tryCatch(fromJSON(data, flatten = TRUE), error = function(e) NULL)
  if (is.null(json_data) || length(json_data$results) == 0) {
    warning("⚠️ No UniProt mapping found for ", gene_symbol)
    return(NA)
  }

  refs <- json_data$results$uniProtKBCrossReferences[[1]]
  if (is.null(refs) || !"database" %in% names(refs)) {
    warning("⚠️ No valid cross-reference data for ", gene_symbol)
    return(NA)
  }

  # Safely subset for EnsemblPlants entries, ignoring NAs
  ensembl_refs <- refs[!is.na(refs$database) & refs$database == "EnsemblPlants", ]
  if (nrow(ensembl_refs) == 0) {
    warning("⚠️ No EnsemblPlants cross-reference found for ", gene_symbol)
    return(NA)
  }

  gene_id <- NA
  for (props in ensembl_refs$properties) {
    if (!is.null(props$key) && !is.null(props$value)) {
      match <- props$value[props$key == "GeneId"]
      if (length(match) > 0) {
        gene_id <- match[1]
        break
      }
    }
  }

  if (is.na(gene_id)) {
    warning("⚠️ No GeneId property found for ", gene_symbol)
  } else {
    gene_id <- sub("\\..*$", "", gene_id) # strip version suffix
  }

  gene_id
}

# Check expression for genes, returning counts and percentages
check_expression <- function(seurat_obj, genes) {
  assay_name <- DefaultAssay(seurat_obj)
  expr_data <- GetAssayData(seurat_obj, assay = assay_name, layer = "counts")[genes, , drop = FALSE]
  total_cells <- ncol(expr_data)

  expressed_counts <- rowSums(expr_data > 0)
  expressed_perc <- (expressed_counts / total_cells) * 100

  data.frame(
    GeneID = names(expressed_counts),
    ExpressedCells = as.integer(expressed_counts),
    TotalCells = total_cells,
    PercentExpressed = round(expressed_perc, 2),
    stringsAsFactors = FALSE
  )
}

# Strip version suffix like ".3" from input query gene symbols
strip_version_suffix <- function(gene_symbol) {
  sub("\\..*$", "", gene_symbol)
}

# ------------------------------------------------------------
# ---- CLI Parsing ----
# ------------------------------------------------------------

option_list <- list(
  make_option(c("-i", "--input"),
    type = "character", default = NULL,
    help = "Path to Seurat RDS file or 10X directory"
  ),
  make_option(c("-g", "--genes"),
    type = "character", default = NULL,
    help = "Comma-separated gene symbols (e.g., HY5,COP1,SPA1)"
  ),
  make_option(c("-p", "--percent"),
    type = "numeric", default = 10,
    help = "Percentage threshold (e.g., 10 means at least 10%% of cells must express the gene)"
  ),
  make_option(c("-t", "--tsv"),
    type = "character", default = NULL,
    help = "Optional: output results to TSV file (e.g., results.tsv)"
  )
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# ---- Argument validation ----
if (is.null(opt$input)) {
  print_help(opt_parser)
  stop("❌ Provide input path with -i")
}
if (is.null(opt$genes)) {
  print_help(opt_parser)
  stop("❌ Provide comma-separated gene symbols with -g")
}

# ------------------------------------------------------------
# ---- Main Script ----
# ------------------------------------------------------------

main <- function(opt) {
  target_symbols <- strsplit(opt$genes, ",")[[1]] %>% trimws()
  target_symbols <- sapply(target_symbols, strip_version_suffix)
  message("🧬 Querying genes (version-stripped): ", paste(target_symbols, collapse = ", "))
  message("📊 Expression threshold: ", opt$percent, "%")

  # Map gene symbols to EnsemblPlants GeneIDs
  mapped_genes <- sapply(target_symbols, map_gene_to_ensemblplants)
  names(mapped_genes) <- target_symbols

  if (any(is.na(mapped_genes))) {
    warning("⚠️ Could not map: ", paste(names(mapped_genes)[is.na(mapped_genes)], collapse = ", "))
  }

  solyc_ids <- na.omit(mapped_genes)
  if (length(solyc_ids) == 0) stop("❌ No genes could be mapped to EnsemblPlants IDs.")

  # Load Seurat object
  seurat_obj <- load_10x_data(opt$input)
  all_genes <- rownames(seurat_obj)

  # ---- Fuzzy matching ----
  present_genes <- unlist(lapply(solyc_ids, function(id) {
    matches <- grep(paste0("^", id, "(\\.|$)"), all_genes, value = TRUE)
    if (length(matches) > 0) matches[1] else NA
  }))
  names(present_genes) <- solyc_ids
  present_genes <- na.omit(present_genes)

  missing_genes <- solyc_ids[!solyc_ids %in% names(present_genes)]
  if (length(missing_genes) > 0) {
    message(
      "ℹ️ Some mapped genes were not found in the dataset: ",
      paste(missing_genes, collapse = ", ")
    )
  }

  results <- data.frame(
    Symbol = character(),
    GeneID = character(),
    ExpressedCells = numeric(),
    PercentExpressed = numeric(),
    Status = character(),
    stringsAsFactors = FALSE
  )

  if (length(present_genes) > 0) {
    expr_stats <- check_expression(seurat_obj, unname(present_genes))

    for (symbol in names(mapped_genes)) {
      gene_id <- mapped_genes[symbol]

      if (is.na(gene_id)) {
        results <- rbind(results, data.frame(
          Symbol = symbol, GeneID = NA,
          ExpressedCells = NA, PercentExpressed = NA,
          Status = "MAPPING FAILED"
        ))
      } else if (!(gene_id %in% names(present_genes))) {
        results <- rbind(results, data.frame(
          Symbol = symbol, GeneID = gene_id,
          ExpressedCells = 0, PercentExpressed = 0,
          Status = "NOT IN DATASET"
        ))
      } else {
        matched_name <- present_genes[[gene_id]]
        row <- expr_stats[expr_stats$GeneID == matched_name, ]
        status <- ifelse(row$PercentExpressed >= opt$percent,
          paste0("EXPRESSED ≥ ", opt$percent, "%"),
          paste0("NOT EXPRESSED (< ", opt$percent, "%)")
        )

        results <- rbind(results, data.frame(
          Symbol = symbol,
          GeneID = matched_name,
          ExpressedCells = row$ExpressedCells,
          PercentExpressed = row$PercentExpressed,
          Status = status,
          stringsAsFactors = FALSE
        ))
      }
    }

    # Print results
    apply(results, 1, function(row) {
      cat(sprintf(
        " %s (%s): %s [%.2f%% of cells]\n",
        row["Symbol"], row["GeneID"], row["Status"],
        as.numeric(row["PercentExpressed"])
      ))
    })

    # Optional TSV export
    if (!is.null(opt$tsv)) {
      write_tsv(results, opt$tsv)
      message("💾 Results written to ", opt$tsv)
    }
  } else {
    cat("❌ None of the mapped genes are present in the dataset.\n")
  }
}

# Run main with error handling
tryCatch(
  {
    main(opt)
  },
  error = function(e) {
    message(e$message)
    quit(status = 1)
  }
)
