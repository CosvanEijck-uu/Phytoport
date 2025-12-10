# ============================================================
# Script: check_target_expression_uniprot_integrated.R
# Purpose: Combine UniProt → EnsemblPlants gene mapping with
#          expression quantification, dynamic thresholds,
#          histogram plotting (Catppuccin colors), raw count export,
#          and summary TSV output.
# ============================================================

suppressPackageStartupMessages({
  library(optparse)
  library(Seurat)
  library(dplyr)
  library(httr)
  library(jsonlite)
  library(readr)
  library(Matrix)
  library(ggplot2)
  library(gridExtra)
})

# ------------------------------------------------------------
# Load 10X dataset and create Seurat object
# ------------------------------------------------------------
load_10x_data <- function(input_path) {
  if (!dir.exists(input_path)) {
    stop("❌ Input path not found or not a directory: ", input_path)
  }
  message("📂 Reading 10X data...")
  counts <- Read10X(data.dir = input_path)
  CreateSeuratObject(counts = counts, project = basename(input_path))
}

# ------------------------------------------------------------
# Strip version suffix
# ------------------------------------------------------------
strip_version_suffix <- function(x) sub("\\..*$", "", x)

# ------------------------------------------------------------
# UniProt → EnsemblPlants gene mapping
# ------------------------------------------------------------
map_gene_to_ensemblplants <- function(symbol) {
  query <- URLencode(symbol)
  url <- paste0(
    "https://rest.uniprot.org/uniprotkb/search?query=", query,
    "&fields=xref_ensemblplants&format=json&size=1"
  )
  message("🔍 Querying UniProt for: ", symbol)

  resp <- tryCatch(GET(url, add_headers(Accept = "application/json")), error = function(e) NULL)
  if (is.null(resp) || status_code(resp) != 200) {
    return(NA)
  }

  json_data <- tryCatch(fromJSON(content(resp, as = "text", encoding = "UTF-8"), flatten = TRUE), error = function(e) NULL)
  if (is.null(json_data) || length(json_data$results) == 0) {
    return(NA)
  }

  refs <- json_data$results$uniProtKBCrossReferences[[1]]
  if (is.null(refs)) {
    return(NA)
  }

  ensembl_refs <- refs[refs$database == "EnsemblPlants", ]
  if (nrow(ensembl_refs) == 0) {
    return(NA)
  }

  gene_id <- NA
  for (props in ensembl_refs$properties) {
    if (!is.null(props$key) && !is.null(props$value)) {
      match <- props$value[props$key == "GeneId"]
      if (length(match) > 0) {
        gene_id <- strip_version_suffix(match[1])
        break
      }
    }
  }

  gene_id
}

# ------------------------------------------------------------
# Compute expression statistics with dynamic threshold
# ------------------------------------------------------------
check_expression <- function(seurat_obj, genes, min_counts = 1) {
  expr <- GetAssayData(seurat_obj, assay = DefaultAssay(seurat_obj), layer = "counts")[genes, , drop = FALSE]
  total <- ncol(expr)

  expressed_cells <- rowSums(expr >= min_counts)
  expressed_perc <- (expressed_cells / total) * 100

  avg_expr <- apply(expr, 1, function(x) {
    pos <- x[x >= min_counts]
    if (length(pos) == 0) {
      return(0)
    }
    mean(pos)
  })

  list(
    stats = data.frame(
      GeneID = rownames(expr),
      ExpressedCells = expressed_cells,
      TotalCells = total,
      PercentExpressed = round(expressed_perc, 2),
      AvgExpression = avg_expr,
      Expressed = expressed_cells > 0,
      stringsAsFactors = FALSE
    ),
    counts = expr
  )
}

# ------------------------------------------------------------
# Histogram
# ------------------------------------------------------------
plot_expression_histogram <- function(gene_id, vector) {
  vector <- as.numeric(vector)
  nz <- vector[vector > 0]

  if (length(nz) == 0) {
    return(
      ggplot() +
        theme_minimal() +
        labs(title = paste0(gene_id, " — no non-zero counts"))
    )
  }

  # Viridis palette hex codes (standard viridis)
  viridis_hex <- c(
    "#440154", "#482475", "#414487", "#355F8D", "#2A788E",
    "#21918C", "#22A884", "#44BF70", "#7AD151", "#BDDF26", "#FDE725"
  )

  df <- data.frame(Expression = nz)

  ggplot(df, aes(x = Expression, fill = after_stat(count))) +
    geom_histogram(binwidth = 1, color = "white") +
    scale_fill_gradientn(colors = viridis_hex) +
    theme_minimal(base_size = 12) +
    theme(
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),
      text = element_text(color = "black"),
      axis.text = element_text(color = "black")
    ) +
    labs(title = paste0(gene_id), fill = "Count")
}

# ------------------------------------------------------------
# CLI options
# ------------------------------------------------------------
option_list <- list(
  make_option(c("-i", "--input"), type = "character", help = "10X directory"),
  make_option(c("-g", "--genes"), type = "character", help = "Comma-separated UniProt symbols"),
  make_option(c("-t", "--tsv"), type = "character", default = NULL, help = "Write summary TSV"),
  make_option(c("-m", "--min_counts"), type = "integer", default = 1, help = "Minimum transcripts to count as expressed")
)
opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$input) || is.null(opt$genes)) stop("❌ Provide -i and -g")

# ------------------------------------------------------------
# Main
# ------------------------------------------------------------
main <- function(opt) {
  symbols <- strsplit(opt$genes, ",")[[1]] %>%
    trimws() %>%
    sapply(strip_version_suffix)
  mapped <- sapply(symbols, map_gene_to_ensemblplants)
  names(mapped) <- symbols

  message("🧬 Mapped IDs: ")
  print(mapped)

  valid_ids <- na.omit(mapped)
  if (length(valid_ids) == 0) stop("❌ No valid EnsemblPlants mappings found.")

  seurat_obj <- load_10x_data(opt$input)
  all_genes <- rownames(seurat_obj)

  present <- unlist(lapply(valid_ids, function(id) {
    m <- grep(paste0("^", id, "(\\.|$)"), all_genes, value = TRUE)
    if (length(m) > 0) m[1] else NA
  }))
  names(present) <- valid_ids
  present <- na.omit(present)

  expr <- check_expression(seurat_obj, unname(present), min_counts = opt$min_counts)
  expr_stats <- expr$stats
  expr_data <- expr$counts

  # ------------------------------------------------------------
  # histogram
  # ------------------------------------------------------------
  pdf("expression_histograms.pdf", width = 8, height = 6)
  plots <- lapply(rownames(expr_data), function(g) plot_expression_histogram(g, expr_data[g, ]))
  do.call(grid.arrange, c(plots, ncol = 2))
  dev.off()

  # ------------------------------------------------------------
  # TSV export
  # ------------------------------------------------------------
  if (!is.null(opt$tsv)) {
    write_tsv(expr_stats, opt$tsv)
    write_tsv(as.data.frame(expr_data, row.names = "GeneID"), sub(".tsv$", "_raw_counts.tsv", opt$tsv))
  }
  message("✅ Finished dataset → count distribution histogram PDF saved")
}

main(opt)
