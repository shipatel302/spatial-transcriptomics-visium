#!/usr/bin/env Rscript
# R/run_rctd.R
# RCTD cell type deconvolution for 10x Visium
# Called from Python via rpy2 or standalone
# Author: Shivani Patel

suppressPackageStartupMessages({
  library(spacexr)   # RCTD
  library(Matrix)
  library(dplyr)
  library(ggplot2)
})


run_rctd_pipeline <- function(spaceranger_path, reference_path,
                               sample_id, outdir) {
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

  # ── Load Visium spatial data ───────────────────────────────────────────────
  counts <- load_visium_counts(spaceranger_path)
  coords <- load_visium_coords(spaceranger_path)

  # ── Build SpatialRNA object ────────────────────────────────────────────────
  puck <- SpatialRNA(coords, counts, require.int = FALSE)

  # ── Load single-cell reference ─────────────────────────────────────────────
  reference <- load_sc_reference(reference_path)

  # ── Run RCTD ──────────────────────────────────────────────────────────────
  message("  Running RCTD deconvolution...")
  myRCTD <- create.RCTD(puck, reference, max_cores = 4)
  myRCTD <- run.RCTD(myRCTD, doublet_mode = "full")

  # ── Extract proportions ────────────────────────────────────────────────────
  results      <- myRCTD@results
  norm_weights <- normalize_weights(results$weights)
  proportions  <- as.data.frame(as.matrix(norm_weights))

  # ── Save ──────────────────────────────────────────────────────────────────
  write.table(proportions,
    file.path(outdir, paste0(sample_id, "_rctd_proportions.tsv")),
    sep = "\t", quote = FALSE)

  # ── Plot spatial proportions ───────────────────────────────────────────────
  plot_rctd_spatial(proportions, coords, sample_id, outdir)

  return(proportions)
}


load_visium_counts <- function(path) {
  h5_file <- file.path(path, "filtered_feature_bc_matrix.h5")
  if (file.exists(h5_file)) {
    return(Seurat::Read10X_h5(h5_file))
  }
  # Try MEX format
  mex_path <- file.path(path, "filtered_feature_bc_matrix")
  if (dir.exists(mex_path)) {
    return(Seurat::Read10X(mex_path))
  }
  stop("Cannot find count matrix in: ", path)
}


load_visium_coords <- function(path) {
  coords_file <- file.path(path, "spatial", "tissue_positions_list.csv")
  if (!file.exists(coords_file)) {
    coords_file <- file.path(path, "spatial", "tissue_positions.csv")
  }
  coords <- read.csv(coords_file, header = FALSE)
  colnames(coords) <- c("barcode", "in_tissue", "row", "col", "y", "x")
  coords <- coords[coords$in_tissue == 1, ]
  rownames(coords) <- coords$barcode
  return(coords[, c("x", "y")])
}


load_sc_reference <- function(ref_path) {
  if (!file.exists(ref_path)) {
    message("  Reference not found — using demo reference")
    return(make_demo_reference())
  }
  # Load h5ad via reticulate/SeuratDisk
  tryCatch({
    library(SeuratDisk)
    Convert(ref_path, dest = "h5seurat", overwrite = TRUE)
    seurat_obj <- LoadH5Seurat(gsub("\\.h5ad$", ".h5seurat", ref_path))
    cell_types <- as.factor(seurat_obj$cell_type)
    counts     <- seurat_obj@assays$RNA@counts
    return(Reference(counts, cell_types))
  }, error = function(e) {
    message("  Reference load failed: ", e$message, " — using demo reference")
    return(make_demo_reference())
  })
}


make_demo_reference <- function() {
  set.seed(42)
  n_cells  <- 500
  n_genes  <- 2000
  cell_types <- factor(sample(
    c("Tumor", "T_cell_CD8", "T_cell_CD4", "Macrophage_M1",
      "Macrophage_M2", "CAF", "Endothelial", "B_cell"),
    n_cells, replace = TRUE
  ))
  counts <- matrix(
    rnbinom(n_cells * n_genes, mu = 2, size = 0.5),
    nrow = n_genes, ncol = n_cells,
    dimnames = list(paste0("GENE_", seq_len(n_genes)),
                    paste0("CELL_", seq_len(n_cells)))
  )
  Reference(Matrix::Matrix(counts, sparse = TRUE), cell_types)
}


plot_rctd_spatial <- function(proportions, coords, sample_id, outdir) {
  cell_types <- colnames(proportions)
  shared     <- intersect(rownames(proportions), rownames(coords))
  if (length(shared) == 0) return(invisible(NULL))

  plot_df <- cbind(
    coords[shared, ],
    proportions[shared, ]
  )

  pdf(file.path(outdir, paste0(sample_id, "_rctd_spatial.pdf")),
      width = 4 * min(4, length(cell_types)),
      height = 4 * ceiling(length(cell_types) / 4))

  par(mfrow = c(ceiling(length(cell_types) / 4), min(4, length(cell_types))),
      mar   = c(2, 2, 2, 2))

  for (ct in cell_types) {
    vals <- plot_df[[ct]]
    col_scale <- colorRampPalette(c("#f8f9fa", "#e76f51"))(100)
    breaks    <- seq(0, max(vals, na.rm = TRUE), length.out = 101)
    col_idx   <- findInterval(vals, breaks, rightmost.closed = TRUE)

    plot(plot_df$x, plot_df$y,
         col  = col_scale[pmax(1, pmin(100, col_idx))],
         pch  = 16, cex = 0.5,
         main = gsub("_", " ", ct),
         xlab = "", ylab = "", axes = FALSE)
  }
  dev.off()
}
