# ==============================================================================
# BEER: Batch EffEct Remover for Single-Cell Data
# ==============================================================================
# Version: 0.2.0 (Seurat v5 compatible)
# Author: Feng Zhang (original), Refactored for v5
# Date: October 14, 2025
# License: MIT
#
# Description:
#   BEER removes batch effects from single-cell RNA-seq and ATAC-seq data
#   by detecting batch-specific principal components using mutual nearest
#   neighbors and selecting biologically relevant PCs for downstream analysis.
#
# Reference:
#   Zhang, F., Wu, Y., & Tian, W. (2019). A novel approach to remove the
#   batch effect of single-cell data. Cell Discovery, 5, 46.
#   https://doi.org/10.1038/s41421-019-0114-x
#
# Usage:
#   source('BEER_v5.R')
#   mybeer <- BEER(DATA, BATCH, GNUM=30, PCNUM=50, ROUND=1)
# ==============================================================================

# Global Settings --------------------------------------------------------------
BEER_VERSION <- "0.2.0"
CORMETHOD <- "spearman"

# Package Dependencies ---------------------------------------------------------
#' @description Load required packages with version checking
#' @importFrom Seurat CreateSeuratObject NormalizeData FindVariableFeatures
#' @importFrom sva ComBat
#' @importFrom limma makeContrasts
.load_dependencies <- function() {
  required_packages <- c("Seurat", "sva", "limma", "stringi")

  for (pkg in required_packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop(sprintf("Package '%s' is required but not installed. Install with:\n  install.packages('%s')",
                   pkg, pkg))
    }
  }

  suppressPackageStartupMessages({
    library(Seurat)
    library(sva)
    library(limma)
    library(stringi)
  })

  # Check Seurat version
  seurat_version <- packageVersion("Seurat")
  if (seurat_version < "5.0.0") {
    warning(sprintf("BEER v%s is optimized for Seurat v5+. Current version: %s",
                    BEER_VERSION, seurat_version))
  }

  invisible(TRUE)
}

# Utility Functions: Matrix Operations ----------------------------------------

#' Convert sparse matrix to dense matrix
#'
#' @param mat A sparse matrix (dgCMatrix or similar)
#' @return A dense matrix
#' @details Handles large sparse matrices efficiently
as_matrix <- function(mat) {
  if (!inherits(mat, "Matrix")) {
    return(as.matrix(mat))
  }

  tmp <- matrix(data = 0L, nrow = mat@Dim[1], ncol = mat@Dim[2])
  row_pos <- mat@i + 1
  col_pos <- findInterval(seq(mat@x) - 1, mat@p[-1]) + 1
  val <- mat@x

  for (i in seq_along(val)) {
    tmp[row_pos[i], col_pos[i]] <- val[i]
  }

  rownames(tmp) <- mat@Dimnames[[1]]
  colnames(tmp) <- mat@Dimnames[[2]]

  return(tmp)
}

#' Remove duplicate gene names from expression matrix
#'
#' @param mat Expression matrix with genes as rows
#' @return Matrix with unique gene names only
.check_duplicate_genes <- function(mat) {
  mat <- as.matrix(mat)
  gene_names <- rownames(mat)

  unique_genes <- names(which(table(gene_names) == 1))
  mat <- mat[rownames(mat) %in% unique_genes, , drop = FALSE]

  n_removed <- length(gene_names) - nrow(mat)
  if (n_removed > 0) {
    message(sprintf("Removed %d duplicate gene names", n_removed))
  }

  return(mat)
}

#' Get Seurat assay data (v5 compatible)
#'
#' @param object Seurat object
#' @param assay Assay name (default: "RNA")
#' @param layer Layer name (default: "data")
#' @return Expression matrix
.get_assay_data <- function(object, assay = "RNA", layer = "data") {
  # Try v5 method first
  tryCatch({
    return(LayerData(object, assay = assay, layer = layer))
  }, error = function(e) {
    # Fall back to v4 method
    slot_name <- switch(layer,
                        "counts" = "counts",
                        "data" = "data",
                        "scale.data" = "scale.data",
                        layer)
    return(GetAssayData(object, assay = assay, slot = slot_name))
  })
}

#' Set Seurat assay data (v5 compatible)
#'
#' @param object Seurat object
#' @param layer Layer name
#' @param new_data New expression matrix
#' @param assay Assay name (default: "RNA")
#' @return Updated Seurat object
.set_assay_data <- function(object, layer, new_data, assay = "RNA") {
  # Try v5 method first
  tryCatch({
    object <- SetAssayData(object, assay = assay, layer = layer, new.data = new_data)
    return(object)
  }, error = function(e) {
    # Fall back to v4 method
    slot_name <- switch(layer,
                        "counts" = "counts",
                        "data" = "data",
                        "scale.data" = "scale.data",
                        layer)
    object <- SetAssayData(object, assay = assay, slot = slot_name, new.data = new_data)
    return(object)
  })
}

# Utility Functions: Data Combination -----------------------------------------

#' Combine two expression matrices
#'
#' @param exp_sc_mat1 First expression matrix
#' @param exp_sc_mat2 Second expression matrix
#' @param fill Keep genes expressed in only one condition (default: FALSE)
#' @return List with combined matrix and individual matrices
#' @export
.simple_combine <- function(exp_sc_mat1, exp_sc_mat2, fill = FALSE) {
  exp_sc_mat <- exp_sc_mat1
  exp_ref_mat <- exp_sc_mat2

  # Optional: fill missing genes with zeros
  if (fill) {
    gene1 <- rownames(exp_sc_mat)
    gene2 <- rownames(exp_ref_mat)
    gene12 <- setdiff(gene2, gene1)
    gene21 <- setdiff(gene1, gene2)

    if (length(gene12) > 0) {
      exp_sc_mat_add <- matrix(0, ncol = ncol(exp_sc_mat), nrow = length(gene12))
      rownames(exp_sc_mat_add) <- gene12
      colnames(exp_sc_mat_add) <- colnames(exp_sc_mat)
      exp_sc_mat <- rbind(exp_sc_mat, exp_sc_mat_add)
    }

    if (length(gene21) > 0) {
      exp_ref_mat_add <- matrix(0, ncol = ncol(exp_ref_mat), nrow = length(gene21))
      rownames(exp_ref_mat_add) <- gene21
      colnames(exp_ref_mat_add) <- colnames(exp_ref_mat)
      exp_ref_mat <- rbind(exp_ref_mat, exp_ref_mat_add)
    }
  }

  # Sort and find overlapping genes
  exp_sc_mat <- exp_sc_mat[order(rownames(exp_sc_mat)), , drop = FALSE]
  exp_ref_mat <- exp_ref_mat[order(rownames(exp_ref_mat)), , drop = FALSE]

  gene_overlap <- intersect(rownames(exp_sc_mat), rownames(exp_ref_mat))

  if (length(gene_overlap) == 0) {
    stop("No overlapping genes found between matrices")
  }

  exp_sc_mat <- exp_sc_mat[gene_overlap, , drop = FALSE]
  exp_ref_mat <- exp_ref_mat[gene_overlap, , drop = FALSE]

  list(
    exp_sc_mat1 = exp_sc_mat,
    exp_sc_mat2 = exp_ref_mat,
    combine = cbind(exp_sc_mat, exp_ref_mat)
  )
}

# Utility Functions: Aggregation ----------------------------------------------

#' Generate aggregate expression by group (sum)
#'
#' @param exp_sc_mat Expression matrix
#' @param tag Group labels (character vector)
#' @param print_step Print progress every N groups (default: 100)
#' @return Aggregated expression matrix
.generate_agg <- function(exp_sc_mat, tag, print_step = 100) {
  tag <- as.character(tag)
  refnames <- unique(tag)
  total_num <- length(refnames)

  new_ref <- matrix(0, ncol = total_num, nrow = nrow(exp_sc_mat))
  outnames <- character(total_num)

  for (i in seq_along(refnames)) {
    one <- refnames[i]
    this_col <- which(tag == one)
    outnames[i] <- one

    if (length(this_col) > 1) {
      new_ref[, i] <- rowSums(exp_sc_mat[, this_col, drop = FALSE])
    } else {
      new_ref[, i] <- exp_sc_mat[, this_col]
    }

    if (i %% print_step == 1) {
      message(sprintf("Processing group %d / %d", i, total_num))
    }
  }

  rownames(new_ref) <- rownames(exp_sc_mat)
  colnames(new_ref) <- outnames

  # Handle single group case
  if (ncol(new_ref) == 1) {
    new_ref <- cbind(new_ref[, 1], new_ref[, 1])
    colnames(new_ref) <- c(outnames, outnames)
  }

  return(new_ref)
}

#' Generate aggregate expression by group (mean)
#'
#' @param exp_sc_mat Expression matrix
#' @param tag Group labels
#' @param print_step Print progress every N groups
#' @return Aggregated expression matrix
.generate_mean <- function(exp_sc_mat, tag, print_step = 100) {
  tag <- as.character(tag)
  refnames <- unique(tag)
  total_num <- length(refnames)

  new_ref <- matrix(0, ncol = total_num, nrow = nrow(exp_sc_mat))
  outnames <- character(total_num)

  for (i in seq_along(refnames)) {
    one <- refnames[i]
    this_col <- which(tag == one)
    outnames[i] <- one

    if (length(this_col) > 1) {
      new_ref[, i] <- rowMeans(exp_sc_mat[, this_col, drop = FALSE])
    } else {
      new_ref[, i] <- exp_sc_mat[, this_col]
    }

    if (i %% print_step == 1) {
      message(sprintf("Processing group %d / %d", i, total_num))
    }
  }

  rownames(new_ref) <- rownames(exp_sc_mat)
  colnames(new_ref) <- outnames

  if (ncol(new_ref) == 1) {
    new_ref <- cbind(new_ref[, 1], new_ref[, 1])
    colnames(new_ref) <- c(outnames, outnames)
  }

  return(new_ref)
}

# Core Functions: Grouping and MN Pairs ---------------------------------------

#' Group cells within each batch using K-means
#'
#' @param x Cell embeddings (e.g., UMAP coordinates)
#' @param tag Batch label
#' @param gnum Number of groups
#' @return Group labels for cells
.get_group <- function(x, tag, gnum) {
  message(sprintf("Creating %d groups for batch: %s", gnum, tag))

  clust <- kmeans(x, centers = gnum, iter.max = 100, nstart = 10)$cluster
  group <- paste0(tag, "_", as.character(clust))

  message(sprintf("Generated %d groups", length(unique(group))))

  return(group)
}

#' Find mutual nearest neighbor (MN) pairs across batches
#'
#' @param pbmc Seurat object with group metadata
#' @param round Number of MN detection rounds
#' @return Matrix of MN pairs
.get_vp_all <- function(pbmc, round) {
  message("Finding mutual nearest neighbor (MN) pairs...")

  # Get aggregated expression
  group_data <- .get_assay_data(pbmc, layer = "data")
  ref <- .generate_agg(group_data, pbmc@meta.data$group)

  # Calculate correlation matrix
  cvref <- cor(ref, method = CORMETHOD)
  ubatch <- unique(pbmc@meta.data$batch)

  # Extract batch from group name
  .get_batch <- function(x) {
    strsplit(x, "_")[[1]][1]
  }
  group_batch <- sapply(colnames(cvref), .get_batch)

  # Find MN pairs in correlation matrix
  .get_mn <- function(cor_mat) {
    vp <- NULL
    for (i in seq_len(nrow(cor_mat))) {
      for (j in seq_len(ncol(cor_mat))) {
        this_cor <- cor_mat[i, j]
        if (this_cor != -99999 &&
            this_cor == max(cor_mat[i, ]) &&
            this_cor == max(cor_mat[, j])) {
          vp <- cbind(vp, c(rownames(cor_mat)[i], colnames(cor_mat)[j]))
        }
      }
    }
    return(vp)
  }

  # Iteratively find MN pairs
  vp <- NULL
  for (round_idx in seq_len(round)) {
    # Mask already found pairs
    if (!is.null(vp)) {
      for (i in seq_len(ncol(vp))) {
        p1 <- vp[1, i]
        p2 <- vp[2, i]
        b1 <- .get_batch(p1)
        b2 <- .get_batch(p2)
        b1_idx <- which(group_batch == b1)
        b2_idx <- which(group_batch == b2)

        cvref[p1, b2_idx] <- -99999
        cvref[p2, b1_idx] <- -99999
        cvref[b2_idx, p1] <- -99999
        cvref[b1_idx, p2] <- -99999
      }
    }

    # Find new pairs between all batch combinations
    for (i in seq_len(length(ubatch) - 1)) {
      for (j in (i + 1):length(ubatch)) {
        b1 <- ubatch[i]
        b2 <- ubatch[j]
        b1_idx <- which(group_batch == b1)
        b2_idx <- which(group_batch == b2)

        this_cor_mat <- cvref[b1_idx, b2_idx, drop = FALSE]
        this_vp <- .get_mn(this_cor_mat)
        vp <- cbind(vp, this_vp)
      }
    }

    message(sprintf("Round %d completed", round_idx))
  }

  vp <- t(unique(t(vp)))
  message(sprintf("Found %d MN pairs", nrow(vp)))

  return(vp)
}

# Core Functions: PC Evaluation ------------------------------------------------

#' Evaluate principal components for batch correlation
#'
#' @param dr Dimension reduction matrix (PCA embeddings)
#' @param group Cell group labels
#' @param vp Valid pairs (MN pairs)
#' @return List with correlation statistics for each PC
.evaluate_pro_beer <- function(dr, group, vp) {
  message("Evaluating principal components for batch effects...")
  message("This may take a few minutes...")

  valid_pair <- vp
  all_cor <- numeric(ncol(dr))
  all_pv <- numeric(ncol(dr))
  all_lcor <- numeric(ncol(dr))
  all_lpv <- numeric(ncol(dr))
  all_lc1 <- numeric(ncol(dr))
  all_lc2 <- numeric(ncol(dr))

  for (this_dr in seq_len(ncol(dr))) {
    lst1_quantile <- numeric()
    lst2_quantile <- numeric()

    for (i in seq_len(nrow(valid_pair))) {
      this_pair <- valid_pair[i, ]
      idx1 <- which(group == this_pair[1])
      idx2 <- which(group == this_pair[2])

      lst1_quantile <- c(lst1_quantile, quantile(dr[idx1, this_dr]))
      lst2_quantile <- c(lst2_quantile, quantile(dr[idx2, this_dr]))
    }

    # Spearman correlation
    test1 <- cor.test(lst1_quantile, lst2_quantile, method = CORMETHOD)
    all_cor[this_dr] <- test1$estimate
    all_pv[this_dr] <- test1$p.value

    # Pearson correlation
    test2 <- cor.test(lst1_quantile, lst2_quantile, method = "pearson")
    all_lcor[this_dr] <- test2$estimate
    all_lpv[this_dr] <- test2$p.value

    # Linear regression
    fit <- lm(lst1_quantile ~ lst2_quantile)
    all_lc1[this_dr] <- summary(fit)$coefficients[1, 4]
    all_lc2[this_dr] <- summary(fit)$coefficients[2, 4]

    if (this_dr %% 10 == 0) {
      message(sprintf("Evaluated PC %d / %d", this_dr, ncol(dr)))
    }
  }

  message("PC evaluation completed!")

  list(
    cor = all_cor,
    pv = all_pv,
    fdr = p.adjust(all_pv, method = "fdr"),
    lc1 = all_lc1,
    lc2 = all_lc2,
    lcor = all_lcor,
    lpv = all_lpv,
    lfdr = p.adjust(all_lpv, method = "fdr")
  )
}

# Main Function: BEER ----------------------------------------------------------

#' BEER: Batch EffEct Remover for single-cell data
#'
#' @param data Expression matrix (genes x cells)
#' @param batch Character vector of batch labels
#' @param gnum Number of groups per batch (default: 30)
#' @param pcnum Number of PCs to compute (default: 50)
#' @param gn Number of variable genes (default: 2000)
#' @param combat Use ComBat correction (default: TRUE)
#' @param seed Random seed (default: 123)
#' @param n_components Number of UMAP components (default: 2)
#' @param round Number of MN detection rounds (default: 1)
#' @param remove_genes Genes to remove from analysis (default: NULL)
#' @return List containing Seurat object and batch effect statistics
#' @export
#'
#' @examples
#' \dontrun{
#' mybeer <- BEER(DATA, BATCH, gnum=30, pcnum=50, round=1)
#' pbmc <- mybeer$seurat
#' pbmc <- RunUMAP(pbmc, dims=mybeer$select)
#' }
BEER <- function(data,
                 batch,
                 gnum = 30,
                 pcnum = 50,
                 gn = 2000,
                 combat = TRUE,
                 seed = 123,
                 n_components = 2,
                 round = 1,
                 remove_genes = NULL) {

  # Initialize
  set.seed(seed)
  .load_dependencies()

  message("================================================================================")
  message(sprintf("BEER v%s - Batch EffEct Remover", BEER_VERSION))
  message("================================================================================")
  message(sprintf("Start time: %s", Sys.time()))

  # Validate inputs
  if (ncol(data) != length(batch)) {
    stop("Number of cells in data must match length of batch vector")
  }

  # Sanitize batch names (replace dots with underscores)
  batch <- stringi::stri_replace_all(batch, "_", fixed = ".")
  ubatch <- unique(batch)

  message(sprintf("Number of batches: %d", length(ubatch)))
  message(sprintf("Number of cells: %d", ncol(data)))
  message(sprintf("Number of genes: %d", nrow(data)))
  message(sprintf("Groups per batch (GNUM): %d", gnum))
  message(sprintf("Variable genes per batch (GN): %d", gn))
  message(sprintf("MN detection rounds: %d", round))

  # Find variable genes per batch
  message("\n--- Step 1: Finding variable genes ---")
  varg <- character()

  for (idx in seq_along(ubatch)) {
    this_batch <- ubatch[idx]
    message(sprintf("Processing batch %d/%d: %s", idx, length(ubatch), this_batch))

    batch_cells <- which(batch == this_batch)
    this_pbmc <- CreateSeuratObject(
      counts = data[, batch_cells, drop = FALSE],
      min.cells = 0,
      min.features = 0,
      project = this_batch
    )

    this_pbmc <- NormalizeData(
      this_pbmc,
      normalization.method = "LogNormalize",
      scale.factor = 10000,
      verbose = FALSE
    )

    this_pbmc <- FindVariableFeatures(
      this_pbmc,
      selection.method = "vst",
      nfeatures = gn,
      verbose = FALSE
    )

    varg <- c(varg, VariableFeatures(this_pbmc))
  }

  varg <- unique(varg)
  message(sprintf("Total variable genes: %d", length(varg)))

  # Create combined Seurat object
  message("\n--- Step 2: Creating Seurat object ---")
  pbmc <- CreateSeuratObject(
    counts = data,
    min.cells = 0,
    min.features = 0,
    project = "BEER"
  )
  pbmc@meta.data$batch <- batch
  VariableFeatures(pbmc) <- varg

  # Remove specified genes
  if (!is.null(remove_genes)) {
    n_remove <- sum(VariableFeatures(pbmc) %in% remove_genes)
    message(sprintf("Removing %d specified genes", n_remove))
    VariableFeatures(pbmc) <- setdiff(VariableFeatures(pbmc), remove_genes)
    message(sprintf("Genes after removal: %d", length(VariableFeatures(pbmc))))
  }

  # Normalize
  pbmc <- NormalizeData(
    pbmc,
    normalization.method = "LogNormalize",
    scale.factor = 10000,
    verbose = FALSE
  )

  # Optional: ComBat correction before scaling
  combat_exp <- NULL
  if (combat) {
    message("\n--- Step 3: ComBat batch correction ---")

    pheno <- data.frame(batch = batch)
    orig_data <- .get_assay_data(pbmc, layer = "data")
    used_gene_idx <- which(rownames(orig_data) %in% varg)
    edata <- as_matrix(orig_data)[used_gene_idx, , drop = FALSE]

    modcombat <- model.matrix(~ 1, data = pheno)
    combat_edata <- ComBat(
      dat = edata,
      batch = pheno$batch,
      mod = modcombat,
      par.prior = TRUE,
      prior.plots = FALSE
    )

    rownames(combat_edata) <- rownames(edata)
    colnames(combat_edata) <- colnames(edata)
    combat_edata <- as.matrix(combat_edata)
    combat_edata[combat_edata < 0] <- 0
    combat_edata[is.na(combat_edata)] <- 0

    pbmc <- .set_assay_data(pbmc, layer = "data", new_data = combat_edata)
    pbmc <- ScaleData(pbmc, features = VariableFeatures(pbmc), verbose = FALSE)
    pbmc <- .set_assay_data(pbmc, layer = "data", new_data = orig_data)

    combat_exp <- combat_edata
    rm(edata, combat_edata, orig_data)
    gc()
  } else {
    pbmc <- ScaleData(pbmc, features = VariableFeatures(pbmc), verbose = FALSE)
  }

  # Calculate PCs and UMAP
  message("\n--- Step 4: Calculating PCA ---")
  pbmc <- RunPCA(
    pbmc,
    features = VariableFeatures(pbmc),
    npcs = pcnum,
    seed.use = seed,
    verbose = FALSE
  )

  pbmc <- RunUMAP(
    pbmc,
    dims = 1:pcnum,
    seed.use = seed,
    n.components = n_components,
    verbose = FALSE
  )

  # Group cells within batches
  message("\n--- Step 5: Grouping cells within batches ---")
  dr <- pbmc@reductions$umap@cell.embeddings
  group <- character(length(batch))

  for (this_batch in ubatch) {
    batch_idx <- which(batch == this_batch)
    this_dr <- dr[batch_idx, , drop = FALSE]

    # Adjust gnum if needed
    this_gnum <- min(gnum, length(batch_idx) - 1)
    if (this_gnum != gnum) {
      message(sprintf("Adjusted GNUM to %d for batch %s (n=%d cells)",
                      this_gnum, this_batch, length(batch_idx)))
    }

    this_group <- .get_group(this_dr, this_batch, this_gnum)
    group[batch_idx] <- this_group
  }

  pbmc@meta.data$group <- group

  # Find MN pairs
  message("\n--- Step 6: Finding mutual nearest neighbors ---")
  vp <- .get_vp_all(pbmc, round)

  # Evaluate PCs
  message("\n--- Step 7: Evaluating PCs for batch correlation ---")
  pca_embeddings <- pbmc@reductions$pca@cell.embeddings
  out <- .evaluate_pro_beer(pca_embeddings, group, vp)

  # Mark cells in MN pairs
  map <- rep("NA", length(group))
  map[group %in% vp[, 1]] <- "V1"
  map[group %in% vp[, 2]] <- "V2"
  pbmc@meta.data$map <- map

  # Select PCs with low batch correlation
  pcuse <- which(
    (rank(out$cor) >= length(out$cor) / 2 | out$cor > 0.7) &
    (rank(out$lcor) >= length(out$lcor) / 2 | out$lcor > 0.7)
  )

  message(sprintf("\nSelected %d / %d PCs with low batch correlation",
                  length(pcuse), pcnum))

  # Prepare results
  result <- list(
    seurat = pbmc,
    vp = vp,
    cor = out$cor,
    pv = out$pv,
    fdr = out$fdr,
    lcor = out$lcor,
    lpv = out$lpv,
    lc1 = out$lc1,
    lc2 = out$lc2,
    lfdr = out$lfdr,
    select = pcuse,
    round = round,
    combat = combat,
    combat_exp = combat_exp,
    remove_genes = remove_genes,
    gnum = gnum,
    gn = gn,
    pcnum = pcnum,
    seed = seed,
    n_components = n_components,
    version = BEER_VERSION,
    app = "BEER"
  )

  message("================================================================================")
  message("BEER completed successfully!")
  message(sprintf("End time: %s", Sys.time()))
  message("================================================================================")

  return(result)
}

# Alias for backwards compatibility
MBEER <- BEER

# Helper Functions: PC Selection ----------------------------------------------

#' Select PCs based on correlation cutoffs
#'
#' @param result BEER result object
#' @param cutr Rank correlation cutoff (default: 0.7)
#' @param cutl Linear correlation cutoff (default: 0.7)
#' @return Vector of selected PC indices
#' @export
.get_use <- function(result, cutr = 0.7, cutl = 0.7) {
  pcuse <- which(result$cor > cutr & result$lcor > cutl)
  return(pcuse)
}

#' Select PCs based on rank and correlation thresholds
#'
#' @param result BEER result object
#' @param cutr Rank correlation cutoff
#' @param cutl Linear correlation cutoff
#' @param rr Rank ratio (default: 0.5)
#' @param rl Linear rank ratio (default: 0.5)
#' @return Vector of selected PC indices
#' @export
.select_use <- function(result, cutr = 0.7, cutl = 0.7, rr = 0.5, rl = 0.5) {
  n_pcs <- length(result$cor)
  pcuse <- which(
    (rank(result$cor) >= n_pcs * rr | result$cor > cutr) &
    (rank(result$lcor) >= n_pcs * rl | result$lcor > cutl)
  )
  return(pcuse)
}

# Main Function: ReBEER --------------------------------------------------------

#' Re-run BEER with adjusted parameters (faster than full BEER)
#'
#' @param mybeer Previous BEER result object
#' @param gnum Number of groups per batch
#' @param pcnum Number of PCs
#' @param seed Random seed
#' @param n_components Number of UMAP components
#' @param round Number of MN detection rounds
#' @param remove_genes Genes to remove
#' @return Updated BEER result object
#' @export
ReBEER <- function(mybeer,
                   gnum = 30,
                   pcnum = 50,
                   seed = 123,
                   n_components = 2,
                   round = 1,
                   remove_genes = NULL) {

  set.seed(seed)
  .load_dependencies()

  message("================================================================================")
  message(sprintf("ReBEER v%s - Re-running with adjusted parameters", BEER_VERSION))
  message("================================================================================")
  message(sprintf("Start time: %s", Sys.time()))

  # Extract data from previous run
  data <- .get_assay_data(mybeer$seurat, layer = "counts")
  batch <- mybeer$seurat@meta.data$batch
  ubatch <- unique(batch)

  pbmc <- mybeer$seurat

  message(sprintf("Groups per batch (GNUM): %d", gnum))
  message(sprintf("Total variable genes: %d", length(VariableFeatures(pbmc))))
  message(sprintf("MN detection rounds: %d", round))

  # Remove specified genes
  if (!is.null(remove_genes)) {
    n_remove <- sum(VariableFeatures(pbmc) %in% remove_genes)
    message(sprintf("Removing %d specified genes", n_remove))
    VariableFeatures(pbmc) <- setdiff(VariableFeatures(pbmc), remove_genes)
    message(sprintf("Genes after removal: %d", length(VariableFeatures(pbmc))))
  }

  # Recalculate PCA and UMAP
  message("Recalculating PCA...")
  pbmc <- RunPCA(
    pbmc,
    features = VariableFeatures(pbmc),
    npcs = pcnum,
    seed.use = seed,
    verbose = FALSE
  )

  pbmc <- RunUMAP(
    pbmc,
    dims = 1:pcnum,
    seed.use = seed,
    n.components = n_components,
    verbose = FALSE
  )

  # Re-group cells
  message("Re-grouping cells...")
  dr <- pbmc@reductions$umap@cell.embeddings
  group <- character(length(batch))

  for (this_batch in ubatch) {
    batch_idx <- which(batch == this_batch)
    this_dr <- dr[batch_idx, , drop = FALSE]
    this_gnum <- min(gnum, length(batch_idx) - 1)
    this_group <- .get_group(this_dr, this_batch, this_gnum)
    group[batch_idx] <- this_group
  }

  pbmc@meta.data$group <- group

  # Find MN pairs
  message("Finding MN pairs...")
  vp <- .get_vp_all(pbmc, round)

  # Evaluate PCs
  message("Evaluating PCs...")
  pca_embeddings <- pbmc@reductions$pca@cell.embeddings
  out <- .evaluate_pro_beer(pca_embeddings, group, vp)

  # Mark cells
  map <- rep("NA", length(group))
  map[group %in% vp[, 1]] <- "V1"
  map[group %in% vp[, 2]] <- "V2"
  pbmc@meta.data$map <- map

  # Select PCs
  pcuse <- which(
    (rank(out$cor) >= length(out$cor) / 2 | out$cor > 0.7) &
    (rank(out$lcor) >= length(out$lcor) / 2 | out$lcor > 0.7)
  )

  message(sprintf("Selected %d / %d PCs", length(pcuse), pcnum))

  # Prepare results
  result <- list(
    seurat = pbmc,
    vp = vp,
    cor = out$cor,
    pv = out$pv,
    fdr = out$fdr,
    lcor = out$lcor,
    lpv = out$lpv,
    lc1 = out$lc1,
    lc2 = out$lc2,
    lfdr = out$lfdr,
    select = pcuse,
    round = round,
    gnum = gnum,
    pcnum = pcnum,
    seed = seed,
    n_components = n_components,
    version = BEER_VERSION,
    app = "ReBEER"
  )

  message("================================================================================")
  message("ReBEER completed!")
  message(sprintf("End time: %s", Sys.time()))
  message("================================================================================")

  return(result)
}

# Enhancement Functions --------------------------------------------------------

#' Apply ComBat correction to PCA embeddings
#'
#' @param pbmc Seurat object
#' @return Seurat object with ComBat-corrected PCA
#' @export
BEER.combat <- function(pbmc) {
  message("Applying ComBat correction to PCA embeddings...")

  batch <- as.character(pbmc@meta.data$batch)
  pca <- pbmc@reductions$pca@cell.embeddings

  pheno <- data.frame(batch = batch)
  edata <- t(pca)
  modcombat <- model.matrix(~ 1, data = pheno)

  combat_edata <- ComBat(
    dat = edata,
    batch = pheno$batch,
    mod = modcombat,
    par.prior = TRUE,
    prior.plots = FALSE
  )

  corrected_pca <- t(combat_edata)
  colnames(corrected_pca) <- colnames(pca)
  rownames(corrected_pca) <- rownames(pca)

  pbmc@reductions$pca@cell.embeddings <- corrected_pca

  message("ComBat correction completed")
  return(pbmc)
}

#' Apply BBKNN enhancement for better batch mixing
#'
#' @param pbmc Seurat object
#' @param pcuse Selected PCs to use
#' @param nb Neighbors within batch (default: 3)
#' @param nt Number of trees for Annoy (default: 10)
#' @param dm Number of UMAP dimensions (default: 2)
#' @return UMAP coordinates
#' @export
#' @importFrom reticulate import py_to_r
BEER.bbknn <- function(pbmc, pcuse, nb = 3, nt = 10, dm = 2) {
  message("Applying BBKNN enhancement...")
  message("This requires Python packages: scanpy, anndata, bbknn")

  if (!requireNamespace("reticulate", quietly = TRUE)) {
    stop("Package 'reticulate' is required for BBKNN. Install with: install.packages('reticulate')")
  }

  library(reticulate)

  # Import Python modules
  tryCatch({
    anndata <- import("anndata", convert = FALSE)
    bbknn <- import("bbknn", convert = FALSE)
    sc <- import("scanpy", convert = FALSE)
  }, error = function(e) {
    stop("Failed to import Python modules. Ensure scanpy, anndata, and bbknn are installed:\n  pip install scanpy anndata bbknn")
  })

  batch <- as.character(pbmc@meta.data$batch)
  pca_all <- pbmc@reductions$pca@cell.embeddings
  pca_use <- pbmc@reductions$pca@cell.embeddings[, pcuse, drop = FALSE]

  # Create AnnData object
  adata <- anndata$AnnData(X = pca_all, obs = batch)
  pcnum <- ncol(pca_use)

  # Run BBKNN
  sc$tl$pca(adata, n_comps = as.integer(pcnum))
  adata$obsm$X_pca <- pca_use

  bbknn$bbknn(
    adata,
    batch_key = 0L,
    neighbors_within_batch = as.integer(nb),
    n_pcs = as.integer(pcnum),
    annoy_n_trees = as.integer(nt)
  )

  sc$tl$umap(adata, n_components = as.integer(dm))

  # Extract UMAP
  umap <- py_to_r(adata$obsm['X_umap'])
  rownames(umap) <- rownames(pbmc@reductions$umap@cell.embeddings)
  colnames(umap) <- paste0("UMAP_", seq_len(ncol(umap)))

  message("BBKNN enhancement completed")
  return(umap)
}

# Utility Functions: I/O -------------------------------------------------------

#' Read table with proper formatting
#'
#' @param path File path
#' @param sep Separator (default: tab)
#' @param upper Convert gene names to uppercase (default: FALSE)
#' @return Expression matrix
#' @export
.read_table <- function(path, sep = "\t", upper = FALSE) {
  data <- read.table(
    file = path,
    sep = sep,
    header = TRUE,
    row.names = NULL
  )
  data <- apply(data, 2, as.character)

  if (upper) {
    data[, 1] <- toupper(data[, 1])
  }

  # Keep only unique genes
  tab <- table(data[, 1])
  uniq <- names(tab)[tab == 1]
  data <- data[data[, 1] %in% uniq, , drop = FALSE]

  rn <- data[, 1]
  data <- data[, -1, drop = FALSE]
  data <- apply(data, 2, as.numeric)
  rownames(data) <- rn

  return(data)
}

#' Write table with proper formatting
#'
#' @param data Data matrix
#' @param path Output file path
#' @param sep Separator (default: tab)
#' @param title First column title (default: "OUTPUT")
#' @export
.write_table <- function(data, path, sep = "\t", title = "OUTPUT") {
  out <- cbind(rownames(data), data)
  colnames(out)[1] <- title

  write.table(
    out,
    file = path,
    sep = sep,
    row.names = FALSE,
    col.names = TRUE,
    quote = FALSE
  )

  message(sprintf("Data written to: %s", path))
}

#' Count positive values in vector
#'
#' @param x Numeric vector
#' @return Number of positive values
.get_pos <- function(x) {
  sum(x > 0)
}

# Startup Message --------------------------------------------------------------
message(sprintf("BEER v%s loaded successfully!", BEER_VERSION))
message("Compatible with Seurat v5+")
message("Source: https://github.com/jumphone/BEER")
message("Citation: Zhang et al. (2019) Cell Discovery, 5, 46")
message("")
message("Quick start:")
message("  mybeer <- BEER(DATA, BATCH, gnum=30, pcnum=50, round=1)")
message("  pbmc <- mybeer$seurat")
message("  pbmc <- RunUMAP(pbmc, dims=mybeer$select)")
message("================================================================================")
