validate_and_extract_sc_data <- function(sc_object) {
  if (!inherits(sc_object, "SingleCellExperiment") && !inherits(sc_object, "Seurat")) {
    stop("Single-cell reference must be SingleCellExperiment or Seurat object")
  }

  if (inherits(sc_object, "SingleCellExperiment")) {
    list(
      expression = SingleCellExperiment::counts(sc_object),
      metadata = SummarizedExperiment::colData(sc_object)
    )
  } else {
    list(
      expression = sc_object[["RNA"]]$counts,
      metadata = sc_object@meta.data
    )
  }
}

sim_bulk_seurat_sce <- function(
    sc_object,
    celltype_col,
    n_mixtures,
    n_cells_in_mix,
    n_genes = NULL,
    seed = 1,
    coverage_threshold = 5
) {
  extracted <- validate_and_extract_sc_data(sc_object)

  expr <- extracted$expression
  meta <- as.data.frame(extracted$metadata)

  if (!(celltype_col %in% colnames(meta))) {
    stop(sprintf("'%s' not found in metadata.", celltype_col))
  }

  if (ncol(expr) != nrow(meta)) {
    stop("Number of columns in expression matrix must match number of rows in metadata.")
  }

  ct_labels <- as.character(meta[[celltype_col]])

  if (any(is.na(ct_labels))) {
    stop("Cell type labels contain NA values.")
  }

  if (!is.null(n_genes) && nrow(expr) < n_genes) {
    stop(
      paste0(
        "Number of genes to select needs to be equal or lower than the amount of genes provided. ",
        "Found ", nrow(expr), " genes but n_genes is ", n_genes, "."
      )
    )
  }

  celltype_counts <- table(ct_labels)
  for (celltype in names(celltype_counts)) {
    count <- as.integer(celltype_counts[[celltype]])
    if (count == 1) {
      warning(
        sprintf(
          "Warning: Only one example provided for cell type '%s'. The mean will be the same as the single example.",
          celltype
        )
      )
    } else if (count <= coverage_threshold) {
      warning(
        sprintf(
          "Warning: Low coverage for cell type '%s' (%d examples).",
          celltype, count
        )
      )
    }
  }

  average_by_celltype <- function(expr_mat, labels) {
    celltypes <- unique(labels)

    avg_list <- lapply(celltypes, function(ct) {
      cols <- which(labels == ct)
      if (length(cols) == 1) {
        as.numeric(expr_mat[, cols, drop = TRUE])
      } else {
        Matrix::rowMeans(expr_mat[, cols, drop = FALSE])
      }
    })

    X_ref <- do.call(cbind, avg_list)
    rownames(X_ref) <- rownames(expr_mat)
    colnames(X_ref) <- celltypes
    as.data.frame(X_ref, check.names = FALSE)
  }

  gene_filter <- function(expr_mat, labels, n_genes) {
    X_ref <- average_by_celltype(expr_mat, labels)
    gene_variances <- apply(X_ref, 1, stats::var)
    names(sort(gene_variances, decreasing = TRUE))[seq_len(n_genes)]
  }

  if (!is.null(n_genes)) {
    selected_genes <- gene_filter(expr, ct_labels, n_genes)
    expr <- expr[selected_genes, , drop = FALSE]
  }

  X_ref <- average_by_celltype(expr, ct_labels)

  set.seed(seed)
  bulk_indices <- matrix(
    sample.int(ncol(expr), size = n_mixtures * n_cells_in_mix, replace = TRUE),
    nrow = n_mixtures,
    byrow = TRUE
  )

  Y_list <- lapply(seq_len(n_mixtures), function(i) {
    Matrix::rowSums(expr[, bulk_indices[i, ], drop = FALSE]) / n_cells_in_mix
  })

  Y_mat <- as.data.frame(do.call(cbind, Y_list), check.names = FALSE)
  rownames(Y_mat) <- rownames(expr)
  colnames(Y_mat) <- as.character(seq_len(n_mixtures))

  celltypes_unique <- unique(ct_labels)
  C_mat <- matrix(
    0,
    nrow = length(celltypes_unique),
    ncol = n_mixtures,
    dimnames = list(celltypes_unique, as.character(seq_len(n_mixtures)))
  )

  for (i in seq_len(n_mixtures)) {
    sampled_celltypes <- ct_labels[bulk_indices[i, ]]
    sampled_counts <- table(sampled_celltypes)
    C_mat[names(sampled_counts), i] <- as.numeric(sampled_counts)
  }

  C_mat <- C_mat / n_cells_in_mix
  C_mat <- C_mat[colnames(X_ref), , drop = FALSE]
  C_mat <- as.data.frame(C_mat, check.names = FALSE)

  list(
    X_ref = X_ref,
    Y_mat = Y_mat,
    C_mat = C_mat
  )
}