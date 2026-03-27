

to_numeric_df <- function(df) {
  mat <- as.matrix(df)
  storage.mode(mat) <- "numeric"
  out <- as.data.frame(mat, check.names = FALSE)
  rownames(out) <- rownames(df)
  colnames(out) <- colnames(df)
  out
}

average_by_celltype <- function(scRNA_df) {
  ct_unique <- unique(colnames(scRNA_df))

  averaged <- lapply(ct_unique, function(ct) {
    sub <- scRNA_df[, colnames(scRNA_df) == ct, drop = FALSE]
    sub_mat <- as.matrix(sub)
    storage.mode(sub_mat) <- "numeric"
    rowMeans(sub_mat, na.rm = TRUE)
  })

  X_ref <- do.call(cbind, averaged)

  if (is.null(dim(X_ref))) {
    X_ref <- matrix(X_ref, ncol = 1)
  }

  rownames(X_ref) <- rownames(scRNA_df)
  colnames(X_ref) <- ct_unique

  as.data.frame(X_ref, check.names = FALSE)
}



simulate_data <- function(
    scRNA_df,
    n_mixtures,
    n_cells_in_mix,
    ct_labels = NULL,
    n_genes = NULL,
    seed = 1,
    coverage_threshold = 5
) {
  # Input Handling
  if (!is.null(n_genes)) {
    if (nrow(scRNA_df) < n_genes) {
      stop(
        paste0(
          "Number of genes to select needs to be equal or lower than the amount of genes provided in scRNA_df. ",
          "Found ", nrow(scRNA_df), " genes in scRNA_df but n_genes is ", n_genes
        )
      )
    }
  }


  ct_unique <- unique(colnames(scRNA_df))
  n_celltypes <- length(ct_unique)
  n_examples <- ncol(scRNA_df)

 
  celltype_counts <- table(colnames(scRNA_df))
  for (celltype in names(celltype_counts)) {
    count <- as.integer(celltype_counts[[celltype]])

    if (count == 1) {
      warning(
        paste0(
          "Warning: Only one example provided for cell type '",
          celltype,
          "'. The mean will be the same as the single example."
        )
      )
    } else if (count <= coverage_threshold) {
      warning(
        paste0(
          "Warning: Low coverage for cell type '",
          celltype,
          "' (",
          count,
          " examples)."
        )
      )
    }
  }

  # Attach ct_labels to scRNA_df if provided
  if (!is.null(ct_labels)) {
    if (length(ct_labels) != ncol(scRNA_df)) {
      stop(
        paste0(
          "Length of ct_labels (", length(ct_labels),
          ") does not match number of samples (columns) in scRNA_df (",
          ncol(scRNA_df), ")."
        )
      )
    }
    colnames(scRNA_df) <- ct_labels
  }

  # Gene filter function
  gene_filter <- function(scRNA_df, n_genes = 1000) {
    X_ref <- average_by_celltype(scRNA_df)
    gene_variances <- apply(X_ref, 1, var)
    selected_genes <- names(sort(gene_variances, decreasing = TRUE))[seq_len(n_genes)]
    selected_genes
  }

  # Filter genes by variance across cell types if n_genes is specified
  if (!is.null(n_genes)) {
    gene_selection <- gene_filter(scRNA_df, n_genes)
    scRNA_df <- scRNA_df[gene_selection, , drop = FALSE]
  }

  # Creating X_ref as average profiles of the unique cell types
  X_ref <- average_by_celltype(scRNA_df)

  # Generate random samples all at once via random index choice
  set.seed(seed)
  bulk_indices <- matrix(
    sample.int(ncol(scRNA_df), size = n_mixtures * n_cells_in_mix, replace = TRUE),
    nrow = n_mixtures,
    byrow = TRUE
  )

  sc_mat <- as.matrix(scRNA_df)
  storage.mode(sc_mat) <- "numeric"

  # Create Y_mat
  Y_list <- lapply(seq_len(n_mixtures), function(i) {
    rowSums(sc_mat[, bulk_indices[i, ], drop = FALSE]) / n_cells_in_mix
  })

  Y_mat <- as.data.frame(do.call(cbind, Y_list), check.names = FALSE)
  rownames(Y_mat) <- rownames(scRNA_df)
  colnames(Y_mat) <- as.character(0:(n_mixtures - 1))

  # Create C_mat by counting occurrences of each cell type in the sampled indices
  C_mat <- matrix(
    0,
    nrow = length(unique(colnames(scRNA_df))),
    ncol = n_mixtures,
    dimnames = list(unique(colnames(scRNA_df)), as.character(0:(n_mixtures - 1)))
  )

  for (i in seq_len(n_mixtures)) {
    sampled_celltypes <- colnames(scRNA_df)[bulk_indices[i, ]]
    counts <- table(sampled_celltypes)
    C_mat[names(counts), i] <- as.numeric(counts)
  }

  C_mat <- C_mat / n_cells_in_mix
  C_mat <- C_mat[colnames(X_ref), , drop = FALSE]
  C_mat <- as.data.frame(C_mat, check.names = FALSE)

  list(X_ref = X_ref, Y_mat = Y_mat, C_mat = C_mat)
}