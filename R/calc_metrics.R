quasar_prop_metrics <- function(estimated_proportions, ground_truth_proportions, epsilon = 1e-12) {
  
 
  ground_truth_proportions <- as.matrix(ground_truth_proportions)
  estimated_proportions <- as.matrix(estimated_proportions)
  
  if (!is.numeric(ground_truth_proportions)) stop("Ground truth must be numeric")
  if (!is.numeric(estimated_proportions)) stop("Estimated proportions must be numeric")
  

  common_samples <- intersect(rownames(estimated_proportions), rownames(ground_truth_proportions))
  common_celltypes <- intersect(colnames(estimated_proportions), colnames(ground_truth_proportions))
  
  if (length(common_celltypes) == 0) {
    stop("No common cell types found.\n",
         "Estimated has: ", paste(colnames(estimated_proportions), collapse=", "), "\n",
         "Ground truth has: ", paste(colnames(ground_truth_proportions), collapse=", "))
  }
  if (length(common_samples) == 0) {
    stop("No common samples between estimated and ground truth proportions")
  }
  

  estimated_proportions <- estimated_proportions[common_samples, common_celltypes, drop = FALSE]
  ground_truth_proportions <- ground_truth_proportions[common_samples, common_celltypes, drop = FALSE]
  

  q_pm_check_input_type(estimated_proportions, colnames(ground_truth_proportions), 2)
  q_pm_check_input_type(ground_truth_proportions, colnames(estimated_proportions), 2)
  

  
  # RMSE per cell type
  cell_type_rmse <- sqrt(colMeans((estimated_proportions - ground_truth_proportions)^2))
  

  cell_type_mad <- apply(abs(estimated_proportions - ground_truth_proportions), 2,
                         function(x) median(abs(x - median(x))))
  

  mean_ground_truth <- colMeans(ground_truth_proportions)
  
  # Normalized MAE per cell type
  compute_nmae <- function(Y, X) {
    if (length(Y) != length(X)) stop("Vectors must be of equal length.")
    na_mask <- !(is.na(X) | is.na(Y))
    X <- X[na_mask]
    Y <- Y[na_mask]
    n <- length(X)
    if (n == 0) stop("No valid (non-NA) pairs remain for NMAE calculation.")
    range_X <- max(X) - min(X)
    if (range_X == 0) {
      if (all(X == Y)) {
        return(0)
      } else {
        return(mean(abs(Y - X))) # fallback when GT is constant
      }
    } else {
      return(mean(abs(Y - X) / range_X))
    }
  }
  
  cell_type_nmae <- sapply(common_celltypes, function(ct) {
    compute_nmae(estimated_proportions[, ct], ground_truth_proportions[, ct])
  })
  
  # Pearson correlation per cell type
  pearson_celltype_cor <- sapply(common_celltypes, function(ct) {
    est <- estimated_proportions[, ct]
    gt  <- ground_truth_proportions[, ct]
    if (sd(est) < epsilon) est <- est + rnorm(length(est), mean = 0, sd = epsilon)
    if (sd(gt)  < epsilon) gt  <- gt  + rnorm(length(gt),  mean = 0, sd = epsilon)
    cor(est, gt, method = "pearson")
  })
  
  # Spearman correlation per cell type
  spearman_celltype_cor <- sapply(common_celltypes, function(ct) {
    est <- estimated_proportions[, ct]
    gt  <- ground_truth_proportions[, ct]
    if (sd(est) < epsilon) est <- est + rnorm(length(est), mean = 0, sd = epsilon)
    if (sd(gt)  < epsilon) gt  <- gt  + rnorm(length(gt),  mean = 0, sd = epsilon)
    cor(est, gt, method = "spearman")
  })
  

  
  KL_div <- function(P, Q, eps = epsilon) {

    P <- P + eps
    Q <- Q + eps
    sum(P * log(P / Q))
  }
  
  JSD_func <- function(P, Q, eps = epsilon) {

    P <- P / sum(P)
    Q <- Q / sum(Q)
    M <- 0.5 * (P + Q)
    sqrt(0.5 * KL_div(P, M, eps) + 0.5 * KL_div(Q, M, eps))
  }
  
  # Per cell type JSD: distribution = proportions across samples for that cell type
  per_celltype_jsd <- sapply(common_celltypes, function(ct) {
    true_vec <- ground_truth_proportions[, ct]
    est_vec  <- estimated_proportions[, ct]
    JSD_func(true_vec, est_vec, epsilon)
  })
  
  # Per sample JSD: distribution = proportions across cell types for that sample
  per_sample_jsd <- sapply(common_samples, function(sm) {
    true_vec <- ground_truth_proportions[sm, ]
    est_vec  <- estimated_proportions[sm, ]
    JSD_func(true_vec, est_vec, epsilon)
  })
  
  message("\n[SUCCESS]: All metrics calculated!\n")
  
  return(list(
    cell_type_rmse = cell_type_rmse,
    cell_type_mad = cell_type_mad,
    cell_type_nmae = cell_type_nmae,
    mean_ground_truth_for_nmae = mean_ground_truth,
    pearson_celltype_cor = pearson_celltype_cor,
    spearman_celltype_cor = spearman_celltype_cor,
    per_celltype_jsd = per_celltype_jsd,
    per_sample_jsd = per_sample_jsd
  ))
}

#' internal: check if datatype is correct
#'
#' @param file file
#' @param columns object single or list specifying column
#' @param option 1 for character, 2 for numeric
#'
#' @author Sergej Ruff
#' @noRd
q_pm_check_input_type <- function(file, columns, option) {
  if (!option %in% c(1, 2)) {
    stop("Invalid option. Please choose 1 for character or 2 for numeric.")
  }
  expected_class <- if (option == 1) "character" else c("numeric", "integer")
  all_names <- colnames(file)
  for (col in columns) {
    if (!(col %in% all_names)) {
      stop("Column '", col, "' not found in file. Available column names: ", paste(all_names, collapse = ", "))
    }
    if (!inherits(file[, col], expected_class)) {
      stop("Error: Column '", col, "' must be of type ", expected_class, ".")
    }
  }
}