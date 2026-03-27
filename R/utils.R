


#'@noRd
calculate_estimated_composition <- function(X_mat, Y_mat, gamma) {
 
  if (nrow(X_mat) != nrow(Y_mat)) {
    stop(
      paste0(
        "Number of genes (rows) must match between X_mat and Y_mat. ",
        "Found ", nrow(X_mat), " in X_mat and ", nrow(Y_mat), " in Y_mat."
      )
    )
  }

  if (nrow(X_mat) != nrow(gamma)) {
    stop(
      paste0(
        "Number of genes (rows) must match between X_mat and gamma. ",
        "Found ", nrow(X_mat), " in X_mat and ", nrow(gamma), " in gamma."
      )
    )
  }

  if (!identical(rownames(X_mat), rownames(Y_mat))) {
    stop("Gene labels (index) must be identical between X_mat and Y_mat.")
  }

  if (!identical(rownames(X_mat), rownames(gamma))) {
    stop("Gene labels (index) must be identical between X_mat and gamma.")
  }


  X <- as.matrix(X_mat)
  Y <- as.matrix(Y_mat)
  gamma_mat <- as.matrix(gamma)

  storage.mode(X) <- "numeric"
  storage.mode(Y) <- "numeric"
  storage.mode(gamma_mat) <- "numeric"

  gamma_vec <- as.numeric(gamma_mat[, 1])
  Gamma <- diag(gamma_vec)

  C_e <- solve(t(X) %*% Gamma %*% X) %*% t(X) %*% Gamma %*% Y
  C_e[C_e < 0] <- 0

  C_estimated <- as.data.frame(C_e, check.names = FALSE)
  rownames(C_estimated) <- colnames(X_mat)
  colnames(C_estimated) <- colnames(Y_mat)

  C_estimated
}



#'@noRd
XYCtoTorch <- function(X_mat, Y_mat, C_mat) {


  X_torch <- torch::torch_tensor(as.matrix(X_mat))
  Y_torch <- torch::torch_tensor(as.matrix(Y_mat))
  C_torch <- torch::torch_tensor(as.matrix(C_mat))

  list(X_torch = X_torch, Y_torch = Y_torch, C_torch = C_torch)
}


#' @noRd
plot_corr <- function(C_true, C_est, title = "", color = "grey", hidden_ct = NULL, c_est = NULL, path = NULL) {
  correlations <- numeric()

  rows <- ceiling(length(rownames(C_true)) / 3)
  empty <- rows * 3 - length(rownames(C_true))

  if (!is.null(path)) {
    grDevices::pdf(paste0(path, ".pdf"))
    on.exit(grDevices::dev.off(), add = TRUE)
  }

  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar), add = TRUE)

  par(mfrow = c(rows, 3))
  title(main = title)

  for (celltype in rownames(C_est)) {
    r_val <- suppressWarnings(
      cor(
        as.numeric(C_true[celltype, ]),
        as.numeric(C_est[celltype, ]),
        method = "spearman",
        use = "complete.obs"
      )
    )

    correlations[celltype] <- r_val

    plot(
      as.numeric(C_true[celltype, ]),
      as.numeric(C_est[celltype, ]),
      main = paste0(celltype, " (R: ", sprintf("%1.2f", r_val), ")"),
      xlab = "true abundance",
      ylab = "estimated abundance",
      col = color,
      pch = 16
    )
  }

  if (!is.null(hidden_ct)) {
    r_val <- suppressWarnings(
      cor(
        as.numeric(C_true[hidden_ct, ]),
        as.numeric(c_est["hidden", ]),
        method = "spearman",
        use = "complete.obs"
      )
    )

    correlations["hidden"] <- r_val

    plot(
      as.numeric(C_true[hidden_ct, ]),
      as.numeric(c_est["hidden", ]),
      main = paste0("Hidden (R: ", sprintf("%1.2f", r_val), ")"),
      xlab = "true abundance",
      ylab = "estimated abundance",
      col = "black",
      pch = 16
    )
  }


  if (empty == 2) {
    plot.new()
    plot.new()
  } else if (empty == 1) {
    plot.new()
  }

  correlations
}


#' @noRd
calculate_corr <- function(C_true, C_est, hidden_ct = NULL, c_est = NULL) {
  correlations <- numeric()

  for (celltype in rownames(C_est)) {
    r_val <- suppressWarnings(
      cor(
        as.numeric(C_true[celltype, ]),
        as.numeric(C_est[celltype, ]),
        method = "spearman",
        use = "complete.obs"
      )
    )
    correlations[celltype] <- r_val
  }

  if (!is.null(hidden_ct)) {
    r_val <- suppressWarnings(
      cor(
        as.numeric(C_true[hidden_ct, ]),
        as.numeric(c_est["hidden", ]),
        method = "spearman",
        use = "complete.obs"
      )
    )
    correlations["hidden"] <- r_val
    return(correlations)
  }


  invisible(NULL)
}
