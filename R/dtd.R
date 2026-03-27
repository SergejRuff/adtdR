


XYCtoTorch <- function(X_mat, Y_mat, C_mat) {
  X_torch <- torch::torch_tensor(as.matrix(X_mat), dtype = torch_float())
  Y_torch <- torch::torch_tensor(as.matrix(Y_mat), dtype = torch_float())
  C_torch <- torch::torch_tensor(as.matrix(C_mat), dtype = torch_float())
  list(X = X_torch, Y = Y_torch, C = C_torch)
}

C_est <- function(X, Y, gamma) {
  Gamma <- torch::torch_diag(gamma)
  A <- torch::torch_mm(torch::torch_mm(torch::torch_t(X), Gamma), X)
  B <- torch::torch_mm(torch::torch_mm(torch::torch_t(X), Gamma), Y)
  C_e <- torch::linalg_solve(A, B)
  C_e <- torch::torch_clamp(C_e, min = 0)
  C_e
}

DTD_loss <- function(C, C_e) {
  C_mean <- torch::torch_mean(C, dim = 2, keepdim = TRUE)
  C_e_mean <- torch::torch_mean(C_e, dim = 2, keepdim = TRUE)

  r_num <- torch::torch_sum((C - C_mean) * (C_e - C_e_mean), dim = 2)
  r_den <- torch::torch_sqrt(
    torch::torch_sum((C - C_mean)^2, dim = 2) *
      torch::torch_sum((C_e - C_e_mean)^2, dim = 2)
  )

  r <- r_num / r_den

  list(
    loss = -torch::torch_sum(r),
    r = r
  )
}

Cosine_loss <- function(C, C_e) {
  loss <- torch::torch_tensor(0, dtype = C$dtype, device = C$device)

  for (i in seq_len(C$size(1))) {
    r_num <- C[i, ]$dot(C_e[i, ])
    r_den <- torch::torch_sqrt(torch::torch_sum(C[i, ]^2)) * torch::torch_sqrt(torch::torch_sum(C_e[i, ]^2))
    loss <- loss + r_num / r_den
  }

  C_mean <- torch::torch_mean(C, dim = 2, keepdim = TRUE)
  C_e_mean <- torch::torch_mean(C_e, dim = 2, keepdim = TRUE)

  r_num <- torch::torch_sum((C - C_mean) * (C_e - C_e_mean), dim = 2)
  r_den <- torch::torch_sqrt(
    torch::torch_sum((C - C_mean)^2, dim = 2) *
      torch::torch_sum((C_e - C_e_mean)^2, dim = 2)
  )

  r <- r_num / r_den

  list(
    loss = -loss,
    r = r
  )
}

get_loss_function <- function(func) {
  if (func == "pearson") {
    return(DTD_loss)
  }
  if (func == "cosine") {
    return(Cosine_loss)
  }
  stop("func must be either 'pearson' or 'cosine'")
}



dtd_run <- function(
    X_mat,
    Y_mat,
    C_mat,
    iterations = 1000,
    plot = FALSE,
    path_plot = NULL,
    func = "pearson"
) {


  if (nrow(X_mat) != nrow(Y_mat)) {
    stop(
      paste0(
        "Number of genes (rows) must match between X_mat and Y_mat. ",
        "Found ", nrow(X_mat), " in X_mat and ", nrow(Y_mat), " in Y_mat."
      )
    )
  }

  if (ncol(X_mat) != nrow(C_mat)) {
    stop(
      paste0(
        "Number of celltypes (columns) in X_mat must match the number of celltypes (rows) in C_mat. ",
        "Found ", ncol(X_mat), " in X_mat and ", nrow(C_mat), " in C_mat."
      )
    )
  }

  if (ncol(Y_mat) != ncol(C_mat)) {
    stop(
      paste0(
        "Number of mixtures (columns) must match between Y_mat and C_mat. ",
        "Found ", ncol(Y_mat), " in Y_mat and ", ncol(C_mat), " in C_mat."
      )
    )
  }

  if (!identical(rownames(X_mat), rownames(Y_mat))) {
    stop("Gene labels (index) must be identical between X_mat and Y_mat.")
  }

  if (!identical(colnames(X_mat), rownames(C_mat))) {
    stop("Cell type labels must be identical between X_mat (columns) and C_mat (rows).")
  }

  genes <- rownames(X_mat)
  p <- length(genes)
  celltypes <- colnames(X_mat)


  torch_inputs <- XYCtoTorch(X_mat, Y_mat, C_mat)
  X <- torch_inputs$X
  Y <- torch_inputs$Y
  C <- torch_inputs$C


  torch::torch_manual_seed(42)
  weights <- torch::torch_tensor(
    stats::runif(p, min = 0.001, max = 0.1),
    dtype = torch::torch_float(),
    requires_grad = TRUE
  )

  opt <- torch::optim_adam(params = list(weights), lr = 0.001)

  losses <- numeric(iterations)
  mean_corr <- numeric(iterations)
  all_corr <- matrix(NA_real_, nrow = iterations, ncol = ncol(X_mat))
  colnames(all_corr) <- celltypes

  loss_func <- get_loss_function(func)

  pb <- utils::txtProgressBar(min = 0, max = iterations, style = 3)

  for (i in seq_len(iterations)) {
    preds <- C_est(X, Y, weights^2)

    loss_out <- loss_func(C, preds)
    loss <- loss_out$loss
    all_r <- loss_out$r

    loss$backward()
    opt$step()
    opt$zero_grad()

    losses[i] <- as.numeric(loss$item())
    all_corr[i, ] <- as.numeric(as_array(all_r$detach()))
    mean_corr[i] <- mean(all_corr[i, ], na.rm = TRUE)

    utils::setTxtProgressBar(pb, i)
  }

  close(pb)

  # Format returns
  gamma <- as.numeric(as_array(weights$detach()^2))
  gamma <- length(gamma) * gamma / sum(gamma)

  gamma <- data.frame(
    "gene weights" = gamma,
    row.names = genes,
    check.names = FALSE
  )


  if (plot) {
    if (!is.null(path_plot)) {
      grDevices::pdf(path_plot)
      on.exit(grDevices::dev.off(), add = TRUE)
    }

    matplot(
      x = seq_len(iterations),
      y = all_corr,
      type = "l",
      lty = 1,
      xlab = "Iterations",
      ylab = "corr(C_hat, C)",
      col = seq_len(ncol(all_corr))
    )

    lines(seq_len(iterations), mean_corr, col = "black", lty = 3, lwd = 3)

    legend(
      "bottomright",
      legend = c(celltypes, "AVG"),
      col = c(seq_len(ncol(all_corr)), "black"),
      lty = c(rep(1, ncol(all_corr)), 3),
      lwd = c(rep(1, ncol(all_corr)), 3),
      bty = "n"
    )
  }


  list(
    gamma = gamma,
    mean_corr = mean_corr,
    all_corr = all_corr,
    losses = losses
  )
}