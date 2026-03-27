##############################################################

hps_setup <- function(
    X_ref,
    Y_test,
    k_folds = 5,
    lambdas = 10^seq(-20, 0, length.out = 21)
) {

  X_df <- as.data.frame(X_ref, check.names = FALSE)
  Y_df <- as.data.frame(Y_test, check.names = FALSE)

  X_df[] <- lapply(X_df, as.numeric)
  Y_df[] <- lapply(Y_df, as.numeric)

  X_mat <- as.matrix(X_df)
  Y_mat <- as.matrix(Y_df)

  gamma_vec <- (1 / ncol(Y_df)) * rep(1, nrow(Y_df)) / (rowMeans(Y_df)^2)
  gamma <- data.frame(gamma_vec, row.names = rownames(X_df), check.names = FALSE)
  G_mat <- diag(sqrt(as.numeric(gamma[, 1])))


  colnames(Y_df) <- as.character(seq_len(ncol(Y_df)) - 1)
  colnames(Y_mat) <- colnames(Y_df)


  message("Preparing Job List")
  jobs <- list()

  set.seed(42)
  fold_ids <- sample(rep(seq_len(k_folds), length.out = ncol(Y_df)))

  fold_num <- 0
  for (fold in seq_len(k_folds)) {
    test_index <- which(fold_ids == fold)
    train_index <- setdiff(seq_len(ncol(Y_df)), test_index)

    for (lmbda in lambdas) {
      jobs[[length(jobs) + 1]] <- list(
        fold = fold_num,
        train_ids = train_index,
        test_ids = test_index,
        lambda = lmbda
      )
    }

    fold_num <- fold_num + 1
  }


  message("Preparing Baseline Model")
  model_baseline <- adtd_run(
    X_mat = X_df,
    Y_mat = Y_df,
    gamma = gamma,
    C_static = TRUE,
    Delta_static = TRUE,
    max_iterations = 1000,
    verbose = FALSE
  )

  Cc_mat <- rbind(
    as.matrix(model_baseline$C_est),
    as.matrix(model_baseline$c_est)
  )

  x_df <- model_baseline$x_est

  list(
    X_df = X_df,
    X_mat = X_mat,
    Y_df = Y_df,
    Y_mat = Y_mat,
    gamma = gamma,
    G_mat = G_mat,
    k_folds = k_folds,
    lambdas = lambdas,
    jobs = jobs,
    Cc_mat = Cc_mat,
    x_df = x_df,
    validation_losses = data.frame(
      fold = numeric(),
      lambda = numeric(),
      loss = numeric()
    ),
    results_raw = NULL,
    results = NULL
  )
}



hps_generalization_loss <- function(state, Delta, test_ids) {
  Y_test_mat <- as.matrix(state$Y_df[, test_ids, drop = FALSE])
  Cc_test_mat <- state$Cc_mat[, test_ids, drop = FALSE]
  Xx_mat <- cbind(as.matrix(Delta) * state$X_mat, as.matrix(state$x_df))
  recon_error <- Y_test_mat - (Xx_mat %*% Cc_test_mat)

  loss_raw <- norm(recon_error, type = "F")^2
  loss_weighted <- norm(state$G_mat %*% recon_error, type = "F")^2

  c(loss_raw = loss_raw, loss_weighted = loss_weighted)
}

hps_run_job <- function(job, state) {
  fold <- job$fold
  lmbda <- job$lambda
  train_ids <- job$train_ids
  test_ids <- job$test_ids

  model <- adtd_run(
    X_mat = state$X_df,
    Y_mat = state$Y_df[, train_ids, drop = FALSE],
    gamma = state$gamma,
    max_iterations = 1000,
    C_static = TRUE,
    Delta_static = FALSE,
    lambda2 = lmbda,
    verbose = FALSE
  )

  losses <- hps_generalization_loss(state, model$Delta_est, test_ids)

  list(
    fold = fold,
    lambda = lmbda,
    loss = unname(losses["loss_raw"])
    # loss_weighted = unname(losses["loss_weighted"])
  )
}



hps_run <- function(state, n_workers = 10) {
  jobs <- state$jobs

  if (is.null(jobs)) {
    stop("No jobs found. Please assign a list of jobs before running.")
  }

  results <- vector("list", length(jobs))

  if (n_workers <= 1) {
    pb <- txtProgressBar(min = 0, max = length(jobs), style = 3)
    for (i in seq_along(jobs)) {
      results[[i]] <- hps_run_job(jobs[[i]], state)
      setTxtProgressBar(pb, i)
    }
    close(pb)
  } else {
    cl <- parallel::makeCluster(n_workers)
    on.exit(parallel::stopCluster(cl), add = TRUE)

    parallel::clusterExport(
      cl,
      varlist = c(
        "state",
        "jobs",
        "hps_run_job",
        "hps_generalization_loss",
        "adtd_run",
        "adtd_setup",
        "adtd_update_C0",
        "adtd_update_C",
        "adtd_update_Delta",
        "adtd_update_x",
        "adtd_loss",
        "adtd_loss_static",
        "solve_qp_nonneg_eq",
        "solve_qp_nonneg_eq_clarabel",
        "solve_qp_nonneg_eq_quadprog"
      ),
      envir = environment()
    )

    parallel::clusterEvalQ(cl, {
      library(Matrix)
      library(nnls)
      library(quadprog)
      NULL
    })

    pb <- txtProgressBar(min = 0, max = length(jobs), style = 3)

    for (i in seq_along(jobs)) {
      # parLapplyLB one job at a time to preserve progress reporting
      results[[i]] <- parallel::parLapplyLB(cl, list(jobs[[i]]), function(job) {
        hps_run_job(job, state)
      })[[1]]
      setTxtProgressBar(pb, i)
    }

    close(pb)
  }

  # Write results to log file with timestamp
  log_filename <- paste0(
    "deconomix_run_results_",
    format(Sys.time(), "%Y%m%d_%H%M%S"),
    ".log"
  )

  tryCatch({
    con <- file(log_filename, open = "wt")
    writeLines(
      c(
        paste("Run timestamp:", format(Sys.time(), "%Y-%m-%dT%H:%M:%S")),
        "Results:",
        vapply(results, function(x) paste(capture.output(str(x)), collapse = " "), character(1))
      ),
      con = con
    )
    close(con)
  }, error = function(e) {
    message("Could not write log file: ", e$message)
  })

  results_raw <- do.call(
    rbind,
    lapply(results, function(x) {
      data.frame(
        fold = x$fold,
        lambda = x$lambda,
        loss = x$loss,
        stringsAsFactors = FALSE
      )
    })
  )

  results_summary <- aggregate(
    loss ~ lambda,
    data = results_raw,
    FUN = function(z) c(mean = mean(z), sd = stats::sd(z))
  )

  results <- data.frame(
    lambda = results_summary$lambda,
    loss.mean = results_summary$loss[, "mean"],
    loss.sd = results_summary$loss[, "sd"],
    check.names = FALSE
  )
  rownames(results) <- results$lambda

  state$results_raw <- results_raw
  state$results <- results

  state
}


get_lambda_min <- function(results) {
  lambdas <- as.numeric(results$lambda)
  lambdas[which.min(results$loss.mean)]
}

get_lambda_1se <- function(results) {
  avgLoss <- results$loss.mean
  stdLoss <- results$loss.sd
  lambdas <- as.numeric(results$lambda)

  idx_min <- which.min(avgLoss)
  minMean <- lambdas[idx_min]
  minMeanValue <- avgLoss[idx_min]
  std_at_min <- stdLoss[idx_min]

  threshold <- minMeanValue + std_at_min

  eligible <- lambdas[avgLoss <= threshold]
  lambda_1se <- max(eligible) * 10

  if (identical(lambda_1se, tail(lambdas, 1))) {
    message("Warning: No index within 1se. Returning minimum.")
    return(minMean)
  }

  lambda_1se
}



plot_results <- function(results, title = NULL, path = NULL) {
  lambdas <- as.numeric(results$lambda)
  loss_raw_mean <- results$loss.mean
  loss_raw_sd <- results$loss.sd

  if (!is.null(path)) {
    grDevices::pdf(path)
    on.exit(grDevices::dev.off(), add = TRUE)
  }

  plot(
    lambdas, loss_raw_mean,
    log = "x",
    type = "b",
    pch = 16,
    lty = 1,
    col = "black",
    xlab = expression(lambda[2]),
    ylab = "Average Out-of-Distribution Error",
    main = title
  )

  polygon(
    c(lambdas, rev(lambdas)),
    c(loss_raw_mean - loss_raw_sd, rev(loss_raw_mean + loss_raw_sd)),
    col = adjustcolor("steelblue", alpha.f = 0.2),
    border = NA
  )

  lines(lambdas, loss_raw_mean, type = "b", pch = 16)

  lambda_min <- get_lambda_min(results)
  lambda_1se <- get_lambda_1se(results)

  idx_min <- which.min(abs(lambdas - lambda_min))
  idx_1se <- which.min(abs(lambdas - lambda_1se))

  min_mean <- loss_raw_mean[idx_min]
  one_se_mean <- loss_raw_mean[idx_1se]

  points(lambda_min, min_mean, col = "black", pch = 16, cex = 1.2)
  points(lambda_1se, one_se_mean, col = "green3", pch = 16, cex = 1.2)

  legend(
    "topright",
    legend = c("Test Loss", "Std. deviation", expression(lambda[2]~min), expression(lambda[2]~"1se")),
    col = c("black", "steelblue", "black", "green3"),
    pch = c(16, 15, 16, 16),
    lty = c(1, NA, NA, NA),
    pt.cex = c(1, 2, 1, 1),
    bty = "n"
  )
}