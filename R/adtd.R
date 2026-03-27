
solve_qp_nonneg_eq_quadprog <- function(P, q, Aeq, b, ridge = 1e-10) {
  P <- as.matrix(P)
  P <- (P + t(P)) / 2
  P <- P + diag(ridge, nrow(P))

  q <- as.numeric(q)

  if (is.null(dim(Aeq))) {
    Aeq <- matrix(Aeq, nrow = 1)
  } else {
    Aeq <- as.matrix(Aeq)
  }

  nvars <- ncol(P)

  if (ncol(Aeq) != nvars) {
    stop("Aeq must have ncol equal to length of optimization variable.")
  }


  Amat <- cbind(t(Aeq), diag(nvars))
  bvec <- c(as.numeric(b), rep(0, nvars))

  sol <- quadprog::solve.QP(
    Dmat = P,
    dvec = -q,
    Amat = Amat,
    bvec = bvec,
    meq = nrow(Aeq)
  )

  sol$solution
}


solve_qp_nonneg_eq_clarabel <- function(P, q, Aeq, b) {
  if (!requireNamespace("clarabel", quietly = TRUE)) {
    return(NULL)
  }

  q <- as.numeric(q)
  nvars <- length(q)

  if (is.null(dim(Aeq))) {
    Aeq <- matrix(Aeq, nrow = 1)
  }

  Aeq <- Matrix::Matrix(Aeq, sparse = TRUE)
  Gmat <- -Matrix::Diagonal(nvars)  # -x <= 0
  h <- rep(0, nvars)

  A <- rbind(Aeq, Gmat)
  rhs <- c(as.numeric(b), h)

  cones <- list(
    z = nrow(Aeq),  
    l = nvars       
  )

  P <- Matrix::forceSymmetric(Matrix::Matrix(P, sparse = TRUE), uplo = "U")

  sol <- tryCatch(
    clarabel::clarabel(
      A = A,
      b = rhs,
      q = q,
      P = P,
      cones = cones
    ),
    error = function(e) NULL
  )

  if (is.null(sol) || is.null(sol$x)) {
    return(NULL)
  }

  as.numeric(sol$x)
}

solve_qp_nonneg_eq <- function(P, q, Aeq, b) {
  sol <- solve_qp_nonneg_eq_clarabel(P, q, Aeq, b)
  if (!is.null(sol)) {
    return(sol)
  }
  solve_qp_nonneg_eq_quadprog(P, q, Aeq, b)
}



adtd_setup <- function(state) {
  A_nnls <- as.matrix(state$G %*% state$X)
  Y_nnls <- as.matrix(state$G %*% state$Y)

  C0_init <- matrix(0, nrow = state$q, ncol = state$n)

  for (i in seq_len(state$n)) {
    fit <- nnls::nnls(A_nnls, Y_nnls[, i])
    C0_init[, i] <- coef(fit)
  }

  x_base <- rowMeans(state$Y - state$X %*% C0_init)
  x_base <- pmax(x_base, 0)
  x <- matrix(x_base / sum(x_base), ncol = 1)

  state$C0_init <- C0_init
  state$x <- x
  state
}

adtd_update_C0 <- function(state) {
  A <- as.matrix(state$G %*% state$Y)                    
  B <- cbind(as.matrix(state$G %*% state$X), as.matrix(state$G %*% state$x))

  P <- 2.0 * crossprod(B)                                 
  Q <- -2.0 * crossprod(B, A)                              

  A_eq <- matrix(c(rep(1, state$q), 1), nrow = 1)
  b <- 1

  C_tilde <- rbind(state$C, state$c)

  for (i in seq_len(state$n)) {
    sol <- solve_qp_nonneg_eq(P, Q[, i], A_eq, b)
    C_tilde[, i] <- sol
  }

  state$C0_init <- C_tilde[seq_len(state$q), , drop = FALSE]
  state$c <- C_tilde[state$q + 1, , drop = FALSE]

  state$C0_init[state$C0_init < 0] <- 0
  state$c[state$c < 0] <- 0

  state
}

adtd_update_C <- function(state) {
  if (state$C_static) {
    state$C <- state$C0_init
    state$c <- matrix(1 - colSums(state$C), nrow = 1)
    state$c[state$c < 0] <- 0
    return(state)
  }

  a1 <- as.matrix(state$G %*% state$Y)
  a2 <- sqrt(state$lambda1) * state$C0_init
  A <- rbind(a1, a2)

  b1 <- cbind(
    as.matrix(state$G %*% (state$Delta * state$X)),
    as.matrix(state$G %*% state$x)
  )
  b2 <- cbind(
    sqrt(state$lambda1) * diag(state$q),
    matrix(0, nrow = state$q, ncol = 1)
  )
  B <- rbind(b1, b2)

  P <- 2.0 * crossprod(B)
  Q <- -2.0 * crossprod(B, A)

  A_eq <- matrix(c(rep(1, state$q), 1), nrow = 1)
  b <- 1

  C_tilde <- rbind(state$C, state$c)

  for (i in seq_len(state$n)) {
    sol <- solve_qp_nonneg_eq(P, Q[, i], A_eq, b)
    C_tilde[, i] <- sol
  }

  state$C <- C_tilde[seq_len(state$q), , drop = FALSE]
  state$c <- C_tilde[state$q + 1, , drop = FALSE]

  state$C[state$C < 0] <- 0
  state$c[state$c < 0] <- 0

  state
}

adtd_update_Delta <- function(state) {
  if (state$Delta_static) {
    return(state)
  }


  Y_tilde <- state$G %*% (state$Y - state$x %*% state$c)  
  CCT <- state$C %*% t(state$C)                            


  a1 <- kronecker(Matrix::Matrix(CCT, sparse = TRUE), state$Gamma)

  
  pq <- state$p * state$q
  a2 <- Matrix::Matrix(0, nrow = pq, ncol = pq, sparse = TRUE)

  for (i in seq_len(state$q)) {
    idx_i <- ((i - 1) * state$p + 1):(i * state$p)

    for (j in i:state$q) {
      idx_j <- ((j - 1) * state$p + 1):(j * state$p)
      block <- Matrix::Diagonal(x = state$X[, i] * state$X[, j])

      a2[idx_i, idx_j] <- block
      if (i != j) {
        a2[idx_j, idx_i] <- block
      }
    }
  }


  a12 <- a1 * a2


  a3 <- state$lambda2 * Matrix::Diagonal(pq)

  P <- 2.0 * (a12 + a3)
  P <- Matrix::forceSymmetric(P, uplo = "U")


  kron_eye <- kronecker(Matrix::Diagonal(state$q), state$G)
  X_flat <- as.vector(state$X)                          
  Y_tilde_C <- state$C %*% t(Y_tilde)                      


  ytc_flat_row_major <- as.vector(t(Y_tilde_C))

  b1 <- as.numeric(
    matrix(ytc_flat_row_major, nrow = 1) %*%
      kron_eye %*%
      Matrix::Diagonal(x = X_flat)
  )

  Q <- -2.0 * (b1 + state$lambda2 * rep(1, pq))


  Aeq <- Matrix::Matrix(0, nrow = state$q, ncol = pq, sparse = TRUE)
  for (i in seq_len(state$q)) {
    idx <- ((i - 1) * state$p + 1):(i * state$p)
    Aeq[i, idx] <- state$X[, i]
  }
  b <- rep(1, state$q)

  sol <- solve_qp_nonneg_eq(P, Q, Aeq, b)


  if (is.null(sol)) {
    stop("No solution returned in adtd_update_Delta().")
  }

  state$Delta <- matrix(sol, nrow = state$p, ncol = state$q, byrow = FALSE)
  state
}

adtd_update_x <- function(state) {
  if (all(state$c == 0)) {
    state$x <- matrix(rep(1 / state$p, state$p), ncol = 1)
    return(state)
  }

  Z <- state$Y - (state$Delta * state$X) %*% state$C

  P <- 2.0 * state$Gamma * sum(state$c^2)
  q <- -2.0 * as.numeric(state$c %*% t(Z) %*% state$Gamma)

  Aeq <- matrix(rep(1, state$p), nrow = 1)
  b <- 1

  sol <- solve_qp_nonneg_eq(P, q, Aeq, b)
  if (is.null(sol)) {
    stop("No solution returned in adtd_update_x().")
  }

  state$x <- matrix(sol, ncol = 1)
  state$x[state$x < 0] <- 0

  state
}



adtd_loss <- function(state) {
  ltemp <- state$G %*% (state$Y - (state$Delta * state$X) %*% state$C - (state$x %*% state$c))
  term1 <- sum(ltemp^2)
  term2 <- state$lambda1 * sum((state$C - state$C0_init)^2)
  term3 <- state$lambda2 * sum((matrix(1, nrow = state$p, ncol = state$q) - state$Delta)^2)

  c(total = term1 + term2 + term3, rss = term1, bias = term2 + term3)
}

adtd_loss_static <- function(state) {
  sum((state$G %*% (state$Y - state$X %*% state$C - state$x %*% state$c))^2)
}



adtd_run <- function(
    X_mat,
    Y_mat,
    gamma,
    lambda1 = 1.0,
    lambda2 = 1.0,
    max_iterations = 200,
    eps = 1e-8,
    C_static = FALSE,
    Delta_static = FALSE,
    gamma_offset = TRUE,
    delta_stepsize = 1,
    verbose = TRUE
) {
  # Input handling
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

  genes <- rownames(X_mat)
  celltypes <- colnames(X_mat)
  mixture_names <- colnames(Y_mat)

  X <- as.matrix(X_mat)
  Y <- as.matrix(Y_mat)
  gamma_vec <- as.numeric(as.matrix(gamma)[, 1])

  if (gamma_offset) {
    gamma_vec <- gamma_vec + 1 / nrow(X)
  }

  # Normalize gamma
  gamma_vec <- gamma_vec / sum(gamma_vec)

  # Normalize expression profiles
  X <- sweep(X, 2, colSums(X), "/")
  Y <- sweep(Y, 2, colSums(Y), "/")

  p <- nrow(X)
  q <- ncol(X)
  n <- ncol(Y)

  state <- list(
    X = X,
    Y = Y,
    gamma = gamma_vec,
    Gamma = Matrix::Diagonal(x = gamma_vec),
    G = Matrix::Diagonal(x = sqrt(gamma_vec)),
    p = p,
    q = q,
    n = n,
    lambda1 = lambda1,
    lambda2 = lambda2,
    max_iterations = max_iterations,
    eps = eps,
    C_static = C_static,
    Delta_static = Delta_static,
    delta_stepsize = delta_stepsize,  
    genes = genes,
    celltypes = celltypes,
    Delta = matrix(1, nrow = p, ncol = q),
    C0_init = NULL,
    x = NULL,
    C = matrix(0, nrow = q, ncol = n),
    c = matrix(stats::runif(n), nrow = 1)
  )

  state <- adtd_setup(state)
  state <- adtd_update_C0(state)

  loss <- numeric(0)

  if (verbose) {
    pb <- utils::txtProgressBar(min = 0, max = max_iterations, style = 3)
    on.exit(close(pb), add = TRUE)
  }

  # Case 1: C_static = TRUE
  if (state$C_static) {
    loss_old <- Inf

    for (i in seq_len(state$max_iterations)) {
      ltemp <- adtd_loss(state)
      loss_new <- ltemp[1]
      loss <- c(loss, loss_new)

      err <- loss_old - loss_new

      if (verbose) {
        utils::setTxtProgressBar(pb, i)
      }

      if (err <= state$eps) {
        if (verbose) {
          message(sprintf("i = %i, err = %1.2e - Convergence reached!", i - 1, err))
        }
        break
      }

      if (i == 1) {
        old_static_loss <- Inf
        new_static_loss <- 0

        for (j in seq_len(1000)) {
          static_err <- old_static_loss - new_static_loss

          if (static_err <= state$eps) {
            state <- adtd_update_C(state)
            state <- adtd_update_x(state)
            state <- adtd_update_Delta(state)
            loss_old <- loss_new
            break
          } else {
            old_static_loss <- adtd_loss_static(state)
            state <- adtd_update_C0(state)
            state$C <- state$C0_init
            state <- adtd_update_x(state)
            new_static_loss <- adtd_loss_static(state)
          }
        }
      } else {
        loss_old <- loss_new
        state <- adtd_update_C(state)
        state <- adtd_update_x(state)
        state <- adtd_update_Delta(state)
      }
    }

  } else {
    # Case 2: C_static = FALSE
    for (i in seq_len(state$max_iterations)) {
      # Save old
      C_copy <- state$C
      c_copy <- state$c
      x_copy <- state$x
      Delta_copy <- state$Delta

      loss_old <- adtd_loss(state)[1]
      loss <- c(loss, loss_old)

      # Update all optimization variables
      state <- adtd_update_C(state)
      state <- adtd_update_x(state)
      state <- adtd_update_Delta(state)

      # Calculate loss and error
      loss_new <- adtd_loss(state)[1]
      err <- loss_old - loss_new

      if (verbose) {
        utils::setTxtProgressBar(pb, i)
      }


      if (i > 2 && err < state$eps) {
        if (err < 0) {
          state$C <- C_copy
          state$c <- c_copy
          state$x <- x_copy
          state$Delta <- Delta_copy
        }

        if (verbose) {
          message(sprintf("i = %i, err = %1.2e, Convergence Reached!", i, err))
        }
        break
      }
    }
  }


  C_est <- as.data.frame(state$C, check.names = FALSE)
  rownames(C_est) <- state$celltypes
  colnames(C_est) <- mixture_names

  c_est <- as.data.frame(state$c, check.names = FALSE)
  rownames(c_est) <- "hidden"
  colnames(c_est) <- mixture_names

  x_est <- data.frame(hidden = as.numeric(state$x), row.names = state$genes, check.names = FALSE)

  Delta_est <- as.data.frame(state$Delta, check.names = FALSE)
  rownames(Delta_est) <- state$genes
  colnames(Delta_est) <- state$celltypes

  list(
    C_est = C_est,
    c_est = c_est,
    x_est = x_est,
    Delta_est = Delta_est,
    loss = loss
  )
}