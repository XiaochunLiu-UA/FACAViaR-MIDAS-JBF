###############################################################################
#  Joint VaR-ES Estimation: CAViaR-RMIDAS Models (Rcpp-accelerated)
#  --------------------------------------------------------------------------
#  R wrapper for CAViaR_MIDAS_cpp.cpp fused C++ functions.
#
#  Key argument:  fix_w1 = TRUE  =>  w1 = 1 fixed (restricted, default)
#                 fix_w1 = FALSE =>  w1 free    (unrestricted)
#
#  Convention: VaR and ES are negative for tau < 0.5; mu = 0
###############################################################################

# ===========================================================================
# 0. Dependencies
# ===========================================================================
suppressPackageStartupMessages({
  if (!requireNamespace("numDeriv", quietly = TRUE))
    install.packages("numDeriv", repos = "https://cloud.r-project.org")
  if (!requireNamespace("Rcpp", quietly = TRUE))
    install.packages("Rcpp", repos = "https://cloud.r-project.org")
  library(numDeriv)
  library(Rcpp)
})

# Source the C++ file (compile once per session)
#Rcpp::sourceCpp("/home/xliu121/Research/CAViaR_MIDAS/est_caviar_midas/CAViaR_MIDAS_cpp.cpp")

Rcpp::sourceCpp("CAViaR_MIDAS_cpp.cpp")

# ===========================================================================
# 1. Parameter helpers
# ===========================================================================

n_params <- function(model, K, fix_w1 = TRUE) {
  n_w <- ifelse(fix_w1, 1L, 2L)
  1L + K + n_w + ifelse(model == "SAV", 2L, 3L) + 1L
}

unpack_params <- function(theta, model, K, fix_w1 = TRUE) {
  idx <- 1L
  phi <- theta[idx]; idx <- idx + 1L
  psi <- theta[idx:(idx + K - 1L)]; idx <- idx + K

  if (fix_w1) {
    w1 <- 1.0; w2 <- theta[idx]; idx <- idx + 1L
  } else {
    w1 <- theta[idx]; idx <- idx + 1L
    w2 <- theta[idx]; idx <- idx + 1L
  }

  if (model == "SAV") {
    alpha   <- theta[idx]; idx <- idx + 1L
    beta_ar <- theta[idx]; idx <- idx + 1L
    sr <- list(alpha = alpha, beta_ar = beta_ar)
  } else {
    alpha_pos <- theta[idx]; idx <- idx + 1L
    alpha_neg <- theta[idx]; idx <- idx + 1L
    beta_ar   <- theta[idx]; idx <- idx + 1L
    sr <- list(alpha_pos = alpha_pos, alpha_neg = alpha_neg,
               beta_ar = beta_ar)
  }
  gamma0 <- theta[idx]

  list(phi = phi, psi = psi, w1 = w1, w2 = w2, sr = sr, gamma0 = gamma0)
}

param_names <- function(model, K, fix_w1 = TRUE) {
  nm <- c("phi", paste0("psi_", 1:K))
  if (fix_w1) nm <- c(nm, "w2")
  else        nm <- c(nm, "w1", "w2")
  if (model == "SAV") nm <- c(nm, "alpha", "beta_ar")
  else                nm <- c(nm, "alpha_pos", "alpha_neg", "beta_ar")
  c(nm, "gamma0_ES")
}

# ===========================================================================
# 2. Dispatcher: select the C++ fused function
# ===========================================================================

get_obj_fun_cpp <- function(model, lr_type) {
  key <- paste0(lr_type, "_", tolower(model))
  switch(key,
    linear_sav      = neg_al_loglik_midas_linear_sav_cpp,
    linear_as       = neg_al_loglik_midas_linear_as_cpp,
    exponential_sav = neg_al_loglik_midas_exp_sav_cpp,
    exponential_as  = neg_al_loglik_midas_exp_as_cpp,
    stop("Unknown model/lr_type combination: ", key)
  )
}

       
                                
# ===========================================================================
# 3. Starting vector generation
# ===========================================================================

#' Generate random starting vectors: n_init is how many initials;
## emp_q=q_init is the empirical quantile at tau
generate_starts_midas <- function(n_init, model, lr_type, K, emp_q, fix_w1 = TRUE) {
  np <- n_params(model, K, fix_w1)
  starts <- matrix(NA_real_, nrow = n_init, ncol = np)

  for (i in seq_len(n_init)) {
    idx <- 1L

    # phi
    if (lr_type == "exponential") {
      starts[i, idx] <- runif(1, 2 * emp_q, 0.01 * emp_q)
    } else {
      starts[i, idx] <- runif(1, 2 * emp_q, 0)
    }
    idx <- idx + 1L

    # psi_k
    for (k in 1:K) {
      if (lr_type == "exponential") starts[i, idx] <- runif(1, -1, 1)
      else                          starts[i, idx] <- runif(1, -2, 0)
      idx <- idx + 1L
    }

    # w1 (only if unrestricted)
    if (!fix_w1) {
      starts[i, idx] <- runif(1, 0.5, 5); idx <- idx + 1L
    }

    # w2
    starts[i, idx] <- runif(1, 2,10 ); idx <- idx + 1L

    # Short-run
    if (model == "SAV") {
      starts[i, idx]     <- runif(1, -0.5, 0)
      starts[i, idx + 1] <- runif(1, 0.5, 0.99)
      idx <- idx + 2L
    } else {
      starts[i, idx]     <- runif(1, -0.3, 0.3)
      starts[i, idx + 1] <- runif(1, -0.5, 0)
      starts[i, idx + 2] <- runif(1, 0.5, 0.99)
      idx <- idx + 3L
    }

    # gamma0
    starts[i, idx] <- runif(1, -2, 2)
  }
  starts
}

# ===========================================================================
# 4. Bounds for L-BFGS-B
# ===========================================================================

get_bounds_midas <- function(model, lr_type, K, fix_w1 = TRUE) {
  np <- n_params(model, K, fix_w1)
  lower <- rep(-Inf, np)
  upper <- rep( Inf, np)

  idx <- K + 2L     # skip phi(1) + psi(K); now at first weight param

  if (!fix_w1) {
    lower[idx] <- 0.01; upper[idx] <- 100; idx <- idx + 1L  # w1
  }
  lower[idx] <- 1; upper[idx] <- 100; idx <- idx + 1L       # w2

  if (model == "SAV") {
    idx <- idx + 1L                                          # alpha
    lower[idx] <- 0; upper[idx] <- 0.9999                   # beta_ar
  } else {
    idx <- idx + 2L                                          # alpha_pos, alpha_neg
    lower[idx] <- 0; upper[idx] <- 0.9999                   # beta_ar
  }

  list(lower = lower, upper = upper)
}

# ===========================================================================
# 5. Main estimation function
# ===========================================================================
 
      
estimate_caviar_midas_es_cpp <- function(y, factors, month_id,
                                     tau = 0.05,
                                     model = c("SAV", "AS"),
                       lr_type = c("linear", "exponential"),
                                     trim,
                                     fix_w1 = TRUE,
                                     S = 12L, n0 = 300L,
                                     n_init = 5000L, m_best = 10L,
                                     n_burn = 0L, maxit = 5000L,
                                     seed = 42L,compute_se = TRUE) {

  model   <- match.arg(model)
  lr_type <- match.arg(lr_type)

  if (is.null(dim(factors))) factors <- matrix(factors, ncol = 1)
  K  <- ncol(factors)
  TT <- length(y)
  stopifnot(TT > n0 + 10)
  stopifnot(length(month_id) == TT)

  factors  <- as.matrix(factors)
  month_id <- as.integer(month_id)

  w_label <- ifelse(fix_w1, "w1=1 fixed", "w1 free")
#  cat(sprintf("Estimating CAViaR-MIDAS-%s (%s LR, %s) | K=%d | tau=%.3f\n",    model, lr_type, w_label, K, tau))

  # ---- initial quantile ---------------------------------------------------
  q_init <- as.numeric(quantile(y[1:n0], tau))
  emp_q  <- q_init
  np     <- n_params(model, K, fix_w1)

 # cat(sprintf("  Parameters: %d | q_init: %.6f\n", np, q_init))

  # ---- C++ fused function --------------------------------------------------
  cpp_fun <- get_obj_fun_cpp(model, lr_type)

  obj_fun <- function(theta) {
    tryCatch(
      cpp_fun(theta, y, factors, month_id, tau, q_init, S, n_burn, fix_w1, trim),
      error = function(e) 1e12
    )
  }

  # ---- bounds --------------------------------------------------------------
  bnds <- get_bounds_midas(model, lr_type, K, fix_w1)

  # ---- random starts -------------------------------------------------------
  set.seed(seed)
  starts <- generate_starts_midas(n_init, model, lr_type, K, emp_q, fix_w1)

  # ---- evaluate objective at all starts ------------------------------------
#  cat("  Evaluating starts...")
  t0 <- proc.time()
  obj_vals <- apply(starts, 1, obj_fun)
  n_finite <- sum(is.finite(obj_vals) & obj_vals < 1e11)
 # cat(sprintf(" %d/%d finite (%.1fs)\n", n_finite, n_init,
 #             (proc.time() - t0)[3]))

  if (n_finite == 0)
    stop("No feasible starting point found. Check data and parameter ranges.")

  # ---- keep m_best ---------------------------------------------------------
  best_idx <- order(obj_vals)[1:min(m_best, n_finite)]

  # ---- refine: Nelder-Mead -> L-BFGS-B -> Nelder-Mead ---------------------
 # cat(sprintf("  Refining top %d...\n", length(best_idx)))
  t0 <- proc.time()

  results <- lapply(seq_along(best_idx), function(j) {
    th <- starts[best_idx[j], ]
  opt <- tryCatch(
      optim(th, obj_fun, method = "L-BFGS-B",
            lower = bnds$lower, upper = bnds$upper,
            control = list(maxit = maxit, factr = 1e-6)),
      error = function(e) list(par = th, value = 1e12, convergence = 1))
      opt
      
 #   opt1 <- tryCatch(
#      optim(th, obj_fun, method = "Nelder-Mead",
#            control = list(maxit = maxit, reltol = 1e-10)),
#      error = function(e) list(par = th, value = 1e12, convergence = 1))
#    opt2 <- tryCatch(
#      optim(opt1$par, obj_fun, method = "L-BFGS-B",
#            lower = bnds$lower, upper = bnds$upper,
#            control = list(maxit = maxit, factr = 1e4)),
#      error = function(e) opt1)
#    opt3 <- tryCatch(
#      optim(opt2$par, obj_fun, method = "Nelder-Mead",
#            control = list(maxit = maxit, reltol = 1e-10)),
#      error = function(e) opt2)
  #  cat(sprintf("    Start %d: nll = %.4f\n", j, opt3$value))
#    opt3
  })

 # cat(sprintf("  Optimisation: %.1fs\n", (proc.time() - t0)[3]))

  # ---- global optimum ------------------------------------------------------
  obj_final <- sapply(results, function(r) r$value)
  best      <- results[[which.min(obj_final)]]
  theta_hat <- best$par

  #cat(sprintf("  Best neg-loglik: %.4f\n", best$value))

  # ---- extract fitted values -----------------------------------------------
  p <- unpack_params(theta_hat, model, K, fix_w1)

  LR_monthly <- as.numeric(
    compute_LR_monthly_cpp(factors, p$psi, p$w1, p$w2, S, p$phi,
                           lr_type == "exponential", trim)
  )
  LR_daily <- LR_monthly[month_id]

  if (model == "SAV") {
    rec <- caviar_midas_sav_cpp(p$sr$alpha, p$sr$beta_ar,
                                y, LR_daily, q_init)
  } else {
    rec <- caviar_midas_as_cpp(p$sr$alpha_pos, p$sr$alpha_neg,
                               p$sr$beta_ar, y, LR_daily, q_init)
  }

  Q       <- as.numeric(rec$Q)
  SR      <- as.numeric(rec$SR)
  es_mult <- 1 + exp(p$gamma0)
  ES      <- es_mult * Q

  # ---- standard errors -------
  
  if (compute_se) {
  se <- tryCatch({
    H <- numDeriv::hessian(obj_fun, theta_hat)
    Hinv <- solve(H)
    sqrt(abs(diag(Hinv)))
  }, error = function(e) rep(NA_real_, np))
  } else {
    se <- rep(NA_real_, np)
  }

  # ---- parameter names -----------------------------------------------------
  pnames <- param_names(model, K, fix_w1)
  names(theta_hat) <- pnames
  names(se)        <- pnames

  # ---- weights at optimum --------------------------------------------------
  midas_weights <- as.numeric(beta_weights_cpp(S, p$w1, p$w2))

  list(
    model         = model,
    lr_type       = lr_type,
    tau           = tau,
    K             = K,
    S             = S,
    fix_w1        = fix_w1,
    theta         = theta_hat,
    se            = se,
    phi           = p$phi,
    psi           = p$psi,
    w1            = p$w1,
    w2            = p$w2,
    sr_params     = p$sr,
    gamma0        = p$gamma0,
    es_mult       = es_mult,
    Q             = Q,
    ES            = ES,
    SR            = SR,
    LR_daily      = LR_daily,
    LR_monthly    = LR_monthly,
    midas_weights = midas_weights,
    neg_loglik    = best$value,
    q_init        = q_init,
    converged     = (best$convergence == 0)
  )
}

# ===========================================================================
# 6. One-step-ahead forecast
# ===========================================================================

forecast_caviar_midas_es <- function(fit, y_last, factors_new, trim) {
  if (is.null(dim(factors_new))) factors_new <- matrix(factors_new, ncol = 1)

  LR_monthly_new <- as.numeric(
    compute_LR_monthly_cpp(as.matrix(factors_new), fit$psi,
                           fit$w1, fit$w2, fit$S, fit$phi,
                           fit$lr_type == "exponential", trim)
  )
  LR_new  <- tail(LR_monthly_new, 1)
  SR_last <- tail(fit$SR, 1)

  if (fit$model == "SAV") {
    SR_new <- fit$sr_params$alpha * abs(y_last) +
              fit$sr_params$beta_ar * SR_last
  } else {
    SR_new <- fit$sr_params$alpha_pos * max(y_last, 0) +
              fit$sr_params$alpha_neg * (-min(y_last, 0)) +
              fit$sr_params$beta_ar * SR_last
  }

  Q_new  <- LR_new + SR_new
  ES_new <- fit$es_mult * Q_new

  list(VaR_forecast = Q_new, ES_forecast = ES_new, LR_forecast = LR_new)
}

# ===========================================================================
# 7. Quantile loss (check function)
# ===========================================================================

#' Quantile check-loss (tick loss)
#'
#' Computes the quantile loss  rho_tau(u) = u * (tau - I{u < 0})
#' where u = y - Q.
#'
#' @param y     numeric vector of realised returns
#' @param Q     numeric vector of VaR forecasts (same length as y)
#' @param tau   quantile level (e.g. 0.05)
#' @return list with:
#'   \item{total}{sum of tick losses}
#'   \item{mean}{average tick loss}
#'   \item{n}{number of observations}
quantile_loss <- function(y, Q, tau) {
  u    <- y - Q
  loss <- u * (tau - (u < 0))
  list(total = sum(loss), mean = mean(loss), n = length(y))
}

# ===========================================================================
# 8. Summary display
# ===========================================================================

print_summary <- function(fit) {
  w_label <- ifelse(fit$fix_w1, "w1=1 fixed", "w1 free")
  cat(sprintf("\n=== CAViaR-MIDAS-%s (%s, %s) | tau=%.3f | K=%d ===\n",
              fit$model, fit$lr_type, w_label, fit$tau, fit$K))

  est <- fit$theta; se <- fit$se
  tab <- data.frame(Estimate = round(est, 6), Std.Err = round(se, 6),
                    t.stat = round(est / se, 3), row.names = names(est))
  print(tab)

  cat(sprintf("\nNeg-loglik: %.4f | ES/VaR: %.4f | Converged: %s\n",
              fit$neg_loglik, fit$es_mult, fit$converged))
  cat(sprintf("Beta weights (S=%d, w1=%.3f, w2=%.3f):\n",
              fit$S, fit$w1, fit$w2))
  cat(sprintf("  First 6: %s\n",
              paste(round(head(fit$midas_weights, 6), 4), collapse = ", ")))
  cat("\n")
}
