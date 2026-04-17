library(FKF)
library(foreign)
library(rjags)
library(coda)

# ============================================================
# Helper Function: Fit AR(1)+WN and Extract Residuals
# With multiple starting values to handle convergence issues
# ============================================================
fit_ar1wn_residuals <- function(ts, verbose = FALSE) {
  
  # ========== Log-Likelihood Function ==========
  loglik_ar1wn <- function(parm, y) {
    phi <- parm[1]
    mu <- parm[2]
    log_sigma2_eps <- parm[3]
    log_sigma2_omega <- parm[4]
    
    sigma2_eps <- exp(log_sigma2_eps)
    sigma2_omega <- exp(log_sigma2_omega)
    
    n <- length(y)
    
    # Bounds check
    if (abs(phi) >= 0.99) {
      return(1e10)
    }
    
    # State-space matrices
    m0 <- 0
    C0 <- matrix(sigma2_eps / (1 - phi^2), nrow = 1, ncol = 1)
    d <- mu
    F <- matrix(1, nrow = 1, ncol = 1)
    c <- 0
    A <- matrix(phi, nrow = 1, ncol = 1)
    Sigma_omega <- matrix(sigma2_omega, nrow = 1, ncol = 1)
    Sigma_eps <- matrix(sigma2_eps, nrow = 1, ncol = 1)
    
    # Run Kalman filter
    tryCatch({
      kf <- fkf(
        a0 = m0,
        P0 = C0,
        dt = d,
        ct = c,
        Tt = A,
        Zt = F,
        HHt = Sigma_omega,
        GGt = Sigma_eps,
        yt = matrix(y, nrow = 1, ncol = n)
      )
      
      if (is.na(kf$logLik) || is.infinite(kf$logLik)) {
        return(1e10)
      }
      
      return(-kf$logLik)
    }, error = function(e) {
      return(1e10)
    })
  }
  
  # ========== Generate Multiple Starting Values ==========
  ts_mean <- mean(ts, na.rm = TRUE)
  ts_var <- var(ts, na.rm = TRUE)
  
  starting_values <- list(
    c(phi = 0.1, mu = ts_mean, log_sigma2_eps = log(ts_var/4), log_sigma2_omega = log(ts_var/4)),
    c(phi = 0.3, mu = ts_mean, log_sigma2_eps = log(ts_var/3), log_sigma2_omega = log(ts_var/3)),
    c(phi = 0.5, mu = ts_mean, log_sigma2_eps = log(ts_var/2), log_sigma2_omega = log(ts_var/2)),
    c(phi = 0.7, mu = ts_mean, log_sigma2_eps = log(ts_var/1.5), log_sigma2_omega = log(ts_var/1.5)),
    c(phi = 0.2, mu = ts_mean, log_sigma2_eps = log(ts_var/2.5), log_sigma2_omega = log(ts_var/2.5))
  )
  
  # ========== Try Each Starting Value ==========
  best_fit <- NULL
  best_loglik <- Inf
  
  for (i in seq_along(starting_values)) {
    parm0 <- starting_values[[i]]
    
    if (verbose) {
      cat("Trying starting value", i, "of", length(starting_values), "...\n")
    }
    
    fit <- optim(
      parm0,
      loglik_ar1wn,
      y = ts,
      method = "L-BFGS-B",
      lower = c(-1, -Inf, -Inf, -Inf),
      upper = c(1, Inf, Inf, Inf),
      hessian = FALSE,
      control = list(maxit = 1000)
    )
    
    # Keep best fit
    if (fit$value < best_loglik) {
      best_loglik <- fit$value
      best_fit <- fit
    }
  }
  
  if (verbose) {
    cat("Convergence code:", best_fit$convergence, "\n")
    cat("Log-likelihood:", -best_fit$value, "\n\n")
  }
  
  # ========== Extract Parameter Estimates ==========
  phi_est <- best_fit$par[1]
  mu_est <- best_fit$par[2]
  ivar_est <- as.numeric(exp(best_fit$par[3])) 
  evar_est <- as.numeric(exp(best_fit$par[4]))  
  
  # ========== Check for Heywood Cases ==========
  heywood_ivar <- ivar_est < 1e-6
  heywood_evar <- evar_est < 1e-6
  
  if (verbose) {
    cat("Parameter Estimates:\n")
    cat("====================\n")
    cat("phi (AR coefficient):", phi_est, "\n")
    cat("mu (Mean):", mu_est, "\n")
    cat("Ivar (Innovation variance):", ivar_est, 
        if (heywood_ivar) " [HEYWOOD CASE]" else "", "\n")
    cat("Evar (Measurement error variance):", evar_est, 
        if (heywood_evar) " [HEYWOOD CASE]" else "", "\n\n")
  }
  
  # ========== Run Kalman Filter and Smoother ==========
  n <- length(ts)
  
  kf_final <- fkf(
    a0 = 0,
    P0 = matrix(ivar_est / (1 - phi_est^2), nrow = 1, ncol = 1),
    dt = mu_est,
    ct = 0,
    Tt = matrix(phi_est, nrow = 1, ncol = 1),
    Zt = matrix(1, nrow = 1, ncol = 1),
    HHt = matrix(evar_est, nrow = 1, ncol = 1),
    GGt = matrix(ivar_est, nrow = 1, ncol = 1),
    yt = matrix(ts, nrow = 1, ncol = n)
  )
  
  # Run smoother to get smoothed states
  ks_final <- fks(kf_final)
  smoothed_states <- as.numeric(ks_final$ahatt)
  
  # ========== Compute Measurement-Error-Free Residuals ==========
  measurement_residuals <- ts - mu_est - smoothed_states
  
  # ========== Return Results (SLIM VERSION) ==========
  return(list(
    residuals = measurement_residuals,
    mu = mu_est,
    ivar = ivar_est,
    evar = evar_est,
    heywood_ivar = heywood_ivar,
    heywood_evar = heywood_evar
  ))
}



fit_ar1wn_residuals_bayesian <- function(ts, verbose = FALSE, 
                                         n_burnin = 40000, 
                                         n_samples = 40000,
                                         n_chains = 3) {
  
  nt <- length(ts)
  
  model_code <- "
  model {
    for (t in 2:nt) {
      y[t] ~ dnorm(muy[t], Epre)
      muy[t] <- mu + ytilde[t]
      ytilde[t] ~ dnorm(muytilde[t], Ipre)
      muytilde[t] <- phi * ytilde[t-1]
    }
    ytilde[1] <- y[1] - mu
    Epre <- 1 / Evar
    Evar ~ dunif(0, 500)
    Ipre <- 1 / Ivar
    Ivar ~ dunif(0, 500)
    phi ~ dunif(-1, 1)
    mu ~ dnorm(50, 0.001)
  }
  "
  
  model_file <- tempfile(fileext = ".txt")
  writeLines(model_code, model_file)
  
  inits <- list(
    list(phi = 0.6, mu = 50, Ivar = 30, Evar = 20),
    list(phi = 0.3, mu = 70, Ivar = 60, Evar = 40),
    list(phi = -0.2, mu = 80, Ivar = 20, Evar = 30)
  )
  
  tryCatch({
    jags_model <- jags.model(
      file = model_file,
      data = list(y = ts, nt = nt),
      inits = inits,
      n.chains = n_chains,
      quiet = !verbose
    )
    
    if (verbose) cat("Running burnin...\n")
    update(jags_model, n.iter = n_burnin, progress.bar = ifelse(verbose, "text", "none"))
    
    # Sample mu, Ivar, and Evar (needed for residuals and Heywood check)
    if (verbose) cat("Sampling from posterior...\n")
    samples_params <- coda.samples(
      jags_model,
      variable.names = c("mu", "Ivar", "Evar", "ytilde"),
      n.iter = n_samples,
      progress.bar = ifelse(verbose, "text", "none")
    )
    
    # ========== Extract Only What You Need ==========
    samples_matrix <- as.matrix(samples_params)
    
    # Get posterior means
    mu_est <- mean(samples_matrix[, "mu"])
    ivar_est <- mean(samples_matrix[, "Ivar"])
    evar_est <- mean(samples_matrix[, "Evar"])
    
    # Check Heywood cases
    heywood_ivar <- ivar_est < 1e-6
    heywood_evar <- evar_est < 1e-6
    
    # Extract latent states
    ytilde_cols <- grep("^ytilde\\[", colnames(samples_matrix))
    ytilde_posterior_mean <- colMeans(samples_matrix[, ytilde_cols])
    
    # Compute residuals
    measurement_residuals <- ts - mu_est - ytilde_posterior_mean
    
    # Remove names to match frequentist format
    names(measurement_residuals) <- NULL
    
    # Clean up
    rm(samples_matrix)
    gc()
    
    # ========== Return Results ==========
    return(list(
      residuals = measurement_residuals,
      mu = mu_est,
      ivar = ivar_est,
      evar = evar_est,
      heywood_ivar = heywood_ivar,
      heywood_evar = heywood_evar
    ))
    
  }, error = function(e) {
    cat("Error in Bayesian AR(1)+WN fitting:", e$message, "\n")
    return(NULL)
  })
}
