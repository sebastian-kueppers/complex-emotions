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
  convergence_info <- list()
  
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
      hessian = TRUE,
      control = list(maxit = 1000)
    )
    
    convergence_info[[i]] <- list(
      start = i,
      convergence = fit$convergence,
      loglik = -fit$value
    )
    
    # Keep best fit
    if (fit$value < best_loglik) {
      best_loglik <- fit$value
      best_fit <- fit
    }
  }
  
  if (verbose) {
    best_start <- which.min(sapply(convergence_info, function(x) x$loglik))
    cat("\nBest fit from starting value:", best_start, "\n")
    cat("Convergence code:", best_fit$convergence, "\n")
    cat("Log-likelihood:", -best_fit$value, "\n\n")
  }
  
  # ========== Extract Parameter Estimates ==========
  phi_est <- best_fit$par[1]
  mu_est <- best_fit$par[2]
  ivar_est <- exp(best_fit$par[3])
  evar_est <- exp(best_fit$par[4])
  
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
  
  # ========== Compute Measurement Residuals ==========
  # Residuals = Observed - (Mean + Smoothed Latent State)
  # These represent the measurement error component
  measurement_residuals <- ts - mu_est - smoothed_states
  
  # ========== Return Results ==========
  return(list(
    residuals = measurement_residuals,
    phi = phi_est,
    mu = mu_est,
    ivar = ivar_est,
    evar = evar_est,
    loglik = -best_fit$value,
    convergence = best_fit$convergence,
    heywood_ivar = heywood_ivar,
    heywood_evar = heywood_evar,
    smoothed_states = smoothed_states,
    fitted_values = mu_est + smoothed_states,
    best_starting_value = which.min(sapply(convergence_info, function(x) x$loglik))
  ))
}



fit_ar1wn_residuals_bayesian <- function(ts, verbose = FALSE, 
                                         n_burnin = 40000, 
                                         n_samples = 40000,
                                         n_chains = 3) {
  
  # ========== Prepare Data ==========
  nt <- length(ts)
  
  # ========== JAGS Model ==========
  model_code <- "
  model {
    # Likelihood
    for (t in 2:nt) {
      y[t] ~ dnorm(muy[t], Epre)
      muy[t] <- mu + ytilde[t]
      ytilde[t] ~ dnorm(muytilde[t], Ipre)
      muytilde[t] <- phi * ytilde[t-1]
    }
    
    # First time point
    ytilde[1] <- y[1] - mu
    
    # Priors (matching the paper)
    Epre <- 1 / Evar
    Evar ~ dunif(0, 500)
    
    Ipre <- 1 / Ivar
    Ivar ~ dunif(0, 500)
    
    phi ~ dunif(-1, 1)
    mu ~ dnorm(50, 0.001)
  }
  "
  
  # Write model to temporary file
  model_file <- tempfile(fileext = ".txt")
  writeLines(model_code, model_file)
  
  # ========== Initial Values ==========
  ts_mean <- mean(ts, na.rm = TRUE)
  ts_var <- var(ts, na.rm = TRUE)
  
  inits <- list(
    list(phi = 0.2, mu = ts_mean, Ivar = ts_var/3, Evar = ts_var/3),
    list(phi = 0.5, mu = ts_mean, Ivar = ts_var/2, Evar = ts_var/2),
    list(phi = 0.8, mu = ts_mean, Ivar = ts_var/1.5, Evar = ts_var/1.5)
  )
  
  # ========== Fit Model ==========
  tryCatch({
    jags_model <- jags.model(
      file = model_file,
      data = list(y = ts, nt = nt),
      inits = inits,
      n.chains = n_chains,
      quiet = !verbose
    )
    
    # Burnin
    if (verbose) cat("Running burnin...\n")
    update(jags_model, n.iter = n_burnin, progress.bar = ifelse(verbose, "text", "none"))
    
    # Sample from posterior (including latent states ytilde)
    if (verbose) cat("Sampling from posterior...\n")
    samples <- coda.samples(
      jags_model,
      variable.names = c("phi", "mu", "Ivar", "Evar", "ytilde"),
      n.iter = n_samples,
      progress.bar = ifelse(verbose, "text", "none")
    )
    
    # ========== Extract Posterior Summaries ==========
    # Convert to matrix for easier extraction
    samples_matrix <- as.matrix(samples)
    
    # Posterior means 
    phi_est <- mean(samples_matrix[, "phi"])
    mu_est <- mean(samples_matrix[, "mu"])
    ivar_est <- mean(samples_matrix[, "Ivar"])
    evar_est <- mean(samples_matrix[, "Evar"])
    
    # Posterior standard deviations
    phi_sd <- sd(samples_matrix[, "phi"])
    mu_sd <- sd(samples_matrix[, "mu"])
    ivar_sd <- sd(samples_matrix[, "Ivar"])
    evar_sd <- sd(samples_matrix[, "Evar"])
    
    # Credible intervals
    phi_ci <- quantile(samples_matrix[, "phi"], c(0.025, 0.975))
    mu_ci <- quantile(samples_matrix[, "mu"], c(0.025, 0.975))
    ivar_ci <- quantile(samples_matrix[, "Ivar"], c(0.025, 0.975))
    evar_ci <- quantile(samples_matrix[, "Evar"], c(0.025, 0.975))
    
    # ========== Check for Heywood Cases ==========
    heywood_ivar <- ivar_est < 1e-6
    heywood_evar <- evar_est < 1e-6
    
    # ========== Extract Latent States (ytilde) ==========
    # Find columns corresponding to ytilde
    ytilde_cols <- grep("^ytilde\\[", colnames(samples_matrix))
    
    # Posterior mean of latent states
    ytilde_posterior_mean <- colMeans(samples_matrix[, ytilde_cols])
    
    # Posterior SD of latent states (uncertainty in estimates)
    ytilde_posterior_sd <- apply(samples_matrix[, ytilde_cols], 2, sd)
    
    # ========== Compute Measurement-Error-Free Residuals ==========
    # Residuals = Observed - (Mean + Latent State)
    measurement_residuals <- ts - mu_est - ytilde_posterior_mean
    
    # ========== Fitted Values ==========
    fitted_values <- mu_est + ytilde_posterior_mean
    
    # ========== Smoothed States ==========
    smoothed_states <- ytilde_posterior_mean
    
    # ========== Convergence Diagnostics ==========
    # Gelman-Rubin statistic (should be < 1.05 for convergence)
    params_only <- samples[, c("phi", "mu", "Ivar", "Evar"), drop = FALSE]
    gelman_diag <- gelman.diag(params_only, confidence = 0.95)
    gelman_psrf <- gelman_diag$psrf[, "Point est."]
    
    # Check if all parameters converged
    all_converged <- all(gelman_psrf < 1.05, na.rm = TRUE)
    
    # ========== Return Results ==========
    return(list(
      # Measurement-error-free residuals (main output)
      residuals = measurement_residuals,
      
      # Parameter estimates (posterior means)
      phi = phi_est,
      mu = mu_est,
      ivar = ivar_est,
      evar = evar_est,
      
      # Parameter uncertainty (posterior SDs)
      phi_sd = phi_sd,
      mu_sd = mu_sd,
      ivar_sd = ivar_sd,
      evar_sd = evar_sd,
      
      # Credible intervals
      phi_ci = phi_ci,
      mu_ci = mu_ci,
      ivar_ci = ivar_ci,
      evar_ci = evar_ci,
      
      # Heywood cases
      heywood_ivar = heywood_ivar,
      heywood_evar = heywood_evar,
      
      # Latent states and fitted values
      smoothed_states = smoothed_states,
      fitted_values = fitted_values,
      ytilde_posterior_mean = ytilde_posterior_mean,
      ytilde_posterior_sd = ytilde_posterior_sd,
      
      # Convergence information
      gelman_psrf = gelman_psrf,
      all_converged = all_converged,
      
      # MCMC samples (for further analysis if needed)
      samples = samples,
      samples_matrix = samples_matrix,
      
      # Metadata
      n_burnin = n_burnin,
      n_samples = n_samples,
      n_chains = n_chains,
      n_obs = nt
    ))
    
  }, error = function(e) {
    cat("Error in Bayesian AR(1)+WN fitting:", e$message, "\n")
    return(NULL)
  })
}

