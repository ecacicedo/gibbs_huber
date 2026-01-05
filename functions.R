### Functions for Weighted Particle-Based Gibbs Posterior learning rate 
### calibration based on arXiv:2405.04845v4 and
### Huber Regression Gibbs posterior inference with RWMH and pCN

## Auxiliary functions
# Huber Loss
huber_loss <- function(r, delta = 1.345) {
  ifelse(abs(r) <= delta, 0.5 * r^2,
         delta * abs(r) - 0.5 * delta^2)
}

# Cumulative Loss
loss_cum <- function(beta, X, y) {
  r <- as.vector(y - X %*% beta)
  sum(huber_loss(r))
}

# Log Prior distribution
log_prior <- function(theta)
  sum(dnorm(theta, 0, 10, log = TRUE))

# Log Target distribution
log_target <- function(beta, eta, X, y) {
  log_prior(beta) - eta * loss_cum(beta, X, y)
}

# Random Walk Metropolis-Hastings sampler

mh_sampler <- function(
    X, 
    y, 
    eta,
    beta_init,
    prop_cov,
    n_iter,
    burnin
) {
  p <- length(beta_init)
  samples <- matrix(NA, n_iter, p)
  weights <- matrix(NA, n_iter, 1)
  
  beta_curr <- beta_init
  log_curr  <- log_target(beta_curr, eta, X, y)
  
  accept <- 0
  n_total <- n_iter + burnin
  
  for (t in 1:n_total) {
    
    beta_prop <- MASS::mvrnorm(1, beta_curr, prop_cov)
    log_prop  <- log_target(beta_prop, eta, X, y)
    
    accepted <- FALSE
    if (log(runif(1)) < log_prop - log_curr) {
      beta_curr <- beta_prop
      log_curr  <- log_prop
      accepted <- TRUE
    }
    
    if (t > burnin) {
      samples[t - burnin, ] <- beta_curr
      weights[t - burnin, ] <- log_curr
      if (accepted) accept <- accept + 1
    }
  }
  
  list(
    weights = weights/sum(weights),
    samples = samples,
    accept_ratio = accept / n_iter,
    mean = colMeans(samples),
    cov = cov(samples)
  )
}

# Empirical coverage for MCMC samples
coverage_estimator <- function(samples, beta_true, alpha) {
  ci <- apply(samples, 2, quantile,
              probs = c(alpha/2, 1 - alpha/2))
  as.numeric(all(ci[1, ] <= beta_true & beta_true <= ci[2, ]))
}

## Calibration Functions
# GPC (Figure 1)

gpc <- function(
    X, 
    y, 
    n, 
    p, 
    eta_init = 1, 
    alpha, 
    epsilon = 0.005, 
    B = 500, 
    n_iter = 2000, 
    burnin = 1000,
    max_wp_iter = 50,
    method
){
  # Bootstrap datasets
  boot_idx <- replicate(B, sample(1:n, n, replace = TRUE), simplify = FALSE)
  
  # Initializer
  s <- 1
  l <- 1
  eta <- eta_init
  eta_count <- rep(eta_init, 2)
  beta_init <- numeric(p)
  
  while (TRUE){
    cover <- numeric(B)
    particles_list <- vector("list", B)
    weights_list <- vector("list", B)
    
    for (b in 1:B) {
      
      Xb <- X[boot_idx[[b]], ]
      yb <- y[boot_idx[[b]]]
      
      mcmc <- mh_sampler(
        Xb, 
        yb, 
        eta,
        beta_init,
        prop_cov, 
        n_iter,
        burnin
      )
      
      particles_list[[b]] <- mcmc$samples
      weights_list[[b]] <- mcmc$weights
      cover[b] <- coverage_estimator(mcmc$samples, beta_true, alpha)
    }
    
    coverage <- mean(cover)
    
    if (abs(coverage - (1 - alpha)) < epsilon){
      break
    }
    
    else if (method == "SA"){
      eta_count <- base::append(eta_count, eta)
      optimize <- robbins_munro(eta_count, coverage, alpha, l)
      eta_vec <- optimize$eta_vec
      eta <- eta_vec[s + 3]
      l <- optimize$l
      s <- s + 1
    }
    
    else if (method == "WP"){
      wb <- wb_optimize(
        eta_count,
        eta,
        coverage,
        epsilon = 0.005,
        s,
        B,
        particles_list,
        weights_list,
        n_iter,
        max_wp_iter,
      )
      eta <- wb$eta
      coverage <- wb$coverage
      eta_list <- wb$eta_list
      coverage_list <- wb$coverage_list
      min_ess_list <- wb$min_ess_list
      break
    }
  }
  
  list(
    eta = eta,
    coverage = coverage,
    s = s,
    eta_list = eta_list,
    coverage_list = coverage_list,
    min_ess_list = min_ess_list
  )
  
}
  
# Stochastic Approximation with Keston's rule (Figure 2)
robbins_munro <- function(eta_vec, coverage, alpha, l){
  i = length(eta_vec)
  if ((eta_vec[i-1] - eta_vec[i-2]) * (eta_vec[i] - eta_vec[i-1]) <0){
      l <- l+1
      }
      eta_new <- eta_vec[i] + (l^(-0.51)) * (coverage - (1 - alpha))
      eta_vec <- base::append(eta_vec, eta_new)
      
    list(
      eta_vec = eta_vec,
      l = l
    )
}

# Weighted particle-based optimization (Figure 3)
wb_optimize <- function(
    eta_count, 
    eta,
    coverage, 
    epsilon = 0.005, 
    s, 
    B, 
    mcmc_samples, 
    mcmc_weights, 
    n_iter,
    max_wp_iter,
    eta_list,
    coverage_list,
    min_ess_list
){
  u <- 1
  l <- 1
  tau <- 0
  ess_threshold = 0.25 * n_iter
  
  wp_iter <- 1
  
  while (wp_iter <= max_wp_iter){
    eta_count <- base::append(eta_count, eta)
    optimize <- robbins_munro(eta_count, coverage, alpha, l)
    eta_vec <- optimize$eta_vec
    eta <- eta_vec[s + 3]
    eta_prev <- eta_vec[s + 2]
    l <- optimize$l
    
    cover <- numeric(B)
    ess <- numeric(B)
    
    for (b in 1:B) {  
    
    y_matrix <- matrix(y, nrow = 100, ncol = 2000)

    resid <- y_matrix - (X %*% t(mcmc_samples[[b]]))
    
    out <- ifelse(abs(resid) <= 1.345, 0.5 * resid^2, 
                  1.345 * abs(resid) - 0.5 * 1.345^2)
    attributes(out) <- attributes(resid) 
    
    loss_matrix <- matrix(out, nrow = 100, ncol = 2000)
    loss_vals <- colSums(loss_matrix)
    
    logw <- log(mcmc_weights[[b]]) - (eta - eta_prev) * loss_vals
    logw <- logw - max(logw)
    weights <- exp(logw)
    weights <- weights/sum(weights)
    
    ess[b] <- 1 / sum(weights^2)
    
    indices <- sample(1:2000, size = 2000, replace = TRUE, prob = weights)
    resample <- mcmc_samples[[b]][indices, ]
    
    
    cover[b] <- coverage_estimator(resample, beta_true, alpha)
    }
    
    coverage <- mean(cover)
    min_ess <- min(ess)
    
    if (abs(coverage - (1 - alpha)) < epsilon
        && min_ess >= ess_threshold){
      tau <- 1
      break
    }   
    
    else if (min_ess < ess_threshold){
      eta_count <- base::append(eta_count, eta)
      optimize <- robbins_munro(eta_count, coverage, alpha, l)
      eta_vec <- optimize$eta_vec
      eta <- eta_vec[s + 3]
    }
    
    else {u <- u + 1}
    
    if (wp_iter %% 5 == 0) {
      cat("WP iter:", wp_iter,
          "eta:", round(eta, 4),
          "coverage:", round(coverage, 4),
          "min ESS:", round(min_ess), "\n")
    }
    
    wp_iter <- wp_iter + 1
    
  }
  
  warning("WP optimizer stopped due to max_wp_iter")
  
  list(
    eta = eta,
    tau = tau,
    coverage = coverage,
    coverage_list = coverage_list,
    eta_list = eta_list,
    min_ess_list = min_ess_list
  )
}

## Inference Functions
# pCN sampler
posterior_huber_pcn <- function(
    X, 
    y, 
    delta, 
    eta, 
    theta_init,
    prior_mean, 
    prior_cov, 
    rho,
    n_iter = 5000, 
    burnin = 1000
) {
  
  p <- length(theta_init)
  
  samples <- matrix(NA, n_iter, p)
  beta <- theta_init
  accept <- 0
  n_total <- n_iter + burnin
  
  for (t in 1:n_total) {
    
    xi <- MASS::mvrnorm(1, mu = rep(0, p), Sigma = prior_cov)
    
    beta_prop <- prior_mean +
      sqrt(1 - rho^2) * (beta - prior_mean) +
      rho * xi
    
    log_alpha <- -eta * (loss_cum(beta_prop, X, y) - loss_cum(beta, X, y))
    
    if (log(runif(1)) < log_alpha) {
      beta <- beta_prop
      if (t > burnin) accept <- accept + 1
    }
    
    if (t > burnin) {
      samples[t - burnin, ] <- beta
    }
  }
  
  list(
    samples = samples,
    mean = colMeans(samples),
    cov = cov(samples),
    accept_ratio = accept / n_iter
  )
}


## Diagnostic Functions
ess_function <- function(x, max_lag = 50) {
  acf_vals <- acf(x, plot = FALSE, lag.max = max_lag)$acf[-1]
  pos <- which(acf_vals > 0)
  if (length(pos) == 0) return(length(x))
  length(x) / (1 + 2 * sum(acf_vals[pos]))
}

mcse_function <- function(samples) {
  p <- ncol(samples)
  se <- numeric(p)
  
  for (j in 1:p) {
    ess_j <- ess_function(samples[, j])
    var_j <- var(samples[, j])
    se[j] <- sqrt(var_j / ess_j)
  }
  se
}