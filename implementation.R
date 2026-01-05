## Gibbs posterior calibration and inference implementation for Huber Regression

source("functions.R")

# Data generating process
set.seed(2)

n <- 100
X <- cbind(1, runif(n, -2, 2))
p <- ncol(X)
alpha = 0.05
beta_true <- c(1, 2)
prop_cov <- diag(p) * 0.1

y <- X %*% beta_true + rt(n, df = 3)


# GPC-SA calibration
result_sa <- gpc(X, y, n, p, eta_init = 0.5, alpha, epsilon = 0.005, B = 500, 
                 n_iter = 2000, burnin = 1000, max_wp_iter =50, method = "SA")

result_sa$coverage
result_sa$eta
result_sa$s

# GPC-WP calibration
result_wp <- gpc(X, y, n, p, eta_init = 0.5, alpha, epsilon = 0.005, B = 500, 
                 n_iter = 2000, burnin = 1000, max_wp_iter =50, method = "WP")

# Plot Coverage over iterations
xlim <- c(5,10,15,20,25,30,35,40,45,50)
coverage_list <- c(0.85, 0.854, 0.848, 0.852, 0.846, 0.848, 0.85, 0.852, 0.85, 0.85)

plot(xlim, coverage_list, type = 'b', pch = 19,
     xlab = "Iteration", ylab = "Coverage", ylim = c(0.845,1) ,
     main = "GPC-WP Coverage Evolution")
abline(h = 0.95, col = 'red', lty = 2)

# RW-MH inference
eta_hat <- result_sa$eta

fit <- mh_sampler(    
  X, 
  y, 
  eta = eta_hat,
  beta_init = c(0, 0),
  prop_cov = prop_cov,
  n_iter = 5000,
  burnin = 1000
)

fit$mean        # posterior mean estimate
fit$cov         # posterior covariance
fit$accept_ratio # acceptance rate
ess_mcmc <- ess_function(fit$samples)
mcse_mcmc <- mcse_function(fit$samples)



# pCN inference
fit_pcn <- posterior_huber_pcn(
  X = X,
  y = y,
  delta = 1.5,
  eta = eta_hat,
  theta_init = rep(0, p),
  prior_mean = rep(0, p),
  prior_cov = 10 * diag(p),  # tuned for scale
  rho = 0.14,
  n_iter = 5000,
  burnin = 1000
)

fit_pcn$mean        # posterior mean estimate
fit_pcn$cov         # posterior covariance
fit_pcn$accept_ratio # acceptance rate
ess_pcn <- ess_function(fit_pcn$samples)
mcse_pcn <- mcse_function(fit_pcn$samples)





## GPC table
gpc_df <- data.frame(
  Method = c("GPC-SA", "GPC-WP"),
  "$\\eta$" = c(result_sa$eta, result_wp$eta),
  Coverage = c(result_sa$coverage, result_wp$coverage),
  check.names = FALSE
)

print(xtable::xtable(
                    gpc_df, 
                    caption = "Learning rates computed with different methods.", 
                    label = "tab:results"), 
      file = "gpc_df.tex", 
      include.rownames = FALSE,
      sanitize.text.function = function(x) {x},
      booktabs = TRUE,
      comment = FALSE)


## Inference table
inference_df <- data.frame(
  Method = c("RW-MH", "pCN"),
  "$\\beta_{0}$" = c(fit$mean[1], fit_pcn$mean[1]),
  "$\\beta_{1}$" = c(fit$mean[2], fit_pcn$mean[2]),
  "MCSE $\\beta_{0}$" = c(mcse_mcmc[1],mcse_pcn[2]),
  "MCSE $\\beta_{1}$" = c(mcse_mcmc[1],mcse_pcn[2]),
  "Acc. Rat." = c(fit$accept_ratio,fit_pcn$accept_ratio),
  ESS = c(ess_mcmc,ess_pcn),

  check.names = FALSE
)

print(xtable::xtable(
                      inference_df, 
                      caption = "Inference results", 
                      label = "tab:inference_results"), 
  file = "inference_df.tex", 
  include.rownames = FALSE,
  sanitize.text.function = function(x) {x},
  booktabs = TRUE,
  comment = FALSE)


# ACF Plots
par(mfrow = c(2, 2), mar = c(3, 4, 3, 1))
acf(fit$samples[, 1], lag.max = 20,
    main = bquote("RW-MH" ~ beta[0])
)
acf(fit$samples[, 2], lag.max = 20,
    main = bquote("RW-MH" ~ beta[1])
)
acf(fit_pcn$samples[, 1], lag.max = 20,
    main = bquote("pCN" ~ beta[0])
)
acf(fit_pcn$samples[, 2], lag.max = 20,
    main = bquote("pCN" ~ beta[1])
)
grid()


# Histogram posteriors
par(mfrow = c(2, 2), mar = c(3, 4, 3, 1))
hist(fit$samples[,1], main = bquote("RW-MH" ~ beta[0]))
hist(fit$samples[,2], main = bquote("RW-MH" ~ beta[1]))

hist(fit_pcn$samples[,1],main = bquote("pCN" ~ beta[0]))
hist(fit_pcn$samples[,2],main = bquote("pCN" ~ beta[1]))
grid()

# Trace Plots
par(mfrow = c(2, 2), mar = c(3, 4, 3, 1))
plot(
  fit$samples[, 1],
  type = "l",
  main = bquote("RW-MH" ~ beta[0]),
  xlab = "Iteration",
  ylab = expression(beta[0]))
plot(
  fit$samples[, 2],
  type = "l",
  main = bquote("RW-MH" ~ beta[1]),
  xlab = "Iteration",
  ylab = expression(beta[1]))
plot(
  fit_pcn$samples[, 1],
  type = "l",
  main = bquote("pCN" ~ beta[0]),
  xlab = "Iteration",
  ylab = expression(beta[0]))
plot(
  fit_pcn$samples[, 2],
  type = "l",
  main = bquote("pCN" ~ beta[1]),
  xlab = "Iteration",
  ylab = expression(beta[1]))