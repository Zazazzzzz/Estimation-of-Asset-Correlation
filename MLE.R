###############################################################
# Monte Carlo Simulation of Default Clustering with MLE Estimation
#
# This script simulates 1,000 datasets of credit portfolios
# under a one-factor Gaussian latent variable model.
#
# Key parameters (can be modified to explore different scenarios):
#   - Unconditional default probability: 0.002 (via threshold)
#   - Asset correlation (w): 0.15
#   - Time periods (t): 20
#   - Number of obligors (n): 500
#
# The script:
#   1. Simulates latent returns
#   2. Determines defaults via a threshold
#   3. Aggregates default counts over time
#   4. Estimates parameters (w, gamma) via MLE
#   5. Analyzes distribution of estimators
#
# Changing parameters (w, n, t, threshold) will lead to
# different simulation and estimation results.
###############################################################

library(ggplot2)

############################
# Model Parameters
############################
w <- 0.15                  # Asset correlation parameter
n <- 500                  # Number of obligors
t <- 20                   # Number of time periods
n_sim <- 1000             # Number of simulation runs

# Default threshold corresponding to unconditional PD = 0.002
gamma_true <- -2.878162   # qnorm(0.002)

############################
# Storage Initialization
############################

# Matrix to store simulated latent returns (n x t)
y_a <- matrix(, nrow = n, ncol = t)

# Matrix to store default indicators (1 = default, 0 = no default)
default_a <- matrix(, nrow = n, ncol = t)

# Vector to store number of defaults per time period
default_time_a <- c()

# Matrix to store MLE estimates: [w_hat, gamma_hat]
result_MLE <- matrix(, nrow = n_sim, ncol = 2)

############################
# Monte Carlo Simulation Loop
############################

times <- 1

while (times < (n_sim + 1)) {
  
  set.seed(times)  
  
  # Systematic risk factor (common across obligors)
  x <- rnorm(t)
  
  ############################
  # Step 1: Simulate Latent Returns
  ############################
  for (i in 1:n) {
    epsilon <- rnorm(t)  # Idiosyncratic risk
    y <- x * w + epsilon * sqrt(1 - w^2)
    y_a[i, ] <- y
  }
  
  ############################
  # Step 2: Determine Defaults
  ############################
  for (i in 1:nrow(y_a)) {
    for (j in 1:ncol(y_a)) {
      if (y_a[i, j] < gamma_true) {
        default_a[i, j] <- 1
      } else {
        default_a[i, j] <- 0
      }
    }
  }
  
  ############################
  # Step 3: Aggregate Defaults per Time Period
  ############################
  for (i in 1:ncol(default_a)) {
    default_time_a[i] <- sum(default_a[, i])
  }
  
  ############################
  # Step 4: Maximum Likelihood Estimation
  ############################
  logll <- function(hat_w) {
    
    F <- c()
    w <- hat_w[1]
    gamma <- hat_w[2]
    
    for (i in 1:length(default_time_a)) {
      
      f <- function(u) {
        # Conditional default probability given systematic factor u
        p <- pnorm((gamma - w * u) / sqrt(1 - w^2))
        
        # Likelihood contribution (Binomial × Normal density)
        dbinom(default_time_a[i], nrow(default_a), p) * dnorm(u)
      }
      
      # Integrate over systematic factor
      F[i] <- integrate(f, lower = -5, upper = 10)$value
    }
    
    # Negative log-likelihood
    return(-sum(log(F)))
  }
  
  # Optimization (bounded MLE)
  re_MLE <- optim(
    c(0.1, -3.5),
    logll,
    method = "L-BFGS-B",
    lower = c(0.1, -3.5),
    upper = c(0.99, -1)
  )$par
  
  # Store results
  result_MLE[times, ] <- re_MLE
  
  times <- times + 1
}

############################
# Results Analysis
############################

result1 <- result_MLE

# Summary statistics
mean(result1[, 2])   # Mean of gamma estimates
sd(result1[, 2])     # Std. dev. of gamma estimates
mean(result1[, 1])   # Mean of w estimates

############################
# Visualization: w Estimates
############################

simu1_of <- data.frame(result1[, 1])
colnames(simu1_of) <- c("estimators")

p <- ggplot(simu1_of, aes(x = estimators)) +
  geom_density() +
  xlim(0, 0.5) +
  geom_vline(aes(xintercept = 0.1516684), linetype = "dashed")

p + theme_minimal()

############################
# Visualization: gamma Estimates
############################

simu1_of_gamma <- data.frame(result1[, 2])
colnames(simu1_of_gamma) <- c("estimators")

p <- ggplot(simu1_of_gamma, aes(x = estimators)) +
  geom_density() +
  xlim(-3.2, -2.5) +
  geom_vline(aes(xintercept = -2.886564), linetype = "dashed")

p + theme_minimal()

