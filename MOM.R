###############################################################
# Monte Carlo Simulation of Method-of-Moments Estimators
# for Asset Correlation in a One-Factor Credit Risk Model
#
# This script simulates credit defaults under a one-factor
# Gaussian latent variable model and estimates the asset
# correlation parameter using:
#   1. The original Method-of-Moments (MoM) estimator
#   2. Several bias-adjusted MoM estimators
#
# Key model inputs are defined at the top of the script so they
# can be changed easily:
#   - Unconditional default probability (pd)
#   - Asset correlation loading (w)
#   - Number of time periods (t)
#   - Number of obligors (n)
#   - Number of simulation runs (n_sim)
#
# For each simulated dataset, the script:
#   1. Simulates latent asset returns
#   2. Determines defaults using a threshold implied by pd
#   3. Computes period-specific default rates
#   4. Estimates the correlation parameter
#   5. Computes adjusted estimators
#   6. Summarizes and plots their distributions
###############################################################

############################
# Model Parameters
############################
pd <- 0.002              # Unconditional default probability
w <- 0.20                # Asset correlation loading
t <- 20                  # Number of time periods
n <- 2000               # Number of obligors
n_sim <- 1000            # Number of Monte Carlo simulation runs

# Default threshold implied by the unconditional default probability
gamma_true <- qnorm(pd)

############################
# Required Packages
############################
library(mvtnorm)
library(ggplot2)
library(plyr)

############################
# Functions for MoM Estimation
############################

# Equation used to recover rho from the bivariate normal model.
# This is solved numerically by uniroot().
binorm <- function(rho) {
  p <- E_Z
  sigma <- matrix(c(1, rho, rho, 1), nrow = 2)
  return(pmvnorm(lower = -Inf, upper = c(a, a), mean = c(0, 0), sigma) - p)
}

# First derivative of the bivariate normal CDF with respect to rho.
# This enters the bias-adjustment formula.
phi_prime <- function(s, t_val, rho) {
  y <- (1 / ((2 * pi) * sqrt(1 - rho^2))) *
    exp(-(s^2 / 2 - rho * s * t_val + t_val^2 / 2) / (1 - rho^2))
  return(y)
}

# Second derivative of the bivariate normal CDF with respect to rho.
# This is also needed for the bias-adjustment formula.
phi_prime_prime <- function(s, t_val, rho) {
  y <- ((s * t_val + rho * (1 - s^2 - t_val^2) + s * t_val * rho^2 - rho^3) /
          (2 * pi * (1 - rho^2)^(5 / 2))) *
    exp(-(s^2 / 2 - rho * s * t_val + t_val^2 / 2) / (1 - rho^2))
  return(y)
}

# Compute several bias-adjusted estimators of rho using
# autocovariances up to lag 5.
adjusted_rho <- function(rho, g_prime, g_prime_prime,
                         alpha_0, alpha_1, alpha_2,
                         alpha_3, alpha_4, alpha_5, t) {
  
  adj_rho <- c()
  
  adj_rho[1] <- rho + (g_prime_prime / (2 * t * (g_prime)^3)) * alpha_0 +
    (g_prime_prime / (t * (g_prime)^3) * (1 - 1 / t) * alpha_1)
  
  adj_rho[2] <- rho + g_prime_prime / (t * (g_prime)^3) *
    (alpha_0 / 2 + (1 - 1 / t) * alpha_1 + (1 - 2 / t) * alpha_2)
  
  adj_rho[3] <- rho + g_prime_prime / (t * (g_prime)^3) *
    (alpha_0 / 2 + (1 - 1 / t) * alpha_1 + (1 - 2 / t) * alpha_2 +
       (1 - 3 / t) * alpha_3)
  
  adj_rho[4] <- rho + g_prime_prime / (t * (g_prime)^3) *
    (alpha_0 / 2 + (1 - 1 / t) * alpha_1 + (1 - 2 / t) * alpha_2 +
       (1 - 3 / t) * alpha_3 + (1 - 4 / t) * alpha_4)
  
  adj_rho[5] <- rho + g_prime_prime / (t * (g_prime)^3) *
    (alpha_0 / 2 + (1 - 1 / t) * alpha_1 + (1 - 2 / t) * alpha_2 +
       (1 - 3 / t) * alpha_3 + (1 - 4 / t) * alpha_4 +
       (1 - 5 / t) * alpha_5)
  
  return(adj_rho)
}

############################
# Storage Objects
############################

# Matrix of latent returns: rows = obligors, columns = time periods
y_a <- matrix(, nrow = n, ncol = t)

# Matrix of default indicators:
# 1 = default, 0 = no default
default_a <- matrix(, nrow = n, ncol = t)

# Number of defaults in each time period
default_time_a <- numeric(t)

# Default rate in each time period
p_default_a <- numeric(t)

# Storage for estimators:
# Column 1 = original MoM estimator
# Columns 2 to 6 = adjusted estimators for k = 1, ..., 5
result_MoM <- matrix(, nrow = n_sim, ncol = 6)

############################
# Monte Carlo Simulation
############################

times <- 1

while (times < (n_sim + 1)) {
  
  # Uncomment for reproducible results
  # set.seed(times)
  
  ############################
  # Step 1: Simulate systematic factor
  ############################
  # One common risk factor for each time period
  x <- rnorm(t)
  
  ############################
  # Step 2: Simulate latent returns
  ############################
  # For each obligor:
  #   y_it = w * x_t + sqrt(1 - w^2) * epsilon_it
  # where x_t is the systematic factor and epsilon_it is
  # the idiosyncratic risk term.
  for (i in 1:n) {
    epsilon <- rnorm(t)
    y <- x * w + epsilon * sqrt(1 - w^2)
    y_a[i, ] <- y
  }
  
  ############################
  # Step 3: Determine defaults
  ############################
  # A default occurs if the latent return falls below
  # the threshold implied by the unconditional PD.
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
  # Step 4: Aggregate defaults by time period
  ############################
  for (i in 1:ncol(default_a)) {
    default_time_a[i] <- sum(default_a[, i])
    p_default_a[i] <- mean(default_a[, i])
  }
  
  ############################
  # Step 5: Construct MoM estimator
  ############################
  
  # Auxiliary vector of ones, used for averaging
  aux <- rep(1, t)
  
  # Z_t = p_t^2, where p_t is the observed default rate
  Z <- p_default_a * p_default_a
  
  # Sample mean of Z_t
  E_Z <- mean(Z)
  
  # Threshold implied by the average default rate
  a <- qnorm(crossprod(p_default_a, aux) / t)
  
  ############################
  # Step 6: Compute variance and autocovariances of Z_t
  ############################
  
  # Variance (lag 0 covariance)
  var_1 <- var(Z)
  
  # Sample autocovariances up to lag 5
  cov_1 <- sum((Z[1:(t - 1)] - E_Z) * (Z[2:t] - E_Z)) / t
  cov_2 <- sum((Z[1:(t - 2)] - E_Z) * (Z[3:t] - E_Z)) / t
  cov_3 <- sum((Z[1:(t - 3)] - E_Z) * (Z[4:t] - E_Z)) / t
  cov_4 <- sum((Z[1:(t - 4)] - E_Z) * (Z[5:t] - E_Z)) / t
  cov_5 <- sum((Z[1:(t - 5)] - E_Z) * (Z[6:t] - E_Z)) / t
  
  ############################
  # Step 7: Estimate rho and convert to w
  ############################
  
  # Recover rho from the bivariate normal equation
  hat_rho <- uniroot(binorm, c(-1, 1))$root
  
  # Since rho = w^2 in this setup, recover w as sqrt(rho)
  hat_w <- sqrt(hat_rho)
  
  ############################
  # Step 8: Compute adjusted estimators
  ############################
  
  ad_rho <- adjusted_rho(
    hat_rho,
    phi_prime(a, a, hat_rho),
    phi_prime_prime(a, a, hat_rho),
    var_1, cov_1, cov_2, cov_3, cov_4, cov_5, t
  )
  
  # Convert adjusted rho values into adjusted w values
  ad_w <- sqrt(ad_rho)
  
  ############################
  # Step 9: Store results
  ############################
  result_MoM[times, ] <- c(hat_w, ad_w)
  
  times <- times + 1
}

############################
# Summary Statistics
############################

# Mean of each estimator across simulations
mean3 <- c(mean(result_MoM[, 1]), mean(result_MoM[, 2]), mean(result_MoM[, 3]),
           mean(result_MoM[, 4]), mean(result_MoM[, 5]), mean(result_MoM[, 6]))

# Standard deviation of each estimator across simulations
sd3 <- c(sd(result_MoM[, 1]), sd(result_MoM[, 2]), sd(result_MoM[, 3]),
         sd(result_MoM[, 4]), sd(result_MoM[, 5]), sd(result_MoM[, 6]))

############################
# Prepare Data for Plotting
############################

simu0 <- as.data.frame(cbind(result_MoM[, 1], rep("original", n_sim)))
simu1 <- as.data.frame(cbind(result_MoM[, 2], rep("k=1", n_sim)))
simu2 <- as.data.frame(cbind(result_MoM[, 3], rep("k=2", n_sim)))
simu3 <- as.data.frame(cbind(result_MoM[, 4], rep("k=3", n_sim)))
simu4 <- as.data.frame(cbind(result_MoM[, 5], rep("k=4", n_sim)))
simu5 <- as.data.frame(cbind(result_MoM[, 6], rep("k=5", n_sim)))

simu <- rbind(simu0, simu1, simu2, simu3, simu4, simu5)
colnames(simu) <- c("estimators", "indicators")

# Convert to appropriate types
simu$estimators <- as.numeric(simu$estimators)
simu$indicator <- factor(simu$indicators)

############################
# Plot Density of Estimators
############################

# Mean of each estimator group, used for dashed vertical lines
mu <- ddply(simu, "indicator", summarise, grp.mean = mean(estimators))

p <- ggplot(simu, aes(x = estimators, color = indicator)) +
  geom_density() +
  xlim(0, 0.4) +
  geom_vline(data = mu,
             aes(xintercept = grp.mean, color = indicator),
             linetype = "dashed")

p + theme_minimal()