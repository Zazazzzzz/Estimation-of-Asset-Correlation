###############################################################
# Estimation of Asset Correlation and Capital Requirements
# Using IBRD Default Data
#
# This script reads default data from an Excel file and applies
# several approaches to estimate asset correlation in a
# one-factor credit risk framework.
#
# The script includes:
#   1. Data preparation from the IBRD Excel file
#   2. Hand-written Method-of-Moments (MM) estimation
#   3. Hand-written Maximum Likelihood Estimation (MLE)
#   4. A closed-form MLE approximation
#   5. Estimation using functions from the AssetCorr package
#   6. Basel-style correlation calibration
#   7. Basel-style unexpected loss calculations
#
# Main outputs:
#   - Original and adjusted MM estimates
#   - MLE estimate
#   - Alternative package-based estimates
#   - Correlation estimates implied by Basel calibration
#   - Required capital / unexpected loss for different LGD values
###############################################################

############################
# Required Packages
############################
library(openxlsx)
library(AssetCorr)
library(mvtnorm)

############################
# Step 1: Read and Prepare IBRD Data
############################

# Read sheet 2 from the Excel file.
# This sheet appears to contain default indicators over time.
IBRD <- read.xlsx("IBRD.xlsx", sheet = 2)

# Initialize vectors:
#   DT = number of defaults per period
#   DF = default rate per period
DF <- c()
DT <- c()

# These vectors are later used to prepend periods with zero defaults.
N1 <- rep(0, 23)
N2 <- rep(22, 23)

# For each time column (excluding the first column),
# compute total defaults and default rate.
for (i in 1:(ncol(IBRD) - 1)) {
  DT[i] <- sum(IBRD[, (i + 1)], na.rm = TRUE)
  DF[i] <- mean(IBRD[, (i + 1)], na.rm = TRUE)
}

# Recover number of observations from defaults / default rate
N <- DT / DF

# Add 23 periods with zero defaults
DF <- c(N1, DF)
DT <- c(N1, DT)
N  <- c(N2, N)

# Read sheet 3, which seems to contain the cleaned / final time series:
#   defaults = number of defaults
#   total    = total obligor count
IBRD <- read.xlsx("IBRD.xlsx", sheet = 3)
DT <- IBRD$defaults
N  <- IBRD$total
DF <- DT / N

############################
# Step 2: Hand-Written Method-of-Moments (MM) Estimation
############################

# Equation used to recover rho from the bivariate normal model.
# It is solved numerically by uniroot().
binorm <- function(rho, p) {
  p <- E_Z
  sigma <- matrix(c(1, rho, rho, 1), nrow = 2)
  return(pmvnorm(lower = -Inf, upper = c(a, a), mean = c(0, 0), sigma) - p)
}

# First derivative of the bivariate normal CDF with respect to rho.
phi_prime <- function(s, t_val, rho) {
  y <- (1 / ((2 * pi) * sqrt(1 - rho^2))) *
    exp(-(s^2 / 2 - rho * s * t_val + t_val^2 / 2) / (1 - rho^2))
  return(y)
}

# Second derivative of the bivariate normal CDF with respect to rho.
phi_prime_prime <- function(s, t_val, rho) {
  y <- ((s * t_val + rho * (1 - s^2 - t_val^2) + s * t_val * rho^2 - rho^3) /
          (2 * pi * (1 - rho^2)^(5 / 2))) *
    exp(-(s^2 / 2 - rho * s * t_val + t_val^2 / 2) / (1 - rho^2))
  return(y)
}

# Compute several adjusted MM estimators of rho.
# The adjustment uses variance and autocovariances up to lag 5.
# Note: 60 is hard-coded because the time dimension here is 60.
adjusted_rho <- function(rho, g_prime, g_prime_prime,
                         alpha_0, alpha_1, alpha_2,
                         alpha_3, alpha_4, alpha_5) {
  
  adj_rho <- c()
  
  adj_rho[1] <- rho + (g_prime_prime / (2 * 60 * (g_prime)^3)) *
    alpha_0 + (g_prime_prime / (60 * (g_prime)^3) * (1 - 1 / 60) * alpha_1)
  
  adj_rho[2] <- rho + g_prime_prime / (60 * (g_prime)^3) *
    (alpha_0 / 2 + (1 - 1 / 60) * alpha_1 + (1 - 2 / 60) * alpha_2)
  
  adj_rho[3] <- rho + g_prime_prime / (60 * (g_prime)^3) *
    (alpha_0 / 2 + (1 - 1 / 60) * alpha_1 + (1 - 2 / 60) * alpha_2 +
       (1 - 3 / 60) * alpha_3)
  
  adj_rho[4] <- rho + g_prime_prime / (60 * (g_prime)^3) *
    (alpha_0 / 2 + (1 - 1 / 60) * alpha_1 + (1 - 2 / 60) * alpha_2 +
       (1 - 3 / 60) * alpha_3 + (1 - 4 / 60) * alpha_4)
  
  adj_rho[5] <- rho + g_prime_prime / (60 * (g_prime)^3) *
    (alpha_0 / 2 + (1 - 1 / 60) * alpha_1 + (1 - 2 / 60) * alpha_2 +
       (1 - 3 / 60) * alpha_3 + (1 - 4 / 60) * alpha_4 +
       (1 - 5 / 60) * alpha_5)
  
  return(adj_rho)
}

# Construct the MM inputs:
#   Z_t = DF_t^2
#   E_Z = average of Z_t
#   a   = normal-threshold corresponding to average default rate
Z <- DF * DF
E_Z <- mean(DF * DF)
a <- qnorm(mean(DF))

# Variance and autocovariances of Z_t
var_1 <- var(Z)
cov_1 <- sum((Z[1:59] - E_Z) * (Z[2:60] - E_Z)) / 60
cov_2 <- sum((Z[1:58] - E_Z) * (Z[3:60] - E_Z)) / 60
cov_3 <- sum((Z[1:57] - E_Z) * (Z[4:60] - E_Z)) / 60
cov_4 <- sum((Z[1:56] - E_Z) * (Z[5:60] - E_Z)) / 60
cov_5 <- sum((Z[1:55] - E_Z) * (Z[6:60] - E_Z)) / 60

# Recover rho from the bivariate normal equation
hat_rho <- uniroot(binorm, c(-1, 1))$root

# Compute adjusted rho estimates
ad_rho <- adjusted_rho(
  hat_rho,
  phi_prime(a, a, hat_rho),
  phi_prime_prime(a, a, hat_rho),
  var_1, cov_1, cov_2, cov_3, cov_4, cov_5
)

# Since rho = w^2 in the one-factor setup, report sqrt(rho)
ad_MM <- sqrt(ad_rho)
or_MM <- sqrt(hat_rho)

# Store MM results:
# first entry = original MM estimate
# next five entries = adjusted MM estimates
MM <- c(or_MM, ad_MM)
MM

############################
# Step 3: Hand-Written Maximum Likelihood Estimation (MLE)
############################

# Negative log-likelihood for the one-factor Gaussian default model.
# Parameters:
#   hat_w[1] = asset correlation loading w
#   hat_w[2] = default threshold gamma
logll <- function(hat_w) {
  F <- c()
  w <- hat_w[1]
  gamma <- hat_w[2]
  
  for (i in 1:length(DT)) {
    f <- function(u) {
      p <- pnorm((gamma - w * u) / sqrt(1 - w^2))
      obj <- dbinom(DT[i], N[i], p) * dnorm(u)
    }
    F[i] <- integrate(f, lower = -5, upper = 10)$value
  }
  
  return(sum(log(F)) * (-1))
}

# Estimate (w, gamma) by bounded optimization
MLE <- optim(
  c(0.01, -3.5),
  logll,
  method = "L-BFGS-B",
  lower = c(0.01, -3.5),
  upper = c(0.99, -1)
)$par

MLE

############################
# Step 4: Closed-Form MLE Approximation
############################

# Transform observed default rates into probit scale
sigma_a <- qnorm(DF)

# Replace -Inf values caused by zero default rates.
# The replacement value is user-chosen and acts as a practical fix.
for (i in 1:length(sigma_a)) {
  if (sigma_a[i] == -Inf) {
    sigma_a[i] <- -1.232
  }
}

# Estimate variance of transformed default rates
Var_sigma_a <- mean(sigma_a * sigma_a) - (mean(sigma_a))^2

# Closed-form estimate of w
c_MLE <- sqrt(Var_sigma_a / (1 + Var_sigma_a))
c_MLE

# Package-based original MLE for comparison
sqrt(intraMLE(DT, N)$Original)

############################
# Step 5: AssetCorr Package Estimates
############################

# Conditional method-of-moments estimators with different lag choices
MM_Asset <- c(
  sqrt(intraCMM(DT, N, l = 0)$Original),
  sqrt(intraCMM(DT, N, l = 1)$Original),
  sqrt(intraCMM(DT, N, l = 2)$Original),
  sqrt(intraCMM(DT, N, l = 3)$Original),
  sqrt(intraCMM(DT, N, l = 4)$Original),
  sqrt(intraCMM(DT, N, l = 5)$Original)
)
MM_Asset

# Further package-based estimators
sqrt(intraFMM(DT, N)$Original)
sqrt(intraAMLE(DT, N, Adjust = 0.001)$Original)
sqrt(intraAMM(DT, N)$Original)
sqrt(intraJDP2(DT, N)$Original)
sqrt(intraBeta(DT, N)$Original)
sqrt(intraMode(DT, N)$Original)

############################
# Step 6: Basel-Style Correlation Calibration
############################

# This function backs out an implied asset correlation from:
#   - mean loss (mu)
#   - standard deviation of loss (sigma)
#   - loss given default (LGD)
#
# The idea is:
#   1. Infer PD from mean loss and LGD
#   2. Fit a beta distribution to portfolio loss
#   3. Take the 99.9% quantile
#   4. Solve for the correlation consistent with the Basel formula
correlation_basel <- function(mu, sigma, LGD) {
  
  # Infer default probability from average loss
  PD <- mu / LGD
  
  # Parameters of the beta distribution
  alfa <- mu * (mu * (1 - mu) / sigma^2 - 1)
  beta <- (1 - mu) * (mu * (1 - mu) / sigma^2 - 1)
  
  # 99.9% quantile of total loss
  x <- qbeta(0.999, alfa, beta)
  
  # Unexpected loss = tail loss - mean loss
  UL <- x - mu
  
  # Quantities used in the Basel capital formula
  a <- qnorm(PD)
  b <- qnorm(0.999)
  c <- x / LGD
  d <- qnorm(c)
  
  # Solve for correlation x
  f <- function(x) {
    return(d - (sqrt(1 / (1 - x)) * (a + sqrt(x) * b)))
  }
  
  r <- uniroot(f, c(0, 0.9999), tol = 0.001)
  return(r$root)
}

# Compute implied correlations under different LGD assumptions
Lr <- 0.7 * DF
C7 <- sqrt(correlation_basel(mean(Lr), sd(Lr), 0.7))

Lr <- 0.8 * DF
C8 <- sqrt(correlation_basel(mean(Lr), sd(Lr), 0.8))

Lr <- 0.4 * DF
C4 <- sqrt(correlation_basel(mean(Lr), sd(Lr), 0.4))

Lr <- 0.6 * DF
C6 <- sqrt(correlation_basel(mean(Lr), sd(Lr), 0.6))

Lr <- 0.45 * DF
C45 <- sqrt(correlation_basel(mean(Lr), sd(Lr), 0.45))

Lr <- 0.1 * DF
C1 <- sqrt(correlation_basel(mean(Lr), sd(Lr), 0.1))

############################
# Step 7: Basel-Style Unexpected Loss / Capital Requirement
############################

# Basel-style unexpected loss formula:
#   UL = LGD * Phi((Phi^{-1}(PD) + sqrt(rho) * Phi^{-1}(0.999)) / sqrt(1-rho))
#        - PD * LGD
UL_basel <- function(df, LGD, rho) {
  PD <- mean(df)
  o <- sqrt(1 / (1 - rho))
  u <- qnorm(PD) + sqrt(rho) * qnorm(0.999)
  UL <- LGD * pnorm(o * u) - PD * LGD
  return(UL)
}

############################
# Capital Calculations Using MLE Estimate
############################
UL_basel(DF, 0.6, MLE[1]^2)
UL_basel(DF, 0.7, MLE[1]^2)
UL_basel(DF, 0.8, MLE[1]^2)
UL_basel(DF, 0.45, MLE[1]^2)
UL_basel(DF, 0.1, MLE[1]^2)

############################
# Capital Calculations Using Original MM Estimate
############################
UL_basel(DF, 0.6, MM[1]^2)
UL_basel(DF, 0.7, MM[1]^2)
UL_basel(DF, 0.8, MM[1]^2)
UL_basel(DF, 0.45, MM[1]^2)
UL_basel(DF, 0.1, MM[1]^2)

############################
# Capital Calculations Using Basel-Implied Correlations
############################
UL_basel(DF, 0.6, C6^2)
UL_basel(DF, 0.7, C7^2)
UL_basel(DF, 0.8, C8^2)