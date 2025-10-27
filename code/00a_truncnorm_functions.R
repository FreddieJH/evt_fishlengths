# Truncated normal distribution functions
# The default lower bound is set at zero for this analysis (postive body lengths only)
# these functions saves loading in the 'truncnorm' package

rtnorm <- function(n, mean, sd, a = 0, b = Inf) {
  if (a >= b) stop("Lower bound must be less than upper bound")
  
  # Handle special cases
  if (a == -Inf && b == Inf) return(rnorm(n, mean, sd))
  if (sd <= 0) stop("Standard deviation must be positive")
  
  # Standardise bounds
  alpha <- (a - mean) / sd
  beta <- (b - mean) / sd
  
  # Calculate standardised CDF values at bounds
  p_alpha <- pnorm(alpha)
  p_beta <- pnorm(beta)
  
  # Generate uniform random variables in the range (p_alpha, p_beta)
  u <- runif(n, p_alpha, p_beta)
  
  # Transform to truncated normal via inverse CDF
  x <- qnorm(u)
  
  # Rescale to original scale
  return(mean + sd * x)
}

dtnorm <- function(x, mean, sd, a = 0, b = Inf) {
  if (a >= b) stop("Lower bound must be less than upper bound")
  if (sd <= 0) stop("Standard deviation must be positive")
  
  # Return 0 for x outside bounds
  result <- rep(0, length(x))
  valid <- (x >= a) & (x <= b)
  
  if (any(valid)) {
    # For valid x, calculate normalized density
    norm_const <- pnorm((b - mean) / sd) - pnorm((a - mean) / sd)
    result[valid] <- dnorm(x[valid], mean, sd) / norm_const
  }
  
  return(result)
}

ptnorm <- function(q, mean, sd, a = 0, b = Inf) {
  if (a >= b) stop("Lower bound must be less than upper bound")
  if (sd <= 0) stop("Standard deviation must be positive")
  
  # Handle vector inputs
  result <- rep(0, length(q))
  
  # Set to 0 for q < a, 1 for q > b
  result[q < a] <- 0
  result[q > b] <- 1
  
  # For a <= q <= b, compute normalized CDF
  valid <- (q >= a) & (q <= b)
  if (any(valid)) {
    norm_const <- pnorm((b - mean) / sd) - pnorm((a - mean) / sd)
    result[valid] <- (pnorm((q[valid] - mean) / sd) - pnorm((a - mean) / sd)) / norm_const
  }
  
  return(result)
}

qtnorm <- function(p, mean, sd, a = 0, b = Inf) {
  if (a >= b) stop("Lower bound must be less than upper bound")
  if (sd <= 0) stop("Standard deviation must be positive")
  if (any(p < 0 | p > 1)) stop("Probabilities must be between 0 and 1")
  
  # Handle standard normal case
  if (a == -Inf && b == Inf) return(qnorm(p, mean, sd))
  
  # Standardise bounds
  alpha <- (a - mean) / sd
  beta <- (b - mean) / sd
  
  # Calculate standardised CDF values at bounds
  p_alpha <- pnorm(alpha)
  p_beta <- pnorm(beta)
  
  # Transform p to the standard normal scale
  p_std <- p_alpha + p * (p_beta - p_alpha)
  
  # Apply inverse normal CDF and rescale
  return(mean + sd * qnorm(p_std))
}
