# About functions --------------------------------------------------------------------------------

# f_x: PDF of a distribution (e.g., lognormal, truncated normal etc)
# F_x: CDF of a distribution (e.g., lognormal, truncated normal etc)
# g_max: Calculates the PDF of the maxima
# G_max: Calculates the CDF of the maxima
# inverse_G_x: Finds the quantile (inverse CDF) of the maxima distribution at probability p (e.g. value at 95% percentile)
# expected_max: Calculates the expected value of the maxima distribution (using integration)
# expected_max_fromsim: Uses expected_max function, but first calculates the distribution parameters based on mean and variance (useful for simulation)
# expected_max_evt: Calculates expected maximum based on a given quantile
# expected_max_brute: Estimates expected maximum using Monte Carlo simulation for validation purposes

# Script requirements  --------------------------------------------------------------

library(evd)

# Functions --------------------------------------------------------------

# PDF(x) = f(x) = Pr(X = x)
f_x <- function(x, distr, par1, par2) {
  get(paste0("d", distr))(x, par1, par2)
}

# CDF(x) = F(x) = Pr(X <= x)
F_x <- function(x, distr, par1, par2) {
  get(paste0("p", distr))(x, par1, par2)
}

# the PDF of the maxima = g(x_max) = Pr(X_n = x) = n * (F(x)^n-1) * f(x)
g_max <- function(x, distr, n, par1, par2) {
  f_x_clean <- function(x) f_x(x, distr = distr, par1 = par1, par2 = par2)
  F_x_clean <- function(x) F_x(x, distr = distr, par1 = par1, par2 = par2)

  # g_max = function(x) n * (F_x(x)^(n - 1)) * f_x(x)

  # using log to avoid problems with very small pdf and cdf values
  log_f_max <- log(n) + (n - 1) * log(F_x_clean(x)) + log(f_x_clean(x))
  return(exp(log_f_max))
}

# the CDF of the maxima = G(x_max) = F(x)^n
G_max <- function(x, distr, n, par1, par2) {
  F_x_clean <- function(x) F_x(x, distr = distr, par1 = par1, par2 = par2)
  return(F_x_clean(x)^n)
}

# used to calculate the value at a given perentile of the G_max()
inverse_G_x <- function(
  distr,
  n,
  par1,
  par2,
  p,
  interval_lwr = 1,
  interval_upr = 1000
) {
  uniroot(
    function(x) G_max(x, distr = distr, n = n, par1 = par1, par2 = par2) - p,
    lower = interval_lwr,
    upper = interval_upr
  )$root
}

# expected value is just integral of x*g(x) where g(x) is the pdf of the max values
expected_max <- function(distr, n, par1, par2) {
  integrand <- function(x) {
    x * g_max(x, distr = distr, n = n, par1 = par1, par2 = par2)
  }
  upper_bound <- get(paste0("q", distr))(0.9999999999999, par1, par2)
  integrate(integrand, lower = 0, upper = upper_bound, rel.tol = 1e-6)$value
}


get_dist_pars <- function(distr, mean, variance) {
  if (distr == "gamma") {
    par1 <- (mean^2) / variance # shape
    par2 <- mean / variance # rate (= 1/scale)
  } else if (distr == "tnorm") {
    par1 <- mean # mu
    par2 <- sqrt(variance) # sigma
  } else if (distr == "lnorm") {
    par1 <- log(mean) - log(1 + variance / (mean^2)) / 2 # meanlog
    par2 <- sqrt(log(1 + variance / (mean^2)))
  }
  return(c(par1, par2))
}

expected_max_fromsim <- function(distr, n, mean, variance) {
  pars <- get_dist_pars(distr = distr, mean = mean, variance = variance)

  est_max <- expected_max(dist = distr, n = n, par1 = pars[1], par2 = pars[2])
  return(est_max)
}

# is this function needed??
expected_max_evt <- function(loc, scale, shape, k) {
  evd::qgev(1 - (1 / k), loc = loc, scale = scale, shape = shape)
}

# this is 'brute force' approach to estimate expected maximum for any of the distributions
# mainly just for testing the expected_max function
expected_max_brute <- function(distr, par1, par2, n, iterations = 1000) {
  sample_func <- function(n, par1, par2) get(paste0("r", distr))(n, par1, par2)
  mean(replicate(iterations, max(sample_func(n = n, par1, par2))))
}
