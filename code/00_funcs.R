dtnorm <- function(x, par1, par2) {
  truncnorm::dtruncnorm(x = x, a = 0, mean = par1, sd = par2)
}
ptnorm <- function(x, par1, par2) {
  truncnorm::ptruncnorm(q = x, a = 0, mean = par1, sd = par2)
}
qtnorm <- function(x, par1, par2) {
  truncnorm::qtruncnorm(p = x, a = 0, mean = par1, sd = par2)
}
rtnorm <- function(x, par1, par2) {
  truncnorm::rtruncnorm(n = x, a = 0, mean = par1, sd = par2)
}

# Function: given distribution mean and variance, what are the parameter values
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

# Function: PDF of maxima given the PDF and CDF of x
g_max <- function(x, distr, n, par1, par2) {
  f_x <- function(x) get(paste0("d", distr))(x, par1, par2)
  F_x <- function(x) get(paste0("p", distr))(x, par1, par2)

  # using log to avoid problems with very small pdf and cdf values
  # g_max = function(x) n * (F_x(x)^(n - 1)) * f_x(x)
  log_g_max <- log(n) + (n - 1) * log(F_x(x)) + log(f_x(x))
  return(exp(log_g_max))
}

# expected value is just integral of x*g(x) where g(x) is the pdf of the max values
expected_max <- function(distr, n, par1, par2) {
  integrand <- function(x) {
    x * g_max(x, distr = distr, n = n, par1 = par1, par2 = par2)
  }
  upper_bound <- get(paste0("q", distr))(0.9999999999999, par1, par2)
  integrate(integrand, lower = 0, upper = upper_bound, rel.tol = 1e-6)$value
}
