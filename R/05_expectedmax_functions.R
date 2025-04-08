# X_n represent the maximum value of a sample size of n
# for that sample of size n, to be all greater than Pr(X_n <= x) = F(x)^n
# the PDF of X_n is the derivation of f(X_n) = n[F(x)^n]^(n-1)f(x)
# the expected value of X_n is the integration of x*f(X_n)
# therefore E(X_n) = integrate(x*(n[F(x)^n]^(n-1)f(x))) wrt x

cdf_max <- function(x, distr, n, par1, par2) {
  cdf <- function(x) get(paste0("p", distr))(x, par1, par2)
  cdf(x)^n
}

pdf_max <- function(x, distr, n, par1, par2) {
  pdf <- function(x) get(paste0("d", distr))(x, par1, par2)
  cdf <- function(x) get(paste0("p", distr))(x, par1, par2)
  # avoiding very small pdf and cdf values
  log_pdf_max <- log(n) + (n - 1) * log(cdf(x)) + log(pdf(x))
  exp(log_pdf_max)
}

expected_max <- function(distr, n, par1, par2) {
  integrand <- function(x) x * pdf_max(x, distr, n, par1, par2)
  upper_bound <- get(paste0("q", distr))(0.9999999999999, par1, par2)
  integrate(integrand, lower = 0, upper = upper_bound, rel.tol = 1e-6)$value
}

expected_max_fromsim <- function(distr, n, mean, variance) {
  dist_pars <- calc_dist_pars(dist = distr, mean = mean, var = variance)
  est_max <- expected_max(dist = distr, n = n, par1 = dist_pars[[1]], par2 = dist_pars[[2]])
  return(est_max)
}

expected_max_evt <- function(loc, scale, shape, k) {
  if (shape == 0) {
    expected_max <- loc + scale * log(k)
  } else {
    expected_max <- loc + (scale / shape) * ((k^shape) - 1)
  }
}



expected_max_brute <- function(distr, par1, par2, n, iterations = 1000) {
  sample_func <- function(n, par1, par2) get(paste0("r", distr))(n, par1, par2)
  mean(replicate(iterations, max(sample_func(n = n, par1, par2))))
}

