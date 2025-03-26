# Helper functions

gev_fun <- function(model_fit, pc_or_q, type = c("qgev", "pgev"), se_factor = 0) {
    type <- match.arg(type)

    gev_loc <- model_fit$estimate["loc"]
    gev_scale <- model_fit$estimate["scale"]
    gev_shape <- model_fit$estimate["shape"]
    gev_loc_se <- model_fit$std.err["loc"]
    gev_scale_se <- model_fit$std.err["scale"]
    gev_shape_se <- model_fit$ std.err["shape"]


    loc <- gev_loc + se_factor * gev_loc_se
    scale <- gev_scale + se_factor * gev_scale_se
    shape <- gev_shape + se_factor * gev_shape_se

    do.call(type, list(pc_or_q, loc, scale, shape))
}

# Simulate data
sim_pois_truncnorm <- function(k, n_lambda, mu, sigma, min_size) {
    n <- rpois(k, lambda = n_lambda) # number of fish per competition
    maxima <- numeric(k)
    for (i in seq_len(k)) {
        maxima[i] <- max(
            truncnorm::rtruncnorm(
                n = n[i],
                a = min_size,
                mean = mu,
                sd = sigma
            )
        )
    }
    return(tibble(maxima, n))
}

find_max_x <- function(mu, sigma, n) {
    optimize(function(x) calc_max(x, n, mu, sigma),
        interval = c(0, mu+sigma*5),
        maximum = TRUE
    )$maximum
}

# f_{max}(x) 
# pdf of the maximum values
f_n_x <- function(x, n, mu, sigma) {
  fx <- dnorm(x, mu, sigma) 
  Fx <- pnorm(x, mu, sigma) 
  (n * (Fx^(n - 1))) * fx
}

f_max_x_log <- function(x, n, mu, sigma) {
  fx_log <- dnorm(x, mu, sigma, log = TRUE)
  Fx <- pnorm(x, mu, sigma)
  
  # Avoid log(0) issues
  Fx[Fx == 0] <- .Machine$double.eps
  
  log_val <- log(n) + (n - 1) * log(Fx) + fx_log
  return(exp(log_val))
}


find_max_mode <- function(n, mu, sigma) {
optimize(
  f = function(x) f_n_x(x, n = n, mu = mu, sigma = sigma), 
  interval = c(0, mu+(5*sigma)), 
  maximum = TRUE
)$maximum
}

f_max_x_log <- function(x, n, mu, sigma) {
  fx_log <- dnorm(x, mu, sigma, log = TRUE)
  Fx <- pnorm(x, mu, sigma)
  
  # Avoid log(0) issues
  Fx[Fx < .Machine$double.eps] <- .Machine$double.eps
  
  log_val <- log(n) + (n - 1) * log(Fx) + fx_log
  return(exp(log_val))
}

# quicker approach to estimate max given parameters and sample size
# expected_max <- function(n, mu, sigma) {
#   integrand <- function(x) x * f_n_x(x, n, mu, sigma)
#   integrate(integrand, lower = -Inf, upper = mu + 5 * sigma)$value
# }

# quicker approach to estimate max given parameters and sample size
# this uses the log function to avoid very large exponents 
expected_max <- function(n, mu, sigma) {
  integrand <- function(x) {
    f_max_x_log(x, n, mu, sigma) * x
  }
  integrate(integrand, lower = mu - 4 * sigma, upper = mu + 6 * sigma, rel.tol = 1e-6)$value
}


# an slower sim approach to test the theory
simulate_expected_max <- function(sample_size, mu = 50, sigma = 16, nreps = 1e5){
    replicate(nreps, max(rnorm(sample_size, mean = mu, sd = sigma))) |> mean()
}


# F_{max}(x)^n
# cdf of the maximum values
F_n_x <- function(x, n, mu, sigma) {
  Fx <- pnorm(x, mu, sigma) 
  Fx^n
}


# an slower sim approach to test the theory
simulate_expected_max_gev <- function(n, loc, scale, shape, nreps = 1e5) {
   replicate(nreps, max(rgev(n, loc = loc, scale = scale, shape = shape))) |> mean()
}

f_max_x_log_gev <- function(x, k, loc, scale, shape) {
  fx_log <- dgev(x = x, loc = loc, scale = scale, shape = shape, log = TRUE)
  Fx <- pgev(q = x, loc = loc, scale = scale, shape = shape)
  
  # Avoid log(0) issues
  Fx[Fx < .Machine$double.eps] <- .Machine$double.eps
  
  log_val <- log(k) + (k - 1) * log(Fx) + fx_log
  return(exp(log_val))
}

expected_max_gev <- function(k, loc, scale, shape) {
  integrand <- function(x) {
    f_max_x_log_gev(x = x, k = k, loc = loc, scale = scale, shape = shape) * x
  }
  integrate(integrand, lower = 0, upper = 200, rel.tol = 1e-6)$value
}
