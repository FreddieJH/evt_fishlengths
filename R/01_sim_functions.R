
calc_dist_pars <- function(dist = c("gamma", "tnorm", "lnorm"), mean, var){
  # check for the dist argument
  if(length(dist) > 1) stop("Please provide one of the following for the 'dist' argument: 'gamma', 'tnorm', 'lnorm'")
  dist <- match.arg(dist)

  if (dist == "gamma") {
    shape <- (mean^2) / var
    rate <- mean / var # rate is 1/scale
    return(list(shape, rate))
    } else if (dist == "tnorm") {
    mu <- mean
    sigma <- sqrt(var)
    return(list(mu, sigma))
    } else if (dist == "lnorm") {
    meanlog <- log(mean) - log(1 + var/(mean^2))/2
    sdlog <- sqrt(log(1 + var/(mean^2)))
    return(list(meanlog, sdlog))
    }
}

# rtnorm func with default limit vals
rtnorm <- function(x, mean, sd, a = 0, b = Inf) {
  truncnorm::rtruncnorm(n = x, mean = mean, sd = sd, a = a, b = b)
}
ptnorm <- function(x, mean, sd, a = 0, b = Inf) {
  truncnorm::ptruncnorm(q = x, mean = mean, sd = sd, a = a, b = b)
}
qtnorm <- function(x, mean, sd, a = 0, b = Inf) {
  truncnorm::qtruncnorm(p = x, mean = mean, sd = sd, a = a, b = b)
}
dtnorm <- function(x, mean, sd, a = 0, b = Inf) {
  truncnorm::dtruncnorm(x = x, mean = mean, sd = sd, a = a, b = b)
}

simulate_max_anydist <- function(k, lambda, mean, variance, dist = c("gamma", "tnorm", "lnorm")) {
  
  if(length(dist) > 1) stop("Please provide one of the following for the 'dist' argument: 'gamma', 'tnorm', 'lnorm'")
  dist <- match.arg(dist)
  
  # sample size per sample
  # each sample has a given sample size that is derived by the poisson distribution with a mean of lambda
  n_k <- rpois(n = k, lambda = lambda)
  
  dist_pars <- calc_dist_pars(dist = dist, mean = mean, var = variance)

  # the order of the parameters is important
  # e.g. the rgamma function takes the shape argument before the rate (= 1/scale) argument 
  # therefore we give it the rate argument for par2 instead of the scale argument
  get_max <- function(x) max(get(paste0("r", dist))(x, dist_pars[[1]], dist_pars[[2]]))

  maxvals <- lapply(X = n_k, FUN = get_max)

  return(list(n = n_k, max = maxvals))
}
simulate_multimax_anydist <- function(k, lambda, mean, variance, dist = c("gamma", "tnorm", "lnorm")) {
  # Select one distribution if multiple options are provided
  dist <- match.arg(dist)
  
  # Generate random sample sizes
  n_k <- rpois(n = k, lambda = lambda)
  
  # Calculate distribution parameters
  dist_pars <- calc_dist_pars(dist = dist, mean = mean, var = variance)
  
  # Generate random number of maxima to extract from each sample
  m <- runif(n = k, min = 1, max = 5) %>% round() %>% as.integer()
  
  # Make sure m doesn't exceed n_k for any sample
  m <- pmin(m, n_k)
  
  # Function to get m_i maxima from a sample of size n_i
  get_m_max <- function(i) {
    sample_size <- n_k[i]
    num_maxima <- m[i]
    sample_values <- get(paste0("r", dist))(sample_size, dist_pars[[1]], dist_pars[[2]])
    sort(sample_values, decreasing = TRUE)[1:num_maxima]
  }
  
  # Generate maxima for each sample
  maxvals_list <- lapply(1:k, get_m_max)
  
  # Create output tibble with the requested columns
  output_tibble <- tibble::tibble(
    sample_id = 1:k,
    n_k = n_k,
    m_k = m,
    maxvals = maxvals_list
  )
  
  # Add simulation parameters as attributes
  attr(output_tibble, "simulation_params") <- list(
    k = k,
    lambda = lambda,
    mean = mean,
    variance = variance,
    dist = dist,
    dist_pars = dist_pars
  )
  
  return(output_tibble)
}


prepare_stan_data <- function(sim_results) {
  # Count unique samples
  k <- length(unique(sim_results$sample_id))
  
  # Count total fish observations
  total_fish <- nrow(sim_results)
  
  # Count number of maxima in each sample
  m <- sim_results %>%
    dplyr::group_by(sample_id) %>%
    dplyr::summarise(count = n()) %>%
    dplyr::pull(count)
  
  # Extract sample_id for each observation
  # Note: Stan uses 1-indexed arrays, so we ensure sample_ids start from 1
  sample_id <- sim_results$sample_id
  
  # Extract fish sizes
  x <- sim_results$value
  
  # Create list in Stan format
  stan_data <- list(
    k = k,
    total_fish = total_fish,
    m = m,
    sample_id = sample_id,
    x = x
  )
  
  return(stan_data)
}


