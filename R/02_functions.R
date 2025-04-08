# get a set of maxima values from a given population
# k = number of samples (therefore number of maxima)
# n = sample size per sample
# x is a vector of values distributed according to some underlying population
get_sample_maxima <- function(x, k, n) {
  maxima <- replicate(k, {
    sample <- sample(x, n)
    max(sample)
  })
}

# fit GEV distribution to maxima values and return a sample of the fitted distribution
# maxima = a set of maxima values (of length k), derived from the get_sample_maxima() function
sample_gev <- function(maxima, n_out = 10000) {
    gev_fit <- evd::fgev(maxima)
    rgev(n_out,
        loc = gev_fit$estimate["loc"], scale = gev_fit$estimate["scale"], shape = gev_fit$estimate["shape"]
    )
}

# generate maxima based on a truncated normal distributon
# where the sample size is definite by a poisson distribution (lambda)
# k = number of samples
# lambda = poisson parameter used to define sample size of each k
# mu, sigma = location and variance parameters of the truncated normal
get_max_truncnorm <- function(k, lambda, mu, sigma, min_size = 0) {
    n_k <- rpois(n = k, lambda = lambda) # sample size per k
    set.seed(1)
    maxvals <- 
    map(
      .x = n_k,
      .f = ~ rtruncnorm(
        n = .x,
        a = min_size,
        mean = mu,
        sd = sigma
      )
    ) %>% 
    map(max) %>% unlist()

  return(list(n = n_k, max = maxvals))
}



# F_{max}(x)^n
# cdf of the maximum values
F_n_x <- function(x, n, mu, sigma) {
  Fx <- pnorm(x, mu, sigma) 
  Fx^n
}

evt_expected_max <- function(loc, scale, shape, k){
if (shape == 0) {
    expected_max <- loc + scale * log(k)
} else {
    expected_max <- loc + (scale / shape) * ((k^shape) - 1)
}
}




multiple_efs_fit <- function(maxima, model_path = "models/efs.stan") {

  mod_results <- fit_maxima_model(maxima, model_path)

  mod_results$draws %>%
    mutate(
        max = pmap_dbl(
            .l = list(n = length(maxima)*lambda, mu = mu, sigma = sigma),
            .f = efs_expected_max
        ),
        max100 = pmap_dbl(
            .l = list(n = 100*lambda, mu = mu, sigma = sigma),
            .f = efs_expected_max
        )
    ) %>% 
    rename_with(~ paste0("est_", .x), setdiff(names(.), c(".chain", ".iteration", ".draw"))) 
}

multiple_mod_fits <- function(type){
simulation_output_filepath <- paste0("data/model_output/", type,"_posterior_scenarios_multdist.parquet")
if (!file.exists(simulation_output_filepath)) {

  plan(multisession)
  max_posterior_multiple <-
    scenarios_maxima |>
    select(-k_i) %>% 
    nest(maxima_tbl = c(maxvals, sim_n)) %>% 
    mutate(maxima = map(maxima_tbl, ~ .x$maxvals)) %>% 
    mutate(out = future_pmap(
      .l = list(maxima = maxima),
      .f = get(paste0("multiple_",type,"_fit")), 
      .options = furrr_options(seed = TRUE)
    )) |>
    unnest(cols = out) %>%
    select(filename, k, contains("est_"))  
  plan(sequential)
  
  arrow::write_parquet(max_posterior_multiple, simulation_output_filepath)
} else {
  max_posterior_multiple <- arrow::read_parquet(simulation_output_filepath)
}
return(max_posterior_multiple)
}

summarise_posterior <- function(posterior_tbl){
    posterior_tbl %>% 
    summarise(
        max = quantile(est_max, 0.5),
        max_lwr = quantile(est_max, 0.025),
        max_upr = quantile(est_max, 0.975),
        max100 = quantile(est_max100, 0.5),
        max100_lwr = quantile(est_max100, 0.025),
        max100_upr = quantile(est_max100, 0.975), 
        .by = filename
    )
}
