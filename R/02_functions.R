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



simulate_maxima <- function(scenario_tbl) {
  simulation_maxima_filepath <- "data/simulation/scenario_maxima.csv"
  if (file.exists(simulation_maxima_filepath)) {
    # check if all the simulations have been done (maximia calculated)
    sim_maxima_tbl <- read_csv(simulation_maxima_filepath)
    completed_scenarios <- unique(sim_maxima_tbl$filename) %>% head(8900)
    
    scenarios_todo <-
      scenario_tbl %>%
      filter(!filename %in% completed_scenarios)
  } else {
    sim_maxima_tbl <- tibble()
    scenarios_todo <- scenario_tbl
  }
  
  if(nrow(scenarios_todo)) {
    
    new_rows <- 
      scenarios_todo %>%
      mutate(max_tbl = pmap(
        .l = list(
          k = k,
          lambda = n_lambda,
          mu = mu,
          sigma = sigma
        ),
        .f = get_max_truncnorm
      )) %>%
      mutate(
        sim_n = map(.x = max_tbl, .f = function(.x) .x$n),
        maxvals = map(.x = max_tbl, .f = function(.x) .x$max),
        k_i = map(k, seq_len)
      ) %>%
      select(-max_tbl) %>%
      unnest(cols = c(sim_n, maxvals, k_i))
    
    sim_maxima_tbl <- bind_rows(sim_maxima_tbl, new_rows)
    
    write_csv(sim_maxima_tbl, simulation_maxima_filepath)
    
  } 
  return(sim_maxima_tbl)
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
efs_expected_max <- function(n, mu, sigma) {
  integrand <- function(x) {
    f_max_x_log(x, n, mu, sigma) * x
  }
  integrate(integrand, lower = mu - 4 * sigma, upper = mu + 6 * sigma, rel.tol = 1e-6)$value
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


# Function to fit Stan models for maxima analysis
fit_maxima_model <- function(maxima, model_path, 
                           iter_warmup = 1000, 
                           iter_sampling = 2000,
                           chains = 4,
                           parallel_chains = 4,
                           refresh = 500,
                           init = NULL) {
  
  # Determine model type from path
  model_type <- if(grepl("evt", model_path)) "evt" else "efs"
  
  # Set default initialization based on model type
  if (is.null(init)) {
    init <- if(model_type == "evt") {
      function() {
        list(
          mu = median(maxima),
          sigma = sd(maxima),
          xi = 0
        )
      }
    } else {
      function() {
        list(
          mu = median(maxima),
          sigma = sd(maxima),
          lambda = 100
        )
      }
    }
  }
  
  # Compile model
  mod <- cmdstan_model(model_path, stanc_options = list("O1"))
  
  # Fit model
  fit <- 
    mod$sample(
    data = list(
      k = length(maxima),
      x = maxima
    ),
    iter_warmup = iter_warmup,
    iter_sampling = iter_sampling,
    chains = chains,
    parallel_chains = parallel_chains,
    refresh = refresh,
    init = init
  )
  
  # Process draws
  draws <- 
    fit %>%
    as_draws_df() %>%
    as_tibble() %>%
    dplyr::select(-lp__)

  # Return results
  list(
    fit = fit,
    draws = draws
  )
}
3

multiple_evt_fit <- function(maxima, model_path = "models/evt.stan") {

  mod_results <- fit_maxima_model(maxima, model_path)

  mod_results$draws %>%
    mutate(
        max = pmap_dbl(
            .l = list(k = length(maxima), loc = mu, scale = sigma, shape = xi),
            .f = evt_expected_max
        ),
        max100 = pmap_dbl(
            .l = list(k = 100, loc = mu, scale = sigma, shape = xi),
            .f = evt_expected_max
        )
    ) %>% 
    rename_with(~ paste0("est_", .x), setdiff(names(.), c(".chain", ".iteration", ".draw"))) 
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
simulation_output_filepath <- paste0("data/model_output/", type,"_posterior_scenarios.parquet")
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
