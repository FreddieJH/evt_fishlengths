

fit_maxima_model <- function(maxima,type,
                             iter_warmup = 3000, 
                             iter_sampling = 4000,
                             chains = 4,
                             parallel_chains = 4,
                             refresh = 500,
                             init = NULL, 
                             gamma_shape = 0,  
                             gamma_rate = 0) {
  
  
  # model_type <- if(grepl("evt", model_path)) "evt" else "efs"
  
  #   init is required for the evt to avoid wild estimates
  if (is.null(init)) {
    init <- if(type == "evt") {
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
          lambda = 1000
        )
      }
    }
  }
  
  
  stan_code_efs <- paste0("
  data {
    int<lower=1> k;  // Number of competitions
    vector[k] x;  // Recorded largest fish sizes
  }
  
  parameters {
    real<lower=0> mu;    // Mean of the fish size distribution
    real<lower=0.001> sigma; // Standard deviation
    real<lower=1> lambda; // Estimated mean number of fish caught per competition
  }
  
  model {
    // Priors
    mu ~ normal(30, 20) T[0.001,];
    sigma ~ normal(10, 5) T[0.001,];
    lambda ~ gamma(", gamma_shape,",", gamma_rate,");
    
    
    target += k*log(lambda) - lambda*k;
    
    for (i in 1:k) {
    // Using truncated normal instead of normal
    real norm_const = 1 - normal_cdf(0 | mu, sigma);  // Normalisation constant
    real trunc_cdf = (normal_cdf(x[i] | mu, sigma) - normal_cdf(0 | mu, sigma)) / norm_const;
    real trunc_lpdf = normal_lpdf(x[i] | mu, sigma) - log(norm_const);
    target += lambda * trunc_cdf + trunc_lpdf;
    }
  }")
  
  if(type == "evt"){
    mod <- cmdstan_model("models/evt.stan", stanc_options = list("O1"))
  } else {
    stan_file <- write_stan_file(stan_code_efs)
    mod <- cmdstan_model(stan_file, stanc_options = list("O1"))
  }
  
  
  
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
  
  draws <- 
    fit %>%
    as_draws_df() %>%
    as_tibble() %>%
    dplyr::select(-lp__)
  
  list(
    fit = fit,
    draws = draws
  )
}

multiple_efs_fit <- function(maxima, gamma_shape, gamma_rate) {
  
  mod_results <- fit_maxima_model(maxima, type = "efs", gamma_shape = gamma_shape, gamma_rate = gamma_rate)
  
  mod_results$draws %>%
    mutate(
      max = pmap_dbl(
        .l = list(distr = "tnorm", n = length(maxima)*lambda, par1 = mu, par2 = sigma),
        .f = expected_max
      ),
      max100 = pmap_dbl(
        .l = list(distr = "tnorm", n = 100*lambda, par1 = mu, par2 = sigma),
        .f = expected_max
      )
    ) %>% 
    rename_with(~ paste0("est_", .x), setdiff(names(.), c(".chain", ".iteration", ".draw"))) 
}


multiple_evt_fit <- function(maxima) {
  
  mod_results <- fit_maxima_model(maxima, type = "evt")
  
  mod_results$draws %>%
    mutate(
      max = pmap_dbl(
        .l = list(k = length(maxima), loc = mu, scale = sigma, shape = xi),
        .f = expected_max_evt
      ),
      max100 = pmap_dbl(
        .l = list(k = 100, loc = mu, scale = sigma, shape = xi),
        .f = expected_max_evt
      )
    ) %>% 
    rename_with(~ paste0("est_", .x), setdiff(names(.), c(".chain", ".iteration", ".draw"))) 
}

multiple_mod_fits <- function(type, gamma_shape = 0, gamma_rate = 0){
  simulation_output_filepath <- paste0("data/model_output/", type,"_posterior_scenarios_shape",gamma_shape ,".parquet")
  if (!file.exists(simulation_output_filepath)) {
    
    plan(multisession)
    max_posterior_multiple <-
      scenarios_maxima |>
      select(-k_i) %>% 
      nest(maxima_tbl = c(maxvals, sim_n)) %>% 
      mutate(maxima = map(maxima_tbl, ~ .x$maxvals)) %>% 
      mutate(out = future_pmap(
        .l = list(maxima = maxima),
        .f = if(type == "efs") {
          function(maxima) multiple_efs_fit(maxima, gamma_shape = gamma_shape, gamma_rate = gamma_rate)
        } else {
          multiple_evt_fit
        },
        .options = furrr_options(seed = TRUE,  globals = TRUE)
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
