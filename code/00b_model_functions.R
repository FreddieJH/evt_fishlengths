library(cmdstanr)
library(posterior)
library(arrow)

# This is when we have multiple maxima per sample
# for this we need to data in slightly different format
#
stan_code_efsmult <- function(
  gamma_shape,
  gamma_rate,
  mu_prior_mean,
  mu_prior_sd
) {
  paste0(
    "
    data {
      int<lower=1> k;                    // Number of competitions/samples
      int<lower=1> n_obs;                // Total number of observations
      array[n_obs] real x;                   // All recorded fish sizes
      array[k] int n_per_sample; // Number of fish per sample
      array[k] int start_idx;   // Starting index for each sample
    }
    
    parameters {
      real<lower=0> mu;                  // Mean of the fish size distribution
      real<lower=0.001> sigma;           // Standard deviation
      real<lower=1> lambda;              // Estimated mean number of fish caught per competition
    }
    
    model {
      // Priors
      mu ~ normal(",
    mu_prior_mean,
    ", ",
    mu_prior_sd,
    ") T[0.001,];
       sigma ~ normal(",
    mu_prior_mean / 3,
    ", ",
    mu_prior_sd / 3,
    ") T[0.001,];
      lambda ~ gamma(",
    gamma_shape,
    ",",
    gamma_rate,
    ");
      
      target += -lambda * k;
      
      {
        
      // area of normal CDF above zero (normalisation constant)
      real norm_const = 1 - normal_cdf(0 | mu, sigma);  
        
      // looping over each sample
        for(j in 1:k) {
          
        target += n_per_sample[j] * log(lambda);

          // we use the smallest value since we know all other values are smaller than x_m
          // we could use the start_indx if we know the x values are in increasing order
          // but this is the safest way
          real smallest_in_m = min(x[start_idx[j]:(start_idx[j] + n_per_sample[j] - 1)]);

          // Truncated normal CDF
          real trunc_cdf = (normal_cdf(smallest_in_m | mu, sigma) - normal_cdf(0 | mu, sigma)) / norm_const;

          target += lambda * trunc_cdf;
          
          // looping over each maxima within the sample j
          for(i in 1:n_per_sample[j]) {
            int obs_idx = start_idx[j] + i - 1;
            real trunc_lpdf = normal_lpdf(x[obs_idx] | mu, sigma) - log(norm_const);
            target += trunc_lpdf;
          }
        }
      }
    }"
  )
}

stan_code_efs <- function(gamma_shape, gamma_rate, mu_prior_mean, mu_prior_sd) {
  paste0(
    "

  data {
    int<lower=1> k;  // Number of competitions
    vector[k] x;  // Recorded largest fish sizes at each competition
  }
  
  parameters {
    real<lower=0> mu;    // Mean of the fish size distribution
    real<lower=0.001> sigma; // Standard deviation
    real<lower=1> lambda; // Estimated mean number of fish caught per competition
  }
  
  model {
    // Priors
    mu ~ normal(",
    mu_prior_mean,
    ", ",
    mu_prior_sd,
    ") T[0.001,];
    sigma ~ normal(",
    mu_prior_mean / 3,
    ", ",
    mu_prior_sd / 3,
    ") T[0.001,];
    lambda ~ gamma(",
    gamma_shape,
    ",",
    gamma_rate,
    ");
    
    
    target += k*log(lambda) - lambda*k;
    
    for (i in 1:k) {

    // Using truncated normal instead of normal
    real norm_const = 1 - normal_cdf(0 | mu, sigma);  // Normalisation constant
    real trunc_cdf = (normal_cdf(x[i] | mu, sigma) - normal_cdf(0 | mu, sigma)) / norm_const;
    real trunc_lpdf = normal_lpdf(x[i] | mu, sigma) - log(norm_const);
    target += lambda * trunc_cdf + trunc_lpdf;

    }


  }"
  )
}

stan_code_evt <- {
  "data {
  int<lower=0> k;           // number of observations
  vector[k] x;              // observed maxima values
}

parameters {
  real loc;                  // location parameter
  real<lower=0> scale;      // scale parameter (must be positive)
  real shape;                  // shape parameter
}

model {
  // Priors
  loc ~ normal(mean(x), 10);     // centered around sample mean with wide variance
  scale ~ lognormal(0, 1);      // ensures positivity, reasonably diffuse
  shape ~ normal(0, 0.5);          // shape parameter typically small, centered at 0
  
  // GEV likelihood
  for (i in 1:k) {
    if (shape == 0) {
      target += -(log(scale) + (x[i] - loc)/scale + exp(-(x[i] - loc)/scale));
    } else {
      real z = (x[i] - loc)/scale;
      if (1 + shape * z > 0) {
        target += -(log(scale) + (1 + 1/shape) * log1p(shape * z) + pow(1 + shape * z, -1/shape));
      } else {
        target += negative_infinity();  // outside support
      }
    }
  }
}"
}

fit_maxima_model <- function(
  maxima,
  type,
  n_obs = NULL,
  n_per_sample = NULL,
  start_id = NULL,
  iter_warmup = 2000,
  iter_sampling = 1000,
  chains = 4,
  refresh = 500,
  init = NULL,
  gamma_shape = NULL,
  gamma_rate = NULL,
  mu_prior_mean = NULL,
  mu_prior_sd = NULL,
  k = NULL
) {
  #   init is required for the evt to avoid wild estimates
  if (is.null(init)) {
    init <- if (type == "evt") {
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

  stan_file <- switch(
    type,
    "evt" = cmdstanr::write_stan_file(stan_code_evt, dir = "models"),
    "efs" = cmdstanr::write_stan_file(
      stan_code_efs(gamma_shape, gamma_rate, mu_prior_mean, mu_prior_sd),
      dir = "models"
    ),
    "efsmult" = cmdstanr::write_stan_file(
      stan_code_efsmult(gamma_shape, gamma_rate, mu_prior_mean, mu_prior_sd),
      dir = "models"
    )
  )

  k <- switch(
    type,
    "evt" = length(maxima),
    "efs" = length(maxima),
    "efsmult" = k
  )

  mod <- cmdstanr::cmdstan_model(stan_file, stanc_options = list("O1"))

  fit <- mod$sample(
    data = list(
      k = k,
      x = maxima,
      n_obs = n_obs,
      n_per_sample = n_per_sample,
      start_idx = start_id
    ),
    iter_warmup = iter_warmup,
    iter_sampling = iter_sampling,
    chains = chains,
    refresh = refresh,
    init = init
  )

  return(fit)
}

get_posterior <- function(fit) {
  fit %>%
    posterior::as_draws_df() %>%
    dplyr::as_tibble() %>%
    dplyr::select(-lp__)
}

get_summary <- function(fit) {
  fit$summary()
}

# concept_evt_fit %>%
#     posterior::as_draws_df() %>%
#     dplyr::as_tibble() %>%
#     dplyr::select(-lp__)

# multiple_efs_fit <- function(maxima, gamma_shape, gamma_rate) {

#   mod_results <- fit_maxima_model(maxima, type = "efs", gamma_shape = gamma_shape, gamma_rate = gamma_rate)

#   mod_results$draws %>%
#     mutate(
#       max = pmap_dbl(
#         .l = list(distr = "tnorm", n = length(maxima)*lambda, par1 = mu, par2 = sigma, p = 1-(1/length(maxima))),
#         .f = inverse_G_x
#       ),
#       max20 = pmap_dbl(
#         .l = list(distr = "tnorm", n = 20*lambda, par1 = mu, par2 = sigma, p = 1-(1/length(maxima))),
#         .f = inverse_G_x
#       )
#     ) %>%
#     rename_with(~ paste0("est_", .x), setdiff(names(.), c(".chain", ".iteration", ".draw")))
# }

# multiple_evt_fit <- function(maxima) {

#   mod_results <- fit_maxima_model(maxima, type = "evt")

#   mod_results$draws %>%
#     mutate(
#       max = pmap_dbl(
#         .l = list(k = length(maxima), loc = mu, scale = sigma, shape = xi),
#         .f = expected_max_evt
#       ),
#       max20 = pmap_dbl(
#         .l = list(k = 20, loc = mu, scale = sigma, shape = xi),
#         .f = expected_max_evt
#       )
#     ) %>%
#     rename_with(~ paste0("est_", .x), setdiff(names(.), c(".chain", ".iteration", ".draw")))
# }

# multiple_mod_fits <- function(type, gamma_shape = NULL, gamma_rate = NULL){
#   simulation_output_filepath <-
#     if(type == "efs"){
#       paste0("data/model_output/efs_posterior_scenarios_shape",gamma_shape ,".parquet")
#     } else {
#       paste0("data/model_output/evt_posterior_scenarios.parquet")
#     }

#   if (!file.exists(simulation_output_filepath)) {

#     plan(multisession)
#     max_posterior_multiple <-
#       scenarios_maxima |>
#       nest(maxima_tbl = c(max, n)) %>%
#       mutate(maxima = map(maxima_tbl, ~ .x$max)) %>%
#       mutate(out = future_pmap(
#         .l = list(maxima = maxima),
#         .f = if(type == "efs") {
#           function(maxima) multiple_efs_fit(maxima, gamma_shape = gamma_shape, gamma_rate = gamma_rate)
#         } else {
#           multiple_evt_fit
#         },
#         .options = furrr_options(seed = TRUE,
#                                  globals = c("ptnorm", "multiple_efs_fit", "fit_maxima_model", "write_stan_file", "stan_code_efs", gamma_rate = gamma_rate,
#                                                            gamma_shape = gamma_shape, "as_tibble", "as_draws_df", "rename_with", "inverse_G_x", "G_max", "F_x"),
#                                  packages = c("cmdstanr"))
#       )) |>
#       unnest(cols = out) %>%
#       select(filename, k, contains("est_"))
#     plan(sequential)

#     arrow::write_parquet(max_posterior_multiple, simulation_output_filepath)
#   } else {
#     max_posterior_multiple <- arrow::read_parquet(simulation_output_filepath)
#   }
#   return(max_posterior_multiple)
# }

# summarise_posterior <- function(posterior_tbl){
#   posterior_tbl %>%
#     summarise(
#       max = quantile(est_max, 0.5),
#       max_lwr = quantile(est_max, 0.025),
#       max_upr = quantile(est_max, 0.975),
#       max20 = quantile(est_max20, 0.5),
#       max20_lwr = quantile(est_max20, 0.025),
#       max20_upr = quantile(est_max20, 0.975),
#       .by = filename
#     )
# }
