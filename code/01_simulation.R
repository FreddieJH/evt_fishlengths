source("code/00a_truncnorm_functions.R")
library(dplyr)
library(tidyr)
library(tibble)
library(purrr)
library(readr)

set.seed(123)

# Scenerios used in the analysis:
# 1. k = number of samples (i.e., this could a fishing competiton)
# 2. lambda = poisson parameter used to define sample size (n) of each k (100, 1000, 10000)
# 3. mean = mean of the underlying distribution (10, 50, 100)
# 4. dist = distribution used to generate the maxima (truncated normal, lognormal, gamma)

# Other:
# 5. variance = variance of the underlying distribution (calculated from mean and cv of 0.34)
# 6. min_size = minimum body length of the fish (postive values only)

# Vectorised sampling function
sample_fish_data <- function(dist_name, dist_mean, n) {
  dist_var <- (dist_mean * 0.34)^2
  dist_pars <-
    switch(
      dist_name,
      "gamma" = list(
        shape = dist_mean^2 / dist_var,
        rate = dist_mean / dist_var
      ),
      "tnorm" = list(mean = dist_mean, sd = sqrt(dist_var)),
      "lnorm" = {
        meanlog <- log(dist_mean) - 0.5 * log(1 + dist_var / dist_mean^2)
        sdlog <- sqrt(log(1 + dist_var / dist_mean^2))
        list(meanlog = meanlog, sdlog = sdlog)
      }
    )

  switch(
    dist_name,
    "gamma" = rgamma(n, shape = dist_pars$shape, rate = dist_pars$rate),
    "tnorm" = rtnorm(
      n,
      a = 0,
      mean = dist_pars$mean,
      sd = dist_pars$sd
    ),
    "lnorm" = rlnorm(n, meanlog = dist_pars$meanlog, sdlog = dist_pars$sdlog)
  )
}

# Underlying distribution scenarios (n=9)
# f(x, \theta)
# f = c(lnorm, tnorm, gamma)
# \theta = (\mu, \sigma^2)
# \mu = c(10, 50, 100) cm
# \sigma^2 = (0.34\mu^2)
# underlyingdist_scenarios <-
#   expand_grid(dist_name = c("lnorm", "tnorm", "gamma"),
# dist_mean = c(10, 50, 100))

# # How we sample from the underlying dist (n=15)
# sampling_scenarios <-
#   expand_grid(
#     k = c(5, 10, 50, 100, 200),
#     lambda = c(100, 1000, 10000)
#   )

# # 135 input scenarios
# input_scenarios <-
#   expand_grid(underlyingdist_scenarios,
#   sampling_scenarios)

# Create scenarios
scenarios <- expand_grid(
  k = c(5, 10, 50, 100, 200),
  lambda = c(100, 1000, 10000),
  dist_name = c("tnorm", "lnorm", "gamma"),
  dist_mean = c(10, 50, 100)
) |>
  mutate(scenario_id = row_number())

model_scenarios <-
  tibble(
    gamma_mean = c(8000, 2000, 8000, 2000),
    gamma_variance = c(4000000, 500000, 1000000, 4000000)
  ) |>
  mutate(
    gamma_shape = (gamma_mean^2) / gamma_variance,
    gamma_rate = gamma_mean / gamma_variance
  ) |>
  expand_grid(
    type = c("efs", "efsmult"),
    mu_prior_mean = c(30, 50, 100, 150)
  ) |>
  mutate(mu_prior_sd = mu_prior_mean / 1.5) |>
  bind_rows(tibble(type = "evt")) |>
  mutate(
    scenario_model = ifelse(
      type == "evt",
      "evt",
      paste0(type, "_shape", gamma_shape, "_mu", mu_prior_mean)
    )
  )

scenarios_filename <- "data/simulation/scenarios_multimaxima.csv"

if (file.exists(scenarios_filename)) {
  scenario_data <-
    read_csv(scenarios_filename, show_col_types = FALSE) |>
    summarise(top_m = list(top_m), .by = -top_m)
} else {
  simulated_data <-
    scenarios |>
    mutate(n = map2(k, lambda, ~ rpois(.x, .y))) %>%
    unnest(n) %>%
    mutate(j = row_number(), .by = scenario_id) %>%
    mutate(x = pmap(list(dist_name, dist_mean, n), sample_fish_data))

  scenario_data <-
    simulated_data |>
    expand_grid(model_scenarios) |>
    mutate(filename = paste0(scenario_model, "_", scenario_id)) |>
    mutate(top5 = map(.x = x, ~ sort(.x, decreasing = TRUE)[1:5])) |>
    select(-x) |>
    mutate(
      m = case_when(
        type == "efsmult" ~ sample(1:4, n(), replace = TRUE),
        .default = 1
      )
    ) |>
    mutate(
      top_m = map2(top5, m, .f = \(x, y) {
        head(sort(x, decreasing = TRUE), n = y)
      })
    )

  scenario_data |>
    select(-top5) |>
    unnest(top_m) |>
    write_csv(file = scenarios_filename)

  scenario_data |>
    select(-top_m) |>
    unnest(top5) |>
    write_csv(file = "data/simulation/scenarios_top5.csv")
}
