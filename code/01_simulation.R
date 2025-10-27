library(tidyr)
library(dplyr)
library(readr)

# Scenerios used in the analysis:
# 1. k = number of samples (i.e., this could a fisher in a competition, 5, 10, 50, 100, and 200 fishers)
# 2. lambda = poisson parameter used to define sample size (n) of each k (100, 1000, 10000)
# 3. dist = distribution used to generate the maxima (truncated normal, lognormal, gamma)
# 4. mean = mean of the underlying distribution (10, 50, 100)
# 5. variance = variance of the underlying distribution (calculated from mean and cv of 0.34)

# Model predictions
# 6. Three models are used to estimate LMAX (EVT, EFS and EFSMULT)

# Other:
# 7. For the two EFS methods, the priors on the parameters (lambda, mean, variance) we varied

# 135 underlying distribution scenarios
underlying_scenarios <- expand_grid(
  k = c(5, 10, 50, 100, 200),
  lambda = c(100, 1000, 10000),
  dist_name = c("tnorm", "lnorm", "gamma"),
  dist_mean = c(10, 50, 100)
) |>
  mutate(dist_scen_id = row_number())

# 33 scenarios for the priors on the model
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
  bind_rows(tibble(type = "evt"))

# A total of 4455 scenario combinations
scenarios <-
  expand_grid(underlying_scenarios, model_scenarios) |>
  mutate(scenario_id = row_number())


if (!file.exists("data/simulated/scenarios_maxima.csv")) {
  sample_topm <- function(distr, dist_mean, n, m) {
    pars <- get_dist_pars(distr, dist_mean, (dist_mean * 0.34)^2)
    get(paste0("r", distr))(n, pars[1], pars[2]) |>
      sort(decreasing = TRUE) |>
      head(m)
  }

  set.seed(1)
  simulated_data <-
    underlying_scenarios |>
    mutate(n = map2(k, lambda, ~ rpois(.x, .y))) %>%
    unnest(n) %>%
    mutate(j = row_number(), .by = dist_scen_id) %>%
    mutate(x = pmap(list(dist_name, dist_mean, n, m = 5), sample_topm))

  scenarios_maxima <-
    simulated_data |>
    expand_grid(model_scenarios) |>
    mutate(m = ifelse(type == "efsmult", sample(1:4, replace = TRUE), 1)) %>%
    mutate(
      top_m = map2(.x = x, .y = m, ~ sort(.x, decreasing = TRUE)[1:.y])
    ) |>
    select(-x) |>
    unnest(top_m)

  write_csv(scenarios_maxima, file = "data/simulated/scenarios_maxima.csv")
} else {
  scenarios_maxima <- read_csv("data/simulated/scenarios_maxima.csv")
}
