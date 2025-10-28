source("R/01_funcs.R")

library(readr)
library(tidyr)
library(dplyr)
library(purrr)


set.seed(1)

# Scenerios used in the analysis:
# 1. k = number of samples (i.e., this could a fisher in a competition, 5, 10, 50, 100, and 200 fishers)
# 2. lambda = poisson parameter used to define sample size (n) of each k (100, 1000, 10000)
# 3. dist = distribution used to generate the maxima (truncated normal, lognormal, gamma)
# 4. mean = mean of the underlying distribution (10, 50, 100) - note variance is based on mean using a fixed CV of 0.34

# 135 underlying distribution scenarios
scenarios_unnested <-
  expand_grid(
    k = c(5, 10, 50, 100, 200),
    lambda = c(100, 1000, 10000),
    dist_name = c("tnorm", "lnorm", "gamma"),
    dist_mean = c(10, 50, 100)
  ) |>
  mutate(n_k = map2(k, lambda, ~ rpois(.x, .y))) |>
  unnest(n_k) |>
  mutate(m = sample(size = n(), 1:4, replace = TRUE)) |>
  mutate(
    rsample = pmap(
      .l = list(dist_name, dist_mean, n_k),
      .f = ~ {
        pars <- get_dist_pars(..1, ..2)
        rsample <- get(paste0("r", ..1))(..3, pars[1], pars[2])
      }
    ),
    top1 = map_dbl(rsample, max),
    topm = map2(
      .x = rsample,
      .y = m,
      ~ sort(.x, decreasing = TRUE)[seq_len(.y)]
    )
  ) |>
  nest(.by = c(k, lambda, dist_name, dist_mean), .key = "samples") |>
  mutate(scenario_id = row_number()) |>
  unnest(samples) |>
  unnest(topm) |>
  select(-rsample)

write_csv(scenarios_unnested, "results/data/scenarios.csv")
