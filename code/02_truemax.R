source("code/00_funcs.R")
source("code/01_simulation.R")
library(purrr)
library(truncnorm)


scenarios_truemax <-
  underlying_scenarios |>
  mutate(
    pars = pmap(
      .l = list(
        distr = dist_name,
        mean = dist_mean,
        variance = (dist_mean * 0.34)^2
      ),
      .f = get_dist_pars
    )
  ) |>
  unnest_wider(pars, names_sep = "") |>
  mutate(
    true_max = pmap_dbl(
      .l = list(
        distr = dist_name,
        n = k * lambda,
        par1 = pars1,
        par2 = pars2
      ),
      .f = expected_max
    ),
    true_max20 = pmap_dbl(
      .l = list(
        distr = dist_name,
        n = 20 * lambda,
        par1 = pars1,
        par2 = pars2
      ),
      .f = expected_max
    )
  )
