source("R/01_funcs.R")

library(readr)
library(tidyr)
library(dplyr)
library(purrr)

if (!exists("scenarios")) {
  if (!file.exists("results/data/scenarios.csv")) {
    source("R/02_simulation.R")
  }
  scenarios <-
    read_csv("results/data/scenarios.csv", show_col_types = FALSE) |>
    nest(data = topm) |>
    mutate(topm = map(data, \(x) as.numeric(x$topm))) |>
    select(-data) |>
    nest(
      .by = c(k, lambda, dist_name, dist_mean, scenario_id),
      .key = "samples"
    )
}

scenarios_truemax <-
  scenarios |>
  mutate(
    true_max = pmap_dbl(
      .l = list(k, lambda, dist_name, dist_mean),
      .f = \(k, lambda, dist_name, dist_mean) {
        pars <- get_dist_pars(distr = dist_name, mean = dist_mean)
        pdf <- \(x) get(paste0("d", dist_name))(x, pars[1], pars[2])
        cdf <- \(x) get(paste0("p", dist_name))(x, pars[1], pars[2])
        mode <- mode_f(\(x) g(x, n = k * lambda, cdf = cdf, pdf = pdf))
        return(mode)
      }
    ),
    true_max20 = pmap_dbl(
      .l = list(lambda, dist_name, dist_mean),
      .f = \(lambda, dist_name, dist_mean) {
        pars <- get_dist_pars(distr = dist_name, mean = dist_mean)
        pdf <- \(x) get(paste0("d", dist_name))(x, pars[1], pars[2])
        cdf <- \(x) get(paste0("p", dist_name))(x, pars[1], pars[2])
        mode <- mode_f(\(x) g(x, n = 20 * lambda, cdf = cdf, pdf = pdf))
        return(mode)
      }
    )
  )

write_csv(scenarios_truemax, "results/data/scenarios_truemax.csv")
