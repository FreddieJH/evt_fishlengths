source("code/02_simulation.R")

get_start_id <- function(vec) lag(c(0, cumsum(vec)))[-1] + 1

get_stan_data <- function(k, data, multiple = FALSE) {
  if (multiple) {
    list(
      k = k,
      x = unlist(data$topm),
      n_obs = sum(data$m),
      n_per_sample = data$m,
      start_idx = get_start_id(data$m)
    )
  } else {
    list(
      k = k,
      x = unlist(data$top1),
      n_obs = k,
      n_per_sample = rep(1, k),
      start_idx = get_start_id(rep(1, k))
    )
  }
}

scenarios_stan <-
  scenarios %>%
  mutate(maxima_median = map_dbl(samples, ~ median(.x$top1))) |>
  mutate(stan_list_single = map2(.x = k, .y = samples, .f = get_stan_data)) %>%
  mutate(
    stan_list_multpl = map2(
      .x = k,
      .y = samples,
      .f = ~ get_stan_data(.x, .y, multiple = TRUE)
    )
  )
