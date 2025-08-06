scenario_data_top5 <-
  read_csv("data/simulation/scenarios_top5.csv", show_col_types = FALSE) |>
  summarise(top5 = list(top5), .by = -top5) |>
  select(scenario_id, k, lambda, top5) |>
  distinct() |>
  mutate(maxima = map_dbl(top5, max)) |>
  mutate(rank_maxima_byscen = n() - rank(maxima) + 1, .by = scenario_id)

max_median5 <-
  scenario_data_top5 |>
  filter(rank_maxima_byscen <= 5) |>
  mutate(median_5maxima = median(maxima), .by = scenario_id)


max_max1 <-
  scenario_data_top5 |>
  filter(rank_maxima_byscen == 1) |>
  mutate(median_1maxima = median(maxima), .by = scenario_id)


scenarios_truemax <-
  scenarios |>
  mutate(
    true_max = pmap_dbl(
      .l = list(
        distr = dist_name,
        n = k * lambda,
        mean = dist_mean,
        variance = (dist_mean * 0.34)^2
      ),
      .f = expected_max_fromsim
    ),
    true_max20 = pmap_dbl(
      .l = list(
        distr = dist_name,
        n = 20 * lambda,
        mean = dist_mean,
        variance = (dist_mean * 0.34)^2
      ),
      .f = expected_max_fromsim
    )
  ) |>
  select(scenario_id, true_max, true_max20)


scenarios_truemax |>
  left_join(max_median5) |>
  left_join(max_max1)


scenarios_repeated <-
  scenarios |>
  expand_grid(rep = 1:1000)


get_pars <- function(distr, mean) {
  sd <- 0.34 * mean
  variance <- sd^2
  if (distr == "gamma") {
    par1 <- (mean^2) / variance # shape
    par2 <- mean / variance # rate (= 1/scale)
  } else if (distr == "tnorm") {
    par1 <- mean # mu
    par2 <- sqrt(variance) # sigma
  } else if (distr == "lnorm") {
    par1 <- log(mean) - log(1 + variance / (mean^2)) / 2 # meanlog
    par2 <- sqrt(log(1 + variance / (mean^2)))
  }
  return(c(par1, par2))
}

sim_dist <- function(distr, pars, n) {
  rdist <- get(paste0("r", distr))
  x <- rdist(n, pars[1], pars[2])
  return(x)
}

get_top5 <- function(x) sort(x) |> tail(5)

extract_top <- function(distr, mean, n) {
  sim_dist(distr, get_pars(distr, mean), n) |>
    max()
}

sim_out <-
  scenarios_repeated |>
  uncount(k) |>
  mutate(
    maxima = pmap_dbl(
      .l = list(distr = dist_name, mean = dist_mean, n = lambda),
      .f = extract_top
    )
  ) |>
  nest(.by = c(rep, scenario_id)) |>
  mutate(
    top5_median = map_dbl(data, function(x) median(get_top5(x$maxima))),
    max = map_dbl(data, function(x) max(x$maxima))
  ) |>
  select(-data)

scenarios_repeated |>
  uncount(k) |>
  head()
mutate(pars = map2(dist_name, dist_mean, get_pars))

xx2 <-
  scenarios_repeated |>
  uncount(k) |>
  head() |>
  mutate(pars = map2(dist_name, dist_mean, get_pars))


sim_dist("tnorm", xx2[1, ]$pars, 100) |> get_top5()
