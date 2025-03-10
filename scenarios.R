library(tidyverse)

# example simulation
# sim_pois_truncnorm(seed = 1, k = 10, n_lambda = 50, mu = 66, sigma = 22, min_size = 10)

# simulation over all scenarios
sim_tbl <- 
  expand_grid(
    seed = 1:10, 
    k = c(3, 5, 10, 20, 30, 100), 
    n_lambda = c(30, 50, 200),
    mu = c(20, 50, 100),
  ) |> 
  mutate(sigma = mu*0.34,
         min_size = 10) |> 
  mutate(maxima = purrr::pmap(
    .l = list(seed, k, n_lambda, mu, sigma, min_size), 
    .f = sim_pois_truncnorm))




