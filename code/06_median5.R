set.seed(123)
distr = "tnorm"
mean = 50 #cm
sd = 50 * 0.34
total_n = 1e6

k = 10 # number of samples
n_k = rep(1000, k) # sample size of each sample

rtnorm <- function(n, mean, sd, a = 0, b = Inf) {
  if (a >= b) {
    stop("Lower bound must be less than upper bound")
  }
  if (a == -Inf && b == Inf) {
    return(rnorm(n, mean, sd))
  }
  if (sd <= 0) {
    stop("Standard deviation must be positive")
  }
  alpha <- (a - mean) / sd
  beta <- (b - mean) / sd
  p_alpha <- pnorm(alpha)
  p_beta <- pnorm(beta)
  u <- runif(n, p_alpha, p_beta)
  x <- qnorm(u)
  return(mean + sd * x)
}

distr_r <- get(paste0("r", distr))
x = distr_r(n = total_n, mean = mean, sd = sd)

split(sample(x, size = sum(n_k)), rep(1:k, n_k)) |>
  lapply(max) |>
  unlist() |>
  sort(decreasing = TRUE) |>
  head(5) |>
  median()


sim_top5 <- function(distr = "tnorm", mean, sd, total_n, k, n_k_fixed = 1000) {
  n_k = rep(n_k_fixed, k) # sample size of each sample

  distr_r <- get(paste0("r", distr))
  x = distr_r(n = total_n, mean = mean, sd = sd)

  split(sample(x, size = sum(n_k)), rep(1:k, n_k)) |>
    lapply(max) |>
    unlist() |>
    sort(decreasing = TRUE) |>
    head(5) |>
    median()
}

sims =
  expand.grid(
    distr = "tnorm",
    mean = 80,
    sd = 80 * 0.34,
    total_n = seq(1e6, 5e6, by = 1e6),
    k = seq(10, 50, by = 5),
    rep = 1:200
  )

# sims$res <- purrr::pmap_dbl(sims, sim_top5)

if (!file.exists("data/simulation/median5_mean80.csv")) {
  # 40 mins to run 200 repeats
  t1 <- Sys.time()
  sims$res <- purrr::pmap_dbl(
    sims[, c("distr", "mean", "sd", "total_n", "k")],
    sim_top5
  )
  Sys.time() - t1
  write.csv(x = sims, "data/simulation/median5_mean80.csv")
}

library(dplyr)
library(ggplot2)
library(readr)
read_csv("data/simulation/median5_mean80.csv") |>
  bind_rows(read_csv("data/simulation/median5.csv")) |>
  dplyr::summarise(
    median = median(res),
    sd = sd(res),
    n = n(),
    .by = c(distr, mean, sd, total_n, k)
  ) |>
  mutate(lwr = median - (1.96 * sd), upr = median + (1.96 * sd)) |>
  ggplot(aes(total_n, median)) +
  geom_point() +
  geom_errorbar(aes(ymin = lwr, ymax = upr)) +
  facet_grid(mean ~ k, scales = "free_y")

ggsave("results/figures/median5_mean50_mean70_mean80.png")
