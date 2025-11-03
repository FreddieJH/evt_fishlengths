# Following the example from Ullman et al. (2022)
library(dplyr)
library(purrr)

source("R/01_funcs.R")

toadfish_maxima <- c(78.5, 80, 83, 120)

fit_mod <- function(maxima, model_type) {
  maxima_list <- as.list(maxima)

  mod_dat <-
    list(
      x = maxima_list |> unlist(),
      n_obs = length(maxima_list |> unlist()),
      n_per_sample = lapply(maxima_list, length) |> unlist(),
      start_idx = cumsum(lapply(maxima_list, length) |> unlist()),
      k = length(maxima_list)
    )

  init_func <- function(type, maxima_median) {
    if (type == "evt") {
      return_func <- function(chain_id) {
        list(
          loc = maxima_median,
          scale = 10,
          shape = 0
        )
      }
    } else if (type == "evtg") {
      return_func <- function(chain_id) {
        list(
          loc = maxima_median,
          scale = 10
        )
      }
    } else {
      return_func <- function(chain_id) {
        list(
          mu = maxima_median,
          sigma = 10,
          lambda = 100
        )
      }
    }
  }

  mod <- switch(
    model_type,
    "evt" = cmdstanr::cmdstan_model("models/evt.stan"),
    "evtg" = cmdstanr::cmdstan_model("models/evt_gumbel.stan"),
    "efs" = cmdstanr::cmdstan_model("models/efs.stan"),
    "efsmult" = cmdstanr::cmdstan_model("models/efs.stan")
  )

  fit <- mod$sample(
    data = mod_dat,
    chains = 4,
    init = init_func(model_type, median(maxima)),
    iter_warmup = 9000,
    iter_sampling = 1000
  )
  return(fit)
}

get_posterior <- function(fit) {
  posterior <-
    fit |>
    posterior::as_draws_df() %>%
    dplyr::as_tibble()
  return(posterior)
}

est_max <- function(fit, ci = 0.8, k) {
  is_evt <- sum(c("loc", "scale", "shape") %in% colnames(get_posterior(fit))) ==
    3
  is_evtg <- sum(c("loc", "scale") %in% colnames(get_posterior(fit))) == 2
  get_posterior(fit) |>
    mutate(
      pdf = if (is_evt) {
        pmap_dbl(
          .l = list(
            p = 1 - (1 / k),
            loc = loc,
            scale = scale,
            shape = shape
          ),
          .f = qgev
        )
      } else if (is_evtg) {
        pmap_dbl(
          .l = list(
            p = 1 - (1 / k),
            loc = loc,
            scale = scale
          ),
          .f = qgumbel
        )
      } else {
        pmap_dbl(
          .l = list(mu, sigma, lambda, k),
          .f = \(mu, sigma, lambda, k) {
            cdf <- \(x) {
              ptnorm(q = x, mean = mu, sd = sigma)
            }
            pdf <- \(x) {
              dtnorm(x = x, mean = mu, sd = sigma)
            }
            gmax <- \(x) {
              g(x = x, n = lambda * k, cdf = cdf, pdf = pdf)
            }
            mode_f(gmax)
          }
        )
      }
    ) |>
    summarise(
      max_fit = mean(pdf),
      max_lwr = quantile(pdf, (1 - ci) / 2),
      max_upr = quantile(pdf, 1 - ((1 - ci) / 2)),
    ) %>%
    as.numeric()
}


est_max_posterior <- function(fit, ci = 0.8, k) {
  is_evt <- sum(c("loc", "scale", "shape") %in% colnames(get_posterior(fit))) ==
    3
  is_evtg <- sum(c("loc", "scale") %in% colnames(get_posterior(fit))) == 2
  get_posterior(fit) |>
    mutate(
      estmax = if (is_evt) {
        pmap_dbl(
          .l = list(
            p = 1 - (1 / k),
            loc = loc,
            scale = scale,
            shape = shape
          ),
          .f = qgev
        )
      } else if (is_evtg) {
        pmap_dbl(
          .l = list(
            p = 1 - (1 / k),
            loc = loc,
            scale = scale
          ),
          .f = qgumbel
        )
      } else {
        pmap_dbl(
          .l = list(mu, sigma, lambda, k),
          .f = \(mu, sigma, lambda, k) {
            cdf <- \(x) {
              ptnorm(q = x, mean = mu, sd = sigma)
            }
            pdf <- \(x) {
              dtnorm(x = x, mean = mu, sd = sigma)
            }
            gmax <- \(x) {
              g(x = x, n = lambda * k, cdf = cdf, pdf = pdf)
            }
            mode_f(gmax)
          }
        )
      }
    )
}

fit_evt <- fit_mod(toadfish_maxima, "evt")
fit_evtg <- fit_mod(toadfish_maxima, "evtg")
fit_efs <- fit_mod(toadfish_maxima, "efs")

est_max(fit_evt, k = 4) |> round()
# est_max(fit_evt, k = 20)
est_max(fit_evtg, k = 4) |> round()
# est_max(fit_evtg, k = 20)
est_max(fit_efs, k = 4) |> round()
# est_max(fit_efs, k = 20)

ecdf(est_max_posterior(fit_evt, k = 4)$estmax)(120)
ecdf(est_max_posterior(fit_evtg, k = 4)$estmax)(120)
ecdf(est_max_posterior(fit_efs, k = 4)$estmax)(120)
