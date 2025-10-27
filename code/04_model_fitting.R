# Script requirements  --------------------------------------------------------------
source("code/01_funcs.R")
source("data-raw/R/calc_dist_pars.R")
source("code/04_model_prep.R")
init_func <- function(model_id, maxima_median) {
  if (model_id == "evt") {
    return_func <- function(chain_id) {
      list(
        loc = maxima_median,
        scale = 10,
        shape = 0
      )
    }
  } else if (model_id == "evtg") {
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

# Note: If your models are not running try a complete reinstall of "cmdstanr"
if (!"cmdstanr" %in% installed.packages()) {
  install.packages(
    "cmdstanr",
    repos = c('https://stan-dev.r-universe.dev', getOption("repos"))
  )
}

if (!file.exists("data/posterior.parquet")) {
  efs_mod <- cmdstanr::cmdstan_model("models/efs.stan")
  evt_mod <- cmdstanr::cmdstan_model("models/evt.stan")
  evt_gumbel_mod <- cmdstanr::cmdstan_model("models/evt_gumbel.stan")

  if (!file.exists("data/results_single.parquet")) {
    future::plan(multisession)
    results_single <-
      scenarios_stan |>
      head(3) |>
      select(scenario_id, stan_list_single) |>
      crossing(
        tibble(
          model_id = c("efs", "evt", "evtg"),
          model = list(efs_mod, evt_mod, evt_gumbel_mod)
        )
      ) %>%
      mutate(
        fit = future_pmap(
          .l = list(model, model_id, stan_list_single),
          .f = ~ ..1$sample(
            data = ..3,
            iter_warmup = 2000,
            iter_sampling = 1000,
            chains = 4,
            refresh = 1000,
            parallel_chains = 1,
            init = init_func(..2, median(..3$x)), # avoiding looking in complete wrong place
            adapt_delta = 0.999, # to avoid the small number of divergent transitions
            max_treedepth = 12 # increased from 10 to 12 after increasing the adapt_delta
          ) %>%
            posterior::as_draws_df() %>%
            tidyr::as_tibble() |>
            dplyr::select(-lp__) |>
            pivot_longer(-c(.chain, .iteration, .draw), names_to = "par"),
          .options = furrr_options(seed = TRUE)
        )
      ) %>%
      select(scenario_id, model_id, fit) %>%
      unnest(fit)
    future::plan(sequential)

    write_parquet(results_single, "data/results_single.parquet")
  } else {
    results_single <- read_parquet("data/results_single.parquet")
  }

  if (!file.exists("data/results_multpl.parquet")) {
    future::plan(multisession)

    results_multpl <-
      scenarios_stan |>
      head(3) |>
      select(scenario_id, stan_list_multpl) |>
      crossing(tibble(
        model_id = c("efsm"),
        model = list(efs_mod)
      )) %>%
      mutate(
        fit = future_pmap(
          .l = list(model, model_id, stan_list_multpl),
          .f = ~ ..1$sample(
            data = ..3,
            iter_warmup = 2000,
            iter_sampling = 1000,
            chains = 4,
            refresh = 1000,
            parallel_chains = 1,
            init = init_func(..2, median(..3$x)), # avoiding looking in complete wrong place
            adapt_delta = 0.999, # to avoid the small number of divergent transitions
            max_treedepth = 12 # increased from 10 to 12 after increasing the adapt_delta
          ) %>%
            posterior::as_draws_df() %>%
            tidyr::as_tibble() |>
            dplyr::select(-lp__) |>
            pivot_longer(-c(.chain, .iteration, .draw), names_to = "par"),
          .options = furrr_options(seed = TRUE)
        )
      ) %>%
      select(scenario_id, model_id, fit) %>%
      unnest(fit)

    future::plan(sequential)
    write_parquet(results_multpl, "data/results_multpl.parquet")
  } else {
    results_multpl <- read_parquet("data/results_multpl.parquet")
  }

  posterior <-
    bind_rows(
      read_parquet("data/results_single.parquet"),
      read_parquet("data/results_multpl.parquet")
    )
  write_parquet(posterior, "data/posterior.parquet")
} else {
  posterior <- read_parquet("data/posterior.parquet")
}

# fit_stan_mod <- function(stan_mod, data_list, init_f) {
#   fit <- stan_mod$sample(
#     data = data_list,
#     iter_warmup = 500,
#     iter_sampling = 1000,
#     chains = 4,
#     refresh = 100,
#     parallel_chains = 4,
#     init = init_f, # avoiding looking in complete wrong place
#     adapt_delta = 0.999, # to avoid the small number of divergent transitions
#     max_treedepth = 12 # increased from 10 to 12 after increasing the adapt_delta
#   )

#   posterior <-
#     fit %>%
#     as_draws_df() %>%
#     as_tibble()

#   summary <- fit$summary()

#   return(list(posterior, summary))
# }

# if (!file.exists("data/model_fits.rds")) {
#   if (!"cmdstanr" %in% installed.packages()) {
#     install.packages(
#       "cmdstanr",
#       repos = c('https://stan-dev.r-universe.dev', getOption("repos"))
#     )
#   }

#   # Compiling the two models -----------------------------------------------

#   efs_mod <- cmdstanr::cmdstan_model("models/efs.stan")
#   evt_mod <- cmdstanr::cmdstan_model("models/evt.stan")
#   evt_gumbel_mod <- cmdstanr::cmdstan_model("models/evt_gumbel.stan")

#   # Running the models in parallel --------------------------
#   future::plan(future::multisession)
#   model_fits <-
#     underlying_scenarios %>%
#     mutate(maxima_median = map_dbl(samples, ~ median(.x$top1))) |>
#     mutate(
#       stan_list_single = map2(.x = k, .y = samples, .f = get_stan_data)
#     ) %>%
#     mutate(
#       stan_list_multpl = map2(
#         .x = k,
#         .y = samples,
#         .f = ~ get_stan_data(.x, .y, multiple = TRUE)
#       )
#     ) %>%
#     mutate(
#       # stan_fit_evt = future_map2(
#       #   .x = stan_list_single,
#       #   .y = maxima_median,
#       #   .f = ~ fit_stan_mod(evt_mod, .x, init_func("evt", .y)),
#       #   .options = furrr_options(seed = TRUE)
#       # ),
#       stan_fit_evtg = future_map2(
#         .x = stan_list_single,
#         .y = maxima_median,
#         .f = ~ fit_stan_mod(evt_gumbel_mod, .x, init_func("evt", .y)),
#         .options = furrr_options(seed = TRUE)
#       ),
#       # stan_fit_efs = future_map2(
#       #   .x = stan_list_single,
#       #   .y = maxima_median,
#       #   .f = ~ fit_stan_mod(efs_mod, .x, init_func("efs", .y)),
#       #   .options = furrr_options(seed = TRUE)
#       # ),
#       # stan_fit_efsm = future_map2(
#       #   .x = stan_list_multpl,
#       #   .y = maxima_median,
#       #   .f = ~ fit_stan_mod(efs_mod, .x, init_func("efs", .y)),
#       #   .options = furrr_options(seed = TRUE)
#       # ) |>
#     ) |>
#     mutate(
#       evtg_posterior = map(stan_fit_evtg, ~ .x[[1]]),
#       evtg_summary = map(stan_fit_evtg, ~ .x[[2]]),
#       # evt_posterior = map(stan_fit_evt, ~ .x[[1]]),
#       # efs_posterior = map(stan_fit_efs, ~ .x[[1]]),
#       # efsm_posterior = map(stan_fit_efsm, ~ .x[[1]]),
#       # evt_summary = map(stan_fit_evt, ~ .x[[2]]),
#       # efs_summary = map(stan_fit_efs, ~ .x[[2]]),
#       # efsm_summary = map(stan_fit_efsm, ~ .x[[2]])
#     ) |>
#     select(!dplyr::contains("stan_fit_"))
#   future::plan(future::sequential)
# } else {
#   model_fits <- read_rds("data/model_fits.rds")
# }
