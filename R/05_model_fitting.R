# Script requirements  --------------------------------------------------------------
source("R/01_funcs.R")
source("R/04_model_prep.R")

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
      .f = function(model, model_id, data_list) {
        model$sample(
          data = data_list,
          iter_warmup = 2000,
          iter_sampling = 1000,
          chains = 4,
          parallel_chains = 4,  # use all chains in parallel per model
          refresh = 1000,
          init = init_func(model_id, median(data_list$x)),
          adapt_delta = 0.999,
          max_treedepth = 12
        ) |> 
        posterior::as_draws_df()
      },
      .options = furrr_options(seed = TRUE)
    )
  ) %>%
  select(scenario_id, model_id, fit) %>%
  mutate(
    fit = map(fit, ~ .x %>%
      as_tibble() %>%
      select(-lp__) %>%
      pivot_longer(-c(.chain, .iteration, .draw), names_to = "par"))
  ) %>%
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
      select(scenario_id, stan_list_multpl) |>
      crossing(tibble(
        model_id = c("efsm"),
        model = list(efs_mod)
      )) %>%
      mutate(
        fit = future_pmap(
          .l = list(model, model_id, stan_list_multpl),
          .f = function(model, model_id, data_list) {
            model$sample(
            data = data_list,
            iter_warmup = 2000,
            iter_sampling = 1000,
            chains = 4,
            refresh = 1000,
            parallel_chains = 1,
            init = init_func(model_id, median(data_list$x)), # avoiding looking in complete wrong place
            adapt_delta = 0.999, # to avoid the small number of divergent transitions
            max_treedepth = 12 # increased from 10 to 12 after increasing the adapt_delta
          ) |> 
            posterior::as_draws_df()
        },
      .options = furrr_options(seed = TRUE)
          ) 
        )|> 
            select(scenario_id, model_id, fit) %>%
  mutate(
    fit = map(fit, ~ .x %>%
      as_tibble() %>%
      select(-lp__) %>%
      pivot_longer(-c(.chain, .iteration, .draw), names_to = "par"))
  ) %>%
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
