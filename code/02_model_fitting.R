# Script requirements  --------------------------------------------------------------

source("code/00a_truncnorm_functions.R")
source("code/00b_model_functions.R")
source("code/00c_expectedmax_functions.R")
source("code/01_simulation.R")
library(future)
library(furrr)
library(stringr)

# Note: If your models are not running try a complete reinstall of "cmdstanr"

# Expected maximum value of simulated samples -------------------------------------------------

# Expected  max is calcuated using the expected_max_fromsim function from "code/00c_expectedmax_functions.R"
# it is calculated as the expected value from the pdf of the maximum values
# Expected values are estimated as integral of x*g(x), where g(x) if the pdf of x

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


posterior_filename <- "data/model_output/posteriors7.parquet"
posterior_summary_filename <- "data/model_output/posteriors7_summary.parquet"

if (!file.exists(posterior_filename)) {
  future::plan(multisession)

  processed_data <-
    scenario_data |>
    nest(
      .by = c(
        filename,
        type,
        gamma_shape,
        gamma_rate,
        mu_prior_mean,
        mu_prior_sd
      )
    ) |>
    mutate(
      k = map_int(data, nrow),
      m = map(data, ~ .x$m),
      x = map(data, ~ pull(unnest(.x, cols = top_m), top_m)),
      n_obs = map_int(x, length),
      start_indx = map(m, ~ lag(c(0, cumsum(.x)))[-1] + 1)
    ) |>
    mutate(
      model_posterior = furrr::future_pmap(
        .l = list(
          maxima = x,
          k = k,
          type = type,
          gamma_shape = gamma_shape,
          gamma_rate = gamma_rate,
          mu_prior_mean = mu_prior_mean,
          mu_prior_sd = mu_prior_sd,
          n_obs = n_obs,
          n_per_sample = m,
          start_id = start_indx
        ),
        .f = fit_maxima_model,
        .options = furrr_options(seed = TRUE)
      )
    ) |>
    mutate(
      posterior_draws = map(model_posterior, get_posterior),
      posterior_summary = map(model_posterior, get_summary)
    )

  processed_data |>
    select(filename, posterior_draws) %>%
    unnest(cols = posterior_draws) |>
    arrow::write_parquet(posterior_filename)

  processed_data |>
    select(filename, posterior_summary) %>%
    unnest(cols = posterior_summary) |>
    arrow::write_parquet(posterior_summary_filename)
  future::plan(sequential)
}
posteriors <- arrow::read_parquet(posterior_filename)
posteriors_summary <- arrow::read_parquet(posterior_summary_filename)

posteriors_max_summary_path <- "data/model_summaries/posteriors7_max_summary.csv"
if (!file.exists(posteriors_max_summary_path)) {
  plan(multisession, workers = availableCores() - 1)

  posteriors_max_summary_evt <-
    posteriors |>
    left_join(scenario_data |> select(filename, k) |> distinct()) |>
    filter(stringr::str_detect(filename, "evt")) |>
    mutate(
      est_max = future_pmap_dbl(
        .l = list(p = 1 - (1 / k), loc = loc, scale = scale, shape = shape),
        .f = evd::qgev,
        .options = furrr_options(seed = TRUE)
      ),
      est_max20 = future_pmap_dbl(
        .l = list(p = 1 - (1 / 20), loc = loc, scale = scale, shape = shape),
        .f = evd::qgev,
        .options = furrr_options(seed = TRUE)
      )
    ) |>
    summarise(
      across(
        c(est_max, est_max20),
        list(
          fit = median,
          lwr = ~ quantile(.x, 0.1),
          upr = ~ quantile(.x, 0.9)
        )
      ),
      .by = filename
    )

  custom_funcs <- list(
    f_x = f_x,
    F_x = F_x,
    g_max = g_max,
    G_max = G_max,
    # inverse_G_x = inverse_G_x,
    ptnorm = ptnorm,
    qtnorm = qtnorm,
    dtnorm = dtnorm
  )

  posteriors_max_summary_efs <-
    posteriors |>
    filter(stringr::str_detect(filename, "efs_")) |>
    left_join(scenario_data |> select(filename, k) |> distinct()) |>
    mutate(
      est_max = future_pmap_dbl(
        .l = list(
          distr = "tnorm",
          par1 = mu,
          par2 = sigma,
          n = lambda * k
        ),
        .f = expected_max,
        .options = furrr_options(
          seed = TRUE,
          globals = custom_funcs
        )
      ),
      est_max20 = future_pmap_dbl(
        .l = list(
          distr = "tnorm",
          par1 = mu,
          par2 = sigma,
          n = lambda * 20
        ),
        .f = expected_max,
        .options = furrr_options(
          seed = TRUE,
          globals = custom_funcs
        )
      )
    ) |>
    summarise(
      across(
        c(est_max, est_max20),
        list(
          fit = median,
          lwr = ~ quantile(.x, 0.1),
          upr = ~ quantile(.x, 0.9)
        )
      ),
      .by = filename
    )

  posteriors_max_summary_efsmult <-
    posteriors |>
    filter(stringr::str_detect(filename, "efsmult_")) |>
    left_join(scenario_data |> select(filename, k) |> distinct()) |>
    mutate(
      est_max = future_pmap_dbl(
        .l = list(
          distr = "tnorm",
          par1 = mu,
          par2 = sigma,
          n = lambda * k
        ),
        .f = expected_max,
        .options = furrr_options(
          seed = TRUE,
          globals = custom_funcs
        )
      ),
      est_max20 = future_pmap_dbl(
        .l = list(
          distr = "tnorm",
          par1 = mu,
          par2 = sigma,
          n = lambda * 20
        ),
        .f = expected_max,
        .options = furrr_options(
          seed = TRUE,
          globals = custom_funcs
        )
      )
    ) |>
    summarise(
      across(
        c(est_max, est_max20),
        list(
          fit = median,
          lwr = ~ quantile(.x, 0.1),
          upr = ~ quantile(.x, 0.9)
        )
      ),
      .by = filename
    )

  posteriors_max_summary <- bind_rows(
    posteriors_max_summary_evt,
    posteriors_max_summary_efs,
    posteriors_max_summary_efsmult
  )

  readr::write_csv(posteriors_max_summary, posteriors_max_summary_path)
} else {
  posteriors_max_summary <- readr::read_csv(
    posteriors_max_summary_path,
    show_col_types = FALSE
  )
}
