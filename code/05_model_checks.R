source("code/04_model_fitting.R")

# check the Rhat and ESS for non-convergence

posterior_summary <-
  posterior |>
  nest(.by = c(scenario_id, model_id)) |>
  mutate(
    data_wide = map(
      data,
      ~ pivot_wider(.x, names_from = par, values_from = value)
    )
  ) |>
  mutate(posterior_summary = map(data_wide, posterior::summarise_draws)) |>
  select(scenario_id, model_id, posterior_summary) |>
  unnest(posterior_summary)

if (!file.exists("results/figures/bayes_check.png")) {
  posterior_summary <-
    bind_rows(
      read_parquet("data/results_single.parquet"),
      read_parquet("data/results_multpl.parquet")
    ) |>
    nest(.by = c(scenario_id, model_id)) |>
    mutate(
      data_wide = map(
        data,
        ~ pivot_wider(.x, names_from = par, values_from = value)
      )
    ) |>
    mutate(posterior_summary = map(data_wide, posterior::summarise_draws)) |>
    select(scenario_id, model_id, posterior_summary) |>
    unnest(posterior_summary)

  p_bayes_check <-
    posterior_summary |>
    pivot_longer(cols = c(rhat, ess_bulk, ess_tail)) |>
    ggplot(aes(x = value, col = name)) +
    geom_density() +
    facet_wrap(~ model_id + name, scales = "free", ncol = 3) +
    theme_classic(20)

  ggsave(
    p_bayes_check,
    filename = "results/figures/bayes_check.png",
    height = 16,
    width = 10
  )
}

# Plotting individual traceplots for each model run
future::plan(multisession)

plot_combinations <-
  posterior |>
  distinct(model_id, scenario_id) |>
  filter(
    !file.exists(paste0(
      "results/figures/traceplots/",
      model_id,
      scenario_id,
      ".png"
    ))
  )

future_walk2(
  plot_combinations$model_id,
  plot_combinations$scenario_id,
  \(m, i) {
    p <- posterior |>
      filter(scenario_id == i, model_id == m) |>
      ggplot(aes(
        .iteration,
        value,
        colour = as.factor(.chain),
        group = .chain
      )) +
      geom_path(alpha = 0.4) +
      facet_wrap(~par, scales = "free") +
      theme_classic() +
      theme(legend.position = "none")

    ggsave(
      filename = paste0("results/figures/traceplots/", m, i, ".png"),
      plot = p,
      height = 2,
      width = 6
    )
  },
  .options = furrr_options(seed = TRUE)
)
future::plan(sequential)
