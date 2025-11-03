library(arrow)
library(ggplot2)
library(future)
library(furrr)

if (!exists("posterior")) {
  if (!file.exists("results/data/posterior.parquet")) {
    source("R/05_model_fitting.R")
  }
  posterior <- read_parquet("results/data/posterior.parquet")
}

#  the Rhat and ESS for non-convergence
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

p_bayes_check <-
  posterior_summary |>
  pivot_longer(cols = c(rhat, ess_bulk, ess_tail)) |>
  ggplot(aes(x = value, col = name)) +
  geom_density() +
  facet_wrap(~ model_id + name, scales = "free", ncol = 3) +
  theme_classic(20)

ggsave(
  p_bayes_check,
  filename = "results/figures/model_checks/bayes_check.png",
  height = 10,
  width = 16
)


# Plotting individual traceplots for each model run
plot_combinations <-
  posterior |>
  distinct(model_id, scenario_id) |>
  filter(
    !file.exists(paste0(
      "results/figures/model_checks/traceplots/",
      model_id,
      scenario_id,
      ".png"
    ))
  )

future::plan(multisession)

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
      filename = paste0(
        "results/figures/model_checks/traceplots/",
        m,
        i,
        ".png"
      ),
      plot = p,
      height = 2,
      width = 6
    )
  },
  .options = furrr_options(seed = TRUE)
)
future::plan(sequential)
