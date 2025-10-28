library(dplyr)
library(readr)
library(ggplot2)
library(forcats)
library(scales)

if (!exists("estmax_posterior")) {
  estmax_posterior <- read_csv(
    "results/data/estmax_posterior.csv",
    show_col_types = FALSE
  )
}

if (!exists("scenarios_truemax")) {
  scenarios_truemax <- read_csv(
    "results/data/scenarios_truemax.csv",
    show_col_types = FALSE
  )
}


bias_summary <-
  estmax_posterior |>
  left_join(scenarios_truemax) |>
  mutate(
    est_max_fit_err = (est_max_fit - true_max) / true_max,
    est_max_lwr_err = (est_max_lwr - true_max) / true_max,
    est_max_upr_err = (est_max_upr - true_max) / true_max
  ) |>
  mutate(
    model_id = forcats::fct_relevel(model_id, c("evt", "evtg", "efs", "efsm"))
  ) |>
  summarise(
    med_err = median(est_max_fit_err),
    .by = c(model_id, dist_name, lambda)
  )

bias_summary |>
  filter(dist_name == "lnorm") |>
  mutate(med_err_pc = med_err * 100)


sensitivity_data <-
  estmax_posterior |>
  left_join(scenarios_truemax) |>
  filter(lambda == 1000) |>
  # mutate(est_max_fit = coalesce(q5_fit, modal_fit),est_max_lwr = coalesce(q5_lwr, modal_lwr), est_max_upr = coalesce(q5_upr, modal_upr)) |>
  mutate(
    est_max_fit_err = (est_max_fit - true_max) / true_max,
    est_max_lwr_err = (est_max_lwr - true_max) / true_max,
    est_max_upr_err = (est_max_upr - true_max) / true_max
  ) |>
  mutate(
    model_id = forcats::fct_relevel(model_id, c("evt", "evtg", "efs", "efsm"))
  )


sensitivity_data |>
  summarise(med_err = median(est_max_fit_err), .by = c(model_id, dist_name))

sensitivity_data |>
  mutate(ci = q95_upr - q95_lwr) |>
  summarise(med_ci = median(ci), .by = c(model_id, dist_name))

p_comparison <-
  sensitivity_data |>
  ggplot(aes(true_max, est_max_fit_err, col = dist_name)) +
  geom_abline(slope = 0, lty = 2) +
  geom_errorbar(aes(
    ymin = est_max_lwr_err,
    ymax = est_max_upr_err,
    width = true_max / 30
  )) +
  geom_point(
    aes(shape = dist_name, fill = dist_name, size = as.factor(k)),
    alpha = 1
  ) +

  facet_grid(
    model_id ~ dist_mean,
    scales = "free",
    labeller = labeller(
      dist_mean = function(x) paste0("True mean = ", x, "cm"),
      model_id = function(x) {
        case_when(
          x == "evt" ~ "EVT (GEV)",
          x == "evtg" ~ "EVT (Gumbel)",
          x == "efs" ~ "EFS",
          x == "efsm" ~ "EFSMM"
        )
      }
    )
  ) +
  labs(
    x = "True maximum length (cm)",
    y = "Estimation Error",
    shape = NULL,
    col = NULL,
    fill = NULL,
    size = "k ="
  ) +
  scale_y_continuous(label = label_percent(), limits = c(-0.25, 0.25)) +
  theme_classic(20) +
  theme(
    panel.border = element_rect(
      colour = "black",
      fill = NA,
      linewidth = 0.5
    ),
    axis.line = element_blank(),
    legend.position = "bottom"
  ) +
  scale_color_manual(values = brewer.pal(name = "Dark2", n = 3)) +
  scale_fill_manual(values = brewer.pal(name = "Dark2", n = 3)) +
  scale_shape_manual(values = c(21, 22, 24)) +
  guides(
    col = guide_legend(
      nrow = 1,
      byrow = TRUE,
      override.aes = list(size = 5, linewidth = 0)
    ),
    size = guide_legend(nrow = 1, byrow = TRUE),
    shape = guide_legend(
      nrow = 2,
      byrow = TRUE,
      override.aes = list(size = 5)
    )
  )

ggsave(
  filename = "results/figures/manuscript_figures/sensitivity.png",
  plot = p_comparison,
  height = 15,
  width = 15
)
