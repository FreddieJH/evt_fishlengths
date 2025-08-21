source("code/01_simulation.R")
source("code/02_model_fitting.R")
library(scales)
library(ggplot2)
library(RColorBrewer)

problematic_fits <-
    posteriors_summary |>
    filter(rhat > 1.1 | ess_bulk < 100 | ess_tail < 100) |>
    pull(filename)


for (selected_dist in c("tnorm", "gamma", "lnorm")) {
    for (model_type in c("efs_", "efsmult_")) {
        for (mu_prior in c(30, 50, 100, 150)) {
            plot_scenario_filename <-
                paste0(
                    "results/figures/p_simulation_bias/",
                    selected_dist,
                    "/",
                    model_type,
                    "/mu",
                    mu_prior,
                    ".png"
                )

            if (file.exists(plot_scenario_filename)) {
                next
            }

            gamma_priors <-
                scenario_data |>
                select(gamma_shape, gamma_rate) |>
                distinct() |>
                drop_na()

            p_simulation_bias <-
                posteriors_max_summary |>
                filter(str_detect(filename, model_type)) %>%
                left_join(
                    scenario_data |>
                        select(
                            filename,
                            dist_name,
                            mu_prior_mean,
                            gamma_mean,
                            gamma_variance
                        ) |>
                        distinct()
                ) |>
                filter(!(filename %in% problematic_fits)) |>
                filter(dist_name == selected_dist) %>%
                filter(mu_prior_mean == mu_prior) %>%
                mutate(
                    gamma_shape = as.numeric(str_extract(
                        filename,
                        "(?<=_shape)[0-9]+"
                    ))
                ) %>%
                left_join(gamma_priors) %>%
                mutate(
                    scenario_id = as.numeric(str_extract(filename, "\\d+$"))
                ) %>%
                left_join(scenarios_truemax) %>%
                left_join(
                    scenario_data |>
                        select(scenario_id, lambda, k) |>
                        distinct()
                ) |>
                mutate(
                    max_error = (est_max_fit - true_max) / true_max,
                    max_error_lwr = (est_max_lwr - true_max) / true_max,
                    max_error_upr = (est_max_upr - true_max) / true_max
                ) |>
                mutate(
                    gamma_prior = paste(
                        paste0(as.character(gamma_mean / 1000), "k"),
                        paste0(as.character(sqrt(gamma_variance) / 1000), "k"),
                        sep = ","
                    )
                ) |>
                # ggplot(aes(x = est_max_fit, y = max_error, col = as.factor(k))) +
                ggplot(aes(x = true_max, y = est_max_fit, col = as.factor(k))) +
                geom_abline(slope = 1) +
                geom_point() +
                # geom_errorbar(aes(ymin = max_error_lwr, ymax = max_error_upr)) +
                geom_errorbar(aes(ymin = est_max_lwr, ymax = est_max_upr)) +
                facet_grid(
                    lambda ~ gamma_prior,
                    scales = "free",
                    labeller = labeller(
                        lambda = function(x) paste0("n = ", x),
                        gamma_prior = function(x) paste0("Gamma(", x, ")")
                    )
                ) +
                theme_minimal(20) +
                theme(
                    panel.border = element_rect(
                        color = "black",
                        fill = NA,
                        linewidth = 1
                    ),
                    legend.position = "bottom",
                    legend.direction = "horizontal",
                    plot.background = element_rect(fill = "white")
                ) +
                guides(
                    color = guide_legend(
                        override.aes = list(size = 5),
                        nrow = 1
                    )
                ) +
                labs(
                    y = expression("Estimated " ~ L[max] ~ "(cm)"),
                    x = expression("True " ~ L[max] ~ "(cm)"),
                    color = "# samples (k):"
                )

            ggsave(
                filename = plot_scenario_filename,
                plot = p_simulation_bias,
                height = 10,
                width = 15
            )
        }
    }
}


# It appears that the choice of priors has not much difference on the outcome.
# We will therefore just take a single choice of priors

selected_scenarios <-
    scenario_data |>
    filter(gamma_mean == 8000, gamma_variance == 1000^2, mu_prior_mean == 50) |>
    pull(filename)


sel_pal <- brewer.pal(name = "Dark2", n = 3)

p_comparison <-
    posteriors_max_summary |>
    filter((filename %in% selected_scenarios) | str_detect(filename, "evt")) |>
    mutate(problematic = filename %in% problematic_fits) |>
    # filter(!(filename %in% problematic_fits)) |>
    left_join(
        scenario_data |>
            select(filename, dist_name, k, lambda, dist_mean, scenario_model) |>
            distinct()
    ) |>
    mutate(scenario_id = as.numeric(str_extract(filename, "\\d+$"))) |>
    left_join(scenarios_truemax) %>%
    mutate(
        pc_error = (true_max - est_max_fit) / true_max,
        pc_error_lwr = (true_max - est_max_lwr) / true_max,
        pc_error_upr = (true_max - est_max_upr) / true_max
    ) |>
    filter(lambda == 1000) |>
    mutate(
        fill_col = case_when(
            problematic ~ "white",
            dist_name == "gamma" ~ sel_pal[1],
            dist_name == "lnorm" ~ sel_pal[2],
            dist_name == "tnorm" ~ sel_pal[3]
        )
    ) |>
    ggplot(aes(true_max, pc_error, col = dist_name)) +
    geom_errorbar(aes(ymin = pc_error_lwr, ymax = pc_error_upr)) +
    geom_point(
        # fill = "white",
        aes(shape = dist_name, fill = fill_col, size = as.factor(k)),
        alpha = 1
    ) +
    geom_abline(slope = 0, lty = 2) +
    facet_grid(
        scenario_model ~ dist_mean,
        scales = "free",
        labeller = labeller(
            dist_mean = function(x) paste0("True mean = ", x, "cm"),
            scenario_model = function(x) {
                case_when(
                    x == "evt" ~ "EVT",
                    str_detect(x, "efs_") ~ "EFS",
                    str_detect(x, "efsmult_") ~ "EFSMM"
                )
            }
        )
    ) +
    labs(
        x = "True maximum length (cm)",
        y = "Estimation Error",
        shape = NULL,
        col = NULL,
        size = "k ="
    ) +
    # scale_shape_manual(
    #     values = c("TRUE" = 21, "FALSE" = 16),
    #     label = c("TRUE" = "Not converged", "FALSE" = "Converged")
    # ) +
    scale_y_continuous(label = label_percent()) +
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
    scale_fill_identity() +
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
p_comparison

ggsave(
    filename = "results/figures/performance.png",
    plot = p_comparison,
    height = 10,
    width = 10
)
