# Figure 1: Conceptual figure comparing the two methods of estimating maximum length of a fish

# Script requirements  --------------------------------------------------------------

source("code/01_simulation.R")
source("code/00b_model_functions.R")
source("code/00c_expectedmax_functions.R")
source("code/00d_colours.R")
source("code/02_model_fitting.R")
library(tidyverse)
library(evd)
library(geomtextpath)
library(patchwork)

# Script parameters  --------------------------------------------------------------

# Choosing scenario parameters for the example
scenario1_k <- 5
scenario1_lambda <- 1000
scenario1_mean <- 50
scenario1_variance <- (scenario1_mean * 0.34)^2

# Scenario calculations  --------------------------------------------------------------

# getting the single scenario based on the pars set above
scenario1 <-
    scenario_data %>%
    filter(
        k == scenario1_k,
        lambda == scenario1_lambda,
        dist_mean == scenario1_mean,
        dist_name == "tnorm",
    ) %>%
    filter(str_detect(filename, "efsmult_shape16_mu50"))


# scenario maxima, list elements for each sample
# more than one maxima per sample in some cases
scenario1_maxima <-
    scenario1 %>%
    pull(top_m)

# the n for each sample (k samples)
scenario1_n <-
    scenario1 |>
    pull(n)

scenario1_posterior_max_summary <-
    posteriors_max_summary %>%
    mutate(scenario_id = as.numeric(str_extract(filename, "\\d+$"))) %>%
    filter(scenario_id == unique(scenario1$scenario_id)) |>
    filter(str_detect(filename, "evt_|shape16_mu50"))

scenario1_posterior_max_summary_evt <-
    scenario1_posterior_max_summary %>%
    filter(str_detect(filename, "evt_"))

scenario1_posterior_max_summary_efs <-
    scenario1_posterior_max_summary %>%
    filter(str_detect(filename, "efs_shape16_mu50"))

scenario1_posterior_max_summary_efsmult <-
    scenario1_posterior_max_summary %>%
    filter(str_detect(filename, "efsmult_shape16_mu50"))

evt_posteriors <-
    posteriors %>%
    mutate(scenario_id = as.numeric(str_extract(filename, "\\d+$"))) %>%
    filter(scenario_id == unique(scenario1$scenario_id)) |>
    filter(str_detect(filename, "evt_"))

evt_posteriors |>
    summarise(
        med_loc = quantile(loc, 0.5),
        med_scale = quantile(scale, 0.5),
        med_shape = quantile(shape, 0.5),
        lwr_loc = quantile(loc, 0.1),
        lwr_scale = quantile(scale, 0.1),
        lwr_shape = quantile(shape, 0.1),
        upr_loc = quantile(loc, 0.9),
        upr_scale = quantile(scale, 0.9),
        upr_shape = quantile(shape, 0.9)
    )

evt_posteriors |>
    mutate(
        qgev = pmap_dbl(
            .l = list(p = 0.8, loc = loc, scale = scale, shape = shape),
            .f = evd::qgev
        )
    ) |>
    summarise(
        qgev_fit = quantile(qgev, 0.5),
        qgev_lwr = quantile(qgev, 0.1),
        qgev_upr = quantile(qgev, 0.9)
    )

evt_posteriors |>
    mutate(
        qgev = pmap_dbl(
            .l = list(p = 0.95, loc = loc, scale = scale, shape = shape),
            .f = evd::qgev
        )
    ) |>
    summarise(
        qgev_fit = quantile(qgev, 0.5),
        qgev_lwr = quantile(qgev, 0.1),
        qgev_upr = quantile(qgev, 0.9)
    )

expected_max_fromsim(
    "tnorm",
    n = scenario1_k * scenario1_lambda,
    mean = scenario1_mean,
    variance = (scenario1_mean * 0.34)^2
)
expected_max_fromsim(
    "tnorm",
    n = 20 * scenario1_lambda,
    mean = scenario1_mean,
    variance = (scenario1_mean * 0.34)^2
)


efs_posteriors <-
    posteriors %>%
    mutate(scenario_id = as.numeric(str_extract(filename, "\\d+$"))) %>%
    filter(scenario_id == unique(scenario1$scenario_id)) |>
    filter(str_detect(filename, "efs_shape16_mu50_11"))

efs_posteriors |>
    summarise(
        med_mu = quantile(mu, 0.5),
        med_sigma = quantile(sigma, 0.5),
        med_lambda = quantile(lambda, 0.5),
        lwr_mu = quantile(mu, 0.1),
        lwr_sigma = quantile(sigma, 0.1),
        lwr_lambda = quantile(lambda, 0.1),
        upr_mu = quantile(mu, 0.9),
        upr_sigma = quantile(sigma, 0.9),
        upr_lambda = quantile(lambda, 0.9)
    )

efs_posteriors |>
    mutate(
        exp_max = pmap_dbl(
            .l = list(
                distr = "tnorm",
                n = scenario1_k * lambda,
                par1 = mu,
                par2 = sigma
            ),
            .f = expected_max
        )
    ) |>
    summarise(
        max_fit = quantile(exp_max, 0.5),
        max_lwr = quantile(exp_max, 0.1),
        max_upr = quantile(exp_max, 0.9)
    )

efs_posteriors |>
    mutate(
        exp_max = pmap_dbl(
            .l = list(
                distr = "tnorm",
                n = 20 * lambda,
                par1 = mu,
                par2 = sigma
            ),
            .f = expected_max
        )
    ) |>
    summarise(
        max_fit = quantile(exp_max, 0.5),
        max_lwr = quantile(exp_max, 0.1),
        max_upr = quantile(exp_max, 0.9)
    )


efsmult_posteriors <-
    posteriors %>%
    mutate(scenario_id = as.numeric(str_extract(filename, "\\d+$"))) %>%
    filter(scenario_id == unique(scenario1$scenario_id)) |>
    filter(str_detect(filename, "efsmult_shape16_mu50_11"))

efsmult_posteriors |>
    summarise(
        med_mu = quantile(mu, 0.5),
        med_sigma = quantile(sigma, 0.5),
        med_lambda = quantile(lambda, 0.5),
        lwr_mu = quantile(mu, 0.1),
        lwr_sigma = quantile(sigma, 0.1),
        lwr_lambda = quantile(lambda, 0.1),
        upr_mu = quantile(mu, 0.9),
        upr_sigma = quantile(sigma, 0.9),
        upr_lambda = quantile(lambda, 0.9)
    )

efsmult_posteriors |>
    mutate(
        exp_max = pmap_dbl(
            .l = list(
                distr = "tnorm",
                n = scenario1_k * lambda,
                par1 = mu,
                par2 = sigma
            ),
            .f = expected_max
        )
    ) |>
    summarise(
        max_fit = quantile(exp_max, 0.5),
        max_lwr = quantile(exp_max, 0.1),
        max_upr = quantile(exp_max, 0.9)
    )

efsmult_posteriors |>
    mutate(
        exp_max = pmap_dbl(
            .l = list(
                distr = "tnorm",
                n = 20 * lambda,
                par1 = mu,
                par2 = sigma
            ),
            .f = expected_max
        )
    ) |>
    summarise(
        max_fit = quantile(exp_max, 0.5),
        max_lwr = quantile(exp_max, 0.1),
        max_upr = quantile(exp_max, 0.9)
    )

scenario1_evt_cdf_summary_path <- "data/model_summaries/scenario1_evt_cdf_summary3.csv"
if (!file.exists(scenario1_evt_cdf_summary_path)) {
    scenario1_evt_cdf_summary <-
        posteriors %>%
        mutate(scenario_id = as.numeric(str_extract(filename, "\\d+$"))) %>%
        filter(scenario_id == unique(scenario1$scenario_id)) |>
        filter(str_detect(filename, "evt_")) |>
        expand_grid(size = 0:180) |>
        mutate(
            cdf = pmap_dbl(
                .l = list(q = size, loc = loc, scale = scale, shape = shape),
                .f = evd::pgev
            )
        ) |>
        summarise(
            p_fit = quantile(cdf, 0.5),
            p_lwr = quantile(cdf, 0.1),
            p_upr = quantile(cdf, 0.9),
            .by = size
        )

    write_csv(scenario1_evt_cdf_summary, scenario1_evt_cdf_summary_path)
} else {
    scenario1_evt_cdf_summary <- read_csv(
        scenario1_evt_cdf_summary_path,
        show_col_types = FALSE
    )
}

scenario1_efs_summary_path <- "data/model_summaries/scenario1_efs_cdf_summary3.csv"
if (!file.exists(scenario1_efs_summary_path)) {
    scenario1_efs_cdf_summary <-
        posteriors %>%
        mutate(scenario_id = as.numeric(str_extract(filename, "\\d+$"))) %>%
        filter(scenario_id == unique(scenario1$scenario_id)) |>
        filter(str_detect(filename, "efs_shape16_mu50")) |>
        expand_grid(size = 0:180) |>
        mutate(
            cdf = pmap_dbl(
                .l = list(
                    x = size,
                    distr = "tnorm",
                    n = lambda * scenario1_k,
                    par1 = mu,
                    par2 = sigma
                ),
                .f = G_max
            ),
            cdf20 = pmap_dbl(
                .l = list(
                    x = size,
                    distr = "tnorm",
                    n = lambda * 20,
                    par1 = mu,
                    par2 = sigma
                ),
                .f = G_max
            )
        ) |>
        summarise(
            p_fit = quantile(cdf, 0.5),
            p_lwr = quantile(cdf, 0.1),
            p_upr = quantile(cdf, 0.9),
            p_fit20 = quantile(cdf20, 0.5),
            p_lwr20 = quantile(cdf20, 0.1),
            p_upr20 = quantile(cdf20, 0.9),
            .by = size
        )

    write_csv(scenario1_efs_cdf_summary, scenario1_efs_summary_path)
} else {
    scenario1_efs_cdf_summary <- read_csv(
        scenario1_efs_summary_path,
        show_col_types = FALSE
    )
}

scenario1_efsmult_summary_path <- "data/model_summaries/scenario1_efsmult_cdf_summary3.csv"
if (!file.exists(scenario1_efsmult_summary_path)) {
    scenario1_efsmult_cdf_summary <-
        posteriors %>%
        mutate(scenario_id = as.numeric(str_extract(filename, "\\d+$"))) %>%
        filter(scenario_id == unique(scenario1$scenario_id)) |>
        filter(str_detect(filename, "efsmult_shape16_mu50")) |>
        expand_grid(size = 0:180) |>
        mutate(
            cdf = pmap_dbl(
                .l = list(
                    x = size,
                    distr = "tnorm",
                    n = lambda * scenario1_k,
                    par1 = mu,
                    par2 = sigma
                ),
                .f = G_max
            ),
            cdf20 = pmap_dbl(
                .l = list(
                    x = size,
                    distr = "tnorm",
                    n = lambda * 20,
                    par1 = mu,
                    par2 = sigma
                ),
                .f = G_max
            )
        ) |>
        summarise(
            p_fit = quantile(cdf, 0.5),
            p_lwr = quantile(cdf, 0.1),
            p_upr = quantile(cdf, 0.9),
            p_fit20 = quantile(cdf20, 0.5),
            p_lwr20 = quantile(cdf20, 0.1),
            p_upr20 = quantile(cdf20, 0.9),
            .by = size
        )

    write_csv(scenario1_efsmult_cdf_summary, scenario1_efsmult_summary_path)
} else {
    scenario1_efsmult_cdf_summary <- read_csv(
        scenario1_efsmult_summary_path,
        show_col_types = FALSE
    )
}

# expected maximum based on 20 samples, each of length 200 (total is 20*1000 samples)
scenario1_max <- expected_max_fromsim(
    dist = "tnorm",
    n = scenario1_k * scenario1_lambda,
    mean = scenario1_mean,
    variance = scenario1_variance
)
scenario1_max20 <- expected_max_fromsim(
    dist = "tnorm",
    n = (20 * scenario1_lambda),
    mean = scenario1_mean,
    variance = scenario1_variance
)

# Figure 2  ------------------------------------------------------------------------------

p_scenario1_underlying <-
    ggplot() +
    geom_textdensity(
        aes(x),
        data = tibble(
            x = rtnorm(
                1e6,
                a = 0,
                mean = scenario1_mean,
                sd = sqrt(scenario1_variance)
            )
        ),
        label = "Observable body size distribution",
        hjust = 0.4
    ) +
    geom_rect(
        aes(
            xmin = scenario1_posterior_max_summary_evt$est_max20_lwr,
            xmax = scenario1_posterior_max_summary_evt$est_max20_upr,
            ymin = -Inf,
            ymax = Inf
        ),
        fill = evt_colour,
        alpha = 0.3
    ) +
    geom_rect(
        aes(
            xmin = scenario1_posterior_max_summary_efs$est_max20_lwr,
            xmax = scenario1_posterior_max_summary_efs$est_max20_upr,
            ymin = -Inf,
            ymax = Inf
        ),
        fill = efs_colour,
        alpha = 0.3
    ) +
    geom_rect(
        aes(
            xmin = scenario1_posterior_max_summary_efsmult$est_max20_lwr,
            xmax = scenario1_posterior_max_summary_efsmult$est_max20_upr,
            ymin = -Inf,
            ymax = Inf
        ),
        fill = efsmult_colour,
        alpha = 0.5
    ) +
    geom_textvline(
        aes(xintercept = scenario1_posterior_max_summary_efs$est_max20_fit),
        color = efs_colour_dark,
        linetype = "dashed",
        label = as.character(expression(L[paste(max, ", ", 20)])),
        parse = TRUE,
        size = 5
    ) +
    geom_textvline(
        aes(xintercept = scenario1_posterior_max_summary_evt$est_max20_fit),
        color = evt_colour,
        linetype = "dashed",
        label = as.character(expression(L[paste(max, ", ", 20)])),
        parse = TRUE,
        size = 5
    ) +
    geom_textvline(
        aes(xintercept = scenario1_posterior_max_summary_efsmult$est_max20_fit),
        color = efsmult_colour_dark,
        linetype = "dashed",
        label = as.character(expression(L[paste(max, ", ", 20)])),
        parse = TRUE,
        size = 5
    ) +
    geom_rug(
        aes(x = unlist(scenario1_maxima)),
        color = data_colour_nearmax,
        length = unit(0.5, units = "cm")
    ) +
    geom_rug(
        aes(x = unlist(map(scenario1_maxima, max))),
        color = data_colour_ismax,
        length = unit(0.5, units = "cm")
    ) +
    geom_point(aes(x = scenario1_max20, y = 0), size = 4) +
    labs(x = "Fish body length (cm)", y = NULL) +
    theme_minimal(20) +
    theme(
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
    ) +
    scale_x_continuous(breaks = seq(40, 150, 20), expand = c(0, 1)) +
    coord_cartesian(xlim = c(40, 150))

p1_lims <- layer_scales(p_scenario1_underlying)$x$range$range

scenario1_plotting <-
    scenario1 |>
    mutate(maxima = unlist(map(top_m, max))) |>
    mutate(plotting_pos = rank(maxima) / (max(rank(maxima)) + 1))

scenario1_unnest <-
    scenario1_plotting |>
    unnest(top_m)

# scenario1_unnest |>
#     mutate(
#         max = max(top_m),
#         near_max = rank(-top_m, ties.method = "max") <= m,
#         .by = j
#     )

# tibble(
#     sample_id = scenario1_unnest$j,
#     value = scenario1_unnest$top_m,
#     is_max = sample_data == max(sample_data),
#     near_max = rank(-value, ties.method = "max") <= threshold,
#     n = nn, # Store sample size for reference
#     rank = case_when(
#         !is_max & !near_max ~ "low",
#         !is_max & near_max ~ "high",
#         is_max & near_max ~ "highest"
#     )
# )

plot_modfit <- function(type, colour, k = 20) {
    posterior <- get(paste0("scenario1_", type, "_cdf_summary"))
    est_summary <- get(paste0("scenario1_posterior_max_summary_", type))

    est_max20_fit <- est_summary$est_max20_fit
    est_max20_lwr <- est_summary$est_max20_lwr
    est_max20_upr <- est_summary$est_max20_upr

    pc <- if (type == "evt") {
        1 - (1 / k)
    } else {
        posterior |>
            mutate(diff = abs(size - est_max20_fit)) |>
            arrange(diff) |>
            head(1) |>
            pull(p_fit20)
    }

    posterior_out <- if (type == "evt") {
        posterior
    } else {
        posterior |>
            select(size, p_fit = p_fit20, p_lwr = p_lwr20, p_upr = p_upr20)
    }

    ggplot(posterior_out) +
        {
            if (type %in% c("efs", "efsmult")) {
                geom_ribbon(
                    aes(x = size, ymin = p_lwr, ymax = p_upr),
                    fill = colour,
                    alpha = 0.2,
                    data = posterior
                )
            }
        } +
        {
            if (type %in% c("efs", "efsmult")) {
                geom_line(
                    aes(x = size, y = p_fit),
                    col = colour,
                    lty = 2,
                    alpha = 0.5,
                    data = posterior
                )
            }
        } +
        geom_ribbon(
            aes(x = size, ymin = p_lwr, ymax = p_upr),
            fill = colour,
            alpha = 0.5
        ) +
        {
            if (type == "efsmult") {
                geom_point(
                    aes(x = top_m, y = 0, size = n_scaled),
                    data = scenario1_plotting |>
                        mutate(n_scaled = (n / 600)^3) |>
                        unnest(top_m) |>
                        filter(top_m != maxima),
                    color = data_colour_nearmax
                )
            }
        } +
        {
            if (type == "evt") {
                geom_point(
                    aes(x = maxima, y = plotting_pos, size = n_scaled),
                    data = scenario1_plotting |> mutate(n_scaled = (n / 600)^3),
                    color = data_colour_ismax
                )
            }
        } +
        {
            if (type != "evt") {
                geom_point(
                    aes(x = maxima, y = 0, size = n_scaled),
                    data = scenario1_plotting |> mutate(n_scaled = (n / 600)^3),
                    color = data_colour_ismax
                )
            }
        } +
        # geom_rug(
        #     aes(x = unlist(map(scenario1_maxima, max))),
        #     color = "purple",
        #     length = unit(1, units = "cm"),
        #     linewidth = 2,
        #     data = tibble()
        # ) +
        geom_line(aes(x = size, y = p_fit), col = colour, linewidth = 1.5) +
        geom_errorbarh(
            aes(xmin = est_max20_lwr, xmax = est_max20_upr, y = pc),
            height = 0.03,
            col = "black",
            lty = "solid",
            data = tibble(n = 1)
        ) +
        geom_vline(xintercept = est_max20_fit, col = colour, lty = "dashed") +
        labs(
            x = "Fish body length (cm)",
            y = "Pr(Lmax < size)",
            size = "Sample size"
        ) +
        theme_minimal(20) +
        theme(
            legend.position = "inside",
            legend.position.inside = c(0.05, 1),
            legend.justification = c(0, 1.2),
            legend.background = element_rect(fill = "white", color = "black"),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank()
        ) +
        scale_x_continuous(breaks = seq(40, 150, 20), expand = c(0, 1)) +
        scale_size_identity() +
        coord_cartesian(xlim = c(40, 150))
}

p_scenario1_evt <- plot_modfit(type = "evt", colour = evt_colour)
p_scenario1_efs <- plot_modfit(type = "efs", colour = efs_colour)
p_scenario1_efsmult <- plot_modfit(type = "efsmult", colour = efsmult_colour)

p_examplefit <-
    p_scenario1_underlying +
    p_scenario1_evt +
    p_scenario1_efs +
    theme(legend.position = "none") +
    p_scenario1_efsmult +
    theme(legend.position = "none") +
    plot_layout(ncol = 1, axis_titles = "collect") +
    plot_annotation(tag_levels = "A")

ggsave(
    filename = "results/figures/p_examplefit.png",
    plot = p_examplefit,
    height = 12,
    width = 10
)
