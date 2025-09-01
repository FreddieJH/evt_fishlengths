# Figure 1: Conceptual figure comparing the two methods of estimating maximum length of a fish

# Script requirements  --------------------------------------------------------------

# source("code/00a_truncnorm_functions.R")
# source("code/00b_model_functions.R")
# source("code/00c_expectedmax_functions.R")

source("code/00_funcs.R")
source("code/00d_colours.R")

library(ggplot2)
library(patchwork)
library(evd)

# Script parameters  --------------------------------------------------------------

concept_popln_mean <- 50
concept_popln_sd <- 16
concept_k <- 6 # number of samples (i.e., maxima values)
concept_n <- 50 # sample size per sample
save_figures <- TRUE


# Underlying population  -------------------------------------------------

# Random sampling of truncated normal distribution
set.seed(1)
concept_popln <- rtnorm(10000, concept_popln_mean, concept_popln_sd)

# PDF of the underlying population
concept_popln_pdf <-
    tibble(size = 0:150) |>
    mutate(p = dtnorm(size, concept_popln_mean, concept_popln_sd))

# get k maxima values from the underlying distribution
concept_maxima <- replicate(
    concept_k,
    max(sample(x = concept_popln, size = concept_n))
)

# Sampling method --------------------------------------------------------

samples_df <- map_dfr(1:concept_k, function(i) {
    n_k <- sample(10:100, 1)
    x <- rtnorm(n_k, concept_popln_mean, concept_popln_sd)
    m <- sample(1:4, 1)

    # Create and return data frame for this sample
    tibble(
        sample_id = i,
        value = x,
        is_max = x == max(x),
        near_max = rank(-value, ties.method = "max") <= m,
        n = n_k, # Store sample size for reference
        rank = case_when(
            !is_max & !near_max ~ "low",
            !is_max & near_max ~ "high",
            is_max & near_max ~ "highest"
        )
    )
})

concept_samples <-
    ggplot(samples_df, aes(x = value)) +
    geom_dotplot(
        aes(fill = rank),
        method = "histodot",
        # binpositions = "all",
        dotsize = 1,
        binwidth = 5
    ) +
    scale_fill_manual(
        values = c(
            "low" = "transparent",
            "high" = data_colour_nearmax,
            "highest" = data_colour_ismax
        )
    ) +
    facet_wrap(
        ~sample_id,
        nrow = 1,
        labeller = labeller(sample_id = function(value) {
            paste0("Sample #", value)
        })
    ) +
    labs(x = "Body Size (cm)", y = "Count") +
    theme_classic(base_size = 20) +
    theme(
        legend.position = "none",
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_blank()
    )

# cleaning data for modelling
get_model_data <- function(model_type) {
    filter_var <- ifelse(model_type == "efsmult", "near_max", "is_max")

    filtered_data <- samples_df %>% filter(.data[[filter_var]])
    sample_counts <- filtered_data %>% count(sample_id, name = "n")

    list(
        maxima = filtered_data %>% pull(value),
        type = model_type,
        n_obs = nrow(filtered_data),
        n_per_sample = sample_counts %>% pull(n),
        start_id = c(1, cumsum(sample_counts$n)[-nrow(sample_counts)] + 1),
        gamma_shape = 5,
        gamma_rate = 0.1,
        mu_prior_mean = 30,
        mu_prior_sd = 30 / 1.5,
        k = length(unique(samples_df$sample_id))
    )
}

get_concept_expmax <- function(model_type, k) {
    get(paste0("concept_", model_type, "_fit")) |>
        get_posterior() |>
        mutate(
            pdf = if (model_type == "evt") {
                pmap_dbl(
                    list(loc = loc, scale = scale, shape = shape, k = k),
                    expected_max_evt
                )
            } else {
                pmap_dbl(
                    .l = list(
                        distr = "tnorm",
                        n = k * lambda,
                        mean = mu,
                        variance = sigma^2
                    ),
                    .f = expected_max_fromsim
                )
            }
        ) |>
        summarise(
            p_fit = mean(pdf)
        ) %>%
        as.numeric()
}

# Model fitting
for (model_type in c("evt", "efs", "efsmult")) {
    # model fitting
    assign(
        x = paste0("concept_", model_type, "_fit"),
        value = do.call(fit_maxima_model, get_model_data(model_type))
    )

    # pdf summary
    assign(
        x = paste0("concept_", model_type, "_pdf_summary"),
        value = get_posterior(get(paste0("concept_", model_type, "_fit"))) |>
            expand_grid(size = 0:180) |>
            mutate(
                pdf = if (model_type == "evt") {
                    pmap_dbl(
                        list(x = size, loc = loc, scale = scale, shape = shape),
                        evd::dgev
                    )
                } else {
                    pmap_dbl(list(x = size, mean = mu, sd = sigma), dtnorm)
                }
            ) |>
            summarise(
                p_fit = mean(pdf),
                p_lwr = quantile(pdf, 0.1),
                p_upr = quantile(pdf, 0.9),
                .by = size
            )
    )

    # pdf of max
    assign(
        x = paste0("concept_", model_type, "_pdfmax_summary"),
        if (model_type != "evt") {
            get_posterior(get(paste0("concept_", model_type, "_fit"))) |>
                expand_grid(size = 0:180) |>
                mutate(
                    pdf = pmap_dbl(
                        .l = list(
                            x = size,
                            distr = "tnorm",
                            n = concept_k * lambda,
                            par1 = mu,
                            par2 = sigma
                        ),
                        .f = g_max
                    )
                ) |>
                summarise(
                    p_fit = mean(pdf),
                    p_lwr = quantile(pdf, 0.1),
                    p_upr = quantile(pdf, 0.9),
                    .by = size
                )
        }
    )
}

p_ordering <- c(
    "evt" = "Extreme Value Theory",
    "efs" = "Exact Finite Sample",
    "efsmult" = "Exact Finite Sample Multiple Maxima"
)

expmax_vlines <-
    tibble(model_type = names(p_ordering), model_type_full = p_ordering) |>
    mutate(
        max = map2_dbl(
            .x = model_type,
            .y = concept_k,
            .f = get_concept_expmax
        ),
        max20 = map2_dbl(.x = model_type, .y = 20, .f = get_concept_expmax),
        max100 = map2_dbl(.x = model_type, .y = 100, .f = get_concept_expmax)
    ) |>
    mutate(model_type_full = factor(model_type_full, levels = p_ordering))

combined_plot <-
    bind_rows(
        concept_evt_pdf_summary |> mutate(model_type = "evt"),
        concept_efs_pdf_summary |> mutate(model_type = "efs"),
        concept_efsmult_pdf_summary |> mutate(model_type = "efsmult")
    ) |>
    mutate(
        model_type_full = case_when(
            model_type == "evt" ~ "Extreme Value Theory",
            model_type == "efs" ~ "Exact Finite Sample",
            model_type == "efsmult" ~ "Exact Finite Sample Multiple Maxima"
        ),
        light_fill = case_when(
            model_type == "evt" ~ evt_colour,
            model_type == "efs" ~ efs_colour,
            model_type == "efsmult" ~ efsmult_colour
        ),
        dark_fill = case_when(
            model_type == "evt" ~ evt_colour_dark,
            model_type == "efs" ~ efs_colour_dark,
            model_type == "efsmult" ~ efsmult_colour_dark
        ),
        main_lty = case_when(
            model_type == "evt" ~ "solid",
            model_type == "efs" ~ "dashed",
            model_type == "efsmult" ~ "dashed"
        )
    ) |>
    mutate(model_type_full = factor(model_type_full, levels = p_ordering)) |>
    filter(size <= 140) |>
    ggplot(aes(x = size, y = p_fit)) +
    geom_line(aes(y = p), data = concept_popln_pdf) +
    geom_ribbon(
        aes(ymin = p_lwr, ymax = p_upr, fill = light_fill),
        alpha = 0.4
    ) +
    geom_line(
        aes(y = p_fit, col = light_fill, lty = main_lty),
        linewidth = 2
    ) +
    geom_point(
        aes(x = value, y = 0),
        col = "black",
        fill = data_colour_ismax,
        size = 4,
        pch = 21,
        alpha = 0.7,
        data = samples_df |> filter(is_max)
    ) +
    geom_point(
        aes(x = value, y = 0),
        col = "black",
        fill = data_colour_nearmax,
        size = 4,
        pch = 21,
        alpha = 0.7,
        data = samples_df |>
            filter(near_max, !is_max) |>
            mutate(
                model_type_full = factor(
                    "Exact Finite Sample Multiple Maxima",
                    levels = p_ordering
                )
            )
    ) +
    geom_ribbon(
        aes(ymin = p_lwr, ymax = p_upr),
        fill = efs_colour_dark,
        alpha = 0.5,
        data = concept_efs_pdfmax_summary |>
            mutate(
                model_type_full = factor(
                    "Exact Finite Sample",
                    levels = p_ordering
                )
            )
    ) +
    geom_line(
        linewidth = 1.5,
        col = efs_colour_dark,
        data = concept_efs_pdfmax_summary |>
            mutate(
                model_type_full = factor(
                    "Exact Finite Sample",
                    levels = p_ordering
                )
            )
    ) +
    geom_ribbon(
        aes(ymin = p_lwr, ymax = p_upr),
        fill = efsmult_colour_dark,
        alpha = 0.5,
        data = concept_efsmult_pdfmax_summary |>
            mutate(
                model_type_full = factor(
                    "Exact Finite Sample Multiple Maxima",
                    levels = p_ordering
                )
            )
    ) +
    geom_line(
        linewidth = 1.5,
        col = efsmult_colour_dark,
        data = concept_efsmult_pdfmax_summary |>
            mutate(
                model_type_full = factor(
                    "Exact Finite Sample Multiple Maxima",
                    levels = p_ordering
                )
            )
    ) +
    geom_vline(aes(xintercept = max), lty = "dotted", data = expmax_vlines) +
    geom_vline(aes(xintercept = max20), data = expmax_vlines, lty = "dashed") +
    geom_vline(aes(xintercept = max100), data = expmax_vlines, lty = "solid") +

    labs(x = "Body Size (cm)", y = "Density") +
    scale_fill_identity() +
    scale_colour_identity() +
    scale_linetype_identity() +
    facet_wrap(~model_type_full) +
    theme_classic(20) +
    theme(
        legend.position = "none",
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_blank()
    )

# Figure 1 ---------------------------------------------------------------------

# ylim_max <-
#     max(c(
#         layer_scales(plot_fit_summary("evt"))$y$range$range,
#         layer_scales(plot_fit_summary("efs"))$y$range$range,
#         layer_scales(plot_fit_summary("efsmult"))$y$range$range
#     ))

concept_plot <-
    # concept_underlying +
    concept_samples +
    plot_spacer() +
    combined_plot +
    plot_layout(design = "A\nB\nC", heights = c(2, -0.5, 5)) +
    plot_annotation(tag_levels = "A")


ggsave(
    filename = "results/figures/p_concept_MULT.png",
    plot = concept_plot,
    height = 10,
    width = 14
)
