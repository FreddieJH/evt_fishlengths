# Figure 1: Conceptual figure comparing the two methods of estimating maximum length of a fish

# Script requirements  --------------------------------------------------------------

source("code/00a_truncnorm_functions.R")
source("code/00b_model_functions.R")
source("code/00c_expectedmax_functions.R")
source("code/00d_colours.R")
library(tidyverse)
library(patchwork)
library(evd)

# Script parameters  --------------------------------------------------------------

concept_popln_mean <- 50
concept_popln_sd <- 16
concept_popln_truncation <- 0 # positive body lengths only
concept_k <- 6 # number of samples (i.e., maxima values)
concept_n <- 50 # sample size per sample
save_figures <- TRUE


# install.packages("rcartocolor")
# rcartocolor::carto_pal(8, "Antique")
# rcartocolor::display_carto_pal(7, "Antique")
# scale_color_brewer(palette = "Antique") # use this one? need a similar dark and light version for each colour

# Underlying population  -------------------------------------------------

# Random sampling of truncated normal distribution
set.seed(1)
concept_popln <- rtnorm(
    n = 10000,
    mean = concept_popln_mean,
    sd = concept_popln_sd,
    a = concept_popln_truncation
)

# PDF of the underlying population
concept_popln_pdf <-
    tibble(size = 0:150) |>
    mutate(
        p = dtnorm(
            size,
            mean = concept_popln_mean,
            sd = concept_popln_sd,
            a = concept_popln_truncation
        )
    )

# get k maxima values from the underlying distribution
concept_maxima <- replicate(
    concept_k,
    max(sample(x = concept_popln, size = concept_n))
)

# concept_underlying <-
#     ggplot(concept_popln_pdf) +
#     geom_line(aes(x = size, y = p), linewidth = 1) +
#     labs(x = "Body Size (cm)", y = "Density") +
#     theme_minimal(20) +
#     theme(legend.position = "none", axis.ticks.y = element_blank(), axis.text.y = element_blank(), axis.title.y = element_blank())

# Sampling method --------------------------------------------------------

samples_df <- map_dfr(1:concept_k, function(i) {
    # Generate random sample sizes between 30 and 200
    nn <- sample(10:100, 1)

    # Generate truncated normal data
    sample_data <- rtnorm(
        nn,
        a = 0,
        mean = concept_popln_mean,
        sd = concept_popln_sd
    )

    threshold <- sample(1:4, 1)

    # Create and return data frame for this sample
    tibble(
        sample_id = i,
        value = sample_data,
        is_max = sample_data == max(sample_data),
        near_max = rank(-value, ties.method = "max") <= threshold,
        n = nn, # Store sample size for reference
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

# Model fitting + plotting
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
                p_lwr = quantile(pdf, 0.025),
                p_upr = quantile(pdf, 0.975),
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
                            n = lambda,
                            par1 = mu,
                            par2 = sigma
                        ),
                        .f = g_max
                    )
                ) |>
                summarise(
                    p_fit = mean(pdf),
                    p_lwr = quantile(pdf, 0.025),
                    p_upr = quantile(pdf, 0.975),
                    .by = size
                )
        }
    )

    #  plotting

    model_labels <- c(
        evt = "Extreme Value Theory",
        efs = "Exact Finite Sample",
        efsmult = "Exact Finite Sample Multiple Maxima"
    )

    plot <-
        get(paste0("concept_", model_type, "_pdf_summary")) |>
        ggplot() +
        geom_line(aes(size, p), data = concept_popln_pdf) +
        geom_ribbon(
            aes(x = size, ymin = p_lwr, ymax = p_upr),
            fill = get(paste0(model_type, "_colour")),
            alpha = 0.3
        ) +
        geom_line(
            aes(x = size, y = p_fit),
            linewidth = 2,
            lty = ifelse(model_type == "evt", "solid", "dashed"),
            col = get(paste0(model_type, "_colour"))
        ) +
        {
            if (model_type != "evt") {
                list(
                    geom_ribbon(
                        aes(x = size, ymin = p_lwr, ymax = p_upr),
                        fill = get(paste0(model_type, "_colour_dark")),
                        alpha = 0.3,
                        data = get(paste0(
                            "concept_",
                            model_type,
                            "_pdfmax_summary"
                        ))
                    ),
                    geom_line(
                        aes(x = size, y = p_fit),
                        linewidth = 1,
                        lty = "solid",
                        col = get(paste0(model_type, "_colour_dark")),
                        data = get(paste0(
                            "concept_",
                            model_type,
                            "_pdfmax_summary"
                        ))
                    )
                )
            }
        } +
        geom_point(
            aes(x = value, y = 0),
            col = "black",
            fill = data_colour_ismax,
            size = 4,
            pch = 21,
            alpha = 0.7,
            data = samples_df |> filter(is_max)
        ) +
        {
            if (model_type == "efsmult") {
                geom_point(
                    aes(x = value, y = 0),
                    col = "black",
                    fill = data_colour_nearmax,
                    size = 4,
                    pch = 21,
                    alpha = 0.7,
                    data = samples_df |> filter(near_max, !is_max)
                )
            }
        } +
        geom_vline(
            xintercept = get_concept_expmax(model_type, concept_k),
            col = "black",
            lty = "dotted"
        ) +
        geom_vline(
            xintercept = get_concept_expmax(model_type, k = 20),
            col = "black",
            lty = "dashed"
        ) +
        geom_vline(
            xintercept = get_concept_expmax(model_type, k = 100),
            col = "black",
            lty = "solid"
        ) +
        labs(x = "Body Size (cm)", y = "Density") +
        theme_classic(20) +
        facet_grid(
            cols = vars(1),
            labeller = labeller(.cols = function(x) model_labels[model_type])
        ) +
        theme(
            legend.position = "none",
            axis.ticks.y = element_blank(),
            axis.text.y = element_blank(),
            axis.title.y = element_blank()
        )

    assign(
        x = paste0("concept_", model_type, "_plot"),
        value = plot
    )
}

# Figure 1 ---------------------------------------------------------------------

ylim_max <-
    max(c(
        layer_scales(concept_evt_plot)$y$range$range,
        layer_scales(concept_efs_plot)$y$range$range,
        layer_scales(concept_efsmult_plot)$y$range$range
    ))

concept_plot <-
    # concept_underlying +
    concept_samples +
    concept_evt_plot +
    concept_efs_plot +
    concept_efsmult_plot +
    plot_layout(design = "AAA\nBCD", heights = c(2, 5)) +
    plot_annotation(tag_levels = "A") &
    ylim(0, ylim_max)


ggsave(
    filename = "results/figures/p_concept_MULT.png",
    plot = concept_plot,
    height = 10,
    width = 14
)
