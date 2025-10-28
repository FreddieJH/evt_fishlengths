# Figure 1: Conceptual figure comparing the two methods of estimating maximum length of a fish

library(readr)
library(dplyr)
library(purrr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(arrow)

source("R/01_funcs.R")

if (!exists("estmax_posterior")) {
    estmax_posterior <- read_csv(
        "results/data/estmax_posterior.csv",
        show_col_types = FALSE
    )
}

if (!exists("posterior")) {
    posterior <- read_parquet("results/data/posterior.parquet")
}

if (!exists("scenarios")) {
    scenarios <-
        read_csv("results/data/scenarios.csv", show_col_types = FALSE) |>
        nest(data = topm) |>
        mutate(topm = map(data, \(x) as.numeric(x$topm))) |>
        select(-data) |>
        nest(
            .by = c(k, lambda, dist_name, dist_mean, scenario_id),
            .key = "samples"
        )
}

concept_popln_mean <- 10
concept_k <- 5 # number of samples (i.e., maxima values)
concept_n <- 1000 # sample size per sample

selected_scenario <-
    scenarios |>
    filter(
        k == concept_k,
        lambda == concept_n,
        dist_name == "tnorm",
        dist_mean == concept_popln_mean
    )

set.seed(1)
concept_data <-
    selected_scenario |>
    pull(samples) |>
    pluck(1) |>
    mutate(sample_id = row_number()) |>
    mutate(
        rsample = pmap(
            .l = list("tnorm", concept_popln_mean, n_k),
            .f = ~ {
                pars <- get_dist_pars(..1, ..2)
                rsample <- get(paste0("r", ..1))(..3, pars[1], pars[2])
            }
        )
    ) |>
    unnest(rsample) |>
    mutate(is_max = top1 == rsample) |>
    mutate(near_max = map2_lgl(rsample, topm, ~ .x %in% .y)) |>
    mutate(
        rank = case_when(
            !is_max & !near_max ~ "low",
            !is_max & near_max ~ "high",
            is_max & near_max ~ "highest"
        )
    )

concept_samples_plot <-
    concept_data |>
    ggplot(aes(x = rsample)) +
    geom_dotplot(
        aes(fill = rank),
        col = "grey70",
        method = "histodot",
        # binpositions = "all",
        dotsize = 4,
        binwidth = 0.2,
        data = concept_data |> filter(rank == "low")
    ) +
    geom_dotplot(
        aes(fill = rank),
        method = "histodot",
        # binpositions = "all",
        dotsize = 4,
        binwidth = 0.2,
        data = concept_data |> filter(rank != "low")
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


selected_posterior_evt <-
    posterior |>
    filter(
        scenario_id == unique(selected_scenario$scenario_id),
        model_id == "evt"
    ) |>
    pivot_wider(names_from = par, values_from = value)

selected_posterior_evtg <-
    posterior |>
    filter(
        scenario_id == unique(selected_scenario$scenario_id),
        model_id == "evtg"
    ) |>
    pivot_wider(names_from = par, values_from = value)

selected_posterior_efs <-
    posterior |>
    filter(
        scenario_id == unique(selected_scenario$scenario_id),
        model_id == "efs"
    ) |>
    pivot_wider(names_from = par, values_from = value)

selected_posterior_efsm <-
    posterior |>
    filter(
        scenario_id == unique(selected_scenario$scenario_id),
        model_id == "efs"
    ) |>
    pivot_wider(names_from = par, values_from = value)


gmax <- \(x, lambda, k, mu, sigma) {
    cdf <- \(y) {
        ptnorm(q = y, mean = mu, sd = sigma)
    }
    pdf <- \(y) {
        dtnorm(x = y, mean = mu, sd = sigma)
    }
    gmax <- g(x = x, n = lambda * k, cdf = cdf, pdf = pdf)

    return(gmax)
}

Gmax <- \(x, lambda, k, mu, sigma) {
    cdf <- \(y) {
        ptnorm(q = y, mean = mu, sd = sigma)
    }
    Gmax <- G(x = x, n = lambda * k, cdf = cdf)
    return(Gmax)
}

x_seq <- seq(concept_popln_mean * 0.8, concept_popln_mean * 3, by = 0.1)

pdf_matrix_evt <- vapply(
    seq_along(selected_posterior_evt$loc),
    \(i) {
        dgev(
            x_seq,
            selected_posterior_evt$loc[i],
            selected_posterior_evt$scale[i],
            selected_posterior_evt$shape[i]
        )
    },
    numeric(length(x_seq))
)
cdf_matrix_evt <- vapply(
    seq_along(selected_posterior_evt$loc),
    \(i) {
        pgev(
            x_seq,
            selected_posterior_evt$loc[i],
            selected_posterior_evt$scale[i],
            selected_posterior_evt$shape[i]
        )
    },
    numeric(length(x_seq))
)
pdf_matrix_evtg <- vapply(
    seq_along(selected_posterior_evt$loc),
    \(i) {
        dgumbel(
            x_seq,
            selected_posterior_evt$loc[i],
            selected_posterior_evt$scale[i]
        )
    },
    numeric(length(x_seq))
)
cdf_matrix_evtg <- vapply(
    seq_along(selected_posterior_evt$loc),
    \(i) {
        pgumbel(
            x_seq,
            selected_posterior_evt$loc[i],
            selected_posterior_evt$scale[i]
        )
    },
    numeric(length(x_seq))
)
pdf_matrix_efs <- vapply(
    seq_along(selected_posterior_efs$lambda),
    \(i) {
        gmax(
            x_seq,
            lambda = selected_posterior_efs$lambda[i],
            k = concept_k,
            mu = selected_posterior_efs$mu[i],
            sigma = selected_posterior_efs$sigma[i]
        )
    },
    numeric(length(x_seq))
)
pdf_matrix_efs20 <- vapply(
    seq_along(selected_posterior_efs$lambda),
    \(i) {
        gmax(
            x_seq,
            lambda = selected_posterior_efs$lambda[i],
            k = 20,
            mu = selected_posterior_efs$mu[i],
            sigma = selected_posterior_efs$sigma[i]
        )
    },
    numeric(length(x_seq))
)
underlying_matrix_efs <- vapply(
    seq_along(selected_posterior_efs$lambda),
    \(i) {
        dtnorm(
            x_seq,
            mean = selected_posterior_efs$mu[i],
            sd = selected_posterior_efs$sigma[i]
        )
    },
    numeric(length(x_seq))
)
cdf_matrix_efs <- vapply(
    seq_along(selected_posterior_efs$lambda),
    \(i) {
        Gmax(
            x_seq,
            lambda = selected_posterior_efs$lambda[i],
            k = concept_k,
            mu = selected_posterior_efs$mu[i],
            sigma = selected_posterior_efs$sigma[i]
        )
    },
    numeric(length(x_seq))
)
cdf_matrix_efs20 <- vapply(
    seq_along(selected_posterior_efs$lambda),
    \(i) {
        Gmax(
            x_seq,
            lambda = selected_posterior_efs$lambda[i],
            k = 20,
            mu = selected_posterior_efs$mu[i],
            sigma = selected_posterior_efs$sigma[i]
        )
    },
    numeric(length(x_seq))
)
pdf_matrix_efsm <- vapply(
    seq_along(selected_posterior_efsm$lambda),
    \(i) {
        gmax(
            x_seq,
            lambda = selected_posterior_efsm$lambda[i],
            k = concept_k,
            mu = selected_posterior_efsm$mu[i],
            sigma = selected_posterior_efsm$sigma[i]
        )
    },
    numeric(length(x_seq))
)
pdf_matrix_efsm20 <- vapply(
    seq_along(selected_posterior_efsm$lambda),
    \(i) {
        gmax(
            x_seq,
            lambda = selected_posterior_efsm$lambda[i],
            k = 20,
            mu = selected_posterior_efsm$mu[i],
            sigma = selected_posterior_efsm$sigma[i]
        )
    },
    numeric(length(x_seq))
)
underlying_matrix_efsm <- vapply(
    seq_along(selected_posterior_efsm$lambda),
    \(i) {
        dtnorm(
            x_seq,
            mean = selected_posterior_efsm$mu[i],
            sd = selected_posterior_efsm$sigma[i]
        )
    },
    numeric(length(x_seq))
)
cdf_matrix_efsm <- vapply(
    seq_along(selected_posterior_efsm$lambda),
    \(i) {
        Gmax(
            x_seq,
            lambda = selected_posterior_efsm$lambda[i],
            k = concept_k,
            mu = selected_posterior_efsm$mu[i],
            sigma = selected_posterior_efsm$sigma[i]
        )
    },
    numeric(length(x_seq))
)
cdf_matrix_efsm20 <- vapply(
    seq_along(selected_posterior_efsm$lambda),
    \(i) {
        Gmax(
            x_seq,
            lambda = selected_posterior_efsm$lambda[i],
            k = 20,
            mu = selected_posterior_efsm$mu[i],
            sigma = selected_posterior_efsm$sigma[i]
        )
    },
    numeric(length(x_seq))
)

g_max_evt <- tibble(
    x = x_seq,
    pdf_fit = apply(pdf_matrix_evt, 1, median),
    pdf_lwr = apply(pdf_matrix_evt, 1, quantile, 0.1),
    pdf_upr = apply(pdf_matrix_evt, 1, quantile, 0.9),
    cdf_fit = apply(cdf_matrix_evt, 1, median),
    cdf_lwr = apply(cdf_matrix_evt, 1, quantile, 0.1),
    cdf_upr = apply(cdf_matrix_evt, 1, quantile, 0.9)
)

g_max_evtg <- tibble(
    x = x_seq,
    pdf_fit = apply(pdf_matrix_evtg, 1, median),
    pdf_lwr = apply(pdf_matrix_evtg, 1, quantile, 0.1),
    pdf_upr = apply(pdf_matrix_evtg, 1, quantile, 0.9),
    cdf_fit = apply(cdf_matrix_evtg, 1, median),
    cdf_lwr = apply(cdf_matrix_evtg, 1, quantile, 0.1),
    cdf_upr = apply(cdf_matrix_evtg, 1, quantile, 0.9)
)

g_max_efs <- tibble(
    x = x_seq,
    pdf_fit = apply(pdf_matrix_efs, 1, median),
    pdf_lwr = apply(pdf_matrix_efs, 1, quantile, 0.1),
    pdf_upr = apply(pdf_matrix_efs, 1, quantile, 0.9),
    pdf20_fit = apply(pdf_matrix_efs20, 1, median),
    pdf20_lwr = apply(pdf_matrix_efs20, 1, quantile, 0.1),
    pdf20_upr = apply(pdf_matrix_efs20, 1, quantile, 0.9),
    cdf_fit = apply(cdf_matrix_efs, 1, median),
    cdf_lwr = apply(cdf_matrix_efs, 1, quantile, 0.1),
    cdf_upr = apply(cdf_matrix_efs, 1, quantile, 0.9)
) |>
    left_join(
        tibble(
            x = x_seq,
            underlying_fit = apply(underlying_matrix_efs, 1, median),
            underlying_lwr = apply(underlying_matrix_efs, 1, quantile, 0.1),
            underlying_upr = apply(underlying_matrix_efs, 1, quantile, 0.9),
        )
    )

g_max_efsm <- tibble(
    x = x_seq,
    pdf_fit = apply(pdf_matrix_efsm, 1, median),
    pdf_lwr = apply(pdf_matrix_efsm, 1, quantile, 0.1),
    pdf_upr = apply(pdf_matrix_efsm, 1, quantile, 0.9),
    pdf20_fit = apply(pdf_matrix_efsm20, 1, median),
    pdf20_lwr = apply(pdf_matrix_efsm20, 1, quantile, 0.1),
    pdf20_upr = apply(pdf_matrix_efsm20, 1, quantile, 0.9),
    cdf_fit = apply(cdf_matrix_efsm, 1, median),
    cdf_lwr = apply(cdf_matrix_efsm, 1, quantile, 0.1),
    cdf_upr = apply(cdf_matrix_efsm, 1, quantile, 0.9)
) |>
    left_join(
        tibble(
            x = x_seq,
            underlying_fit = apply(underlying_matrix_efsm, 1, median),
            underlying_lwr = apply(underlying_matrix_efsm, 1, quantile, 0.1),
            underlying_upr = apply(underlying_matrix_efsm, 1, quantile, 0.9),
        )
    )


# g_max_evt <-
#     posterior |>
#     filter(
#         scenario_id == unique(selected_scenario$scenario_id),
#         model_id == "evt"
#     ) |>
#     pivot_wider(names_from = par, values_from = value) |>
#     expand_grid(
#         x = seq(concept_popln_mean * 0.8, concept_popln_mean * 3, by = 0.1)
#     ) |>
#     mutate(
#         pdf = pmap_dbl(
#             .l = list(x = x, loc = loc, scale = scale, shape = shape),
#             .f = dgev
#         ),
#         cdf = pmap_dbl(
#             .l = list(q = x, loc = loc, scale = scale, shape = shape),
#             .f = pgev
#         )
#     ) |>
#     dplyr::summarise(
#         pdf_fit = quantile(pdf, 0.5),
#         pdf_lwr = quantile(pdf, 0.1),
#         pdf_upr = quantile(pdf, 0.9),
#         cdf_fit = quantile(cdf, 0.5),
#         cdf_lwr = quantile(cdf, 0.1),
#         cdf_upr = quantile(cdf, 0.9),
#         .by = x
#     )

# g_max_evtg <-
#     posterior |>
#     filter(
#         scenario_id == unique(selected_scenario$scenario_id),
#         model_id == "evtg"
#     ) |>
#     pivot_wider(names_from = par, values_from = value) |>
#     expand_grid(
#         x = seq(concept_popln_mean * 0.8, concept_popln_mean * 3, by = 0.1)
#     ) |>
#     mutate(
#         pdf = pmap_dbl(
#             .l = list(x = x, loc = loc, scale = scale),
#             .f = dgumbel
#         ),
#         cdf = pmap_dbl(
#             .l = list(q = x, loc = loc, scale = scale),
#             .f = pgumbel
#         )
#     ) |>
#     dplyr::summarise(
#         pdf_fit = quantile(pdf, 0.5),
#         pdf_lwr = quantile(pdf, 0.1),
#         pdf_upr = quantile(pdf, 0.9),
#         cdf_fit = quantile(cdf, 0.5),
#         cdf_lwr = quantile(cdf, 0.1),
#         cdf_upr = quantile(cdf, 0.9),
#         .by = x
#     )

# g_max_efs <-
#     posterior |>
#     filter(
#         scenario_id == unique(selected_scenario$scenario_id),
#         model_id == "efs"
#     ) |>
#     pivot_wider(names_from = par, values_from = value) |>
#     left_join(scenarios |> select(scenario_id, k)) |>
#     expand_grid(
#         x = seq(concept_popln_mean * 0.8, concept_popln_mean * 3, by = 0.1)
#     ) |>
#     mutate(
#         gmax = pmap_dbl(
#             .l = list(x, lambda, k, mu, sigma),
#             .f = \(x, lambda, k, mu, sigma) {
#                 cdf <- \(y) {
#                     ptnorm(q = y, mean = mu, sd = sigma)
#                 }
#                 pdf <- \(y) {
#                     dtnorm(x = y, mean = mu, sd = sigma)
#                 }
#                 gmax <- g(x = x, n = lambda * k, cdf = cdf, pdf = pdf)

#                 return(gmax)
#             }
#         ),
#         gmax20 = pmap_dbl(
#             .l = list(x, lambda, k, mu, sigma),
#             .f = \(x, lambda, k, mu, sigma) {
#                 cdf <- \(y) {
#                     ptnorm(q = y, mean = mu, sd = sigma)
#                 }
#                 pdf <- \(y) {
#                     dtnorm(x = y, mean = mu, sd = sigma)
#                 }
#                 gmax <- g(x = x, n = lambda * 20, cdf = cdf, pdf = pdf)

#                 return(gmax)
#             }
#         ),
#         pdf_under = pmap_dbl(
#             .l = list(
#                 x = x,
#                 mean = mu,
#                 sd = sigma
#             ),
#             .f = dtnorm
#         ),
#         Gmax = pmap_dbl(
#             .l = list(x, lambda, k, mu, sigma),
#             .f = \(x, lambda, k, mu, sigma) {
#                 cdf <- \(y) {
#                     ptnorm(q = y, mean = mu, sd = sigma)
#                 }
#                 Gmax <- G(x = x, n = lambda * k, cdf = cdf)
#                 return(Gmax)
#             }
#         ),
#         Gmax20 = pmap_dbl(
#             .l = list(x, lambda, k, mu, sigma),
#             .f = \(x, lambda, k, mu, sigma) {
#                 cdf <- \(y) {
#                     ptnorm(q = y, mean = mu, sd = sigma)
#                 }
#                 Gmax <- G(x = x, n = lambda * 20, cdf = cdf)
#                 return(Gmax)
#             }
#         )
#     ) |>
#     dplyr::summarise(
#         pdf_fit = quantile(gmax, 0.5),
#         pdf_lwr = quantile(gmax, 0.1),
#         pdf_upr = quantile(gmax, 0.9),
#         pdf20_fit = quantile(gmax20, 0.5),
#         pdf20_lwr = quantile(gmax20, 0.1),
#         pdf20_upr = quantile(gmax20, 0.9),
#         pdf_fit_under = quantile(pdf_under, 0.5),
#         pdf_lwr_under = quantile(pdf_under, 0.1),
#         pdf_upr_under = quantile(pdf_under, 0.9),
#         cdf_fit = quantile(Gmax, 0.5),
#         cdf_lwr = quantile(Gmax, 0.1),
#         cdf_upr = quantile(Gmax, 0.9),
#         cdf20_fit = quantile(Gmax20, 0.5),
#         cdf20_lwr = quantile(Gmax20, 0.1),
#         cdf20_upr = quantile(Gmax20, 0.9),
#         .by = x
#     )

# g_max_efsm <-
#     posterior |>
#     filter(
#         scenario_id == unique(selected_scenario$scenario_id),
#         model_id == "efsm"
#     ) |>
#     pivot_wider(names_from = par, values_from = value) |>
#     left_join(scenarios |> select(scenario_id, k)) |>
#     expand_grid(
#         x = seq(concept_popln_mean * 0.8, concept_popln_mean * 3, by = 0.1)
#     ) |>
#     mutate(
#         pdf = pmap_dbl(
#             .l = list(x, lambda, k, mu, sigma),
#             .f = \(x, lambda, k, mu, sigma) {
#                 cdf <- \(y) {
#                     ptnorm(q = y, mean = mu, sd = sigma)
#                 }
#                 pdf <- \(y) {
#                     dtnorm(x = y, mean = mu, sd = sigma)
#                 }
#                 gmax <- g(x = x, n = lambda * k, cdf = cdf, pdf = pdf)

#                 return(gmax)
#             }
#         ),
#         pdf20 = pmap_dbl(
#             .l = list(x, lambda, k, mu, sigma),
#             .f = \(x, lambda, k, mu, sigma) {
#                 cdf <- \(y) {
#                     ptnorm(q = y, mean = mu, sd = sigma)
#                 }
#                 pdf <- \(y) {
#                     dtnorm(x = y, mean = mu, sd = sigma)
#                 }
#                 gmax <- g(x = x, n = lambda * 20, cdf = cdf, pdf = pdf)

#                 return(gmax)
#             }
#         ),
#         pdf_under = pmap_dbl(
#             .l = list(
#                 x = x,
#                 mean = mu,
#                 sd = sigma
#             ),
#             .f = dtnorm
#         ),
#         cdf = pmap_dbl(
#             .l = list(x, lambda, k, mu, sigma),
#             .f = \(x, lambda, k, mu, sigma) {
#                 cdf <- \(y) {
#                     ptnorm(q = y, mean = mu, sd = sigma)
#                 }
#                 Gmax <- G(x = x, n = lambda * k, cdf = cdf)
#                 return(Gmax)
#             }
#         ),
#         cdf20 = pmap_dbl(
#             .l = list(x, lambda, k, mu, sigma),
#             .f = \(x, lambda, k, mu, sigma) {
#                 cdf <- \(y) {
#                     ptnorm(q = y, mean = mu, sd = sigma)
#                 }
#                 Gmax <- G(x = x, n = lambda * 20, cdf = cdf)
#                 return(Gmax)
#             }
#         )
#     ) |>
#     dplyr::summarise(
#         pdf_fit = quantile(pdf, 0.5),
#         pdf_lwr = quantile(pdf, 0.1),
#         pdf_upr = quantile(pdf, 0.9),
#         pdf20_fit = quantile(pdf20, 0.5),
#         pdf20_lwr = quantile(pdf20, 0.1),
#         pdf20_upr = quantile(pdf20, 0.9),
#         pdf_fit_under = quantile(pdf_under, 0.5),
#         pdf_lwr_under = quantile(pdf_under, 0.1),
#         pdf_upr_under = quantile(pdf_under, 0.9),
#         cdf_fit = quantile(cdf, 0.5),
#         cdf_lwr = quantile(cdf, 0.1),
#         cdf_upr = quantile(cdf, 0.9),
#         cdf20_fit = quantile(cdf20, 0.5),
#         cdf20_lwr = quantile(cdf20, 0.1),
#         cdf20_upr = quantile(cdf20, 0.9),
#         .by = x
#     )

p_ordering <- c(
    "evt" = "Extreme Value Theory",
    "evtg" = "Extreme Value Theory (Gumbel)",
    "efs" = "Exact Finite Sample",
    "efsm" = "Exact Finite Sample Multiple Maxima"
)

mods <- c(
    "evt" = "gev",
    "evtg" = "gumbel",
    "efs" = "tnorm",
    "efsm" = "tnorm"
)


concept_plot_data <-
    bind_rows(
        g_max_evt |> mutate(model_id = "evt"),
        g_max_evtg |> mutate(model_id = "evtg"),
        g_max_efs |> mutate(model_id = "efs"),
        g_max_efsm |> mutate(model_id = "efsm")
    ) |>
    mutate(
        model_id_full = case_when(
            model_id == "evt" ~ "Extreme Value Theory",
            model_id == "evtg" ~ "Extreme Value Theory (Gumbel)",
            model_id == "efs" ~ "Exact Finite Sample",
            model_id == "efsm" ~ "Exact Finite Sample Multiple Maxima"
        )
    ) |>
    mutate(model_id_full = factor(model_id_full, levels = p_ordering))


est_vals <-
    estmax_posterior |>
    left_join(scenarios) |>
    filter(
        k == concept_k,
        lambda == concept_n,
        dist_name == "tnorm",
        dist_mean == concept_popln_mean
    ) |>
    select(contains("est_max"), model_id) |>
    # pivot_longer(-model_id) |>
    mutate(
        model_id_full = case_when(
            model_id == "evt" ~ "Extreme Value Theory",
            model_id == "evtg" ~ "Extreme Value Theory (Gumbel)",
            model_id == "efs" ~ "Exact Finite Sample",
            model_id == "efsm" ~ "Exact Finite Sample Multiple Maxima"
        )
    ) |>
    mutate(model_id_full = factor(model_id_full, levels = p_ordering))


pdf_plot <-
    concept_plot_data |>
    ggplot(aes(x = x, y = pdf_fit)) +
    geom_line(
        data = tibble(
            x = seq(concept_popln_mean * 0.8, concept_popln_mean * 3, by = 0.1),
            pdf_fit = dtnorm(x, concept_popln_mean, concept_popln_mean * 0.34)
        )
    ) +
    geom_ribbon(
        aes(ymin = underlying_lwr, ymax = underlying_upr, fill = model_id),
        alpha = 0.3,
    ) +
    geom_line(
        aes(y = underlying_fit, col = model_id),
        linewidth = 2,
        lty = "31",
        alpha = 0.5
    ) +
    geom_ribbon(
        aes(ymin = pdf_lwr, ymax = pdf_upr, fill = model_id),
        alpha = 0.5
    ) +
    geom_line(aes(y = pdf_fit, col = model_id), linewidth = 2) +
    geom_ribbon(
        aes(ymin = pdf20_lwr, ymax = pdf20_upr, fill = model_id),
        alpha = 0.5
    ) +
    geom_line(aes(y = pdf20_fit, col = model_id), linewidth = 2) +

    # geom_ribbon(aes(ymin = pdf_lwr, ymax = pdf_upr, fill = model_id), alpha = 0.5) +
    # geom_line(aes(y = pdf_fit, col = model_id), linewidth = 2) +
    # geom_vline(aes(xintercept = max), lty = "dotted", data = expmax_vlines) +
    # geom_vline(aes(xintercept = max20), data = expmax_vlines, lty = "dashed") +
    # geom_vline(aes(xintercept = max100), data = expmax_vlines, lty = "solid") +
    geom_point(
        aes(x = topm, y = 0),
        col = "black",
        fill = data_colour_nearmax,
        size = 4,
        pch = 21,
        alpha = 0.7,
        data = concept_data |>
            filter(is_max) |>
            mutate(pos = rank(top1) / max(rank(top1) + 1)) |>
            unnest(topm) |>
            filter(topm != top1) |>
            mutate(
                model_id_full = factor(
                    "Exact Finite Sample Multiple Maxima",
                    levels = p_ordering
                )
            )
    ) +
    geom_point(
        aes(x = rsample, y = 0),
        col = "black",
        fill = data_colour_ismax,
        size = 4,
        pch = 21,
        alpha = 0.7,
        data = concept_data |>
            filter(is_max) |>
            mutate(pos = rank(top1) / max(rank(top1) + 1))
    ) +
    geom_errorbarh(
        aes(y = -0.1, xmin = est_max20_lwr, xmax = est_max20_upr),
        data = est_vals,
        inherit.aes = FALSE,
        height = 0.01
    ) +
    geom_point(
        aes(y = -0.1, x = est_max20_fit),
        data = est_vals,
        inherit.aes = FALSE,
        fill = "white",
        pch = 21,
        size = 4
    ) +
    geom_errorbar(
        aes(y = -0.05, xmin = est_max_lwr, xmax = est_max_upr),
        data = est_vals,
        inherit.aes = FALSE,
        height = 0.01,
        orientation = "y"
    ) +
    geom_point(
        aes(y = -0.05, x = est_max_fit),
        data = est_vals,
        inherit.aes = FALSE,
        fill = "white",
        pch = 21,
        size = 4
    ) +
    labs(x = "Body Size (cm)", y = "Density") +
    scale_fill_manual(
        values = c("evt" = evt_colour, "efs" = efs_colour, "efsm" = efsm_colour)
    ) +
    scale_colour_manual(
        values = c("evt" = evt_colour, "efs" = efs_colour, "efsm" = efsm_colour)
    ) +
    facet_wrap(~model_id_full) +
    theme_classic(20) +
    theme(
        legend.position = "none" #,
        # axis.ticks.y = element_blank(),
        # axis.text.y = element_blank(),
        # axis.title.y = element_blank()
    )
pdf_plot

# cdf_plot <-
#     concept_plot_data |>
#     ggplot(aes(x = x, y = pdf_fit)) +
#     geom_ribbon(
#         aes(ymin = cdf_lwr, ymax = cdf_upr, fill = model_id),
#         alpha = 0.5
#     ) +
#     geom_ribbon(
#         aes(ymin = cdf20_lwr, ymax = cdf20_upr, fill = model_id),
#         alpha = 0.8
#     ) +
#     geom_line(aes(y = cdf20_fit, col = model_id), linewidth = 2) +
#     geom_line(aes(y = cdf_fit, col = model_id), linewidth = 2) +
#     scale_alpha_manual(values = c("TRUE" = 0.2, "FALSE" = 1), guide = "none") +
#     geom_point(
#         aes(
#             x = topm,
#             y = pos,
#             alpha = model_id %in% c("efs", "efsm") & pos < 0.8
#         ),
#         col = "black",
#         fill = data_colour_nearmax,
#         size = 4,
#         pch = 21,
#         data = concept_data |>
#             filter(is_max) |>
#             mutate(pos = rank(top1) / max(rank(top1) + 1)) |>
#             unnest(topm) |>
#             filter(topm != top1) |>
#             mutate(
#                 model_id_full = factor(
#                     "Exact Finite Sample Multiple Maxima",
#                     levels = p_ordering
#                 )
#             ) |>
#             expand_grid(model_id = c("evt", "evtg", "efs", "efsm")) |>
#             mutate(
#                 model_id_full = case_when(
#                     model_id == "evt" ~ "Extreme Value Theory",
#                     model_id == "evtg" ~ "Extreme Value Theory (Gumbel)",
#                     model_id == "efs" ~ "Exact Finite Sample",
#                     model_id == "efsm" ~ "Exact Finite Sample Multiple Maxima"
#                 )
#             )
#     ) +
#     geom_point(
#         aes(
#             x = rsample,
#             y = pos,
#             alpha = model_id %in% c("efs", "efsm") & pos < 0.8
#         ),
#         col = "black",
#         fill = data_colour_ismax,
#         size = 4,
#         pch = 21,
#         data = concept_data |>
#             filter(is_max) |>
#             mutate(pos = rank(top1) / max(rank(top1) + 1)) |>
#             expand_grid(model_id = c("evt", "evtg", "efs", "efsm")) |>
#             mutate(
#                 model_id_full = case_when(
#                     model_id == "evt" ~ "Extreme Value Theory",
#                     model_id == "evtg" ~ "Extreme Value Theory (Gumbel)",
#                     model_id == "efs" ~ "Exact Finite Sample",
#                     model_id == "efsm" ~ "Exact Finite Sample Multiple Maxima"
#                 )
#             )
#     ) +
#     # horizontal at 0.95
#     geom_segment(
#         aes(group = model_id_full, xend = est_max20_upr),
#         y = 0.95,
#         yend = 0.95,
#         x = -Inf,
#         col = "grey70",
#         lty = 2,
#         data = est_vals |> filter(model_id %in% c("evt", "evtg"))
#     ) +
#     geom_errorbar(
#         aes(y = 0.95, xmin = est_max20_lwr, xmax = est_max20_upr),
#         data = est_vals |> filter(model_id %in% c("evt", "evtg")),
#         inherit.aes = FALSE,
#         width = 0.01,
#         orientation = "y"
#     ) +
#     geom_point(
#         aes(y = 0.95, x = est_max20_fit),
#         data = est_vals |> filter(model_id %in% c("evt", "evtg")),
#         inherit.aes = FALSE,
#         fill = "white",
#         pch = 21,
#         size = 4
#     ) +
#     geom_segment(
#         aes(group = model_id_full, xend = est_max20_lwr, x = est_max20_lwr),
#         y = -Inf,
#         yend = 0.95,
#         col = "grey60",
#         lty = 2,
#         data = est_vals |> filter(model_id %in% c("evt", "evtg"))
#     ) +
#     geom_segment(
#         aes(group = model_id_full, x = est_max20_fit, xend = est_max20_fit),
#         y = -Inf,
#         yend = 0.95,
#         col = "grey60",
#         lty = 2,
#         data = est_vals |> filter(model_id %in% c("evt", "evtg"))
#     ) +
#     geom_segment(
#         aes(group = model_id_full, x = est_max20_upr, xend = est_max20_upr),
#         y = -Inf,
#         yend = 0.95,
#         col = "grey60",
#         lty = 2,
#         data = est_vals |> filter(model_id %in% c("evt", "evtg"))
#     ) +
#     geom_vline(
#         aes(
#             xintercept = scenarios_truemax |>
#                 filter(
#                     k == concept_k,
#                     lambda == concept_n,
#                     dist_mean == concept_popln_mean,
#                     dist_name == "tnorm"
#                 ) |>
#                 pull(true_max20)
#         ),
#     ) +
#     labs(x = "Body Size (cm)", y = "Cumulative density") +
#     scale_fill_manual(
#         values = c("evt" = evt_colour, "efs" = efs_colour, "efsm" = efsm_colour)
#     ) +
#     scale_colour_manual(
#         values = c("evt" = evt_colour, "efs" = efs_colour, "efsm" = efsm_colour)
#     ) +
#     facet_wrap(~model_id_full) +
#     theme_classic(20) +
#     theme(
#         legend.position = "none" #,
#         # axis.ticks.y = element_blank(),
#         # axis.text.y = element_blank(),
#         # axis.title.y = element_blank()
#     )
# cdf_plot
# Figure 1 ---------------------------------------------------------------------

concept_plot <-
    concept_samples_plot +
    plot_spacer() +
    pdf_plot +
    # cdf_plot +
    plot_layout(design = "A\nB\nC", heights = c(2, -0.5, 5)) +
    plot_annotation(tag_levels = "A")

ggsave(
    filename = "results/figures/manuscript_figures/concept.png",
    plot = concept_plot,
    height = 16,
    width = 16
)

# selected_scenario |>
#     pull(samples) |>
#     pluck(1) |>
#     pull(top1) |>
#     sort() |>
#     signif(3)
# selected_scenario |>
#     pull(samples) |>
#     pluck(1) |>
#     pull(top1) |>
#     sort() |>
#     tail(3) |>
#     median()
# mode_f(\(x) {
#     g(x, 1000, cdf = \(x) pnorm(x, 10, 3.4), pdf = \(x) dnorm(x, 10, 3.4))
# })
# mode_f(\(x) {
#     g(x, 5 * 1000, cdf = \(x) pnorm(x, 10, 3.4), pdf = \(x) dnorm(x, 10, 3.4))
# })
# mode_f(\(x) {
#     g(x, 20 * 1000, cdf = \(x) pnorm(x, 10, 3.4), pdf = \(x) dnorm(x, 10, 3.4))
# })

# estmax_posterior |> filter(scenario_id == unique(selected_scenario$scenario_id))
# posterior |>
#     filter(scenario_id == unique(selected_scenario$scenario_id)) |>
#     summarise(
#         median = quantile(value, 0.5),
#         lwr = quantile(value, 0.1),
#         upr = quantile(value, 0.9),
#         .by = c(model_id, par)
#     )
