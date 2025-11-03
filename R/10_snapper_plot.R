# make underlying dist fainter in the figure
# test sensitivity figure with extremely high lambdas (way above prior)

library(dplyr)
library(readr)
library(ggplot2)
library(tidyr)
library(purrr)
library(magick)
library(scales)
library(patchwork)

source("R/01_funcs.R")

kg2cm <- function(w, a = 0.04478, b = 2.673) ((w * 1000) / a)^(1 / b)

snapper_maxima <-
  read_csv("data/snapper_maxima.csv") |>
  mutate(
    max = case_when(
      !is.na(length_cm) ~ length_cm,
      is.na(length_cm) ~ kg2cm(weight_kg, a = 0.01, b = 3)
    )
  )
snapper_maxima |> pull(max) |> signif(digits = 4)

# maxima is either a vector of maxmia (for EVT and EFS) list(x1, x2, x3) or c(x1, x2, x3)
# or it is a list of vectors (for EFSMM) list(c(x1, x2, x3), x(y1, y2), c(z1, z2, z3, z4))
fit_mod <- function(maxima, model_type) {
  maxima_list <- as.list(maxima)

  mod_dat <-
    list(
      x = maxima_list |> unlist(),
      n_obs = length(maxima_list |> unlist()),
      n_per_sample = lapply(maxima_list, length) |> unlist(),
      start_idx = cumsum(lapply(maxima_list, length) |> unlist()),
      k = length(maxima_list)
    )

  init_func <- function(type, maxima_median) {
    if (type == "evt") {
      return_func <- function(chain_id) {
        list(
          loc = maxima_median,
          scale = 10,
          shape = 0
        )
      }
    } else if (type == "evtg") {
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

  mod <- switch(
    model_type,
    "evt" = cmdstanr::cmdstan_model("models/evt.stan"),
    "evtg" = cmdstanr::cmdstan_model("models/evt_gumbel.stan"),
    "efs" = cmdstanr::cmdstan_model("models/efs.stan"),
    "efsmult" = cmdstanr::cmdstan_model("models/efs.stan")
  )

  fit <- mod$sample(
    data = mod_dat,
    chains = 4,
    init = init_func(model_type, median(maxima)),
    iter_warmup = 9000,
    iter_sampling = 1000
  )
  return(fit)
}

fit_evt <- fit_mod(snapper_maxima$max, "evt")
fit_evtg <- fit_mod(snapper_maxima$max, "evtg")
fit_efs <- fit_mod(snapper_maxima$max, "efs")

get_posterior <- function(fit) {
  posterior <-
    fit |>
    posterior::as_draws_df() %>%
    dplyr::as_tibble()
  return(posterior)
}

traceplot <- function(fit) {
  posterior <- get_posterior(fit)
  plot <-
    posterior %>%
    pivot_longer(cols = -c(.chain, .iteration, .draw)) |>
    ggplot(aes(
      .iteration,
      value,
      col = as.factor(.chain),
      .group = .chain
    )) +
    geom_path(alpha = 0.4) +
    facet_wrap(
      ~name,
      scales = "free",
      ncol = ceiling(sqrt(ncol(posterior) - 3))
    ) +
    theme_classic(20) +
    theme(legend.position = "none")

  return(plot)
}

est_max_posterior <- function(fit, ci = 0.8, k) {
  is_evt <- sum(c("loc", "scale", "shape") %in% colnames(get_posterior(fit))) ==
    3
  is_evtg <- sum(c("loc", "scale") %in% colnames(get_posterior(fit))) == 2
  get_posterior(fit) |>
    mutate(
      estmax = if (is_evt) {
        pmap_dbl(
          .l = list(
            p = 1 - (1 / k),
            loc = loc,
            scale = scale,
            shape = shape
          ),
          .f = qgev
        )
      } else if (is_evtg) {
        pmap_dbl(
          .l = list(
            p = 1 - (1 / k),
            loc = loc,
            scale = scale
          ),
          .f = qgumbel
        )
      } else {
        pmap_dbl(
          .l = list(mu, sigma, lambda, k),
          .f = \(mu, sigma, lambda, k) {
            cdf <- \(x) {
              ptnorm(q = x, mean = mu, sd = sigma)
            }
            pdf <- \(x) {
              dtnorm(x = x, mean = mu, sd = sigma)
            }
            gmax <- \(x) {
              g(x = x, n = lambda * k, cdf = cdf, pdf = pdf)
            }
            mode_f(gmax)
          }
        )
      }
    )
}

est_max <- function(fit, ci = 0.8, k) {
  est_max_posterior(fit = fit, ci = ci, k = k) |>
    summarise(
      max_fit = mean(estmax),
      max_lwr = quantile(estmax, (1 - ci) / 2),
      max_upr = quantile(estmax, 1 - ((1 - ci) / 2)),
    ) %>%
    as.numeric()
}


get_pdf <- function(
  fit,
  k = length(snapper_maxima$max),
  xmin = 0,
  xmax = 300,
  xstep = 1,
  ci = 0.8
) {
  is_evt <- sum(c("loc", "scale", "shape") %in% colnames(get_posterior(fit))) ==
    3
  is_evtg <- sum(c("loc", "scale") %in% colnames(get_posterior(fit))) == 2
  get_posterior(fit) |>
    expand_grid(size = seq(xmin, xmax, xstep)) |>
    mutate(
      pdf = if (is_evt) {
        pmap_dbl(
          .l = list(
            x = size,
            loc = loc,
            scale = scale,
            shape = shape
          ),
          .f = dgev
        )
      } else if (is_evtg) {
        pmap_dbl(
          .l = list(
            x = size,
            loc = loc,
            scale = scale
          ),
          .f = dgumbel
        )
      } else {
        pmap_dbl(
          .l = list(mu, sigma, lambda, k, size),
          .f = \(mu, sigma, lambda, k, size) {
            cdf <- \(x) {
              ptnorm(q = x, mean = mu, sd = sigma)
            }
            pdf <- \(x) {
              dtnorm(x = x, mean = mu, sd = sigma)
            }
            gmax <- \(x) {
              g(x = x, n = lambda * k, cdf = cdf, pdf = pdf)
            }
            gmax(size)
          }
        )
      }
    ) |>
    summarise(
      pdf_fit = quantile(pdf, 0.5),
      pdf_lwr = quantile(pdf, (1 - ci) / 2),
      pdf_upr = quantile(pdf, 1 - ((1 - ci) / 2)),
      .by = size
    )
}


get_underlying <- function(fit, xmin = 0, xmax = 300, xstep = 1, ci = 0.8) {
  get_posterior(fit) |>
    expand_grid(size = seq(xmin, xmax, xstep)) |>
    mutate(
      pdf = pmap_dbl(
        .l = list(
          x = size,
          mean = mu,
          sd = sigma
        ),
        .f = dtnorm
      )
    ) |>
    summarise(
      pdf_fit = quantile(pdf, 0.5),
      pdf_lwr = quantile(pdf, (1 - ci) / 2),
      pdf_upr = quantile(pdf, 1 - ((1 - ci) / 2)),
      .by = size
    )
}

get_cdf <- function(fit, xmin = 0, xmax = 300, xstep = 1, ci = 0.8) {
  is_evt <- sum(c("loc", "scale", "shape") %in% colnames(get_posterior(fit))) ==
    3
  is_evtg <- sum(c("loc", "scale") %in% colnames(get_posterior(fit))) == 2
  get_posterior(fit) |>
    expand_grid(size = seq(xmin, xmax, xstep)) |>
    mutate(
      cdf = if (is_evt) {
        pmap_dbl(
          .l = list(
            q = size,
            loc = loc,
            scale = scale,
            shape = shape
          ),
          .f = pgev
        )
      } else if (is_evtg) {
        pmap_dbl(
          .l = list(
            q = size,
            loc = loc,
            scale = scale
          ),
          .f = pgumbel
        )
      } else {
        pmap_dbl(
          .l = list(x, lambda, k, mu, sigma),
          .f = \(x, lambda, k, mu, sigma) {
            cdf <- \(y) {
              ptnorm(q = y, mean = mu, sd = sigma)
            }
            Gmax <- G(x = x, n = lambda * k, cdf = cdf)
            return(Gmax)
          }
        )
      }
    ) |>
    summarise(
      cdf_fit = quantile(cdf, 0.5),
      cdf_lwr = quantile(cdf, (1 - ci) / 2),
      cdf_upr = quantile(cdf, 1 - ((1 - ci) / 2)),
      .by = size
    )
}

# get_posterior(fit1)
# get_posterior(fit2)
# traceplot(fit1)
# est_max(fit1, k = 20)
# est_max(fit2, k = 20)
# get_pdf(fit1)
# get_pdf(fit2)

# get_posterior(fit2) |>
#   summarise(
#     mu_fit = quantile(mu, 0.5),
#     mu_lwr = quantile(mu, 0.1),
#     mu_upr = quantile(mu, 0.9),
#     sigma_fit = quantile(sigma, 0.5),
#     sigma_lwr = quantile(sigma, 0.1),
#     sigma_upr = quantile(sigma, 0.9),
#     lambda_fit = quantile(lambda, 0.5),
#     lambda_lwr = quantile(lambda, 0.1),
#     lambda_upr = quantile(lambda, 0.9)
#   )

# get_pdf(fit1) |>
#   ggplot(aes(x = size, y = pdf_fit)) +
#   geom_ribbon(
#     aes(ymin = pdf_lwr, ymax = pdf_upr),
#     alpha = 0.8,
#     fill = "grey50"
#   ) +
#   geom_line() +
#   theme_classic(20) +
#   labs(x = "Body size", y = "Probability density")

# get_cdf(fit1) |>
#   ggplot(aes(x = size, y = cdf_fit)) +
#   geom_ribbon(
#     aes(ymin = cdf_lwr, ymax = cdf_upr),
#     alpha = 0.8,
#     fill = "grey50"
#   ) +
#   geom_line() +
#   theme_classic(20) +
#   labs(x = "Body size", y = "Cumulative density")

evt_pdf <- get_pdf(fit_evt)
evtg_pdf <- get_pdf(fit_evtg)
efs_pdf <- get_pdf(fit_efs)
efs_underlying <- get_underlying(fit_efs)

p_snapper <-
  evt_pdf |>
  ggplot(aes(size, pdf_fit)) +
  geom_ribbon(
    aes(ymin = pdf_lwr, ymax = pdf_upr),
    alpha = 0.2,
    fill = "#A23B72",
    data = efs_underlying
  ) +
  geom_line(
    aes(y = pdf_fit),
    col = "#A23B72",
    lty = 2,
    data = efs_underlying,
    alpha = 0.3
  ) +
  geom_ribbon(
    aes(ymin = pdf_lwr, ymax = pdf_upr),
    alpha = 0.5,
    fill = "#7DB3D3"
  ) +
  geom_line(col = "#2E86AB", linewidth = 2) +
  geom_ribbon(
    aes(ymin = pdf_lwr, ymax = pdf_upr),
    alpha = 0.5,
    data = efs_pdf,
    fill = "#C77BA0"
  ) +
  geom_line(col = "#C77BA0", data = efs_pdf, linewidth = 2) +
  geom_rug(
    aes(x = x),
    data = tibble(x = snapper_maxima$max),
    inherit.aes = FALSE
  ) +
  scale_x_continuous(
    label = scales::label_number(suffix = "cm"),
    limits = c(0, 200)
  ) +
  labs(x = "Body size", y = "Probability density") +
  theme_classic(20)


p2_data <-
  tibble(k = c(length(as.list(snapper_maxima$max)), 20)) |>
  mutate(
    evt = map(.x = k, .f = ~ est_max_posterior(fit_evtg, k = .x)$estmax),
    efs = map(.x = k, .f = ~ est_max_posterior(fit_efs, k = .x)$estmax)
  )

p2_data_summary <-
  tibble(k = c(length(as.list(snapper_maxima$max)), 20)) |>
  mutate(
    evt = map(.x = k, .f = ~ est_max(fit_evtg, k = .x)),
    efs = map(.x = k, .f = ~ est_max(fit_efs, k = .x))
  ) |>
  unnest(cols = c(evt, efs)) |>
  mutate(pred = rep(c("fit", "lwr", "upr"), 2)) |>
  pivot_longer(cols = evt:efs) |>
  pivot_wider(values_from = value, names_from = pred) |>
  mutate(fact = paste0(name, k)) |>
  mutate(
    fact = factor(
      fact,
      levels = c("efs20", "efs14", "evt20", "evt14"),
      ordered = TRUE,
    )
  )


p_partb <-
  p2_data |>
  unnest(cols = c(evt, efs)) |>
  pivot_longer(cols = evt:efs) |>
  mutate(fact = paste0(name, k)) |>
  mutate(
    fact = factor(
      fact,
      levels = c("efs20", "efs14", "evt20", "evt14"),
      ordered = TRUE,
    )
  ) |>
  ggplot(aes(x = fact, y = value, fill = name)) +
  geom_hline(yintercept = 130, lty = 2, col = "grey70") +
  coord_flip() +
  see::geom_violinhalf() +
  geom_errorbar(
    aes(x = fact, ymin = lwr, ymax = upr),
    width = 0,
    data = p2_data_summary,
    inherit.aes = FALSE,
    linewidth = 2
  ) +
  geom_point(aes(y = fit), data = p2_data_summary, size = 4) +
  labs(x = NULL, y = expression(paste("Estimated ", L[max]))) +
  scale_linewidth_manual(values = c(`14` = 1, `20` = 3)) +
  scale_y_continuous(
    labels = label_number(suffix = "cm"),
    limits = layer_scales(p_snapper)$x$range$range
  ) +
  scale_x_discrete(
    labels = c(
      "evt14" = "EVT (k = 14)",
      "evt20" = "EVT (k = 20)",
      "efs14" = "EFS (k = 14)",
      "efs20" = "EFS (k = 20)"
    )
  ) +
  scale_fill_manual(values = c("evt" = evt_colour, "efs" = efs_colour)) +
  theme_classic(20) +
  theme(legend.position = "none")

snapper_img <- magick::image_read("data/snapper.png")

img_grob <- grid::rasterGrob(snapper_img, interpolate = TRUE)
text_grob <- grid::textGrob(
  "Chrysophrys auratus",
  gp = grid::gpar(fontsize = 15, fontface = "italic")
)

p_snapper_img <-
  p_snapper +
  annotation_custom(
    img_grob,
    xmin = 0,
    xmax = 70,
    ymin = layer_scales(p_snapper)$y$range$range[2] * 0.7, # Start at 70% of max y-value
    ymax = layer_scales(p_snapper)$y$range$range[2] # End at max y-value (top)
  ) +
  annotation_custom(
    text_grob,
    xmin = 0,
    xmax = 70,
    ymin = layer_scales(p_snapper)$y$range$range[2] * 0.68, # Just below image
    ymax = layer_scales(p_snapper)$y$range$range[2] * 0.7 # Up to image start
  )

p_snapper_final <-
  p_snapper_img +
  p_partb +
  patchwork::plot_layout(ncol = 1, heights = c(4, 1)) +
  patchwork::plot_annotation(tag_levels = "A")


ggsave(
  filename = "results/figures/manuscript_figures/snapper.png",
  plot = p_snapper_final,
  height = 12,
  width = 12
)

get_posterior(fit_evt)
get_posterior(fit_evt) |>
  pivot_longer(cols = c(mu, sigma, lambda)) |>
  summarise(
    fit = quantile(value, 0.5),
    lwr = quantile(value, 0.1),
    upr = quantile(value, 0.9),
    .by = name
  )
est_max(fit_evt, k = 14)
est_max(fit_evt, k = 20)
est_max(fit_evtg, k = 14)
est_max(fit_evtg, k = 20)
est_max(fit_efs, k = 14)
est_max(fit_efs, k = 20)
