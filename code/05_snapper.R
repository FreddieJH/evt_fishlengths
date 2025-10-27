library(tidyverse)
library(scales)
library(grid)
library(patchwork)

source("code/00a_truncnorm_functions.R")
source("code/00b_model_functions.R")
source("code/00c_expectedmax_functions.R")
source("code/00d_colours.R")

kg2cm <- function(w, a = 0.04478, b = 2.673) ((w * 1000) / a)^(1 / b)
# cm2kg = function(l, a = 0.04478, b = 2.673) (a*l^b)/1000
# cm2kg = function(l, a = 0.02239, b = 2.97) (a*l^b)/1000
# cm2kg = function(l, a = 0.01, b = 3) (a*l^b)/1000

# snapper_maxima <- tibble(
#   type = c(
#     rep("length", 8),
#     rep("weight", 4)
#   ),
#   max = c(
#     91.3,
#     102,
#     112,
#     107,
#     107,
#     99.2,
#     95,
#     82.2,
#     kg2cm(c(11.8, 18.4, 16.5, 17.2))
#   )
# )

# lm1 <-
# xx |>
#   mutate(w_g = weight_kg/1000) %>%
# lm(log(w_g) ~ log(length_cm), data = .)

# xx <- read_csv("data/snapper_maxima.csv")
# xx |> ggplot(aes(length_cm, weight_kg)) +
#   geom_point() +
#   geom_line(data = tibble(length_cm = 0:130) |> mutate(weight_kg = cm2kg(length_cm)))

snapper_maxima <-
  read_csv("data/snapper_maxima.csv") |>
  mutate(
    max = case_when(
      !is.na(length_cm) ~ length_cm,
      is.na(length_cm) ~ kg2cm(weight_kg, a = 0.01, b = 3)
    )
  )

# maxima is either a vector of maxmia (for EVT and EFS) list(x1, x2, x3) or c(x1, x2, x3)
# or it is a list of vectors (for EFSMM) list(c(x1, x2, x3), x(y1, y2), c(z1, z2, z3, z4))
get_maxima <- function(maxima, model_type) {
  maxima_list <- as.list(maxima)
  gamma_shape = 5
  gamma_rate = 0.1
  mu_prior_mean = 30
  mu_prior_sd = 30 / 1.5

  mod_dat <-
    list(
      x = maxima_list |> unlist(),
      n_obs = length(maxima_list |> unlist()),
      n_per_sample = lapply(maxima_list, length) |> unlist(),
      start_idx = cumsum(lapply(maxima_list, length) |> unlist()),
      k = length(maxima_list)
    )

  init <- if (model_type == "evt") {
    function() {
      list(
        mu = median(maxima_list |> unlist()),
        sigma = sd(maxima_list |> unlist()),
        xi = 0
      )
    }
  } else {
    function() {
      list(
        mu = median(maxima_list |> unlist()),
        sigma = sd(maxima_list |> unlist()),
        lambda = 1000
      )
    }
  }

  stan_file <- switch(
    model_type,
    "evt" = cmdstanr::write_stan_file(stan_code_evt, dir = "models"),
    "efs" = cmdstanr::write_stan_file(
      stan_code_efs(gamma_shape, gamma_rate, mu_prior_mean, mu_prior_sd),
      dir = "models"
    ),
    "efsmult" = cmdstanr::write_stan_file(
      stan_code_efsmult(gamma_shape, gamma_rate, mu_prior_mean, mu_prior_sd),
      dir = "models"
    )
  )

  mod <- cmdstanr::cmdstan_model(stan_file, stanc_options = list("O1"))

  fit <- mod$sample(
    data = mod_dat,
    chains = 4,
    init = init,
    iter_warmup = 2000,
    iter_sampling = 2000
  )

  traceplot <- bayesplot::mcmc_trace(fit$draws())

  summary_pdf <-
    fit |>
    posterior::as_draws_df() %>%
    dplyr::as_tibble() %>%
    dplyr::select(-lp__) |>
    expand_grid(size = 0:round(((max(snapper_maxima$max)) * 1.5))) |>
    mutate(
      pdf = if (model_type == "evt") {
        0
      } else {
        pmap_dbl(list(x = size, mean = mu, sd = sigma), dtnorm)
      },
      pdf_max = if (model_type == "evt") {
        pmap_dbl(
          list(x = size, loc = loc, scale = scale, shape = shape),
          evd::dgev
        )
      } else {
        pmap_dbl(
          .l = list(
            x = size,
            distr = "tnorm",
            n = length(maxima_list) * lambda,
            par1 = mu,
            par2 = sigma
          ),
          .f = g_max
        )
      }
    ) |>
    summarise(
      p_fit = mean(pdf),
      p_lwr = quantile(pdf, 0.025),
      p_upr = quantile(pdf, 0.975),
      pmax_fit = mean(pdf_max),
      pmax_lwr = quantile(pdf_max, 0.025),
      pmax_upr = quantile(pdf_max, 0.975),
      .by = size
    )

  return(list(fit = fit, traceplot = traceplot, output_table = summary_pdf))
}

fit1 <- get_maxima(snapper_maxima$max, "evt")
fit2 <- get_maxima(snapper_maxima$max, "efs")


evt_maxest <- function(k) {
  fit1$fit |>
    posterior::as_draws_df() %>%
    dplyr::as_tibble() %>%
    dplyr::select(-lp__) |>
    mutate(
      pdf = pmap_dbl(
        list(loc = loc, scale = scale, shape = shape, k = k),
        expected_max_evt
      )
    ) |>
    summarise(
      max_fit = mean(pdf),
      max_lwr = quantile(pdf, 0.025),
      max_upr = quantile(pdf, 0.975),
    ) %>%
    as.numeric()
}

efs_maxest <- function(k) {
  fit2$fit |>
    posterior::as_draws_df() %>%
    dplyr::as_tibble() %>%
    dplyr::select(-lp__) |>
    mutate(
      pdf = pmap_dbl(
        .l = list(
          distr = "tnorm",
          n = k * lambda,
          mean = mu,
          variance = sigma^2
        ),
        .f = expected_max_fromsim
      )
    ) |>
    summarise(
      max_fit = mean(pdf),
      max_lwr = quantile(pdf, 0.1),
      max_upr = quantile(pdf, 0.9),
    ) %>%
    as.numeric()
}


p_snapper <-
  fit1[[3]] |>
  ggplot(aes(size, pmax_fit)) +
  geom_ribbon(
    aes(ymin = p_lwr, ymax = p_upr),
    alpha = 0.3,
    fill = "#A23B72",
    data = fit2[[3]]
  ) +
  geom_line(aes(y = p_fit), col = "#A23B72", lty = 2, data = fit2[[3]]) +
  geom_ribbon(
    aes(ymin = pmax_lwr, ymax = pmax_upr),
    alpha = 0.5,
    fill = "#7DB3D3"
  ) +
  geom_line(col = "#2E86AB", linewidth = 2) +
  geom_ribbon(
    aes(ymin = pmax_lwr, ymax = pmax_upr),
    alpha = 0.5,
    data = fit2[[3]],
    fill = "#C77BA0"
  ) +
  # geom_line(col = "#C77BA0", data = fit2[[3]], linewidth = 2) +
  # geom_ribbon(
  #   aes(ymin = pmax20_lwr, ymax = pmax20_upr),
  #   alpha = 0.5,
  #   data = fit2[[3]],
  #   fill = "#C77BA0"
  # ) +
  geom_line(col = "#C77BA0", data = fit2[[3]], linewidth = 2) +
  geom_rug(
    aes(x = x),
    data = tibble(x = snapper_maxima$max),
    inherit.aes = FALSE
  ) +
  scale_x_continuous(label = label_number(suffix = "cm")) +
  labs(x = "Body size", y = "Probability density") +
  theme_classic(20)

p2_data <-
  tibble(k = c(length(as.list(snapper_maxima$max)), 20)) |>
  mutate(evt = map(.x = k, .f = evt_maxest), efs = map(.x = k, .f = efs_maxest))

p_partb <-
  p2_data |>
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
  ) |>
  ggplot(aes(x = fit, y = fact, col = name, linewidth = as.factor(k))) +
  geom_vline(xintercept = 130, lty = 2, col = "grey70") +
  geom_point(size = 5) +
  geom_errorbarh(aes(xmin = lwr, xmax = upr), height = 0) +

  labs(y = NULL, x = expression(paste("Estimated ", L[max]))) +
  scale_linewidth_manual(values = c(`14` = 1, `20` = 3)) +
  scale_x_continuous(
    labels = label_number(suffix = "cm"),
    limits = layer_scales(p_snapper)$x$range$range
  ) +
  scale_y_discrete(
    labels = c(
      "evt14" = "EVT (k = 14)",
      "evt20" = "EVT (k = 20)",
      "efs14" = "EFS (k = 14)",
      "efs20" = "EFS (k = 20)"
    )
  ) +
  scale_color_manual(values = c("evt" = evt_colour, "efs" = efs_colour)) +
  theme_classic(20) +
  theme(legend.position = "none")

snapper_img <- magick::image_read("docs/pagrus_aurata_RSS.jpg")

img_grob <- grid::rasterGrob(snapper_img, interpolate = TRUE)
text_grob <- grid::textGrob(
  "Chrysophrys aurata",
  gp = grid::gpar(fontsize = 15, fontface = "italic")
)

p_snapper_img <-
  p_snapper +
  annotation_custom(
    img_grob,
    xmin = 0,
    xmax = 50,
    ymin = layer_scales(p_snapper)$y$range$range[2] * 0.7, # Start at 70% of max y-value
    ymax = layer_scales(p_snapper)$y$range$range[2] # End at max y-value (top)
  ) +
  annotation_custom(
    text_grob,
    xmin = 0,
    xmax = 50,
    ymin = layer_scales(p_snapper)$y$range$range[2] * 0.68, # Just below image
    ymax = layer_scales(p_snapper)$y$range$range[2] * 0.7 # Up to image start
  )

p_snapper_final <-
  p_snapper_img +
  p_partb +
  patchwork::plot_layout(ncol = 1, heights = c(4, 1))


ggsave(
  filename = "results/figures/snapper.png",
  plot = p_snapper_final,
  height = 10,
  width = 10
)
