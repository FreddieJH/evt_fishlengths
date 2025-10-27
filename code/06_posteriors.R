# Script requirements  --------------------------------------------------------------

# source("code/00_pkgs.R")
# source("code/01_funcs.R")
# source("code/03_truemax.R")
source("code/04_model_fitting.R")


posterior |>
  filter(model_id == "evt") |>
  pivot_wider(names_from = par, values_from = value) |>
  left_join(scenarios_truemax |> select(-samples), by = "scenario_id") |>
  filter(scenario_id == 10) |>
  expand_grid(
    x = seq(concept_popln_mean * 0.8, concept_popln_mean * 3, by = 0.1)
  ) |>
  mutate(
    pdf = pmap_dbl(
      .l = list(x = x, loc = loc, scale = scale, shape = shape),
      .f = dgev
    ),
    cdf = pmap_dbl(
      .l = list(q = x, loc = loc, scale = scale, shape = shape),
      .f = pgev
    )
  )

xx |>
  summarise(
    pdf_lwr = quantile(pdf, 0.1),
    pdf_fit = quantile(pdf, 0.5),
    pdf_upr = quantile(pdf, 0.9),
    cdf_lwr = quantile(cdf, 0.1),
    cdf_fit = quantile(cdf, 0.5),
    cdf_upr = quantile(cdf, 0.9),
    .by = c(x)
  )


if (!file.exists("data/estmax_posterior.csv")) {
  # EVT METHOD
  plan(multisession)
  evt_ests <-
    posterior |>
    filter(model_id == "evt") |>
    pivot_wider(names_from = par, values_from = value) |>
    left_join(scenarios_truemax |> select(-samples), by = "scenario_id") |>
    mutate(
      est_max20 = future_pmap_dbl(
        .l = list(loc, scale, shape),
        .f = \(loc, scale, shape) {
          qgev(p = 0.95, loc = loc, scale = scale, shape = shape)
        }
      ),
      est_max = future_pmap_dbl(
        .l = list(loc, scale, shape, k),
        .f = \(loc, scale, shape, k) {
          qgev(p = 1 - (1 / k), loc = loc, scale = scale, shape = shape)
        }
      )
    ) |>
    summarise(
      est_max20_lwr = quantile(est_max20, 0.1),
      est_max20_fit = quantile(est_max20, 0.5),
      est_max20_upr = quantile(est_max20, 0.9),
      est_max_lwr = quantile(est_max, 0.1),
      est_max_fit = quantile(est_max, 0.5),
      est_max_upr = quantile(est_max, 0.9),
      .by = c(scenario_id, model_id)
    )
  plan(sequential)
  # EVT (GUMBEL) METHOD
  plan(multisession)
  evtg_ests <-
    posterior |>
    filter(model_id == "evtg") |>
    pivot_wider(names_from = par, values_from = value) |>
    left_join(scenarios_truemax |> select(-samples), by = "scenario_id") |>
    mutate(
      est_max20 = future_pmap_dbl(
        .l = list(loc, scale),
        .f = \(loc, scale) {
          qgumbel(p = 0.95, loc = loc, scale = scale)
        }
      ),
      est_max = future_pmap_dbl(
        .l = list(loc, scale, k),
        .f = \(loc, scale, k) {
          qgumbel(p = 1 - (1 / k), loc = loc, scale = scale)
        }
      )
    ) |>
    summarise(
      est_max20_lwr = quantile(est_max20, 0.1),
      est_max20_fit = quantile(est_max20, 0.5),
      est_max20_upr = quantile(est_max20, 0.9),
      est_max_lwr = quantile(est_max, 0.1),
      est_max_fit = quantile(est_max, 0.5),
      est_max_upr = quantile(est_max, 0.9),
      .by = c(scenario_id, model_id)
    )
  plan(sequential)

  # EFS METHOD
  plan(multisession)
  efs_ests <-
    posterior |>
    filter(model_id %in% c("efs", "efsm")) |>
    pivot_wider(names_from = par, values_from = value) |>
    left_join(
      scenarios_truemax |> select(-c(samples, lambda)),
      by = "scenario_id"
    ) |>
    mutate(
      est_max = future_pmap_dbl(
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
      ),
      est_max20 = future_pmap_dbl(
        .l = list(mu, sigma, lambda, k),
        .f = \(mu, sigma, lambda, k) {
          cdf <- \(x) {
            ptnorm(q = x, mean = mu, sd = sigma)
          }
          pdf <- \(x) {
            dtnorm(x = x, mean = mu, sd = sigma)
          }
          gmax <- \(x) {
            g(x = x, n = lambda * 20, cdf = cdf, pdf = pdf)
          }
          mode_f(gmax)
        }
      )
    ) |>
    summarise(
      est_max_lwr = quantile(est_max, 0.1),
      est_max_fit = quantile(est_max, 0.5),
      est_max_upr = quantile(est_max, 0.9),
      est_max20_lwr = quantile(est_max20, 0.1),
      est_max20_fit = quantile(est_max20, 0.5),
      est_max20_upr = quantile(est_max20, 0.9),
      .by = c(scenario_id, model_id)
    )
  plan(sequential)

  estmax_posterior <-
    bind_rows(evt_ests, evtg_ests, efs_ests)

  write_csv(estmax_posterior, "data/estmax_posterior.csv")
} else {
  estmax_posterior <- read_csv("data/estmax_posterior.csv")
}

p_estmax <-
  estmax_posterior |>
  pivot_longer(
    cols = starts_with("est_"),
    names_to = c("source", ".value"),
    names_pattern = "est_(max20|max)_(.*)"
  ) |>
  ggplot(aes(x, fit)) +
  geom_ribbon(aes(ymin = lwr, ymax = upr, fill = model_id), alpha = 0.3) +
  geom_path(aes(col = type)) +
  facet_wrap(
    . ~ type_full,
    scales = "free",
  ) +
  scale_fill_manual(
    values = c("evt" = evt_colour, "efs" = efs_colour, "efsm" = efsm_colour),
    labels = c("evt" = "EVT", "efs" = "EFS", "efsm" = "EFSMM")
  ) +
  scale_colour_manual(
    values = c("evt" = evt_colour, "efs" = efs_colour, "efsm" = efsm_colour),
    labels = c("evt" = "EVT", "efs" = "EFS", "efsm" = "EFSMM")
  ) +
  labs(x = "Size (cm)", y = "Density", fill = NULL, col = NULL) +
  theme_classic(20) +
  theme(legend.position = "bottom")

# # Calculate the PDF and CDF of maxima at size X
# if (!file.exists("data/estmax_posteriors.rds")) {
#   future::plan(future::multisession)

#   model_fits_temp <-
#     posterior |>
#     mutate(
#       maxima_min = map_dbl(data, ~ min(.x$top1)),
#       maxima_max = map_dbl(data, ~ max(.x$top1)),
#       xmin = maxima_min * 0.8,
#       xmax = maxima_max * 1.2
#     )

#   # Vectorise the furrr_options creation
#   opts <- furrr_options(
#     globals = c(
#       "g_max",
#       "G_max",
#       "ptnorm",
#       "dtnorm",
#       "dgev",
#       "pgev",
#       "dgumbel",
#       "pgumbel",
#       "g_max_evt_posterior",
#       "g_max_evtg_posterior",
#       "g_max_efs_posterior"
#     ),
#     packages = c("dplyr", "purrr", "tidyr"),
#     seed = TRUE
#   )

#   estmax_posteriors <-
#     model_fits_temp |>
#     head(1) |>
#     mutate(
#       posteriors = furrr::future_pmap(
#         list(
#           evt_posterior,
#           evtg_posterior,
#           efs_posterior,
#           efsm_posterior,
#           dist_mean,
#           xmin,
#           xmax
#         ),
#         \(evt, evtg, efs, efsm, mean, xmin, xmax) {
#           list(
#             g_max_evt = g_max_evt_posterior(evt, mean, xmin, xmax),
#             g_max_evtg = g_max_evtg_posterior(evtg, mean, xmin, xmax),
#             g_max_efs = g_max_efs_posterior(efs, mean, xmin, xmax),
#             g_max_efsm = g_max_efs_posterior(efsm, mean, xmin, xmax)
#           )
#         },
#         .options = opts
#       )
#     ) |>
#     unnest_wider(posteriors)

#   future::plan(future::sequential())
#   write_rds(estmax_posteriors, "data/estmax_posteriors.rds")
# } else {
#   estmax_posteriors <- read_rds("data/estmax_posteriors.rds")
# }

# quantile_evt <- function(evt_posterior, q = 0.95) {
#   evt_posterior |>
#     mutate(
#       q = pmap_dbl(
#         .l = list(p = q, loc = loc, scale = scale, shape = shape),
#         .f = qgev
#       )
#     ) |>
#     summarise(
#       fit = quantile(q, 0.5),
#       lwr = quantile(q, 0.1),
#       upr = quantile(q, 0.9)
#     )
# }

# quantile_evtg <- function(evt_posterior, q = 0.95) {
#   evt_posterior |>
#     mutate(
#       q = pmap_dbl(
#         .l = list(p = q, loc = loc, scale = scale),
#         .f = qgumbel
#       )
#     ) |>
#     summarise(
#       fit = quantile(q, 0.5),
#       lwr = quantile(q, 0.1),
#       upr = quantile(q, 0.9)
#     )
# }

# quantile_efs <- function(posterior, q = 0.95) {
#   posterior |>
#     mutate(
#       q = pmap_dbl(
#         .l = list(
#           distr = "tnorm",
#           n = lambda,
#           par1 = mu,
#           par2 = sigma,
#           p = q
#         ),
#         .f = inverse_G_x
#       )
#     ) |>
#     summarise(
#       fit = quantile(q, 0.5),
#       lwr = quantile(q, 0.1),
#       upr = quantile(q, 0.9)
#     )
# }

# modal_max <- function(distr, nk, par1, par2) {
#   g_max_new <- function(x) {
#     g_max(x, distr = distr, n = nk, par1 = par1, par2 = par2)
#   }
#   qfunc <- function(p) get(paste0("q", distr))(p, par1, par2)
#   optimise(
#     f = g_max_new,
#     interval = c(qfunc(0.9), qfunc(1 - 1 / (2 * nk))),
#     maximum = TRUE
#   )$maximum
# }

# mode_efs <- function(posterior, k) {
#   posterior |>
#     mutate(
#       q = pmap_dbl(
#         .l = list(
#           distr = "tnorm",
#           n = lambda * k,
#           par1 = mu,
#           par2 = sigma
#         ),
#         .f = expected_max
#       )
#     ) |>
#     summarise(
#       fit = quantile(q, 0.5),
#       lwr = quantile(q, 0.1),
#       upr = quantile(q, 0.9)
#     )
# }

# get_max_error <- function(posterior_name, quanile_func) {
#   estmax_posteriors |>
#     select(all_of(c(
#       "k",
#       "lambda",
#       "dist_name",
#       "dist_mean",
#       posterior_name
#     ))) |>
#     expand_grid(q = c(1 - (1 / k), 0.95, 0.5)) |>
#     mutate(value = map2(!!sym(posterior_name), q, quanile_func)) |>
#     select(-!!sym(posterior_name)) |>
#     unnest(value)
# }

# get_max_error_EFS <- function(posterior_name, modal_func) {
#   estmax_posteriors |>
#     select(all_of(c(
#       "k",
#       "lambda",
#       "dist_name",
#       "dist_mean",
#       posterior_name
#     ))) |>
#     mutate(value = map2(!!sym(posterior_name), k, modal_func)) |>
#     select(-!!sym(posterior_name)) |>
#     unnest(value)
# }

# if (!file.exists("data/scenarios_estmax.csv")) {
#   scenarios_estmax_quantile <-
#     bind_rows(
#       get_max_error("evt_posterior", quantile_evt) |>
#         mutate(type = "evt"),
#       get_max_error("evtg_posterior", quantile_evtg) |>
#         mutate(type = "evtg"),
#       get_max_error("efs_posterior", quantile_efs) |>
#         mutate(type = "efs"),
#       get_max_error("efsm_posterior", quantile_efs) |>
#         mutate(type = "efsm")
#     )

#   scenarios_estmax_mode <-
#     bind_rows(
#       get_max_error_EFS("efs_posterior", mode_efs) |>
#         mutate(type = "efs_mode"),
#       get_max_error_EFS("efsm_posterior", mode_efs) |>
#         mutate(type = "efsm_mode")
#     )

#   scenarios_estmax <-
#     bind_rows(scenarios_estmax_quantile, scenarios_estmax_mode) |>
#     left_join(
#       scenarios_truemax |>
#         select(k, lambda, dist_name, dist_mean, true_max, true_max20)
#     )
#   write_csv(scenarios_estmax, "data/scenarios_estmax.csv")
# } else {
#   scenarios_estmax <- read_csv("data/scenarios_estmax.csv")
# }

# estmax_posteriors %>%
#   select(k, lambda, dist_name, dist_mean, g_max_evt, g_max_efs, g_max_efsm) %>%
#   unnest(cols = contains("g_max_"), names_sep = "_") %>%
#   select(-c(g_max_efs_x, g_max_efsm_x)) %>%
#   rename(x = g_max_evt_x) %>%
#   pivot_longer(cols = contains("g_max_")) %>%
#   mutate(
#     type = str_extract(name, "(?<=g_max_).*(?=_)"),
#     pred_name = str_extract(name, "[^_]*$")
#   ) %>%
#   select(-name) %>%
#   filter(
#     pred_name != "under",
#     lambda == 1000,
#     k == 10,
#     str_detect(type, "_pdf")
#   ) %>%
#   mutate(type = str_remove(type, "_pdf")) |>
#   pivot_wider(names_from = pred_name, values_from = value) %>%
#   mutate(
#     type_full = paste0(
#       case_when(
#         dist_name == "lnorm" ~ "Lognormal",
#         dist_name == "gamma" ~ "Gamma",
#         dist_name == "tnorm" ~ "Truncated-Normal"
#       ),
#       " (mean = ",
#       dist_mean,
#       "cm)"
#     ),
#     # Convert to factor with desired order
#     type_full = factor(
#       type_full,
#       levels = c(
#         paste0("Gamma (mean = ", c("10cm)", "50cm)", "100cm)")),
#         paste0("Lognormal (mean = ", c("10cm)", "50cm)", "100cm)")),
#         paste0("Truncated-Normal (mean = ", c("10cm)", "50cm)", "100cm)"))
#       )
#     )
#   ) |>
#   ggplot(aes(x, fit)) +
#   geom_ribbon(aes(ymin = lwr, ymax = upr, fill = type), alpha = 0.3) +
#   geom_path(aes(col = type)) +
#   facet_wrap(
#     . ~ type_full,
#     scales = "free",
#   ) +
#   scale_fill_manual(
#     values = c("evt" = evt_colour, "efs" = efs_colour, "efsm" = efsm_colour),
#     labels = c("evt" = "EVT", "efs" = "EFS", "efsm" = "EFSMM")
#   ) +
#   scale_colour_manual(
#     values = c("evt" = evt_colour, "efs" = efs_colour, "efsm" = efsm_colour),
#     labels = c("evt" = "EVT", "efs" = "EFS", "efsm" = "EFSMM")
#   ) +
#   labs(x = "Size (cm)", y = "Density", fill = NULL, col = NULL) +
#   theme_classic(20) +
#   theme(legend.position = "bottom")

# estmax_posteriors %>%
#   select(k, lambda, dist_name, dist_mean, g_max_evt, g_max_efs, g_max_efsm) %>%
#   unnest(cols = contains("g_max_"), names_sep = "_") %>%
#   select(-c(g_max_efs_x, g_max_efsm_x)) %>%
#   rename(x = g_max_evt_x) %>%
#   pivot_longer(cols = contains("g_max_")) %>%
#   mutate(
#     type = str_extract(name, "(?<=g_max_).*(?=_)"),
#     pred_name = str_extract(name, "[^_]*$")
#   ) %>%
#   select(-name) %>%
#   filter(
#     pred_name != "under",
#     lambda == 1000,
#     k == 10,
#     str_detect(type, "_pdf")
#   ) %>%
#   mutate(type = str_remove(type, "_pdf")) |>
#   pivot_wider(names_from = pred_name, values_from = value) %>%
#   mutate(
#     type_full = paste0(
#       ifelse(
#         dist_name == "lnorm",
#         "Lognormal",
#         ifelse(
#           dist_name == "gamma",
#           "Gamma",
#           ifelse(dist_name == "tnorm", "Truncated-Normal", NA)
#         )
#       ),
#       " (mean = ",
#       dist_mean,
#       "cm)"
#     )
#   ) |>
#   ggplot(aes(x, fit)) +
#   geom_ribbon(aes(ymin = lwr, ymax = upr, fill = type), alpha = 0.3) +
#   geom_path(aes(col = type)) +
#   facet_grid(
#     dist_name ~ dist_mean,
#     scales = "free",
#     labeller = labeller(
#       dist_mean = function(value) {
#         paste0("Mean = ", value, "cm")
#       },
#       dist_name = function(value) {
#         x <- ifelse(
#           value == "lnorm",
#           "Lognormal",
#           ifelse(
#             value == "gamma",
#             "Gamma",
#             ifelse(value == "tnorm", "Truncated-Normal", NA)
#           )
#         )
#         return(x)
#       }
#     )
#   ) +
#   scale_fill_manual(
#     values = c("evt" = evt_colour, "efs" = efs_colour, "efsm" = efsm_colour),
#     labels = c("evt" = "EVT", "efs" = "EFS", "efsm" = "EFSMM")
#   ) +
#   scale_colour_manual(
#     values = c("evt" = evt_colour, "efs" = efs_colour, "efsm" = efsm_colour),
#     labels = c("evt" = "EVT", "efs" = "EFS", "efsm" = "EFSMM")
#   ) +
#   labs(x = "Size (cm)", y = "Density", fill = NULL, col = NULL) +
#   theme_classic(20) +
#   theme(legend.position = "bottom")
