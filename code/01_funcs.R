list.files(
  "data-raw/R",
  pattern = "\\.R$",
  full.names = TRUE,
  ignore.case = TRUE
) |>
  purrr::walk(source)

dtnorm <- function(x, mean, sd) {
  truncnorm::dtruncnorm(x = x, mean = mean, sd = sd, a = 0)
}
ptnorm <- function(q, mean, sd) {
  truncnorm::ptruncnorm(q = q, mean = mean, sd = sd, a = 0)
}
qtnorm <- function(p, mean, sd) {
  truncnorm::qtruncnorm(p = p, mean = mean, sd = sd, a = 0)
}
rtnorm <- function(n, mean, sd) {
  truncnorm::rtruncnorm(n = n, mean = mean, sd = sd, a = 0)
}
# modified from the evd package but made vectorised
dgev <- Vectorize(evd::dgev)
pdgev <- Vectorize(evd::pgev)
dgumbel <- Vectorize(evd::dgumbel)
pgumbel <- Vectorize(evd::pgumbel)

# # Function: given distribution mean and variance, what are the parameter values
# get_dist_pars <- function(distr, mean, variance) {
#   switch(
#     distr,
#     gamma = c((mean^2) / variance, mean / variance),
#     tnorm = c(mean, sqrt(variance)),
#     lnorm = {
#       log_var_ratio <- log(1 + variance / (mean^2))
#       c(log(mean) - log_var_ratio / 2, sqrt(log_var_ratio))
#     },
#     stop("Unsupported distribution: ", distr)
#   )
# }

# # PDF of maxima given the PDF and CDF of x
# g_max <- function(x, distr, n, par1, par2) {
#   f_x <- get(paste0("d", distr), mode = "function")
#   F_x <- get(paste0("p", distr), mode = "function")

#   # log to avoid very small values
#   log_g_max <- log(n) +
#     (n - 1) * log(F_x(x, par1, par2)) +
#     log(f_x(x, par1, par2))
#   exp(log_g_max)
# }

# # CDF for maxima
# G_max <- function(x, distr, n, par1, par2) {
#   F_x <- function(x) get(paste0("p", distr))(x, par1, par2)
#   return(F_x(x)^n)
# }

# # quantile function for maxima
# inverse_G_x <- function(
#   distr,
#   n,
#   par1,
#   par2,
#   p,
#   interval_lwr = 1,
#   interval_upr = 1000
# ) {
#   uniroot(
#     function(x) G_max(x, distr = distr, n = n, par1 = par1, par2 = par2) - p,
#     lower = interval_lwr,
#     upper = interval_upr
#   )$root
# }

# # expected value is just integral of x*g(x) where g(x) is the pdf of the max values
# expected_max <- function(distr, n, par1, par2) {
#   integrand <- function(x) {
#     x * g_max(x, distr = distr, n = n, par1 = par1, par2 = par2)
#   }
#   upper_bound <- get(paste0("q", distr))(0.9999999999999, par1, par2)
#   integrate(integrand, lower = 0, upper = upper_bound, rel.tol = 1e-6)$value
# }

# g_max_evt_posterior <- function(evt_posterior, mean, xmin, xmax, nsteps = 100) {
#   evt_posterior |>
#     select(-lp__) |>
#     expand_grid(x = seq(xmin, xmax, length.out = nsteps)) |>
#     mutate(
#       pdf = dgev(x = x, loc = loc, scale = scale, shape = shape),
#       cdf = pgev(q = x, loc = loc, scale = scale, shape = shape)
#     ) |>
#     summarise(
#       pdf_fit = quantile(pdf, 0.5),
#       pdf_lwr = quantile(pdf, 0.1),
#       pdf_upr = quantile(pdf, 0.9),
#       cdf_fit = quantile(cdf, 0.5),
#       cdf_lwr = quantile(cdf, 0.1),
#       cdf_upr = quantile(cdf, 0.9),
#       .by = x
#     )
# }

# g_max_evtg_posterior <- function(
#   evtg_posterior,
#   mean,
#   xmin,
#   xmax,
#   nsteps = 100
# ) {
#   evtg_posterior |>
#     select(-lp__) |>
#     expand_grid(x = seq(xmin, xmax, length.out = nsteps)) |>
#     mutate(
#       pdf = dgumbel(x, loc = loc, scale = scale),
#       cdf = pgumbel(x, loc = loc, scale = scale)
#     ) |>
#     summarise(
#       pdf_fit = quantile(pdf, 0.5),
#       pdf_lwr = quantile(pdf, 0.1),
#       pdf_upr = quantile(pdf, 0.9),
#       cdf_fit = quantile(cdf, 0.5),
#       cdf_lwr = quantile(cdf, 0.1),
#       cdf_upr = quantile(cdf, 0.9),
#       .by = x
#     )
# }

# g_max_efs_posterior <- function(efs_posterior, mean, xmin, xmax, nsteps = 100) {
#   efs_posterior |>
#     select(-lp__) |>
#     expand_grid(x = seq(xmin, xmax, length.out = nsteps)) |>
#     mutate(
#       pdf = g_max(x, "tnorm", lambda, mu, sigma),
#       pdf_under = dtnorm(x, mu, sigma),
#       cdf = G_max(x, "tnorm", lambda, mu, sigma)
#     ) |>
#     summarise(
#       pdf_fit = quantile(pdf, 0.5),
#       pdf_lwr = quantile(pdf, 0.1),
#       pdf_upr = quantile(pdf, 0.9),
#       pdf_fit_under = quantile(pdf_under, 0.5),
#       pdf_lwr_under = quantile(pdf_under, 0.1),
#       pdf_upr_under = quantile(pdf_under, 0.9),
#       cdf_fit = quantile(cdf, 0.5),
#       cdf_lwr = quantile(cdf, 0.1),
#       cdf_upr = quantile(cdf, 0.9),
#       .by = x
#     )
# }

evt_colour_dark <- "#2E86AB"
efs_colour_dark <- "#A23B72"
efsm_colour_dark <- "#F18F01"

evt_colour <- "#7DB3D3"
efs_colour <- "#C77BA0"
efsm_colour <- "#F4B942"

data_colour_nearmax <- "#89e4dcff"
data_colour_ismax <- "#000000"

sel_pal <- RColorBrewer::brewer.pal(name = "Dark2", n = 3)
