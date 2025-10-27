source("R/01_funcs.R")
source("R/02_simulation.R")


scenarios_truemax <-
  scenarios |>
  mutate(
    true_max = pmap_dbl(
      .l = list(k, lambda, dist_name, dist_mean),
      .f = \(k, lambda, dist_name, dist_mean) {
        pars <- get_dist_pars(distr = dist_name, mean = dist_mean)
        pdf <- \(x) get(paste0("d", dist_name))(x, pars[1], pars[2])
        cdf <- \(x) get(paste0("p", dist_name))(x, pars[1], pars[2])
        mode <- mode_f(\(x) g(x, n = k * lambda, cdf = cdf, pdf = pdf))
        return(mode)
      }
    ),
    true_max20 = pmap_dbl(
      .l = list(lambda, dist_name, dist_mean),
      .f = \(lambda, dist_name, dist_mean) {
        pars <- get_dist_pars(distr = dist_name, mean = dist_mean)
        pdf <- \(x) get(paste0("d", dist_name))(x, pars[1], pars[2])
        cdf <- \(x) get(paste0("p", dist_name))(x, pars[1], pars[2])
        mode <- mode_f(\(x) g(x, n = 20 * lambda, cdf = cdf, pdf = pdf))
        return(mode)
      }
    )
  )
