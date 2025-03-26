
# calculating the maxima for each scenario
est_sim_max <- function(scenario_input){ 

    future::plan(multisession)
    sim_tbl_max <-
        scenario_input %>%
        dplyr::mutate(maxima = future_pmap(
            list(k, n_lambda, mu, sigma, min_size),
            sim_pois_truncnorm,
            .options = furrr::furrr_options(seed = TRUE)
        ))
    future::plan(sequential)

 return(sim_tbl_max)
}

