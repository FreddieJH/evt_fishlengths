
# calculating the maxima for each scenario
est_sim_max <- function(scenario_input){ 

    maxima_filepath <- "data/simulation/scenario_maxima.csv"
    # get stored maxima values
    if(file.exists(maxima_filepath)) {
        read_csv(maxima_filepath)
    }
    
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

