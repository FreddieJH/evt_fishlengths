
fit_mod <- function(maxima) {
    mod <- cmdstan_model("models/max_est.stan", stanc_options = list("O1"))

    fit <- mod$sample(
        data = list(k = length(maxima), x = maxima),
        iter_warmup = 1000,
        iter_sampling = 1000,
        chains = 4,
        parallel_chains = 4,
        refresh = 500
    )
}
# calc_max <- function(x, n, mu, sigma) {
#     fx <- dnorm(x, mu, sigma)
#     Fx <- pnorm(x, mu, sigma)
#     max <- (n * (Fx^(n - 1))) * fx
#     return(max)
# }

# find_max_x <- function(mu, sigma, n) {
#     optimize(function(x) calc_max(x, n, mu, sigma),
#         interval = c(mu - (5 * sigma), mu + (5 * sigma)),
#         maximum = TRUE
#     )$maximum
# }

# Stan fitting function
fit_stan_model <- function(rep, k, lambda, mu, maxima) {
mod <- cmdstan_model("models/max_est.stan", stanc_options = list("O1"))

    fit <- mod$sample(
        data = list(k = length(maxima), x = maxima),
        iter_warmup = 1000,
        iter_sampling = 1000,
        chains = 4,
        parallel_chains = 4,
        refresh = 500
    )

    filename <-
        paste0(
            "rep", rep,
            "_k", k,
            "_lambda", lambda,
            "_mu", mu
        )

    est_max_tbl <- fit %>% posterior::as_draws_df() 

    data.table::fwrite(fit$summary(), file = paste0("models/output_data/stan_outputs/", filename, ".csv"))
    data.table::fwrite(est_max_tbl, file = paste0("models/output_data/max_posterior/", filename, ".csv"))

    return(NULL)
}

run_scenarios <- function(scenarios){
processed_files <-
    list.files(
        "models/output_data/stan_outputs",
        full.names = FALSE
    ) %>%
    str_remove(".csv")

unprocessed_tbl <-
    scenarios %>% 
    mutate(maxvals = map(maxima, function(.x) .x$max)) %>% 
    filter(!filename %in% processed_files) 

if(nrow(unprocessed_tbl)){
    plan(multisession)
    future_pmap(
        .l = list(
            rep = unprocessed_tbl$rep,
            k = unprocessed_tbl$k,
            lambda = unprocessed_tbl$n_lambda,
            mu = unprocessed_tbl$mu,
            maxima = unprocessed_tbl$maxvals
        ),
        .f = fit_stan_model,
        .options = furrr_options(seed = TRUE),
        .progress = TRUE
    )
    plan(sequential)
}
}


combine_stan_results <- function() {
    filename_pattern <- "rep[0-9]+_k[0-9]+_lambda[0-9]+_mu[0-9]+"
    vroom::vroom(
        list.files("models/output_data/stan_outputs",
            pattern = filename_pattern,
            full.names = TRUE
        ),
        id = "filename",
        num_threads = 1
    ) %>%
        mutate(filename = str_extract(filename, filename_pattern)) %>%
        left_join(sim_tbl) %>%
        fwrite(file = "output/data/stan_output.csv")
}

if(file.exists("output/data/stan_output.csv")) {
    saved_files <-
        vroom::vroom("output/data/stan_output.csv") %>%
        pull(filename) %>%
        unique() %>%
        str_remove(".csv")

    unsaved_files <-
        scenarios %>%
        filter(!filename %in% saved_files) %>%
        pull(filename)

    if (length(unsaved_files)) {
        combine_stan_results()
    }
} else {
    combine_stan_results()
}
