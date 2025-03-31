fit_stan_model <- function(rep, k, lambda, mu, maxima, type) {
  
  mod <- 
    paste0("models/", type, ".stan") %>% 
    cmdstan_model(stanc_options = list("O1"))
  
  if(type == "evt"){
    init_func = function() {
      list(
        mu = median(maxima),
        sigma = sd(maxima),
        xi = 0
      )
    }
  } else {
    init_func = function() {
      list(
        mu = median(maxima),
        sigma = sd(maxima), 
        lambda = 1000
      )
    }
  }
  
  fit <- mod$sample(
    data = list(k = length(maxima), x = maxima),
    iter_warmup = 1000,
    iter_sampling = 4000,
    chains = 4,
    parallel_chains = 4,
    refresh = 500, 
    init = init_func
  )
  
  filename <-
    paste0(
      "rep", rep,
      "_k", k,
      "_lambda", lambda,
      "_mu", mu
    )
  
  draws_obj <- posterior::as_draws_df(fit)
  summary_obj <- fit$summary()
  
  summary_path <- paste0("models/output_data/", type,"/summary/", filename, ".csv")
  draws_path <- paste0("models/output_data/", type,"/draws/", filename, ".csv")
  
  data.table::fwrite(summary_obj, file = summary_path)
  data.table::fwrite(draws_obj, file = draws_path)
  
  return(NULL)
}

run_scenarios <- function(scenarios, type = c("evt", "finite_sample")){
  processed_files <-
    list.files(
      paste0("models/output_data/", type,"/summary"),
      full.names = FALSE
    ) %>%
    str_remove(".csv")
  
  unprocessed_tbl <-
    scenarios %>% 
    mutate(maxvals = map(maxima, function(.x) .x$max)) %>% 
    filter(!filename %in% processed_files) 
  
  if(nrow(unprocessed_tbl)){
    
    plan(multisession)
    if(type == "evt"){
      furrr::future_pmap(
      .l = list(
        rep = unprocessed_tbl$rep,
        k = unprocessed_tbl$k,
        lambda = unprocessed_tbl$n_lambda,
        mu = unprocessed_tbl$mu,
        maxima = unprocessed_tbl$maxvals,
        type = "evt"
      ),
      .f = fit_stan_model,
      .options = furrr_options(seed = TRUE),
      .progress = TRUE)
      } else {
        future_pmap(
          .l = list(
            rep = unprocessed_tbl$rep,
            k = unprocessed_tbl$k,
            lambda = unprocessed_tbl$n_lambda,
            mu = unprocessed_tbl$mu,
            maxima = unprocessed_tbl$maxvals, 
            type = "finite_sample"
          ),
          .f = fit_stan_model,
          .options = furrr_options(seed = TRUE),
          .progress = TRUE
        )
    }
    plan(sequential)
  }
}


combine_stan_results <- function(type = c("evt", "finite_sample")) {
  filename_pattern <- "rep[0-9]+_k[0-9]+_lambda[0-9]+_mu[0-9]+"
  vroom::vroom(
    list.files(paste0("models/output_data/", type,"/summary"),
               pattern = filename_pattern,
               full.names = TRUE
    ),
    id = "filename",
    num_threads = 1
  ) %>%
    mutate(filename = str_extract(filename, filename_pattern)) %>%
    left_join(sim_tbl) %>%
    fwrite(file = paste0("output/data/", type,"-summary.csv"))
}

combine_csv <- function(type = c("evt", "finite_sample")){
  if(file.exists(paste0("output/data/", type,"-summary.csv"))) {
    saved_files <-
      vroom::vroom(paste0("output/data/", type,"-summary.csv")) %>%
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
}

if (!file.exists("output/simulated_data/simulation_scenarios.csv")) {
  scenarios <-
    expand_grid(
      rep = 1:100,
      k = c(3, 5, 10, 20, 30, 100),
      n_lambda = c(30, 50, 200, 2000, 10000),
      mu = c(20, 50, 100),
    ) |>
    mutate(
      sigma = mu * 0.34,
      min_size = 10,
      filename = paste0(
        "rep", rep,
        "_k", k,
        "_lambda", n_lambda,
        "_mu", mu
      )
    ) |>
    est_sim_max()
  
  if (!dir.exists("output/simulated_data/")) dir.create("output/simulated_data/")
  
  scenarios %>%
    tidyr::unnest(cols = maxima) %>%
    data.table::fwrite("output/simulated_data/simulation_scenarios.csv")
} else {
  scenarios <- 
    data.table::fread("output/simulated_data/simulation_scenarios.csv") %>% 
    nest(.by =-c(maxima, n), .key = "maxima")
}

run_scenarios(scenarios = scenarios, type = "evt")
run_scenarios(scenarios = scenarios, type = "finite_samples")
combine_csv(type = "evt")
combine_csv(type = "finite_samples")



