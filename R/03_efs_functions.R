


multiple_efs_fit <- function(maxima, model_path = "models/efs.stan") {

  mod_results <- fit_maxima_model(maxima, model_path)

  mod_results$draws %>%
    mutate(
        max = pmap_dbl(
            .l = list(n = length(maxima)*lambda, mu = mu, sigma = sigma),
            .f = efs_expected_max
        ),
        max100 = pmap_dbl(
            .l = list(n = 100*lambda, mu = mu, sigma = sigma),
            .f = efs_expected_max
        )
    ) %>% 
    rename_with(~ paste0("est_", .x), setdiff(names(.), c(".chain", ".iteration", ".draw"))) 
}