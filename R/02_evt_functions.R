evt_expected_max <- function(loc, scale, shape, k){
if (shape == 0) {
    expected_max <- loc + scale * log(k)
} else {
    expected_max <- loc + (scale / shape) * ((k^shape) - 1)
}
}



multiple_evt_fit <- function(maxima, model_path = "models/evt.stan") {

  mod_results <- fit_maxima_model(maxima, model_path)

  mod_results$draws %>%
    mutate(
        max = pmap_dbl(
            .l = list(k = length(maxima), loc = mu, scale = sigma, shape = xi),
            .f = expected_max_evt
        ),
        max100 = pmap_dbl(
            .l = list(k = 100, loc = mu, scale = sigma, shape = xi),
            .f = expected_max_evt
        )
    ) %>% 
    rename_with(~ paste0("est_", .x), setdiff(names(.), c(".chain", ".iteration", ".draw"))) 
}
