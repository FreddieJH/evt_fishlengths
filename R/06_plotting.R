# Plotting functions

source("R/02_functions.R")
source("R/04_evt_analysis.R")
source("R/05_numerical.R")

# Compare EVT and numerical estimation methods
compare_methods <- function(snapper_maxima) {
    est1 <- estimate_lmax_evt(snapper_maxima$max)
    est2 <- estimate_lmax_numerical(snapper_maxima$max)

    combined_est <- est1 %>%
        left_join(est2) %>%
        pivot_longer(cols = -x, names_to = c("model", ".value"), names_pattern = "(pgev|pmax)_(.*)") %>%
        mutate(model_clean = case_when(model == "pgev" ~ "EVT", model == "pmax" ~ "Numerical estimation"))

    p <- combined_est %>%
        ggplot(aes(x, y = fit, col = model_clean)) +
        geom_line() +
        geom_ribbon(aes(ymin = lwr, ymax = upr, fill = model_clean), alpha = 0.1) +
        geom_rug(aes(x = snapper_maxima$max), inherit.aes = FALSE) +
        scale_x_continuous(label = scales::label_number(suffix = "cm")) +
        labs(x = "Body size", y = "Pr(x>Lmax)") +
        theme_classic(20) +
        theme(legend.position = "inside", legend.position.inside = c(0, 1), legend.justification = c(-0.1, 1.5), legend.title = element_blank())

    ggsave("output/figures/snapper_maxlength_estimation.png", p, width = 10, height = 7)
}
