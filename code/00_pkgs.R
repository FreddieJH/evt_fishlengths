pkgs <- c(
  # general
  "dplyr",
  "readr",
  "purrr",
  "arrow",
  "tibble",
  "forcats",
  "stringr",
  "tidyr",
  # models
  "truncnorm",
  "evd",
  "posterior",
  # parralellisation
  "future",
  "furrr",
  # visualisation
  "ggplot2",
  "patchwork",
  "scales",
  "RColorBrewer"
)

new_pkgs <- pkgs[!(pkgs %in% installed.packages()[, "Package"])]
if (length(new_pkgs)) {
  install.packages(new_pkgs)
}

purrr::walk(pkgs, ~ library(.x, character.only = TRUE))
rm(pkgs, new_pkgs)
