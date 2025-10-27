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




# Function: given distribution mean and variance, what are the parameter values
get_dist_pars <- function(distr, mean) {
  variance = (mean*0.34)^2
  switch(
    distr,
    gamma = c((mean^2) / variance, mean / variance),
    tnorm = c(mean, sqrt(variance)),
    lnorm = {
      log_var_ratio <- log(1 + variance / (mean^2))
      c(log(mean) - log_var_ratio / 2, sqrt(log_var_ratio))
    },
    stop("Unsupported distribution: ", distr)
  )
}

# get the mode of a function
mode_f <- function(f, lwr = -500, upr = 500) {
  # Evaluate at several points to find best starting region
  x_grid <- seq(lwr, upr, length.out = 100)
  f_vals <- sapply(x_grid, f)
  best_idx <- which.max(f_vals)

  # Narrow the search around the best point
  if (best_idx == 1) {
    search_lwr <- lwr
    search_upr <- x_grid[min(best_idx + 10, length(x_grid))]
  } else if (best_idx == length(x_grid)) {
    search_lwr <- x_grid[max(best_idx - 10, 1)]
    search_upr <- upr
  } else {
    search_lwr <- x_grid[max(best_idx - 10, 1)]
    search_upr <- x_grid[min(best_idx + 10, length(x_grid))]
  }

  optimise(f, interval = c(search_lwr, search_upr), maximum = TRUE)$maximum
}

# PDF of the maxima
g <- function(x, n, cdf, pdf) {
  n * cdf(x)^(n - 1) * pdf(x)
}


validate_fit <- function(fit) {
  valid_names <- c("evt", "evt_gumbel", "efs", "efsmm", "maxima")

  # Check fit is a list
  if (!is.list(fit)) {
    stop("'fit' must be a named list of CmdStanMCMC objects", call. = FALSE)
  }

  # Check fit has names
  if (is.null(names(fit)) || any(names(fit) == "")) {
    stop("'fit' must be a named list", call. = FALSE)
  }

  # Check all names are valid
  invalid_names <- setdiff(names(fit), valid_names)
  if (length(invalid_names) > 0) {
    stop(
      sprintf(
        "Invalid model names: %s. Must be one of: %s",
        paste(invalid_names, collapse = ", "),
        paste(valid_names, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  # Check all elements (except maxima) are CmdStanMCMC objects
  is_cmdstan <- vapply(
    fit[names(fit) != "maxima"],
    function(x) {
      inherits(x, "CmdStanMCMC")
    },
    logical(1)
  )

  if (!all(is_cmdstan)) {
    stop("All elements of 'fit' must be CmdStanMCMC objects", call. = FALSE)
  }

  invisible(TRUE)
}

check_cmdstan <- function() {
  if (!requireNamespace("cmdstanr", quietly = TRUE)) {
    stop(
      "Package 'cmdstanr' is required. Install with: 
         install.packages('cmdstanr', repos = c('https://stan-dev.r-universe.dev', getOption('repos')))"
    )
  }

  if (is.null(cmdstanr::cmdstan_version(error_on_NA = FALSE))) {
    stop("CmdStan is not installed. Install with: cmdstanr::install_cmdstan()")
  }
}


get_posterior <- function(fit) {
  fit_slim <- fit[names(fit) != "maxima"]
  validate_fit(fit_slim)

  # Extract posteriors for all models
  output_list <-
    purrr::map(fit_slim, \(model_fit) {
      tryCatch(
        {
          posterior::as_draws_df(model_fit) |>
            dplyr::as_tibble()
        },
        error = function(e) {
          stop(
            sprintf("Failed to extract posterior samples: %s", e$message),
            call. = FALSE
          )
        }
      )
    })
  names(output_list) <- names(fit_slim)
  output_list[["maxima"]] <- fit[["maxima"]]
  return(output_list)
}

fit_max_model <- function(
  length_maxima,
  model_type = c("evt", "evt_gumbel", "efs", "efsmm"),
  chains = 4,
  iter_warmup = 2000,
  iter_sampling = 1000,
  adapt_delta = 0.999,
  max_treedepth = 12
) {
  check_cmdstan()
  # Input validation
  checkmate::assert(
    checkmate::test_numeric(
      length_maxima,
      finite = TRUE,
      any.missing = FALSE,
      min.len = 1
    ),
    checkmate::test_list(length_maxima, types = "numeric", min.len = 1),
    combine = "or",
    .var.name = "length_maxima"
  )

  # If list, check at least one vector has length > 1 (for EFSMM)
  if (is.list(length_maxima)) {
    if (!any(lengths(length_maxima) > 1)) {
      stop(
        "When length_maxima is a list, at least one vector must have length > 1",
        call. = FALSE
      )
    }
    # Check all elements are numeric and finite
    if (
      !all(vapply(
        length_maxima,
        function(x) {
          is.numeric(x) && all(is.finite(x)) && !anyNA(x)
        },
        logical(1)
      ))
    ) {
      stop(
        "All elements in length_maxima list must be numeric, finite, and non-missing",
        call. = FALSE
      )
    }
  }

  # Validate model_type selection
  all_models <- c("evt", "evt_gumbel", "efs", "efsmm")
  if (missing(model_type)) {
    if (is.list(length_maxima)) {
      model_type <- all_models
    } else {
      model_type <- all_models[which(all_models != "efsmm")]
    }
  } else {
    model_type <- match.arg(model_type, several.ok = TRUE)
  }

  # Check EFSMM only used with list input
  if ("efsmm" %in% model_type && !is.list(length_maxima)) {
    stop(
      "model_type 'efsmm' requires length_maxima to be a list",
      call. = FALSE
    )
  }

  checkmate::assert_int(chains, lower = 1)
  checkmate::assert_int(iter_warmup, lower = 100)
  checkmate::assert_int(iter_sampling, lower = 100)

  # Convert to list format
  maxima_list <- if (is.list(length_maxima)) {
    length_maxima
  } else {
    as.list(length_maxima)
  }

  # Fit each model
  fits <- lapply(model_type, function(mtype) {
    fit_single_model(
      maxima_list = maxima_list,
      model_type = mtype,
      chains = chains,
      iter_warmup = iter_warmup,
      iter_sampling = iter_sampling,
      adapt_delta = adapt_delta,
      max_treedepth = max_treedepth
    )
  })

  names(fits) <- model_type
  fits[["maxima"]] <- length_maxima
  return(fits)
}

fit_single_model <- function(
  maxima_list,
  model_type,
  chains,
  iter_warmup,
  iter_sampling,
  adapt_delta,
  max_treedepth
) {
  if (model_type != "efsmm" & is.list(maxima_list)) {
    maxima_list <- unlist(lapply(maxima_list, FUN = max))
  }
  mod_dat <- list(
    x = unlist(maxima_list),
    n_obs = length(unlist(maxima_list)),
    n_per_sample = lengths(maxima_list),
    start_idx = cumsum(c(0, lengths(maxima_list)[-length(maxima_list)])) + 1,
    k = length(maxima_list)
  )

  init_func <- function(type, maxima_median) {
    if (type %in% c("evt", "evt_gumbel")) {
      function(chain_id) {
        list(loc = maxima_median, scale = 10, shape = 0)
      }
    } else {
      function(chain_id) {
        list(mu = maxima_median, sigma = 10, lambda = 100)
      }
    }
  }

  model_file <- system.file(
    "stan",
    paste0(ifelse(model_type == "efsmm", "efs", model_type), ".stan"),
    package = "fishmax"
  )

  if (!file.exists(model_file) || model_file == "") {
    stop(
      "Stan model file not found. Available files: ",
      paste(
        list.files(system.file("stan", package = "fishmax")),
        collapse = ", "
      ),
      "\nLooking for: ",
      paste0(ifelse(model_type == "efsmm", "efs", model_type), ".stan"),
      call. = FALSE
    )
  }

  mod <- cmdstanr::cmdstan_model(model_file)

  fit <- mod$sample(
    data = mod_dat,
    chains = chains,
    init = init_func(model_type, median(unlist(maxima_list))),
    iter_warmup = iter_warmup,
    iter_sampling = iter_sampling,
    adapt_delta = adapt_delta,
    max_treedepth = max_treedepth
  )

  return(fit)
}


evt_colour_dark <- "#2E86AB"
efs_colour_dark <- "#A23B72"
efsm_colour_dark <- "#F18F01"

evt_colour <- "#7DB3D3"
efs_colour <- "#C77BA0"
efsm_colour <- "#F4B942"

data_colour_nearmax <- "#89e4dcff"
data_colour_ismax <- "#000000"

sel_pal <- RColorBrewer::brewer.pal(name = "Dark2", n = 3)
