# fit the GEV distribution
fit_gev <- function(maxima) {
  # Check for sufficient data
  if (length(maxima) < 3 || all(is.na(maxima))) {
    warning("Insufficient data for fitting GEV.")
    return(tibble(loc = NA_real_, scale = NA_real_, shape = NA_real_))
  }

  # Try fitting GEV
  gev_fit <- suppressWarnings(try(evd::fgev(x = maxima), silent = TRUE))

  # Check for errors in GEV fit
  if (inherits(gev_fit, "try-error")) {
    warning("GEV fitting failed, returning NA values.")
    return(tibble(loc = NA_real_, scale = NA_real_, shape = NA_real_))
  }

  # Extract estimates
  est <- gev_fit$estimate
  se <- gev_fit$std.err
  se <- gev_fit$std.err

  # Validate parameter estimates
  if (est["scale"] <= 0 || abs(est["shape"]) > 2 || any(is.na(se)) || any(se > 10)) {
    warning("Unreliable GEV parameter estimates, returning NA values.")
    return(tibble(loc = NA_real_, scale = NA_real_, shape = NA_real_))
  }

  # Return results
  return(tibble(loc = est["loc"], scale = est["scale"], shape = est["shape"],
  loc_se= se["loc"], scale_se = se["scale"], shape_se = se["shape"]))
}


gev_exp_max <- function(loc, scale, shape, k){
if (shape == 0) {
    expected_max <- loc + scale * log(k)
} else {
    expected_max <- loc + (scale / shape) * ((k^shape) - 1)
}
}

