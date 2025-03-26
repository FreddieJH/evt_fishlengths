# Maximum Length Estimation for Fish Species

This repository contains code for estimating maximum body lengths of fish species using Extreme Value Theory (EVT) and numerical methods.

## Overview

We use two main approaches to estimate maximum fish lengths:
1. Extreme Value Theory (EVT) using the Generalized Extreme Value (GEV) distribution
2. Numerical estimation using Bayesian inference

## Setup

### Required R Packages
```r
install.packages(c(
  "tidyverse",  # For data manipulation and visualization
  "evd",        # For extreme value distributions
  "patchwork",  # For combining plots
  "geomtextpath", # For text on paths in plots
  "multidplyr", # For parallel processing
  "furrr",      # For parallel processing
  "cmdstanr",   # For Bayesian inference
  "posterior",  # For handling MCMC output
  "truncnorm"   # For truncated normal distributions
))
```

## Analysis Steps

### 1. Simulated data
- Simulate maxima from truncated normal distributions using the `sim_pois_truncnorm()` function

### 2. EVT Approach
- Fit GEV distribution to sample maxima
- Estimate maximum length using EVT quantiles
- Code location: See `estimate-lmax` chunk

### 3. Numerical Estimation
- Use Bayesian model to estimate underlying distribution parameters
- Calculate maximum length using numerical optimization
- Code location: See Stan model in `models/max_est.stan`

### 4. Real Data Application
- Apply both methods to snapper (Pagrus auratus) data
- Compare EVT and numerical approaches
- Code location: See sections "Applying EVT to real data" and "Numerical estimation - real data"

## File Structure
```
.
├── README.md           # This file
├── main.qmd           # Main analysis document
├── models/
│   └── max_est.stan   # Stan model for numerical estimation
├── output/
│   ├── data/          # Processed data
│   ├── figures/       # Generated plots
│   └── simulated_data/# Simulation results
├── R/
│   ├── simulation/    # simulate data
│   ├── evt/            # functions to fit evt to data
│   ├── numerical/    # functions to fit numerical model to data
│   └── plotting/    # functions for plotting data
```

## Instructions for AI

1. Load required packages as specified in the setup section
2. Follow the analysis steps in order
3. Each step builds on previous results
4. Use provided functions in main.qmd:
   - `gev_fun()` for GEV calculations
   - `sim_pois_truncnorm()` for simulations
   - `find_max_x()` for numerical optimization


## Contributing

1. Fork the repository
2. Create a feature branch
3. Submit a pull request

## License

[Insert chosen license]
