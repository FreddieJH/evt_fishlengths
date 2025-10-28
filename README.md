# Maximum Length Estimation for Fish Species

This repository contains code for estimating maximum body lengths of fish species using Extreme Value Theory (EVT) and and Exact-finite Sampling (EFS) approach.

## Overview

We use two main approaches to estimate maximum fish lengths:
1. Extreme Value Theory (EVT) using either the Generalized Extreme Value (GEV) or Gumbel distribution
2. Using Bayesian inference with some knowledge about the underlying distribution, which we term Exact Finite Sample (EFS) approach


The benefit of the EFS approach is that we can easily incorporate more than 

## File Structure

The code to produce the results and figures presented in the manuscript are arranged as follows:

```
.
├── models/
|   ├── efs.stan               # Exact finite sample Stan model
|   ├── evt.stan               # Extreme value theory (GEV distribution) Stan model
|   └── evt_gumbel.stan        # Extreme value theory (Gumbel distribution) Stan model
├── R/
│   ├── 00_pkgs.R              # packages to be used in the analysis
│   ├── 01_funcs.R             # custom functions used in the analysis
│   ├── 02_simulation.R        # simulate data (to be used in the sensitivity and concept plots)
│   ├── 03_truemax.R           # calculate the 'true' maximum for each simulation
|   ├── 04_model_prep.R        # prepare simulation data for model fitting
|   ├── 05_model_fitting.R     # model fitting of simulated data
|   ├── 06_model_checks.R      # check the fitting proceedure has performed correctly
|   ├── 07_posteriors.R        # extract and summarise fitted model posteriors
|   ├── 08_concept_plot.R      # plot showing the concept of the approaches
|   ├── 09_sensitivity_plot.R  # plot showing the sensitivity to various simulation parameters
|   └── 10_snapper_plot.R      # plot of case-study example (Australasian Snapper)
├── results/
│   └── data/    
|        ├── estmax_posterior.csv  # estimated max and max20 (20-sample maximum) for each scenario
|        ├── scenarios_truemax.csv  # 
|        ├── scenarios.csv          # 
|        └── posterior.parquet     # raw posterior draws (dataframe)
│   └── figures/    
│       └── manuscript_figures/ 
|           ├── concept.png
|           ├── sensitivity.png
|           └── snapper.png
│       └── supplementary_figures/ 
|           ├── XX.png
|           ├── XX.png
|           └── XX.png  
│       └── model_checks/  
|           ├── bayes_check.png
|           └── traceplots/
|         
```
