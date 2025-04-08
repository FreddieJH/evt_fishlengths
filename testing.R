

 bind_rows(evt_posteriors_summary |> mutate(type = "evt"), 
        efs_posteriors_summary |> mutate(type = "efs")) |> 
    left_join(expand_grid(scenarios_truemax, type = c("evt", "efs")))


# Load necessary libraries
library(ggplot2)
library(truncnorm)
library(dplyr)
library(tidyr)
library(patchwork)

# Function to calculate distribution parameters based on mean and CV
calc_dist_pars <- function(dist = c("gamma", "tnorm", "lnorm"), mean, cv = 0.34) {
  var <- (cv * mean)^2  # Calculate variance from CV
  
  if (dist == "gamma") {
    shape <- (mean^2) / var
    rate <- mean / var  # rate is 1/scale
    return(list(shape = shape, rate = rate))
  } else if (dist == "tnorm") {
    mu <- mean
    sigma <- sqrt(var)
    return(list(mu = mu, sigma = sigma))
  } else if (dist == "lnorm") {
    meanlog <- log(mean) - log(1 + var/(mean^2))/2
    sdlog <- sqrt(log(1 + var/(mean^2)))
    return(list(meanlog = meanlog, sdlog = sdlog))
  }
}

# Create a function to generate density curves
get_density_data <- function(mean_val, cv = 0.34) {
  # Define x range (from 0 to mean + 3*SD)
  var <- (cv * mean_val)^2
  sd_val <- sqrt(var)
  x_max <- mean_val + 3*sd_val
  x <- seq(0, x_max, length.out = 1000)
  
  # Calculate parameters for each distribution
  gamma_pars <- calc_dist_pars("gamma", mean_val, cv)
  tnorm_pars <- calc_dist_pars("tnorm", mean_val, cv)
  lnorm_pars <- calc_dist_pars("lnorm", mean_val, cv)
  
  # Calculate densities
  gamma_dens <- dgamma(x, shape = gamma_pars$shape, rate = gamma_pars$rate)
  tnorm_dens <- dtruncnorm(x, a = 0, b = Inf, mean = tnorm_pars$mu, sd = tnorm_pars$sigma)
  lnorm_dens <- dlnorm(x, meanlog = lnorm_pars$meanlog, sdlog = lnorm_pars$sdlog)
  
  # Create data frame
  data.frame(
    x = rep(x, 3),
    density = c(gamma_dens, tnorm_dens, lnorm_dens),
    distribution = rep(c("Gamma", "Truncated Normal", "Lognormal"), each = length(x)),
    mean = mean_val
  )
}

# Create a function to generate CDF curves
get_cdf_data <- function(mean_val, cv = 0.34) {
  # Define x range (from 0 to mean + 3*SD)
  var <- (cv * mean_val)^2
  sd_val <- sqrt(var)
  x_max <- mean_val + 3*sd_val
  x <- seq(0, x_max, length.out = 1000)
  
  # Calculate parameters for each distribution
  gamma_pars <- calc_dist_pars("gamma", mean_val, cv)
  tnorm_pars <- calc_dist_pars("tnorm", mean_val, cv)
  lnorm_pars <- calc_dist_pars("lnorm", mean_val, cv)
  
  # Calculate CDFs
  gamma_cdf <- pgamma(x, shape = gamma_pars$shape, rate = gamma_pars$rate)
  tnorm_cdf <- ptruncnorm(x, a = 0, b = Inf, mean = tnorm_pars$mu, sd = tnorm_pars$sigma)
  lnorm_cdf <- plnorm(x, meanlog = lnorm_pars$meanlog, sdlog = lnorm_pars$sdlog)
  
  # Create data frame
  data.frame(
    x = rep(x, 3),
    cdf = c(gamma_cdf, tnorm_cdf, lnorm_cdf),
    distribution = rep(c("Gamma", "Truncated Normal", "Lognormal"), each = length(x)),
    mean = mean_val
  )
}

# Function to generate upper tail comparison
get_tail_comparison <- function(mean_val, cv = 0.34, prob_threshold = 0.95) {
  # Define x range focusing on upper tail
  var <- (cv * mean_val)^2
  sd_val <- sqrt(var)
  
  # Calculate parameters for each distribution
  gamma_pars <- calc_dist_pars("gamma", mean_val, cv)
  tnorm_pars <- calc_dist_pars("tnorm", mean_val, cv)
  lnorm_pars <- calc_dist_pars("lnorm", mean_val, cv)
  
  # Calculate quantiles for prob_threshold
  gamma_q <- qgamma(prob_threshold, shape = gamma_pars$shape, rate = gamma_pars$rate)
  tnorm_q <- qtruncnorm(prob_threshold, a = 0, b = Inf, mean = tnorm_pars$mu, sd = tnorm_pars$sigma)
  lnorm_q <- qlnorm(prob_threshold, meanlog = lnorm_pars$meanlog, sdlog = lnorm_pars$sdlog)
  
  # Find maximum quantile for x range
  x_start <- min(gamma_q, tnorm_q, lnorm_q) * 0.95
  x_end <- max(gamma_q, tnorm_q, lnorm_q) * 1.2
  
  x <- seq(x_start, x_end, length.out = 500)
  
  # Calculate densities
  gamma_dens <- dgamma(x, shape = gamma_pars$shape, rate = gamma_pars$rate)
  tnorm_dens <- dtruncnorm(x, a = 0, b = Inf, mean = tnorm_pars$mu, sd = tnorm_pars$sigma)
  lnorm_dens <- dlnorm(x, meanlog = lnorm_pars$meanlog, sdlog = lnorm_pars$sdlog)
  
  # Create data frame
  data.frame(
    x = rep(x, 3),
    density = c(gamma_dens, tnorm_dens, lnorm_dens),
    distribution = rep(c("Gamma", "Truncated Normal", "Lognormal"), each = length(x)),
    mean = mean_val
  )
}

# Function to simulate expected max values
simulate_max_values <- function(mean_values, cv = 0.34, sample_sizes = c(10, 20, 50, 100)) {
  results <- data.frame()
  
  for (mean_val in mean_values) {
    for (n in sample_sizes) {
      # Calculate parameters
      gamma_pars <- calc_dist_pars("gamma", mean_val, cv)
      tnorm_pars <- calc_dist_pars("tnorm", mean_val, cv)
      lnorm_pars <- calc_dist_pars("lnorm", mean_val, cv)
      
      # Function to estimate expected max through simulation
      estimate_max <- function(dist, n_samples = 10000) {
        if (dist == "gamma") {
          samples <- replicate(n_samples, max(rgamma(n, shape = gamma_pars$shape, rate = gamma_pars$rate)))
        } else if (dist == "tnorm") {
          samples <- replicate(n_samples, max(rtruncnorm(n, a = 0, b = Inf, mean = tnorm_pars$mu, sd = tnorm_pars$sigma)))
        } else if (dist == "lnorm") {
          samples <- replicate(n_samples, max(rlnorm(n, meanlog = lnorm_pars$meanlog, sdlog = lnorm_pars$sdlog)))
        }
        mean(samples)
      }
      
      # Estimate expected max for each distribution
      gamma_max <- estimate_max("gamma")
      tnorm_max <- estimate_max("tnorm")
      lnorm_max <- estimate_max("lnorm")
      
      # Add to results
      results <- rbind(results, data.frame(
        mean = mean_val,
        sample_size = n,
        gamma_max = gamma_max,
        tnorm_max = tnorm_max,
        lnorm_max = lnorm_max,
        gamma_ratio = gamma_max / tnorm_max,
        lnorm_ratio = lnorm_max / tnorm_max
      ))
    }
  }
  
  return(results)
}

# Define mean values to test
mean_values <- c(5, 10, 25, 50)
cv <- 0.34

# Generate density data for each mean value
density_data <- do.call(rbind, lapply(mean_values, get_density_data, cv = cv))

# Generate CDF data for each mean value
cdf_data <- do.call(rbind, lapply(mean_values, get_cdf_data, cv = cv))

# Generate tail comparison data
tail_data <- do.call(rbind, lapply(mean_values, get_tail_comparison, cv = cv))

# Simulate expected max values
max_results <- simulate_max_values(mean_values, cv = cv)

# Reshape max results for plotting
max_long <- max_results %>%
  pivot_longer(cols = c(gamma_ratio, lnorm_ratio),
               names_to = "distribution", 
               values_to = "ratio")

# Create plots
p1 <- ggplot(density_data, aes(x = x, y = density, color = distribution)) +
  geom_line(size = 1) +
  facet_wrap(~ mean, scales = "free", 
             labeller = labeller(mean = function(x) paste("Mean =", x))) +
  labs(title = "PDF Comparison (CV = 0.34)",
       x = "x", y = "Density") +
  theme_minimal() +
  theme(legend.position = "bottom")

p2 <- ggplot(cdf_data, aes(x = x, y = cdf, color = distribution)) +
  geom_line(size = 1) +
  facet_wrap(~ mean, scales = "free_x", 
             labeller = labeller(mean = function(x) paste("Mean =", x))) +
  labs(title = "CDF Comparison (CV = 0.34)",
       x = "x", y = "Cumulative Probability") +
  theme_minimal() +
  theme(legend.position = "bottom")

p3 <- ggplot(tail_data, aes(x = x, y = density, color = distribution)) +
  geom_line(size = 1) +
  facet_wrap(~ mean, scales = "free", 
             labeller = labeller(mean = function(x) paste("Mean =", x))) +
  labs(title = "Upper Tail Comparison (CV = 0.34)",
       x = "x", y = "Density") +
  theme_minimal() +
  theme(legend.position = "bottom")

p4 <- ggplot(max_long, aes(x = mean, y = ratio, color = distribution, shape = factor(sample_size))) +
  geom_point(size = 3) +
  geom_line(aes(linetype = factor(sample_size))) +
  labs(title = "Ratio of Expected Max Values (relative to Truncated Normal)",
       x = "Mean", y = "Ratio to Truncated Normal",
       shape = "Sample Size", linetype = "Sample Size") +
  scale_color_discrete(labels = c("Gamma/TNorm", "Lognormal/TNorm")) +
  theme_minimal() +
  theme(legend.position = "bottom")

# Arrange plots
p5 <- (p1 | p2) / (p3 | p4)

ggsave(
    filename = paste0("output/figures/understanding_dist_tails.png"),
    plot = p5,
    height = 10,
    width = 15
)
