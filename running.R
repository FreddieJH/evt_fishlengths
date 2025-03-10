source("scenarios.R")
library(rstan)

sim_pois_truncnorm <-
  function(seed, 
           k, # number of comps (therefore number of maxima)
           n_lambda, # poisson lambda = mean number of fish per competition
           mu, #cm = mean population fish size (truncated normal)
           sigma, 
           min_size #cm the smallest 'detectable' or 'catchable' fish
           ){
    
    set.seed(seed)
    n = rpois(k, lambda = n_lambda) # sample size of each competition
    
    truncnorm_max <- function(n, mu, sigma, lower) {
      max(truncnorm::rtruncnorm(n = n, 
                                a = lower, 
                                mean = mu, 
                                sd = sigma))
    }
    
    
    sapply(X = n, 
           FUN = truncnorm_max, 
           mu = mu, 
           sigma = sigma, 
           lower = min_size)
  }

# simulation over all scenarios
sim_tbl <- 
  expand_grid(
    seed = 1:10, 
    k = c(3, 5, 10, 20, 30, 100), 
    n_lambda = c(30, 50, 200),
    mu = c(20, 50, 100),
  ) |> 
  mutate(sigma = mu*0.34,
         min_size = 10) |> 
  mutate(maxima = purrr::pmap(
    .l = list(seed, k, n_lambda, mu, sigma, min_size), 
    .f = sim_pois_truncnorm))


stan_code <- "
data {
  int<lower=1> k;  // Number of competitions
  vector[k] x;  // Recorded largest fish sizes
}

parameters {
  real<lower=0> mu;    // Mean of the fish size distribution
  real<lower=0> sigma; // Standard deviation
  real<lower=1> lambda;     // Estimated mean number of fish caught per competition
}

model {
  // Priors
  mu ~ normal(30, 20);
  sigma ~ normal(10, 5);
  lambda ~ gamma(50, 0.1);  


  target += k*log(lambda) - lambda*k;

  for (i in 1:k) {

  target += lambda*normal_cdf(x[i] | mu, sigma) + normal_lpdf(x[i]|mu, sigma);

}

}
"

for(i in 1:nrow(sim_tbl)){
  
  if(file.exists(paste0("stan_outputs/i", i, ".csv"))) next
  
  maxima = unlist(sim_tbl$maxima[i])
  
  stan_data <- list(
    k = length(maxima),
    x = maxima
  )
  
  fit3 <- 
    stan(model_code = stan_code, 
       data = stan_data, 
       iter = 2000, 
       chains = 4, 
       cores = 4)
  
  
  write_csv(x = summarise_draws(fit), 
            file = paste0("stan_outputs/i", i, ".csv"))
  
  # sim_tbl$fit[i] <- summarise_draws(fit)
  
}


if(!file.exists("stan_output.csv")){
  sim_outputs <- 
    vroom::vroom(list.files("stan_outputs/", pattern = "\\.csv$", full.names = TRUE), 
                 id = "filename") %>%
    mutate(filename = basename(filename)) 
  
  
  sim_output_join <- 
    sim_outputs |> 
    mutate(sim_id = as.numeric(str_extract(filename, "\\d+"))) |> 
    left_join(sim_tbl |> 
                select(-maxima) |> 
                mutate(sim_id = row_number())) 
  
  write_csv(sim_output_join, file = "stan_output.csv")
} else {
  sim_output_join <- read_csv("stan_output.csv")
}


sim_output_join |> 
  filter(variable == "mu") |> 
  ggplot(aes(median, mu, col = n_lambda, shape = as.factor(k))) +
    geom_point(position = "jitter") +
    geom_abline(slope = 1) +
  labs(y = "Population mu", 
       x = "Estimated mu (posterior median)",
       title = "3 true means (20cm, 50cm and 100cm) - y axis is jittered")

ggsave("mu_posterior.png")

sim_output_join |> 
  filter(variable == "sigma") |> 
  ggplot(aes(median, sigma, col = n_lambda, shape = as.factor(k))) +
  geom_point(position = "jitter") +
  geom_abline(slope = 1)+
  labs(y = "Population sigma", 
       x = "Estimated sigma (posterior median)",
       title = "3 true sd (6.6cm, 17cm and 34cm) - y axis is jittered")

ggsave("sigma_posterior.png")

sim_output_join |> 
  filter(variable == "lambda") |> 
  ggplot(aes(median, n_lambda, col = k)) +
  geom_point() +
  geom_abline(slope = 1)+
  labs(y = "True poisson lambda", 
       x = "Estimated lambda (posterior median)",
       title = "3 true pois lambda (20, 50 and 100)")
ggsave("lambda_posterior.png")


sim_output_join |> 
  select(sim_id, variable, median, real_mu = mu) |> 
  pivot_wider(values_from = median, 
              names_from = variable) |> 
  ggplot(aes(mu, lambda, col = as.factor(real_mu))) +
  geom_point() +
  geom_abline(slope = 1)


