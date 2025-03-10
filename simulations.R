library(truncnorm)  

# Simulating fishing competitions
sim_pois_truncnorm <-
  function(seed, k, n_lambda, mu, sigma, min_size){
    
    # k = 5 # number of comps (therefore number of maxima)
    # n_lambda = 30 # poisson lambda = mean number of fish per competition
    
    # mu = 60 #cm = mean population fish size (truncated normal)
    # sigma = 60*0.34
    # min_size = 10 #cm the smallest 'detectable' or 'catchable' fish
    
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

# example simulation
sim_pois_truncnorm(seed = 1, k = 10, n_lambda = 50, mu = 66, sigma = 22, min_size = 10)

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



