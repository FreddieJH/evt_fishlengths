data {
  int<lower=1> k;  // Number of competitions
  vector[k] x;  // Recorded largest fish sizes
}

parameters {
  real<lower=0> mu;    // Mean of the fish size distribution
  real<lower=0> sigma; // Standard deviation
  real<lower=1> lambda; // Estimated mean number of fish caught per competition
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




