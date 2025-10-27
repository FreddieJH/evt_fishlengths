data {
  int<lower=0> k;           // number of observations
  vector[k] x;              // observed maxima values
}

parameters {
  real loc;                  // location parameter
  real<lower=0.00001> scale;      // scale parameter (must be positive)
}

model {
  // Priors
  loc ~ normal(100, 30);
  scale ~ lognormal(0, 2); 

  // GUMBEL likelihood
x ~ gumbel(loc, scale);
}
