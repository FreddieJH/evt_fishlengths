data {
  int<lower=0> k;
  vector[k] x;
}

parameters {
  real loc;
  real<lower=0> scale;
  real shape;
}

model {
  // Priors
  loc ~ normal(mean(x), 10);
  scale ~ lognormal(0, 1);
  shape ~ normal(0, 0.5);
  
  // GEV likelihood
  for (i in 1:k) {
    real z = (x[i] - loc) / scale;
    
    if (abs(shape) < 1e-8) {
      // Gumbel case (shape ≈ 0)
      target += -log(scale) - z - exp(-z);
    } else {
      // General GEV case
      real t = 1 + shape * z;
      if (t > 0) {
        target += -log(scale) - (1 + 1/shape) * log(t) - pow(t, -1/shape);
      } else {
        target += negative_infinity();  // Reject invalid parameter values
      }
    }
  }
}
