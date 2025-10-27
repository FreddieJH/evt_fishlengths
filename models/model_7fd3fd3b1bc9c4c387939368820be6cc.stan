data {
  int<lower=0> k;           // number of observations
  vector[k] x;              // observed maxima values
}

parameters {
  real loc;                  // location parameter
  real<lower=0> scale;      // scale parameter (must be positive)
  real shape;                  // shape parameter
}

model {
  // Priors
  loc ~ normal(mean(x), 10);     // centered around sample mean with wide variance
  scale ~ lognormal(0, 1);      // ensures positivity, reasonably diffuse
  shape ~ normal(0, 0.5);          // shape parameter typically small, centered at 0
  
  // GEV likelihood
  for (i in 1:k) {
    if (shape == 0) {
      target += -(log(scale) + (x[i] - loc)/scale + exp(-(x[i] - loc)/scale));
    } else {
      real z = (x[i] - loc)/scale;
      if (1 + shape * z > 0) {
        target += -(log(scale) + (1 + 1/shape) * log1p(shape * z) + pow(1 + shape * z, -1/shape));
      } else {
        target += negative_infinity();  // outside support
      }
    }
  }
}
