// Bayesian GEV model for estimating maximum length
data {
  int<lower=0> k;           // number of observations
  vector[k] x;              // observed maxima values
}

parameters {
  real mu;                  // location parameter
  real<lower=0> sigma;      // scale parameter (must be positive)
  real xi;                  // shape parameter
}

model {
  // Priors
  mu ~ normal(mean(x), 10);     // centered around sample mean with wide variance
  sigma ~ lognormal(0, 1);      // ensures positivity, reasonably diffuse
  xi ~ normal(0, 0.5);          // shape parameter typically small, centered at 0
  
  // GEV likelihood
  for (i in 1:k) {
    if (xi == 0) {
      target += -(log(sigma) + (x[i] - mu)/sigma + exp(-(x[i] - mu)/sigma));
    } else {
      real z = (x[i] - mu)/sigma;
      if (1 + xi * z > 0) {
        target += -(log(sigma) + (1 + 1/xi) * log1p(xi * z) + pow(1 + xi * z, -1/xi));
      } else {
        target += negative_infinity();  // outside support
      }
    }
  }
}
