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
   shape ~ normal(0, 0.1);
  scale ~ normal(0, 10); 
  loc ~ normal(0, 50);  

  // GEV likelihood
  for (i in 1:k) {
    if (shape == 0) {
    real t_i = exp(-(x[i]-loc)/scale);
    target += log(1/scale) + (shape+1)*(log(t_i)) - t_i;
    } else {
      real t_i = pow(1 + (shape * ((x[i] - loc)/scale)), -1/shape);
      target += log(1/scale) + (shape+1)*(log(t_i)) - t_i;
    }
  }
}
