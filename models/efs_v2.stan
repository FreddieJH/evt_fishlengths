data {
  int<lower=1> k;                    // Number of sample maxima
  int<lower=1> total_fish;           // Total number of observed maxima
  int<lower=1> m[k];                 // Number of maxima in each sample
  int<lower=1> sample_id[total_fish]; // Sample ID for each fish observation
  vector[total_fish] x;              // Recorded fish sizes (largest m fish from each sample)
}

parameters {
  real<lower=0> mu;     // Mean of the fish size distribution
  real<lower=0> sigma;  // Standard deviation
  real<lower=1> lambda; // Estimated mean number of fish caught per sample
}

transformed parameters {
  // Calculate the minimum size for each sample (smallest of the m largest fish)
  vector[k] x_min;
  
  for (j in 1:k) {
    real min_val = positive_infinity();
    for (i in 1:total_fish) {
      if (sample_id[i] == j) {
        min_val = min(min_val, x[i]);
      }
    }
    x_min[j] = min_val;
  }
}

model {
  // Priors
  mu ~ normal(30, 20);
  sigma ~ normal(10, 5);
  lambda ~ gamma(50, 0.1);
  
  // Sum of m_j across all samples
  int sum_m = sum(m);
  
  // First part of the likelihood: λ^(sum of m_j)
  target += sum_m * log(lambda);
  
  // Second part: exp(-λ[k-sum(F(x_m_j))])
  real sum_cdf = 0;
  for (j in 1:k) {
    sum_cdf += normal_cdf(x_min[j], mu, sigma);
  }
  target += -lambda * (k - sum_cdf);
  
  // Third part: product of f(x_i,j) for all observed fish
  for (i in 1:total_fish) {
    target += normal_lpdf(x[i] | mu, sigma);
  }
}