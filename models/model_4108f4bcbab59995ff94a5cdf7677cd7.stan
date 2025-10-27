
    data {
      int<lower=1> k;                    // Number of competitions/samples
      int<lower=1> n_obs;                // Total number of observations
      vector[n_obs] x;                   // All recorded fish sizes
      array[k] int<lower=1> n_per_sample; // Number of fish per sample
      array[k] int<lower=1> start_idx;   // Starting index for each sample
    }
    
    parameters {
      real<lower=0> mu;                  // Mean of the fish size distribution
      real<lower=0.001> sigma;           // Standard deviation
      real<lower=1> lambda;              // Estimated mean number of fish caught per competition
    }
    
    model {
      // Priors
      mu ~ normal(30, 20) T[0.001,];
      sigma ~ normal(10, 5) T[0.001,];
      lambda ~ gamma(16,0.002);
      
      target += -lambda * k;
      
      {
        real norm_const = 1 - normal_cdf(0 | mu, sigma);  // Normalisation constant
        
        for(j in 1:k) {
          target += n_per_sample[j] * log(lambda);
          int end_idx = start_idx[j] + n_per_sample[j] - 1;
          real trunc_cdf = (normal_cdf(min(x[start_idx[j]:end_idx]) | mu, sigma) - normal_cdf(0 | mu, sigma)) / norm_const;
          target += lambda * trunc_cdf;
          
          for(i in 1:n_per_sample[j]) {
            int obs_idx = start_idx[j] + i - 1;
            real trunc_lpdf = normal_lpdf(x[obs_idx] | mu, sigma) - log(norm_const);
            target += trunc_lpdf;
          }
        }
      }
    }
