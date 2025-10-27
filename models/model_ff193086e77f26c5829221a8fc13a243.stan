
    data {
      int<lower=1> k;                    // Number of competitions/samples
      int<lower=1> n_obs;                // Total number of observations
      array[n_obs] real x;                   // All recorded fish sizes
      array[k] rintal n_per_sample; // Number of fish per sample
      array[k] rintal start_idx;   // Starting index for each sample
    }
    
    parameters {
      real<lower=0> mu;                  // Mean of the fish size distribution
      real<lower=0.001> sigma;           // Standard deviation
      real<lower=1, upper=10000> lambda;              // Estimated mean number of fish caught per competition
    }
    
    model {
      // Priors
      mu ~ normal(30, 20) T[0.001,];
      sigma ~ normal(10, 5) T[0.001,];
      lambda ~ gamma(16,0.002);
      
      target += -lambda * k;
      
      {
        
      // area of normal CDF above zero (normalisation constant)
      real norm_const = 1 - normal_cdf(0 | mu, sigma);  
        
      // looping over each sample
        for(j in 1:k) {
          
        target += n_per_sample[j] * log(lambda);

          // we use the smallest value since we know all other values are smaller than x_m
          // we could use the start_indx if we know the x values are in increasing order
          // but this is the safest way
          real smallest_in_m = min(x[start_idx[j]:(start_idx[j] + n_per_sample[j] - 1)]);

          // Truncated normal CDF
          real trunc_cdf = (normal_cdf(smallest_in_m | mu, sigma) - normal_cdf(0 | mu, sigma)) / norm_const;

          target += lambda * trunc_cdf;

          // Add safeguard against overflow
          if (trunc_cdf > 0 && trunc_cdf < 1) {
            target += lambda * trunc_cdf;
          } else {
            target += negative_infinity();  // Reject this sample
          }
          
          // looping over each maxima within the sample j
          for(i in 1:n_per_sample[j]) {
            int obs_idx = start_idx[j] + i - 1;
            real trunc_lpdf = normal_lpdf(x[obs_idx] | mu, sigma) - log(norm_const);
            target += trunc_lpdf;
          }
        }
      }
    }
