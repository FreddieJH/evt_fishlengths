data {
    int<lower=0> N;             // number of observations
    vector[N] x;                // input data (maxima)
}

parameters {
    real<lower=0> mu;           // location parameter
    real<lower=0> sigma;        // scale parameter
    real xi;                    // shape parameter (can be negative, zero, or positive)
}

model {
    // Priors
    mu ~ normal(0, 100);        // Weakly informative prior for location
    sigma ~ cauchy(0, 5);       // Half-Cauchy prior for scale (positive)
    xi ~ normal(0, 10);         // Weakly informative prior for shape

    // Likelihood using GEV distribution
    for (i in 1:N) {
        // GEV CDF and log-likelihood calculation
        real z = 1 + xi * (x[i] - mu) / sigma;
        
        if (z > 0) {
            target += log(1/sigma) - (1/xi + 1) * log(z) - pow(z, -1/xi);
        }
    }
}

generated quantities {
    vector[N] log_lik;          // Log-likelihood for each observation
    array[N] real ypred;              // Posterior predictive samples

    for (i in 1:N) {
        // Log-likelihood calculation
        real z = 1 + xi * (x[i] - mu) / sigma;
        
        if (z > 0) {
            log_lik[i] = log(1/sigma) - (1/xi + 1) * log(z) - pow(z, -1/xi);
        } else {
            log_lik[i] = negative_infinity();
        }

        // Posterior predictive sampling
        ypred[i] = mu + (sigma / xi) * (pow(uniform_rng(0, 1), -xi) - 1);
    }
}
