
import scipy.stats as stats

def expected_max_gamma(mean, variance, n):
    k = mean**2 / variance
    theta = variance / mean
    p = n / (n + 1)  # Quantile method
    return stats.gamma.ppf(p, a=k, scale=theta)
