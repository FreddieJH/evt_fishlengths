The likelihood for the first EFS method, where we have a single maximum per sample is:

$$L=\lambda^ke^{-\lambda[k-\sum_{i=1}^kF(x_i,\theta)]}\prod_{i=1}^k f(x_i,\theta)$$

To get the log-likelihood we need to take the log of each side: 

$$LL=log(\lambda^ke^{-\lambda[k-\sum_{i=1}^kF(x_i,\theta)]}\prod_{i=1}^k f(x_i,\theta))$$

Following the product rule:

$$LL=log(\lambda^k) + log(e^{-\lambda[k-\sum_{i=1}^kF(x_i,\theta)]}) + log(\prod_{i=1}^k f(x_i,\theta))$$

For the first part of the RHS, can be re-written according to the power rule. 

$$LL=k \cdot log(\lambda) + log(e^{-\lambda[k-\sum_{i=1}^kF(x_i,\theta)]}) + log(\prod_{i=1}^k f(x_i,\theta))$$

For the second part of the RHS, can be re-written according to the log-of-base rule. 

$$LL=k \cdot log(\lambda) -\lambda[k-\sum_{i=1}^kF(x_i,\theta)] + log(\prod_{i=1}^k f(x_i,\theta))$$

For the third part of the RHS, can be re-written according to the product rule. 

$$LL=k \cdot log(\lambda) -\lambda[k-\sum_{i=1}^kF(x_i,\theta)] + \sum_{i=1}^k log(f(x_i,\theta))$$

We can now expand the second part of the RHS:

$$LL=k \cdot log(\lambda) -\lambda k+ \lambda\sum_{i=1}^kF(x_i,\theta) + \sum_{i=1}^k log(f(x_i,\theta))$$

Re-writing to make it easier to code:

$$LL=k \cdot log(\lambda) -\lambda k+ \sum_{i=1}^k (\lambda F(x_i,\theta) + log(f(x_i,\theta)))$$