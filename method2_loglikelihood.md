For the second EFS method, that is with multiple maxima per sample. The likelihood is: 

$$L=\lambda^{\sum_{j=1}^km_j}e^{-\lambda[k-\sum_{j=1}^kF(x_{m_j})]}\prod_{j=1}^k\prod_{i=1}^{m_j}f(x_{i,j})$$

To get the log-likelihood, we log both sides. 

$$LL=log(\lambda^{\sum_{j=1}^km_j}e^{-\lambda[k-\sum_{j=1}^kF(x_{m_j})]}\prod_{j=1}^k\prod_{i=1}^{m_j}f(x_{i,j}))$$

Using the product rule:

$$LL=log(\lambda^{\sum_{j=1}^km_j}) + log(e^{-\lambda[k-\sum_{j=1}^kF(x_{m_j})]}) + log(\prod_{j=1}^k\prod_{i=1}^{m_j}f(x_{i,j}))$$

Focusing on the first component on the RHS, we use the power rule:

$$LL={log(\lambda)\sum_{j=1}^km_j} + log(e^{-\lambda[k-\sum_{j=1}^kF(x_{m_j})]}) + log(\prod_{j=1}^k\prod_{i=1}^{m_j}f(x_{i,j}))$$

Focusing on the second component on the RHS, we use the log-of-base rule:

$$LL={log(\lambda)\sum_{j=1}^km_j} + -\lambda[k-\sum_{j=1}^kF(x_{m_j})] + log(\prod_{j=1}^k\prod_{i=1}^{m_j}f(x_{i,j}))$$

We can them expand the second component:

$$LL={log(\lambda)\sum_{j=1}^km_j} -\lambda k + \lambda \sum_{j=1}^kF(x_{m_j}) + log(\prod_{j=1}^k\prod_{i=1}^{m_j}f(x_{i,j}))$$

Focusing on the third component on the RHS, we use the product rule:

$$LL={log(\lambda)\sum_{j=1}^km_j} -\lambda k + \lambda \sum_{j=1}^kF(x_{m_j}) + \sum_{j=1}^k\sum_{i=1}^{m_j}log(f(x_{i,j}))$$

We can re-arrange this to make it easier to code up:

$$LL=-\lambda k + log(\lambda)\sum_{j=1}^km_j + \lambda \sum_{j=1}^kF(x_{m_j}) + \sum_{j=1}^k\sum_{i=1}^{m_j}log(f(x_{i,j}))$$

and move all the summations together:

$$LL=-\lambda k + \sum_{j=1}^k(log(\lambda)m_j + \lambda F(x_{m_j}) + \sum_{i=1}^{m_j}log(f(x_{i,j})))$$
