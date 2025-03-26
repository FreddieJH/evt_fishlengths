The goal of these notes is to think through inference based on the largest observation in a finite sample of fish.  

Let $x$ represent the size of an individual and $f(x)$ denote the pdf.  The cumulative distribution function is $F(x)$ which gives us the probability that the size is less than or equal to $x$.  This distribution depends on some parameters that we'd like to infer, but for now, to keep the notation simple, we'll ignore that and introduce these parameters later.

If we sample $n$ fish, the probability that all of them are smaller than $x$ is $F(x)^n$.  Note that if all of them are smaller than $x$, then the maximum size is also smaller, i.e. $P(\text{largest fish has size }<x|n \text{ fish sampled})=F(x)^n$.  To obtain the pdf for the max, say $f_n(x)$, we differentiate and find that $f_n(x)=nF(x)^{n-1}f(x)$.

Extreme value distributions are derived by asking what happens as $n\to\infty$ and asserting that $F(x)^n = F(a_nx+b_n)$ for some sequence $a_n, b_n$.  This is not necessarily the case - apparently this isn't true for Poisson, Geometric, etc.  But it works for normal, lognormal, Gamma...However, as Fisher and Tippett point out if $f$ is Gaussian, the convergence to the Gumbell distribution is quite slow in $n$. I see as a rationale for doing something somewhat different...

Instead, we might imagine that $n$ is a parameter that we can assign a prior and integrate out.  That is, we could write $f_{max}(x)=P(\text{max size is } x)=\sum_{n=0}^\infty P(\text{max size is } x|n)p(n)=\sum_{n=0}^\infty p(n)f_n(x) = \sum_{n=0}^\infty p(n) n F(x)^{n-1} f(x)$

For example, if we assume that $n$ follows a Poisson distribution with mean $\lambda$ then $f_{max}(x)=\sum_{n=0}^\infty \frac{e^{-\lambda}\lambda^n}{n!} n F(x)^{n-1} f(x) = \lambda e^{-\lambda (1-F(x))}f(x)$.  

If we have several, say k, independent samples for which each maximum was recorded, the likelihood is $L=\prod_{i=1}^k f_{max}(x_i)$ which in the Poisson-n case is 

\begin{equation}
L=\lambda^ke^{-\lambda[k-\sum_{i=1}^kF(x_i,\theta)]}\prod_{i=1}^k f(x_i,\theta)
\end{equation}
\noindent where $\theta$ is the collection of parameters governing the underlying size distribution.  For inference, we can use maximum likelihood or Bayes, to obtain estimates of $\theta$ conditional on $\lambda$. If we want to allow for a bit more uncertainty in $n$, we can assign $\lambda$ a prior as well, e.g. Gamma($\alpha, \beta$).  If we assume the sizes are normally distributed, then $\theta$ includes the mean and variance and we probably want $k>3.$
