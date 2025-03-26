Prompting Copilot inside VScode (using Claude 3.5 sonnet Preview) to generate a my readme.md that will be used as instructions for the Claude model inside cline: 
<general-info>
I am a expert in R programming, I write predominantly using the tidyverse. I want to make clear and concise code, but that can be easily read. I like to keep my files organised. I want to write in British English spelling.
</general-info>

<aim>
I want to develop new methods for estimating the maximum body length of a fish population. I want this estimate to be more statistically robust than just taking the maximum observed length, or, for example, the median of the five largest fish.

I will use two methods: 

1. Extreme value theory (EVT)
2. Numerical estimation (NE)

The EVT method will develop from the idea of extreme value theory, that says that the maxima of a set of samples will follow a generalized extreme value distribution. This distribution has three parameters: location, scale, and shape. The location parameter is the mode of the distribution, the scale parameter is the spread, and the shape parameter is the tail behaviour. The tail behaviour is important because it determines how quickly the probability of observing a very large value decreases as the value increases.

EVT is similar to the central limit theorem, where the means of a set of samples will follow a normal distribution. The normal distribution has two parameters: mean and standard deviation. The mean is the location parameter, and the standard deviation is the scale parameter. The normal distribution has a fixed shape, so there is no shape parameter.

The second method will use a numerical estimation method. This method will try to estimate the parameters of the underlying body size distribution (we assume a normal distribution, truncated positive values only), and the sample size that would give rise to the observed maxima values. The theory is outlined more in the file @munch_notes.md.

I want to be able to give a set of maxima values (say 12 values), which the two methods will estimate the cdf of the maxima values. 
</aim>


<methods>
I will use gev() from the evd package to fit the generalized extreme value distribution to the maxima values. I will use the stan model outlined in @max_est.stan (which is written based on @munch_notes.md) to estimate the parameters of the underlying body size distribution and the sample size.

I want to use a simulation approach to first test the two methods. To see how well they perform estimating the maximum length of simulated data. I want the simulation data to simulate a fishing competition with $k$ fishers, each only reporting their maximum fish. The number of fish caught by each fisher, $n$, will be derived from a Poisson distribution with mean $\lambda$. The body size of the fish will be normally distributed with mean $\mu$ and standard deviation $\sigma$, and truncated to only positive values. 

The parameter combinations will be as follows: 

sim_tbl <-
    expand_grid(
        rep = 1:10,
        k = c(3, 5, 10, 20, 30, 100),
        n_lambda = c(30, 50, 200, 2000, 10000),
        mu = c(20, 50, 100),
    ) |>
    mutate(
        sigma = mu * 0.34,
        min_size = 10,
        filename = paste0(
            "rep", rep,
            "_k", k,
            "_lambda", n_lambda,
            "_mu", mu
        )
    )

I want both methods to use the same simulated data to test each of the approaches. 

</methods>


<output>

I want to produce the four figures in the paper: 
1. Conceptual figure showing the concept of trying to estimate the maximum length of a fish population, and how the methods estimate the distribution of the maxima values.
2. A figure with three panels. Similar to what is shown in the plot saved as "output/figures/simulation_evt.png" in the script @main.qmd. But I want a third panel in the plot that includes the numerical estimation method and well as the EVT method that is currently shown.
3. I want to produce a figure showing how the estimate of maximum length relates to the true maximum length for each of the methods, and for each of the simulation scenarios. Similar to the plot @output/figures/maxima_vs_estmax.png.
4. I want to apply the models to a case-study in this case it will be using maximum lengths derived from snapper. here are the snapper values:
kg2cm <- function(w, a = 0.04478, b = 2.673) ((w * 1000) / a)^(1 / b)

snapper_maxima <- tibble(
    type = c(
        rep("length", 8),
        rep("weight", 4)
    ),
    max = c(
        91.3, 102, 112, 107, 107, 99.2, 95, 82.2,
        kg2cm(c(11.8, 18.4, 16.5, 17.2))
    )
)

</output>


<format>
I want the code to be written in simple modular fashion. With a separate script for each component of the analysis that is easy to follow. I want these scripts to be in a new folder and labelled sequentially. 

I want the figures to be saved in a new folder called "output/figures". I want the figures to be saved in a high-quality format, such as .png. I want the figures to be labelled with a clear and concise title.

I want a qmd file that will run the R scripts and then perform the plotting explaining each step. I want the qmd file to be saved in the root directory.

</format>