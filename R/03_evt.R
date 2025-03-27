fit_mod(scenario1_maxima)

evd::fgev(x = scenario1_maxima)


evt_stan <- cmdstan_model("models/evt.stan", stanc_options = list("O1"))

    fit <- evt_stan$sample(
        data = list(k = length(maxima), x = maxima),
        iter_warmup = 1000,
        iter_sampling = 1000,
        chains = 4,
        parallel_chains = 4,
        refresh = 500
    )