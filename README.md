
## Uncertainty in AZ results

AstraZeneca [released some
details](https://www.astrazeneca.com/media-centre/press-releases/2020/azd1222hlr.html)
of results from ongoing clinical trials of their COVID-19 vaccine. Below
are calculations of what 95% confidence intervals might look like for
these results based on the anticipated case splits given the reported
VE. I wrote about it [on
Twitter](https://twitter.com/biosbenk/status/1331673833934368771?s=20).

**tldr**: UK trial = 90.0% (67.3%, 96.9%), Brazil = 62.0% (40.9%, 75.5%)

``` r
uk_trial_param <- list(n = 2741, 
                       n_cases_vax = 3,
                       n_cases_placebo = 30)
brazil_trial_param <- list(n = 8895, 
                           n_cases_vax = 27,
                           n_cases_placebo = 71)
library(sandwich)

get_ve_ci <- function(param){
    n <- param$n 
    n_vax <- round(n / 2)
    n_placebo <- n - n_vax
    n_cases_vax <- param$n_cases_vax
    n_cases_placebo <- param$n_cases_placebo

    # a simple data set
    y_vax <- rep(0, n_vax)
    y_vax[seq_len(n_cases_vax)] <- 1
    y_placebo <- rep(0, n_placebo)
    y_placebo[seq_len(n_cases_placebo)] <- 1
    y <- c(y_vax, y_placebo)
    x <- c(rep(1, n_vax), rep(0, n_placebo))

    # fit regression and get results based on sandwich covariance matrix
    fit <- glm(y ~ x, family = poisson())
    beta_hat <- fit$coef[2]
    se_beta_hat <- sqrt(sandwich::vcovHC(fit)[2,2])
    ve <- 1 - exp(beta_hat)
    ve_cil <- 1 - exp(beta_hat + 1.96 * se_beta_hat)
    ve_ciu <- 1 - exp(beta_hat - 1.96 * se_beta_hat)
    out <- paste0(
      sprintf("%2.1f", ve * 100), " (", 
      sprintf("%2.1f", ve_cil * 100), ", ",
      sprintf("%2.1f", ve_ciu * 100), ")"
    )
    return(out)
}

get_ve_ci(uk_trial_param)
```

    ## [1] "90.0 (67.3, 96.9)"

``` r
get_ve_ci(brazil_trial_param)
```

    ## [1] "62.0 (40.9, 75.5)"
