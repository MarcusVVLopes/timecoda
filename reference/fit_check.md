# Check the goodness of fit of a composition model

Runs residual and goodness-of-fit diagnostics on a fitted
[`comp_iv()`](https://marcusvvlopes.github.io/timecoda/reference/comp_iv.md)
model, with the diagnostic backend chosen from the model class. Returns
a `fit_check` list holding a tidy tibble of test statistics and the raw
backend object. With `plot = TRUE` it also draws the backend's
diagnostic plot and prints the report.

## Usage

``` r
fit_check(object, plot = TRUE, ...)

# S3 method for class 'comp_iv'
fit_check(object, plot = TRUE, ...)
```

## Arguments

- object:

  A `comp_iv` object.

- plot:

  Logical; if `TRUE` (default), draw the backend's diagnostic plot(s)
  and print the report. The tidy results are returned either way. The
  `clmm`, `mblogit` and `coxme` fits draw nothing.

- ...:

  Passed to
  [`DHARMa::simulateResiduals()`](https://rdrr.io/pkg/DHARMa/man/simulateResiduals.html)
  for the DHARMa backend (e.g., `n`, the number of simulations, or
  `seed`); ignored by the other backends. DHARMa simulates under its own
  seed (123 by default) and then restores the session's random state.
  The results are therefore the same in a fresh session and after any
  [`set.seed()`](https://rdrr.io/r/base/Random.html); only the arguments
  passed here, such as `seed` or `n`, change them.

## Value

Invisibly, a `fit_check` object: a list with

- `tests`:

  a tibble with one row per check and columns `check`, `statistic`, `df`
  and `p.value` (`NA` where a check has none).

- `raw`:

  the backend object: the DHARMa simulation, the
  [`survival::cox.zph()`](https://rdrr.io/pkg/survival/man/cox.zph.html)
  result, a list of the surrogate residuals and the
  [`ordinal::nominal_test()`](https://rdrr.io/pkg/ordinal/man/nominal.test.html)
  table (`clm`), or the robustness weights (`lmrob`, `rlmer`). It is
  `NULL` for `clmm`, `mblogit` and `coxme` fits. The diagnostic plot is
  drawn as a side effect of `plot = TRUE` and is not returned;
  `plot(x$raw)` re-draws it for the DHARMa and `cox.zph` backends only.

- `note`:

  an explanatory note, or `NULL`.

- `backend`:

  the backend used: `"dharma"`, `"coxzph"`, `"ordinal"`, `"robust"`,
  `"robust_mixed"` or `"limited"`.

- `model_type`:

  the model class.

- `convergence`:

  the fit's convergence state as
  [`comp_iv()`](https://marcusvvlopes.github.io/timecoda/reference/comp_iv.md)
  stores it: a one-row tibble with `converged` (`TRUE`, `FALSE`, or `NA`
  when the backend reports none), `singular` and `n_warnings` (the
  number of backend warnings captured at fit time).

## Details

The backend follows the model class:

- `lm`, `glm`, `lmer`, `glmer`:

  DHARMa simulation-based quantile residuals (KS uniformity, dispersion,
  and outlier tests).

- `cox`:

  [`survival::cox.zph()`](https://rdrr.io/pkg/survival/man/cox.zph.html)
  (proportional-hazards test and Schoenfeld residual plots).

- ordinal `clm`:

  sure surrogate residuals, drawn by sure's own plot method as a Q-Q
  plot and a residuals-against-fitted plot, plus
  [`ordinal::nominal_test()`](https://rdrr.io/pkg/ordinal/man/nominal.test.html)
  (a likelihood-ratio test of proportional odds per term; Christensen,
  ordinal). The surrogate residuals are random draws made under a fixed
  seed (123) with the session's random state restored, so they repeat
  from call to call.

- robust `lmrob` / `rlmer`:

  the robustness weights of the MM-estimate (Yohai 1987) or of the
  robust mixed fit (Koller 2016): the share of observations downweighted
  below 0.5 (half the full weight) and the smallest weight. A
  Shapiro-Wilk test of the residuals follows.

- mixed ordinal `clmm`, multinomial `mblogit`, mixed-effects Cox
  `coxme`:

  fit statistics (AIC and log-likelihood) and a note only. No residual
  diagnostic is run for these classes, nor a proportional-odds test for
  `clmm` or a proportional-hazards test for `coxme`, and nothing is
  drawn. For `coxme` the log-likelihood is the penalized one that
  [`stats::logLik()`](https://rdrr.io/r/stats/logLik.html) returns, with
  the fit's effective degrees of freedom in the AIC. For an `mblogit`
  fit with a random effect, the convergence line shows mclogit's own
  flag, which can read as converged although the inner iterations
  warned. The log-likelihood row is `NA`.

The header of the printed report carries one line on the fit's
convergence (`Convergence: converged` or `Convergence: NOT converged`)
with the number of backend warnings captured at fit time. The line is
left out when the backend reports no convergence state (`lm`, `cox` and
`coxme` fits). The same state is returned in the `convergence` element
and has no row in `tests`.

DHARMa, sure (with tinyplot for its plots) and robustlmm are optional
(Suggests) and checked at call time.

## References

Christensen, R. H. B. ordinal: Regression models for ordinal data. R
package. <https://CRAN.R-project.org/package=ordinal>
[doi:10.32614/CRAN.package.ordinal](https://doi.org/10.32614/CRAN.package.ordinal)

Dunn, P. K., & Smyth, G. K. (1996). Randomized quantile residuals.
*Journal of Computational and Graphical Statistics*, 5(3), 236-244.
[doi:10.1080/10618600.1996.10474708](https://doi.org/10.1080/10618600.1996.10474708)

Grambsch, P. M., & Therneau, T. M. (1994). Proportional hazards tests
and diagnostics based on weighted residuals. *Biometrika*, 81(3),
515-526.
[doi:10.1093/biomet/81.3.515](https://doi.org/10.1093/biomet/81.3.515)

Greenwell, B. M., McCarthy, A. J., Boehmke, B. C., & Liu, D. (2018).
Residuals and diagnostics for binary and ordinal regression models: An
introduction to the sure package. *The R Journal*, 10(1).
[doi:10.32614/RJ-2018-004](https://doi.org/10.32614/RJ-2018-004)

Hartig, F. DHARMa: Residual diagnostics for hierarchical (multi-level /
mixed) regression models. R package.
<https://CRAN.R-project.org/package=DHARMa>
[doi:10.32614/CRAN.package.DHARMa](https://doi.org/10.32614/CRAN.package.DHARMa)

Koller, M. (2016). robustlmm: An R package for robust estimation of
linear mixed-effects models. *Journal of Statistical Software*, 75(6),
1-24. [doi:10.18637/jss.v075.i06](https://doi.org/10.18637/jss.v075.i06)

Liu, D., & Zhang, H. (2018). Residuals and diagnostics for ordinal
regression models: A surrogate approach. *Journal of the American
Statistical Association*, 113(522), 845-854.
[doi:10.1080/01621459.2017.1292915](https://doi.org/10.1080/01621459.2017.1292915)

Yohai, V. J. (1987). High breakdown-point and high efficiency robust
estimates for regression. *The Annals of Statistics*, 15(2), 642-656.
[doi:10.1214/aos/1176350366](https://doi.org/10.1214/aos/1176350366)

## See also

[`comp_iv()`](https://marcusvvlopes.github.io/timecoda/reference/comp_iv.md),
the model object this diagnoses.

## Examples

``` r
cd <- declare_composition(timeuse, c("MVPA", "LPA", "SB", "SPT"),
                          closure = 1440)
m <- comp_iv(cd, outcome = "obesity", model = "glm", family = "binomial")
if (requireNamespace("DHARMa", quietly = TRUE)) {
  fc <- fit_check(m, plot = FALSE)
  fc$tests
}
#> # A tibble: 3 × 4
#>   check         statistic    df p.value
#>   <chr>             <dbl> <dbl>   <dbl>
#> 1 KS uniformity    0.0381    NA   0.777
#> 2 Dispersion       1.01      NA   0.976
#> 3 Outliers         0         NA   1    
```
