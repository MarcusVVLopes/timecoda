# Advanced model engines and the formula interface

Beyond single-level linear and generalized linear models
([`vignette("composition-predictor", package = "timecoda")`](https://marcusvvlopes.github.io/timecoda/articles/composition-predictor.md)),
[`comp_iv()`](https://marcusvvlopes.github.io/timecoda/reference/comp_iv.md)
supports multilevel, survival, ordinal, multinomial, and outlier-robust
designs, along with custom
[`ilr()`](https://marcusvvlopes.github.io/timecoda/reference/ilr.md)
formulas, moderation, varying-total compositions, and engine-specific
diagnostics.

In addition to `zbmi`, `obesity`, and `sex`, `timeuse` includes a
clustering variable (`school`, ten schools), an ordered academic band
(`grade`: `below`, `expected`, `above`), and follow-up time in months to
obesity onset (`time`, `event`).

``` r

library(timecoda)
```

``` r

labs <- c(
  MVPA = "Moderate-to-vigorous PA",
  LPA  = "Light PA",
  SB   = "Sedentary behaviour",
  SPT  = "Sleep"
)

cd <- declare_composition(
  data    = timeuse,
  parts   = c("MVPA", "LPA", "SB", "SPT"),
  closure = 1440,
  units   = "min/day",
  labels  = labs
)
```

Several backends (`coxme`, `ordinal`, `mclogit`, `robustbase`,
`robustlmm`) are suggested packages; examples requiring them run only
when installed.

## Multilevel models

For clustered or repeated-measures data, pass the grouping variable to
`random` with `model = "lmer"` or `"glmer"`:

``` r

fit_lmer <- comp_iv(comp = cd, outcome = "zbmi", model = "lmer",
                    random = "school")
summary(object = fit_lmer)
#> <comp_iv> compositional linear mixed model
#>   Parts: Moderate-to-vigorous PA, Light PA, Sedentary behaviour, Sleep (total = 1440 min/day)
#>   Omnibus composition test (Wald): chi-squared(3) = 343.561, p = <2e-16
#>   Intraclass correlation (school): 0.264
#>   Convergence: OK
#> 
#> Compositional effects (Mean difference)
#> # A tibble: 5 × 6
#>   parameter          label                  estimate conf.low conf.high  p.value
#>   <chr>              <chr>                     <dbl>    <dbl>     <dbl>    <dbl>
#> 1 Global Composition Global Composition       NA       NA       NA      3.70e-74
#> 2 MVPA vs remaining  Moderate-to-vigorous …   -0.691   -0.778   -0.604  2.02e-40
#> 3 LPA vs remaining   Light PA                 -0.200   -0.426    0.0255 8.19e- 2
#> 4 SB vs remaining    Sedentary behaviour       1.59     1.28     1.89   3.52e-21
#> 5 SPT vs remaining   Sleep                    -0.697   -1.00    -0.391  1.08e- 5
```

The summary adds the intraclass correlation (ICC) for each grouping term
(here, 26.4% of total variance lies between schools). When `pbkrtest` is
installed, both [`summary()`](https://rdrr.io/r/base/summary.html) and
[`predict_realloc()`](https://marcusvvlopes.github.io/timecoda/reference/predict_realloc.md)
use Kenward–Roger degrees of freedom (Kenward & Roger, 1997; see
[`?comp_iv`](https://marcusvvlopes.github.io/timecoda/reference/comp_iv.md)
for alternatives):

``` r

predict_realloc(object = fit_lmer, from = "SB", to = "MVPA", delta = 10)
#> # A tibble: 1 × 10
#>   from  to    delta estimate std.error    df conf.low conf.high statistic
#>   <chr> <chr> <dbl>    <dbl>     <dbl> <dbl>    <dbl>     <dbl>     <dbl>
#> 1 SB    MVPA     10   -0.203    0.0119  288.   -0.226    -0.179     -17.0
#> # ℹ 1 more variable: p.value <dbl>
```

For binary or count clustered outcomes, switch to `"glmer"` with
`family`:

``` r

fit_glmer <- comp_iv(comp = cd, outcome = "obesity", model = "glmer",
                     family = "binomial", random = "school")
summary(object = fit_glmer)
#> <comp_iv> compositional generalized linear mixed model (binomial)
#>   Parts: Moderate-to-vigorous PA, Light PA, Sedentary behaviour, Sleep (total = 1440 min/day)
#>   Omnibus composition test (Wald): chi-squared(3) = 28.757, p = 2.52e-06
#>   Latent intraclass correlation (school): 0.001
#>   Convergence: OK
#> 
#> Compositional effects (Odds ratio)
#> # A tibble: 5 × 6
#>   parameter          label                   estimate conf.low conf.high p.value
#>   <chr>              <chr>                      <dbl>    <dbl>     <dbl>   <dbl>
#> 1 Global Composition Global Composition        NA       NA        NA     2.52e-6
#> 2 MVPA vs remaining  Moderate-to-vigorous PA    0.501    0.382     0.656 5.06e-7
#> 3 LPA vs remaining   Light PA                   0.984    0.506     1.91  9.61e-1
#> 4 SB vs remaining    Sedentary behaviour        2.98     1.16      7.68  2.37e-2
#> 5 SPT vs remaining   Sleep                      0.681    0.271     1.71  4.15e-1
```

On logistic mixed models, the between-school share is reported as a
latent ICC using \pi^2 / 3 as the level-1 variance (Goldstein et al.,
2002)—near zero (0.001) here.

## Survival outcomes

For time-to-event data, supply `time` and `event` with `model = "cox"`:

``` r

fit_cox <- comp_iv(comp = cd, model = "cox", time = "time", event = "event",
                   covariates = "sex")
summary(object = fit_cox)
#> <comp_iv> compositional Cox proportional hazards model
#>   Parts: Moderate-to-vigorous PA, Light PA, Sedentary behaviour, Sleep (total = 1440 min/day)
#>   Omnibus composition test (Wald): chi-squared(3) = 81.235, p = <2e-16
#> 
#> Compositional effects (Hazard ratio)
#> # A tibble: 5 × 6
#>   parameter          label                  estimate conf.low conf.high  p.value
#>   <chr>              <chr>                     <dbl>    <dbl>     <dbl>    <dbl>
#> 1 Global Composition Global Composition       NA       NA        NA     1.67e-17
#> 2 MVPA vs remaining  Moderate-to-vigorous …    0.589    0.516     0.672 4.22e-15
#> 3 LPA vs remaining   Light PA                  1.04     0.763     1.41  8.17e- 1
#> 4 SB vs remaining    Sedentary behaviour       3.44     2.17      5.45  1.36e- 7
#> 5 SPT vs remaining   Sleep                     0.476    0.314     0.720 4.45e- 4
```

Both [`summary()`](https://rdrr.io/r/base/summary.html) and
[`predict_realloc()`](https://marcusvvlopes.github.io/timecoda/reference/predict_realloc.md)
report hazard ratios (no convergence line is printed because
[`comp_iv()`](https://marcusvvlopes.github.io/timecoda/reference/comp_iv.md)
does not probe Cox fits for convergence):

``` r

predict_realloc(object = fit_cox, from = "SB", to = "MVPA", delta = c(10, 30))
#> # A tibble: 2 × 10
#>   from  to    delta estimate std.error    df conf.low conf.high statistic
#>   <chr> <chr> <dbl>    <dbl>     <dbl> <dbl>    <dbl>     <dbl>     <dbl>
#> 1 SB    MVPA     10    0.856    0.0158   Inf    0.826     0.888     -8.44
#> 2 SB    MVPA     30    0.683    0.0304   Inf    0.626     0.745     -8.56
#> # ℹ 1 more variable: p.value <dbl>
```

Moving 10 minutes from `SB` to `MVPA` is associated with a hazard of
obesity onset 0.856 times that of the reference day.

To add a school-level random intercept (`(1 | school)`), switch to
`"coxme"`:

``` r

fit_coxme <- comp_iv(comp = cd, model = "coxme", time = "time",
                     event = "event", random = "school")
summary(object = fit_coxme)
#> <comp_iv> compositional mixed-effects Cox proportional hazards model
#>   Parts: Moderate-to-vigorous PA, Light PA, Sedentary behaviour, Sleep (total = 1440 min/day)
#>   Omnibus composition test (Wald): chi-squared(3) = 78.846, p = <2e-16
#>   Median hazard ratio (school): 1.231
#> 
#> Compositional effects (Hazard ratio)
#> # A tibble: 5 × 6
#>   parameter          label                  estimate conf.low conf.high  p.value
#>   <chr>              <chr>                     <dbl>    <dbl>     <dbl>    <dbl>
#> 1 Global Composition Global Composition       NA       NA        NA     5.43e-17
#> 2 MVPA vs remaining  Moderate-to-vigorous …    0.600    0.527     0.683 1.21e-14
#> 3 LPA vs remaining   Light PA                  1.06     0.778     1.44  7.22e- 1
#> 4 SB vs remaining    Sedentary behaviour       3.16     2.01      4.98  7.00e- 7
#> 5 SPT vs remaining   Sleep                     0.499    0.330     0.754 9.80e- 4
```

The extra line reports the median hazard ratio for school, \exp(\sqrt{2
v}\\\Phi^{-1}(0.75)) (Austin et al., 2017)—always \ge 1, with 1
indicating no between-school variation in hazard.

## Ordinal and unordered categorical outcomes

For an ordered factor like `grade`, `model = "clm"` fits a
cumulative-link model via
[`ordinal::clm()`](https://rdrr.io/pkg/ordinal/man/clm.html)
(Christensen, 2026), where odds ratios above 1 indicate higher odds of
being at or above any cut-point:

``` r

fit_clm <- comp_iv(comp = cd, outcome = "grade", model = "clm")
summary(object = fit_clm)
#> <comp_iv> compositional cumulative-link (ordinal) model
#>   Parts: Moderate-to-vigorous PA, Light PA, Sedentary behaviour, Sleep (total = 1440 min/day)
#>   Omnibus composition test (Wald): chi-squared(3) = 23.941, p = 2.57e-05
#>   Convergence: OK
#> 
#> Compositional effects (Cumulative odds ratio)
#> # A tibble: 5 × 6
#>   parameter          label                   estimate conf.low conf.high p.value
#>   <chr>              <chr>                      <dbl>    <dbl>     <dbl>   <dbl>
#> 1 Global Composition Global Composition        NA       NA        NA     2.57e-5
#> 2 MVPA vs remaining  Moderate-to-vigorous PA    1.58     1.27      1.97  3.91e-5
#> 3 LPA vs remaining   Light PA                   1.05     0.627     1.76  8.52e-1
#> 4 SB vs remaining    Sedentary behaviour        0.343    0.167     0.705 3.57e-3
#> 5 SPT vs remaining   Sleep                      1.75     0.855     3.59  1.26e-1
```

Clustered ordinal outcomes use `"clmm"` with `random`, reporting the
same latent ICC as `"glmer"`:

``` r

fit_clmm <- comp_iv(comp = cd, outcome = "grade", model = "clmm",
                    random = "school")
summary(object = fit_clmm)
#> <comp_iv> compositional cumulative-link mixed (ordinal) model
#>   Parts: Moderate-to-vigorous PA, Light PA, Sedentary behaviour, Sleep (total = 1440 min/day)
#>   Omnibus composition test (Wald): chi-squared(3) = 21.830, p = 7.08e-05
#>   Latent intraclass correlation (school): 0.069
#>   Convergence: OK
#> 
#> Compositional effects (Cumulative odds ratio)
#> # A tibble: 5 × 6
#>   parameter          label                   estimate conf.low conf.high p.value
#>   <chr>              <chr>                      <dbl>    <dbl>     <dbl>   <dbl>
#> 1 Global Composition Global Composition        NA       NA        NA     7.08e-5
#> 2 MVPA vs remaining  Moderate-to-vigorous PA    1.58     1.26      1.98  6.19e-5
#> 3 LPA vs remaining   Light PA                   0.981    0.574     1.68  9.44e-1
#> 4 SB vs remaining    Sedentary behaviour        0.386    0.184     0.809 1.17e-2
#> 5 SPT vs remaining   Sleep                      1.67     0.796     3.48  1.76e-1
```

For unordered categories, `"mblogit"` rejects ordered factors, so we
first add an unordered copy of `grade` and re-declare the composition:

``` r

tu <- timeuse
tu$grade_cat <- factor(as.character(tu$grade), levels = levels(tu$grade))
cd_cat <- declare_composition(data = tu, parts = c("MVPA", "LPA", "SB", "SPT"),
                              closure = 1440, units = "min/day",
                              labels = labs)

fit_mb <- comp_iv(comp = cd_cat, outcome = "grade_cat", model = "mblogit")
summary(object = fit_mb)
#> <comp_iv> compositional multinomial baseline-logit model
#>   Parts: Moderate-to-vigorous PA, Light PA, Sedentary behaviour, Sleep (total = 1440 min/day)
#>   Omnibus composition test (Wald): chi-squared(6) = 24.174, p = 0.000485
#>   Convergence: OK
#> 
#> Compositional effects (Relative risk ratio)
#> # A tibble: 9 × 7
#>   parameter          label    estimate conf.low conf.high p.value response_level
#>   <chr>              <chr>       <dbl>    <dbl>     <dbl>   <dbl> <chr>         
#> 1 Global Composition Global …   NA      NA         NA     4.85e-4 NA            
#> 2 MVPA vs remaining  Moderat…    1.57    1.18       2.08  1.76e-3 expected      
#> 3 MVPA vs remaining  Moderat…    1.76    1.31       2.37  1.73e-4 above         
#> 4 LPA vs remaining   Light PA    1.28    0.634      2.59  4.90e-1 expected      
#> 5 LPA vs remaining   Light PA    1.08    0.529      2.19  8.39e-1 above         
#> 6 SB vs remaining    Sedenta…    0.552   0.206      1.48  2.37e-1 expected      
#> 7 SB vs remaining    Sedenta…    0.245   0.0902     0.666 5.85e-3 above         
#> 8 SPT vs remaining   Sleep       0.902   0.344      2.37  8.35e-1 expected      
#> 9 SPT vs remaining   Sleep       2.15    0.795      5.82  1.31e-1 above
```

Each non-reference category (`expected`, `above`) is contrasted with the
first level (`below`) as a relative-risk ratio, giving a 6-df omnibus
test across the two equations and one row per category in reallocations:

``` r

predict_realloc(object = fit_mb, from = "SB", to = "MVPA", delta = 10)
#> # A tibble: 2 × 11
#>   from  to    delta response_level estimate std.error    df conf.low conf.high
#>   <chr> <chr> <dbl> <chr>             <dbl>     <dbl> <dbl>    <dbl>     <dbl>
#> 1 SB    MVPA     10 expected           1.13    0.0441   Inf     1.05      1.22
#> 2 SB    MVPA     10 above              1.18    0.0484   Inf     1.09      1.28
#> # ℹ 2 more variables: statistic <dbl>, p.value <dbl>
```

Clustered multinomial models also accept `random` (either a grouping
name such as `"school"` or a random-slope term like `"(ilr1 | school)"`;
see
[`?comp_iv`](https://marcusvvlopes.github.io/timecoda/reference/comp_iv.md)).

Note that survival (`"cox"`, `"coxme"`), ordinal (`"clm"`, `"clmm"`),
and multinomial (`"mblogit"`) engines have no single response scale for
absolute predictions:
[`predict_realloc()`](https://marcusvvlopes.github.io/timecoda/reference/predict_realloc.md)
and [`plot()`](https://rdrr.io/r/graphics/plot.default.html) give
relative estimates only (pass `category` to
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) for a specific
multinomial level).

## Robust linear models

When continuous outcomes contain influential outliers, `"lmrob"` fits an
MM-estimator (Yohai, 1987) and `"rlmer"` fits a robust linear mixed
model (Koller, 2016):

``` r

fit_lmrob <- comp_iv(comp = cd, outcome = "zbmi", model = "lmrob")
summary(object = fit_lmrob)
#> <comp_iv> compositional robust linear model
#>   Parts: Moderate-to-vigorous PA, Light PA, Sedentary behaviour, Sleep (total = 1440 min/day)
#>   Omnibus composition test (Wald): chi-squared(3) = 309.625, p = <2e-16
#>   Convergence: OK
#> 
#> Compositional effects (Mean difference)
#> # A tibble: 5 × 6
#>   parameter          label                  estimate conf.low conf.high  p.value
#>   <chr>              <chr>                     <dbl>    <dbl>     <dbl>    <dbl>
#> 1 Global Composition Global Composition       NA       NA       NA      8.21e-67
#> 2 MVPA vs remaining  Moderate-to-vigorous …   -0.718   -0.830   -0.607  9.75e-30
#> 3 LPA vs remaining   Light PA                 -0.223   -0.505    0.0580 1.19e- 1
#> 4 SB vs remaining    Sedentary behaviour       1.74     1.39     2.09   6.19e-20
#> 5 SPT vs remaining   Sleep                    -0.799   -1.13    -0.466  3.54e- 6
```

``` r

fit_rlmer <- comp_iv(comp = cd, outcome = "zbmi", model = "rlmer",
                     random = "school")
summary(object = fit_rlmer)
#> <comp_iv> compositional robust linear mixed model
#>   Parts: Moderate-to-vigorous PA, Light PA, Sedentary behaviour, Sleep (total = 1440 min/day)
#>   Omnibus composition test (Wald): chi-squared(3) = 324.018, p = <2e-16
#>   Intraclass correlation (school): 0.313
#>   Convergence: OK
#> 
#> Compositional effects (Mean difference)
#> # A tibble: 5 × 6
#>   parameter          label                  estimate conf.low conf.high  p.value
#>   <chr>              <chr>                     <dbl>    <dbl>     <dbl>    <dbl>
#> 1 Global Composition Global Composition       NA       NA      NA       6.29e-70
#> 2 MVPA vs remaining  Moderate-to-vigorous …   -0.702   -0.792  -0.612   1.33e-52
#> 3 LPA vs remaining   Light PA                 -0.230   -0.464   0.00392 5.40e- 2
#> 4 SB vs remaining    Sedentary behaviour       1.61     1.29    1.93    2.91e-23
#> 5 SPT vs remaining   Sleep                    -0.675   -0.994  -0.357   3.18e- 5
```

Because `robustlmm` does not estimate residual degrees of freedom,
[`summary()`](https://rdrr.io/r/base/summary.html) and
[`predict_realloc()`](https://marcusvvlopes.github.io/timecoda/reference/predict_realloc.md)
compute confidence intervals and p-values from a standard normal
approximation (`df = Inf`) rather than a t-distribution (for example,
`df = 288` for `"lmer"`):

``` r

predict_realloc(object = fit_rlmer, from = "SB", to = "MVPA", delta = 10)
#> # A tibble: 1 × 10
#>   from  to    delta estimate std.error    df conf.low conf.high statistic
#>   <chr> <chr> <dbl>    <dbl>     <dbl> <dbl>    <dbl>     <dbl>     <dbl>
#> 1 SB    MVPA     10   -0.206    0.0124   Inf   -0.230    -0.182     -16.6
#> # ℹ 1 more variable: p.value <dbl>
```

## The formula interface and `ilr()`

Instead of passing `outcome`, `covariates`, `random`, `time`, `event`,
or `moderator` separately, you can supply a standard R `formula` using
[`ilr()`](https://marcusvvlopes.github.io/timecoda/reference/ilr.md) as
a placeholder for the composition:

``` r

fit_main <- comp_iv(comp = cd, formula = zbmi ~ ilr() + sex, model = "lm")
fit_main
#> <timecoda_fit> compositional linear model
#>   Parts: MVPA, LPA, SB, SPT (total = 1440 min/day)
#>   Reference composition: MVPA = 28.1, LPA = 235.5, SB = 677.2, SPT = 499.2
#> 
#> 
#> Call:
#> stats::lm(formula = zbmi ~ (ilr1 + ilr2 + ilr3) + sex, data = data)
#> 
#> Coefficients:
#> (Intercept)         ilr1         ilr2         ilr3      sexmale  
#>    -2.34267     -0.70652     -0.47025      1.54824     -0.06241
```

All other terms pass through to the backend unchanged, while
[`ilr()`](https://marcusvvlopes.github.io/timecoda/reference/ilr.md)
expands by two rules
([`?ilr`](https://marcusvvlopes.github.io/timecoda/reference/ilr.md)).
Inside formula operators (`+`, `-`, `*`, `:`, `/`, `^`, `%in%`, `|`,
`||`) and parentheses,
[`ilr()`](https://marcusvvlopes.github.io/timecoda/reference/ilr.md)
becomes `(ilr1 + ... + ilr(D-1))`, so `(ilr() | school)` adds random
slopes on all coordinates and `ilr() * sex` creates the full
composition-by-sex interaction:

``` r

fit_int <- comp_iv(comp = cd, formula = zbmi ~ ilr() * sex, model = "lm")
fit_int
#> <timecoda_fit> compositional linear model
#>   Parts: MVPA, LPA, SB, SPT (total = 1440 min/day)
#>   Reference composition: MVPA = 28.1, LPA = 235.5, SB = 677.2, SPT = 499.2
#> 
#> 
#> Call:
#> stats::lm(formula = zbmi ~ (ilr1 + ilr2 + ilr3) * sex, data = data)
#> 
#> Coefficients:
#>  (Intercept)          ilr1          ilr2          ilr3       sexmale  
#>      -2.2024       -0.6580       -0.3483        1.8404       -0.4386  
#> ilr1:sexmale  ilr2:sexmale  ilr3:sexmale  
#>      -0.1255       -0.3135       -0.7065
```

Inside any other function call,
[`ilr()`](https://marcusvvlopes.github.io/timecoda/reference/ilr.md)
distributes across coordinates—so `I(ilr()^2)` expands to
`(I(ilr1^2) + I(ilr2^2) + I(ilr3^2))`, and
[`poly()`](https://rdrr.io/r/stats/poly.html) or
[`splines::ns()`](https://rdrr.io/r/splines/ns.html) distribute the same
way. Combining both rules gives the full quadratic composition model:

``` r

fit_quad <- comp_iv(comp = cd, formula = zbmi ~ ilr()^2 + I(ilr()^2) + sex,
                    model = "lm")
fit_quad
#> <timecoda_fit> compositional linear model
#>   Parts: MVPA, LPA, SB, SPT (total = 1440 min/day)
#>   Reference composition: MVPA = 28.1, LPA = 235.5, SB = 677.2, SPT = 499.2
#> 
#> 
#> Call:
#> stats::lm(formula = zbmi ~ (ilr1 + ilr2 + ilr3)^2 + (I(ilr1^2) + 
#>     I(ilr2^2) + I(ilr3^2)) + sex, data = data)
#> 
#> Coefficients:
#> (Intercept)         ilr1         ilr2         ilr3    I(ilr1^2)    I(ilr2^2)  
#>    -2.42867     -0.98660      0.04431      1.59720     -0.04117      0.53721  
#>   I(ilr3^2)      sexmale    ilr1:ilr2    ilr1:ilr3    ilr2:ilr3  
#>    -1.06062     -0.07917     -0.15457     -0.22287      0.18520
```

If a formula includes
[`ilr()`](https://marcusvvlopes.github.io/timecoda/reference/ilr.md)
without any coordinate appearing as a standalone term or first argument
to [`poly()`](https://rdrr.io/r/stats/poly.html), `ns()`, or `bs()`,
[`comp_iv()`](https://marcusvvlopes.github.io/timecoda/reference/comp_iv.md)
warns that estimates may be biased.

## Moderation and non-linear composition terms

To test whether the compositional association differs across a
covariate, naming it in `moderator` fits the same interaction as
`zbmi ~ ilr() * sex` and adds a joint Wald moderation test plus per-part
interaction columns:

``` r

fit_mod <- comp_iv(comp = cd, outcome = "zbmi", model = "lm",
                   moderator = "sex")
summary(object = fit_mod)
#> <comp_iv> compositional linear model
#>   Parts: Moderate-to-vigorous PA, Light PA, Sedentary behaviour, Sleep (total = 1440 min/day)
#>   Omnibus composition test (Wald): chi-squared(3) = 156.595, p = <2e-16
#>   Moderation by sex (Wald): chi-squared(3) = 5.209, p = 0.157
#> 
#> Compositional effects (Mean difference)
#> # A tibble: 5 × 9
#>   parameter      label estimate conf.low conf.high  p.value mod_statistic mod_df
#>   <chr>          <chr>    <dbl>    <dbl>     <dbl>    <dbl>         <dbl>  <int>
#> 1 Global Compos… Glob…   NA       NA        NA     9.95e-34        NA         NA
#> 2 MVPA vs remai… Mode…   -0.658   -0.790    -0.527 6.12e-20         1.53       1
#> 3 LPA vs remain… Ligh…   -0.109   -0.473     0.255 5.56e- 1         0.952      1
#> 4 SB vs remaini… Sede…    1.89     1.40      2.38  5.02e-13         1.21       1
#> 5 SPT vs remain… Sleep   -1.12    -1.58     -0.663 2.19e- 6         4.59       1
#> # ℹ 1 more variable: mod_p.value <dbl>
```

Here, the per-part test for sleep yields p = 0.0322 (unadjusted across
the four parts). On a `moderator` fit,
[`predict_realloc()`](https://marcusvvlopes.github.io/timecoda/reference/predict_realloc.md)
automatically splits results by a factor moderator unless pinned with
`at`; on a raw formula fit like `fit_int`, you must supply `by` or `at`
explicitly (see
[`vignette("predicting-reallocations", package = "timecoda")`](https://marcusvvlopes.github.io/timecoda/articles/predicting-reallocations.md)).

For a quadratic fit like `fit_quad`,
[`summary()`](https://rdrr.io/r/base/summary.html) reports a 9-df
omnibus test and a joint curvature Wald test across all quadratic and
cross-product terms. Because no single `ilr1` slope exists when
curvature is present, each per-part row gives the local slope at the
reference composition via
[`emmeans::emtrends()`](https://rvlenth.github.io/emmeans/reference/emtrends.html):

``` r

summary(object = fit_quad)
#> <comp_iv> compositional linear model
#>   Parts: Moderate-to-vigorous PA, Light PA, Sedentary behaviour, Sleep (total = 1440 min/day)
#>   Omnibus composition test (Wald): chi-squared(9) = 308.015, p = <2e-16
#>   Curvature (Wald): chi-squared(6) = 12.629, p = 0.0493
#> 
#> Compositional effects (Mean difference)
#> # A tibble: 5 × 6
#>   parameter          label                  estimate conf.low conf.high  p.value
#>   <chr>              <chr>                     <dbl>    <dbl>     <dbl>    <dbl>
#> 1 Global Composition Global Composition      NA        NA        NA     5.20e-61
#> 2 MVPA vs remaining  Moderate-to-vigorous …  -0.726    -0.835    -0.617 5.52e-31
#> 3 LPA vs remaining   Light PA                -0.0818   -0.351     0.188 5.51e- 1
#> 4 SB vs remaining    Sedentary behaviour      1.65      1.30      2.01  1.51e-17
#> 5 SPT vs remaining   Sleep                   -0.845    -1.20     -0.487 5.28e- 6
#> Per-part effects are slopes at the reference composition, from emmeans::emtrends(), as the
#> model's composition terms are non-linear or carry no main effect.
```

Reallocations and surface plots work identically on quadratic
fits—notice below that the 30-minute `zbmi` difference is less than
three times the 10-minute difference:

``` r

predict_realloc(object = fit_quad, from = "SB", to = "MVPA", delta = c(10, 30))
#> # A tibble: 2 × 10
#>   from  to    delta estimate std.error    df conf.low conf.high statistic
#>   <chr> <chr> <dbl>    <dbl>     <dbl> <dbl>    <dbl>     <dbl>     <dbl>
#> 1 SB    MVPA     10   -0.215    0.0160   289   -0.247    -0.184     -13.5
#> 2 SB    MVPA     30   -0.537    0.0429   289   -0.621    -0.452     -12.5
#> # ℹ 1 more variable: p.value <dbl>
```

``` r

plot(x = fit_quad)
```

![One of four ternary panels, one for each group of three parts, of the
quadratic linear model's surface of the mean difference in z-BMI from
the reference composition as tiles. Several dozen small hexagonal tiles
form a patch inside the triangle, coloured from dark purple for the
panel's lowest mean difference to yellow for its highest, with a black
point at the reference composition; the rest of the triangle is blank
grey, and a legend titled Mean difference gives the
scale.](advanced-models_files/figure-html/quadratic-plot-1.png)![One of
four ternary panels, one for each group of three parts, of the quadratic
linear model's surface of the mean difference in z-BMI from the
reference composition as tiles. Several dozen small hexagonal tiles form
a patch inside the triangle, coloured from dark purple for the panel's
lowest mean difference to yellow for its highest, with a black point at
the reference composition; the rest of the triangle is blank grey, and a
legend titled Mean difference gives the
scale.](advanced-models_files/figure-html/quadratic-plot-2.png)![One of
four ternary panels, one for each group of three parts, of the quadratic
linear model's surface of the mean difference in z-BMI from the
reference composition as tiles. Several dozen small hexagonal tiles form
a patch inside the triangle, coloured from dark purple for the panel's
lowest mean difference to yellow for its highest, with a black point at
the reference composition; the rest of the triangle is blank grey, and a
legend titled Mean difference gives the
scale.](advanced-models_files/figure-html/quadratic-plot-3.png)![One of
four ternary panels, one for each group of three parts, of the quadratic
linear model's surface of the mean difference in z-BMI from the
reference composition as tiles. Several dozen small hexagonal tiles form
a patch inside the triangle, coloured from dark purple for the panel's
lowest mean difference to yellow for its highest, with a black point at
the reference composition; the rest of the triangle is blank grey, and a
legend titled Mean difference gives the
scale.](advanced-models_files/figure-html/quadratic-plot-4.png)

The full quadratic (`ilr()^2 + I(ilr()^2)`, or `curvature = TRUE`) is
invariant to part order. In contrast, squares-only
(`ilr() + I(ilr()^2)`) and `poly(ilr(), 2)` depend on the order in which
parts were declared (see
[`?ilr`](https://marcusvvlopes.github.io/timecoda/reference/ilr.md)).

## Varying-total compositions

When a composition covers only part of the day—such as waking behaviours
(`MVPA`, `LPA`, `SB`)—row totals vary across participants. Declaring the
composition with `closure = "mean"` keeps it on the average waking-time
scale, and `include_total = TRUE` adds the total covariate
t\_{\mathrm{total}} = \sqrt{D}\ln(g) of Coenders et al. (2017), where g
is each person’s unclosed geometric mean:

``` r

cd_wake <- declare_composition(
  data    = timeuse,
  parts   = c("MVPA", "LPA", "SB"),
  closure = "mean",
  units   = "min/day",
  labels  = labs[c("MVPA", "LPA", "SB")]
)
summary(object = cd_wake, desc = FALSE, var_matrix = FALSE)
#> <compD> compositional variable summary
#>   Parts: Moderate-to-vigorous PA, Light PA, Sedentary behaviour
#>   Closure: 949.1823 min/day
#>   Observations: 300
#> 
#>   Row totals (before closure): mean 949.2 (SD 106.2), range 624.8-1200.1
#>   No zeroes detected.
```

``` r

fit_tot <- comp_iv(comp = cd_wake, outcome = "zbmi", model = "lm",
                   include_total = TRUE)
summary(object = fit_tot)
#> <comp_iv> compositional linear model
#>   Parts: Moderate-to-vigorous PA, Light PA, Sedentary behaviour (total = 949.1823 min/day)
#>   Omnibus composition test (Wald): chi-squared(2) = 94.057, p = <2e-16
#>   Total (log geometric mean, Wald): beta = 1.116, p = <0.001
#> 
#> Compositional effects (Mean difference)
#> # A tibble: 4 × 6
#>   parameter          label                  estimate conf.low conf.high  p.value
#>   <chr>              <chr>                     <dbl>    <dbl>     <dbl>    <dbl>
#> 1 Global Composition Global Composition       NA       NA        NA     3.76e-21
#> 2 MVPA vs remaining  Moderate-to-vigorous …   -1.49    -1.79     -1.18  7.75e-19
#> 3 LPA vs remaining   Light PA                 -0.628   -0.898    -0.358 6.98e- 6
#> 4 SB vs remaining    Sedentary behaviour       2.11     1.65      2.58  3.28e-17
```

The ILR rows now reflect waking-time shares at constant waking duration,
while the positive `t_total` coefficient indicates that a longer waking
day (less sleep, holding waking proportions fixed) is associated with
higher `zbmi`. During reallocations, `t_total` is held at its sample
mean:

``` r

predict_realloc(object = fit_tot, from = "SB", to = "MVPA", delta = 10)
#> # A tibble: 1 × 10
#>   from  to    delta estimate std.error    df conf.low conf.high statistic
#>   <chr> <chr> <dbl>    <dbl>     <dbl> <dbl>    <dbl>     <dbl>     <dbl>
#> 1 SB    MVPA     10   -0.392    0.0410   296   -0.473    -0.312     -9.56
#> # ℹ 1 more variable: p.value <dbl>
```

Using `include_total = TRUE` requires varying row totals; on a
fixed-total composition like `cd`,
[`comp_iv()`](https://marcusvvlopes.github.io/timecoda/reference/comp_iv.md)
stops with an error.

## Engine-specific diagnostics

Each model class uses diagnostic checks suited to its likelihood and
residual structure. `"lm"`, `"glm"`, `"lmer"`, and `"glmer"` fits use
`DHARMa` simulated quantile residuals (Hartig, 2026; see
[`vignette("composition-predictor", package = "timecoda")`](https://marcusvvlopes.github.io/timecoda/articles/composition-predictor.md)).

For `"cox"` fits,
[`fit_check()`](https://marcusvvlopes.github.io/timecoda/reference/fit_check.md)
runs
[`survival::cox.zph()`](https://rdrr.io/pkg/survival/man/cox.zph.html)
to test proportional hazards via scaled Schoenfeld residuals (Grambsch &
Therneau, 1994); pass `plot = TRUE` to inspect the residual curves over
time:

``` r

fc_cox <- fit_check(object = fit_cox, plot = FALSE)
fc_cox
#> <fit_check> Cox proportional hazards model 
#>   Backend: coxzph 
#> 
#> # A tibble: 5 × 4
#>   check                        statistic    df p.value
#>   <chr>                            <dbl> <dbl>   <dbl>
#> 1 Proportional hazards: ilr1       0.504     1   0.478
#> 2 Proportional hazards: ilr2       0.296     1   0.587
#> 3 Proportional hazards: ilr3       0.114     1   0.735
#> 4 Proportional hazards: sex        0.141     1   0.707
#> 5 Proportional hazards: GLOBAL     1.03      4   0.906
```

For `"clm"` fits, it combines
[`ordinal::nominal_test()`](https://rdrr.io/pkg/ordinal/man/nominal.test.html)
(Christensen, 2026) with `sure` surrogate residuals under a fixed seed
(Greenwell et al., 2018; Liu & Zhang, 2018):

``` r

fc_clm <- fit_check(object = fit_clm, plot = FALSE)
fc_clm
#> <fit_check> cumulative-link (ordinal) model 
#>   Backend: ordinal 
#>   Convergence: converged (0 backend warnings)
#>   Note: Proportional-odds rows: likelihood-ratio test of each term against a model letting its effect vary across the outcome thresholds; a small p-value is evidence that the proportional-odds assumption fails for that term. 
#> 
#> # A tibble: 5 × 4
#>   check                   statistic    df p.value
#>   <chr>                       <dbl> <dbl>   <dbl>
#> 1 Proportional odds: ilr1     1.26      1   0.262
#> 2 Proportional odds: ilr2     0.858     1   0.354
#> 3 Proportional odds: ilr3     0.294     1   0.588
#> 4 AIC                       644.       NA  NA    
#> 5 log-likelihood           -317.       NA  NA
```

For `"lmrob"` and `"rlmer"` fits, it summarizes the robustness weights
(here, three of 300 observations are downweighted below 0.5):

``` r

fc_lmrob <- fit_check(object = fit_lmrob, plot = FALSE)
fc_lmrob
#> <fit_check> robust linear model 
#>   Backend: robust 
#>   Convergence: converged (0 backend warnings)
#> 
#> # A tibble: 3 × 4
#>   check                             statistic    df p.value
#>   <chr>                                 <dbl> <dbl>   <dbl>
#> 1 Downweighted obs (weight < 0.5)      0.01      NA  NA    
#> 2 Minimum robustness weight            0.0825    NA  NA    
#> 3 Residual normality (Shapiro-Wilk)    0.994     NA   0.340
```

For `"clmm"`, `"mblogit"`, and `"coxme"` fits, no residual diagnostic is
run;
[`fit_check()`](https://marcusvvlopes.github.io/timecoda/reference/fit_check.md)
reports AIC and log-likelihood only (for `"coxme"`, the penalized
log-likelihood and effective degrees of freedom).

## Engine reference

All eleven engines share the same
[`comp_iv()`](https://marcusvvlopes.github.io/timecoda/reference/comp_iv.md)
interface. `family` is required for `"glm"` and `"glmer"` (and rejected
elsewhere); `random` is required for `"lmer"`, `"glmer"`, `"coxme"`,
`"clmm"`, and `"rlmer"`, optional for `"mblogit"`, and rejected by
single-level engines.

| `model` | Backend | Use case | Effect measure |
|----|----|----|----|
| `"lm"` | [`stats::lm()`](https://rdrr.io/r/stats/lm.html) | Continuous outcome, independent observations | Mean difference |
| `"glm"` | [`stats::glm()`](https://rdrr.io/r/stats/glm.html) | GLM family, independent observations | Odds, rate, risk, or mean ratio; mean difference |
| `"cox"` | [`survival::coxph()`](https://rdrr.io/pkg/survival/man/coxph.html) | Time-to-event, independent observations | Hazard ratio (relative only) |
| `"coxme"` | [`coxme::coxme()`](https://rdrr.io/pkg/coxme/man/coxme.html) | Time-to-event, clustered or repeated | Hazard ratio (relative only) |
| `"lmer"` | [`lme4::lmer()`](https://rdrr.io/pkg/lme4/man/lmer.html) | Continuous outcome, clustered or repeated | Mean difference |
| `"glmer"` | [`lme4::glmer()`](https://rdrr.io/pkg/lme4/man/glmer.html) | GLM family, clustered or repeated | Set by link, as for `"glm"` |
| `"clm"` | [`ordinal::clm()`](https://rdrr.io/pkg/ordinal/man/clm.html) | Ordered-factor outcome | Cumulative odds ratio (relative only) |
| `"clmm"` | [`ordinal::clmm()`](https://rdrr.io/pkg/ordinal/man/clmm.html) | Ordered-factor outcome, clustered or repeated | Cumulative odds ratio (relative only) |
| `"mblogit"` | [`mclogit::mblogit()`](https://melff.github.io/mclogit/reference/mblogit.html) | Unordered-factor outcome (3+ levels), optional `random` | Relative-risk ratio per category (relative only) |
| `"lmrob"` | [`robustbase::lmrob()`](https://rdrr.io/pkg/robustbase/man/lmrob.html) | Continuous outcome, outlier-robust MM-estimation | Mean difference |
| `"rlmer"` | [`robustlmm::rlmer()`](https://rdrr.io/pkg/robustlmm/man/rlmer.html) | Continuous outcome, outlier-robust and clustered | Mean difference |

## References

Austin, P. C., Wagner, P., & Merlo, J. (2017). The median hazard ratio:
A useful measure of variance and general contextual effects in
multilevel survival analysis. *Statistics in Medicine*, 36, 928-938.

Christensen, R. H. B. (2026). *ordinal: Regression models for ordinal
data*. R package version 2026.7-26.

Coenders, G., Martín-Fernández, J. A., & Ferrer-Rosell, B. (2017). When
relative and absolute information matter: Compositional predictor with a
total in generalized linear models. *Statistical Modelling*, 17,
494-512.

Goldstein, H., Browne, W., & Rasbash, J. (2002). Partitioning variation
in multilevel models. *Understanding Statistics*, 1, 223-231.

Grambsch, P. M., & Therneau, T. M. (1994). Proportional hazards tests
and diagnostics based on weighted residuals. *Biometrika*, 81, 515-526.

Greenwell, B. M., McCarthy, A. J., Boehmke, B. C., & Liu, D. (2018).
Residuals and diagnostics for binary and ordinal regression models: An
introduction to the sure package. *The R Journal*, 10(1).

Hartig, F. (2026). *DHARMa: Residual diagnostics for hierarchical
(multi-level / mixed) regression models*. R package version 0.5.0.

Kenward, M. G., & Roger, J. H. (1997). Small sample inference for fixed
effects from restricted maximum likelihood. *Biometrics*, 53, 983-997.

Koller, M. (2016). robustlmm: An R package for robust estimation of
linear mixed-effects models. *Journal of Statistical Software*, 75(6),
1-24.

Liu, D., & Zhang, H. (2018). Residuals and diagnostics for ordinal
regression models: A surrogate approach. *Journal of the American
Statistical Association*, 113, 845-854.

Yohai, V. J. (1987). High breakdown-point and high efficiency robust
estimates for regression. *The Annals of Statistics*, 15, 642-656.

## See also

- [`vignette("composition-predictor", package = "timecoda")`](https://marcusvvlopes.github.io/timecoda/articles/composition-predictor.md):
  linear and logistic fits, summary interpretation, and fitted surfaces.
- [`vignette("predicting-reallocations", package = "timecoda")`](https://marcusvvlopes.github.io/timecoda/articles/predicting-reallocations.md):
  reallocation estimates, schemes, reference days, publication tables,
  and transfer curves.
- [`vignette("optimal-composition", package = "timecoda")`](https://marcusvvlopes.github.io/timecoda/articles/optimal-composition.md):
  searching for the Goldilocks zone.
- [`vignette("timecoda", package = "timecoda")`](https://marcusvvlopes.github.io/timecoda/articles/timecoda.md):
  the whole workflow in one analysis.
