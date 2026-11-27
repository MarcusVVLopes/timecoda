# Advanced model engines and the formula interface

For the model engines of
[`comp_iv()`](https://marcusvvlopes.github.io/timecoda/reference/comp_iv.md)
beyond `"lm"` and `"glm"`, the formula interface and the diagnostics of
each engine, see the four parts below.

``` r

library(timecoda)
```

## The example data

To model an outcome as a function of a time-use composition, which is
the independent variable, use
[`comp_iv()`](https://marcusvvlopes.github.io/timecoda/reference/comp_iv.md).
For its linear and generalized linear models, see
[`vignette("composition-predictor", package = "timecoda")`](https://marcusvvlopes.github.io/timecoda/articles/composition-predictor.md).

The examples use the bundled `timeuse` data. For a description of the
data, see
[`vignette("composition-variable", package = "timecoda")`](https://marcusvvlopes.github.io/timecoda/articles/composition-variable.md).
The data hold the body-mass-index z-score (`zbmi`), the obesity
indicator (`obesity`) and sex (`sex`) used in the other vignettes,
together with each person’s school (`school`, ten schools) and an
ordered academic-achievement band (`grade`, with the levels `below`,
`expected` and `above`). They also hold the follow-up time in months to
the onset of obesity or censoring (`time`), with its event indicator
(`event`, 1 for onset).

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

## Part 1. One builder, eleven engines

Every engine is fitted by
[`comp_iv()`](https://marcusvvlopes.github.io/timecoda/reference/comp_iv.md)
through the same procedure. The engines differ in the effect measure
reported by [`summary()`](https://rdrr.io/r/base/summary.html), in
whether `family` or `random` must be given, and in whether a prediction
on the scale of the outcome exists. You must give `family` for `"glm"`
and `"glmer"`, and the other engines reject it. `random` is required for
the clustered engines (`"coxme"`, `"lmer"`, `"glmer"`, `"clmm"` and
`"rlmer"`), optional for `"mblogit"` and rejected by the rest.

Reallocation estimates are contrasts. For details, see
[`vignette("predicting-reallocations", package = "timecoda")`](https://marcusvvlopes.github.io/timecoda/articles/predicting-reallocations.md).

| `model` | Backend | Use case | Effect measure |
|----|----|----|----|
| `"lm"` | [`stats::lm()`](https://rdrr.io/r/stats/lm.html) | A continuous outcome with independent observations | mean difference |
| `"glm"` | [`stats::glm()`](https://rdrr.io/r/stats/glm.html) | A binary, count or continuous outcome with a GLM family (e.g., logistic or Poisson regression) and independent observations | set by the link: odds ratio, rate, risk or mean ratio, or mean difference |
| `"cox"` | [`survival::coxph()`](https://rdrr.io/pkg/survival/man/coxph.html) | A time-to-event outcome with independent observations | hazard ratio |
| `"coxme"` | [`coxme::coxme()`](https://rdrr.io/pkg/coxme/man/coxme.html) | A time-to-event outcome with random effects for clustered or repeated data | hazard ratio |
| `"lmer"` | [`lme4::lmer()`](https://rdrr.io/pkg/lme4/man/lmer.html) | A continuous outcome with random effects for clustered or repeated data | mean difference |
| `"glmer"` | [`lme4::glmer()`](https://rdrr.io/pkg/lme4/man/glmer.html) | A binary, count or continuous outcome with a GLM family and random effects for clustered or repeated data | set by the link, as for `"glm"` |
| `"clm"` | [`ordinal::clm()`](https://rdrr.io/pkg/ordinal/man/clm.html) | An ordered-factor outcome (e.g., a graded scale) | cumulative odds ratio |
| `"clmm"` | [`ordinal::clmm()`](https://rdrr.io/pkg/ordinal/man/clmm.html) | An ordered-factor outcome with random effects for clustered or repeated data | cumulative odds ratio |
| `"mblogit"` | [`mclogit::mblogit()`](https://melff.github.io/mclogit/reference/mblogit.html) | An unordered-factor outcome with three or more categories, with optional random effects | relative-risk ratio, one per non-reference category |
| `"lmrob"` | [`robustbase::lmrob()`](https://rdrr.io/pkg/robustbase/man/lmrob.html) | A continuous outcome with outlying observations, down-weighted by MM-estimation | mean difference |
| `"rlmer"` | [`robustlmm::rlmer()`](https://rdrr.io/pkg/robustlmm/man/rlmer.html) | A continuous outcome with outlying observations and random effects for clustered or repeated data | mean difference |

Five engines (`"cox"`, `"coxme"`, `"clm"`, `"clmm"` and `"mblogit"`)
offer no single scale on which to predict the outcome. On these fits,
[`predict_realloc()`](https://marcusvvlopes.github.io/timecoda/reference/predict_realloc.md)
accepts only `type = "relative"`, its default, and
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) draws only the
relative surface. To draw one outcome level of a multinomial fit, give
that level to [`plot()`](https://rdrr.io/r/graphics/plot.default.html)
as `category`.

The backends in the coxme, ordinal, mclogit, robustbase and robustlmm
packages are suggested packages, and
[`comp_iv()`](https://marcusvvlopes.github.io/timecoda/reference/comp_iv.md)
checks for them when it is called. An example that uses one of them runs
only when that package is installed.

## Part 2. The engines in turn

For each family of engines, only the differences from the linear and
logistic summaries of
[`vignette("composition-predictor", package = "timecoda")`](https://marcusvvlopes.github.io/timecoda/articles/composition-predictor.md)
are described.

### Survival outcomes: `"cox"` and `"coxme"`

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

The effects are hazard ratios, and a reallocation on a Cox fit is a
hazard ratio as well. The report has no convergence line, because
[`comp_iv()`](https://marcusvvlopes.github.io/timecoda/reference/comp_iv.md)
does not probe Cox fits for convergence.

``` r

predict_realloc(object = fit_cox, from = "SB", to = "MVPA", delta = c(10, 30))
#> # A tibble: 2 × 10
#>   from  to    delta estimate std.error    df conf.low conf.high statistic
#>   <chr> <chr> <dbl>    <dbl>     <dbl> <dbl>    <dbl>     <dbl>     <dbl>
#> 1 SB    MVPA     10    0.856    0.0158   Inf    0.826     0.888     -8.44
#> 2 SB    MVPA     30    0.683    0.0304   Inf    0.626     0.745     -8.56
#> # ℹ 1 more variable: p.value <dbl>
```

We observed that moving 10 minutes from sedentary behaviour to
moderate-to-vigorous activity was associated with a hazard of obesity
onset 0.856 times that at the reference composition.

`"coxme"` extends the Cox model with random effects. To give each school
a random intercept (the term `(1 | school)`), set `random = "school"`.

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

One line is added to the summary, the median hazard ratio for school,
which is computed from the school variance `v` as
`exp(sqrt(2 * v) * qnorm(0.75))` (Austin et al., 2017). It is the median
hazard ratio between two people with identical covariates drawn from
different schools, taken as the person at higher risk over the other,
and it is therefore never below 1. A value of 1 would mean no variation
in hazard between schools.

### Multilevel models: `"lmer"` and `"glmer"`

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

We read the intraclass correlation for school, 0.264, as the share of
the total variance that lies between schools. With several grouping
terms (e.g., `random = c("region", "school")`), one such share is
reported per term.

``` r

predict_realloc(object = fit_lmer, from = "SB", to = "MVPA", delta = 10)
#> # A tibble: 1 × 10
#>   from  to    delta estimate std.error    df conf.low conf.high statistic
#>   <chr> <chr> <dbl>    <dbl>     <dbl> <dbl>    <dbl>     <dbl>     <dbl>
#> 1 SB    MVPA     10   -0.203    0.0119  288.   -0.226    -0.179     -17.0
#> # ℹ 1 more variable: p.value <dbl>
```

The `df` column is computed with the Kenward-Roger approximation
(Kenward & Roger, 1997), which `emmeans` applies through the pbkrtest
package when that package is installed. The same method is used for the
per-part rows of the summary. For the alternatives, see the help page
[`?comp_iv`](https://marcusvvlopes.github.io/timecoda/reference/comp_iv.md).

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

The between-school line is now a latent intraclass correlation. It is
the same share, computed on the latent scale of the logistic model with
`pi^2 / 3` in place of the residual variance (Goldstein et al., 2002).
With a value of 0.001, almost none of the variation in the odds of
obesity lies between schools.

### Ordinal outcomes: `"clm"` and `"clmm"`

To fit a cumulative-link model to an ordered factor, use `"clm"`. Here
the factor is the three grade bands (`below`, `expected` and `above`, in
that order).

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

With the parameterization of
[`ordinal::clm()`](https://rdrr.io/pkg/ordinal/man/clm.html)
(Christensen, 2026), we read an odds ratio above 1 as higher odds of a
higher grade band (i.e., of being at or above any given cut-point).

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

As for `"glmer"`, the latent intraclass correlation for school is
computed with `pi^2 / 3` in place of the residual variance.

### Unordered categories: `"mblogit"`

[`comp_iv()`](https://marcusvvlopes.github.io/timecoda/reference/comp_iv.md)
stops with an error when an ordered `grade` is given to this engine.

To fit `"mblogit"`, add an unordered copy of `grade` to the data, with
the same levels in the same order, and declare the composition again
from it.

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

`below`, the first level, is the reference category. Each effect is a
relative-risk ratio comparing one of the other two categories with the
reference, and the new `response_level` column names the category. The
omnibus test has 6 degrees of freedom, which are three coordinates for
each non-reference category. On this fit, a reallocation also gives one
row per non-reference category.

``` r

predict_realloc(object = fit_mb, from = "SB", to = "MVPA", delta = 10)
#> # A tibble: 2 × 11
#>   from  to    delta response_level estimate std.error    df conf.low conf.high
#>   <chr> <chr> <dbl> <chr>             <dbl>     <dbl> <dbl>    <dbl>     <dbl>
#> 1 SB    MVPA     10 expected           1.13    0.0441   Inf     1.05      1.22
#> 2 SB    MVPA     10 above              1.18    0.0484   Inf     1.09      1.28
#> # ℹ 2 more variables: statistic <dbl>, p.value <dbl>
```

To fit a clustered model with `"mblogit"`, add `random`.
`random = "school"` adds a random intercept for each school, and a bar
term such as `random = "(ilr1 | school)"` adds a random slope on a
coordinate. For details, see the help page
[`?comp_iv`](https://marcusvvlopes.github.io/timecoda/reference/comp_iv.md).

### Robust linear models: `"lmrob"` and `"rlmer"`

To fit a linear model by MM-estimation (Yohai, 1987), use `"lmrob"`.

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

To fit the robust linear mixed model (Koller, 2016), use `"rlmer"`.

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

Because the robustlmm package reports no degrees of freedom, every
standard error, interval and p-value of an `"rlmer"` fit is asymptotic,
in the summary and in the reallocations alike. The `df` column of a
reallocation is `Inf`, compared with 288 for the same reallocation on
the `"lmer"` fit.

``` r

predict_realloc(object = fit_rlmer, from = "SB", to = "MVPA", delta = 10)
#> # A tibble: 1 × 10
#>   from  to    delta estimate std.error    df conf.low conf.high statistic
#>   <chr> <chr> <dbl>    <dbl>     <dbl> <dbl>    <dbl>     <dbl>     <dbl>
#> 1 SB    MVPA     10   -0.206    0.0124   Inf   -0.230    -0.182     -16.6
#> # ℹ 1 more variable: p.value <dbl>
```

## Part 3. The formula interface

### The `ilr()` marker

To specify the model as a formula, give `formula` in place of `outcome`,
`covariates`, `random`, `time`, `event` and `moderator`. Within it, the
marker
[`ilr()`](https://marcusvvlopes.github.io/timecoda/reference/ilr.md)
stands for the whole composition. The discrete call
`outcome = "zbmi", covariates = "sex"` can be written as the formula
`zbmi ~ ilr() + sex`, and the printed fit shows the call with
[`ilr()`](https://marcusvvlopes.github.io/timecoda/reference/ilr.md)
expanded.

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

Every other term is passed to the backend unchanged, so a formula can
hold interactions, random terms and transformed variables.
[`ilr()`](https://marcusvvlopes.github.io/timecoda/reference/ilr.md) is
expanded by two rules
([`?ilr`](https://marcusvvlopes.github.io/timecoda/reference/ilr.md)).

The first rule applies inside the formula operators `+`, `-`, `*`, `:`,
`/`, `^`, `%in%`, `|` and `||` and inside parentheses. There,
[`ilr()`](https://marcusvvlopes.github.io/timecoda/reference/ilr.md) is
replaced by the sum `(ilr1 + ... + ilr(D-1))`, and R’s formula algebra
applies to it. In an lme4-style engine, `(ilr() | school)` is a random
slope on every coordinate. `ilr()^2` expands to the main effects and
every pairwise product of the coordinates, with no squares.

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

The second rule applies inside any other call, where
[`ilr()`](https://marcusvvlopes.github.io/timecoda/reference/ilr.md) is
distributed. The call is repeated once per coordinate and the copies are
added, and `I(ilr()^2)` becomes `(I(ilr1^2) + I(ilr2^2) + I(ilr3^2))`.
Combining both rules, the formula `ilr()^2 + I(ilr()^2)` is the full
quadratic in the composition, with main effects, pairwise products and
squares.

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

The second rule distributes `poly(ilr(), 2)` and `splines::ns(ilr(), 3)`
in the same way.
[`comp_iv()`](https://marcusvvlopes.github.io/timecoda/reference/comp_iv.md)
warns that the estimates may be biased when the formula holds
[`ilr()`](https://marcusvvlopes.github.io/timecoda/reference/ilr.md) but
no coordinate appears as a term of its own or as the first argument of
[`poly()`](https://rdrr.io/r/stats/poly.html), `ns()` or `bs()`.

### Moderation

To test the interaction of one covariate with the composition, name that
covariate with `moderator`. The fit is the model `zbmi ~ ilr() * sex` of
`fit_int` with the moderation test added.

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

The moderation line reports a joint Wald test of the null hypothesis
that every composition-by-sex coefficient is zero (i.e., that the
association between the composition and z-BMI does not differ between
the sexes). For each part, the three new columns (`mod_statistic`,
`mod_df` and `mod_p.value`) hold the test of the null hypothesis that
the part’s own interaction with sex is zero. For sleep, the test gives p
= 0.0322, with no adjustment for the four tests.

For reallocation estimates on a fit with a composition-by-covariate
interaction, see
[`vignette("predicting-reallocations", package = "timecoda")`](https://marcusvvlopes.github.io/timecoda/articles/predicting-reallocations.md).

For a factor `moderator`, `by` defaults to the moderator when neither
`by` nor `at` names it. `fit_int`, the formula fit, is the same model,
although a formula carries no default `by`. A reallocation on `fit_int`
depends on sex, and
[`predict_realloc()`](https://marcusvvlopes.github.io/timecoda/reference/predict_realloc.md)
therefore stops with an error until you give a value for every covariate
that interacts with the composition, through `at` or `by`.

### Non-linear composition terms

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

The omnibus test now has 9 degrees of freedom, one for each coefficient
that involves only the coordinates. The curvature line reports a joint
Wald test of the null hypothesis that the coefficients of the squares
and products are all zero (i.e., that the association between the
composition and z-BMI has no curvature). Because a rotated fit with
squares and products has no single `ilr1` coefficient, each per-part row
is the slope of the linear predictor in `ilr1` at the reference
composition, taken from
[`emmeans::emtrends()`](https://rvlenth.github.io/emmeans/reference/emtrends.html).
[`predict_realloc()`](https://marcusvvlopes.github.io/timecoda/reference/predict_realloc.md)
reads the quadratic fit in the same way as a linear one.

``` r

predict_realloc(object = fit_quad, from = "SB", to = "MVPA", delta = c(10, 30))
#> # A tibble: 2 × 10
#>   from  to    delta estimate std.error    df conf.low conf.high statistic
#>   <chr> <chr> <dbl>    <dbl>     <dbl> <dbl>    <dbl>     <dbl>     <dbl>
#> 1 SB    MVPA     10   -0.215    0.0160   289   -0.247    -0.184     -13.5
#> 2 SB    MVPA     30   -0.537    0.0429   289   -0.621    -0.452     -12.5
#> # ℹ 1 more variable: p.value <dbl>
```

We observed that the mean z-BMI difference associated with moving 30
minutes is less than three times that for 10 minutes.

[`plot()`](https://rdrr.io/r/graphics/plot.default.html) draws the
fitted surface as one ternary panel for each group of three parts. For
this figure on a logistic fit, see
[`vignette("composition-predictor", package = "timecoda")`](https://marcusvvlopes.github.io/timecoda/articles/composition-predictor.md).

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

The full quadratic is the same model in every ILR basis, so its fit and
reallocations do not depend on the order of the parts. A squares-only
specification, `ilr() + I(ilr()^2)`, and `poly(ilr(), 2)` lack this
invariance, and their fits depend on the order of the parts given to
[`declare_composition()`](https://marcusvvlopes.github.io/timecoda/reference/declare_composition.md).
For an example on `timeuse`, see the help page
[`?ilr`](https://marcusvvlopes.github.io/timecoda/reference/ilr.md).

### The total as a covariate

The total of a composition of the waking behaviours alone
(moderate-to-vigorous activity, light activity and sedentary behaviour)
varies between people. To add this total to the model, set
`include_total = TRUE`, which enters it as the covariate `t_total` of
Coenders et al. (2017), `t = sqrt(D) * log(g)`, where `g` is the
geometric mean of a person’s parts on their original scale. The ILR
coefficients then give the association between the outcome and the
shares of waking time, at a constant waking time. The coefficient of
`t_total` gives the difference in the outcome associated with a longer
waking day, in which every waking behaviour is larger in the same
proportion and sleep is shorter.

The waking composition is declared from the same data with three parts.
Declared with `closure = "mean"`, it is reported on the scale of the
average waking time, and the row totals before closure, shown in its
summary, are the waking times.

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

The null hypothesis of the one-degree-of-freedom Wald test of `t_total`
is that its coefficient is zero (i.e., that waking time and z-BMI are
not associated at constant shares of waking time). At constant shares, a
longer waking day, and so a shorter time asleep, is associated with a
higher z-BMI in these data, as the positive coefficient shows.

``` r

predict_realloc(object = fit_tot, from = "SB", to = "MVPA", delta = 10)
#> # A tibble: 1 × 10
#>   from  to    delta estimate std.error    df conf.low conf.high statistic
#>   <chr> <chr> <dbl>    <dbl>     <dbl> <dbl>    <dbl>     <dbl>     <dbl>
#> 1 SB    MVPA     10   -0.392    0.0410   296   -0.473    -0.312     -9.56
#> # ℹ 1 more variable: p.value <dbl>
```

In a reallocation, every composition is evaluated with `t_total` held at
its reference value, the sample mean of `t_total`, and the estimate is
therefore the difference associated with the move at a constant waking
time. On the fixed-total `cd`, `include_total = TRUE` ends with an
error.

## Part 4. Diagnostics per engine

The diagnostic backend used by
[`fit_check()`](https://marcusvvlopes.github.io/timecoda/reference/fit_check.md)
depends on the model class.

### Models checked with DHARMa

For `"lm"`, `"glm"`, `"lmer"` and `"glmer"` fits,
[`fit_check()`](https://marcusvvlopes.github.io/timecoda/reference/fit_check.md)
uses the DHARMa package (Hartig, 2026). For the simulated quantile
residuals, see
[`vignette("composition-predictor", package = "timecoda")`](https://marcusvvlopes.github.io/timecoda/articles/composition-predictor.md).
The rows reported for a `"glmer"` fit are those of a `"glm"` fit.

### Cox models

A `"cox"` fit is checked with
[`survival::cox.zph()`](https://rdrr.io/pkg/survival/man/cox.zph.html),
a test of the proportional-hazards assumption based on the scaled
Schoenfeld residuals (Grambsch & Therneau, 1994).

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

In the global test, we found no evidence of a departure from
proportional hazards (p = 0.906). To draw the Schoenfeld residuals of
each term against time, set `plot = TRUE`.

### Ordinal models

A `"clm"` fit is checked with
[`ordinal::nominal_test()`](https://rdrr.io/pkg/ordinal/man/nominal.test.html)
(Christensen, 2026) and with the sure package (Greenwell et al., 2018),
which provides surrogate residuals (Liu & Zhang, 2018).

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

In the proportional-odds tests, we found no evidence of a departure from
proportional odds for any of the three coordinates (p-values of 0.262
and above). The surrogate residuals are random draws made under a fixed
seed and therefore repeat from call to call. To draw them, set
`plot = TRUE`.

### Robust models

`"lmrob"` and `"rlmer"` fits are checked through the robustness weights
of their estimates.

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

In this fit, the downweighted share of 0.01 corresponds to three of the
300 people, each weighted below 0.5. The rows reported for an `"rlmer"`
fit are the same.

### Engines with fit statistics only

For `"clmm"`, `"mblogit"` and `"coxme"` fits, no residual diagnostic is
run, and
[`fit_check()`](https://marcusvvlopes.github.io/timecoda/reference/fit_check.md)
reports the AIC and the log-likelihood with a note and draws nothing. A
`"clmm"` fit therefore has no proportional-odds test, and a `"coxme"`
fit has no proportional-hazards test. The log-likelihood of a `"coxme"`
fit is the penalized one, and the AIC uses the fit’s effective degrees
of freedom.

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
  the linear and logistic fits, the summary read line by line, and the
  fitted surface.
- [`vignette("predicting-reallocations", package = "timecoda")`](https://marcusvvlopes.github.io/timecoda/articles/predicting-reallocations.md):
  reallocation estimates, their schemes and anchors, the publication
  table and the reallocation figure.
- [`vignette("optimal-composition", package = "timecoda")`](https://marcusvvlopes.github.io/timecoda/articles/optimal-composition.md):
  the search for the Goldilocks zone, the compositions with the best
  predicted outcome.
- [`vignette("timecoda", package = "timecoda")`](https://marcusvvlopes.github.io/timecoda/articles/timecoda.md):
  the whole workflow in one analysis.
