# Getting started with timecoda

`timecoda` supports time-use compositional data analysis (CoDA) from
data declaration through regression modelling, isotemporal time
reallocations, and optimal-composition search. We assume basic
familiarity with R regression workflows and core CoDA concepts (simplex
closure, isometric log-ratio transformations, and pivot coordinates).

``` r

library(timecoda)
```

Every analysis follows a single pipeline built around two objects:
`compD` (a validated, closed composition) and `comp_iv` (a regression
model with the composition as predictor).

![Flowchart of the timecoda workflow. A data frame with the parts and an
outcome enters declare_composition(), which returns a compD object. The
compD is described with summary(), plot() and pub_ready(), and modelled
as a predictor of an outcome with comp_iv(). A comp_iv() fit leads to
fit_check() for diagnostics, to predict_realloc() for reallocation
estimates, which are summarized with summary() and pub_ready() and drawn
with plot_realloc(), and to optimal_composition() for a search of the
fit for the best compositions. A second path, drawn in grey and dashed,
would model the composition as the outcome with comp_dv(), described
with summary() and plot(). A caption at the bottom explains that grey
marks what is planned and not yet implemented.](workflow.svg)

*(Modelling the composition as the outcome via `comp_dv()` is planned
for a future release.)*

## Declaring and describing a composition

We will use the bundled `timeuse` dataset, a simulated sample of 300
participants whose 1,440-minute days are partitioned into
moderate-to-vigorous physical activity (`MVPA`), light physical activity
(`LPA`), sedentary behaviour (`SB`), and sleep period time (`SPT`),
alongside BMI z-score (`zbmi`) and binary obesity status (`obesity`).

``` r

parts <- c("MVPA", "LPA", "SB", "SPT")
head(timeuse[c(parts, "zbmi", "obesity")])
#> # A tibble: 6 × 6
#>    MVPA   LPA    SB   SPT   zbmi obesity
#>   <dbl> <dbl> <dbl> <dbl>  <dbl> <fct>  
#> 1  36.1  320   541.  543. -1.31  yes    
#> 2  14.7  214.  815.  397  -0.472 no     
#> 3 195.   206.  500.  539. -1.45  no     
#> 4  95.5  229.  445.  670. -2.40  no     
#> 5   1.7  504.  479.  455.  2.56  no     
#> 6  11.2  150.  583.  696. -0.259 no
```

To prepare these columns for analysis, we declare the composition once
with
[`declare_composition()`](https://marcusvvlopes.github.io/timecoda/reference/declare_composition.md).
It validates the parts, imputes any zeroes, closes rows to 1,440
minutes, and stores the metadata so downstream functions reuse the same
labels and scale:

``` r

cd <- declare_composition(
  data    = timeuse,
  parts   = c("MVPA", "LPA", "SB", "SPT"),
  closure = 1440,
  labels  = c(
    MVPA = "Moderate-to-vigorous PA",
    LPA  = "Light PA",
    SB   = "Sedentary behaviour",
    SPT  = "Sleep"
  )
)

summary(object = cd)
#> <compD> compositional variable summary
#>   Parts: Moderate-to-vigorous PA, Light PA, Sedentary behaviour, Sleep
#>   Closure: 1440
#>   Observations: 300
#> 
#>   Row totals (before closure): mean 1440.0 (SD 0.0), range 1440.0-1440.0
#>   No zeroes detected.
#> 
#> Descriptive statistics
#>   geom_mean: imputed, closed composition; other columns: original data
#> # A tibble: 4 × 10
#>   part  label               geom_mean  mean    sd median   p25   p75   min   max
#>   <chr> <chr>                   <dbl> <dbl> <dbl>  <dbl> <dbl> <dbl> <dbl> <dbl>
#> 1 MVPA  Moderate-to-vigoro…      28.1  45.9  44.7   35.0  12.7  61.8   0.1  240 
#> 2 LPA   Light PA                235.  243.   86.9  232.  173.  308.   51.8  504.
#> 3 SB    Sedentary behaviour     677.  661.  118.   669.  586.  743.  316.  1037.
#> 4 SPT   Sleep                   499.  491.  106.   489.  426.  559.  240.   815.
#> 
#> Variation matrix
#> # A tibble: 4 × 6
#>   part  label                    MVPA   LPA    SB   SPT
#>   <chr> <chr>                   <dbl> <dbl> <dbl> <dbl>
#> 1 MVPA  Moderate-to-vigorous PA  0    1.67  1.55  1.44 
#> 2 LPA   Light PA                 1.67 0     0.247 0.260
#> 3 SB    Sedentary behaviour      1.55 0.247 0     0.140
#> 4 SPT   Sleep                    1.44 0.260 0.140 0
```

In the descriptive table, `geom_mean` gives the closed compositional
centre summing to 1,440 minutes, while the remaining columns summarize
the raw data. To format the summary for a manuscript with `gt`, pass it
to
[`pub_ready()`](https://marcusvvlopes.github.io/timecoda/reference/pub_ready.md):

``` r

pub_ready(x = summary(object = cd), digits = 1)$desc
```

| Descriptive summary |  |  |  |
|----|----|----|----|
| Spread of the original data and the compositional centre per part |  |  |  |
| Part | Mean ± SD (original data) | Geometric mean (closed to 1440) | % of total |
| Moderate-to-vigorous PA | 45.9 ± 44.7 | 28.1 | 2.0 |
| Light PA | 242.6 ± 86.9 | 235.5 | 16.4 |
| Sedentary behaviour | 660.7 ± 117.6 | 677.2 | 47.0 |
| Sleep | 490.8 ± 106.2 | 499.2 | 34.7 |

We can inspect the marginal distributions, pairwise ternary projections
(with the compositional centre and 95% confidence region), and variation
matrix together using
[`plot()`](https://rdrr.io/r/graphics/plot.default.html):

``` r

plot(x = cd)
```

![Density curves and boxplots of the four parts in minutes per day; four
ternary plots, one for each set of three parts, each with the geometric
mean marked by dashed crosshairs and a 95% confidence region drawn as a
closed curve; and a heatmap of the variation matrix between the four
parts.](timecoda_files/figure-html/plot-cd-1.png)![Density curves and
boxplots of the four parts in minutes per day; four ternary plots, one
for each set of three parts, each with the geometric mean marked by
dashed crosshairs and a 95% confidence region drawn as a closed curve;
and a heatmap of the variation matrix between the four
parts.](timecoda_files/figure-html/plot-cd-2.png)![Density curves and
boxplots of the four parts in minutes per day; four ternary plots, one
for each set of three parts, each with the geometric mean marked by
dashed crosshairs and a 95% confidence region drawn as a closed curve;
and a heatmap of the variation matrix between the four
parts.](timecoda_files/figure-html/plot-cd-3.png)![Density curves and
boxplots of the four parts in minutes per day; four ternary plots, one
for each set of three parts, each with the geometric mean marked by
dashed crosshairs and a 95% confidence region drawn as a closed curve;
and a heatmap of the variation matrix between the four
parts.](timecoda_files/figure-html/plot-cd-4.png)![Density curves and
boxplots of the four parts in minutes per day; four ternary plots, one
for each set of three parts, each with the geometric mean marked by
dashed crosshairs and a 95% confidence region drawn as a closed curve;
and a heatmap of the variation matrix between the four
parts.](timecoda_files/figure-html/plot-cd-5.png)![Density curves and
boxplots of the four parts in minutes per day; four ternary plots, one
for each set of three parts, each with the geometric mean marked by
dashed crosshairs and a 95% confidence region drawn as a closed curve;
and a heatmap of the variation matrix between the four
parts.](timecoda_files/figure-html/plot-cd-6.png)![Density curves and
boxplots of the four parts in minutes per day; four ternary plots, one
for each set of three parts, each with the geometric mean marked by
dashed crosshairs and a 95% confidence region drawn as a closed curve;
and a heatmap of the variation matrix between the four
parts.](timecoda_files/figure-html/plot-cd-7.png)

For zero imputation, multivariate outlier screening, and varying-total
closures, see
[`vignette("composition-variable", package = "timecoda")`](https://marcusvvlopes.github.io/timecoda/articles/composition-variable.md).

## Modelling an outcome from the composition

With `cd` declared, we can test whether daily time use is associated
with obesity. Here we fit a logistic regression with
[`comp_iv()`](https://marcusvvlopes.github.io/timecoda/reference/comp_iv.md),
which handles the pivot ILR transformation behind the scenes (linear,
multilevel, survival, ordinal, and robust engines use the same
interface; see
[`vignette("advanced-models", package = "timecoda")`](https://marcusvvlopes.github.io/timecoda/articles/advanced-models.md)):

``` r

fit_ob <- comp_iv(comp = cd, outcome = "obesity", model = "glm",
                  family = "binomial")
fit_ob
#> <timecoda_fit> compositional generalized linear model (binomial)
#>   Parts: MVPA, LPA, SB, SPT (total = 1440)
#>   Reference composition: MVPA = 28.1, LPA = 235.5, SB = 677.2, SPT = 499.2
#> 
#> 
#> Call:  stats::glm(formula = obesity ~ ilr1 + ilr2 + ilr3, family = "binomial", 
#>     data = data)
#> 
#> Coefficients:
#> (Intercept)         ilr1         ilr2         ilr3  
#>     -3.0796      -0.6916      -0.2635       0.9014  
#> 
#> Degrees of Freedom: 299 Total (i.e. Null);  296 Residual
#> Null Deviance:       355.8 
#> Residual Deviance: 322   AIC: 330
```

Printing the fit shows the underlying model call in the initial pivot
basis. The summary runs an order-invariant omnibus Wald test across all
coordinates and rotates the pivot basis D times so each part appears
first against the remaining D - 1 parts:

``` r

summary(object = fit_ob)
#> <comp_iv> compositional generalized linear model (binomial)
#>   Parts: Moderate-to-vigorous PA, Light PA, Sedentary behaviour, Sleep (total = 1440)
#>   Omnibus composition test (Wald): chi-squared(3) = 28.922, p = 2.33e-06
#>   Convergence: OK
#> 
#> Compositional effects (Odds ratio)
#> # A tibble: 5 × 6
#>   parameter          label                   estimate conf.low conf.high p.value
#>   <chr>              <chr>                      <dbl>    <dbl>     <dbl>   <dbl>
#> 1 Global Composition Global Composition        NA       NA        NA     2.33e-6
#> 2 MVPA vs remaining  Moderate-to-vigorous PA    0.501    0.382     0.656 4.89e-7
#> 3 LPA vs remaining   Light PA                   0.982    0.506     1.91  9.58e-1
#> 4 SB vs remaining    Sedentary behaviour        2.98     1.16      7.65  2.34e-2
#> 5 SPT vs remaining   Sleep                      0.683    0.273     1.71  4.15e-1
```

An odds ratio below 1, as for `MVPA`, means that a day with more `MVPA`
relative to the rest of the day is associated with lower odds of
obesity. Because a one-unit change in a log-ratio coordinate has no
direct reading in minutes, we use reallocation estimates next to express
the same model on the minute scale. See
[`vignette("composition-predictor", package = "timecoda")`](https://marcusvvlopes.github.io/timecoda/articles/composition-predictor.md)
for full model interpretation and response surfaces.

## Checking the fit

Before interpreting predictions, we check residual diagnostics with
[`fit_check()`](https://marcusvvlopes.github.io/timecoda/reference/fit_check.md),
which uses simulated quantile residuals from `DHARMa` for linear and
generalized linear models:

``` r

fit_check(object = fit_ob)
#> Warning in newton(lsp = lsp, X = G$X, y = G$y, Eb = G$Eb, UrS = G$UrS, L = G$L,
#> : Fitting terminated with step failure - check results carefully
#> <fit_check> generalized linear model 
#>   Backend: dharma 
#>   Convergence: converged (0 backend warnings)
#> 
#> # A tibble: 3 × 4
#>   check         statistic    df p.value
#>   <chr>             <dbl> <dbl>   <dbl>
#> 1 KS uniformity    0.0381    NA   0.777
#> 2 Dispersion       1.01      NA   0.976
#> 3 Outliers         0         NA   1
```

None of the diagnostic checks flags a lack of fit.

## Estimating time reallocations

To express the logistic model in minutes, we use compositional
isotemporal substitution via
[`predict_realloc()`](https://marcusvvlopes.github.io/timecoda/reference/predict_realloc.md).
By default, reallocations start from the sample’s closed geometric mean
and return relative contrasts on the model’s effect scale (here, an odds
ratio). Moving 15 minutes from sedentary behaviour to
moderate-to-vigorous activity, holding sleep and light activity
constant, is associated with lower odds of obesity:

``` r

predict_realloc(object = fit_ob, from = "SB", to = "MVPA", delta = 15)
#> # A tibble: 1 × 10
#>   from  to    delta estimate std.error    df conf.low conf.high statistic
#>   <chr> <chr> <dbl>    <dbl>     <dbl> <dbl>    <dbl>     <dbl>     <dbl>
#> 1 SB    MVPA     15    0.758    0.0400   Inf    0.683     0.840     -5.26
#> # ℹ 1 more variable: p.value <dbl>
```

To see the predicted probability of obesity at both compositions instead
of the odds ratio, switch to `type = "absolute"`:

``` r

predict_realloc(object = fit_ob, from = "SB", to = "MVPA", delta = 15,
                type = "absolute")
#> # A tibble: 2 × 8
#>   from  to    delta estimate std.error    df conf.low conf.high
#>   <chr> <chr> <dbl>    <dbl>     <dbl> <dbl>    <dbl>     <dbl>
#> 1 NA    NA        0    0.258    0.0271   Inf    0.209     0.315
#> 2 SB    MVPA     15    0.209    0.0269   Inf    0.161     0.266
```

You can also evaluate several reallocation magnitudes in one call by
passing a vector to `delta`, then format the sweep with
[`pub_ready()`](https://marcusvvlopes.github.io/timecoda/reference/pub_ready.md):

``` r

r <- predict_realloc(object = fit_ob, from = "SB", to = "MVPA",
                     delta = c(10, 20, 30))
summary(object = r)
#> Time-reallocation effects: Odds ratio versus the reference composition
#>   3 reallocations of 10, 20, 30; 95% confidence intervals
#> 
#> # A tibble: 3 × 9
#>   to    from  to_label      from_label delta estimate conf.low conf.high p.value
#>   <chr> <chr> <chr>         <chr>      <dbl>    <dbl>    <dbl>     <dbl>   <dbl>
#> 1 MVPA  SB    Moderate-to-… Sedentary…    10    0.822    0.764     0.884 1.51e-7
#> 2 MVPA  SB    Moderate-to-… Sedentary…    20    0.704    0.618     0.802 1.36e-7
#> 3 MVPA  SB    Moderate-to-… Sedentary…    30    0.620    0.519     0.740 1.25e-7
```

``` r

pub_ready(x = summary(object = r), digits = 2)$effects
```

[TABLE]

To visualize reallocation curves across every pair of behaviours up to
60 minutes, pass the fit to
[`plot_realloc()`](https://marcusvvlopes.github.io/timecoda/reference/plot_realloc.md):

``` r

plot_realloc(object = fit_ob, delta_max = 60)
```

![Reallocation plot in a lower-triangular grid of six panels, one per
pair of parts, each showing the odds ratio of obesity with its
confidence ribbon against the minutes reallocated between the two parts,
and a dashed reference line at an odds ratio of
1.](timecoda_files/figure-html/realloc-plot-1.png)

See
[`vignette("predicting-reallocations", package = "timecoda")`](https://marcusvvlopes.github.io/timecoda/articles/predicting-reallocations.md)
for asymmetric reversals, one-to-all schemes, custom reference
compositions, and subgroup moderation.

## Finding the optimal composition

To identify daily compositions with the best predicted outcome (the
Goldilocks zone), we search a feasible grid inside the empirical sample
range with
[`optimal_composition()`](https://marcusvvlopes.github.io/timecoda/reference/optimal_composition.md).
Because `fit_ob` is linear in ILR coordinates, predictions keep
improving in one direction and the zone lands on the edge of the
feasible region; refitting with `curvature = TRUE` allows an interior
optimum (see
[`vignette("optimal-composition", package = "timecoda")`](https://marcusvvlopes.github.io/timecoda/articles/optimal-composition.md)).

``` r

opt <- optimal_composition(object = fit_ob, direction = "min")
#> "obesity" is linear in the composition, so its best predicted outcome lies on
#> the edge of the region.
#> ℹ The zone then marks where the region ends; see on_boundary.
#> ℹ A curved surface can hold an interior optimum: refit with `curvature = TRUE`
#>   (or `ilr()^2 + I(ilr()^2)` in a formula).
opt
#> <optimal_comp> Goldilocks-zone search: obesity
#> 
#> Search
#>   Parts       Moderate-to-vigorous PA, Light PA, Sedentary behaviour, Sleep
#>               (total = 1440)
#>   Region      per-part box of the mean +/- 3 SD
#>   Grid        multiples of 10; 36,611 candidates
#>   Scale       Odds ratio versus the reference day
#>   Zone rule   best 5% of predictions, lower is better
#> 
#> Region limits
#>                              lower     upper
#>   Moderate-to-vigorous PA      0.1*    180.0
#>   Light PA                    51.8*    503.3
#>   Sedentary behaviour        315.6*   1013.4
#>   Sleep                      239.9*    809.4
#>   * the part's observed minimum
#> 
#> Zone for obesity: 1,831 compositions (centre [min; max])
#>   Moderate-to-vigorous PA   151.2  [ 90.0; 170.0]
#>   Light PA                  262.2  [ 60.0; 500.0]
#>   Sedentary behaviour       367.9  [320.0; 450.0]
#>   Sleep                     658.7  [400.0; 800.0]
#>   A fit linear in the composition has no interior optimum: the zone is the edge
#>   of the region in the direction of benefit.
#>   Refit with `curvature = TRUE` for a surface that can hold one.
#> 
#> Gain from the reference day to the zone centre (95% CI)
#>   obesity   Odds ratio   0.187  [0.0842; 0.414]  p = 3.65e-05
```

``` r

plot(x = opt)
```

![Four ternary plots of the search, one for each set of three parts with
the fourth part held at the zone centre. Inside the searched region,
tiles are coloured by the odds ratio of obesity against the reference
day, a red outline surrounds the zone and a red point marks its
centre.](timecoda_files/figure-html/optimal-plot-1.png)![Four ternary
plots of the search, one for each set of three parts with the fourth
part held at the zone centre. Inside the searched region, tiles are
coloured by the odds ratio of obesity against the reference day, a red
outline surrounds the zone and a red point marks its
centre.](timecoda_files/figure-html/optimal-plot-2.png)![Four ternary
plots of the search, one for each set of three parts with the fourth
part held at the zone centre. Inside the searched region, tiles are
coloured by the odds ratio of obesity against the reference day, a red
outline surrounds the zone and a red point marks its
centre.](timecoda_files/figure-html/optimal-plot-3.png)![Four ternary
plots of the search, one for each set of three parts with the fourth
part held at the zone centre. Inside the searched region, tiles are
coloured by the odds ratio of obesity against the reference day, a red
outline surrounds the zone and a red point marks its
centre.](timecoda_files/figure-html/optimal-plot-4.png)

In each ternary panel, the fourth behaviour is held at the zone centre
(red point), and the red outline marks the zone boundary.

## See also

- [`vignette("composition-variable", package = "timecoda")`](https://marcusvvlopes.github.io/timecoda/articles/composition-variable.md):
  validating data, imputing zeroes, flagging multivariate outliers, and
  descriptive plots.
- [`vignette("composition-predictor", package = "timecoda")`](https://marcusvvlopes.github.io/timecoda/articles/composition-predictor.md):
  modelling an outcome from a composition and checking the fit.
- [`vignette("predicting-reallocations", package = "timecoda")`](https://marcusvvlopes.github.io/timecoda/articles/predicting-reallocations.md):
  reallocation estimates, publication tables, and transfer curves.
- [`vignette("optimal-composition", package = "timecoda")`](https://marcusvvlopes.github.io/timecoda/articles/optimal-composition.md):
  single- and multi-outcome Goldilocks-zone optimization.
- [`vignette("advanced-models", package = "timecoda")`](https://marcusvvlopes.github.io/timecoda/articles/advanced-models.md):
  multilevel, survival, ordinal, robust, and formula-based models.
- Check each function’s help page for full argument lists and runnable
  examples.
