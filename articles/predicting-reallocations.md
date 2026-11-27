# Predicting time reallocations

``` r

library(timecoda)
```

## Reallocating time

[`predict_realloc()`](https://marcusvvlopes.github.io/timecoda/reference/predict_realloc.md)
takes a model fitted by
[`comp_iv()`](https://marcusvvlopes.github.io/timecoda/reference/comp_iv.md)
and estimates the difference in an outcome associated with moving a
stated number of minutes from one part of a composition to another.

The examples use the bundled `timeuse` data, with the composition
declared once as `cd`. For a description of these data, see
[`vignette("composition-variable", package = "timecoda")`](https://marcusvvlopes.github.io/timecoda/articles/composition-variable.md).
Two models are fitted to it. For fitting these models, see
[`vignette("composition-predictor", package = "timecoda")`](https://marcusvvlopes.github.io/timecoda/articles/composition-predictor.md).
`fit_bmi` is a linear model of the body-mass-index z-score (`zbmi`), and
`fit_ob` is a logistic regression of an obesity indicator (`obesity`).

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

fit_bmi <- comp_iv(comp = cd, outcome = "zbmi", model = "lm")
fit_ob  <- comp_iv(comp = cd, outcome = "obesity", model = "glm",
                   family = "binomial")
```

## Part 1. What a reallocation estimate is

### Two days compared

A reallocation estimate compares two whole days. One is the reference
composition, a typical day for the sample, and the other is the same day
with `delta` minutes moved from one part, the donor, to another, the
receiver, while the remaining parts keep their reference values. The
estimate is the difference between the outcomes the model predicts for
the two days (Dumuid, Stanford, Pedišić, et al., 2018; Dumuid et al.,
2020). This approach is the compositional isotemporal substitution model
of Dumuid et al. (2019).

The difference is taken on the link scale of the model. For a linear
model, which has the identity link, it is a mean difference in the
outcome itself, here in z-BMI units. For a logistic regression, which
has the logit link, it is a difference in log-odds, and exponentiating
it gives an odds ratio, that is, the odds of the outcome on the
reallocated day divided by the odds on the reference day.

The estimates are associations in observational data. Each compares the
outcome predicted for two different compositions across the sample, and
none measures the change that would follow if a person altered their own
day (Dumuid, Stanford, Pedišić, et al., 2018).

### The reference composition

The reference composition is the closed geometric mean of the sample,
and every estimate in this vignette is anchored at it (Dumuid, Stanford,
Pedišić, et al., 2018; Dumuid et al., 2020) except those in the section
on a chosen reference day. To see the reference composition, print the
fit.

``` r

print(x = fit_ob)
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

Because both fits share the declared composition, `fit_bmi` has the same
reference composition. The reference day holds only 28.1 minutes of
moderate-to-vigorous activity, so a transfer of 28.1 minutes or more out
of that part leaves it at zero or below (see the section on invalid
transfers in Part 2).

### What `predict_realloc()` does

A call to
[`predict_realloc()`](https://marcusvvlopes.github.io/timecoda/reference/predict_realloc.md)
runs six steps.

![Flowchart of the predict_realloc() pipeline: take the reference
composition stored in the fit, build each reallocated composition by
moving delta minutes from one part to another, check that every part
stays positive with an invalid transfer handled by on_invalid, express
every composition in ILR coordinates on an emmeans reference grid with
covariates averaged, contrast each reallocated composition against the
reference on the link scale or, when type is absolute, predict the
response at each composition, back-transform to the model's effect
measure, and return a realloc tibble read by summary(), pub_ready() and
plot_realloc().](realloc-pipeline.svg)

### A complete call

``` r

r <- predict_realloc(
  object = fit_ob,
  from   = "SB",
  to     = "MVPA",
  delta  = c(10, 20, 30)
)
```

### Reading the summary

``` r

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

Each `p.value` comes from a test of the null hypothesis that the odds
ratio in its row equals 1, that is, that the odds of obesity do not
differ between the reallocated day and the reference day. The p-values
are not adjusted for the number of rows.

We read an odds ratio of 0.822 for 10 minutes as odds of obesity about
18% lower on the reallocated day than on the reference day. Each further
10 minutes multiplies the odds ratio by a factor closer to 1, because
the model works on log-ratios.

## Part 2. The arguments in turn

### Direction and size

To reverse the transfer, swap `from` and `to`, here moving time from
moderate-to-vigorous activity to sedentary behaviour.

``` r

predict_realloc(object = fit_ob, from = "MVPA", to = "SB", delta = c(10, 20, 30))
#> # A tibble: 3 × 10
#>   from  to    delta estimate std.error    df conf.low conf.high statistic
#>   <chr> <chr> <dbl>    <dbl>     <dbl> <dbl>    <dbl>     <dbl>     <dbl>
#> 1 MVPA  SB       10     1.32    0.0705   Inf     1.19      1.47      5.20
#> 2 MVPA  SB       20     2.17    0.325    Inf     1.61      2.91      5.16
#> 3 MVPA  SB       30    NA      NA         NA    NA        NA        NA   
#> # ℹ 1 more variable: p.value <dbl>
```

The reversed estimate is not the reciprocal of the original. Moving 10
minutes into moderate-to-vigorous activity gave an odds ratio of 0.822,
whose reciprocal is about 1.22, whereas moving 10 minutes out of it
gives 1.32, because removing 10 of 28.1 minutes is a larger relative
change than adding 10 to them. Dumuid et al. (2020) note that estimated
differences need not be symmetric when the same duration is added to or
taken from a part, and Dumuid, Stanford, Pedišić, et al. (2018) report
such asymmetry for adiposity.

You can also reverse the transfer with a negative `delta`. The estimates
for the original `from` and `to` with negative magnitudes are those of
the swapped call.

### All pairs at once

To use every ordered pair of parts, omit `from` and `to`. This gives
`D (D - 1)` reallocations for a composition of `D` parts.

``` r

predict_realloc(object = fit_ob, delta = 15)
#> # A tibble: 12 × 10
#>    from  to    delta estimate std.error    df conf.low conf.high statistic
#>    <chr> <chr> <dbl>    <dbl>     <dbl> <dbl>    <dbl>     <dbl>     <dbl>
#>  1 LPA   MVPA     15    0.775    0.0422   Inf    0.696     0.862    -4.69 
#>  2 SB    MVPA     15    0.758    0.0400   Inf    0.683     0.840    -5.26 
#>  3 SPT   MVPA     15    0.782    0.0427   Inf    0.702     0.870    -4.51 
#>  4 MVPA  LPA      15    1.58     0.146    Inf    1.32      1.89      4.92 
#>  5 SB    LPA      15    0.978    0.0230   Inf    0.934     1.02     -0.942
#>  6 SPT   LPA      15    1.01     0.0250   Inf    0.961     1.06      0.368
#>  7 MVPA  SB       15    1.61     0.149    Inf    1.35      1.93      5.18 
#>  8 LPA   SB       15    1.02     0.0250   Inf    0.974     1.07      0.888
#>  9 SPT   SB       15    1.03     0.0205   Inf    0.992     1.07      1.55 
#> 10 MVPA  SPT      15    1.56     0.147    Inf    1.30      1.88      4.77 
#> 11 LPA   SPT      15    0.991    0.0254   Inf    0.943     1.04     -0.342
#> 12 SB    SPT      15    0.970    0.0191   Inf    0.933     1.01     -1.57 
#> # ℹ 1 more variable: p.value <dbl>
```

The odds ratios farthest from 1 are in the six transfers that involve
moderate-to-vigorous activity.

### One part against all the others

To move minutes between one focal part and all the others, set
`scheme = "one-to-all"`. The minutes are shared among the other parts in
proportion to their reference values, so that those parts keep their
ratios to one another.

``` r

r_rest <- predict_realloc(object = fit_ob, to = "MVPA", delta = c(10, 20, 30),
                          scheme = "one-to-all")
summary(object = r_rest)
#> Time-reallocation effects: Odds ratio versus the reference composition
#>   3 reallocations of 10, 20, 30; 95% confidence intervals
#> 
#> # A tibble: 3 × 9
#>   to    from   to_label     from_label delta estimate conf.low conf.high p.value
#>   <chr> <chr>  <chr>        <chr>      <dbl>    <dbl>    <dbl>     <dbl>   <dbl>
#> 1 MVPA  (rest) Moderate-to… All other…    10    0.830    0.771     0.892 4.89e-7
#> 2 MVPA  (rest) Moderate-to… All other…    20    0.718    0.632     0.817 4.89e-7
#> 3 MVPA  (rest) Moderate-to… All other…    30    0.639    0.536     0.761 4.89e-7
```

The odds ratio for 10 minutes is close to the 0.822 obtained when all 10
minutes come from sedentary behaviour.

To make the focal part give minutes to the others, use `from = "MVPA"`
in place of `to`. If neither is given, every part in turn is the focal
part receiving minutes.

### Custom reallocations

To reallocate time among three or more parts at once, give `realloc` a
named vector of the minutes each part gains or loses. The changes must
sum to zero, and parts left out keep their reference values. `realloc`
cannot be combined with `from`, `to` or `scheme`.

``` r

r_custom <- predict_realloc(object = fit_ob,
                            realloc = c(SPT = 30, SB = -15, LPA = -15))
summary(object = r_custom)
#> Time-reallocation effects: Odds ratio versus the reference composition
#>   1 reallocation of 30; 95% confidence intervals
#> 
#> # A tibble: 1 × 9
#>   to      from     to_label from_label delta estimate conf.low conf.high p.value
#>   <chr>   <chr>    <chr>    <chr>      <dbl>    <dbl>    <dbl>     <dbl>   <dbl>
#> 1 SPT +30 LPA -15… Sleep +… Light PA …    30    0.961    0.892      1.04   0.302
```

The scenario is named in the row by its signed changes, `SPT +30` under
`to` and `LPA -15, SB -15` under `from`, and `delta` is the 30 minutes
moved in total. To estimate several scenarios in one call, give a list
of such vectors. Each scenario gets one row.

### Invalid transfers

A transfer larger than the donor’s reference minutes drives that part to
zero or below. To choose how such an invalid transfer is handled, set
`on_invalid`. Its default, `"na"`, keeps the transfer without an
estimate, and `"omit"` and `"error"` are the alternatives.

With `"error"`, the call halts at the first invalid transfer and returns
no result.

With `"omit"`, the invalid transfer is dropped and the valid one
returned.

``` r

predict_realloc(object = fit_ob, from = "MVPA", to = "SB", delta = c(15, 60),
                on_invalid = "omit")
#> # A tibble: 1 × 10
#>   from  to    delta estimate std.error    df conf.low conf.high statistic
#>   <chr> <chr> <dbl>    <dbl>     <dbl> <dbl>    <dbl>     <dbl>     <dbl>
#> 1 MVPA  SB       15     1.61     0.149   Inf     1.35      1.93      5.18
#> # ℹ 1 more variable: p.value <dbl>
```

For a sweep over several magnitudes or pairs, use the default or
`"omit"`. A single invalid transfer then does not halt the whole call.

A transfer that leaves a part close to zero passes the check, although
the reallocated composition then lies near the edge of the simplex,
where predicted outcomes tend towards infinity as a part approaches zero
(Dumuid et al., 2020).

### Relative and absolute output

To choose what each row reports, set `type`. With `type = "absolute"`,
each row is instead the model’s prediction at that day, on the scale of
the outcome.

``` r

predict_realloc(object = fit_bmi, from = "SB", to = "MVPA", delta = c(10, 20, 30))
#> # A tibble: 3 × 10
#>   from  to    delta estimate std.error    df conf.low conf.high statistic
#>   <chr> <chr> <dbl>    <dbl>     <dbl> <dbl>    <dbl>     <dbl>     <dbl>
#> 1 SB    MVPA     10   -0.209    0.0135   296   -0.236    -0.182     -15.5
#> 2 SB    MVPA     20   -0.374    0.0240   296   -0.422    -0.327     -15.6
#> 3 SB    MVPA     30   -0.513    0.0326   296   -0.577    -0.449     -15.8
#> # ℹ 1 more variable: p.value <dbl>
```

Moving 10 minutes from sedentary behaviour to moderate-to-vigorous
activity is associated with a z-BMI 0.209 units lower.

``` r

predict_realloc(object = fit_bmi, from = "SB", to = "MVPA",
                delta = c(10, 20, 30), type = "absolute")
#> # A tibble: 4 × 8
#>   from  to    delta estimate std.error    df conf.low conf.high
#>   <chr> <chr> <dbl>    <dbl>     <dbl> <dbl>    <dbl>     <dbl>
#> 1 NA    NA        0  -0.0228    0.0521   296   -0.125    0.0797
#> 2 SB    MVPA     10  -0.232     0.0539   296   -0.338   -0.126 
#> 3 SB    MVPA     20  -0.397     0.0574   296   -0.510   -0.284 
#> 4 SB    MVPA     30  -0.536     0.0615   296   -0.657   -0.415
```

The first row is the reference day, which has no donor or receiver and a
`delta` of 0. Subtracting the reference row from a reallocated row gives
the mean difference of the relative output. For the logistic regression,
the absolute output is the predicted probability of obesity.

``` r

predict_realloc(object = fit_ob, from = "SB", to = "MVPA",
                delta = c(10, 20, 30), type = "absolute")
#> # A tibble: 4 × 8
#>   from  to    delta estimate std.error    df conf.low conf.high
#>   <chr> <chr> <dbl>    <dbl>     <dbl> <dbl>    <dbl>     <dbl>
#> 1 NA    NA        0    0.258    0.0271   Inf    0.209     0.315
#> 2 SB    MVPA     10    0.222    0.0268   Inf    0.174     0.279
#> 3 SB    MVPA     20    0.197    0.0270   Inf    0.149     0.255
#> 4 SB    MVPA     30    0.177    0.0271   Inf    0.130     0.237
```

The odds ratios of the complete call are ratios of the odds at these
probabilities.

### A chosen reference day

To anchor the comparison at a chosen day instead (e.g., a guideline
schedule or a policy target), give the `reference` argument one positive
value for each part. The values must sum to the total of the
composition, 1440 minutes here.

``` r

predict_realloc(object = fit_ob, from = "SB", to = "MVPA",
                delta = c(10, 20, 30),
                reference = c(MVPA = 30, LPA = 240, SB = 660, SPT = 510))
#> # A tibble: 3 × 10
#>   from  to    delta estimate std.error    df conf.low conf.high statistic
#>   <chr> <chr> <dbl>    <dbl>     <dbl> <dbl>    <dbl>     <dbl>     <dbl>
#> 1 SB    MVPA     10    0.830    0.0294   Inf    0.774     0.889     -5.26
#> 2 SB    MVPA     20    0.715    0.0454   Inf    0.632     0.810     -5.28
#> 3 SB    MVPA     30    0.632    0.0548   Inf    0.533     0.749     -5.30
#> # ℹ 1 more variable: p.value <dbl>
```

The odds ratio for 10 minutes is now 0.830, against 0.822 at the closed
geometric mean. Because the same transfer is a different relative change
on a different day, the estimate depends on the day it starts from.
Dumuid, Stanford, Pedišić, et al. (2018) show such differences between
reference compositions for boys and girls.

### Covariates

Neither fit so far includes covariates. When a model adjusts for
covariates without composition-by-covariate interactions, the covariates
cancel in the contrast, so a relative estimate is the same at every
covariate value, and an absolute prediction is made with the covariates
averaged over, as in the emmeans grid step. When the model includes an
interaction between the composition and a covariate, the contrast
differs between values of that covariate, and the call ends with an
error unless you fix those values through `at` or split them through
`by`. `at` takes a list of covariate values such as
`at = list(sex = "female")`. For this case, see the help page
[`?predict_realloc`](https://marcusvvlopes.github.io/timecoda/reference/predict_realloc.md).

In `fit_sex`, sex is added as a `moderator` of the composition, and the
reallocation is estimated once per sex. For moderators, see
[`vignette("advanced-models", package = "timecoda")`](https://marcusvvlopes.github.io/timecoda/articles/advanced-models.md).

``` r

fit_sex <- comp_iv(comp = cd, outcome = "obesity", model = "glm",
                   family = "binomial", moderator = "sex")
r_sex <- predict_realloc(object = fit_sex, from = "SB", to = "MVPA",
                         delta = c(10, 20, 30), by = "sex")
summary(object = r_sex)
#> Time-reallocation effects: Odds ratio versus the reference composition
#>   6 reallocations of 10, 20, 30; 95% confidence intervals
#>   Estimated per level of sex
#> 
#> # A tibble: 6 × 10
#>   to    from  sex    to_label       from_label delta estimate conf.low conf.high
#>   <chr> <chr> <chr>  <chr>          <chr>      <dbl>    <dbl>    <dbl>     <dbl>
#> 1 MVPA  SB    female Moderate-to-v… Sedentary…    10    0.812    0.734     0.899
#> 2 MVPA  SB    female Moderate-to-v… Sedentary…    20    0.689    0.575     0.825
#> 3 MVPA  SB    female Moderate-to-v… Sedentary…    30    0.599    0.468     0.767
#> 4 MVPA  SB    male   Moderate-to-v… Sedentary…    10    0.832    0.748     0.926
#> 5 MVPA  SB    male   Moderate-to-v… Sedentary…    20    0.722    0.597     0.873
#> 6 MVPA  SB    male   Moderate-to-v… Sedentary…    30    0.643    0.496     0.832
#> # ℹ 1 more variable: p.value <dbl>
```

On a fit with a factor `moderator`, `by` defaults to the moderator, so
the same rows are returned without naming it.
[`plot_realloc()`](https://marcusvvlopes.github.io/timecoda/reference/plot_realloc.md)
passes `by` to
[`predict_realloc()`](https://marcusvvlopes.github.io/timecoda/reference/predict_realloc.md)
and stacks one panel grid per level, and the grids share the legend and
the vertical axis, so that you can compare a curve between the sexes
panel by panel.

``` r

plot_realloc(object = fit_sex, delta_max = 60, n_points = 20, by = "sex")
```

![Grouped reallocation figure for obesity with sex as moderator: two
lower-triangular grids of six panels stacked as rows, the first titled
sex: female and the second sex: male, each panel showing the odds ratio
against the minutes reallocated between a pair of parts with a
confidence ribbon and a dashed line at 1, on one vertical axis shared by
both
grids.](predicting-reallocations_files/figure-html/covariates-plot-1.png)

### The ILR basis

A different ILR basis rotates the coordinates and the coefficients
together, and for a model linear in the coordinates, like the two fits
here, or one holding their full quadratic, every prediction is
unchanged.
[`predict_realloc()`](https://marcusvvlopes.github.io/timecoda/reference/predict_realloc.md)
evaluates every contrast in the single basis stored in the fit. A
squares-only or [`poly()`](https://rdrr.io/r/stats/poly.html) term
changes with the order of the parts (see
[`?ilr`](https://marcusvvlopes.github.io/timecoda/reference/ilr.md)).

## Part 3. Publication output

### Publication table

The table returned by [`summary()`](https://rdrr.io/r/base/summary.html)
keeps full numeric precision. To build the publication table from it,
use
[`pub_ready()`](https://marcusvvlopes.github.io/timecoda/reference/pub_ready.md),
which uses the `gt` package. `gt` is a suggested package, and
[`pub_ready()`](https://marcusvvlopes.github.io/timecoda/reference/pub_ready.md)
is the only function that requires it.

``` r

r_all <- predict_realloc(object = fit_ob, delta = c(10, 20, 30))
```

``` r

pub_ready(x = summary(object = r_all), digits = 2)$effects
```

[TABLE]

[`pub_ready()`](https://marcusvvlopes.github.io/timecoda/reference/pub_ready.md)
returns a list of tables, and at the console the bare call
`pub_ready(x = summary(object = r_all), digits = 2)` displays every
table in the list.

The dashes in the 30-minute group are the three invalid transfers out of
moderate-to-vigorous activity, which the default `on_invalid = "na"`
keeps as missing estimates.

To change the threshold for the asterisk, set `p_flag` in
[`pub_ready()`](https://marcusvvlopes.github.io/timecoda/reference/pub_ready.md),
here to `1e-6`. Only the asterisks differ from those of the first table.

``` r

pub_ready(x = summary(object = r_all), digits = 2, p_flag = 1e-6)$effects
```

[TABLE]

### Transfer figure

[`plot_realloc()`](https://marcusvvlopes.github.io/timecoda/reference/plot_realloc.md)
draws the estimates as continuous curves. It sweeps each pair of parts
in both directions up to `delta_max` minutes and evaluates each curve at
`n_points` magnitudes across the two directions. Twenty is used here to
keep this example quick, and the default of 41 draws smoother curves.

Each of the six panels in the lower triangle of the grid holds one
unordered pair of parts. The minutes moved between the two parts are
plotted as a positive magnitude on both sides of zero, with time moving
to the part named first in the panel’s title to the right of zero and to
the second part to the left of zero. The vertical axis shows the odds
ratio on a logarithmic scale shared by every panel, and a dashed line
marks 1 for no difference. The shaded ribbon is the 95% confidence
interval.

``` r

plot_realloc(object = fit_ob, delta_max = 60, n_points = 20)
```

![Lower-triangular grid of six panels, one per pair of parts, each
showing the odds ratio of obesity against the minutes reallocated
between the two parts in either direction, with a confidence ribbon, a
dashed line at 1, and colour marking the part gaining
time.](predicting-reallocations_files/figure-html/plot-1.png)

In the panels with moderate-to-vigorous activity, the half in which that
part gives time climbs more steeply than the other half falls, which is
the asymmetry found when the transfer was reversed. That half stops at
24 minutes, the largest evaluated magnitude that leaves the part
positive at its reference value of 28.1 minutes.

[`plot_realloc()`](https://marcusvvlopes.github.io/timecoda/reference/plot_realloc.md)
also returns the panels, so you can draw one of them on its own or
assemble it with others. The returned list holds `panels`, named by the
parts of each panel, and `figure`, the assembled figure.

``` r

p <- plot_realloc(object = fit_ob, delta_max = 60, n_points = 20)
```

![A single reallocation panel for obesity: the odds ratio against the
minutes reallocated between sedentary behaviour and moderate-to-vigorous
activity, with a confidence ribbon and a dashed line at
1.](predicting-reallocations_files/figure-html/plot-panel-1.png)

``` r

p$panels[["MVPA.SB"]]
```

![A single reallocation panel for obesity: the odds ratio against the
minutes reallocated between sedentary behaviour and moderate-to-vigorous
activity, with a confidence ribbon and a dashed line at
1.](predicting-reallocations_files/figure-html/plot-panel-2.png)

To hide the panels of a part, name it in `exclude`.

``` r

plot_realloc(object = fit_ob, delta_max = 60, n_points = 20, exclude = "MVPA")
```

![Reallocation figure for obesity without the moderate-to-vigorous
activity panels: three panels showing the odds ratio against the minutes
reallocated between light activity, sedentary behaviour and sleep, with
confidence ribbons and a dashed line at
1.](predicting-reallocations_files/figure-html/plot-exclude-1.png)

The excluded part stays in the composition, so the filter changes only
the display. With `scheme = "one-to-all"`,
[`plot_realloc()`](https://marcusvvlopes.github.io/timecoda/reference/plot_realloc.md)
draws one panel per part, each showing that part against all the others,
and an excluded part still takes its share of each reallocation.

For the linear model, the vertical axis shows the mean difference in
z-BMI on an ordinary scale, and the dashed reference line is drawn at 0.

``` r

plot_realloc(object = fit_bmi, delta_max = 60, n_points = 20)
```

![Reallocation figure for the linear model: six panels showing the mean
difference in z-BMI against the minutes reallocated between each pair of
parts, with a confidence ribbon and a dashed reference line at
0.](predicting-reallocations_files/figure-html/plot-bmi-1.png)

For z-BMI, the ribbons of the transfers between sedentary behaviour and
either light activity or sleep lie wholly above or below the dashed line
at reallocations away from zero, whereas the ribbons of the transfers
between light activity and sleep include it.

## References

Dumuid, D., Pedišić, Ž., Palarea-Albaladejo, J., Martín-Fernández, J.
A., Hron, K., & Olds, T. (2020). Compositional data analysis in time-use
epidemiology: What, why, how. *International Journal of Environmental
Research and Public Health*, 17, 2220.

Dumuid, D., Pedišić, Ž., Stanford, T. E., Martín-Fernández, J. A., Hron,
K., Maher, C. A., Lewis, L. K., & Olds, T. (2019). The compositional
isotemporal substitution model: A method for estimating changes in a
health outcome for reallocation of time between sleep, physical activity
and sedentary behaviour. *Statistical Methods in Medical Research*, 28,
846-857.

Dumuid, D., Stanford, T. E., Pedišić, Ž., Maher, C., Lewis, L. K.,
Martín-Fernández, J. A., Katzmarzyk, P. T., Chaput, J. P., Fogelholm,
M., Standage, M., Tremblay, M. S., & Olds, T. (2018). Adiposity and the
isotemporal substitution of physical activity, sedentary time and sleep
among school-aged children: A compositional data analysis approach. *BMC
Public Health*, 18, 311.

## See also

- [`vignette("composition-predictor", package = "timecoda")`](https://marcusvvlopes.github.io/timecoda/articles/composition-predictor.md):
  fitting the models whose reallocations this vignette estimates.
- [`vignette("timecoda", package = "timecoda")`](https://marcusvvlopes.github.io/timecoda/articles/timecoda.md):
  the complete workflow, from declaring a composition to the
  reallocation figure.
