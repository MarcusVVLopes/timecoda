# Predict time-reallocation effects

Estimates the effect of reallocating time between the parts of a
composition, from a
[`comp_iv()`](https://marcusvvlopes.github.io/timecoda/reference/comp_iv.md)
fit. Each estimate is the contrast of the model's linear predictor
between a reallocated composition and a reference composition, computed
by `emmeans`. It is back-transformed to the response scale: an odds,
hazard, rate or risk ratio for a ratio-scale model, and a mean
difference for a linear one.

## Usage

``` r
predict_realloc(
  object,
  delta,
  from = NULL,
  to = NULL,
  scheme = c("one-to-one", "one-to-all"),
  realloc = NULL,
  reference = NULL,
  type = c("relative", "absolute"),
  at = NULL,
  by = NULL,
  conf_level = 0.95,
  on_invalid = c("na", "omit", "error")
)
```

## Arguments

- object:

  A `comp_iv` / `timecoda_fit` from
  [`comp_iv()`](https://marcusvvlopes.github.io/timecoda/reference/comp_iv.md).

- delta:

  Time to reallocate, in the composition's closure units: a number or
  vector of numbers, each applied to every `from`/`to` pair (a negative
  value moves time the other way). Ignored in custom mode (when
  `realloc` is supplied), where each scenario has its own signed
  changes.

- from, to:

  Character vectors of part names defining the reallocations. For
  `scheme = "one-to-one"` (the default), they name the donor (`from`)
  and receiver (`to`) parts and default to all ordered pairs. If given,
  they are matched elementwise (scalars recycled). For
  `scheme = "one-to-all"`, supply only one of them. Each part it names
  is a focal part in turn, and the other side is implicitly all
  remaining parts (see `scheme`).

- scheme:

  How time is reallocated. `"one-to-one"` (default) moves `delta`
  closure units directly between the `from` and `to` parts.
  `"one-to-all"` moves `delta` closure units between a single focal part
  \\f\\ and all the remaining parts, in proportion to their values in
  the reference composition \\r\\. Remaining part \\j\\ changes by \\\mp
  \delta \\ r_j / \sum\_{k \ne f} r_k\\, and the remaining parts keep
  their ratios to one another. With `to = "MVPA"` the focal part gains
  time from the rest, and with `from = "MVPA"` it gives time to the
  rest. With neither, every part is used as the focal receiver in turn.
  The implicit "all remaining parts" side is labelled `"(rest)"` in the
  `from`/`to` output column.

- realloc:

  Optional custom multi-part reallocation, in the style of Dumuid et
  al. (2022) "Your best day". In place of a directed `from`/`to` pair,
  give a signed change per part, in the composition's closure units. The
  changes are added to the reference composition, which is then
  re-closed. For example, `c(SPT = 30, SB = -15, LPA = -15)` moves 30
  minutes into sleep, taken 15 each from `SB` and `LPA`. Accepts (1) a
  named numeric vector (one scenario), (2) a list of named numeric
  vectors (several scenarios), or (3) a data frame with one column per
  part and one row per scenario. Every scenario must name known parts,
  sum to zero (to within `1e-8`; a closed reallocation), and hold at
  least one positive and one negative entry. Parts left unnamed are held
  at their reference value. Custom mode is mutually exclusive with
  `from`, `to` and `scheme`. `type`, `at`, `by` and `on_invalid` apply
  as for a directed reallocation, and a multinomial fit still gives one
  row set per outcome category.

- reference:

  Optional named numeric vector (or one-row data frame) giving the
  composition to anchor the reallocations at: one strictly positive
  value per part, summing to the fit's composition total. The default
  (`NULL`) uses the sample's closed geometric mean, stored in the fit.

- type:

  `"relative"` (default) contrasts each reallocated composition against
  the reference; `"absolute"` returns the predicted response at each
  composition. `"cox"`, `"coxme"`, `"clm"`, `"clmm"` and `"mblogit"`
  fits support `"relative"` only.

- at:

  Optional named list of covariate values at which to evaluate. Each
  covariate that interacts with the composition must be pinned by `at`
  or split by `by`. Covariates not in `at` are averaged over in the
  usual `emmeans` way.

- by:

  Optional name of a factor or character covariate, to compute the
  reallocation per group. The contrast is evaluated once at each level
  present in the fitted data. The result gains a column named after the
  covariate, holding the level (e.g., one odds ratio per `sex` for a
  moderated effect). When the fit was built with a factor `moderator`
  (see
  [`comp_iv()`](https://marcusvvlopes.github.io/timecoda/reference/comp_iv.md))
  and `at` does not pin it, `by` defaults to that moderator. `at` and
  `by` cannot name the same covariate. Pin numeric covariates with `at`
  instead.

- conf_level:

  Confidence level for the intervals.

- on_invalid:

  What to do when a reallocation drives a part to zero or below: `"na"`
  (default) keeps the rows but returns `NA` estimates for them; `"omit"`
  drops those rows; `"error"` aborts.

## Value

A `realloc` object, a tibble with columns `from`, `to`, `delta`,
`estimate`, `std.error`, `df`, `conf.low`, `conf.high`, `statistic` and
`p.value` (`statistic` and `p.value` only for `type = "relative"`).
`std.error` is on the scale of `estimate`. For a ratio effect it is the
standard error of the ratio itself, as `emmeans` reports it on the
response scale; the standard error of the log ratio is not returned. For
`type = "absolute"`, the first row is the reference composition (`from`,
`to` = `NA`, `delta` = 0). A `by` grouping adds a column named after the
covariate, and a multinomial fit adds `response_level`, one row set per
non-reference outcome category. In custom mode, each scenario's `from`
and `to` are compact, comma-separated signed part strings (e.g.,
`from = "LPA -15, SB -15"`, `to = "SPT +30"`), and `delta` is the
scenario's total time moved, in closure units.
[summary()](https://marcusvvlopes.github.io/timecoda/reference/summary.realloc.md)
builds the summary table from it. The fit metadata that method needs is
stored as attributes, among them `units`, the composition's unit label
(`NULL` when none was declared).

## Details

The reference and every reallocated composition are converted to the
fit's ILR coordinates and evaluated as their own `emmeans` reference
grids, averaged over the covariates. Each reallocated row is contrasted
with the reference row on the link scale,
\\\eta(x\_{\mathrm{realloc}}) - \eta(x\_{\mathrm{ref}})\\, with a
delta-method standard error, and back-transformed with
`type = "response"`. Every row is its own Wald test, with no
multiplicity adjustment across rows. The estimates do not depend on the
ILR basis when the model is linear in the coordinates or holds their
full quadratic
([`ilr()`](https://marcusvvlopes.github.io/timecoda/reference/ilr.md)).
A row with `delta = 0` contrasts the reference with itself and is
reported as the reference value (1 for a ratio, 0 for a difference),
with a zero standard error, `NA` degrees of freedom and no test.

Without composition-by-covariate interactions, the covariates cancel in
the contrast and a relative estimate does not depend on them. With such
interactions the contrast is covariate-conditional. Every interacting
covariate must then be pinned with `at`, or the contrast evaluated per
level with `by`; the call aborts naming any covariate left unpinned. On
a fit built with `include_total = TRUE`, every composition is evaluated
with `t_total` held at the fit's `total_ref`. A reallocation then
changes the shares at a constant total, and the total cancels from a
relative contrast.

`type = "absolute"` returns the predicted response at each composition
instead of a contrast. It is not available for `"cox"`, `"coxme"`,
`"clm"`, `"clmm"` or `"mblogit"` fits. These have no single
absolute-response scale: Cox leaves its baseline hazard unestimated, and
ordinal and multinomial predictions are per outcome category.

The `df` column holds the degrees of freedom `emmeans` computes for the
fit, and the per-part table of
[summary()](https://marcusvvlopes.github.io/timecoda/reference/summary.comp_iv.md)
uses the same degrees of freedom. For an `"lmer"` fit, the default is
the Kenward-Roger approximation through pbkrtest. The Satterthwaite
approximation through lmerTest is used when it is selected with
`emmeans::emm_options(lmer.df =)`, or when pbkrtest cannot be used
(absent, disabled, or more observations than
`emmeans::emm_options(pbkrtest.limit =)`). Both packages are optional,
and when neither can be used the degrees of freedom are `Inf`
(asymptotic). For `"glmer"`, `"cox"`, `"coxme"` and `"rlmer"` fits the
degrees of freedom are `Inf` (robustlmm reports none), and for `"lmrob"`
fits they are the residual degrees of freedom.

A reallocation that drives any part to zero or below has no composition
to evaluate. By default (`on_invalid = "na"`), such a row is kept with
`NA` estimates. A sweep over many pairs or `delta`s then still returns
every feasible row, and the printed summary gives the number of rows
with no estimate. `"omit"` drops those rows instead, and `"error"`
aborts the whole call on the first one.

## The reference composition and "isotemporal substitution"

The `reference` argument sets the composition every reallocation starts
from, and the contrast is computed the same way at any reference. The
default, the cohort's geometric-mean day, suits population-typical
questions. A day supplied by the user, such as a guideline schedule or
one person's day, suits scenario questions. In the compositional
time-use literature this prediction is *compositional isotemporal
substitution* (Dumuid et al. 2019). The traditional leave-one-out
isotemporal substitution model (Mekary et al. 2009) is a different
regression parameterization, and no choice of `reference` reproduces it.
[`vignette("predicting-reallocations")`](https://marcusvvlopes.github.io/timecoda/articles/predicting-reallocations.md)
compares the two.

## References

Dumuid, D., Pedišić, Ž., Stanford, T. E., et al. (2019). The
compositional isotemporal substitution model: A method for estimating
changes in a health outcome for reallocation of time between sleep,
physical activity and sedentary behaviour. *Statistical Methods in
Medical Research*, 28(3), 846-857.
[doi:10.1177/0962280217737805](https://doi.org/10.1177/0962280217737805)

Mekary, R. A., Willett, W. C., Hu, F. B., & Ding, E. L. (2009).
Isotemporal substitution paradigm for physical activity epidemiology and
weight change. *American Journal of Epidemiology*, 170(4), 519-527.
[doi:10.1093/aje/kwp163](https://doi.org/10.1093/aje/kwp163)

Dumuid, D., Olds, T., Wake, M., et al. (2022). "Your best day": An
interactive app to translate how time reallocations within a 24-hour day
are associated with health measures. *PLOS ONE*, 17(9), e0272343.
[doi:10.1371/journal.pone.0272343](https://doi.org/10.1371/journal.pone.0272343)

Kenward, M. G., & Roger, J. H. (1997). Small sample inference for fixed
effects from restricted maximum likelihood. *Biometrics*, 53(3),
983-997. [doi:10.2307/2533558](https://doi.org/10.2307/2533558)

## See also

[`comp_iv()`](https://marcusvvlopes.github.io/timecoda/reference/comp_iv.md)
to fit the model this consumes;
[`summary.realloc()`](https://marcusvvlopes.github.io/timecoda/reference/summary.realloc.md)
for the summary table;
[`plot_realloc()`](https://marcusvvlopes.github.io/timecoda/reference/plot_realloc.md)
for the continuous sweep figure;
[`optimal_composition()`](https://marcusvvlopes.github.io/timecoda/reference/optimal_composition.md)
to search for the best-predicted composition instead of a directed
transfer.

Other reallocation and prediction:
[`optimal_composition()`](https://marcusvvlopes.github.io/timecoda/reference/optimal_composition.md),
[`plot_realloc()`](https://marcusvvlopes.github.io/timecoda/reference/plot_realloc.md),
[`summary.realloc()`](https://marcusvvlopes.github.io/timecoda/reference/summary.realloc.md)

## Examples

``` r
cd <- declare_composition(timeuse, c("MVPA", "LPA", "SB", "SPT"),
                          closure = 1440)
fit <- comp_iv(cd, outcome = "obesity", model = "glm", family = "binomial")
predict_realloc(fit, delta = 30, from = "SB", to = "MVPA")
#> # A tibble: 1 × 10
#>   from  to    delta estimate std.error    df conf.low conf.high statistic
#>   <chr> <chr> <dbl>    <dbl>     <dbl> <dbl>    <dbl>     <dbl>     <dbl>
#> 1 SB    MVPA     30    0.620    0.0561   Inf    0.519     0.740     -5.29
#> # ℹ 1 more variable: p.value <dbl>

# Proportional reallocation: add 30 min to MVPA, taken from the rest of the
# day in proportion to its composition.
predict_realloc(fit, delta = 30, scheme = "one-to-all", to = "MVPA")
#> # A tibble: 1 × 10
#>   from   to    delta estimate std.error    df conf.low conf.high statistic
#>   <chr>  <chr> <dbl>    <dbl>     <dbl> <dbl>    <dbl>     <dbl>     <dbl>
#> 1 (rest) MVPA     30    0.639    0.0569   Inf    0.536     0.761     -5.03
#> # ℹ 1 more variable: p.value <dbl>

# Custom multi-part reallocation ("best-day" style): add 30 min to SPT,
# taken 15 each from SB and LPA.
predict_realloc(fit, realloc = c(SPT = 30, SB = -15, LPA = -15))
#> # A tibble: 1 × 10
#>   from         to    delta estimate std.error    df conf.low conf.high statistic
#>   <chr>        <chr> <dbl>    <dbl>     <dbl> <dbl>    <dbl>     <dbl>     <dbl>
#> 1 LPA -15, SB… SPT …    30    0.961    0.0367   Inf    0.892      1.04     -1.03
#> # ℹ 1 more variable: p.value <dbl>

# Anchor the same transfer at a chosen day instead of the cohort mean.
predict_realloc(fit, delta = 30, from = "SB", to = "MVPA",
                reference = c(MVPA = 30, LPA = 240, SB = 660, SPT = 510))
#> # A tibble: 1 × 10
#>   from  to    delta estimate std.error    df conf.low conf.high statistic
#>   <chr> <chr> <dbl>    <dbl>     <dbl> <dbl>    <dbl>     <dbl>     <dbl>
#> 1 SB    MVPA     30    0.632    0.0548   Inf    0.533     0.749     -5.30
#> # ℹ 1 more variable: p.value <dbl>
```
