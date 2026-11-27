# Summarize time-reallocation effects

Tabulates the reallocations of a
[`predict_realloc()`](https://marcusvvlopes.github.io/timecoda/reference/predict_realloc.md)
result as a tidy tibble at full numeric precision: one row per receiving
part, donor part and reallocation magnitude (`delta`), with the parts'
display labels. Printing the result shows a text report, as with
[`summary.lm()`](https://rdrr.io/r/stats/summary.lm.html).
[`pub_ready()`](https://marcusvvlopes.github.io/timecoda/reference/pub_ready.md)
builds the wide publication table from it, with the effect estimate and
its confidence interval under a column spanner for each `delta`.

## Usage

``` r
# S3 method for class 'realloc'
summary(object, labels = attr(object, "labels"), ...)
```

## Arguments

- object:

  A `realloc` object from
  [`predict_realloc()`](https://marcusvvlopes.github.io/timecoda/reference/predict_realloc.md).
  An infeasible transfer (a row with `NA` estimates, the
  [`predict_realloc()`](https://marcusvvlopes.github.io/timecoda/reference/predict_realloc.md)
  default) is kept, and shows as `-` in the
  [`pub_ready()`](https://marcusvvlopes.github.io/timecoda/reference/pub_ready.md)
  table.

- labels:

  Optional named character vector mapping raw part names to display
  labels (e.g., `c(MVPA = "Exercise", SB = "Sedentary")`). Parts not
  named keep their raw names, and the implicit "all remaining parts"
  side of a one-to-all reallocation is shown as "All other parts".
  Defaults to the labels of the
  [`declare_composition()`](https://marcusvvlopes.github.io/timecoda/reference/declare_composition.md)
  object behind the fit (if any); pass a vector here to override them.

- ...:

  Unused; present for S3 consistency. Passing `display`, `digits`,
  `p_flag` or `conf_level` is an error: rounding and the significance
  flag are set in
  [`pub_ready()`](https://marcusvvlopes.github.io/timecoda/reference/pub_ready.md),
  and the confidence level in
  [`predict_realloc()`](https://marcusvvlopes.github.io/timecoda/reference/predict_realloc.md).

## Value

An object of class `summary_realloc`, printed as a text report by
[`print.summary_realloc()`](https://marcusvvlopes.github.io/timecoda/reference/print.summary_realloc.md)
and passed to
[`pub_ready()`](https://marcusvvlopes.github.io/timecoda/reference/pub_ready.md)
for a `gt` table. A list with elements:

- `table`:

  A tibble with one row per reallocation: `to`, `from`, the `by`
  grouping column and `response_level` when present, the display labels
  `to_label` and `from_label`, then `delta`, `estimate`, `conf.low`,
  `conf.high` and, for relative output, `p.value`. The reference row of
  an absolute prediction is dropped; a zero reallocation (`delta = 0`)
  is kept.

- `effect_label`, `type`, `conf_level`, `parts`, `by`:

  Kept for
  [`print.summary_realloc()`](https://marcusvvlopes.github.io/timecoda/reference/print.summary_realloc.md)
  and
  [`pub_ready()`](https://marcusvvlopes.github.io/timecoda/reference/pub_ready.md).

- `units`:

  The composition's unit label, shown beside each `delta` in the
  [`pub_ready()`](https://marcusvvlopes.github.io/timecoda/reference/pub_ready.md)
  table; `NULL` when none was declared.

## See also

[`predict_realloc()`](https://marcusvvlopes.github.io/timecoda/reference/predict_realloc.md)
to build the `realloc` object;
[`plot_realloc()`](https://marcusvvlopes.github.io/timecoda/reference/plot_realloc.md)
for the continuous sweep figure;
[`pub_ready()`](https://marcusvvlopes.github.io/timecoda/reference/pub_ready.md)
for the publication `gt` table built from this summary.

Other reallocation and prediction:
[`optimal_composition()`](https://marcusvvlopes.github.io/timecoda/reference/optimal_composition.md),
[`plot_realloc()`](https://marcusvvlopes.github.io/timecoda/reference/plot_realloc.md),
[`predict_realloc()`](https://marcusvvlopes.github.io/timecoda/reference/predict_realloc.md)

## Examples

``` r
cd <- declare_composition(timeuse, c("MVPA", "LPA", "SB", "SPT"),
                          closure = 1440)
fit <- comp_iv(cd, outcome = "obesity", model = "glm", family = "binomial")
r <- predict_realloc(fit, delta = c(10, 20, 30), on_invalid = "na")
s <- summary(r, labels = c(MVPA = "Exercise"))
s$table
#> # A tibble: 36 × 9
#>    to    from  to_label from_label delta estimate conf.low conf.high     p.value
#>    <chr> <chr> <chr>    <chr>      <dbl>    <dbl>    <dbl>     <dbl>       <dbl>
#>  1 MVPA  LPA   Exercise LPA           10    0.834    0.773     0.899 0.00000225 
#>  2 MVPA  LPA   Exercise LPA           20    0.726    0.634     0.831 0.00000347 
#>  3 MVPA  LPA   Exercise LPA           30    0.648    0.538     0.782 0.00000548 
#>  4 MVPA  SB    Exercise SB            10    0.822    0.764     0.884 0.000000151
#>  5 MVPA  SB    Exercise SB            20    0.704    0.618     0.802 0.000000136
#>  6 MVPA  SB    Exercise SB            30    0.620    0.519     0.740 0.000000125
#>  7 MVPA  SPT   Exercise SPT           10    0.839    0.778     0.905 0.00000543 
#>  8 MVPA  SPT   Exercise SPT           20    0.734    0.641     0.841 0.00000791 
#>  9 MVPA  SPT   Exercise SPT           30    0.660    0.549     0.795 0.0000114  
#> 10 LPA   MVPA  LPA      Exercise      10    1.30     1.17      1.45  0.00000102 
#> # ℹ 26 more rows
```
