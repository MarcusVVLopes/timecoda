# Validate a compositional dataset

Checks that the composition parts are present in the data, numeric,
strictly positive, free of missing values, and (optionally) that every
row sums to the expected total.

## Usage

``` r
check_composition(data, parts, total = NULL, tol = 0.01, allow_zero = FALSE)
```

## Arguments

- data:

  A data frame containing the composition.

- parts:

  Character vector of column names making up the composition (at least
  two parts).

- total:

  Expected row total for the parts (e.g., `1440` for minutes in a
  24-hour day). If `NULL` (default), row totals are not checked.

- tol:

  Numeric tolerance for the row-total check, as a proportion of `total`
  (default `0.01`, i.e., 1%). The margin allows for the imprecision of
  accelerometer- or diary-derived time-use minutes, which rarely sum to
  exactly 1440.

- allow_zero:

  Logical. If `FALSE` (default), parts must be strictly positive. If
  `TRUE`, zeroes are permitted and only negative values are rejected,
  for data whose zeroes are imputed before any log-ratio transformation,
  as in
  [`declare_composition()`](https://marcusvvlopes.github.io/timecoda/reference/declare_composition.md).

## Value

`data`, invisibly, if all checks pass; otherwise an informative error.

## See also

Other composition helpers:
[`geom_mean_comp()`](https://marcusvvlopes.github.io/timecoda/reference/geom_mean_comp.md),
[`ilr_backtransform()`](https://marcusvvlopes.github.io/timecoda/reference/ilr_backtransform.md),
[`ilr_transform()`](https://marcusvvlopes.github.io/timecoda/reference/ilr_transform.md),
[`variation_matrix()`](https://marcusvvlopes.github.io/timecoda/reference/variation_matrix.md)

## Examples

``` r
check_composition(timeuse, c("MVPA", "LPA", "SB", "SPT"), total = 1440)
```
