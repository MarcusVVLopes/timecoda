# Geometric mean composition

Computes the closed geometric mean of a composition on the simplex. The
sample's closed geometric mean, stored in a
[`comp_iv()`](https://marcusvvlopes.github.io/timecoda/reference/comp_iv.md)
fit, is the default reference composition of
[`predict_realloc()`](https://marcusvvlopes.github.io/timecoda/reference/predict_realloc.md).

## Usage

``` r
geom_mean_comp(data, parts, total = 1, robust = FALSE, ...)
```

## Arguments

- data:

  A data frame containing the composition.

- parts:

  Character vector of column names making up the composition (at least
  two parts).

- total:

  Row total to close the result to (e.g., `1440` for minutes per day).
  Defaults to `1`, returning proportions.

- robust:

  Passed to
  [`compositions::mean.acomp()`](https://rdrr.io/pkg/compositions/man/meanAcomp.html).
  Defaults to `FALSE` (the classical geometric mean); pass `TRUE` for
  the MCD-based robust centre.

- ...:

  Passed on to
  [`compositions::mean.acomp()`](https://rdrr.io/pkg/compositions/man/meanAcomp.html).

## Value

A one-row tibble with one column per part, summing to `total`.

## References

Aitchison, J. (1986). *The Statistical Analysis of Compositional Data.*
Chapman & Hall.

## See also

Other composition helpers:
[`check_composition()`](https://marcusvvlopes.github.io/timecoda/reference/check_composition.md),
[`ilr_backtransform()`](https://marcusvvlopes.github.io/timecoda/reference/ilr_backtransform.md),
[`ilr_transform()`](https://marcusvvlopes.github.io/timecoda/reference/ilr_transform.md),
[`variation_matrix()`](https://marcusvvlopes.github.io/timecoda/reference/variation_matrix.md)

## Examples

``` r
geom_mean_comp(timeuse, c("MVPA", "LPA", "SB", "SPT"), total = 1440)
#> # A tibble: 1 × 4
#>    MVPA   LPA    SB   SPT
#>   <dbl> <dbl> <dbl> <dbl>
#> 1  28.1  235.  677.  499.
```
