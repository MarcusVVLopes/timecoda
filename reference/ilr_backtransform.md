# Back-transform ILR coordinates to a composition

Inverts
[`ilr_transform()`](https://marcusvvlopes.github.io/timecoda/reference/ilr_transform.md):
converts `D - 1` pivot ILR coordinates back to the `D`-part composition,
closed to `total`.

## Usage

``` r
ilr_backtransform(coords, parts, total = 1, pivot = 1)
```

## Arguments

- coords:

  A data frame or matrix of ILR coordinates with `D - 1` columns, in the
  column order produced by
  [`ilr_transform()`](https://marcusvvlopes.github.io/timecoda/reference/ilr_transform.md).

- parts:

  Character vector of the `D` composition part names, in the same order
  passed to
  [`ilr_transform()`](https://marcusvvlopes.github.io/timecoda/reference/ilr_transform.md).

- total:

  Row total the composition is closed to (e.g., `1440`). Defaults to
  `1`, returning proportions.

- pivot:

  The pivot used in
  [`ilr_transform()`](https://marcusvvlopes.github.io/timecoda/reference/ilr_transform.md)
  (name or position).

## Value

A tibble with one column per part, rows summing to `total`.

## References

Filzmoser, P., Hron, K., & Templ, M. (2018). *Applied Compositional Data
Analysis: With Worked Examples in R.* Springer.
[doi:10.1007/978-3-319-96422-5](https://doi.org/10.1007/978-3-319-96422-5)

## See also

Other composition helpers:
[`check_composition()`](https://marcusvvlopes.github.io/timecoda/reference/check_composition.md),
[`geom_mean_comp()`](https://marcusvvlopes.github.io/timecoda/reference/geom_mean_comp.md),
[`ilr_transform()`](https://marcusvvlopes.github.io/timecoda/reference/ilr_transform.md),
[`variation_matrix()`](https://marcusvvlopes.github.io/timecoda/reference/variation_matrix.md)

## Examples

``` r
ilr <- ilr_transform(timeuse, c("MVPA", "LPA", "SB", "SPT"), total = 1440)
comp <- ilr_backtransform(
  ilr[c("ilr1", "ilr2", "ilr3")],
  parts = c("MVPA", "LPA", "SB", "SPT"),
  total = 1440
)
all.equal(comp$MVPA, timeuse$MVPA)
#> [1] TRUE
```
