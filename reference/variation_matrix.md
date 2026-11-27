# Compositional variation matrix

Computes the variation matrix `T`, with entries
`T[i, j] = var(log(x_i / x_j))`. Small values indicate parts that are
highly proportional.

## Usage

``` r
variation_matrix(data, parts, robust = FALSE, ...)
```

## Arguments

- data:

  A data frame containing the composition.

- parts:

  Character vector of column names making up the composition (at least
  two parts).

- robust:

  Passed to
  [`compositions::variation.acomp()`](https://rdrr.io/pkg/compositions/man/variation.html).
  Defaults to `FALSE` (the classical variation matrix); pass `TRUE` for
  the MCD-based robust variant.

- ...:

  Passed on to
  [`compositions::variation.acomp()`](https://rdrr.io/pkg/compositions/man/variation.html).

## Value

A symmetric numeric matrix with `parts` as row and column names.

## References

Aitchison, J. (1986). *The Statistical Analysis of Compositional Data.*
Chapman & Hall.

## See also

Other composition helpers:
[`check_composition()`](https://marcusvvlopes.github.io/timecoda/reference/check_composition.md),
[`geom_mean_comp()`](https://marcusvvlopes.github.io/timecoda/reference/geom_mean_comp.md),
[`ilr_backtransform()`](https://marcusvvlopes.github.io/timecoda/reference/ilr_backtransform.md),
[`ilr_transform()`](https://marcusvvlopes.github.io/timecoda/reference/ilr_transform.md)

## Examples

``` r
variation_matrix(timeuse, c("MVPA", "LPA", "SB", "SPT"))
#>          MVPA       LPA        SB       SPT
#> MVPA 0.000000 1.6691328 1.5485439 1.4432409
#> LPA  1.669133 0.0000000 0.2466768 0.2601569
#> SB   1.548544 0.2466768 0.0000000 0.1403835
#> SPT  1.443241 0.2601569 0.1403835 0.0000000
```
