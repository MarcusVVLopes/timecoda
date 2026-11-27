# Isometric log-ratio (pivot) coordinates of a composition

Transforms a `D`-part composition into `D - 1` pivot ILR coordinates via
[`robCompositions::pivotCoord()`](https://rdrr.io/pkg/robCompositions/man/pivotCoord.html),
appended to `data` as plain numeric columns `ilr1 ... ilr{D-1}`. The
first coordinate expresses the pivot part relative to the geometric mean
of the remaining parts. Reallocation predictions are invariant to the
choice of pivot when the model is linear in the coordinates or holds
their full quadratic
([`ilr()`](https://marcusvvlopes.github.io/timecoda/reference/ilr.md)).

## Usage

``` r
ilr_transform(data, parts, total = NULL, pivot = 1)
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

- pivot:

  The part forming the numerator of the first coordinate: either a
  column name from `parts` or its position. Defaults to the first part.

## Value

`data` as a tibble with `length(parts) - 1` additional numeric columns
named `ilr1`, `ilr2`, ....

## References

Filzmoser, P., Hron, K., & Templ, M. (2018). *Applied Compositional Data
Analysis: With Worked Examples in R.* Springer.
[doi:10.1007/978-3-319-96422-5](https://doi.org/10.1007/978-3-319-96422-5)

## See also

Other composition helpers:
[`check_composition()`](https://marcusvvlopes.github.io/timecoda/reference/check_composition.md),
[`geom_mean_comp()`](https://marcusvvlopes.github.io/timecoda/reference/geom_mean_comp.md),
[`ilr_backtransform()`](https://marcusvvlopes.github.io/timecoda/reference/ilr_backtransform.md),
[`variation_matrix()`](https://marcusvvlopes.github.io/timecoda/reference/variation_matrix.md)

## Examples

``` r
ilr_transform(timeuse, c("MVPA", "LPA", "SB", "SPT"), total = 1440)
#> # A tibble: 300 × 16
#>    id    sex    school region  MVPA   LPA    SB   SPT   zbmi obesity  time event
#>    <chr> <fct>  <fct>  <fct>  <dbl> <dbl> <dbl> <dbl>  <dbl> <fct>   <dbl> <int>
#>  1 p001  male   schoo… A       36.1  320   541.  543. -1.31  yes      3.84     1
#>  2 p002  male   schoo… C       14.7  214.  815.  397  -0.472 no      16.2      1
#>  3 p003  male   schoo… C      195.   206.  500.  539. -1.45  no      22.0      1
#>  4 p004  female schoo… A       95.5  229.  445.  670. -2.40  no      19.6      0
#>  5 p005  male   schoo… C        1.7  504.  479.  455.  2.56  no       0.21     1
#>  6 p006  male   schoo… C       11.2  150.  583.  696. -0.259 no      21.4      1
#>  7 p007  male   schoo… B        2    186.  861   391.  2.35  no       1.93     1
#>  8 p008  female schoo… A       68.7  397.  646.  329. -2.11  no      27.9      1
#>  9 p009  female schoo… C       56.6  307   625.  451. -0.709 no      20.2      1
#> 10 p010  female schoo… B       61.1  224.  788.  367.  0.413 no      10.6      1
#> # ℹ 290 more rows
#> # ℹ 4 more variables: grade <ord>, ilr1 <dbl>, ilr2 <dbl>, ilr3 <dbl>
```
