# Detect compositional outliers

Flags observations whose *ratios between parts* are atypical, using the
robust Mahalanobis distance of Filzmoser and Hron (2008) computed in
pivot isometric log-ratio coordinates. An observation is an outlier when
its distance exceeds `sqrt(qchisq(quantile, df = D - 1))` for a `D`-part
composition. The distance is invariant to the closure of the composition
and to the choice of ILR pivot.

## Usage

``` r
comp_outliers(data, parts, quantile = 0.975, ...)
```

## Arguments

- data:

  A data frame containing the composition parts, free of zeroes.

- parts:

  Character vector of column names making up the composition (at least
  two parts).

- quantile:

  Chi-squared probability setting the outlier cutoff. Defaults to
  `0.975`.

- ...:

  Additional arguments forwarded to
  [`robCompositions::outCoDa()`](https://rdrr.io/pkg/robCompositions/man/outCoDa.html).
  The useful one is `alpha`, the proportion of observations MCD retains
  in its subset search (default `0.5`, the most robust; higher values
  trade robustness for efficiency). The `coda` argument of `outCoDa()`
  cannot be overridden, and passing it is an error.

## Value

A tibble with one row per observation and the columns `.row` (position
in `data`), `distance` (robust Mahalanobis distance in ILR coordinates),
`limit` (the cutoff, constant down the column), and `is_outlier`
(logical, `distance > limit`). Rows are never removed; acting on the
flags is left to the caller.

## Details

The centre and covariance come from the minimum covariance determinant
(MCD) estimator, which resists the outliers it is meant to find. Its
classical counterpart, the arithmetic mean and sample covariance, is
itself distorted by them (the masking effect) and is not offered.

MCD estimates the covariance from a randomized subset search, so
distances vary slightly between calls. A call to
[`set.seed()`](https://rdrr.io/r/base/Random.html) beforehand makes the
flags reproducible.

Detection requires a zero-free composition. Zeroes can be imputed first
with
[`comp_impute()`](https://marcusvvlopes.github.io/timecoda/reference/comp_impute.md),
or the composition built with
[`declare_composition()`](https://marcusvvlopes.github.io/timecoda/reference/declare_composition.md),
which imputes before it flags.

## References

Filzmoser, P., & Hron, K. (2008). Outlier detection for compositional
data using robust methods. *Mathematical Geosciences*, 40, 233-248.
[doi:10.1007/s11004-007-9141-5](https://doi.org/10.1007/s11004-007-9141-5)

Rousseeuw, P.J., & Van Driessen, K. (1999). A fast algorithm for the
minimum covariance determinant estimator. *Technometrics*, 41, 212-223.
[doi:10.1080/00401706.1999.10485670](https://doi.org/10.1080/00401706.1999.10485670)

## See also

[`declare_composition()`](https://marcusvvlopes.github.io/timecoda/reference/declare_composition.md),
which calls this when `detect_outliers = TRUE`;
[`comp_impute()`](https://marcusvvlopes.github.io/timecoda/reference/comp_impute.md)
to remove zeroes first.

Other declare and describe composition:
[`comp_impute()`](https://marcusvvlopes.github.io/timecoda/reference/comp_impute.md),
[`declare_composition()`](https://marcusvvlopes.github.io/timecoda/reference/declare_composition.md),
[`plot.compD()`](https://marcusvvlopes.github.io/timecoda/reference/plot.compD.md),
[`summary.compD()`](https://marcusvvlopes.github.io/timecoda/reference/summary.compD.md)

## Examples

``` r
set.seed(1)
out <- comp_outliers(timeuse, c("MVPA", "LPA", "SB", "SPT"))
sum(out$is_outlier)
#> [1] 14
```
