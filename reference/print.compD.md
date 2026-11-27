# Print a compositional variable

Prints the parts, the closure (with its unit label, when one was
declared) and the number of observations of a
[`declare_composition()`](https://marcusvvlopes.github.io/timecoda/reference/declare_composition.md)
object, with a line each for imputed zeroes and flagged outliers when
there are any.
[`summary.compD()`](https://marcusvvlopes.github.io/timecoda/reference/summary.compD.md)
gives the descriptive tables and
[`plot.compD()`](https://marcusvvlopes.github.io/timecoda/reference/plot.compD.md)
the figures.

## Usage

``` r
# S3 method for class 'compD'
print(x, ...)
```

## Arguments

- x:

  A `compD` object from
  [`declare_composition()`](https://marcusvvlopes.github.io/timecoda/reference/declare_composition.md).

- ...:

  Unused; present for S3 consistency.

## Value

`x`, invisibly.

## See also

[`summary.compD()`](https://marcusvvlopes.github.io/timecoda/reference/summary.compD.md)
for the descriptive tables and the fuller provenance report;
[`plot.compD()`](https://marcusvvlopes.github.io/timecoda/reference/plot.compD.md)
for the figures.

## Examples

``` r
declare_composition(timeuse, c("MVPA", "LPA", "SB", "SPT"), closure = 1440)
#> <compD> compositional variable
#>   Parts: MVPA, LPA, SB, SPT
#>   Closure: 1440
#>   Observations: 300
```
