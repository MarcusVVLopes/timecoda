# Print a compositional variable summary

Prints a
[`summary.compD()`](https://marcusvvlopes.github.io/timecoda/reference/summary.compD.md)
object as a text report, automatically at the console, as base R's
`print.summary.lm()` does. The report describes the composition by its
parts, closure (with its unit label, when one was declared), number of
observations and how much the raw row totals varied before imputation or
closure. It states which parts had zeroes and, if any were imputed, the
method used, and how many observations were flagged as compositional
outliers when flagging was run. The requested tables follow as plain
tibble output.

## Usage

``` r
# S3 method for class 'summary_compD'
print(x, ...)
```

## Arguments

- x:

  A `summary_compD` object, as returned by
  [`summary.compD()`](https://marcusvvlopes.github.io/timecoda/reference/summary.compD.md).

- ...:

  Unused; present for S3 consistency.

## Value

`x`, invisibly.

## See also

[`summary.compD()`](https://marcusvvlopes.github.io/timecoda/reference/summary.compD.md)
to build this object;
[`pub_ready()`](https://marcusvvlopes.github.io/timecoda/reference/pub_ready.md)
for a formatted `gt` table.

## Examples

``` r
cd <- declare_composition(timeuse, c("MVPA", "LPA", "SB", "SPT"),
                          closure = 1440)
print(summary(cd))
#> <compD> compositional variable summary
#>   Parts: MVPA, LPA, SB, SPT
#>   Closure: 1440
#>   Observations: 300
#> 
#>   Row totals (before closure): mean 1440.0 (SD 0.0), range 1440.0-1440.0
#>   No zeroes detected.
#> 
#> Descriptive statistics
#>   geom_mean: imputed, closed composition; other columns: original data
#> # A tibble: 4 × 10
#>   part  label geom_mean  mean    sd median   p25   p75   min   max
#>   <chr> <chr>     <dbl> <dbl> <dbl>  <dbl> <dbl> <dbl> <dbl> <dbl>
#> 1 MVPA  MVPA       28.1  45.9  44.7   35.0  12.7  61.8   0.1  240 
#> 2 LPA   LPA       235.  243.   86.9  232.  173.  308.   51.8  504.
#> 3 SB    SB        677.  661.  118.   669.  586.  743.  316.  1037.
#> 4 SPT   SPT       499.  491.  106.   489.  426.  559.  240.   815.
#> 
#> Variation matrix
#> # A tibble: 4 × 6
#>   part  label  MVPA   LPA    SB   SPT
#>   <chr> <chr> <dbl> <dbl> <dbl> <dbl>
#> 1 MVPA  MVPA   0    1.67  1.55  1.44 
#> 2 LPA   LPA    1.67 0     0.247 0.260
#> 3 SB    SB     1.55 0.247 0     0.140
#> 4 SPT   SPT    1.44 0.260 0.140 0    
```
