# Print a time-reallocation summary

Prints a text report of a
[`summary.realloc()`](https://marcusvvlopes.github.io/timecoda/reference/summary.realloc.md)
object: a short header on the effect and its reference, then the summary
table as plain tibble output. The report prints automatically at the
console, as base R's `print.summary.lm()` does.

## Usage

``` r
# S3 method for class 'summary_realloc'
print(x, ...)
```

## Arguments

- x:

  A `summary_realloc` object, as returned by
  [`summary.realloc()`](https://marcusvvlopes.github.io/timecoda/reference/summary.realloc.md).

- ...:

  Unused; present for S3 consistency.

## Value

`x`, invisibly.

## See also

[`summary.realloc()`](https://marcusvvlopes.github.io/timecoda/reference/summary.realloc.md)
to build this object;
[`pub_ready()`](https://marcusvvlopes.github.io/timecoda/reference/pub_ready.md)
for the publication `gt` table.

## Examples

``` r
cd <- declare_composition(timeuse, c("MVPA", "LPA", "SB", "SPT"),
                          closure = 1440)
fit <- comp_iv(cd, outcome = "obesity", model = "glm", family = "binomial")
r <- predict_realloc(fit, delta = c(10, 20), from = "SB", to = "MVPA")
print(summary(r))
#> Time-reallocation effects: Odds ratio versus the reference composition
#>   2 reallocations of 10, 20; 95% confidence intervals
#> 
#> # A tibble: 2 × 9
#>   to    from  to_label from_label delta estimate conf.low conf.high     p.value
#>   <chr> <chr> <chr>    <chr>      <dbl>    <dbl>    <dbl>     <dbl>       <dbl>
#> 1 MVPA  SB    MVPA     SB            10    0.822    0.764     0.884 0.000000151
#> 2 MVPA  SB    MVPA     SB            20    0.704    0.618     0.802 0.000000136
```
