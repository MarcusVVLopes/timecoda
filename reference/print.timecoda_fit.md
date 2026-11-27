# Print a timecoda model fit

Prints a compact summary of a `timecoda_fit` object, followed by the
underlying fitted model object. The summary gives the model class and
(for GLM-family engines) its family, the composition parts and total
(with its unit label, when one was declared) and the reference
composition (the sample's closed geometric mean, from which
reallocations start).

## Usage

``` r
# S3 method for class 'timecoda_fit'
print(x, ...)
```

## Arguments

- x:

  A `timecoda_fit` object (or subclass, e.g., `comp_iv`), as returned by
  [`comp_iv()`](https://marcusvvlopes.github.io/timecoda/reference/comp_iv.md).

- ...:

  Further arguments (unused; present for S3 consistency).

## Value

`x`, invisibly.

## See also

[`comp_iv()`](https://marcusvvlopes.github.io/timecoda/reference/comp_iv.md)
to fit the model;
[`summary.comp_iv()`](https://marcusvvlopes.github.io/timecoda/reference/summary.comp_iv.md)
for the full effect table.

## Examples

``` r
cd <- declare_composition(timeuse, c("MVPA", "LPA", "SB", "SPT"),
                          closure = 1440)
m <- comp_iv(cd, outcome = "obesity", model = "glm", family = "binomial")
print(m)
#> <timecoda_fit> compositional generalized linear model (binomial)
#>   Parts: MVPA, LPA, SB, SPT (total = 1440)
#>   Reference composition: MVPA = 28.1, LPA = 235.5, SB = 677.2, SPT = 499.2
#> 
#> 
#> Call:  stats::glm(formula = obesity ~ ilr1 + ilr2 + ilr3, family = "binomial", 
#>     data = data)
#> 
#> Coefficients:
#> (Intercept)         ilr1         ilr2         ilr3  
#>     -3.0796      -0.6916      -0.2635       0.9014  
#> 
#> Degrees of Freedom: 299 Total (i.e. Null);  296 Residual
#> Null Deviance:       355.8 
#> Residual Deviance: 322   AIC: 330
```
