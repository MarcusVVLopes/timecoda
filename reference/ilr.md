# Composition marker for `comp_iv()` formulas

A marker, valid only inside a `formula` passed to
[`comp_iv()`](https://marcusvvlopes.github.io/timecoda/reference/comp_iv.md),
that stands for the whole composition: the pivot isometric log-ratio
(ILR) coordinates `ilr1, ..., ilr{D-1}`. The formula parser expands it
before the model is fitted, and calling `ilr()` directly is an error.

## Usage

``` r
ilr()
```

## Value

None; `ilr()` signals an error when called outside a formula.

## Details

A
[`comp_iv()`](https://marcusvvlopes.github.io/timecoda/reference/comp_iv.md)
formula is fitted as written once `ilr()` is expanded, by one of two
rules that depend on where the marker sits. No composition term is
refused; only `ilr()` with arguments is an error.

Inside the formula operators `+ - * : / ^ %in% | ||` and parentheses,
`ilr()` is replaced by the sum `(ilr1 + ... + ilr{D-1})`, and R's
formula algebra applies to it. `outcome ~ ilr() * sex` fits the
composition, `sex` and every composition-by-`sex` interaction;
`(ilr() | school)` is an `lme4` random slope on every coordinate.
`ilr()^2` means what `(ilr1 + ilr2 + ilr3)^2` means in R: the main
effects and every pairwise product of coordinates, with no squares.

Inside any other call, `ilr()` is distributed: the call is repeated once
per coordinate and the copies are added. `I(ilr()^2)` becomes
`(I(ilr1^2) + ... + I(ilr{D-1}^2))` and `poly(ilr(), 2)` becomes
`(poly(ilr1, 2) + ... + poly(ilr{D-1}, 2))`;
[`splines::ns()`](https://rdrr.io/r/splines/ns.html),
[`log()`](https://rdrr.io/r/base/Log.html) and
[`scale()`](https://rdrr.io/r/base/scale.html) are treated alike. Every
`ilr()` within one such call takes the same coordinate, so
`I(ilr() * ilr())` gives the squares; the pairwise products come from
`ilr()^2`.

The full quadratic `ilr()^2 + I(ilr()^2)` (main effects, squares and
pairwise products) is the same model in every ILR basis, so its fit and
its reallocations do not depend on the order of the parts. The fit of a
squares-only or `poly(ilr(), 2)` specification depends on the part order
given to
[`declare_composition()`](https://marcusvvlopes.github.io/timecoda/reference/declare_composition.md).
On the `timeuse` data, the odds ratio of obesity for 10 minutes moved
from SB to MVPA under `ilr() + I(ilr()^2)` ranged from 0.828 to 0.834
across three part orders.
[`comp_iv()`](https://marcusvvlopes.github.io/timecoda/reference/comp_iv.md)
fits the full quadratic without a formula through `curvature = TRUE`.

[`comp_iv()`](https://marcusvvlopes.github.io/timecoda/reference/comp_iv.md)
warns when the formula leaves out the composition's main effect, that
is, when a coordinate appears neither as a term of its own nor as the
first argument of [`poly()`](https://rdrr.io/r/stats/poly.html), `ns()`
or `bs()`.

## See also

[`comp_iv()`](https://marcusvvlopes.github.io/timecoda/reference/comp_iv.md),
whose `formula` argument accepts this marker.

## Examples

``` r
cd <- declare_composition(timeuse, c("MVPA", "LPA", "SB", "SPT"),
                          closure = 1440)
# The composition enters via ilr(); sex and its interaction are
# ordinary terms
comp_iv(cd, obesity ~ ilr() * sex, model = "glm", family = "binomial")
#> <timecoda_fit> compositional generalized linear model (binomial)
#>   Parts: MVPA, LPA, SB, SPT (total = 1440)
#>   Reference composition: MVPA = 28.1, LPA = 235.5, SB = 677.2, SPT = 499.2
#> 
#> 
#> Call:  stats::glm(formula = obesity ~ (ilr1 + ilr2 + ilr3) * sex, family = "binomial", 
#>     data = data)
#> 
#> Coefficients:
#>  (Intercept)          ilr1          ilr2          ilr3       sexmale  
#>     -3.68401      -0.69592      -0.87505       1.52183       1.15199  
#> ilr1:sexmale  ilr2:sexmale  ilr3:sexmale  
#>      0.02295       1.14463      -1.06266  
#> 
#> Degrees of Freedom: 299 Total (i.e. Null);  292 Residual
#> Null Deviance:       355.8 
#> Residual Deviance: 318.8     AIC: 334.8
# A full quadratic in the composition, which does not depend on the
# part order
comp_iv(cd, zbmi ~ ilr()^2 + I(ilr()^2) + sex, model = "lm")
#> <timecoda_fit> compositional linear model
#>   Parts: MVPA, LPA, SB, SPT (total = 1440)
#>   Reference composition: MVPA = 28.1, LPA = 235.5, SB = 677.2, SPT = 499.2
#> 
#> 
#> Call:
#> stats::lm(formula = zbmi ~ (ilr1 + ilr2 + ilr3)^2 + (I(ilr1^2) + 
#>     I(ilr2^2) + I(ilr3^2)) + sex, data = data)
#> 
#> Coefficients:
#> (Intercept)         ilr1         ilr2         ilr3    I(ilr1^2)    I(ilr2^2)  
#>    -2.42867     -0.98660      0.04431      1.59720     -0.04117      0.53721  
#>   I(ilr3^2)      sexmale    ilr1:ilr2    ilr1:ilr3    ilr2:ilr3  
#>    -1.06062     -0.07917     -0.15457     -0.22287      0.18520  
#> 
```
