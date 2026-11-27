# Print a composition-as-predictor summary

Prints a
[`summary.comp_iv()`](https://marcusvvlopes.github.io/timecoda/reference/summary.comp_iv.md)
object as a text report, automatically at the console, as base R's
`print.summary.lm()` does. The report gives the model description, the
composition's parts and total (with its unit label, when one was
declared) and the omnibus composition test. When the fit carries them,
it adds the moderation test, the curvature test, the total-as-covariate
test, each grouping term's intraclass correlation or median hazard ratio
and any captured backend convergence messages. The effect table follows
as plain tibble output, with a closing note when the per-part effects
are slopes at the reference composition.

## Usage

``` r
# S3 method for class 'summary_comp_iv'
print(x, ...)
```

## Arguments

- x:

  A `summary_comp_iv` object, as returned by
  [`summary.comp_iv()`](https://marcusvvlopes.github.io/timecoda/reference/summary.comp_iv.md).

- ...:

  Unused; present for S3 consistency.

## Value

`x`, invisibly.

## See also

[`summary.comp_iv()`](https://marcusvvlopes.github.io/timecoda/reference/summary.comp_iv.md)
to build this object;
[`pub_ready()`](https://marcusvvlopes.github.io/timecoda/reference/pub_ready.md)
for a formatted `gt` table.

## Examples

``` r
cd <- declare_composition(timeuse, c("MVPA", "LPA", "SB", "SPT"),
                          closure = 1440)
m <- comp_iv(cd, outcome = "obesity", model = "glm", family = "binomial")
print(summary(m))
#> <comp_iv> compositional generalized linear model (binomial)
#>   Parts: MVPA, LPA, SB, SPT (total = 1440)
#>   Omnibus composition test (Wald): chi-squared(3) = 28.922, p = 2.33e-06
#>   Convergence: OK
#> 
#> Compositional effects (Odds ratio)
#> # A tibble: 5 × 6
#>   parameter          label              estimate conf.low conf.high     p.value
#>   <chr>              <chr>                 <dbl>    <dbl>     <dbl>       <dbl>
#> 1 Global Composition Global Composition   NA       NA        NA     0.00000233 
#> 2 MVPA vs remaining  MVPA                  0.501    0.382     0.656 0.000000489
#> 3 LPA vs remaining   LPA                   0.982    0.506     1.91  0.958      
#> 4 SB vs remaining    SB                    2.98     1.16      7.65  0.0234     
#> 5 SPT vs remaining   SPT                   0.683    0.273     1.71  0.415      
```
