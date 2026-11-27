# Impute compositional zeroes

Replaces zero values in a set of composition parts using one of the
continuous imputation methods of `zCompositions`.
[`declare_composition()`](https://marcusvvlopes.github.io/timecoda/reference/declare_composition.md)
calls this function on the raw (pre-closure) scale, before it closes the
composition.

## Usage

``` r
comp_impute(
  data,
  parts,
  method = c("lrEM", "multRepl", "lrDA", "multLN", "lrSVD"),
  ...
)
```

## Arguments

- data:

  A data frame containing the composition parts.

- parts:

  Character vector of column names making up the composition.

- method:

  Imputation method: `"lrEM"` (default, log-ratio EM,
  [`zCompositions::lrEM()`](https://rdrr.io/pkg/zCompositions/man/lrEM.html)),
  `"multRepl"` (simple multiplicative replacement,
  [`zCompositions::multRepl()`](https://rdrr.io/pkg/zCompositions/man/multRepl.html)),
  `"lrDA"` (log-ratio data augmentation,
  [`zCompositions::lrDA()`](https://rdrr.io/pkg/zCompositions/man/lrDA.html)),
  `"multLN"` (multiplicative lognormal,
  [`zCompositions::multLN()`](https://rdrr.io/pkg/zCompositions/man/multLN.html))
  or `"lrSVD"` (SVD-based,
  [`zCompositions::lrSVD()`](https://rdrr.io/pkg/zCompositions/man/lrSVD.html)).
  `"multLN"` is called with `rob = TRUE` (the ROS-based alternative)
  unless `rob` is supplied through `...`; its classical-MLE default
  (`rob = FALSE`) fails inside `comp_impute()` with "could not find
  function \\survreg\\".
  [`zCompositions::multKM()`](https://rdrr.io/pkg/zCompositions/man/multKM.html)
  (Kaplan-Meier based) fails in the same way and is not offered.

- ...:

  Additional arguments forwarded to the chosen method (most usefully
  `dl`, to override the derived per-part detection limit, or a
  method-specific control such as `lrEM`'s
  `rob`/`tolerance`/`max.iter`). The internal defaults (`label = 0`, the
  derived `dl`, and `suppress.print = TRUE` where the method accepts it)
  fill in only the arguments not supplied.

## Value

A list with elements `parts_imputed` (a tibble: the zero-free parts on
the original scale), `imput_table` (a tibble summarizing zero occurrence
per part: `part`, `n`, `n_zero`, `pct_zero`, `min_nonzero`), and
`impute_method` (the `method` used, or `NA_character_` when no part had
a zero). When no part has a zero, `parts_imputed` equals the input
`parts` columns unchanged and no imputation method is called.

## References

Palarea-Albaladejo, J., & Martín-Fernández, J.A. (2015). zCompositions –
R package for multivariate imputation of left-censored data under a
compositional approach. *Chemometrics and Intelligent Laboratory
Systems*, 143, 85-96.
[doi:10.1016/j.chemolab.2015.02.019](https://doi.org/10.1016/j.chemolab.2015.02.019)

## See also

[`declare_composition()`](https://marcusvvlopes.github.io/timecoda/reference/declare_composition.md),
which calls this internally.

Other declare and describe composition:
[`comp_outliers()`](https://marcusvvlopes.github.io/timecoda/reference/comp_outliers.md),
[`declare_composition()`](https://marcusvvlopes.github.io/timecoda/reference/declare_composition.md),
[`plot.compD()`](https://marcusvvlopes.github.io/timecoda/reference/plot.compD.md),
[`summary.compD()`](https://marcusvvlopes.github.io/timecoda/reference/summary.compD.md)

## Examples

``` r
cd_data <- timeuse
cd_data$MVPA[1:3] <- 0
imp <- comp_impute(cd_data, c("MVPA", "LPA", "SB", "SPT"))
imp$imput_table
#> # A tibble: 4 × 5
#>   part      n n_zero pct_zero min_nonzero
#>   <chr> <int>  <int>    <dbl>       <dbl>
#> 1 MVPA    300      3        1         0.1
#> 2 LPA     300      0        0        51.8
#> 3 SB      300      0        0       316. 
#> 4 SPT     300      0        0       240. 
```
