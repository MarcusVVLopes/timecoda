# Describe a compositional variable

Describes a
[`declare_composition()`](https://marcusvvlopes.github.io/timecoda/reference/declare_composition.md)
object by a per-part descriptive table and the compositional
[`variation_matrix()`](https://marcusvvlopes.github.io/timecoda/reference/variation_matrix.md).
The table leads with the closed geometric mean of the imputed
composition (the compositional centre, from
[`geom_mean_comp()`](https://marcusvvlopes.github.io/timecoda/reference/geom_mean_comp.md));
its mean, SD, quartiles and range describe the original data, before
zero imputation and closure. Printing the result shows a text report, as
with [`summary.lm()`](https://rdrr.io/r/stats/summary.lm.html), and
[`pub_ready()`](https://marcusvvlopes.github.io/timecoda/reference/pub_ready.md)
builds a publication `gt` table from it.

## Usage

``` r
# S3 method for class 'compD'
summary(object, desc = TRUE, var_matrix = TRUE, ...)
```

## Arguments

- object:

  A `compD` object from
  [`declare_composition()`](https://marcusvvlopes.github.io/timecoda/reference/declare_composition.md).

- desc:

  Logical; include the per-part descriptive table. Default `TRUE`.

- var_matrix:

  Logical; include the variation matrix. Default `TRUE`.

- ...:

  Unused; present for S3 consistency.

## Value

An object of class `summary_compD`, printed as a text report by
[`print.summary_compD()`](https://marcusvvlopes.github.io/timecoda/reference/print.summary_compD.md)
and passed to
[`pub_ready()`](https://marcusvvlopes.github.io/timecoda/reference/pub_ready.md)
for a `gt` table. A list with elements:

- `desc`:

  One row per part: `part`, `label` (display label), then `geom_mean`
  (the closed geometric mean of the imputed composition), and `mean`,
  `sd`, `median`, `p25`, `p75`, `min`, `max` of the original data,
  before zero imputation and closure. `NULL` when `desc = FALSE`.

- `var_matrix`:

  The variation matrix of the imputed composition, with leading `part`
  and `label` columns. `NULL` when `var_matrix = FALSE`.

- `parts`, `labels`, `closure`, `units`:

  The part names, resolved display labels and closure of the `compD`,
  carried through for
  [`pub_ready()`](https://marcusvvlopes.github.io/timecoda/reference/pub_ready.md),
  and its unit label (`NULL` when none was declared), shown by
  [`print.summary_compD()`](https://marcusvvlopes.github.io/timecoda/reference/print.summary_compD.md).

- `row_totals`:

  A one-row tibble (`mean`, `sd`, `min`, `max`) of the raw part sums per
  observation, before imputation or closure.
  [`print.summary_compD()`](https://marcusvvlopes.github.io/timecoda/reference/print.summary_compD.md)
  reports it to show how much row totals varied in the original data.

- `n_obs`:

  The number of observations in the composition.

- `imput_table`, `impute_method`:

  Carried through from the `compD` object (see
  [`declare_composition()`](https://marcusvvlopes.github.io/timecoda/reference/declare_composition.md)),
  for
  [`print.summary_compD()`](https://marcusvvlopes.github.io/timecoda/reference/print.summary_compD.md)'s
  zero-detection report.

- `outliers`:

  A one-row tibble (`n_outlier`, `pct`, `limit`) summarizing the
  `.outlier` columns
  [`declare_composition()`](https://marcusvvlopes.github.io/timecoda/reference/declare_composition.md)
  appended to `data_raw`, for
  [`print.summary_compD()`](https://marcusvvlopes.github.io/timecoda/reference/print.summary_compD.md)'s
  outlier line. `NULL` when outliers were not flagged.

## See also

[`declare_composition()`](https://marcusvvlopes.github.io/timecoda/reference/declare_composition.md)
to build the `compD` object;
[`plot.compD()`](https://marcusvvlopes.github.io/timecoda/reference/plot.compD.md)
for the figures;
[`pub_ready()`](https://marcusvvlopes.github.io/timecoda/reference/pub_ready.md)
for a formatted `gt` table built from this summary.

Other declare and describe composition:
[`comp_impute()`](https://marcusvvlopes.github.io/timecoda/reference/comp_impute.md),
[`comp_outliers()`](https://marcusvvlopes.github.io/timecoda/reference/comp_outliers.md),
[`declare_composition()`](https://marcusvvlopes.github.io/timecoda/reference/declare_composition.md),
[`plot.compD()`](https://marcusvvlopes.github.io/timecoda/reference/plot.compD.md)

## Examples

``` r
cd <- declare_composition(timeuse, c("MVPA", "LPA", "SB", "SPT"),
                          closure = 1440)
s <- summary(cd)
s$desc
#> # A tibble: 4 × 10
#>   part  label geom_mean  mean    sd median   p25   p75   min   max
#>   <chr> <chr>     <dbl> <dbl> <dbl>  <dbl> <dbl> <dbl> <dbl> <dbl>
#> 1 MVPA  MVPA       28.1  45.9  44.7   35.0  12.7  61.8   0.1  240 
#> 2 LPA   LPA       235.  243.   86.9  232.  173.  308.   51.8  504.
#> 3 SB    SB        677.  661.  118.   669.  586.  743.  316.  1037.
#> 4 SPT   SPT       499.  491.  106.   489.  426.  559.  240.   815.
```
