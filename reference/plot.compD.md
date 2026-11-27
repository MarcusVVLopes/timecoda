# Visualize a compositional variable

Draws per-part density curves and boxplots of the untransformed parts of
a
[`declare_composition()`](https://marcusvvlopes.github.io/timecoda/reference/declare_composition.md)
object, the ternary plots of the composition and the heatmap of its
variation matrix, all from the object's imputed composition. When the
object carries flagged compositional outliers
([`declare_composition()`](https://marcusvvlopes.github.io/timecoda/reference/declare_composition.md)'s
`detect_outliers`), two further panels show the outlier distances and
the flagged observations in the simplex.

## Usage

``` r
# S3 method for class 'compD'
plot(
  x,
  histogram = TRUE,
  boxplot = TRUE,
  ternary = TRUE,
  variation = TRUE,
  outliers = TRUE,
  ...
)
```

## Arguments

- x:

  A `compD` object from
  [`declare_composition()`](https://marcusvvlopes.github.io/timecoda/reference/declare_composition.md).

- histogram:

  Logical; include per-part density curves
  ([`ggplot2::geom_density()`](https://ggplot2.tidyverse.org/reference/geom_density.html),
  faceted by part). Default `TRUE`.

- boxplot:

  Logical; include per-part boxplots
  ([`ggplot2::geom_boxplot()`](https://ggplot2.tidyverse.org/reference/geom_boxplot.html)).
  Default `TRUE`.

- ternary:

  Logical; include the ternary plot(s)
  ([`plot_ternary()`](https://marcusvvlopes.github.io/timecoda/reference/plot_ternary.md)).
  Default `TRUE`. Set `FALSE` for a two-part composition, which has no
  ternary representation; this also drops the outlier ternary.

- variation:

  Logical; include the variation-matrix heatmap
  ([`plot_variation()`](https://marcusvvlopes.github.io/timecoda/reference/plot_variation.md)).
  Default `TRUE`.

- outliers:

  Logical; include the outlier panels. Default `TRUE`. The panels appear
  only when `x` carries flagged outliers, so a composition declared
  without `detect_outliers` is unaffected.

- ...:

  Further arguments passed to
  [`plot_ternary()`](https://marcusvvlopes.github.io/timecoda/reference/plot_ternary.md)
  (e.g., `group`). They do not reach the outlier panels, which always
  describe the whole sample (see Value).

## Value

Invisibly, a named list of the requested plots, with elements
`histogram`, `boxplot`, `ternary` (a single `ggtern` plot or, for more
than three parts, a named list of `choose(D, 3)` panels) and `variation`
(the variation-matrix heatmap). When `x` carries flagged outliers, two
further elements follow: `outliers`, the Mahalanobis distance of each
observation against the cutoff, and, unless `ternary = FALSE`,
`outlier_ternary`, the ternary panel(s) with the flagged observations
highlighted. Both describe the whole sample and are never split by
`group`. Calling
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) draws the plots
as a side effect.

## See also

[`declare_composition()`](https://marcusvvlopes.github.io/timecoda/reference/declare_composition.md)
to build the `compD` object;
[`summary.compD()`](https://marcusvvlopes.github.io/timecoda/reference/summary.compD.md)
for the descriptive tables;
[`plot_ternary()`](https://marcusvvlopes.github.io/timecoda/reference/plot_ternary.md)
and
[`plot_variation()`](https://marcusvvlopes.github.io/timecoda/reference/plot_variation.md)
for the ternary and variation-matrix figures alone.

Other declare and describe composition:
[`comp_impute()`](https://marcusvvlopes.github.io/timecoda/reference/comp_impute.md),
[`comp_outliers()`](https://marcusvvlopes.github.io/timecoda/reference/comp_outliers.md),
[`declare_composition()`](https://marcusvvlopes.github.io/timecoda/reference/declare_composition.md),
[`summary.compD()`](https://marcusvvlopes.github.io/timecoda/reference/summary.compD.md)

## Examples

``` r
# \donttest{
cd <- declare_composition(timeuse, c("MVPA", "LPA", "SB", "SPT"),
                          closure = 1440)
plots <- plot(cd)







plots$boxplot

# }
```
