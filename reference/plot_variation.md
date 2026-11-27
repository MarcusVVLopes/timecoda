# Heatmap of the compositional variation matrix

Draws the compositional variation matrix as a tile heatmap, shading each
cell by the pairwise log-ratio variance `var(log(x_i / x_j))`, zero by
construction on the diagonal. Small values mark highly proportional
parts (parts that move together); with the default colours they are
drawn darker and large values lighter. The axes carry the part labels of
the
[`declare_composition()`](https://marcusvvlopes.github.io/timecoda/reference/declare_composition.md)
object.

## Usage

``` r
plot_variation(
  x,
  digits = 2,
  show_values = TRUE,
  low = "#132B43",
  high = "#EAF3FB",
  theme = NULL
)
```

## Arguments

- x:

  A `compD` object from
  [`declare_composition()`](https://marcusvvlopes.github.io/timecoda/reference/declare_composition.md).

- digits:

  Decimal places for the value labels drawn in each cell. Default `2`.

- show_values:

  Logical; print the rounded variance in each cell. Default `TRUE`.

- low, high:

  Colours for the low (dark, highly proportional) and high (light) ends
  of the gradient. Defaults give a dark-to-light blue scale.

- theme:

  Optional `ggplot2` theme added over the default heatmap theme (e.g.,
  [`ggplot2::theme_bw()`](https://ggplot2.tidyverse.org/reference/ggtheme.html)
  or `ggplot2::theme(legend.position = "bottom")`). Default `NULL` adds
  none. Complete themes are listed in the `ggplot2` reference under
  "Themes"
  ([`?ggplot2::theme_bw`](https://ggplot2.tidyverse.org/reference/ggtheme.html)).
  [`plot_ternary()`](https://marcusvvlopes.github.io/timecoda/reference/plot_ternary.md)
  takes a separate `ggtern` theme.

## Value

Invisibly, a `ggplot` object. The call draws the heatmap as a side
effect.

## See also

[`summary.compD()`](https://marcusvvlopes.github.io/timecoda/reference/summary.compD.md)
for the variation matrix as a table.

Other visualization:
[`plot_ternary()`](https://marcusvvlopes.github.io/timecoda/reference/plot_ternary.md)

## Examples

``` r
cd <- declare_composition(timeuse, c("MVPA", "LPA", "SB", "SPT"),
                          closure = 1440)
plot_variation(cd)

plot_variation(cd, theme = ggplot2::theme_bw())

```
