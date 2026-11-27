# Build publication-ready tables from a timecoda summary

Turns the object returned by a `timecoda`
[`summary()`](https://rdrr.io/r/base/summary.html) method into a named
list of [`gt::gt()`](https://gt.rstudio.com/reference/gt.html) tables,
with the numbers arranged in the layout of a manuscript or report table.

## Usage

``` r
pub_ready(x, digits = 2, p_flag = 0.05)
```

## Arguments

- x:

  A summary object from a `timecoda`
  [`summary()`](https://rdrr.io/r/base/summary.html) method: a
  `summary_compD` (from
  [`summary.compD()`](https://marcusvvlopes.github.io/timecoda/reference/summary.compD.md)),
  a `summary_comp_iv` (from
  [`summary.comp_iv()`](https://marcusvvlopes.github.io/timecoda/reference/summary.comp_iv.md))
  or a `summary_realloc` (from
  [`summary.realloc()`](https://marcusvvlopes.github.io/timecoda/reference/summary.realloc.md)).

- digits:

  Number of decimal places in the rendered table: a single non-negative
  whole number. Default `2`.

- p_flag:

  Significance threshold below which the `summary_realloc` table flags
  an effect with `*` (relative output only). A single number from 0 to
  1; default `0.05`. At `0` no effect is flagged. The other summary
  types ignore it.

## Value

A named list of [`gt::gt()`](https://gt.rstudio.com/reference/gt.html)
tables, one per table the summary carries, with no entry for a table it
lacks.

- `summary_compD`:

  `desc` (each part's arithmetic Mean \\\pm\\ SD of the original data,
  before zero imputation and closure, alongside its geometric mean on
  the closed composition and, when the composition is not already closed
  to proportions or percentages, its share of the total; each column
  header names its scale) and `var_matrix`.

- `summary_comp_iv`:

  `effects` (the omnibus composition row followed by one "part vs
  remaining" row per part, with the effect, its confidence interval and
  p-value; a per-category model adds an Outcome column and a moderated
  fit a per-part moderation p-value).

- `summary_realloc`:

  `effects` (one row per reallocation pair, with the effect, flagged `*`
  below `p_flag`, and its confidence interval under a column spanner for
  each `delta`; a `by` grouping or a multinomial outcome adds a column).

## See also

[`summary.compD()`](https://marcusvvlopes.github.io/timecoda/reference/summary.compD.md),
[`summary.comp_iv()`](https://marcusvvlopes.github.io/timecoda/reference/summary.comp_iv.md)
and
[`summary.realloc()`](https://marcusvvlopes.github.io/timecoda/reference/summary.realloc.md)
for the underlying summaries.

## Examples

``` r
cd <- declare_composition(timeuse, c("MVPA", "LPA", "SB", "SPT"),
                          closure = 1440, units = "min/day")
s <- summary(cd)
pub_ready(s)$desc


  


Descriptive summary
```
