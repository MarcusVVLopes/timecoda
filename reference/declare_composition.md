# Declare a compositional variable

Builds a `compD` object from a set of composition parts. It validates
the parts (zeroes permitted), summarizes zero occurrence, imputes any
zeroes, optionally flags compositional outliers and stores the closed
composition.

## Usage

``` r
declare_composition(
  data,
  parts,
  closure = 1,
  labels = NULL,
  units = NULL,
  impute_method = "lrEM",
  detect_outliers = FALSE,
  outlier_args = list(),
  ...
)
```

## Arguments

- data:

  A data frame containing the composition parts.

- parts:

  Character vector of column names making up the composition (at least
  two parts). Required; there is no default.

- closure:

  Scale the composition is reported on. Either a positive number (e.g.,
  `1440` for minutes-per-day, `1` for proportions) or one of the strings
  `"percent"` (equivalent to `closure = 100`) or `"mean"` (the average
  of the observed row totals of the original input parts). The resolved
  numeric value is stored as `compD$closure`, and every downstream
  summary, plot and prediction uses it. Defaults to `1` (unit
  proportions).

- labels:

  Optional display labels for the parts, reported in summary and plot
  output in place of the raw part names. Either a named character vector
  mapping part names to labels (e.g., `c(MVPA = "Exercise (min/day)")`;
  parts not named keep their raw name) or an unnamed character vector of
  length `length(parts)` giving a label per part in order. Labels must
  be unique across parts (a label may not repeat another part's label or
  raw name). Stored on the object as a complete named vector keyed by
  every part; `NULL` (default) means the raw part names are used
  throughout.

- units:

  Optional unit label for the closure scale, a single string such as
  `"min/day"` for a 1440-minute day. The print methods show it verbatim
  after the closure and after values reported on that scale. It also
  appears in the
  [`summary.compD()`](https://marcusvvlopes.github.io/timecoda/reference/summary.compD.md)
  and
  [`summary.realloc()`](https://marcusvvlopes.github.io/timecoda/reference/summary.realloc.md)
  publication tables
  ([`pub_ready()`](https://marcusvvlopes.github.io/timecoda/reference/pub_ready.md)).
  `NULL` (default) shows no unit.

- impute_method:

  Zero-imputation method, passed to
  [`comp_impute()`](https://marcusvvlopes.github.io/timecoda/reference/comp_impute.md)
  (ignored when the composition has no zeroes). Defaults to `"lrEM"`
  (log-ratio EM); see
  [`comp_impute()`](https://marcusvvlopes.github.io/timecoda/reference/comp_impute.md)
  for the alternatives.

- detect_outliers:

  Single logical; flag compositional outliers
  ([`comp_outliers()`](https://marcusvvlopes.github.io/timecoda/reference/comp_outliers.md))
  on the zero-free parts. Default `FALSE`.

- outlier_args:

  A named list of further arguments for
  [`comp_outliers()`](https://marcusvvlopes.github.io/timecoda/reference/comp_outliers.md),
  most usefully `quantile` (the chi-squared cutoff probability) or
  `alpha` (MCD's subset proportion). With `detect_outliers = FALSE` it
  is only checked to be a named list. With `TRUE`, a name that
  [`comp_outliers()`](https://marcusvvlopes.github.io/timecoda/reference/comp_outliers.md)
  does not accept is an error.

- ...:

  Additional arguments forwarded to
  [`comp_impute()`](https://marcusvvlopes.github.io/timecoda/reference/comp_impute.md),
  and from there to the chosen method. The most useful is `dl`, the
  per-part detection limit below which an imputed value must fall. By
  default `declare_composition()` derives it from each part's smallest
  observed positive value on the original scale. Override it with
  `dl = <vector or matrix of limits>`, or pass other method-specific
  controls such as `lrEM()`'s `rob`, `tolerance` or `max.iter`.

## Value

An object of class `compD`: a list with elements

- `data_raw`:

  The original input data (a tibble), unchanged, so it may still contain
  zeroes. With `detect_outliers = TRUE` it also carries three appended
  columns: `.outlier_distance` (the Mahalanobis distance in ILR
  coordinates), `.outlier_limit` (the cutoff, constant down the column)
  and `.outlier` (logical).

- `parts`:

  The character vector of part names.

- `closure`:

  The resolved numeric closure (always a positive number, even when the
  user supplied `"percent"` or `"mean"`).

- `composition`:

  The imputed, closed composition, a
  [`compositions::acomp()`](https://rdrr.io/pkg/compositions/man/acomp.html)
  object, from which
  [`summary.compD()`](https://marcusvvlopes.github.io/timecoda/reference/summary.compD.md)
  and
  [`plot.compD()`](https://marcusvvlopes.github.io/timecoda/reference/plot.compD.md)
  compute everything they show.

- `parts_imputed`:

  The zero-free parts on the original (pre-closure) scale, a tibble: the
  input values with any zeroes replaced by
  [`comp_impute()`](https://marcusvvlopes.github.io/timecoda/reference/comp_impute.md).
  It keeps each row's absolute total, which the `include_total`
  (Coenders T-space) covariate of
  [`comp_iv()`](https://marcusvvlopes.github.io/timecoda/reference/comp_iv.md)
  requires.

- `imput_table`:

  A tibble describing zero occurrence per part: `part`, `n`, `n_zero`,
  `pct_zero`, `min_nonzero`.

- `labels`:

  Either `NULL` or a complete named character vector keyed by every
  part, giving the display label used in summaries and plots.

- `units`:

  The `units` argument as given, or `NULL`.

- `impute_method`:

  The `impute_method` argument, or `NA_character_` when no part had a
  zero (the method was requested but never run).
  [`summary.compD()`](https://marcusvvlopes.github.io/timecoda/reference/summary.compD.md)
  reports it as the method used.

## Details

Compositional data carry only relative information, and a zero has no
log-ratio. When zeroes are present,
[`comp_impute()`](https://marcusvvlopes.github.io/timecoda/reference/comp_impute.md)
imputes them on the original (pre-closure) scale before
[`compositions::acomp()`](https://rdrr.io/pkg/compositions/man/acomp.html)
closes the result. Every downstream description, plot and coordinate is
computed from this imputed, closed composition. The zero-free values on
the original scale keep each row's absolute total and are stored as
`parts_imputed` for the `include_total` T-space covariate of
[`comp_iv()`](https://marcusvvlopes.github.io/timecoda/reference/comp_iv.md).
Without zeroes, no imputation is performed and `parts_imputed` equals
the input parts.

## Outlier flagging

With `detect_outliers = TRUE`,
[`comp_outliers()`](https://marcusvvlopes.github.io/timecoda/reference/comp_outliers.md)
runs on the zero-free parts after imputation, and its per-observation
columns are appended to `data_raw` as `.outlier_distance`,
`.outlier_limit` and `.outlier`. Rows are flagged and never dropped.
Every field of the object keeps one entry per input observation, and
`data_raw`, `composition` and `parts_imputed` stay row-aligned with the
input and with each other.

Filtering is left to the caller. Dropping the flagged rows means
declaring the composition again, which recomputes it (and
`closure = "mean"`) on the remaining rows. The `.outlier` columns must
be dropped in the same step. `declare_composition()` writes them itself
and stops with an error on any input that already carries them, whether
or not detection is requested.

    kept <- cd$data_raw |>
      dplyr::filter(!.outlier) |>
      dplyr::select(-dplyr::starts_with(".outlier"))

    cd_kept <- declare_composition(kept, parts, closure = 1440)

## References

Palarea-Albaladejo, J., & Martín-Fernández, J.A. (2015). zCompositions –
R package for multivariate imputation of left-censored data under a
compositional approach. *Chemometrics and Intelligent Laboratory
Systems*, 143, 85-96.
[doi:10.1016/j.chemolab.2015.02.019](https://doi.org/10.1016/j.chemolab.2015.02.019)

## See also

[`summary.compD()`](https://marcusvvlopes.github.io/timecoda/reference/summary.compD.md)
and
[`plot.compD()`](https://marcusvvlopes.github.io/timecoda/reference/plot.compD.md)
to describe and visualize the object;
[`comp_impute()`](https://marcusvvlopes.github.io/timecoda/reference/comp_impute.md)
and
[`comp_outliers()`](https://marcusvvlopes.github.io/timecoda/reference/comp_outliers.md)
for the zero-imputation and outlier-detection steps on their own;
[`comp_iv()`](https://marcusvvlopes.github.io/timecoda/reference/comp_iv.md)
to model it as a predictor of an outcome.

Other declare and describe composition:
[`comp_impute()`](https://marcusvvlopes.github.io/timecoda/reference/comp_impute.md),
[`comp_outliers()`](https://marcusvvlopes.github.io/timecoda/reference/comp_outliers.md),
[`plot.compD()`](https://marcusvvlopes.github.io/timecoda/reference/plot.compD.md),
[`summary.compD()`](https://marcusvvlopes.github.io/timecoda/reference/summary.compD.md)

## Examples

``` r
cd <- declare_composition(
  timeuse,
  parts = c("MVPA", "LPA", "SB", "SPT"),
  closure = 1440,
  units = "min/day"
)
cd$imput_table
#> # A tibble: 4 × 5
#>   part      n n_zero pct_zero min_nonzero
#>   <chr> <int>  <int>    <dbl>       <dbl>
#> 1 MVPA    300      0        0         0.1
#> 2 LPA     300      0        0        51.8
#> 3 SB      300      0        0       316. 
#> 4 SPT     300      0        0       240. 

set.seed(1)
flagged <- declare_composition(
  timeuse,
  parts = c("MVPA", "LPA", "SB", "SPT"),
  closure = 1440,
  detect_outliers = TRUE
)
sum(flagged$data_raw$.outlier)
#> [1] 14
```
