# Package index

## Declare & describe composition

Build the `compD` object, with validation, zero imputation and optional
compositional outlier flagging, and describe it with tidy tables and
plots.

- [`declare_composition()`](https://marcusvvlopes.github.io/timecoda/reference/declare_composition.md)
  : Declare a compositional variable
- [`comp_impute()`](https://marcusvvlopes.github.io/timecoda/reference/comp_impute.md)
  : Impute compositional zeroes
- [`comp_outliers()`](https://marcusvvlopes.github.io/timecoda/reference/comp_outliers.md)
  : Detect compositional outliers
- [`print(`*`<compD>`*`)`](https://marcusvvlopes.github.io/timecoda/reference/print.compD.md)
  : Print a compositional variable
- [`summary(`*`<compD>`*`)`](https://marcusvvlopes.github.io/timecoda/reference/summary.compD.md)
  : Describe a compositional variable
- [`print(`*`<summary_compD>`*`)`](https://marcusvvlopes.github.io/timecoda/reference/print.summary_compD.md)
  : Print a compositional variable summary
- [`plot(`*`<compD>`*`)`](https://marcusvvlopes.github.io/timecoda/reference/plot.compD.md)
  : Visualize a compositional variable
- [`pub_ready()`](https://marcusvvlopes.github.io/timecoda/reference/pub_ready.md)
  : Build publication-ready tables from a timecoda summary

## Composition as predictor

Model the composition as an independent variable across eleven model
engines, with a publication effect table and a predicted-outcome plot.

- [`comp_iv()`](https://marcusvvlopes.github.io/timecoda/reference/comp_iv.md)
  : Model a composition as an independent variable

- [`summary(`*`<comp_iv>`*`)`](https://marcusvvlopes.github.io/timecoda/reference/summary.comp_iv.md)
  : Summarize a composition-as-predictor model

- [`print(`*`<summary_comp_iv>`*`)`](https://marcusvvlopes.github.io/timecoda/reference/print.summary_comp_iv.md)
  : Print a composition-as-predictor summary

- [`plot(`*`<comp_iv>`*`)`](https://marcusvvlopes.github.io/timecoda/reference/plot.comp_iv.md)
  : Visualize a composition-as-predictor model

- [`print(`*`<timecoda_fit>`*`)`](https://marcusvvlopes.github.io/timecoda/reference/print.timecoda_fit.md)
  : Print a timecoda model fit

- [`ilr()`](https://marcusvvlopes.github.io/timecoda/reference/ilr.md) :

  Composition marker for
  [`comp_iv()`](https://marcusvvlopes.github.io/timecoda/reference/comp_iv.md)
  formulas

## Reallocation & prediction

Estimate the effects of time reallocations from a
[`comp_iv()`](https://marcusvvlopes.github.io/timecoda/reference/comp_iv.md)
fit, plot them, tabulate them for publication, and search the fit for an
optimal composition.

- [`predict_realloc()`](https://marcusvvlopes.github.io/timecoda/reference/predict_realloc.md)
  : Predict time-reallocation effects
- [`summary(`*`<realloc>`*`)`](https://marcusvvlopes.github.io/timecoda/reference/summary.realloc.md)
  : Summarize time-reallocation effects
- [`print(`*`<summary_realloc>`*`)`](https://marcusvvlopes.github.io/timecoda/reference/print.summary_realloc.md)
  : Print a time-reallocation summary
- [`optimal_composition()`](https://marcusvvlopes.github.io/timecoda/reference/optimal_composition.md)
  [`print(`*`<optimal_comp>`*`)`](https://marcusvvlopes.github.io/timecoda/reference/optimal_composition.md)
  [`plot(`*`<optimal_comp>`*`)`](https://marcusvvlopes.github.io/timecoda/reference/optimal_composition.md)
  **\[experimental\]** : Search for the optimal composition (Goldilocks
  zone)
- [`plot_realloc()`](https://marcusvvlopes.github.io/timecoda/reference/plot_realloc.md)
  [`print(`*`<realloc_plots>`*`)`](https://marcusvvlopes.github.io/timecoda/reference/plot_realloc.md)
  : Plot a time-reallocation sweep

## Visualization

Standalone ternary and variation-matrix plots.

- [`plot_ternary()`](https://marcusvvlopes.github.io/timecoda/reference/plot_ternary.md)
  : Ternary plot(s) of a composition
- [`plot_variation()`](https://marcusvvlopes.github.io/timecoda/reference/plot_variation.md)
  : Heatmap of the compositional variation matrix

## Diagnostics

Goodness-of-fit and residual diagnostics matched to the fitted model
class.

- [`fit_check()`](https://marcusvvlopes.github.io/timecoda/reference/fit_check.md)
  : Check the goodness of fit of a composition model

## Composition helpers

Validation, the ILR transform and its inverse, and other composition
utilities.

- [`check_composition()`](https://marcusvvlopes.github.io/timecoda/reference/check_composition.md)
  : Validate a compositional dataset
- [`ilr_transform()`](https://marcusvvlopes.github.io/timecoda/reference/ilr_transform.md)
  : Isometric log-ratio (pivot) coordinates of a composition
- [`ilr_backtransform()`](https://marcusvvlopes.github.io/timecoda/reference/ilr_backtransform.md)
  : Back-transform ILR coordinates to a composition
- [`geom_mean_comp()`](https://marcusvvlopes.github.io/timecoda/reference/geom_mean_comp.md)
  : Geometric mean composition
- [`variation_matrix()`](https://marcusvvlopes.github.io/timecoda/reference/variation_matrix.md)
  : Compositional variation matrix

## Data

Example dataset shipped with the package.

- [`timeuse`](https://marcusvvlopes.github.io/timecoda/reference/timeuse.md)
  : Simulated 24-hour time-use dataset

## Package overview

- [`timecoda`](https://marcusvvlopes.github.io/timecoda/reference/timecoda-package.md)
  [`timecoda-package`](https://marcusvvlopes.github.io/timecoda/reference/timecoda-package.md)
  : timecoda: Compositional Data Analysis for Time-Use Behaviour
  Research
