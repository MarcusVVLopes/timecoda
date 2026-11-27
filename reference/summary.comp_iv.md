# Summarize a composition-as-predictor model

Tabulates the omnibus test of the whole composition of a
[`comp_iv()`](https://marcusvvlopes.github.io/timecoda/reference/comp_iv.md)
fit and one "part vs remaining" effect per part. Each effect is the
rotated model's `ilr1` coefficient, or its slope at the reference
composition (see Details), back-transformed to the model's effect scale.
Printing the result shows a text report, as with
[`summary.lm()`](https://rdrr.io/r/stats/summary.lm.html), and
[`pub_ready()`](https://marcusvvlopes.github.io/timecoda/reference/pub_ready.md)
builds a publication `gt` table from it.

## Usage

``` r
# S3 method for class 'comp_iv'
summary(object, conf_level = 0.95, ...)
```

## Arguments

- object:

  A `comp_iv` object from
  [`comp_iv()`](https://marcusvvlopes.github.io/timecoda/reference/comp_iv.md).

- conf_level:

  Confidence level for the per-part Wald intervals.

- ...:

  Unused; present for S3 consistency.

## Value

An object of class `summary_comp_iv`, printed as a text report by
[`print.summary_comp_iv()`](https://marcusvvlopes.github.io/timecoda/reference/print.summary_comp_iv.md)
and passed to
[`pub_ready()`](https://marcusvvlopes.github.io/timecoda/reference/pub_ready.md)
for a `gt` table. A list with elements:

- `global`:

  The one-row omnibus composition test.

- `per_part`:

  The per-part effect table: one row per part, or per part and
  non-reference outcome category for a multinomial fit. It shares the
  columns `parameter`, `label`, `estimate`, `conf.low`, `conf.high`,
  `p.value` with `global`, whose `estimate`/`conf.low`/`conf.high` are
  `NA`. `label` carries the part's display label (or its raw name). A
  per-category model (multinomial) adds `response_level`.

- `moderation`:

  The joint moderation Wald test when the fit was built with a
  `moderator`, else `NULL`. A moderated fit also gives `per_part` three
  further columns (`mod_statistic`, `mod_df`, `mod_p.value`). They hold
  the joint Wald test of that part's own `ilr1 : moderator`
  coefficient(s), which asks whether the part's effect varies by
  moderator level; the `moderation` test pools every part. Like the
  per-part effect rows, they are descriptive and carry no multiplicity
  adjustment. An unmoderated fit omits them.

- `curvature_test`:

  The joint Wald test of the squares and products of the coordinates
  when the fit holds the full quadratic (see the `curvature_test`
  element of
  [`comp_iv()`](https://marcusvvlopes.github.io/timecoda/reference/comp_iv.md)),
  else `NULL`.

- `convergence`:

  A list of the fit's `status` one-row tibble (`converged`, `singular`,
  `n_warnings`) and its `warnings` tibble of captured backend
  conditions.

- `icc`:

  For a clustered fit, a tibble with one row per grouping term (`group`,
  `variance`, `measure`, `value`), read from the fit's own variance
  components. For `"lmer"` and `"rlmer"`, `value` is each term's share
  of the total variance, the intraclass correlation \\v_g / (\sum_h
  v_h + \sigma^2)\\, where the sum runs over all grouping terms \\h\\. A
  logit- or probit-link binomial `"glmer"` and an ordinal `"clmm"` give
  its latent-scale form, with \\\pi^2/3\\ (logit) or 1 (probit) in place
  of \\\sigma^2\\. A `"coxme"` fit gives the median hazard ratio
  \\\exp(\sqrt{2 v_g}\\ \Phi^{-1}(0.75))\\: the median, over pairs of
  clusters, of the hazard ratio between two individuals with identical
  covariates, one from each cluster, taken from the higher-risk cluster
  to the lower-risk one (Austin et al. 2017). `measure` is
  `"not available"` and `value` `NA` for the other clustered engines,
  for a log-link binomial `"glmer"`, which has no latent-scale form, and
  for a `"glmer"` of any other family. `NULL` for a fit without a random
  term.

- `total_test`, `total_beta`:

  With `include_total = TRUE`, the one-row Wald test of the `t_total`
  coefficient and its fitted value on the link scale (the first
  category's for a multinomial fit, whose test is joint); `NULL`
  otherwise.

- `per_part_slopes`:

  Present, and `TRUE`, only when at least one per-part row is a slope at
  the reference composition rather than a coefficient.

- `per_part_no_slope`:

  The reason, as text, why
  [`emmeans::emtrends()`](https://rvlenth.github.io/emmeans/reference/emtrends.html)
  gave no slope and those rows are `NA`. Present only in that case.

- `labels`, `effect_measure`, `conf_level`, `moderator`, `fit_label`,
  `total`, `model_type`, `omnibus`:

  Carried through for
  [`print.summary_comp_iv()`](https://marcusvvlopes.github.io/timecoda/reference/print.summary_comp_iv.md)
  and
  [`pub_ready()`](https://marcusvvlopes.github.io/timecoda/reference/pub_ready.md).

- `units`:

  The composition's unit label, printed after `total`; `NULL` when none
  was declared.

## Details

The effect scale is a hazard ratio for a Cox model, a cumulative odds
ratio for an ordinal model and a relative-risk ratio for a multinomial
model. A GLM's scale follows its link: an odds ratio (logit), a rate,
risk or mean ratio (log) or a mean difference (identity). Every linear
engine reports a mean difference.

For mixed models the per-part effects come from one refit per part
(stored at fit time in `object$models`), so building this summary does
no extra fitting. For an `"lmer"` fit, each per-part standard error,
interval and p-value uses the degrees-of-freedom method that `emmeans`
uses for
[`predict_realloc()`](https://marcusvvlopes.github.io/timecoda/reference/predict_realloc.md)
on the same fit. The default is the Kenward-Roger approximation through
pbkrtest. The Satterthwaite approximation through lmerTest is used when
`emmeans::emm_options(lmer.df =)` selects it, or when pbkrtest is
absent, disabled or the fit has more observations than
`emmeans::emm_options(pbkrtest.limit =)`. Both pbkrtest and lmerTest are
optional, and when neither applies, or
`emmeans::emm_options(lmer.df = "asymptotic")` is set, the degrees of
freedom are asymptotic (normal). `"glmer"`, `"coxme"` and `"rlmer"` rows
are asymptotic, as their `emmeans` predictions are (robustlmm reports no
degrees of freedom).

When the formula gives the composition non-linear terms (see
[`ilr()`](https://marcusvvlopes.github.io/timecoda/reference/ilr.md)), a
product of coordinates or no main effect, a rotated fit has no single
`ilr1` coefficient to read. Each per-part effect is then the slope of
its linear predictor in `ilr1` at the reference composition, from
[`emmeans::emtrends()`](https://rvlenth.github.io/emmeans/reference/emtrends.html),
back-transformed in the same way, and the printed report says so. With a
composition-by-covariate interaction as well, that slope is the
equal-weight average over the covariate's levels. An engine for which
`emtrends()` gives no such slope (`"mblogit"`, `"lmrob"`) reports `NA`
rows.

## References

Austin, P. C., Wagner, P., & Merlo, J. (2017). The median hazard ratio:
A useful measure of variance and general contextual effects in
multilevel survival analysis. *Statistics in Medicine*, 36(6), 928-938.
[doi:10.1002/sim.7188](https://doi.org/10.1002/sim.7188)

Goldstein, H., Browne, W., & Rasbash, J. (2002). Partitioning variation
in multilevel models. *Understanding Statistics*, 1(4), 223-231.
[doi:10.1207/S15328031US0104_02](https://doi.org/10.1207/S15328031US0104_02)

Kenward, M. G., & Roger, J. H. (1997). Small sample inference for fixed
effects from restricted maximum likelihood. *Biometrics*, 53(3),
983-997. [doi:10.2307/2533558](https://doi.org/10.2307/2533558)

## See also

[`comp_iv()`](https://marcusvvlopes.github.io/timecoda/reference/comp_iv.md)
to fit the model;
[`plot.comp_iv()`](https://marcusvvlopes.github.io/timecoda/reference/plot.comp_iv.md)
for the predicted-outcome surface;
[`predict_realloc()`](https://marcusvvlopes.github.io/timecoda/reference/predict_realloc.md)
for reallocation predictions;
[`pub_ready()`](https://marcusvvlopes.github.io/timecoda/reference/pub_ready.md)
for a formatted `gt` table built from this summary.

Other composition as predictor:
[`comp_iv()`](https://marcusvvlopes.github.io/timecoda/reference/comp_iv.md),
[`plot.comp_iv()`](https://marcusvvlopes.github.io/timecoda/reference/plot.comp_iv.md)

## Examples

``` r
cd <- declare_composition(timeuse, c("MVPA", "LPA", "SB", "SPT"),
                          closure = 1440)
m <- comp_iv(cd, outcome = "obesity", model = "glm", family = "binomial")
s <- summary(m)
s$per_part
#> # A tibble: 4 × 6
#>   parameter         label estimate conf.low conf.high     p.value
#>   <chr>             <chr>    <dbl>    <dbl>     <dbl>       <dbl>
#> 1 MVPA vs remaining MVPA     0.501    0.382     0.656 0.000000489
#> 2 LPA vs remaining  LPA      0.982    0.506     1.91  0.958      
#> 3 SB vs remaining   SB       2.98     1.16      7.65  0.0234     
#> 4 SPT vs remaining  SPT      0.683    0.273     1.71  0.415      
```
