# Search for the optimal composition (Goldilocks zone)

**\[experimental\]**

Finds the compositions with the best predicted outcome for one
[`comp_iv()`](https://marcusvvlopes.github.io/timecoda/reference/comp_iv.md)
fit, or one fit per outcome, among candidates on a grid in closure units
(by default every 10 minutes of a 1,440-minute day) inside a feasible
region of the sample's compositions. Each outcome's zone holds the
candidates whose prediction ranks in its best `1 - zone` share or meets
a `threshold`, and the closed geometric mean of the zones' intersection
is the optimal composition. The search method, its arguments and its
defaults may change in a later release.

## Usage

``` r
optimal_composition(
  object,
  direction = "max",
  zone = 0.95,
  region = "box",
  sd = 3,
  coverage = NULL,
  spacing = NULL,
  threshold = NULL,
  widen = FALSE,
  type = c("relative", "absolute"),
  category = NULL,
  at = NULL,
  by = NULL,
  from = NULL,
  conf_level = 0.95
)

# S3 method for class 'optimal_comp'
print(x, ...)

# S3 method for class 'optimal_comp'
plot(x, outcome = NULL, level = NULL, resolution = 40L, smooth = NULL, ...)
```

## Arguments

- object:

  A
  [`comp_iv()`](https://marcusvvlopes.github.io/timecoda/reference/comp_iv.md)
  fit, or a list of them, one per outcome, on the same parts (in the
  same order) and total. A list is named by outcome; an unnamed one
  takes each fit's outcome variable, the left-hand side of its model
  formula. The region is estimated from the first fit's composition.

- direction:

  `"max"` (default) when a higher predicted outcome is better, `"min"`
  when a lower one is. Give one value, or one per outcome. On a Cox fit
  (`"cox"`, `"coxme"`), `"max"` favours the highest hazard, and a search
  for the lowest risk uses `"min"`. On an ordinal fit (`"clm"`,
  `"clmm"`), the cumulative logit \\P(Y \le j) =
  \mathrm{logit}^{-1}(\theta_j - \eta)\\ moves probability to the higher
  outcome categories as the latent predictor \\\eta\\ rises. `"max"`
  therefore favours the higher categories. On a multinomial fit, `"max"`
  favours the `category` level over the reference level.

- zone:

  The rank a candidate needs to be in an outcome's zone, a number in (0,
  1); the default `0.95` keeps the best 5%.

- region:

  The feasible region. `"box"` (default) sets per-part limits (Dumuid et
  al. 2021; Brakenridge et al. 2024). `"ellipsoid"` is the chi-square
  ellipsoid in ILR coordinates of Mellow et al. (2026). `"observed"`
  takes the sample's own compositions (Lund Rasmussen et al. 2023;
  Vindiš et al. 2025). A data frame gives the compositions directly, one
  column per part on the closure scale; every row must be positive and
  sum to the fit's total.

- sd:

  For `"box"`: a positive number of standard deviations, default `3`.
  Each part is limited to its mean \\\pm\\ `sd` standard deviations on
  the closure scale, never beyond the part's observed range. Passing
  `sd` with another region, or with `coverage`, is an error, and
  `sd = NULL` needs `coverage`.

- coverage:

  For `"box"`, `NULL` (default) uses `sd`. Otherwise the box takes
  percentile limits in place of `sd`: either the quantiles `c(lo, hi)`
  or one share for the central quantiles
  `c((1 - coverage) / 2, (1 + coverage) / 2)`. For `"ellipsoid"`, the
  chi-square level of its boundary, a share in (0, 1), default `0.80`;
  the printed result also gives the share of the sample inside.
  `coverage` is ignored, with a message, for the other regions.

- spacing:

  For `"box"` and `"ellipsoid"`: the grid step, a positive number of
  closure units. `NULL` (default) uses 1/144 of the closure declared in
  [`declare_composition()`](https://marcusvvlopes.github.io/timecoda/reference/declare_composition.md),
  10 minutes of a 1,440-minute day. A finer grid holds more candidates,
  and the search stops when more than 1,000,000 lie inside the region,
  or more than 5,000,000 on the grid before the region is applied. The
  number of candidates grows steeply with the number of parts. With five
  or more parts the default search takes about a minute, and a coarser
  grid (e.g., 1/72 of the closure) shortens it.

- threshold:

  `NULL` (default), or a named list with one entry per outcome ruled by
  a value instead of `zone`. An entry is one number (a cut-off, met at
  or above it for `"max"` and at or below it for `"min"`) or two (an
  inclusive band; e.g., `list(zbmi = c(-1, 1))`). The value is on the
  `type` scale. For `"relative"`, it is the effect versus `from`: a
  ratio, and so positive, on a ratio scale such as an odds or hazard
  ratio, and otherwise a difference. For `"absolute"`, it is the
  predicted response.

- widen:

  `FALSE` (default) keeps every zone at `zone`. `TRUE` lowers `zone` in
  steps of 0.05 for the outcomes it rules, when their zones share no
  candidate, until they do. It never goes below 0.5 (the Intersection
  item of Details). With `by`, each level is widened on its own.

- type:

  `"relative"` (default) values each candidate by its effect versus
  `from` on the fit's effect scale, the measure
  [`predict_realloc()`](https://marcusvvlopes.github.io/timecoda/reference/predict_realloc.md)
  reports (e.g., an odds ratio, a hazard ratio or a mean difference).
  `"absolute"` values it by the predicted response, which the survival,
  ordinal and multinomial engines do not provide. The zone ruled by
  `zone` is the same under both. `threshold`, the zone limits, the
  `estimate` at the centre and
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) follow the
  scale.

- category:

  For a multinomial (`"mblogit"`) fit, the non-reference outcome level
  whose log-odds against the reference level is searched: required for
  such a fit, and refused when no fit is multinomial. `NULL` (default).

- at:

  Optional named list of covariate values held in every prediction, as
  in
  [`predict_realloc()`](https://marcusvvlopes.github.io/timecoda/reference/predict_realloc.md).

- by:

  Optional name of a factor covariate of every fit. The search then has
  one region and one zone per level, with the level held in its
  predictions. Defaults to a factor moderator not pinned by `at`, as in
  [`predict_realloc()`](https://marcusvvlopes.github.io/timecoda/reference/predict_realloc.md),
  and every fit must resolve to the same `by`.

- from:

  The composition the gain is measured from, and with it every value of
  a `type = "relative"` search. `NULL` (default) uses each fit's
  reference composition, stored by
  [`comp_iv()`](https://marcusvvlopes.github.io/timecoda/reference/comp_iv.md):
  the closed geometric mean of the composition the fit was given. Every
  level of `by` is measured from it. Otherwise, give a named numeric
  vector with one positive value per part, re-closed to the closure
  declared in
  [`declare_composition()`](https://marcusvvlopes.github.io/timecoda/reference/declare_composition.md).
  When the zone centre equals `from`, the gain is the reference value
  (0, or 1 on a ratio scale) with no test.

- conf_level:

  Confidence level of the intervals, default `0.95`.

- x:

  An `optimal_comp` object.

- ...:

  Unused by [`print()`](https://rdrr.io/r/base/print.html), present for
  S3 consistency.
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) takes no
  further arguments, and one given there is an error.

- outcome:

  The outcome whose surface is drawn; `NULL` (default) for the first.

- level:

  For a search with `by`, the level whose zone is drawn; `NULL`
  (default) for the first. Not used without `by`.

- resolution:

  For [`plot()`](https://rdrr.io/r/graphics/plot.default.html), the
  steps along each edge of a panel's simplex, as in
  [`plot.comp_iv()`](https://marcusvvlopes.github.io/timecoda/reference/plot.comp_iv.md)
  (default 40).

- smooth:

  `NULL` (default) draws the computed cells as tiles for a fit linear in
  the composition and a quadratic surface through them for any other
  fit; `TRUE` or `FALSE` draws the smooth or the tiles for any fit, as
  in
  [`plot.comp_iv()`](https://marcusvvlopes.github.io/timecoda/reference/plot.comp_iv.md).
  A panel whose smooth cannot be fitted is drawn as tiles, with a
  warning.

## Value

An `optimal_comp` object, a list with elements:

- `candidates`:

  A tibble with one row per candidate inside the region, holding the
  parts on the closure scale, the ILR coordinates (`ilr1`, `ilr2`, ...),
  `mahalanobis` (the distance for the ellipsoid, else `NA`) and the `by`
  column. Per outcome follow `estimate_<outcome>`, `rank_<outcome>` and
  `in_zone_<outcome>`, then the overall `in_zone`. `estimate_<outcome>`
  is on the link scale: the prediction, or for `type = "relative"` its
  difference from the prediction at `from`.

- `zones`:

  A tibble with one row per outcome and one `"overall"` row (per `by`
  level). Its columns are `rule` (`"quantile"` or `"threshold"`),
  `threshold_low`, `threshold_high`, `n`, per part `centre_<part>`,
  `min_<part>`, `max_<part>` and `best_<part>`, then `on_boundary`, and
  `q_overlap` on the overall row. The two thresholds are the limits the
  zone's values meet, on the `type` scale (the effect versus `from`, or
  the predicted response). `best_<part>` is the best-ranked candidate;
  for the overall zone, it is the candidate whose lowest rank across
  outcomes is highest.

- `optimum`:

  A tibble with one row per outcome (per `by` level). Its columns are
  the parts of the centre evaluated, `centre_of` (`"overall"`, or
  `"outcome"` when the zones do not overlap), and `estimate`, `conf.low`
  and `conf.high` at that centre on the `type` scale. Under `"relative"`
  the estimate is the effect versus `from` and equals the gain; under
  `"absolute"` it is the predicted response. The gain from `from`
  follows on the fit's effect scale (`gain`, `gain.low`, `gain.high`,
  `p.value`), then the `from` parts as `from_<part>`.

- `fits`:

  The fits, named by outcome, which
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) draws.

- Settings:

  `direction`, `zone` (as given), `threshold`, `widen`, `type`,
  `category`, `region`, the `coverage`, `sd` and `spacing` used (`NULL`
  where they do not apply), `n_candidates` (per level) and `zone_used`
  (per level, the `zone` searched, which is `zone` or the step `widen`
  reached). It also holds `parts`, `labels`, `units`, `total`,
  `outcomes`, `linear` (per outcome, `TRUE` for a fit with no curvature
  term), `measures` (the gain's effect label per outcome), `by`, `at`
  and `conf_level`. `region` holds its `name` and, per level, its
  estimated `parameters`. For the ellipsoid these include the limit `r2`
  on the squared distance and `share`, the share of the level's sample
  inside. For the box they are the per-part limits `lo` and `hi` and the
  flags `lo_clipped` and `hi_clipped`, `TRUE` where a limit was cut to
  the part's observed minimum or maximum;
  [`print()`](https://rdrr.io/r/base/print.html) reports those.

## Details

Every prediction comes from `emmeans`, on the reference grid
[`predict_realloc()`](https://marcusvvlopes.github.io/timecoda/reference/predict_realloc.md)
uses. With `type = "relative"` (default), a candidate's value is its
contrast with `from` (by default the reference composition) on the fit's
effect scale, as for a
[`predict_realloc()`](https://marcusvvlopes.github.io/timecoda/reference/predict_realloc.md)
reallocation, so every engine of
[`comp_iv()`](https://marcusvvlopes.github.io/timecoda/reference/comp_iv.md)
can be searched. With `type = "absolute"`, the value is the predicted
response.

The search runs in this order.

1.  Region. The feasible region is estimated from the sample's
    compositions, per level of `by`. The default `"box"` limits every
    part, on the closure scale, to its mean \\\pm\\ `sd` standard
    deviations, 3 by default (Dumuid et al. 2021). With `coverage`, the
    limits are a pair of the part's quantiles instead (Brakenridge et
    al. 2024 use the 1st and 99th percentiles). No limit extends beyond
    the part's observed range. The `"ellipsoid"` region holds the
    compositions whose squared Mahalanobis distance from the mean of the
    pivot isometric log-ratio (ILR) coordinates is at most the
    `coverage` quantile of the chi-square distribution on \\D - 1\\
    degrees of freedom (Mellow et al. 2026). That distance is \\(z -
    \bar z)^\top S^{-1} (z - \bar z)\\, with \\S\\ the sample
    covariance. `"observed"` takes the sample's own compositions as the
    candidates (Lund Rasmussen et al. 2023), and a data frame takes the
    compositions it holds.

2.  Candidates. For `"box"` and `"ellipsoid"`, the candidates are the
    compositions inside the region whose first `D - 1` parts are
    positive multiples of `spacing`, the last part being the remainder
    of the closure declared in
    [`declare_composition()`](https://marcusvvlopes.github.io/timecoda/reference/declare_composition.md).
    The default `spacing` is 1/144 of that closure, the 10-minute grid
    of Dumuid et al. (2021) on a 1,440-minute day; Mellow et al. (2026)
    use 5 minutes. The search stops with an error when more than
    1,000,000 candidates lie inside the region, or when the grid holds
    more than 5,000,000 compositions before the region is applied.

3.  Prediction. Each outcome is predicted at every candidate on the
    fit's link scale. For `type = "relative"`, each prediction is
    contrasted there with the prediction at `from`. As in
    [`predict_realloc()`](https://marcusvvlopes.github.io/timecoda/reference/predict_realloc.md),
    factor covariates that do not interact with the composition are
    averaged with equal weight, and numeric ones are held at their mean.
    Each `emmeans` reference grid holds many candidates, so a search
    over tens of thousands of them takes a few seconds.

4.  Zone. A candidate's rank is the share of candidates whose link-scale
    prediction is no better than its own (ties averaged), so it is the
    same for both values of `type`. An outcome's zone holds the
    candidates whose rank is at least `zone`. For an outcome named in
    `threshold`, a candidate is in the zone when its value on the `type`
    scale (the effect versus `from`, or the predicted response) meets
    the cut-off or lies in the band (Vindiš et al. 2025). In this
    package, the rule "healthy below the 85th percentile" of Lund
    Rasmussen et al. (2023) is `direction = "min", zone = 0.15`, which
    keeps 85% of the candidates.

5.  Intersection. The overall zone holds the candidates in every
    outcome's zone (Dumuid et al. 2021; Duncan et al. 2024). When the
    zones do not overlap, `q_overlap` gives the largest `zone` at which
    they would share a candidate, and the intersection stays empty
    unless `widen = TRUE`. Widening lowers `zone` for every outcome it
    rules, together, in steps of 0.05 (`zone - 0.05`, `zone - 0.10`,
    ...) until the zones share a candidate, as in Dumuid et al. (2021)
    and Brakenridge et al. (2024). It stops at the first such step and
    never goes below 0.5, the best 50% of each outcome. A `threshold` is
    never changed. The result is that of a search at the level reached,
    which `zone_used` reports. A widened zone may not improve every
    outcome, and each outcome's gain and its interval show whether it
    does.

6.  Centre. Each zone is summarized by its closed geometric mean (the
    inverse ILR of its mean coordinates), its per-part range and its
    best-ranked candidate.

7.  Gain. Each outcome is predicted at the overall zone centre (at its
    own centre when the zones do not overlap) and contrasted with
    `from`, as one
    [`predict_realloc()`](https://marcusvvlopes.github.io/timecoda/reference/predict_realloc.md)
    contrast. For `type = "relative"`, this gain is also the value
    reported at the centre.

A predicted outcome linear in the ILR coordinates has no interior
optimum, and its zone lies on the edge of the region. A fit with no
curvature term (squares, products,
[`poly()`](https://rdrr.io/r/stats/poly.html), `ns()`, `bs()`) is still
searched, with one message stating that the fit is linear in the
composition. The zone of a curved fit also lies on the boundary when its
stationary point is outside the region.

A zone's `on_boundary` value is the share of its candidates in the outer
shell of the region. On the ellipsoid, that shell holds the candidates
whose squared Mahalanobis distance is above 0.9 of the limit; on the
box, those with a part within one `spacing` of a limit.

The confidence intervals of the prediction at the centre and of the gain
are conditional on the centre, and the location of the zone itself has
no interval.

The box, the ellipsoid and the grid do not depend on the order of the
parts when `spacing` divides the total, as the default does.

A covariate that interacts with the composition must be pinned with `at`
or split with `by`; otherwise the call aborts and names it. A factor
moderator is split by default, as in
[`predict_realloc()`](https://marcusvvlopes.github.io/timecoda/reference/predict_realloc.md).
With `by`, each level has its own region, estimated from its own rows,
and its own zone. On a fit built with `include_total = TRUE`, every
prediction holds `t_total` at the fit's `total_ref`.

`type = "absolute"` is refused for the engines without an absolute
prediction (`"cox"`, `"coxme"`, `"clm"`, `"clmm"` and `"mblogit"`; see
[`comp_iv()`](https://marcusvvlopes.github.io/timecoda/reference/comp_iv.md)).
A multinomial (`"mblogit"`) fit is searched on the log-odds of the
outcome level named by `category` against the reference level, so each
`category` has its own ranking and zone.

## Functions

- `print(optimal_comp)`: Prints a text report in sections. The report
  gives the search's settings, each part's region limits, each zone's
  centre and per-part range, and the gain from `from` to the zone
  centre. The limits are the box's, or the candidates' range for another
  region, and a mark shows a box limit that is the part's observed
  minimum or maximum. With `by`, the limits, zones and gain are given
  for each level. A widened search states how far its zones were
  widened, and `type = "absolute"` adds the predicted outcome at each
  zone centre.

- `plot(optimal_comp)`: Ternary panels (`choose(D, 3)`) of one outcome's
  surface on the search's `type` scale (the effect versus `from`, or the
  predicted response), drawn for every engine. An ordinal fit's surface
  is its cumulative odds ratio versus `from`, and a multinomial fit's is
  the relative risk ratio of the searched `category`. Each panel is a
  slice through the zone centre (the overall one, or the outcome's own
  when the zones do not overlap). The surface is drawn as
  [`plot.comp_iv()`](https://marcusvvlopes.github.io/timecoda/reference/plot.comp_iv.md)
  draws one, on the cells inside the search's region. For the box and
  the ellipsoid these are the cells the region itself holds (on the
  parameters of the `level` drawn); otherwise they are the cells within
  the candidates' per-part range. The surface holds the search's `at`
  values and, for a search with `by`, the `level` drawn. The convex hull
  of the computed cells in the outcome's zone is outlined in red, and
  the centre is marked. A panel with fewer than three computed cells in
  the zone has no outline, and a higher `resolution` computes more
  cells. A panel carries no caption: the zone centre, the held
  covariates and the region's limits are in
  [`print()`](https://rdrr.io/r/base/print.html). Returns the plot
  object (one plot, or a list named by triplet), drawn as a side effect.
  An argument the method does not take is an error.

## References

Brakenridge, C. J., Koster, A., de Galan, B. E., Carver, A., Dumuid, D.,
Dzakpasu, F. Q. S., Eussen, S. J. P. M., Savelberg, H. H. C. M., Bosma,
H., Owen, N., Schaper, N. C., Healy, G. N., & Dunstan, D. W. (2024).
Associations of 24 h time-use compositions of sitting, standing,
physical activity and sleeping with optimal cardiometabolic risk and
glycaemic control: The Maastricht Study. *Diabetologia*, 67(7),
1356-1367.
[doi:10.1007/s00125-024-06145-0](https://doi.org/10.1007/s00125-024-06145-0)

Dumuid, D., Wake, M., Burgner, D., Tremblay, M. S., Okely, A. D.,
Edwards, B., Dwyer, T., & Olds, T. (2021). Balancing time use for
children's fitness and adiposity: Evidence to inform 24-hour guidelines
for sleep, sedentary time and physical activity. *PLoS ONE*, 16(1),
e0245501.
[doi:10.1371/journal.pone.0245501](https://doi.org/10.1371/journal.pone.0245501)

Duncan, M. J., Kuzik, N., Silva, D. A. S., Bélanger, R. E., Carson, V.,
Chaput, J.-P., Faulkner, G., Ferro, M. A., Turcotte-Tremblay, A.-M.,
Leatherdale, S. T., Patte, K. A., & Tremblay, M. S. (2024). "Goldilocks
days" for adolescent mental health: Movement behaviour combinations for
well-being, anxiety and depression by gender. *Mental Health and
Physical Activity*, 26, 100572.
[doi:10.1016/j.mhpa.2023.100572](https://doi.org/10.1016/j.mhpa.2023.100572)

Lund Rasmussen, C., Gába, A., Stanford, T., Dygrýn, J., Dumuid, D.,
Janda, D., & Hron, K. (2023). The Goldilocks Day for healthy adiposity
measures among children and adolescents. *Frontiers in Public Health*,
11, 1158634.
[doi:10.3389/fpubh.2023.1158634](https://doi.org/10.3389/fpubh.2023.1158634)

Mellow, M. L., Stanford, T. E., Olds, T., Miatke, A., Smith, A. E., &
Dumuid, D. (2026). An interactive tool to personalise 24-hour activity,
sitting and sleep prescription for optimal health outcomes. *npj Digital
Medicine*.
[doi:10.1038/s41746-026-02542-4](https://doi.org/10.1038/s41746-026-02542-4)

Vindiš, J., Pelclová, J., Jašková, P., & Hron, K. (2025). The Goldilocks
days represent optimal time-use to prevent obesity, low physical
performance, risk and fear of falling in older adults. *Scientific
Reports*.
[doi:10.1038/s41598-025-09645-0](https://doi.org/10.1038/s41598-025-09645-0)

## See also

[`comp_iv()`](https://marcusvvlopes.github.io/timecoda/reference/comp_iv.md)
to fit the models this consumes;
[`predict_realloc()`](https://marcusvvlopes.github.io/timecoda/reference/predict_realloc.md)
for a directed reallocation between named parts.

Other reallocation and prediction:
[`plot_realloc()`](https://marcusvvlopes.github.io/timecoda/reference/plot_realloc.md),
[`predict_realloc()`](https://marcusvvlopes.github.io/timecoda/reference/predict_realloc.md),
[`summary.realloc()`](https://marcusvvlopes.github.io/timecoda/reference/summary.realloc.md)

## Examples

``` r
# \donttest{
cd <- declare_composition(timeuse, c("MVPA", "LPA", "SB", "SPT"),
                          closure = 1440, units = "min/day")
quad <- comp_iv(cd, formula = zbmi ~ ilr()^2 + I(ilr()^2) + sex,
                model = "lm")
# The per-part box of the mean +/- 3 standard deviations, on a 10-minute grid
opt <- optimal_composition(quad, direction = "min")
opt
#> <optimal_comp> Goldilocks-zone search: zbmi
#> 
#> Search
#>   Parts       MVPA, LPA, SB, SPT (total = 1440 min/day)
#>   Region      per-part box of the mean +/- 3 SD
#>   Grid        multiples of 10 min/day; 36,611 candidates
#>   Scale       Mean difference versus the reference day
#>   Zone rule   best 5% of predictions, lower is better
#> 
#> Region limits (min/day)
#>           lower     upper
#>   MVPA      0.1*    180.0
#>   LPA      51.8*    503.3
#>   SB      315.6*   1013.4
#>   SPT     239.9*    809.4
#>   * the part's observed minimum
#> 
#> Zone for zbmi: 1,831 compositions (centre [min; max])
#>   MVPA   125.2  [ 30.0; 170.0]
#>   LPA    267.2  [110.0; 480.0]
#>   SB     353.9  [320.0; 410.0]
#>   SPT    693.6  [470.0; 800.0]
#> 
#> Gain from the reference day to the zone centre (95% CI)
#>   zbmi   Mean difference   -2.54  [-3.2; -1.88]  p = 5.44e-13

# Two outcomes, their zones intersected
risk <- comp_iv(cd, "obesity", model = "glm", family = "binomial")
optimal_composition(list(zbmi = quad, obesity = risk), direction = "min")
#> "obesity" is linear in the composition, so its best predicted outcome lies on
#> the edge of the region.
#> ℹ The zone then marks where the region ends; see on_boundary.
#> ℹ A curved surface can hold an interior optimum: refit with `curvature = TRUE`
#>   (or `ilr()^2 + I(ilr()^2)` in a formula).
#> <optimal_comp> Goldilocks-zone search: zbmi, obesity
#> 
#> Search
#>   Parts       MVPA, LPA, SB, SPT (total = 1440 min/day)
#>   Region      per-part box of the mean +/- 3 SD
#>   Grid        multiples of 10 min/day; 36,611 candidates
#>   Scale       Mean difference (zbmi), Odds ratio (obesity), versus the reference
#>               day
#>   Zone rule   best 5% of predictions, lower is better
#> 
#> Region limits (min/day)
#>           lower     upper
#>   MVPA      0.1*    180.0
#>   LPA      51.8*    503.3
#>   SB      315.6*   1013.4
#>   SPT     239.9*    809.4
#>   * the part's observed minimum
#> 
#> Zone for zbmi: 1,831 compositions (centre [min; max])
#>   MVPA   125.2  [ 30.0; 170.0]
#>   LPA    267.2  [110.0; 480.0]
#>   SB     353.9  [320.0; 410.0]
#>   SPT    693.6  [470.0; 800.0]
#> 
#> Zone for obesity: 1,831 compositions (centre [min; max])
#>   MVPA   151.2  [ 90.0; 170.0]
#>   LPA    262.2  [ 60.0; 500.0]
#>   SB     367.9  [320.0; 450.0]
#>   SPT    658.7  [400.0; 800.0]
#>   A fit linear in the composition has no interior optimum: the zone is the edge
#>   of the region in the direction of benefit.
#>   Refit with `curvature = TRUE` for a surface that can hold one.
#> 
#> Overall zone: 1,278 compositions (centre [min; max])
#>   MVPA   146.0  [ 90.0; 170.0]
#>   LPA    261.9  [110.0; 480.0]
#>   SB     353.9  [320.0; 410.0]
#>   SPT    678.2  [470.0; 800.0]
#>   A fit linear in the composition has no interior optimum: the zone is the edge
#>   of the region in the direction of benefit.
#>   Refit with `curvature = TRUE` for a surface that can hold one.
#> 
#> Gain from the reference day to the zone centre (95% CI)
#>   zbmi      Mean difference    -2.6  [-3.27; -1.93]  p = 3.32e-13
#>   obesity   Odds ratio        0.182  [0.079; 0.42]   p = 6.45e-05

# The chi-square ellipsoid at the 80% level
optimal_composition(quad, direction = "min", region = "ellipsoid")
#> <optimal_comp> Goldilocks-zone search: zbmi
#> 
#> Search
#>   Parts       MVPA, LPA, SB, SPT (total = 1440 min/day)
#>   Region      ILR ellipsoid at the 80% chi-square level, holding 82.7% of the
#>               sample
#>   Grid        multiples of 10 min/day; 18,814 candidates
#>   Scale       Mean difference versus the reference day
#>   Zone rule   best 5% of predictions, lower is better
#> 
#> Candidate range (min/day)
#>          lower   upper
#>   MVPA    10.0   300.0
#>   LPA    100.0   500.0
#>   SB     430.0   930.0
#>   SPT    290.0   750.0
#> 
#> Zone for zbmi: 941 compositions (centre [min; max])
#>   MVPA   207.3  [ 90.0; 300.0]
#>   LPA    225.6  [160.0; 320.0]
#>   SB     486.4  [430.0; 550.0]
#>   SPT    520.7  [390.0; 680.0]
#>   63% of the zone lies near the edge of the region (its outer shell); the best
#>   composition may lie beyond it.
#> 
#> Gain from the reference day to the zone centre (95% CI)
#>   zbmi   Mean difference   -1.89  [-2.27; -1.5]  p < 2e-16

# A box from the 1st to the 99th percentile of each part, on a 20-minute grid
optimal_composition(quad, direction = "min", coverage = c(0.01, 0.99),
                    spacing = 20)
#> <optimal_comp> Goldilocks-zone search: zbmi
#> 
#> Search
#>   Parts       MVPA, LPA, SB, SPT (total = 1440 min/day)
#>   Region      per-part box from the 1st to the 99th percentile
#>   Grid        multiples of 20 min/day; 3,243 candidates
#>   Scale       Mean difference versus the reference day
#>   Zone rule   best 5% of predictions, lower is better
#> 
#> Region limits (min/day)
#>          lower   upper
#>   MVPA     1.2   199.0
#>   LPA     68.0   422.0
#>   SB     390.9   925.0
#>   SPT    255.4   739.9
#> 
#> Zone for zbmi: 163 compositions (centre [min; max])
#>   MVPA   145.7  [ 60.0; 180.0]
#>   LPA    247.6  [120.0; 420.0]
#>   SB     424.4  [400.0; 460.0]
#>   SPT    622.3  [440.0; 720.0]
#>   63% of the zone lies near the edge of the region (its outer shell); the best
#>   composition may lie beyond it.
#> 
#> Gain from the reference day to the zone centre (95% CI)
#>   zbmi   Mean difference   -2.08  [-2.53; -1.63]  p < 2e-16

# A Cox fit, searched for the lowest hazard ratio versus the reference day
cox <- comp_iv(cd, formula = survival::Surv(time, event) ~ ilr()^2 +
                 I(ilr()^2) + sex, model = "cox")
optimal_composition(list(onset = cox), direction = "min")
#> <optimal_comp> Goldilocks-zone search: onset
#> 
#> Search
#>   Parts       MVPA, LPA, SB, SPT (total = 1440 min/day)
#>   Region      per-part box of the mean +/- 3 SD
#>   Grid        multiples of 10 min/day; 36,611 candidates
#>   Scale       Hazard ratio versus the reference day
#>   Zone rule   best 5% of predictions, lower is better
#> 
#> Region limits (min/day)
#>           lower     upper
#>   MVPA      0.1*    180.0
#>   LPA      51.8*    503.3
#>   SB      315.6*   1013.4
#>   SPT     239.9*    809.4
#>   * the part's observed minimum
#> 
#> Zone for onset: 1,831 compositions (centre [min; max])
#>   MVPA   151.4  [ 80.0; 170.0]
#>   LPA    339.9  [ 90.0; 500.0]
#>   SB     375.7  [320.0; 570.0]
#>   SPT    573.0  [240.0; 800.0]
#> 
#> Gain from the reference day to the zone centre (95% CI)
#>   onset   Hazard ratio   0.205  [0.0876; 0.48]  p = 0.000259
# }
```
