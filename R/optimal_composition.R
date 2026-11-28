# Module 2: the Goldilocks-zone search (optimal_composition.R)
#
# The Goldilocks-zone search: a grid of candidate compositions in closure units
# inside a feasible region, predicted on a few emmeans grids per outcome, the
# best share of them per outcome intersected across outcomes, and the zone's
# centre compared with a chosen day.

#' Search for the optimal composition (Goldilocks zone)
#'
#' @description
#' `r lifecycle::badge("experimental")`
#'
#' Finds the compositions with the best predicted outcome for one [comp_iv()]
#' fit, or one fit per outcome, among candidates on a grid in closure units (by
#' default every 10 minutes of a 1,440-minute day) inside a feasible region of
#' the sample's compositions. Each outcome's zone holds the candidates whose
#' prediction ranks in its best `1 - zone` share or meets a `threshold`, and
#' the closed geometric mean of the zones' intersection is the optimal
#' composition. The search method, its arguments and its defaults may change
#' in a later release.
#'
#' @details
#' Every prediction comes from `emmeans`, on the reference grid
#' [predict_realloc()] uses. With `type = "relative"` (default), a candidate's
#' value is its contrast with `from` (by default the reference composition) on
#' the fit's effect scale, as for a [predict_realloc()] reallocation, so every
#' engine of [comp_iv()] can be searched. With `type = "absolute"`, the value is
#' the predicted response.
#'
#' The search runs in this order.
#'
#' 1. Region. The feasible region is estimated from the sample's
#'    compositions, per level of `by`. The default `"box"` limits every part,
#'    on the closure scale, to its mean \eqn{\pm}{+/-} `sd` standard
#'    deviations, 3 by default (Dumuid et al. 2021). With `coverage`, the
#'    limits are a pair of the part's quantiles instead (Brakenridge et al.
#'    2024 use the 1st and 99th percentiles). No limit extends beyond the
#'    part's observed range. The `"ellipsoid"` region holds the compositions
#'    whose squared Mahalanobis distance from the mean of the pivot isometric
#'    log-ratio (ILR) coordinates is at most the `coverage` quantile of the
#'    chi-square distribution on \eqn{D - 1} degrees of freedom (Mellow et al.
#'    2026). That distance is
#'    \eqn{(z - \bar z)^\top S^{-1} (z - \bar z)}{(z - zbar)' S^-1 (z - zbar)},
#'    with \eqn{S} the sample covariance. `"observed"` takes the sample's own
#'    compositions as the candidates (Lund Rasmussen et al. 2023), and a data
#'    frame takes the compositions it holds.
#' 2. Candidates. For `"box"` and `"ellipsoid"`, the candidates are the
#'    compositions inside the region whose first `D - 1` parts are positive
#'    multiples of `spacing`, the last part being the remainder of the closure
#'    declared in [declare_composition()]. The default `spacing` is 1/144 of
#'    that closure, the 10-minute grid of Dumuid et al. (2021) on a
#'    1,440-minute day; Mellow et al. (2026) use 5 minutes. The search stops
#'    with an error when more than 1,000,000 candidates lie inside the region,
#'    or when the grid holds more than 5,000,000 compositions before the
#'    region is applied.
#' 3. Prediction. Each outcome is predicted at every candidate on the fit's
#'    link scale. For `type = "relative"`, each prediction is contrasted there
#'    with the prediction at `from`. As in [predict_realloc()], factor
#'    covariates that do not interact with the composition are averaged with
#'    equal weight, and numeric ones are held at their mean. Each `emmeans`
#'    reference grid holds many candidates, so a search over tens of thousands
#'    of them takes a few seconds.
#' 4. Zone. A candidate's rank is the share of candidates whose link-scale
#'    prediction is no better than its own (ties averaged), so it is the same
#'    for both values of `type`. An outcome's zone holds the candidates whose
#'    rank is at least `zone`. For an outcome named in `threshold`, a candidate
#'    is in the zone when its value on the `type` scale (the effect versus
#'    `from`, or the predicted response) meets the cut-off or lies in the band
#'    (Vindiš et al. 2025). In this package, the rule "healthy below the 85th
#'    percentile" of Lund Rasmussen et al. (2023) is
#'    `direction = "min", zone = 0.15`, which keeps 85% of the candidates.
#' 5. Intersection. The overall zone holds the candidates in every outcome's
#'    zone (Dumuid et al. 2021; Duncan et al. 2024). When the zones do not
#'    overlap, `q_overlap` gives the largest `zone` at which they would share
#'    a candidate, and the intersection stays empty unless `widen = TRUE`.
#'    Widening lowers `zone` for every outcome it rules, together, in steps of
#'    0.05 (`zone - 0.05`, `zone - 0.10`, ...) until the zones share a
#'    candidate, as in Dumuid et al. (2021) and Brakenridge et al. (2024). It
#'    stops at the first such step and never goes below 0.5, the best 50% of
#'    each outcome. A `threshold` is never changed. The result is that of a
#'    search at the level reached, which `zone_used` reports. A widened zone
#'    may not improve every outcome, and each outcome's gain and its interval
#'    show whether it does.
#' 6. Centre. Each zone is summarized by its closed geometric mean (the
#'    inverse ILR of its mean coordinates), its per-part range and its
#'    best-ranked candidate.
#' 7. Gain. Each outcome is predicted at the overall zone centre (at its own
#'    centre when the zones do not overlap) and contrasted with `from`, as one
#'    [predict_realloc()] contrast. For `type = "relative"`, this gain is also
#'    the value reported at the centre.
#'
#' A predicted outcome linear in the ILR coordinates has no interior optimum,
#' and its zone lies on the edge of the region. A fit with no curvature term
#' (squares, products, `poly()`, `ns()`, `bs()`) is still searched, with one
#' message stating that the fit is linear in the composition. The zone of a
#' curved fit also lies on the boundary when its stationary point is outside
#' the region.
#'
#' A zone's `on_boundary` value is the share of its candidates in the outer
#' shell of the region. On the ellipsoid, that shell holds the candidates whose
#' squared Mahalanobis distance is above 0.9 of the limit; on the box, those
#' with a part within one `spacing` of a limit.
#'
#' The confidence intervals of the prediction at the centre and of the gain
#' are conditional on the centre, and the location of the zone itself has no
#' interval.
#'
#' The box, the ellipsoid and the grid do not depend on the order of the parts
#' when `spacing` divides the total, as the default does.
#'
#' A covariate that interacts with the composition must be pinned with `at` or
#' split with `by`; otherwise the call aborts and names it. A factor moderator
#' is split by default, as in [predict_realloc()]. With `by`, each level has
#' its own region, estimated from its own rows, and its own zone. On a fit
#' built with `include_total = TRUE`, every prediction holds `t_total` at the
#' fit's `total_ref`.
#'
#' `type = "absolute"` is refused for the engines without an absolute
#' prediction (`"cox"`, `"coxme"`, `"clm"`, `"clmm"` and `"mblogit"`; see
#' [comp_iv()]). A multinomial (`"mblogit"`) fit is searched on the log-odds of
#' the outcome level named by `category` against the reference level, so each
#' `category` has its own ranking and zone.
#'
#' @param object A [comp_iv()] fit, or a list of them, one per outcome, on the
#'   same parts (in the same order) and total. A list is named by outcome; an
#'   unnamed one takes each fit's outcome variable, the left-hand side of its
#'   model formula. The region is estimated from the first fit's composition.
#' @param direction `"max"` (default) when a higher predicted outcome is
#'   better, `"min"` when a lower one is. Give one value, or one per outcome.
#'   On a Cox fit (`"cox"`, `"coxme"`), `"max"` favours the highest hazard,
#'   and a search for the lowest risk uses `"min"`. On an ordinal fit
#'   (`"clm"`, `"clmm"`), the cumulative logit
#'   \eqn{P(Y \le j) = \mathrm{logit}^{-1}(\theta_j - \eta)}{P(Y <= j) = logit^-1(theta_j - eta)}
#'   moves probability to the higher outcome categories as the latent
#'   predictor \eqn{\eta}{eta} rises. `"max"` therefore favours the higher
#'   categories. On a multinomial fit, `"max"` favours the `category` level
#'   over the reference level.
#' @param zone The rank a candidate needs to be in an outcome's zone, a number
#'   in (0, 1); the default `0.95` keeps the best 5%.
#' @param region The feasible region. `"box"` (default) sets per-part limits
#'   (Dumuid et al. 2021; Brakenridge et al. 2024). `"ellipsoid"` is the
#'   chi-square ellipsoid in ILR coordinates of Mellow et al. (2026).
#'   `"observed"` takes the sample's own compositions (Lund Rasmussen et al.
#'   2023; Vindiš et al. 2025). A data frame gives the compositions directly,
#'   one column per part on the closure scale; every row must be positive and
#'   sum to the fit's total.
#' @param sd For `"box"`: a positive number of standard deviations, default
#'   `3`. Each part is limited to its mean \eqn{\pm}{+/-} `sd` standard
#'   deviations on the closure scale, never beyond the part's observed range.
#'   Passing `sd` with another region, or with `coverage`, is an error, and
#'   `sd = NULL` needs `coverage`.
#' @param coverage For `"box"`, `NULL` (default) uses `sd`. Otherwise the box
#'   takes percentile limits in place of `sd`: either the quantiles
#'   `c(lo, hi)` or one share for the central quantiles
#'   `c((1 - coverage) / 2, (1 + coverage) / 2)`. For `"ellipsoid"`, the
#'   chi-square level of its boundary, a share in (0, 1), default `0.80`; the
#'   printed result also gives the share of the sample inside. `coverage` is
#'   ignored, with a message, for the other regions.
#' @param spacing For `"box"` and `"ellipsoid"`: the grid step, a positive
#'   number of closure units. `NULL` (default) uses 1/144 of the closure
#'   declared in [declare_composition()], 10 minutes of a 1,440-minute day. A
#'   finer grid holds more candidates, and the search stops when more than
#'   1,000,000 lie inside the region, or more than 5,000,000 on the grid
#'   before the region is applied. The number of candidates grows steeply with
#'   the number of parts. With five or more parts the default search takes
#'   about a minute, and a coarser grid (e.g., 1/72 of the closure) shortens
#'   it.
#' @param resolution For [plot()], the steps along each edge of a panel's
#'   simplex, as in [plot.comp_iv()] (default 40).
#' @param threshold `NULL` (default), or a named list with one entry per
#'   outcome ruled by a value instead of `zone`. An entry is one number (a
#'   cut-off, met at or above it for `"max"` and at or below it for `"min"`)
#'   or two (an inclusive band; e.g., `list(zbmi = c(-1, 1))`). The value is
#'   on the `type` scale. For `"relative"`, it is the effect versus `from`: a
#'   ratio, and so positive, on a ratio scale such as an odds or hazard ratio,
#'   and otherwise a difference. For `"absolute"`, it is the predicted
#'   response.
#' @param widen `FALSE` (default) keeps every zone at `zone`. `TRUE` lowers
#'   `zone` in steps of 0.05 for the outcomes it rules, when their zones share
#'   no candidate, until they do. It never goes below 0.5 (the Intersection
#'   item of Details). With `by`, each level is widened on its own.
#' @param type `"relative"` (default) values each candidate by its effect
#'   versus `from` on the fit's effect scale, the measure [predict_realloc()]
#'   reports (e.g., an odds ratio, a hazard ratio or a mean difference).
#'   `"absolute"` values it by the predicted response, which the survival,
#'   ordinal and multinomial engines do not provide. The zone ruled by `zone`
#'   is the same under both. `threshold`, the zone limits, the `estimate` at
#'   the centre and [plot()] follow the scale.
#' @param category For a multinomial (`"mblogit"`) fit, the non-reference
#'   outcome level whose log-odds against the reference level is searched:
#'   required for such a fit, and refused when no fit is multinomial. `NULL`
#'   (default).
#' @param at Optional named list of covariate values held in every
#'   prediction, as in [predict_realloc()].
#' @param by Optional name of a factor covariate of every fit. The search then
#'   has one region and one zone per level, with the level held in its
#'   predictions. Defaults to a factor moderator not pinned by `at`, as in
#'   [predict_realloc()], and every fit must resolve to the same `by`.
#' @param from The composition the gain is measured from, and with it every
#'   value of a `type = "relative"` search. `NULL` (default) uses each fit's
#'   reference composition, stored by [comp_iv()]: the closed geometric mean
#'   of the composition the fit was given. Every level of `by` is measured from
#'   it. Otherwise, give a named numeric vector with one positive value per
#'   part, re-closed to the closure declared in [declare_composition()].
#'   When the zone centre equals `from`, the gain is the reference value (0,
#'   or 1 on a ratio scale) with no test.
#' @param conf_level Confidence level of the intervals, default `0.95`.
#'
#' @return An `optimal_comp` object, a list with elements:
#'   \describe{
#'     \item{`candidates`}{A tibble with one row per candidate inside the
#'       region, holding the parts on the closure scale, the ILR coordinates
#'       (`ilr1`, `ilr2`, ...), `mahalanobis` (the distance for the
#'       ellipsoid, else `NA`) and the `by` column. Per outcome follow
#'       `estimate_<outcome>`, `rank_<outcome>` and `in_zone_<outcome>`, then
#'       the overall `in_zone`. `estimate_<outcome>` is on the link scale: the
#'       prediction, or for `type = "relative"` its difference from the
#'       prediction at `from`.}
#'     \item{`zones`}{A tibble with one row per outcome and one `"overall"`
#'       row (per `by` level). Its columns are `rule` (`"quantile"` or
#'       `"threshold"`), `threshold_low`, `threshold_high`, `n`, per part
#'       `centre_<part>`, `min_<part>`, `max_<part>` and `best_<part>`, then
#'       `on_boundary`, and `q_overlap` on the overall row. The two thresholds
#'       are the limits the zone's values meet, on the `type` scale (the effect
#'       versus `from`, or the predicted response). `best_<part>` is the
#'       best-ranked candidate; for the overall zone, it is the candidate
#'       whose lowest rank across outcomes is highest.}
#'     \item{`optimum`}{A tibble with one row per outcome (per `by` level).
#'       Its columns are the parts of the centre evaluated, `centre_of`
#'       (`"overall"`, or `"outcome"` when the zones do not overlap), and
#'       `estimate`, `conf.low` and `conf.high` at that centre on the `type`
#'       scale. Under `"relative"` the estimate is the effect versus `from`
#'       and equals the gain; under `"absolute"` it is the predicted response.
#'       The gain from `from` follows on the fit's effect scale (`gain`,
#'       `gain.low`, `gain.high`, `p.value`), then the `from` parts as
#'       `from_<part>`.}
#'     \item{`fits`}{The fits, named by outcome, which [plot()] draws.}
#'     \item{Settings}{`direction`, `zone` (as given), `threshold`, `widen`,
#'       `type`, `category`, `region`, the `coverage`, `sd` and `spacing` used
#'       (`NULL` where they do not apply), `n_candidates` (per level) and
#'       `zone_used` (per level, the `zone` searched, which is `zone` or the
#'       step `widen` reached). It also holds `parts`, `labels`, `units`,
#'       `total`, `outcomes`, `linear` (per outcome, `TRUE` for a fit with no
#'       curvature term), `measures` (the gain's effect label per outcome),
#'       `by`, `at` and `conf_level`. `region` holds its `name` and, per
#'       level, its estimated `parameters`. For the ellipsoid these include
#'       the limit `r2` on the squared distance and `share`, the share of the
#'       level's sample inside. For the box they are the per-part limits `lo`
#'       and `hi` and the flags `lo_clipped` and `hi_clipped`, `TRUE` where a
#'       limit was cut to the part's observed minimum or maximum; [print()]
#'       reports those.}
#'   }
#'
#' @references
#' Brakenridge, C. J., Koster, A., de Galan, B. E., Carver, A., Dumuid, D.,
#'   Dzakpasu, F. Q. S., Eussen, S. J. P. M., Savelberg, H. H. C. M., Bosma,
#'   H., Owen, N., Schaper, N. C., Healy, G. N., & Dunstan, D. W. (2024).
#'   Associations of 24 h time-use compositions of sitting, standing, physical
#'   activity and sleeping with optimal cardiometabolic risk and glycaemic
#'   control: The Maastricht Study. *Diabetologia*, 67(7), 1356-1367.
#'   \doi{10.1007/s00125-024-06145-0}
#'
#' Dumuid, D., Wake, M., Burgner, D., Tremblay, M. S., Okely, A. D., Edwards,
#'   B., Dwyer, T., & Olds, T. (2021). Balancing time use for children's
#'   fitness and adiposity: Evidence to inform 24-hour guidelines for sleep,
#'   sedentary time and physical activity. *PLoS ONE*, 16(1), e0245501.
#'   \doi{10.1371/journal.pone.0245501}
#'
#' Duncan, M. J., Kuzik, N., Silva, D. A. S., Bélanger, R. E., Carson, V.,
#'   Chaput, J.-P., Faulkner, G., Ferro, M. A., Turcotte-Tremblay, A.-M.,
#'   Leatherdale, S. T., Patte, K. A., & Tremblay, M. S. (2024). "Goldilocks
#'   days" for adolescent mental health: Movement behaviour combinations for
#'   well-being, anxiety and depression by gender. *Mental Health and Physical
#'   Activity*, 26, 100572. \doi{10.1016/j.mhpa.2023.100572}
#'
#' Lund Rasmussen, C., Gába, A., Stanford, T., Dygrýn, J., Dumuid, D., Janda,
#'   D., & Hron, K. (2023). The Goldilocks Day for healthy adiposity measures
#'   among children and adolescents. *Frontiers in Public Health*, 11,
#'   1158634. \doi{10.3389/fpubh.2023.1158634}
#'
#' Mellow, M. L., Stanford, T. E., Olds, T., Miatke, A., Smith, A. E., &
#'   Dumuid, D. (2026). An interactive tool to personalise 24-hour activity,
#'   sitting and sleep prescription for optimal health outcomes. *npj Digital
#'   Medicine*. \doi{10.1038/s41746-026-02542-4}
#'
#' Vindiš, J., Pelclová, J., Jašková, P., & Hron, K. (2025). The Goldilocks
#'   days represent optimal time-use to prevent obesity, low physical
#'   performance, risk and fear of falling in older adults. *Scientific
#'   Reports*. \doi{10.1038/s41598-025-09645-0}
#'
#' @examples
#' \donttest{
#' cd <- declare_composition(timeuse, c("MVPA", "LPA", "SB", "SPT"),
#'                           closure = 1440, units = "min/day")
#' quad <- comp_iv(cd, formula = zbmi ~ ilr()^2 + I(ilr()^2) + sex,
#'                 model = "lm")
#' # The per-part box of the mean +/- 3 standard deviations, on a 10-minute grid
#' opt <- optimal_composition(quad, direction = "min")
#' opt
#'
#' # Two outcomes, their zones intersected
#' risk <- comp_iv(cd, "obesity", model = "glm", family = "binomial")
#' optimal_composition(list(zbmi = quad, obesity = risk), direction = "min")
#'
#' # The chi-square ellipsoid at the 80% level
#' optimal_composition(quad, direction = "min", region = "ellipsoid")
#'
#' # A box from the 1st to the 99th percentile of each part, on a 20-minute grid
#' optimal_composition(quad, direction = "min", coverage = c(0.01, 0.99),
#'                     spacing = 20)
#'
#' # A Cox fit, searched for the lowest hazard ratio versus the reference day
#' cox <- comp_iv(cd, formula = survival::Surv(time, event) ~ ilr()^2 +
#'                  I(ilr()^2) + sex, model = "cox")
#' optimal_composition(list(onset = cox), direction = "min")
#' }
#'
#' @seealso [comp_iv()] to fit the models this consumes; [predict_realloc()]
#'   for a directed reallocation between named parts.
#' @family reallocation and prediction
#' @export
optimal_composition <- function(object, direction = "max", zone = 0.95,
                                region = "box", sd = 3, coverage = NULL,
                                spacing = NULL, threshold = NULL, widen = FALSE,
                                type = c("relative", "absolute"),
                                category = NULL, at = NULL, by = NULL,
                                from = NULL, conf_level = 0.95) {
  env <- rlang::current_env()
  # `sd` counts as given only when passed: its default serves the box alone.
  sd_given <- !missing(sd) && !is.null(sd)
  type <- match.arg(type)
  # purrr wraps an abort inside a walk in "In index: i"; the checks over
  # the fits rethrow the abort itself, so its message stays the user's.
  walk_fits <- function(.x, .f) {
    tryCatch(purrr::iwalk(.x, .f),
             purrr_error_indexed = function(e) rlang::cnd_signal(e$parent))
  }
  fits <- if (inherits(object, "timecoda_fit") || !is.list(object)) list(object) else object
  if (length(fits) == 0) {
    cli::cli_abort("{.arg object} must be a {.cls comp_iv} fit or a list of them.")
  }
  walk_fits(fits, function(fit, k) {
    tc_check_fit(fit, call = env)
    meta <- tc_model_meta(fit$model_type)
    if (type == "absolute" && !meta$absolute_ok) {
      cli::cli_abort(c(
        "{.val {fit$model_type}} models support {.code type = \"relative\"} only.",
        "i" = "Their {.pkg emmeans} predictions have no single absolute-response scale ({meta$absolute_note})."
      ), call = env)
    }
  })
  # An unnamed fit takes its outcome variable, the model formula's left-hand side.
  outcomes <- names(fits)
  if (is.null(outcomes)) outcomes <- rep("", length(fits))
  unnamed <- is.na(outcomes) | outcomes == ""
  outcomes[unnamed] <- purrr::map_chr(fits[unnamed], function(fit) {
    deparse1(rlang::f_lhs(stats::formula(fit$model)))
  })
  if (anyDuplicated(outcomes)) {
    cli::cli_abort(c(
      "Two fits share the outcome name {.val {outcomes[duplicated(outcomes)]}}.",
      "i" = "Name the list of fits, one name per outcome."
    ))
  }
  names(fits) <- outcomes
  # A multinomial fit is searched on one outcome level's log-odds against the
  # reference level; `category` names it, and every other fit ignores it.
  per_category <- tc_check_category(fits, category, "search", call = env)
  categories <- purrr::map(per_category, function(p) if (p) category)
  report <- purrr::map(fits, function(fit) {
    tc_report_meta(fit$model_type, fit$family)
  })
  first <- fits[[1]]
  parts <- first$parts
  total <- first$total
  d <- length(parts)
  coords <- tc_coord_names(d)
  walk_fits(fits[-1], function(fit, k) {
    if (!identical(fit$parts, parts) || !isTRUE(all.equal(fit$total, total))) {
      cli::cli_abort(c(
        "Every fit must share the first fit's parts, in the same order, and its total.",
        "x" = "{.val {k}} has parts {.val {fit$parts}} and total {fit$total}; {.val {outcomes[1]}} has {.val {parts}} and {total}."
      ), call = env)
    }
  })
  comp_obs <- as.matrix(tc_comp_df(first$comp)[parts])
  same_comp <- purrr::map_lgl(fits, function(fit) {
    other <- as.matrix(tc_comp_df(fit$comp)[parts])
    identical(dim(other), dim(comp_obs)) && all(other == comp_obs)
  })
  if (!all(same_comp)) {
    cli::cli_inform(
      "The fits' compositions differ; the region is estimated from the first fit's ({.val {outcomes[1]}})."
    )
  }

  if (!is.character(direction) || !all(direction %in% c("max", "min")) ||
      !length(direction) %in% c(1, length(fits))) {
    cli::cli_abort(
      "{.arg direction} must be {.val max} or {.val min}: one value, or one per outcome ({length(fits)})."
    )
  }
  direction <- rlang::set_names(rep_len(direction, length(fits)), outcomes)
  if (!is.numeric(zone) || length(zone) != 1 || is.na(zone) ||
      zone <= 0 || zone >= 1) {
    cli::cli_abort("{.arg zone} must be a single number in (0, 1).")
  }
  if (!rlang::is_bool(widen)) {
    cli::cli_abort("{.arg widen} must be {.code TRUE} or {.code FALSE}.")
  }
  if (!is.null(threshold)) {
    ok <- is.list(threshold) && !is.null(names(threshold)) &&
      all(names(threshold) %in% outcomes) && !anyDuplicated(names(threshold)) &&
      all(purrr::map_lgl(threshold, function(t) {
        is.numeric(t) && length(t) %in% 1:2 && !anyNA(t) &&
          (length(t) == 1 || t[1] <= t[2])
      }))
    if (!ok) {
      cli::cli_abort(c(
        "{.arg threshold} must be a named list with one entry per outcome it rules.",
        "i" = "Each entry is one cut-off or two increasing band limits on the {.arg type} scale (the effect versus {.arg from} for {.val relative}, the predicted response for {.val absolute}); the outcomes are {.val {outcomes}}."
      ))
    }
    # A relative threshold on a ratio scale is a ratio versus `from`.
    ratio <- type == "relative" &
      purrr::map_lgl(report[names(threshold)], function(r) isTRUE(r$exponentiate))
    bad <- names(threshold)[ratio &
                              purrr::map_lgl(threshold, function(t) any(t <= 0))]
    if (length(bad) > 0) {
      cli::cli_abort(
        "With {.code type = \"relative\"}, {.arg threshold} for {.val {bad}} is a ratio versus {.arg from} and must be positive."
      )
    }
  }
  if (!is.numeric(conf_level) || length(conf_level) != 1 ||
      is.na(conf_level) || conf_level <= 0 || conf_level >= 1) {
    cli::cli_abort("{.arg conf_level} must be a single number in (0, 1).")
  }

  # The region: a registry row by name, or the "given" row for a data frame,
  # whose rows become the candidates.
  given <- NULL
  if (is.data.frame(region)) {
    missing_parts <- setdiff(parts, names(region))
    if (length(missing_parts) > 0) {
      cli::cli_abort(
        "A data-frame {.arg region} needs a column for every part; missing {.val {missing_parts}}."
      )
    }
    check_composition(region, parts, total = total, tol = 1e-6)
    if (nrow(region) == 0) {
      cli::cli_abort("A data-frame {.arg region} needs at least one row.")
    }
    given <- as.matrix(region[parts]) / rowSums(region[parts]) * total
    region <- "given"
  } else if (!rlang::is_string(region) ||
             !region %in% c("ellipsoid", "box", "observed")) {
    cli::cli_abort(
      "{.arg region} must be {.val ellipsoid}, {.val box}, {.val observed} or a data frame of compositions."
    )
  }
  row <- tc_optimal_regions[[region]]
  if (sd_given && region != "box") {
    cli::cli_abort("{.arg sd} applies to {.code region = \"box\"} only.")
  }
  if (region == "ellipsoid") {
    if (is.null(coverage)) coverage <- 0.80
    if (!is.numeric(coverage) || length(coverage) != 1 || is.na(coverage) ||
        coverage <= 0 || coverage >= 1) {
      cli::cli_abort("{.arg coverage} must be a single share in (0, 1) for the ellipsoid.")
    }
    sd <- NULL
  } else if (region == "box" && !is.null(coverage)) {
    # A percentile box: `coverage` replaces the default `sd`.
    if (sd_given) {
      cli::cli_abort("Give {.arg coverage} or {.arg sd} for the box, not both.")
    }
    if (is.numeric(coverage) && length(coverage) == 1 && !is.na(coverage) &&
        coverage > 0 && coverage < 1) {
      coverage <- c((1 - coverage) / 2, (1 + coverage) / 2)
    }
    if (!is.numeric(coverage) || length(coverage) != 2 || anyNA(coverage) ||
        coverage[1] < 0 || coverage[1] >= coverage[2] || coverage[2] > 1) {
      cli::cli_abort(
        "{.arg coverage} must be one share in (0, 1) or two increasing quantiles in [0, 1] for the box."
      )
    }
    sd <- NULL
  } else if (region == "box") {
    if (is.null(sd)) {
      cli::cli_abort("The box needs {.arg sd} or {.arg coverage}; give one of them.")
    }
    if (!is.numeric(sd) || length(sd) != 1 || is.na(sd) || sd <= 0) {
      cli::cli_abort("{.arg sd} must be a single positive number.")
    }
  }
  if (!row$grid) {
    if (!is.null(coverage)) cli::cli_inform("{.arg coverage} is ignored for this region.")
    coverage <- NULL
    sd <- NULL
    if (!is.null(spacing)) {
      cli::cli_abort(
        "{.arg spacing} applies to the {.val ellipsoid} and {.val box} regions only."
      )
    }
  } else {
    # Ten minutes of a 1,440-minute day, for any closure.
    if (is.null(spacing)) spacing <- total / 144
    if (!is.numeric(spacing) || length(spacing) != 1 || is.na(spacing) ||
        spacing <= 0 || spacing * d > total) {
      cli::cli_abort(
        "{.arg spacing} must be a positive number of closure units, at most {signif(total / d, 4)} for {d} parts summing to {total}."
      )
    }
  }

  # Covariates: `at` as predict_realloc() takes it, `by` one zone per level.
  if (!is.null(at) &&
      (!is.list(at) || is.null(names(at)) || any(names(at) == ""))) {
    cli::cli_abort("{.arg at} must be a named list of covariate values.")
  }
  # `by` defaults to a factor moderator, as in predict_realloc(); every fit
  # must resolve it alike.
  by_first <- tc_resolve_by(first, by, at, call = env)
  walk_fits(fits[-1], function(fit, k) {
    by_k <- tc_resolve_by(fit, by, at, call = env)
    if (!identical(by_k, by_first)) {
      shown <- purrr::map_chr(list(by_first, by_k), function(b) {
        if (is.null(b)) "none" else b
      })
      cli::cli_abort(c(
        "The fits resolve {.arg by} differently.",
        "x" = "{.val {outcomes[1]}} gives {.val {shown[1]}}; {.val {k}} gives {.val {shown[2]}}.",
        "i" = "Give {.arg by}, or pin the moderator with {.arg at}."
      ), call = env)
    }
  })
  by <- by_first
  if (!is.null(by) && nrow(first$data) != nrow(comp_obs)) {
    cli::cli_abort(
      "The fit's data and its composition have different row counts, so {.arg by} cannot split the region."
    )
  }
  walk_fits(fits, function(fit, k) {
    unpinned <- setdiff(tc_interacting_covariates(fit), c(names(at), by))
    if (length(unpinned) > 0) {
      cli::cli_abort(c(
        "The {.val {k}} model includes composition-by-covariate interactions, so its optimal composition is covariate-conditional.",
        "i" = "Pin {.val {unpinned}} with {.arg at}, or find one zone per level with {.arg by}."
      ), call = env)
    }
  })
  # The gain is measured from `from`, by default each fit's own reference.
  if (is.null(from)) {
    from <- purrr::map(fits, function(fit) unlist(fit$ref_comp)[parts])
  } else {
    if (!is.numeric(from) || is.null(names(from)) ||
        !setequal(names(from), parts) || length(from) != d ||
        anyNA(from) || any(from <= 0)) {
      cli::cli_abort(
        "{.arg from} must be a named numeric vector with one positive value per part ({.val {parts}})."
      )
    }
    from <- purrr::map(fits, function(fit) from[parts] / sum(from) * total)
  }

  # A fit without a curvature term is linear in the coordinates: each of its
  # terms in them holds one coordinate, entered bare.
  linear <- purrr::map_lgl(fits, function(fit) {
    terms <- tc_composition_terms(c(coords, fit$covariates), coords)
    all(purrr::map_lgl(c(terms$only, terms$crossed), function(term) {
      cmp <- strsplit(term, "(?<!:):(?!:)", perl = TRUE)[[1]]
      has_coord <- purrr::map_lgl(cmp, function(x) {
        any(all.vars(str2lang(x)) %in% coords)
      })
      sum(has_coord) == 1 && all(cmp[has_coord] %in% coords)
    }))
  })
  if (any(linear)) {
    cli::cli_inform(c(
      "{.val {outcomes[linear]}} {cli::qty(sum(linear))}{?is/are} linear in the composition, so {?its/their} best predicted outcome lies on the edge of the region.",
      "i" = "The zone then marks where the region ends; see {.field on_boundary}.",
      "i" = "A curved surface can hold an interior optimum: refit with {.code curvature = TRUE} (or {.code ilr()^2 + I(ilr()^2)} in a formula)."
    ))
  }

  strata <- if (is.null(by)) {
    list(NULL)
  } else {
    col <- first$data[[by]]
    as.list(if (is.factor(col)) levels(droplevels(col)) else sort(unique(stats::na.omit(col))))
  }
  sign <- ifelse(direction == "max", 1, -1)
  # One region, candidate set and zone per level of `by`. lapply(), not
  # purrr::map(), so an abort inside keeps its own message.
  results <- lapply(strata, function(lev) {
    rows <- if (is.null(lev)) seq_len(nrow(comp_obs)) else which(first$data[[by]] == lev)
    pins <- c(at, if (!is.null(lev)) rlang::set_names(list(lev), by))
    cand <- tc_optimal_candidates(row, comp_obs[rows, , drop = FALSE], given,
                                  coverage, sd, spacing, total, first$pivot,
                                  call = env)
    n <- nrow(cand$comps)
    preds <- purrr::imap(fits, function(fit, k) {
      tc_optimal_predict(fit, cand, pins, type, from[[k]], categories[[k]])
    })
    # A candidate emmeans cannot estimate (NA) is ranked among none and is in
    # no zone.
    purrr::iwalk(preds, function(p, k) {
      n_na <- sum(is.na(p$link))
      if (n_na > 0) {
        cli::cli_warn(
          "{n_na} candidate{?s} could not be predicted for {.val {k}} and {cli::qty(n_na)}{?is/are} left out of its zone."
        )
      }
    })
    ranks <- purrr::map2(preds, sign, function(p, s) {
      rank(s * p$link, na.last = "keep") / sum(!is.na(p$link))
    })
    # Each outcome's zone at the rank `level`; a threshold ignores the level.
    zone_at <- function(level) {
      purrr::imap(preds, function(p, k) {
        t <- threshold[[k]]
        member <- if (is.null(t)) {
          ranks[[k]] >= level
        } else if (length(t) == 2) {
          p$response >= t[1] & p$response <= t[2]
        } else if (direction[[k]] == "max") {
          p$response >= t
        } else {
          p$response <= t
        }
        !is.na(member) & member
      })
    }
    in_zone <- zone_at(zone)
    overall <- purrr::reduce(in_zone, `&`)
    # The largest `zone` at which the quantile-ruled zones would share a
    # candidate meeting every threshold.
    quantile_ruled <- setdiff(outcomes, names(threshold))
    feasible <- purrr::reduce(in_zone[names(threshold)], `&`, .init = rep(TRUE, n))
    shared <- if (length(quantile_ruled) > 0) {
      stats::na.omit(do.call(pmin, unname(ranks[quantile_ruled]))[feasible])
    }
    q_overlap <- if (length(shared) > 0) max(shared) else NA_real_
    # `widen`: an empty overall zone takes the first 5-point step below `zone`,
    # no lower than 0.5, at which the zones share a candidate, i.e. the first
    # step at or below `q_overlap`. Rounding keeps each step the typed number
    # (0.6, not 0.6000000000000001); the search then runs at that level.
    level <- zone
    steps <- numeric()
    if (widen && !any(overall) && !is.na(q_overlap)) {
      n_steps <- max(0, floor((zone - 0.5) / 0.05 + 1e-9))
      steps <- round(zone - 0.05 * seq_len(n_steps), 12)
      reached <- steps[steps <= q_overlap][1]
      if (!is.na(reached)) {
        level <- reached
        in_zone <- zone_at(level)
        overall <- purrr::reduce(in_zone, `&`)
      }
    }
    summarise_zone <- function(member, best) {
      comps <- cand$comps[member, , drop = FALSE]
      out <- list(n = sum(member))
      summ <- if (out$n > 0) {
        list(centre = unlist(geom_mean_comp(as.data.frame(comps), parts, total = total)),
             min = apply(comps, 2, min), max = apply(comps, 2, max),
             best = comps[which.max(best[member]), ])
      } else {
        purrr::map(c(centre = 1, min = 1, max = 1, best = 1), function(i) rep(NA_real_, d))
      }
      for (s in names(summ)) {
        out[paste0(s, "_", parts)] <- as.list(unname(summ[[s]]))
      }
      out$on_boundary <- if (out$n > 0) mean(cand$boundary[member]) else NA_real_
      tibble::as_tibble(out)
    }
    zones <- purrr::imap(in_zone, function(member, k) {
      t <- threshold[[k]]
      resp <- preds[[k]]$response[member]
      limits <- if (length(t) == 2) {
        t
      } else if (length(t) == 1) {
        if (direction[[k]] == "max") c(t, NA) else c(NA, t)
      } else if (!any(member)) {
        c(NA, NA)
      } else if (direction[[k]] == "max") {
        c(min(resp), NA)
      } else {
        c(NA, max(resp))
      }
      dplyr::bind_cols(
        tibble::tibble(outcome = k,
                       rule = if (is.null(t)) "quantile" else "threshold",
                       threshold_low = as.numeric(limits[1]),
                       threshold_high = as.numeric(limits[2])),
        summarise_zone(member, ranks[[k]]),
        tibble::tibble(q_overlap = NA_real_)
      )
    }) |>
      dplyr::bind_rows()
    zones <- dplyr::bind_rows(zones, dplyr::bind_cols(
      tibble::tibble(outcome = "overall", rule = NA_character_,
                     threshold_low = NA_real_, threshold_high = NA_real_),
      summarise_zone(overall, do.call(pmin, unname(ranks))),
      tibble::tibble(q_overlap = q_overlap)
    ))
    at_level <- if (!is.null(lev)) glue::glue(" for {by} = {lev}") else ""
    if (level < zone) {
      best <- format(100 * (1 - level))
      cli::cli_inform(if (is.null(threshold)) {
        "The outcomes' zones{at_level} share no candidate at {.code zone = {zone}}; they were widened to {.code zone = {level}}, the best {best}% of each outcome's predictions."
      } else {
        "The outcomes' zones{at_level} share no candidate at {.code zone = {zone}}; they were widened to {.code zone = {level}}, the best {best}% of the predictions of each outcome ruled by {.arg zone}."
      })
    }
    if (!any(overall)) {
      if (is.na(q_overlap)) {
        cli::cli_inform(
          "The outcomes' zones{at_level} share no candidate, and no candidate meets every threshold."
        )
      } else {
        floor_level <- if (length(steps) > 0) min(steps)
        # Rounded down, so a search at the printed level overlaps.
        meet <- format(floor(q_overlap * 1000) / 1000)
        cli::cli_inform(c(
          "The outcomes' zones{at_level} share no candidate; they would meet at {.code zone = {meet}}.",
          "i" = if (!is.null(floor_level)) {
            "Widening in 5-point steps to {.code zone = {floor_level}}, the best {format(100 * (1 - floor_level))}%, found no overlap."
          }
        ))
      }
    }

    # The prediction at the centre and the gain from `from`.
    centre_of <- if (any(overall)) "overall" else "outcome"
    optimum <- purrr::imap(fits, function(fit, k) {
      z <- zones[zones$outcome == (if (any(overall)) "overall" else k), ]
      centre <- rlang::set_names(unlist(z[paste0("centre_", parts)]), parts)
      values <- rep(NA_real_, 7)
      if (!anyNA(centre)) {
        emm <- tc_emm_at_points(fit, rbind(from[[k]], centre), at = pins,
                                conf_level = conf_level,
                                category = categories[[k]])
        # A centre equal to `from` is the reference against itself: the
        # reference value and no test, never sent to emmeans' contrast.
        gain <- if (isTRUE(all.equal(unname(from[[k]]), unname(centre),
                                     tolerance = 1e-8))) {
          ref_value <- if (isTRUE(report[[k]]$exponentiate)) 1 else 0
          c(ref_value, ref_value, ref_value, NA_real_)
        } else {
          g <- summary(emmeans::contrast(emm, method = "trt.vs.ctrl", ref = 1),
                       infer = c(TRUE, TRUE), level = conf_level,
                       type = "response", adjust = "none")
          c(tc_estimate_col(g), tc_ci_col(g, "lower"), tc_ci_col(g, "upper"),
            g$p.value)
        }
        # A relative search reports the centre's value versus `from`, the gain.
        values <- if (type == "relative") {
          c(gain[1:3], gain)
        } else {
          s <- summary(emm, infer = c(TRUE, FALSE), level = conf_level,
                       type = "response")
          c(tc_estimate_col(s)[2], tc_ci_col(s, "lower")[2],
            tc_ci_col(s, "upper")[2], gain)
        }
      }
      dplyr::bind_cols(
        tibble::tibble(outcome = k),
        tibble::as_tibble(as.list(centre)),
        tibble::tibble(centre_of = centre_of, estimate = values[1],
                       conf.low = values[2], conf.high = values[3],
                       gain = values[4], gain.low = values[5],
                       gain.high = values[6], p.value = values[7]),
        tibble::as_tibble(as.list(rlang::set_names(from[[k]], paste0("from_", parts))))
      )
    }) |>
      dplyr::bind_rows()

    candidates <- dplyr::bind_cols(
      tibble::as_tibble(as.data.frame(cand$comps)),
      tibble::as_tibble(as.data.frame(cand$z)),
      tibble::tibble(mahalanobis = cand$mahalanobis)
    )
    for (k in outcomes) {
      candidates[[paste0("estimate_", k)]] <- preds[[k]]$estimate
      candidates[[paste0("rank_", k)]] <- ranks[[k]]
      candidates[[paste0("in_zone_", k)]] <- in_zone[[k]]
    }
    candidates$in_zone <- overall
    if (!is.null(lev)) {
      candidates[[by]] <- lev
      zones[[by]] <- lev
      optimum[[by]] <- lev
      candidates <- dplyr::relocate(candidates, dplyr::all_of(by), .after = "mahalanobis")
      zones <- dplyr::relocate(zones, dplyr::all_of(by), .after = "outcome")
      optimum <- dplyr::relocate(optimum, dplyr::all_of(by), .after = "outcome")
    }
    list(candidates = candidates, zones = zones, optimum = optimum,
         parameters = cand$parameters, n = n, zone_used = level)
  })
  level_names <- if (is.null(by)) "all" else as.character(unlist(strata))

  structure(
    list(
      candidates = dplyr::bind_rows(purrr::map(results, "candidates")),
      zones = dplyr::bind_rows(purrr::map(results, "zones")),
      optimum = dplyr::bind_rows(purrr::map(results, "optimum")),
      fits = fits,
      direction = direction,
      zone = zone,
      threshold = threshold,
      widen = widen,
      type = type,
      category = category,
      region = list(
        name = region,
        parameters = rlang::set_names(purrr::map(results, "parameters"), level_names)
      ),
      coverage = coverage,
      sd = sd,
      spacing = spacing,
      n_candidates = rlang::set_names(purrr::map_int(results, "n"), level_names),
      zone_used = rlang::set_names(purrr::map_dbl(results, "zone_used"), level_names),
      parts = parts,
      labels = first$comp$labels,
      units = first$comp$units,
      total = total,
      outcomes = outcomes,
      linear = linear,
      measures = purrr::map_chr(report, "measure"),
      by = by,
      at = at,
      conf_level = conf_level
    ),
    class = "optimal_comp"
  )
}

# Feasible regions ------------------------------------------------------------

# One row per `region`. `grid` marks a region searched on a `spacing` grid.
# `estimate` reads its parameters off a level's compositions (closure scale)
# and ILR coordinates; `range` bounds each part of the grid, given `V`, the
# D x (D - 1) matrix with clr(x) = V z in the fit's pivot basis; `distance`,
# `inside` and `boundary` read candidates; `label` names the region in print()
# from every level's parameters.
tc_optimal_regions <- list(
  ellipsoid = list(
    grid = TRUE,
    # The chi-square fence of Mellow et al. (2026); `share` is the part of the
    # level's sample inside it.
    estimate = function(comps, z, coverage, sd) {
      mu <- colMeans(z)
      s <- stats::cov(z)
      r2 <- stats::qchisq(coverage, df = ncol(z))
      list(mu = mu, S = s, coverage = coverage, r2 = r2,
           share = mean(stats::mahalanobis(z, mu, s) <= r2))
    },
    # Conservative per-part bounds: with a = V[j, ] - V[k, ], the largest
    # log(x_j / x_k) on the ellipsoid is u_jk = a' mu + sqrt(r2 a' S a), so
    # x_j <= total / sum_k exp(-u_jk) and x_j >= total / sum_k exp(u_kj), the
    # k = j term being 1. Base-R matrix algebra.
    range = function(par, V, total) {
      d <- nrow(V)
      pairs <- which(diag(d) == 0, arr.ind = TRUE)
      a <- V[pairs[, 1], , drop = FALSE] - V[pairs[, 2], , drop = FALSE]
      u <- matrix(0, d, d)
      u[pairs] <- drop(a %*% par$mu) + sqrt(par$r2 * rowSums((a %*% par$S) * a))
      cbind(total / colSums(exp(u)), total / rowSums(exp(-u)))
    },
    distance = function(z, par) stats::mahalanobis(z, par$mu, par$S),
    inside = function(comps, z, par) {
      stats::mahalanobis(z, par$mu, par$S) <= par$r2
    },
    boundary = function(comps, z, par, spacing) {
      stats::mahalanobis(z, par$mu, par$S) > 0.9 * par$r2
    },
    label = function(pars, sd) {
      share <- round(100 * purrr::map_dbl(pars, "share"), 1)
      pct <- function(v) format(v, nsmall = 1)
      held <- if (length(share) == 1) {
        glue::glue("{pct(share)}% of the sample")
      } else if (min(share) == max(share)) {
        # Levels whose shares round alike are given one share.
        glue::glue("{pct(share[1])}% of each level's sample")
      } else {
        glue::glue("{pct(min(share))}% to {pct(max(share))}% of each level's sample")
      }
      glue::glue("ILR ellipsoid at the {100 * pars[[1]]$coverage}% chi-square level, holding {held}")
    }
  ),
  box = list(
    grid = TRUE,
    estimate = function(comps, z, coverage, sd) {
      if (is.null(sd)) {
        lo <- apply(comps, 2, stats::quantile, probs = coverage[1], names = FALSE)
        hi <- apply(comps, 2, stats::quantile, probs = coverage[2], names = FALSE)
      } else {
        lo <- colMeans(comps) - sd * apply(comps, 2, stats::sd)
        hi <- colMeans(comps) + sd * apply(comps, 2, stats::sd)
      }
      # The limits never leave the observed range of each part; the flags
      # record which limits were cut to it.
      mins <- apply(comps, 2, min)
      maxs <- apply(comps, 2, max)
      list(lo = pmax(lo, mins), hi = pmin(hi, maxs), coverage = coverage,
           lo_clipped = lo < mins, hi_clipped = hi > maxs)
    },
    range = function(par, V, total) cbind(par$lo, par$hi),
    distance = function(z, par) rep(NA_real_, nrow(z)),
    inside = function(comps, z, par) {
      rowSums(sweep(comps, 2, par$lo, ">=") & sweep(comps, 2, par$hi, "<=")) ==
        ncol(comps)
    },
    # A candidate is on the boundary when a part is within one spacing of a
    # limit.
    boundary = function(comps, z, par, spacing) {
      near <- sweep(comps, 2, par$lo + spacing, "<") |
        sweep(comps, 2, par$hi - spacing, ">")
      rowSums(near) > 0
    },
    label = function(pars, sd) {
      # An ordinal percentile: 1st, 2nd, 3rd, 11th, 95th, 2.5th.
      nth <- function(p) {
        n <- 100 * p
        last <- if (n == round(n) && !round(n) %% 100 %in% 11:13) round(n) %% 10 else 0
        paste0(n, if (last %in% 1:3) c("st", "nd", "rd")[last] else "th")
      }
      if (is.null(sd)) {
        cov <- pars[[1]]$coverage
        glue::glue("per-part box from the {nth(cov[1])} to the {nth(cov[2])} percentile")
      } else {
        glue::glue("per-part box of the mean +/- {sd} SD")
      }
    }
  ),
  observed = list(
    grid = FALSE,
    estimate = function(comps, z, coverage, sd) list(n = nrow(comps)),
    label = function(pars, sd) "the sample's own compositions"
  ),
  given = list(
    grid = FALSE,
    estimate = function(comps, z, coverage, sd) list(),
    label = function(pars, sd) "the compositions given as the region"
  )
)

# One level's candidates inside the region: the `spacing` grid or the given
# rows, with their distances and boundary flags.
tc_optimal_candidates <- function(row, comp_obs, given, coverage, sd, spacing,
                                  total, pivot, call = rlang::caller_env()) {
  parts <- colnames(comp_obs)
  d <- length(parts)
  if (nrow(comp_obs) <= d) {
    cli::cli_abort(
      "The region needs more than {d} compositions to estimate; this level has {nrow(comp_obs)}.",
      call = call
    )
  }
  z_obs <- tc_pivot_coords(comp_obs, pivot = pivot)
  par <- row$estimate(comp_obs, z_obs, coverage, sd)
  out <- list(parameters = par)
  if (!row$grid) {
    comps <- if (is.null(given)) comp_obs else given
    out$comps <- comps
    out$z <- tc_pivot_coords(comps, pivot = pivot)
    out$mahalanobis <- rep(NA_real_, nrow(comps))
    out$boundary <- rep(NA, nrow(comps))
    return(out)
  }
  # V from the centred log-ratios of the unit coordinate vectors.
  unit <- log(as.matrix(ilr_backtransform(diag(d - 1), parts, pivot = pivot)))
  lims <- row$range(par, t(unit - rowMeans(unit)), total)
  # The first D - 1 parts on multiples of `spacing` within their limits, the
  # last part the remainder; the tolerance keeps a value at a limit.
  tol <- 1e-9 * total
  values <- purrr::map(seq_len(d - 1), function(j) {
    v <- seq(spacing, total - (d - 1) * spacing, by = spacing)
    v[v >= lims[j, 1] - tol & v <= lims[j, 2] + tol]
  })
  n_raw <- prod(lengths(values))
  if (n_raw > 5e6) {
    cli::cli_abort(c(
      "The {.arg spacing} grid would hold {format(n_raw, big.mark = ',')} compositions before the region filters it, above the limit of 5,000,000.",
      "i" = "Use a coarser {.arg spacing}."
    ), call = call)
  }
  nodes <- as.matrix(expand.grid(values, KEEP.OUT.ATTRS = FALSE))
  comps <- cbind(nodes, total - rowSums(nodes))
  last <- comps[, d]
  comps <- comps[last >= max(spacing, lims[d, 1]) - tol & last <= lims[d, 2] + tol,
                 , drop = FALSE]
  colnames(comps) <- parts
  z <- if (nrow(comps) > 0) tc_pivot_coords(comps, pivot = pivot)
  keep <- if (nrow(comps) > 0) which(row$inside(comps, z, par)) else integer(0)
  if (length(keep) > 1e6) {
    cli::cli_abort(c(
      "The {.arg spacing} grid holds {format(length(keep), big.mark = ',')} compositions inside the region, above the limit of 1,000,000.",
      "i" = "Use a coarser {.arg spacing}."
    ), call = call)
  }
  if (length(keep) == 0) {
    cli::cli_abort(c(
      "No composition of the {.arg spacing} grid lies inside the region.",
      "i" = "Use a finer {.arg spacing}, or widen the region with {.arg sd} or {.arg coverage}."
    ), call = call)
  }
  out$comps <- comps[keep, , drop = FALSE]
  out$z <- z[keep, , drop = FALSE]
  out$mahalanobis <- row$distance(out$z, par)
  out$boundary <- row$boundary(out$comps, out$z, par, spacing)
  out
}

# Methods ---------------------------------------------------------------------

#' @describeIn optimal_composition Prints a text report in sections. The report
#'   gives the search's settings, each part's region limits, each zone's centre
#'   and per-part range, and the gain from `from` to the zone centre. The limits
#'   are the box's, or the candidates' range for another region, and a mark
#'   shows a box limit that is the part's observed minimum or maximum. With
#'   `by`, the limits, zones and gain are given for each level. A widened search
#'   states how far its zones were widened, and `type = "absolute"` adds the
#'   predicted outcome at each zone centre.
#' @param x An `optimal_comp` object.
#' @param ... Unused by [print()], present for S3 consistency. [plot()] takes
#'   no further arguments, and one given there is an error.
#' @export
print.optimal_comp <- function(x, ...) {
  labs <- tc_part_labels(x$parts, x$labels)
  lab_col <- format(labs)
  digits <- max(1, 3 - floor(log10(x$total)))
  num <- function(v) as.character(signif(v, 3))
  count <- function(n) format(n, big.mark = ",")
  # Right-align a column of numbers given as text.
  right <- function(v) formatC(v, width = max(nchar(v)))
  # Wrap `text` to the console width: the first line follows `first`, the
  # others are indented as deep as `first` is long.
  wrap <- function(text, first = "  ") {
    lines <- strwrap(text, width = getOption("width") - nchar(first) + 1)
    lead <- c(first, rep(strrep(" ", nchar(first)), length(lines) - 1))
    glue::glue("{lead}{lines}")
  }
  # Each row is compared with its own fit's reference day.
  refs <- purrr::map(x$optimum$outcome, function(k) {
    unname(unlist(x$fits[[k]]$ref_comp)[x$parts])
  })
  from_ref <- all(purrr::map2_lgl(seq_len(nrow(x$optimum)), refs, function(i, ref) {
    isTRUE(all.equal(unname(unlist(x$optimum[i, paste0("from_", x$parts)])), ref))
  }))
  same_ref <- all(purrr::map_lgl(refs, function(ref) isTRUE(all.equal(ref, refs[[1]]))))
  from_name <- if (!from_ref) {
    "the given day"
  } else if (same_ref) {
    "the reference day"
  } else {
    "each fit's reference day"
  }
  relative <- x$type == "relative"
  # The value a threshold rules: the predicted response, or the effect versus
  # `from` on the fit's effect scale.
  ruled <- function(k) {
    if (relative) glue::glue("{x$measures[[k]]} versus {from_name}") else "predicted response"
  }
  rules <- purrr::map_chr(x$outcomes, function(k) {
    t <- x$threshold[[k]]
    better <- if (x$direction[[k]] == "max") "higher" else "lower"
    if (is.null(t)) {
      glue::glue("best {format(100 * (1 - x$zone))}% of predictions, {better} is better")
    } else if (length(t) == 2) {
      glue::glue("{ruled(k)} in [{num(t[1])}, {num(t[2])}]")
    } else {
      glue::glue("{ruled(k)} {if (better == 'higher') '>=' else '<='} {num(t)}")
    }
  })
  # One rule for every outcome is given once, without their names.
  rules <- if (length(unique(rules)) == 1) rules[1] else glue::glue("{x$outcomes}: {rules}")
  scale <- if (!relative) {
    "predicted response"
  } else if (length(x$outcomes) == 1) {
    glue::glue("{x$measures[[1]]} versus {from_name}")
  } else {
    per_outcome <- glue::glue_collapse(glue::glue("{x$measures} ({x$outcomes})"), sep = ", ")
    glue::glue("{per_outcome}, versus {from_name}")
  }
  # With `by`, each level's candidate count heads its own zones.
  n_all <- if (is.null(x$by)) count(x$n_candidates[[1]])
  grid <- if (!is.null(x$spacing)) {
    n_grid <- if (!is.null(n_all)) glue::glue("; {n_all} candidates") else ""
    glue::glue("multiples of {tc_with_units(signif(x$spacing, 4), x$units)}{n_grid}")
  }
  settings <- purrr::compact(list(
    Parts = glue::glue(
      "{glue::glue_collapse(labs, sep = ', ')} (total = {tc_with_units(x$total, x$units)})"
    ),
    Region = tc_optimal_regions[[x$region$name]]$label(x$region$parameters, x$sd),
    Grid = grid,
    Candidates = if (is.null(grid)) n_all,
    Scale = scale,
    `Zone rule` = rules
  ))
  key_width <- max(nchar(names(settings)))
  search <- c("Search", unlist(purrr::imap(settings, function(values, key) {
    first <- glue::glue("  {formatC(key, width = -key_width)}   ")
    pad <- strrep(" ", nchar(first))
    purrr::map2(values, seq_along(values), function(v, i) {
      wrap(v, if (i == 1) first else pad)
    })
  })))
  # Each part's limits per level: the box's stored limits, with a mark where
  # one was cut to the part's observed minimum or maximum, else the
  # candidates' range.
  box <- x$region$name == "box"
  units <- if (!is.null(x$units)) glue::glue(" ({x$units})") else ""
  limits <- purrr::map(names(x$n_candidates), function(lev) {
    if (box) {
      p <- x$region$parameters[[lev]]
      lo <- p$lo[x$parts]
      hi <- p$hi[x$parts]
      lo_cut <- unname(p$lo_clipped[x$parts])
      hi_cut <- unname(p$hi_clipped[x$parts])
    } else {
      cand <- x$candidates
      if (!is.null(x$by)) cand <- cand[as.character(cand[[x$by]]) == lev, ]
      lo <- purrr::map_dbl(x$parts, function(q) min(cand[[q]]))
      hi <- purrr::map_dbl(x$parts, function(q) max(cand[[q]]))
      lo_cut <- hi_cut <- rep(FALSE, length(x$parts))
    }
    marked <- any(lo_cut) || any(hi_cut)
    # A blank in place of the mark keeps every number on one column.
    mark <- function(cut) if (marked) ifelse(cut, "*", " ") else ""
    lower <- glue::glue("{tc_fmt(lo, digits)}{mark(lo_cut)}")
    upper <- glue::glue("{tc_fmt(hi, digits)}{mark(hi_cut)}")
    head <- glue::glue("{c('lower', 'upper')}{mark(FALSE)}")
    width <- max(nchar(c(lower, upper, head)))
    col <- function(v) formatC(v, width = width)
    kinds <- c(if (any(lo_cut)) "minimum", if (any(hi_cut)) "maximum")
    title <- if (box) "Region limits" else "Candidate range"
    at_level <- if (!is.null(x$by)) glue::glue(", {x$by} = {lev}") else ""
    c(
      glue::glue("{title}{at_level}{units}"),
      glue::glue("  {strrep(' ', nchar(lab_col[1]))}   {col(head[1])}   {col(head[2])}"),
      glue::glue("  {lab_col}   {col(lower)}   {col(upper)}"),
      if (marked) {
        glue::glue("  * the part's observed {glue::glue_collapse(kinds, sep = ' or ')}")
      }
    )
  })
  strata <- if (is.null(x$by)) list(NULL) else as.list(names(x$n_candidates))
  results <- purrr::map(strata, function(lev) {
    zones <- x$zones
    optimum <- x$optimum
    if (!is.null(lev)) {
      zones <- zones[zones[[x$by]] == lev, ]
      optimum <- optimum[optimum[[x$by]] == lev, ]
    }
    n_cand <- x$n_candidates[[if (is.null(lev)) 1 else lev]]
    # A widened level says so once, above the zones it widened.
    used <- x$zone_used[[if (is.null(lev)) 1 else lev]]
    best <- format(100 * (1 - used))
    widened <- if (used < x$zone) {
      whose <- if (is.null(x$threshold)) {
        "each outcome's predictions"
      } else {
        "the predictions of each outcome ruled by zone"
      }
      wrap(glue::glue(
        "No candidate is in every zone at zone = {x$zone}; the zones are widened to zone = {used}, the best {best}% of {whose}"
      ))
    }
    level <- if (!is.null(lev)) {
      c(glue::glue("{x$by} = {lev}: {count(n_cand)} candidates"), widened)
    } else if (!is.null(widened)) {
      c("Zones", widened)
    }
    blocks <- purrr::map(c(x$outcomes, if (length(x$outcomes) > 1) "overall"), function(k) {
      z <- zones[zones$outcome == k, ]
      head_line <- if (k == "overall" && used < x$zone) {
        glue::glue("Overall zone, widened to the best {best}%")
      } else if (k == "overall") {
        "Overall zone"
      } else {
        glue::glue("Zone for {k}")
      }
      if (z$n == 0) {
        line <- if (k == "overall" && !is.na(z$q_overlap)) {
          # Rounded down, as in the search's message.
          meet <- format(floor(z$q_overlap * 1000) / 1000)
          glue::glue("no overlap; the zones meet at zone = {meet}")
        } else {
          "empty"
        }
        return(glue::glue("{head_line}: {line}"))
      }
      value <- function(s) right(tc_fmt(unlist(z[paste0(s, "_", x$parts)]), digits))
      notes <- character()
      # A relative search's value at the centre is its gain, printed below.
      if (k != "overall" && !relative) {
        o <- optimum[optimum$outcome == k, ]
        where <- if (o$centre_of == "outcome") {
          "its own centre"
        } else if (length(x$outcomes) > 1) {
          "the overall centre"
        } else {
          "the zone centre"
        }
        notes <- c(notes, wrap(glue::glue(
          "Predicted {k} at {where}: {num(o$estimate)} [{num(o$conf.low)}; {num(o$conf.high)}]"
        )))
      }
      linear <- if (k == "overall") any(x$linear) else x$linear[[k]]
      if (linear) {
        notes <- c(
          notes,
          wrap("A fit linear in the composition has no interior optimum: the zone is the edge of the region in the direction of benefit."),
          wrap("Refit with `curvature = TRUE` for a surface that can hold one.")
        )
      } else if (!is.na(z$on_boundary) && z$on_boundary > 0.5) {
        notes <- c(notes, wrap(glue::glue(
          "{round(100 * z$on_boundary)}% of the zone lies near the edge of the region (its outer shell); the best composition may lie beyond it."
        )))
      }
      c(
        glue::glue("{head_line}: {count(z$n)} compositions (centre [min; max])"),
        glue::glue("  {lab_col}   {value('centre')}  [{value('min')}; {value('max')}]"),
        notes
      )
    })
    to_name <- if (optimum$centre_of[1] == "overall") {
      "the zone centre"
    } else {
      "each outcome's own centre"
    }
    measure <- purrr::map_chr(optimum$outcome, function(k) x$measures[[k]])
    p <- format.pval(optimum$p.value, digits = 3)
    p <- ifelse(startsWith(p, "<"), glue::glue("p < {sub('^< *', '', p)}"), glue::glue("p = {p}"))
    gain <- c(
      glue::glue("Gain from {from_name} to {to_name} ({100 * x$conf_level}% CI)"),
      glue::glue(
        "  {format(optimum$outcome)}   {format(measure)}   {right(num(optimum$gain))}  ",
        "{format(glue::glue('[{num(optimum$gain.low)}; {num(optimum$gain.high)}]'))}  {p}"
      )
    )
    c(list(level), blocks, list(gain))
  })
  sections <- purrr::compact(c(
    list(glue::glue(
      "<optimal_comp> Goldilocks-zone search: {glue::glue_collapse(x$outcomes, sep = ', ')}"
    )),
    list(search),
    limits,
    purrr::list_flatten(results)
  ))
  # One blank line between sections, and no blank at the end of a line (the
  # space that stands for an absent mark in the last column).
  lines <- unlist(purrr::map(sections, function(s) c("", as.character(s))))[-1]
  print(glue::glue_collapse(sub(" +$", "", lines), sep = "\n"))
  invisible(x)
}

#' @describeIn optimal_composition Ternary panels (`choose(D, 3)`) of one
#'   outcome's surface on the search's `type` scale (the effect versus `from`,
#'   or the predicted response), drawn for every engine. An ordinal fit's
#'   surface is its cumulative odds ratio versus `from`, and a multinomial
#'   fit's is the relative risk ratio of the searched `category`. Each panel is
#'   a slice through the zone centre (the overall one, or the outcome's own
#'   when the zones do not overlap). The surface is drawn as [plot.comp_iv()]
#'   draws one, on the cells inside the search's region. For the box and the
#'   ellipsoid these are the cells the region itself holds (on the parameters
#'   of the `level` drawn); otherwise they are the cells within the
#'   candidates' per-part range. The surface holds the search's `at` values
#'   and, for a search with `by`, the `level` drawn. The convex hull of the
#'   computed cells in the outcome's zone is outlined in red, and the centre is
#'   marked. A panel with fewer than three computed cells in the zone has no
#'   outline, and a higher `resolution` computes more cells. A panel carries no
#'   caption: the zone centre, the held covariates and the region's limits are
#'   in [print()]. Returns the plot object (one plot, or a list named by
#'   triplet), drawn as a side effect. An argument the method does not take is
#'   an error.
#' @param outcome The outcome whose surface is drawn; `NULL` (default) for the
#'   first.
#' @param level For a search with `by`, the level whose zone is drawn; `NULL`
#'   (default) for the first. Not used without `by`.
#' @param smooth `NULL` (default) draws the computed cells as tiles for a fit
#'   linear in the composition and a quadratic surface through them for any
#'   other fit; `TRUE` or `FALSE` draws the smooth or the tiles for any fit,
#'   as in [plot.comp_iv()]. A panel whose smooth cannot be fitted is drawn as
#'   tiles, with a warning.
#' @export
plot.optimal_comp <- function(x, outcome = NULL, level = NULL,
                              resolution = 40L, smooth = NULL, ...) {
  if (!inherits(x, "optimal_comp")) {
    cli::cli_abort("{.arg x} must be an {.cls optimal_comp} object.")
  }
  tc_check_plot_dots(
    "plot.optimal_comp", ...,
    .hint = "The region and its band are set in {.fn optimal_composition} ({.arg region}, {.arg sd}, {.arg coverage}, {.arg spacing})."
  )
  if (!is.null(smooth) && !rlang::is_bool(smooth)) {
    cli::cli_abort("{.arg smooth} must be {.code NULL}, {.code TRUE} or {.code FALSE}.")
  }
  if (is.null(outcome)) outcome <- x$outcomes[1]
  if (!rlang::is_string(outcome) || !outcome %in% x$outcomes) {
    cli::cli_abort("{.arg outcome} must be one of {.val {x$outcomes}}.")
  }
  fit <- x$fits[[outcome]]
  # The Cox, ordinal and multinomial engines are searched with type =
  # "relative" only, and their surface is the effect versus `from`, the same
  # contrast the search ranked. A multinomial fit is drawn for the searched
  # `category`, as plot.comp_iv() draws one for its own `category`.
  category <- if (tc_model_meta(fit$model_type)$per_category) x$category
  levels_by <- names(x$n_candidates)
  zones <- x$zones
  candidates <- x$candidates
  optimum <- x$optimum[x$optimum$outcome == outcome, ]
  pins <- x$at
  if (is.null(x$by)) {
    if (!is.null(level)) {
      cli::cli_abort("{.arg level} applies to a search with {.arg by} only.")
    }
  } else {
    if (is.null(level)) level <- levels_by[1]
    if (!rlang::is_scalar_atomic(level) || !as.character(level) %in% levels_by) {
      cli::cli_abort("{.arg level} must be one of {.val {levels_by}}, the levels of {.val {x$by}}.")
    }
    level <- as.character(level)
    zones <- zones[as.character(zones[[x$by]]) == level, ]
    candidates <- candidates[as.character(candidates[[x$by]]) == level, ]
    optimum <- optimum[as.character(optimum[[x$by]]) == level, ]
    pins <- c(pins, rlang::set_names(list(level), x$by))
  }
  parts <- x$parts
  overall <- zones[zones$outcome == "overall", ]
  own <- zones[zones$outcome == outcome, ]
  zone <- if (overall$n > 0) overall else own
  if (zone$n == 0) {
    cli::cli_abort("The zone for {.val {outcome}} is empty; there is nothing to outline.")
  }
  centre <- rlang::set_names(unlist(zone[paste0("centre_", parts)]), parts)
  region <- tc_optimal_regions[[x$region$name]]
  par <- x$region$parameters[[if (is.null(x$by)) 1 else level]]
  # A grid region keeps the cells its own test accepts, on the level's stored
  # parameters; the observed days and a data frame keep the candidates' range.
  bounds <- if (region$grid) {
    function(comps) region$inside(comps, tc_pivot_coords(comps, fit$pivot), par)
  } else {
    purrr::map(rlang::set_names(parts), function(p) range(candidates[[p]]))
  }
  # A relative surface is the effect versus the search's `from`.
  scale <- x$type
  from_k <- rlang::set_names(unlist(optimum[paste0("from_", parts)]), parts)
  from_name <- if (isTRUE(all.equal(unname(from_k),
                                    unname(unlist(fit$ref_comp)[parts])))) {
    "the reference day"
  } else {
    "the given day"
  }
  report <- tc_report_meta(fit$model_type, fit$family)
  fill_name <- if (scale == "relative" && !is.null(category)) {
    glue::glue("{report$measure} of {category} versus {from_name}")
  } else if (scale == "relative") {
    glue::glue("{report$measure} versus {from_name}")
  } else {
    glue::glue("Predicted {outcome}")
  }
  # The smooth is fitted where the model is linear in the coordinates, as in
  # plot.comp_iv(): the log of a ratio, the family's link for a response.
  fam <- if (scale == "absolute" && inherits(fit$family, "family")) fit$family
  link <- inv <- identity
  if (scale == "relative" && isTRUE(report$exponentiate)) {
    link <- log
    inv <- exp
  } else if (!is.null(fam)) {
    link <- fam$linkfun
    inv <- fam$linkinv
  }
  # By default, a fit linear in the composition has a flat surface and its
  # cells are drawn as tiles; any other fit is smoothed.
  if (is.null(smooth)) smooth <- !x$linear[[outcome]]
  out <- tc_triplet_plots(parts, function(tri) {
    p <- tc_surface_panel(fit, tri, centre, resolution, scale, fill_name,
                          bounds, smooth, link, inv, at = pins,
                          control = if (scale == "relative") from_k,
                          category = category,
                          kept_by = if (region$grid) "region" else "candidates")
    # The cells already lie in the region; the outline keeps those in the zone.
    cells <- p$data
    keep <- !is.na(cells$.pred) &
      (is.na(own$threshold_low) | cells$.pred >= own$threshold_low) &
      (is.na(own$threshold_high) | cells$.pred <= own$threshold_high)
    cells <- cells[keep, , drop = FALSE]
    if (nrow(cells) >= 3) {
      hull <- grDevices::chull(cells$x, cells$y)
      p <- p + ggplot2::geom_path(data = cells[c(hull, hull[1]), ],
                                  colour = "red", linewidth = 0.8)
    } else {
      cli::cli_inform(
        "Fewer than three computed cells of {.val {tri}} lie in the zone; no outline is drawn."
      )
    }
    p +
      ggplot2::geom_point(
        data = tibble::tibble(x = centre[[tri[1]]], y = centre[[tri[2]]],
                              z = centre[[tri[3]]]),
        colour = "red", size = 2.5
      )
  })
  tc_draw_plots(out)
}
