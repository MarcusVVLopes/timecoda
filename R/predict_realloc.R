# Module 2: time reallocations via `emmeans` contrasts.
#
# A reallocation is a contrast between a reallocated composition and a
# reference composition on the link scale, not an absolute prediction. Each
# composition is evaluated as its own emmeans reference grid (averaged over
# covariates), the grids are combined with rbind(), and the reallocated rows
# are contrasted against the reference row; `type = "response"` then yields
# odds, hazard or rate ratios. What each directed `scheme` means is defined
# once, in `tc_realloc_schemes`, which plot_realloc() reads as well.

#' Predict time-reallocation effects
#'
#' Estimates the effect of reallocating time between the parts of a
#' composition, from a [comp_iv()] fit. Each estimate is the contrast of the
#' model's linear predictor between a reallocated composition and a reference
#' composition, computed by `emmeans`. It is back-transformed to the response
#' scale: an odds, hazard, rate or risk ratio for a ratio-scale model, and a
#' mean difference for a linear one.
#'
#' @details
#' The reference and every reallocated composition are converted to the fit's
#' ILR coordinates and evaluated as their own `emmeans` reference grids,
#' averaged over the covariates. Each reallocated row is contrasted with the
#' reference row on the link scale,
#' \eqn{\eta(x_{\mathrm{realloc}}) - \eta(x_{\mathrm{ref}})}{eta(x_realloc) - eta(x_ref)},
#' with a delta-method standard error, and back-transformed with
#' `type = "response"`. Every row is its own Wald test, with no multiplicity
#' adjustment across rows. The estimates do not depend on the ILR basis when
#' the model is linear in the coordinates or holds their full quadratic
#' ([ilr()]). A row with `delta = 0` contrasts the reference with itself and is
#' reported as the reference value (1 for a ratio, 0 for a difference), with a
#' zero standard error, `NA` degrees of freedom and no test.
#'
#' Without composition-by-covariate interactions, the covariates cancel in the
#' contrast and a relative estimate does not depend on them. With such
#' interactions the contrast is covariate-conditional. Every interacting
#' covariate must then be pinned with `at`, or the contrast evaluated per level
#' with `by`; the call aborts naming any covariate left unpinned. On a fit
#' built with `include_total = TRUE`, every composition is evaluated with
#' `t_total` held at the fit's `total_ref`. A reallocation then changes the
#' shares at a constant total, and the total cancels from a relative contrast.
#'
#' `type = "absolute"` returns the predicted response at each composition
#' instead of a contrast. It is not available for `"cox"`, `"coxme"`, `"clm"`,
#' `"clmm"` or `"mblogit"` fits. These have no single absolute-response scale:
#' Cox leaves its baseline hazard unestimated, and ordinal and multinomial
#' predictions are per outcome category.
#'
#' The `df` column holds the degrees of freedom `emmeans` computes for the
#' fit, and the per-part table of [summary()][summary.comp_iv] uses the same
#' degrees of freedom. For an `"lmer"` fit, the default is the Kenward-Roger
#' approximation through \pkg{pbkrtest}. The Satterthwaite approximation
#' through \pkg{lmerTest} is used when it is selected with
#' `emmeans::emm_options(lmer.df =)`, or when \pkg{pbkrtest} cannot be used
#' (absent, disabled, or more observations than
#' `emmeans::emm_options(pbkrtest.limit =)`). Both packages are optional, and
#' when neither can be used the degrees of freedom are `Inf` (asymptotic). For
#' `"glmer"`, `"cox"`, `"coxme"` and `"rlmer"` fits the degrees of freedom are
#' `Inf` (\pkg{robustlmm} reports none), and for `"lmrob"` fits they are the
#' residual degrees of freedom.
#'
#' A reallocation that drives any part to zero or below has no composition to
#' evaluate. By default (`on_invalid = "na"`), such a row is kept with `NA`
#' estimates. A sweep over many pairs or `delta`s then still returns every
#' feasible row, and the printed summary gives the number of rows with no
#' estimate. `"omit"` drops those rows instead, and `"error"` aborts the whole
#' call on the first one.
#'
#' @section The reference composition and "isotemporal substitution":
#' The `reference` argument sets the composition every reallocation starts
#' from, and the contrast is computed the same way at any reference. The
#' default, the cohort's geometric-mean day, suits population-typical
#' questions. A day supplied by the user, such as a guideline schedule or one
#' person's day, suits scenario questions. In the compositional
#' time-use literature this prediction is *compositional isotemporal
#' substitution* (Dumuid et al. 2019). The traditional leave-one-out
#' isotemporal substitution model (Mekary et al. 2009) is a different
#' regression parameterization, and no choice of `reference` reproduces it.
#' `vignette("predicting-reallocations")` compares the two.
#'
#' @param object A `comp_iv` / `timecoda_fit` from [comp_iv()].
#' @param reference Optional named numeric vector (or one-row data frame)
#'   giving the composition to anchor the reallocations at: one strictly
#'   positive value per part, summing to the fit's composition total. The
#'   default (`NULL`) uses the sample's closed geometric mean, stored in the
#'   fit.
#' @param delta Time to reallocate, in the composition's closure units: a
#'   number or vector of numbers, each applied to every `from`/`to` pair (a
#'   negative value moves time the other way). Ignored in custom mode (when
#'   `realloc` is supplied), where each scenario has its own signed changes.
#' @param realloc Optional custom multi-part reallocation, in the style of
#'   Dumuid et al. (2022) "Your best day". In place of a directed `from`/`to`
#'   pair, give a signed change per part, in the composition's closure units.
#'   The changes are added to the reference composition, which is then
#'   re-closed. For example, `c(SPT = 30, SB = -15, LPA = -15)` moves 30
#'   minutes into sleep, taken 15 each from `SB` and `LPA`. Accepts (1) a named
#'   numeric vector (one scenario), (2) a list of named numeric vectors
#'   (several scenarios), or (3) a data frame with one column per part and one
#'   row per scenario. Every scenario must name known parts, sum to zero (to
#'   within `1e-8`; a closed reallocation), and hold at least one positive and
#'   one negative entry. Parts left unnamed are held at their reference value.
#'   Custom mode is mutually exclusive with `from`, `to` and `scheme`. `type`,
#'   `at`, `by` and `on_invalid` apply as for a directed reallocation, and a
#'   multinomial fit still gives one row set per outcome category.
#' @param from,to Character vectors of part names defining the
#'   reallocations. For `scheme = "one-to-one"` (the default), they name the
#'   donor (`from`) and receiver (`to`) parts and default to all ordered pairs.
#'   If given, they are matched elementwise (scalars recycled). For
#'   `scheme = "one-to-all"`, supply only one of them. Each part it names is a
#'   focal part in turn, and the other side is implicitly all remaining parts
#'   (see `scheme`).
#' @param scheme How time is reallocated. `"one-to-one"` (default) moves
#'   `delta` closure units directly between the `from` and `to` parts.
#'   `"one-to-all"` moves `delta` closure units between a single focal part
#'   \eqn{f} and all the remaining parts, in proportion to their values in the
#'   reference composition \eqn{r}. Remaining part \eqn{j} changes by
#'   \eqn{\mp \delta \, r_j / \sum_{k \ne f} r_k}{-/+ delta * r_j / sum(r_k, k != f)},
#'   and the remaining parts keep their ratios to one another. With
#'   `to = "MVPA"` the focal part gains time from the rest, and with
#'   `from = "MVPA"` it gives time to the rest. With neither, every part is
#'   used as the focal receiver in turn. The implicit "all remaining parts"
#'   side is labelled `"(rest)"` in the `from`/`to` output column.
#' @param type `"relative"` (default) contrasts each reallocated composition
#'   against the reference; `"absolute"` returns the predicted response at
#'   each composition. `"cox"`, `"coxme"`, `"clm"`, `"clmm"` and `"mblogit"`
#'   fits support `"relative"` only.
#' @param at Optional named list of covariate values at which to evaluate.
#'   Each covariate that interacts with the composition must be pinned by `at`
#'   or split by `by`. Covariates not in `at` are averaged over in the usual
#'   `emmeans` way.
#' @param by Optional name of a factor or character covariate, to compute the
#'   reallocation per group. The contrast is evaluated once at each level
#'   present in the fitted data. The result gains a column named after the
#'   covariate, holding the level (e.g., one odds ratio per `sex` for a
#'   moderated effect). When the fit was built with a factor `moderator` (see
#'   [comp_iv()]) and `at` does not pin it, `by` defaults to that moderator.
#'   `at` and `by` cannot name the same covariate. Pin numeric covariates with
#'   `at` instead.
#' @param conf_level Confidence level for the intervals.
#' @param on_invalid What to do when a reallocation drives a part to zero or
#'   below: `"na"` (default) keeps the rows but returns `NA` estimates for
#'   them; `"omit"` drops those rows; `"error"` aborts.
#'
#' @return A `realloc` object, a tibble with columns `from`, `to`, `delta`,
#'   `estimate`, `std.error`, `df`, `conf.low`, `conf.high`, `statistic` and
#'   `p.value` (`statistic` and `p.value` only for `type = "relative"`).
#'   `std.error` is on the scale of `estimate`. For a ratio effect it is the
#'   standard error of the ratio itself, as `emmeans` reports it on the
#'   response scale; the standard error of the log ratio is not returned. For
#'   `type = "absolute"`, the first row is the reference composition (`from`,
#'   `to` = `NA`, `delta` = 0). A `by` grouping adds a column named after the
#'   covariate, and a multinomial fit adds `response_level`, one row set per
#'   non-reference outcome category. In custom mode, each scenario's `from` and
#'   `to` are compact, comma-separated signed part strings (e.g.,
#'   `from = "LPA -15, SB -15"`, `to = "SPT +30"`), and `delta` is the
#'   scenario's total time moved, in closure units. [summary()][summary.realloc]
#'   builds the summary table from it. The fit metadata that method needs is
#'   stored as attributes, among them `units`, the composition's unit label
#'   (`NULL` when none was declared).
#'
#' @references
#' Dumuid, D., Pedišić, Ž., Stanford, T. E., et al. (2019). The compositional
#'   isotemporal substitution model: A method for estimating changes in a
#'   health outcome for reallocation of time between sleep, physical activity
#'   and sedentary behaviour. *Statistical Methods in Medical Research*,
#'   28(3), 846-857. \doi{10.1177/0962280217737805}
#'
#' Mekary, R. A., Willett, W. C., Hu, F. B., & Ding, E. L. (2009).
#'   Isotemporal substitution paradigm for physical activity epidemiology and
#'   weight change. *American Journal of Epidemiology*, 170(4), 519-527.
#'   \doi{10.1093/aje/kwp163}
#'
#' Dumuid, D., Olds, T., Wake, M., et al. (2022). "Your best day": An
#'   interactive app to translate how time reallocations within a 24-hour day
#'   are associated with health measures. *PLOS ONE*, 17(9), e0272343.
#'   \doi{10.1371/journal.pone.0272343}
#'
#' Kenward, M. G., & Roger, J. H. (1997). Small sample inference for fixed
#'   effects from restricted maximum likelihood. *Biometrics*, 53(3), 983-997.
#'   \doi{10.2307/2533558}
#'
#' @examples
#' cd <- declare_composition(timeuse, c("MVPA", "LPA", "SB", "SPT"),
#'                           closure = 1440)
#' fit <- comp_iv(cd, outcome = "obesity", model = "glm", family = "binomial")
#' predict_realloc(fit, delta = 30, from = "SB", to = "MVPA")
#'
#' # Proportional reallocation: add 30 min to MVPA, taken from the rest of the
#' # day in proportion to its composition.
#' predict_realloc(fit, delta = 30, scheme = "one-to-all", to = "MVPA")
#'
#' # Custom multi-part reallocation ("best-day" style): add 30 min to SPT,
#' # taken 15 each from SB and LPA.
#' predict_realloc(fit, realloc = c(SPT = 30, SB = -15, LPA = -15))
#'
#' # Anchor the same transfer at a chosen day instead of the cohort mean.
#' predict_realloc(fit, delta = 30, from = "SB", to = "MVPA",
#'                 reference = c(MVPA = 30, LPA = 240, SB = 660, SPT = 510))
#'
#' @seealso [comp_iv()] to fit the model this consumes;
#'   [summary.realloc()] for the summary table; [plot_realloc()] for
#'   the continuous sweep figure; [optimal_composition()] to search for
#'   the best-predicted composition instead of a directed transfer.
#' @family reallocation and prediction
#' @export
predict_realloc <- function(object, delta, from = NULL, to = NULL,
                            scheme = c("one-to-one", "one-to-all"),
                            realloc = NULL, reference = NULL,
                            type = c("relative", "absolute"), at = NULL,
                            by = NULL, conf_level = 0.95,
                            on_invalid = c("na", "omit", "error")) {
  tc_check_fit(object)
  scheme_given <- !missing(scheme)
  scheme <- match.arg(scheme)
  type <- match.arg(type)
  on_invalid <- match.arg(on_invalid)
  if (!is.numeric(conf_level) || length(conf_level) != 1 ||
      is.na(conf_level) || conf_level <= 0 || conf_level >= 1) {
    cli::cli_abort("{.arg conf_level} must be a single number in (0, 1).")
  }
  reference <- tc_check_reference(reference, object)
  by <- tc_resolve_by(object, by, at)

  custom <- NULL
  if (!is.null(realloc)) {
    if (!is.null(from) || !is.null(to) || scheme_given) {
      cli::cli_abort(c(
        "{.arg realloc} (custom reallocation) is mutually exclusive with {.arg from}, {.arg to}, and {.arg scheme}.",
        "i" = "Supply either a directed {.arg from}/{.arg to} (with {.arg scheme}) or a custom {.arg realloc} vector, not both."
      ))
    }
    custom <- tc_realloc_normalise(realloc, object$parts)
    # `delta` is unused in custom mode; give it a concrete value rather than
    # leaving an unevaluated missing promise to travel with the call.
    delta <- NA_real_
  }
  if (is.null(custom) && missing(delta)) {
    cli::cli_abort(c(
      "{.arg delta} is missing.",
      "i" = "Supply {.arg delta} (in the composition's units) with {.arg from}/{.arg to}, or a custom {.arg realloc} vector."
    ))
  }

  # Run the engine once, or once per level of `by`, tagging each block with
  # the level. Its aborts name this call.
  env <- rlang::current_env()
  run <- function(at) {
    tc_realloc_engine(
      object, reference = reference, delta = delta, from = from, to = to,
      scheme = scheme, custom = custom, type = type, at = at,
      conf_level = conf_level, on_invalid = on_invalid, call = env
    )
  }
  out <- if (is.null(by)) {
    run(at)
  } else {
    col <- object$data[[by]]
    # Only the levels present in the data: a declared but unused level has no
    # rows in the model frame, and emmeans refuses it as a new level.
    levs <- if (is.factor(col)) levels(droplevels(col)) else sort(unique(stats::na.omit(col)))
    # lapply(), not purrr::map(), so an abort inside keeps its own message.
    stacked <- lapply(levs, function(lev) {
      res <- run(c(at, stats::setNames(list(lev), by)))
      res[[by]] <- lev
      res
    }) |>
      dplyr::bind_rows()
    # Surface the grouping column just after from/to for readability.
    front <- intersect(c("from", "to", "delta", by), names(stacked))
    stacked[c(front, setdiff(names(stacked), front))]
  }
  # The metadata summary.realloc() needs travels with the result, so the
  # summary does not require the fit.
  attr(out, "effect_label") <- tc_effect_label(object, type)
  attr(out, "realloc_type") <- type
  attr(out, "scheme") <- if (is.null(custom)) scheme else "custom"
  attr(out, "conf_level") <- conf_level
  attr(out, "parts") <- object$parts
  attr(out, "total") <- object$total
  attr(out, "labels") <- object$comp$labels
  attr(out, "units") <- object$comp$units
  attr(out, "by") <- by
  # A tibble subclass, so the result prints and subsets as a tibble.
  class(out) <- c("realloc", class(out))
  out
}

# Validate a user-supplied reference composition against the fit's parts and
# total; NULL falls back to the fit's stored (geometric-mean) reference.
tc_check_reference <- function(reference, object,
                               call = rlang::caller_env()) {
  if (is.null(reference)) {
    return(unlist(object$ref_comp))
  }
  if (is.data.frame(reference)) {
    if (nrow(reference) != 1) {
      cli::cli_abort("{.arg reference} must be a single composition (one row).",
                     call = call)
    }
    reference <- unlist(reference[intersect(object$parts, names(reference))])
  }
  if (!is.numeric(reference) || anyNA(reference) || is.null(names(reference)) ||
      !setequal(names(reference), object$parts)) {
    cli::cli_abort(
      "{.arg reference} must be a named numeric vector with one value per part ({.val {object$parts}}).",
      call = call
    )
  }
  reference <- reference[object$parts]
  if (any(reference <= 0)) {
    cli::cli_abort("All parts of {.arg reference} must be strictly positive.",
                   call = call)
  }
  if (abs(sum(reference) - object$total) > 1e-6 * object$total) {
    cli::cli_abort(
      "{.arg reference} sums to {round(sum(reference), 4)}, but the fit's composition total is {format(object$total)}.",
      call = call
    )
  }
  reference
}

# Resolve the `by` group-wise argument. Auto-defaults to the fit's `moderator`
# (when it is a factor/character and `at` does not pin it), so a moderator
# fit yields per-level estimates by default. An explicit `by` must
# name a factor/character covariate that `at` does not also pin; pin numeric
# covariates with `at` instead.
tc_resolve_by <- function(object, by, at, call = rlang::caller_env()) {
  if (is.null(by)) {
    mod <- object$moderator
    if (!is.null(mod) && !mod %in% names(at)) {
      col <- object$data[[mod]]
      if (is.factor(col) || is.character(col)) {
        return(mod)
      }
    }
    return(NULL)
  }
  if (!is.character(by) || length(by) != 1) {
    cli::cli_abort("{.arg by} must be a single covariate column name.", call = call)
  }
  # The fit's data also carries the outcome and any survival or weight
  # columns, none of which the grid can vary, so `by` must name a variable of a
  # covariate term (`sex`, or inside `ilr1:sex` or `factor(sex)`).
  term_vars <- unlist(purrr::map(object$covariates, function(t) {
    all.vars(str2lang(t))
  }))
  if (!by %in% c(term_vars, object$moderator)) {
    cli::cli_abort("{.arg by} {.val {by}} is not a covariate in the fitted model.",
                   call = call)
  }
  col <- object$data[[by]]
  if (!is.factor(col) && !is.character(col)) {
    cli::cli_abort(
      c("{.arg by} must name a factor or character covariate.",
        "i" = "For a numeric covariate, pin specific values with {.arg at}."),
      call = call
    )
  }
  if (by %in% names(at)) {
    cli::cli_abort(
      c("{.arg at} and {.arg by} both name {.val {by}}.",
        "i" = "Pin one level with {.arg at}, or estimate every level with {.arg by}, not both."),
      call = call
    )
  }
  by
}

# Internal engine --------------------------------------------------------------

# Validate the request against the fit, build the grid of reallocated
# compositions and contrast it: once, or once per outcome category of a
# multinomial fit (stacked with a `response_level` column).
tc_realloc_engine <- function(object, reference, delta, from, to, scheme,
                              custom = NULL, type,
                              at, conf_level, on_invalid,
                              call = rlang::caller_env()) {
  meta <- tc_model_meta(object$model_type)
  if (!meta$absolute_ok && type == "absolute") {
    cli::cli_abort(
      c(
        "{.val {object$model_type}} models support {.code type = \"relative\"} only.",
        "i" = "Their {.pkg emmeans} predictions have no single absolute-response scale ({meta$absolute_note})."
      ),
      call = call
    )
  }
  if (object$has_ilr_interaction) {
    # Every covariate that interacts with the composition must be pinned; a
    # `by` level arrives here already merged into `at`.
    unpinned <- setdiff(tc_interacting_covariates(object), names(at))
    if (length(unpinned) > 0) {
      int_terms <- tc_composition_terms(
        object$covariates, tc_coord_names(length(object$parts))
      )$crossed
      # The example value has the covariate's own type: a logical, the first
      # level (also of a numeric column the formula wraps in factor()), or the
      # mean.
      example <- purrr::map(rlang::set_names(unpinned), function(v) {
        col <- object$data[[v]]
        as_level <- any(grepl(paste0("factor(", v, ")"), int_terms, fixed = TRUE) |
                          grepl(paste0("ordered(", v, ")"), int_terms, fixed = TRUE))
        if (is.logical(col)) {
          sort(unique(col))[1]
        } else if (!is.numeric(col)) {
          levels(droplevels(as.factor(col)))[1]
        } else if (as_level) {
          as.numeric(min(col, na.rm = TRUE))
        } else {
          signif(mean(col, na.rm = TRUE), 3)
        }
      })
      example <- deparse1(rlang::call2("list", !!!example))
      cli::cli_abort(
        c(
          "This model includes composition-by-covariate interactions, so reallocation effects are covariate-conditional.",
          "i" = "Supply {cli::qty(length(unpinned))}{?a value/values} for {.val {unpinned}} via {.arg at} (e.g. {.code at = {example}})."
        ),
        call = call
      )
    }
  }
  if (is.null(custom) &&
      (!is.numeric(delta) || length(delta) < 1 || anyNA(delta) ||
       any(!is.finite(delta)))) {
    cli::cli_abort("{.arg delta} must be a numeric vector of finite values (in the composition's units).", call = call)
  }

  grid <- if (!is.null(custom)) {
    tc_realloc_grid_custom(reference, custom, on_invalid, call = call)
  } else {
    moves <- tc_realloc_schemes[[scheme]]$moves(object$parts, from, to,
                                                call = call)
    tc_realloc_grid(reference, moves, delta, scheme, on_invalid, call = call)
  }

  # Drop invalid rows up front when requested; otherwise they stay in the grid
  # and (for on_invalid = "na") receive NA results.
  if (on_invalid == "omit") {
    keep <- which(grid$valid)
    if (length(keep) == 0) {
      cli::cli_abort(
        "Every reallocation drives a part to zero or below; nothing to predict.",
        call = call
      )
    }
    grid$from <- grid$from[keep]
    grid$to <- grid$to[keep]
    grid$delta <- grid$delta[keep]
    grid$comps <- grid$comps[keep, , drop = FALSE]
    grid$valid <- rep(TRUE, length(keep))
  }

  if (!meta$per_category) {
    return(tc_realloc_result(object, reference, grid, at, conf_level, type))
  }
  # The non-reference outcome categories, read off the coefficient names
  # ("<level>~<term>").
  # lapply(), not purrr::map(), so an abort inside keeps its own message.
  levs <- unique(sub("~.*$", "", names(stats::coef(object$model))))
  out <- lapply(levs, function(lev) {
    res <- tc_realloc_result(object, reference, grid, at, conf_level,
                             type, category = lev)
    res$response_level <- lev
    res
  }) |>
    dplyr::bind_rows()
  front <- intersect(c("from", "to", "delta", "response_level"), names(out))
  out[c(front, setdiff(names(out), front))]
}

# The emmeans contrast + tidy tibble for one reallocation grid. For multinomial
# fits `category` names the non-reference outcome level whose per-category
# linear predictor is used; otherwise the model's single linear predictor.
tc_realloc_result <- function(object, reference, grid, at, conf_level,
                              type, category = NULL) {
  relative <- type == "relative"
  # A relative zero-delta row is the reference against itself: it is filled
  # below with the reference value and no test, never sent to emmeans.
  zero <- relative & grid$valid & grid$delta == 0
  vidx <- which(grid$valid & !zero)

  # emmeans runs on the valid reallocations only (the reference row is always
  # valid). A relative call with no such row has nothing to contrast, so
  # emmeans is skipped and every other estimate is NA.
  s <- NULL
  if (length(vidx) > 0 || !relative) {
    emm <- tc_emm_at_points(
      object,
      comps = rbind(reference, grid$comps[vidx, , drop = FALSE]),
      at = at, conf_level = conf_level, category = category
    )
    s <- if (relative) {
      # `"trt.vs.ctrl"` defaults to a Dunnett-type ("dunnettx") multiplicity
      # adjustment across every reallocation contrasted against the same
      # reference row, which both inflates p-values and widens the CI. Each
      # row here is its own independent Wald test, not a
      # simultaneous-coverage family, so adjustment must be off explicitly.
      summary(emmeans::contrast(emm, method = "trt.vs.ctrl", ref = 1),
              infer = c(TRUE, TRUE), level = conf_level,
              type = "response", adjust = "none")
    } else {
      summary(emm, infer = c(TRUE, FALSE), level = conf_level,
              type = "response")
    }
  }

  # The absolute output leads with the reference row (from/to NA, delta 0).
  keys <- if (relative) {
    tibble::tibble(from = grid$from, to = grid$to, delta = grid$delta)
  } else {
    tibble::tibble(from = c(NA_character_, grid$from),
                   to = c(NA_character_, grid$to),
                   delta = c(0, grid$delta))
  }
  # The rows the emmeans summary fills, in order; every other row stays NA.
  slots <- if (relative) vidx else c(1L, vidx + 1L)
  values <- if (!is.null(s)) {
    # A zero standard error (a zero reallocation) leaves the Wald test
    # undefined: emmeans returns 0/0, reported here as NA rather than NaN.
    undefined <- !is.na(s$SE) & s$SE == 0
    list(
      estimate = tc_estimate_col(s),
      std.error = s$SE,
      df = s$df,
      conf.low = tc_ci_col(s, "lower"),
      conf.high = tc_ci_col(s, "upper"),
      statistic = replace(if (!is.null(s$t.ratio)) s$t.ratio else s$z.ratio,
                          undefined, NA_real_),
      p.value = replace(s$p.value, undefined, NA_real_)
    )
  }
  cols <- c("estimate", "std.error", "df", "conf.low", "conf.high",
            if (relative) c("statistic", "p.value"))
  filled <- purrr::map(rlang::set_names(cols), function(col) {
    replace(rep(NA_real_, nrow(keys)), slots, values[[col]])
  })
  if (any(zero)) {
    ref_value <- if (isTRUE(tc_report_meta(object$model_type,
                                           object$family)$exponentiate)) 1 else 0
    cols <- c("estimate", "conf.low", "conf.high")
    filled[cols] <- purrr::map(filled[cols], function(v) {
      replace(v, zero, ref_value)
    })
    filled$std.error[zero] <- 0
  }
  dplyr::bind_cols(keys, tibble::as_tibble(filled))
}

# Directed reallocation schemes ------------------------------------------------

# The label of the pooled "all other parts" side of a one-to-all reallocation,
# as it appears in the from/to columns.
tc_rest_key <- function() "(rest)"

# The display label of that side, in summaries and figures.
tc_rest_label <- function() "All other parts"

# All ordered from/to pairs, or the validated user-supplied ones.
tc_realloc_pairs <- function(parts, from, to, call = rlang::caller_env()) {
  if (is.null(from) != is.null(to)) {
    cli::cli_abort("Supply both {.arg from} and {.arg to}, or neither.", call = call)
  }
  if (is.null(from)) {
    pairs <- expand.grid(
      from = parts, to = parts,
      stringsAsFactors = FALSE, KEEP.OUT.ATTRS = FALSE
    )
    return(pairs[pairs$from != pairs$to, , drop = FALSE])
  }
  bad <- setdiff(unique(c(from, to)), parts)
  if (length(bad) > 0) {
    cli::cli_abort(
      c(
        "{.arg from}/{.arg to} names {.val {bad}}, which {?is/are} not in the fitted composition.",
        "i" = "Composition parts are {.val {parts}}."
      ),
      call = call
    )
  }
  n <- max(length(from), length(to))
  if (!all(c(length(from), length(to)) %in% c(1, n))) {
    cli::cli_abort(
      "{.arg from} (length {length(from)}) and {.arg to} (length {length(to)}) must have equal lengths, or be scalars.",
      call = call
    )
  }
  pairs <- data.frame(from = rep_len(from, n), to = rep_len(to, n))
  if (any(pairs$from == pairs$to)) {
    cli::cli_abort("{.arg from} and {.arg to} must differ for every reallocation.", call = call)
  }
  pairs
}

# Validate the focal spec for one-to-all reallocation: at most one of from/to
# names the focal part(s) (neither: every part receives in turn); the other
# side is the pooled rest.
tc_realloc_focals <- function(parts, from, to, call = rlang::caller_env()) {
  if (!is.null(from) && !is.null(to)) {
    cli::cli_abort(
      c(
        "In {.code scheme = \"one-to-all\"} supply only one of {.arg from} or {.arg to}.",
        "i" = "The other side is implicitly all remaining parts."
      ),
      call = call
    )
  }
  spec <- if (!is.null(to)) {
    list(parts = to, role = "to")
  } else if (!is.null(from)) {
    list(parts = from, role = "from")
  } else {
    list(parts = parts, role = "to") # default: each part as focal receiver
  }
  bad <- setdiff(spec$parts, parts)
  if (length(bad) > 0) {
    cli::cli_abort(
      c(
        "{.arg from}/{.arg to} names {.val {bad}}, which {?is/are} not in the fitted composition.",
        "i" = "Composition parts are {.val {parts}}."
      ),
      call = call
    )
  }
  gains <- spec$role == "to"
  data.frame(
    focal = spec$parts, role = spec$role,
    from = if (gains) tc_rest_key() else spec$parts,
    to = if (gains) spec$parts else tc_rest_key(),
    stringsAsFactors = FALSE
  )
}

# What each directed `scheme` means, read by predict_realloc() and
# plot_realloc() alike. Fields:
#   moves(parts, from, to, call)  validate `from`/`to` into a data frame with
#                                 one row per directed reallocation, carrying
#                                 its `from`/`to` output labels
#   shift(reference, move)        the change to `reference` per unit moved
#                                 along one such row (it sums to zero, so the
#                                 composition stays closed)
#   infeasible(move, delta, reference, call)  abort for a row that `delta`
#                                 closure units drive non-positive
tc_realloc_schemes <- list(
  "one-to-one" = list(
    moves = tc_realloc_pairs,
    shift = function(reference, move) {
      unit <- stats::setNames(numeric(length(reference)), names(reference))
      unit[[move$from]] <- -1
      unit[[move$to]] <- 1
      unit
    },
    infeasible = function(move, delta, reference, call) {
      cli::cli_abort(
        c(
          "Reallocating {delta} from {.val {move$from}} to {.val {move$to}} makes a part non-positive.",
          "i" = "The reference composition has {.val {move$from}} = {round(reference[[move$from]], 1)} and {.val {move$to}} = {round(reference[[move$to]], 1)}; choose {.arg delta} within those bounds."
        ),
        call = call
      )
    }
  ),
  # The focal part gains (role "to") or gives (role "from") the time moved; the
  # remaining parts give or gain them in proportion to their reference values.
  "one-to-all" = list(
    moves = tc_realloc_focals,
    shift = function(reference, move) {
      gains <- move$role == "to"
      fi <- match(move$focal, names(reference))
      rest <- setdiff(seq_along(reference), fi)
      w <- reference[rest] / sum(reference[rest])
      unit <- stats::setNames(numeric(length(reference)), names(reference))
      unit[fi] <- if (gains) 1 else -1
      unit[rest] <- if (gains) -w else w
      unit
    },
    infeasible = function(move, delta, reference, call) {
      cli::cli_abort(
        c(
          "Reallocating {delta} (proportional, {move$role} {.val {move$focal}}) makes a part non-positive.",
          "i" = "Choose a smaller {.arg delta}, or use {.code on_invalid = \"na\"}."
        ),
        call = call
      )
    }
  )
)

# Expand moves x delta into reallocated compositions (reference + delta x the
# scheme's per-unit shift), flagging any row that leaves a part non-positive.
tc_realloc_grid <- function(reference, moves, delta, scheme,
                            on_invalid = "error",
                            call = rlang::caller_env()) {
  spec <- tc_realloc_schemes[[scheme]]
  grid <- expand.grid(
    move = seq_len(nrow(moves)), delta = delta,
    KEEP.OUT.ATTRS = FALSE
  )
  unit <- purrr::map(seq_len(nrow(moves)), function(i) {
    spec$shift(reference, moves[i, , drop = FALSE])
  })
  # Matrix algebra, kept in base R: one reallocated composition per row.
  unit <- do.call(rbind, unit)[grid$move, , drop = FALSE]
  comps <- matrix(
    reference, nrow = nrow(grid), ncol = length(reference),
    byrow = TRUE, dimnames = list(NULL, names(reference))
  ) + grid$delta * unit
  bad <- rowSums(comps <= 0) > 0
  if (any(bad) && on_invalid == "error") {
    b <- which(bad)[1]
    spec$infeasible(moves[grid$move[b], , drop = FALSE], grid$delta[b],
                    reference, call)
  }
  list(from = moves$from[grid$move], to = moves$to[grid$move],
       delta = grid$delta, comps = comps, valid = !bad)
}

# Custom multi-part reallocation ---------------------------------------------

# Normalize `realloc` (a named vector, a list of them, or a data frame with one
# row per scenario) into a list of validated, parts-ordered change vectors.
tc_realloc_normalise <- function(realloc, parts, call = rlang::caller_env()) {
  # Coerce every accepted input to a list of named numeric vectors.
  scenarios <- if (is.data.frame(realloc)) {
    if (nrow(realloc) < 1) {
      cli::cli_abort("{.arg realloc} data frame must have at least one row (scenario).", call = call)
    }
    purrr::map(seq_len(nrow(realloc)), function(i) {
      stats::setNames(as.numeric(unlist(realloc[i, , drop = FALSE])),
                      names(realloc))
    })
  } else if (is.list(realloc)) {
    if (length(realloc) < 1) {
      cli::cli_abort("{.arg realloc} list must contain at least one scenario.", call = call)
    }
    realloc
  } else {
    list(realloc)
  }

  # Validate each scenario and expand it to a full, parts-ordered vector
  # (unnamed parts 0). lapply(), not purrr::map(): map() re-raises an abort
  # from inside its function wrapped in an "In index: i" message, burying the
  # scenario's own.
  lapply(seq_along(scenarios), function(i) {
    vec <- scenarios[[i]]
    label <- if (length(scenarios) > 1) paste0(" (scenario ", i, ")") else ""
    if (!is.numeric(vec) || is.null(names(vec)) || anyNA(names(vec)) ||
        any(names(vec) == "")) {
      cli::cli_abort(
        "{.arg realloc}{label} must be a named numeric vector (one signed change per part).",
        call = call
      )
    }
    if (anyNA(vec) || any(!is.finite(vec))) {
      cli::cli_abort(
        "{.arg realloc}{label} must contain only finite numeric values.",
        call = call
      )
    }
    if (anyDuplicated(names(vec))) {
      dup <- unique(names(vec)[duplicated(names(vec))])
      cli::cli_abort(
        "{.arg realloc}{label} names {cli::qty(dup)}part{?s} {.val {dup}} more than once.",
        call = call
      )
    }
    bad <- setdiff(names(vec), parts)
    if (length(bad) > 0) {
      cli::cli_abort(
        c("{.arg realloc}{label} names {.val {bad}}, which {?is/are} not in the fitted composition.",
          "i" = "Composition parts are {.val {parts}}."),
        call = call
      )
    }
    imbalance <- sum(vec)
    if (abs(imbalance) > 1e-8) {
      cli::cli_abort(
        c("{.arg realloc}{label} must sum to 0 (a closed reallocation), but sums to {signif(imbalance, 4)}.",
          "i" = "Balance the time moved in and out so the day still totals the composition closure."),
        call = call
      )
    }
    if (!any(vec > 0) || !any(vec < 0)) {
      cli::cli_abort(
        "{.arg realloc}{label} must have at least one positive and one negative entry (time both taken and given).",
        call = call
      )
    }
    full <- stats::setNames(rep(0, length(parts)), parts)
    full[names(vec)] <- vec
    full
  })
}

# Expand the validated scenarios into reallocated compositions (reference +
# scenario, re-closed to the reference total), flagging positivity. Each
# scenario gets a compact label per side (receivers "SPT +30", donors
# "SB -15, LPA -15", largest move first) and a `delta` equal to the total
# time moved, in closure units.
tc_realloc_grid_custom <- function(reference, scenarios, on_invalid = "error",
                                   call = rlang::caller_env()) {
  parts <- names(reference)
  total <- sum(reference)
  vecs <- purrr::map(unname(scenarios), function(s) s[parts])
  comps <- purrr::map(vecs, function(vec) {
    raw <- reference + vec
    # Re-close to the original total (guards tiny numeric drift; a genuine
    # zero-sum scenario keeps the total exactly).
    raw / sum(raw) * total
  })
  comps <- do.call(rbind, comps)
  sides <- purrr::map(vecs, function(vec) {
    purrr::map_chr(c(to = 1, from = -1), function(sgn) {
      x <- vec[sign(vec) == sgn]
      x <- x[order(-abs(x))]
      amount <- purrr::map_chr(abs(x), format, trim = TRUE, scientific = FALSE)
      paste(paste0(names(x), " ", if (sgn > 0) "+" else "-", amount),
            collapse = ", ")
    })
  })
  from <- purrr::map_chr(sides, "from")
  to <- purrr::map_chr(sides, "to")
  bad <- rowSums(comps <= 0) > 0
  if (any(bad) && on_invalid == "error") {
    b <- which(bad)[1]
    neg <- which(comps[b, ] <= 0)
    cli::cli_abort(
      c("Custom reallocation {.val {to[b]}} from {.val {from[b]}} drives {.val {parts[neg]}} non-positive.",
        "i" = "Move less time, or use {.code on_invalid = \"na\"}."),
      call = call
    )
  }
  list(from = from, to = to,
       delta = purrr::map_dbl(vecs, function(vec) sum(vec[vec > 0])),
       comps = comps, valid = !bad)
}

# emmeans ----------------------------------------------------------------------

# timecoda never prints emmeans' interaction note; a path that averages over an
# interacting covariate says so once itself.
tc_quiet_emmeans <- function(frame = rlang::caller_env()) {
  rlang::local_options(
    emmeans = utils::modifyList(
      getOption("emmeans", list()), list(msg.interaction = FALSE)
    ),
    .frame = frame
  )
}

# The engine's emmeans reference grid by its registry row: the one place the
# engine fields (emm_via, emm_mode, emm_df, link_tran, response) are applied.
# `...` is `at` and the caller's other ref_grid()/qdrg() arguments; `category`
# names a multinomial fit's non-reference outcome level. `tc_iv_ref_slope()`
# mirrors this for `emtrends()`.
tc_ref_grid <- function(object, ..., meta = tc_model_meta(object$model_type),
                        category = NULL) {
  if (!is.null(category)) {
    # One multinomial category's baseline-relative log-odds model: its slice
    # of the coefficients ("<level>~<term>") and covariance. Contrasting on
    # this scale and exponentiating (link_tran = "log") gives the relative-risk
    # ratio.
    b <- stats::coef(object$model)
    idx <- which(sub("~.*$", "", names(b)) == category)
    coef_c <- stats::setNames(b[idx], sub("^[^~]*~", "", names(b)[idx]))
    vcov_c <- as.matrix(stats::vcov(object$model))[idx, idx, drop = FALSE]
    dimnames(vcov_c) <- list(names(coef_c), names(coef_c))
    grid <- emmeans::qdrg(formula = stats::formula(object$model),
                          data = object$data, coef = coef_c, vcov = vcov_c, ...)
  } else if (meta$emm_via == "qdrg") {
    # lmrob has no native emmeans support: the grid is built from its
    # coefficients and covariance, and emmeans still owns the contrast maths.
    grid <- emmeans::qdrg(
      formula = stats::formula(object$model), data = object$data,
      coef = stats::coef(object$model),
      vcov = as.matrix(stats::vcov(object$model)),
      df = stats::df.residual(object$model), ...
    )
  } else {
    args <- list(object = object$model, data = object$data, ...)
    if (!is.null(meta$emm_mode)) args$mode <- meta$emm_mode
    grid <- do.call(emmeans::ref_grid, args)
  }
  # The engine's forced df (rlmer: asymptotic, as robustlmm reports none).
  if (!is.null(meta$emm_df)) grid <- stats::update(grid, df = meta$emm_df)
  # Tag the link so `type = "response"` back-transforms the contrast: an
  # ordinal latent-scale difference exponentiates to a cumulative odds ratio.
  if (!is.null(meta$link_tran)) grid <- stats::update(grid, tran = meta$link_tran)
  # A survival response is not a transformation, but emmeans reads a
  # namespaced `survival::Surv()` left-hand side as one (`tran2`) and would
  # then back-transform twice.
  if (meta$response == "surv") grid <- stats::update(grid, tran2 = NULL)
  grid
}

# One emmeans reference grid per composition (covariates averaged), combined
# into a single emmGrid so reallocations can be contrasted against row 1.
tc_emm_at_points <- function(object, comps, at, conf_level,
                             category = NULL) {
  meta <- tc_model_meta(object$model_type)
  ilr <- tc_pivot_coords(as.data.frame(comps), pivot = object$pivot)
  tc_quiet_emmeans()
  # Coenders T-space fit (include_total = TRUE): pin the total covariate at the
  # reference composition's value in every grid so a reallocation stays a
  # constant-total contrast (a pure proportional shift among the parts). Held
  # identical across the reference and reallocated rows, t_total cancels in the
  # relative contrast and fixes the total for the absolute prediction.
  if (!is.null(object$total_ref)) {
    at <- c(at, list(t_total = object$total_ref))
  }
  # lapply(), not purrr::map(), so an emmeans abort keeps its own message.
  grids <- lapply(seq_len(nrow(ilr)), function(i) {
    at_i <- c(as.list(ilr[i, ]), at)
    rg <- tc_ref_grid(object, at = at_i, options = list(level = conf_level),
                      meta = meta, category = category)
    # Collapse covariate factor levels by averaging: one row per composition.
    emmeans::emmeans(rg, stats::as.formula("~ ilr1"))
  })
  em <- do.call(rbind, c(grids, list(adjust = "none")))
  # Each row is labelled by its `ilr1` value, so two compositions sharing it
  # (a zero or repeated `delta`, two equal reference parts, symmetric grid
  # points) would give emmeans::contrast() duplicate row names. The label is
  # display only (each row's prediction is fixed in the grid), so make it
  # unique here, once, for every consumer.
  em@grid[["ilr1"]] <- seq_len(nrow(em@grid))
  em
}

# One fit's predictions at the candidates, on one emmeans grid per pass with
# the coordinates tied: `at` gives a pass's values of one key coordinate, and
# a `cov.reduce` formula per other coordinate looks up the candidate holding
# that key. Returns the raw `link` (ranked by the caller), the `estimate`
# reported (the link, or for `type = "relative"` the link-scale difference
# from `from`) and the `response` thresholds read (the predicted response, or
# the effect versus `from` on the fit's effect scale).
tc_optimal_predict <- function(fit, cand, pins, type, from, category) {
  absolute <- type == "absolute"
  coords <- colnames(cand$z)
  # Keys are rounded to 12 decimals for the split, the `at` values, the lookup
  # and the match back alike: equal coordinates can differ by floating-point
  # noise, and rounding in one place but not in another gives silently wrong
  # values. The key is the coordinate whose rounded values repeat least.
  rounded <- round(cand$z, 12)
  repeats <- purrr::map_int(seq_along(coords), function(j) {
    max(tabulate(match(rounded[, j], unique(rounded[, j]))))
  })
  key <- coords[which.min(repeats)]
  tied <- setdiff(coords, key)
  # A candidate's pass is its occurrence number among equal keys, so no pass
  # holds a key twice: emmeans refuses a duplicated `at` value.
  seen <- match(rounded[, key], unique(rounded[, key]))
  pass <- stats::ave(seq_along(seen), seen, FUN = seq_along)
  # Factor covariates neither pinned nor interacting are averaged as nuisance,
  # with equal weight as the per-point path does, keeping one row per
  # candidate. A variable is a factor when the model's own term classes say
  # so, so `factor(region)` over a numeric column is averaged as `region`.
  vars <- intersect(all.vars(reformulas::nobars(stats::formula(fit$model))[[3]]),
                    names(fit$data))
  cls <- attr(stats::terms(fit$model), "dataClasses")
  if (is.null(cls)) {
    # coxme has no model.frame() method: the default one evaluates its
    # `(1 | g)` term and warns, so it falls back to the data's column types.
    cls <- tryCatch(
      attr(attr(stats::model.frame(fit$model), "terms"), "dataClasses"),
      warning = function(w) NULL, error = function(e) NULL
    )
  }
  factors <- if (is.null(cls)) {
    vars[purrr::map_lgl(vars, function(v) {
      col <- fit$data[[v]]
      is.factor(col) || is.character(col) || is.logical(col)
    })]
  } else {
    labels <- names(cls)[cls %in% c("factor", "ordered", "character", "logical")]
    unlist(purrr::map(labels, function(l) {
      tryCatch(all.vars(str2lang(l)), error = function(e) l)
    }))
  }
  unpinned <- setdiff(intersect(vars, factors), c(coords, "t_total", names(pins)))
  nuisance <- setdiff(unpinned, tc_interacting_covariates(fit))
  # A grid's row count before averaging, an upper bound: the pass's rows plus
  # the two that `nuisance` adds, times each pinned variable's values, times
  # the distinct values of every unpinned factor (one that is not nuisance,
  # such as a factor interacting with another covariate, stays in the grid).
  n_levels <- purrr::map_int(unpinned, function(v) {
    length(unique(stats::na.omit(fit$data[[v]])))
  })
  tc_quiet_emmeans()
  # lapply(), not purrr::map(), so an emmeans abort keeps its own message.
  est <- lapply(seq_len(max(pass)), function(p) {
    rows <- which(pass == p)
    keys <- rounded[rows, key]
    values <- cand$z[rows, , drop = FALSE]
    # The pass's value of coordinate `j` at a matching key, and 0 elsewhere:
    # emmeans also calls it on the model data, where no key matches.
    .tc_tied <- function(v, j) {
      out <- values[match(round(v, 12), keys), j]
      out[is.na(out)] <- 0
      out
    }
    # One formula per tied coordinate, fitted by lm() on the data and
    # predicted on the grid through the formula's environment, this frame.
    # The coordinate is given by number: emmeans counts a string literal as
    # a variable and refuses the formula.
    frame <- rlang::current_env()
    reduce <- purrr::map(tied, function(nm) {
      stats::as.formula(
        paste0(nm, " ~ 0 + offset(.tc_tied(", key, ", ", match(nm, coords), "))"),
        env = frame
      )
    })
    args <- list(
      at = c(rlang::set_names(list(keys), key), pins,
             if (!is.null(fit$total_ref)) list(t_total = fit$total_ref)),
      category = category,
      rg.limit = (length(rows) + 2) * prod(lengths(pins)) * prod(n_levels)
    )
    if (length(tied) > 0) args$cov.reduce <- reduce
    if (length(nuisance) > 0) args$nuisance <- nuisance
    rg <- do.call(tc_ref_grid, c(list(fit), args))
    em <- emmeans::emmeans(rg, stats::reformulate(key))
    hit <- match(keys, round(em@grid[[key]], 12))
    if (anyNA(hit)) {
      cli::cli_abort("The {.pkg emmeans} grid rows do not match the candidates.",
                     .internal = TRUE)
    }
    cbind(rows, stats::predict(em)[hit],
          if (absolute) stats::predict(em, type = "response")[hit])
  })
  est <- do.call(rbind, est)
  est <- est[order(est[, 1]), , drop = FALSE]
  link <- est[, 2]
  response <- if (absolute) est[, 3]
  if (absolute) {
    return(list(link = link, estimate = link, response = response))
  }
  # `pins` carries the `by` level, so `from` is predicted under the same pins.
  link_from <- stats::predict(tc_emm_at_points(fit, t(from), at = pins,
                                               conf_level = 0.95,
                                               category = category))
  estimate <- link - link_from
  exponentiate <- isTRUE(tc_report_meta(fit$model_type, fit$family)$exponentiate)
  list(link = link, estimate = estimate,
       response = if (exponentiate) exp(estimate) else estimate)
}

# Column-name normalization for emmeans summaries -----------------------------

# The estimate column of an emmeans summary, whatever emmeans named it.
tc_estimate_col <- function(s) {
  candidates <- c(
    "ratio", "odds.ratio", "rate", "response", "prob", "emmean", "estimate"
  )
  col <- intersect(candidates, names(s))[1]
  if (is.na(col)) {
    cli::cli_abort("The {.pkg emmeans} summary has no estimate column.",
                   .internal = TRUE)
  }
  s[[col]]
}

# The lower or upper confidence-limit column of an emmeans summary made with
# `infer[1] = TRUE`, whatever emmeans named it.
tc_ci_col <- function(s, side = c("lower", "upper")) {
  side <- match.arg(side)
  candidates <- if (side == "lower") {
    c("lower.CL", "asymp.LCL")
  } else {
    c("upper.CL", "asymp.UCL")
  }
  col <- intersect(candidates, names(s))[1]
  if (is.na(col)) {
    cli::cli_abort("The {.pkg emmeans} summary has no {side} confidence-limit column.",
                   .internal = TRUE)
  }
  s[[col]]
}

# The effect a reallocation reports, by model type and output mode.
tc_effect_label <- function(object, type) {
  if (type == "absolute") {
    return("Predicted response")
  }
  tc_report_meta(object$model_type, object$family)$measure
}

# Abort unless `object` is a timecoda model fit.
tc_check_fit <- function(object, call = rlang::caller_env()) {
  if (!inherits(object, "timecoda_fit")) {
    cli::cli_abort(
      "{.arg object} must be a {.cls timecoda_fit} created by {.fn comp_iv}.",
      call = call
    )
  }
  invisible(object)
}

# Check `category` against `fits`, a list named by outcome: refused when no fit
# is multinomial, and for each multinomial fit required and one of its
# non-reference outcome levels; `action`, the verb ending the message that asks
# for it, is "search" or "draw". Returns which fits are multinomial.
tc_check_category <- function(fits, category, action,
                              call = rlang::caller_env()) {
  per_category <- purrr::map_lgl(fits, function(fit) {
    tc_model_meta(fit$model_type)$per_category
  })
  if (!is.null(category) && !any(per_category)) {
    cli::cli_abort("{.arg category} applies to multinomial fits only.",
                   call = call)
  }
  # A loop, not purrr::iwalk(), so an abort keeps its own message.
  for (k in names(fits)[per_category]) {
    levs <- unique(sub("~.*$", "", names(stats::coef(fits[[k]]$model))))
    if (is.null(category)) {
      cli::cli_abort(c(
        "The {.val {k}} model is multinomial, so {.arg category} must name the outcome level to {action}.",
        "i" = "Its non-reference levels are {.val {levs}}."
      ), call = call)
    }
    if (!rlang::is_string(category) || !category %in% levs) {
      cli::cli_abort(
        "{.arg category} must be one of the {.val {k}} model's non-reference outcome levels, {.val {levs}}.",
        call = call
      )
    }
  }
  per_category
}
