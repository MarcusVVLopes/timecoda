# Module 1: summarizing a comp_iv object.
#
# `summary.comp_iv()` returns a `summary_comp_iv` object: a "Global Composition"
# omnibus row plus one "<part> vs remaining" row per part, each carrying the
# back-transformed effect (a ratio for a log- or logit-link GLM, Cox and the
# ordinal and multinomial engines; a mean difference otherwise), its Wald
# confidence interval and its p-value. Each per-part effect is the `ilr1` coefficient of
# the rotated model in `object$models` (that part pivoted into the first
# coordinate), or, when that fit has no plain linear `ilr1` term, its slope in
# `ilr1` at the reference composition. Printing the object shows a text report;
# `pub_ready()` builds the publication `gt` table from it.

#' Summarize a composition-as-predictor model
#'
#' Tabulates the omnibus test of the whole composition of a [comp_iv()] fit and
#' one "part vs remaining" effect per part. Each effect is the rotated model's
#' `ilr1` coefficient, or its slope at the reference composition (see Details),
#' back-transformed to the model's effect scale. Printing the result shows a
#' text report, as with [summary.lm()], and [pub_ready()] builds a publication
#' `gt` table from it.
#'
#' @param object A `comp_iv` object from [comp_iv()].
#' @param conf_level Confidence level for the per-part Wald intervals.
#' @param ... Unused; present for S3 consistency.
#'
#' @details The effect scale is a hazard ratio for a Cox model, a cumulative
#'   odds ratio for an ordinal model and a relative-risk ratio for a multinomial
#'   model. A GLM's scale follows its link: an odds ratio (logit), a rate, risk
#'   or mean ratio (log) or a mean difference (identity). Every linear engine
#'   reports a mean difference.
#'
#'   For mixed models the per-part effects come from one refit per part (stored
#'   at fit time in `object$models`), so building this summary does no extra
#'   fitting. For an `"lmer"` fit, each per-part standard error, interval and
#'   p-value uses the degrees-of-freedom method that `emmeans` uses for
#'   [predict_realloc()] on the same fit. The default is the Kenward-Roger
#'   approximation through \pkg{pbkrtest}. The Satterthwaite approximation
#'   through \pkg{lmerTest} is used when `emmeans::emm_options(lmer.df =)`
#'   selects it, or when \pkg{pbkrtest} is absent, disabled or the fit has more
#'   observations than `emmeans::emm_options(pbkrtest.limit =)`. Both
#'   \pkg{pbkrtest} and \pkg{lmerTest} are optional, and when neither applies,
#'   or `emmeans::emm_options(lmer.df = "asymptotic")` is set, the degrees of
#'   freedom are asymptotic (normal). `"glmer"`, `"coxme"` and `"rlmer"` rows
#'   are asymptotic, as their `emmeans` predictions are (\pkg{robustlmm} reports
#'   no degrees of freedom).
#'
#'   When the formula gives the composition non-linear terms (see [ilr()]), a
#'   product of coordinates or no main effect, a rotated fit has no single
#'   `ilr1` coefficient to read. Each per-part effect is then the slope of its
#'   linear predictor in `ilr1` at the reference composition, from
#'   [emmeans::emtrends()], back-transformed in the same way, and the printed
#'   report says so. With a composition-by-covariate interaction as well, that
#'   slope is the equal-weight average over the covariate's levels. An engine
#'   for which `emtrends()` gives no such slope (`"mblogit"`, `"lmrob"`) reports
#'   `NA` rows.
#'
#' @return An object of class `summary_comp_iv`, printed as a text report by
#'   [print.summary_comp_iv()] and passed to [pub_ready()] for a `gt` table. A
#'   list with elements:
#'   \describe{
#'     \item{`global`}{The one-row omnibus composition test.}
#'     \item{`per_part`}{The per-part effect table: one row per part, or per
#'       part and non-reference outcome category for a multinomial fit. It
#'       shares the columns `parameter`, `label`, `estimate`, `conf.low`,
#'       `conf.high`, `p.value` with `global`, whose
#'       `estimate`/`conf.low`/`conf.high` are `NA`.
#'       `label` carries the part's display label (or its raw name). A
#'       per-category model (multinomial) adds `response_level`.}
#'     \item{`moderation`}{The joint moderation Wald test when the fit was built
#'       with a `moderator`, else `NULL`. A moderated fit also gives `per_part`
#'       three further columns (`mod_statistic`, `mod_df`, `mod_p.value`). They
#'       hold the joint Wald test of that part's own `ilr1 : moderator`
#'       coefficient(s), which asks whether the part's effect varies by
#'       moderator level; the `moderation` test pools every part. Like the
#'       per-part effect rows, they are descriptive and carry no multiplicity
#'       adjustment. An unmoderated fit omits them.}
#'     \item{`curvature_test`}{The joint Wald test of the squares and products
#'       of the coordinates when the fit holds the full quadratic (see the
#'       `curvature_test` element of [comp_iv()]), else `NULL`.}
#'     \item{`convergence`}{A list of the fit's `status` one-row tibble
#'       (`converged`, `singular`, `n_warnings`) and its `warnings` tibble of
#'       captured backend conditions.}
#'     \item{`icc`}{For a clustered fit, a tibble with one row per grouping term
#'       (`group`, `variance`, `measure`, `value`), read from the fit's own
#'       variance components. For `"lmer"` and `"rlmer"`, `value` is each term's
#'       share of the total variance, the intraclass correlation
#'       \eqn{v_g / (\sum_h v_h + \sigma^2)}, where the sum runs over all
#'       grouping terms \eqn{h}. A logit- or probit-link binomial `"glmer"` and
#'       an ordinal `"clmm"` give its latent-scale form, with \eqn{\pi^2/3}
#'       (logit) or 1 (probit) in place of \eqn{\sigma^2}. A `"coxme"` fit gives
#'       the median hazard ratio \eqn{\exp(\sqrt{2 v_g}\, \Phi^{-1}(0.75))}: the
#'       median, over pairs of clusters, of the hazard ratio between two
#'       individuals with identical covariates, one from each cluster, taken
#'       from the higher-risk cluster to the lower-risk one (Austin et al.
#'       2017). `measure` is `"not available"` and `value`
#'       `NA` for the other clustered engines, for a log-link binomial
#'       `"glmer"`, which has no latent-scale form, and for a `"glmer"` of any
#'       other family. `NULL` for a fit without a random term.}
#'     \item{`total_test`, `total_beta`}{With `include_total = TRUE`, the
#'       one-row Wald test of the `t_total` coefficient and its fitted value on
#'       the link scale (the first category's for a multinomial fit, whose test
#'       is joint); `NULL` otherwise.}
#'     \item{`per_part_slopes`}{Present, and `TRUE`, only when at least one
#'       per-part row is a slope at the reference composition rather than a
#'       coefficient.}
#'     \item{`per_part_no_slope`}{The reason, as text, why
#'       `emmeans::emtrends()` gave no slope and those rows are `NA`. Present
#'       only in that case.}
#'     \item{`labels`, `effect_measure`, `conf_level`, `moderator`,
#'       `fit_label`, `total`, `model_type`, `omnibus`}{Carried through for
#'       [print.summary_comp_iv()] and [pub_ready()].}
#'     \item{`units`}{The composition's unit label, printed after `total`;
#'       `NULL` when none was declared.}
#'   }
#'
#' @examples
#' cd <- declare_composition(timeuse, c("MVPA", "LPA", "SB", "SPT"),
#'                           closure = 1440)
#' m <- comp_iv(cd, outcome = "obesity", model = "glm", family = "binomial")
#' s <- summary(m)
#' s$per_part
#'
#' @references
#' Austin, P. C., Wagner, P., & Merlo, J. (2017). The median hazard ratio: A
#'   useful measure of variance and general contextual effects in multilevel
#'   survival analysis. *Statistics in Medicine*, 36(6), 928-938.
#'   \doi{10.1002/sim.7188}
#'
#' Goldstein, H., Browne, W., & Rasbash, J. (2002). Partitioning variation in
#'   multilevel models. *Understanding Statistics*, 1(4), 223-231.
#'   \doi{10.1207/S15328031US0104_02}
#'
#' Kenward, M. G., & Roger, J. H. (1997). Small sample inference for fixed
#'   effects from restricted maximum likelihood. *Biometrics*, 53(3), 983-997.
#'   \doi{10.2307/2533558}
#'
#' @seealso [comp_iv()] to fit the model; [plot.comp_iv()] for the
#'   predicted-outcome surface; [predict_realloc()] for reallocation
#'   predictions; [pub_ready()] for a formatted `gt` table built from this
#'   summary.
#' @family composition as predictor
#' @export
summary.comp_iv <- function(object, conf_level = 0.95, ...) {
  if (!inherits(object, "comp_iv")) {
    cli::cli_abort("{.arg object} must be a {.cls comp_iv} object.")
  }
  if (!is.numeric(conf_level) || length(conf_level) != 1 ||
      is.na(conf_level) || conf_level <= 0 || conf_level >= 1) {
    cli::cli_abort("{.arg conf_level} must be a single number in (0, 1).")
  }
  labels <- object$comp$labels
  report <- tc_report_meta(object$model_type, object$family)
  slopes <- !purrr::map_lgl(object$models, tc_iv_plain_ilr1,
                            model_type = object$model_type,
                            coords = tc_coord_names(length(object$parts)))
  per_part <- tc_iv_per_part_table(object, slopes, report$exponentiate,
                                   conf_level)
  no_slope <- unique(stats::na.omit(per_part$.no_slope))
  per_part$.no_slope <- NULL

  # Between-cluster share per grouping term by the registry's `icc` measure
  # (intercept variance over the total, or coxme's median hazard ratio); else NA.
  icc <- NULL
  if (!is.null(object$random)) {
    icc_type <- tc_model_meta(object$model_type)$icc
    if (icc_type == "none") {
      # The grouping name of a plain term or of a bar term such as "(1 | school)".
      group <- trimws(sub("^.*\\|", "", gsub("[()]", "", object$random)))
      variance <- NA_real_
      measure <- "not available"
      value <- NA_real_
    } else {
      vc <- if (icc_type == "coxme") {
        coxme::VarCorr(object$model)
      } else {
        lme4::VarCorr(object$model)
      }
      group <- names(vc)
      if (icc_type == "coxme") {
        variance <- purrr::map_dbl(vc, 1L)
        measure <- "median hazard ratio"
        value <- exp(sqrt(2 * variance) * stats::qnorm(0.75))
      } else {
        variance <- purrr::map_dbl(vc, function(x) x[1L, 1L])
        # The residual: sigma^2, or the latent pi^2/3 (logit) or 1 (probit),
        # by the link of a binomial glmer or of a clmm (which has no family).
        link <- if (icc_type == "clmm") {
          object$model$link
        } else {
          fam <- stats::family(object$model)
          if (identical(fam$family, "binomial")) fam$link
        }
        residual <- if (icc_type == "lmm") {
          stats::sigma(object$model)^2
        } else if (!is.null(link)) {
          switch(link, logit = pi^2 / 3, probit = 1, NA_real_)
        } else {
          NA_real_
        }
        measure <- if (icc_type == "lmm") {
          "ICC"
        } else if (is.na(residual)) {
          "not available"
        } else {
          "latent ICC"
        }
        value <- variance / (sum(variance) + residual)
      }
    }
    icc <- tibble::tibble(
      group = group, variance = unname(variance), measure = measure,
      value = unname(value)
    )
  }

  out <- structure(
    list(
      global = tibble::tibble(
        parameter = "Global Composition",
        label = "Global Composition",
        estimate = NA_real_,
        conf.low = NA_real_,
        conf.high = NA_real_,
        p.value = object$omnibus$p.value
      ),
      per_part = per_part,
      moderation = object$moderation,
      curvature_test = object$curvature_test,
      convergence = list(
        status = object$convergence,
        warnings = object$fit_warnings
      ),
      icc = icc,
      total_test = object$total_test,
      # The first t_total coefficient, read here because the summary object
      # does not carry the fitted model.
      total_beta = if (!is.null(object$total_test)) {
        b <- tc_model_coefs(object$model, object$model_type)$b
        unname(b[tc_ilr_coef_names(names(b), object$model_type, "t_total")[1]])
      },
      labels = tc_part_labels(object$parts, labels),
      effect_measure = report$measure,
      conf_level = conf_level,
      moderator = object$moderator,
      fit_label = report$label,
      total = object$total,
      units = object$comp$units,
      model_type = object$model_type,
      omnibus = object$omnibus
    ),
    class = "summary_comp_iv"
  )
  if (any(slopes)) out$per_part_slopes <- TRUE
  if (length(no_slope) > 0) out$per_part_no_slope <- no_slope[1]
  out
}

#' Print a composition-as-predictor summary
#'
#' Prints a [summary.comp_iv()] object as a text report, automatically at the
#' console, as base R's `print.summary.lm()` does. The report gives the model
#' description, the composition's parts and total (with its unit label, when
#' one was declared) and the omnibus composition test. When the fit carries
#' them, it adds the moderation test, the curvature test, the total-as-covariate
#' test, each grouping term's intraclass correlation or median hazard ratio and
#' any captured backend convergence messages. The effect table follows as plain
#' tibble output, with a closing note when the per-part effects are slopes at
#' the reference composition.
#'
#' @param x A `summary_comp_iv` object, as returned by [summary.comp_iv()].
#' @param ... Unused; present for S3 consistency.
#'
#' @return `x`, invisibly.
#'
#' @examples
#' cd <- declare_composition(timeuse, c("MVPA", "LPA", "SB", "SPT"),
#'                           closure = 1440)
#' m <- comp_iv(cd, outcome = "obesity", model = "glm", family = "binomial")
#' print(summary(m))
#'
#' @seealso [summary.comp_iv()] to build this object; [pub_ready()] for a
#'   formatted `gt` table.
#' @export
print.summary_comp_iv <- function(x, ...) {
  om <- x$omnibus
  print(glue::glue("
<comp_iv> compositional {x$fit_label}
  Parts: {glue::glue_collapse(x$labels, sep = ', ')} (total = {tc_with_units(x$total, x$units)})
  Omnibus composition test (Wald): chi-squared({om$df}) = {tc_fmt(om$statistic, 3)}, p = {format.pval(om$p.value, digits = 3)}
"))
  if (!is.null(x$moderation)) {
    md <- x$moderation
    print(glue::glue(
      "  Moderation by {x$moderator} (Wald): chi-squared({md$df}) = {tc_fmt(md$statistic, 3)}, p = {format.pval(md$p.value, digits = 3)}"
    ))
  }
  if (!is.null(x$curvature_test)) {
    cv <- x$curvature_test
    print(glue::glue(
      "  Curvature (Wald): chi-squared({cv$df}) = {tc_fmt(cv$statistic, 3)}, p = {format.pval(cv$p.value, digits = 3)}"
    ))
  }
  if (!is.null(x$total_test)) {
    print(glue::glue(
      "  Total (log geometric mean, Wald): beta = {tc_fmt(x$total_beta, 3)}, p = {format.pval(x$total_test$p.value, digits = 3, eps = 0.001)}"
    ))
  }
  if (!is.null(x$icc)) {
    known <- x$icc[x$icc$measure != "not available", ]
    if (nrow(known) > 0) {
      label <- c(
        "ICC" = "Intraclass correlation",
        "latent ICC" = "Latent intraclass correlation",
        "median hazard ratio" = "Median hazard ratio"
      )[known$measure]
      print(glue::glue("  {label} ({known$group}): {tc_fmt(known$value, 3)}"))
    } else {
      print(glue::glue("  Intraclass correlation: not available for this model"))
    }
  }
  tc_iv_print_convergence(x$convergence, x$model_type)
  cat("\n")
  print(glue::glue("Compositional effects ({x$effect_measure})"))
  print(dplyr::bind_rows(x$global, x$per_part))
  if (isTRUE(x$per_part_slopes)) {
    no_slope <- if (!is.null(x$per_part_no_slope)) {
      glue::glue(" NA rows: {x$per_part_no_slope}")
    } else {
      ""
    }
    print(glue::glue("
Per-part effects are slopes at the reference composition, from emmeans::emtrends(), as the
model's composition terms are non-linear or carry no main effect.{no_slope}
"))
  }
  invisible(x)
}

# Print the convergence section of summary(). For a clean iterative fit
# (converged TRUE) print a plain "Convergence: OK". When conditions were
# captured, list them in a cli block; the optimizer-retry hint is added only for
# the classes whose `control` takes an optimizer. lm/cox (converged NA and no
# warnings) print nothing, matching the header's silence.
tc_iv_print_convergence <- function(convergence, model_type) {
  status <- convergence$status
  warns <- convergence$warnings
  n_warn <- nrow(warns)
  singular <- !is.null(status) && isTRUE(status$singular)
  if (n_warn == 0 && !singular) {
    # Only claim OK when the engine actually reported success: NA means the
    # class has no convergence concept (lm, cox), FALSE means it failed.
    if (!is.null(status) && isTRUE(status$converged)) {
      print(glue::glue("  Convergence: OK"))
    } else if (!is.null(status) && isFALSE(status$converged)) {
      print(glue::glue("  Convergence: reported as NOT converged"))
    }
    return(invisible(NULL))
  }
  # Singularity is stated once: as the header when it is the only news, or as
  # the leading bullet when captured messages follow it.
  bullets <- character(0)
  if (n_warn > 0) {
    if (singular) {
      bullets <- c(bullets, "!" = "The fit is singular (a variance component is at the boundary).")
    }
    msgs <- glue::glue("[{warns$source}] {warns$message}")
    # The messages are backend text handed to cli as bullet templates, so a
    # brace inside one must be doubled or cli tries to evaluate it.
    msgs <- gsub("([{}])", "\\1\\1", as.character(msgs))
    names(msgs) <- rep("*", length(msgs))
    bullets <- c(bullets, msgs)
  }
  head_line <- if (n_warn > 0) {
    "Convergence: {n_warn} captured backend message{?s}:"
  } else {
    "Convergence: the fit is singular (a variance component is at the boundary)."
  }
  # Only the lme4-backed classes take an optimizer through `control`; telling a
  # plain glm to pass lme4::lmerControl() would be wrong advice.
  hint <- if (identical(tc_model_meta(model_type)$convergence, "lme4")) {
    c("i" = "Retry with a different optimizer via {.arg model_args} (e.g. {.code model_args = list(control = lme4::lmerControl(optimizer = \"bobyqa\"))}).")
  } else {
    character(0)
  }
  cli::cli_inform(c(head_line, bullets, hint))
  invisible(NULL)
}

# Per-part "vs remaining" effect table: for each part, pull the ilr1 estimate,
# SE, and p from its rotated model (or, where `slopes` marks it, the slope at
# the reference composition), form a Wald CI, and back-transform. When the fit
# has a moderator, also runs the per-part composition-by-moderator interaction
# test on that same rotated fit and appends its statistic/df/p.
tc_iv_per_part_table <- function(object, slopes, exponentiate,
                                 conf_level = 0.95) {
  alpha <- 1 - conf_level
  model_type <- object$model_type
  parts <- object$parts
  moderator <- object$moderator
  per_cat <- tc_model_meta(model_type)$per_category
  # lapply(), not purrr::map(), so an abort keeps its own message.
  lapply(seq_along(parts), function(k) {
    p <- parts[k]
    st <- if (slopes[[k]]) {
      tc_iv_ref_slope(object, k)
    } else {
      tc_iv_ilr1_stats(object$models[[p]], model_type)
    }
    # qt() with df = Inf is qnorm(), so this one line covers every model class:
    # the interval uses the same reference distribution as the p-value beside it.
    crit <- stats::qt(1 - alpha / 2, df = st$df)
    est <- st$est
    lo <- st$est - crit * st$se
    hi <- st$est + crit * st$se
    if (exponentiate) {
      est <- exp(est)
      lo <- exp(lo)
      hi <- exp(hi)
    }
    tib <- tibble::tibble(
      parameter = paste(p, "vs remaining"),
      label = tc_part_labels(p, object$comp$labels),
      estimate = est, conf.low = lo, conf.high = hi, p.value = st$p
    )
    if (per_cat) tib$response_level <- st$level
    if (!is.null(moderator)) {
      mod <- tc_moderation_test(object$models[[p]], model_type, moderator, parts,
                                per_part = TRUE)
      tib$mod_statistic <- mod$statistic
      tib$mod_df <- mod$df
      tib$mod_p.value <- mod$p.value
    }
    # Why a slope row is NA; summary.comp_iv() moves it off the table.
    tib$.no_slope <- if (slopes[[k]]) st[["no_slope"]][1] else NA_character_
    tib
  }) |> purrr::list_rbind()
}

# TRUE when a rotated fit's ilr1 effect is its plain `ilr1` coefficient: the
# coefficient exists, and every term in ilr1 is ilr1 alone or ilr1 crossed with
# covariates (`ilr1:sex`), not inside a call nor times another coordinate.
tc_iv_plain_ilr1 <- function(model, model_type, coords) {
  b_names <- names(tc_model_coefs(model, model_type)$b)
  labels <- attr(stats::terms(model), "term.labels")
  plain <- purrr::map_lgl(labels, function(label) {
    cmp <- strsplit(label, "(?<!:):(?!:)", perl = TRUE)[[1]]
    !"ilr1" %in% all.vars(str2lang(label)) ||
      (all(purrr::map_lgl(cmp, function(x) is.name(str2lang(x)))) &&
         !any(setdiff(cmp, "ilr1") %in% coords))
  })
  all(plain) && length(tc_ilr_coef_names(b_names, model_type, "ilr1")) > 0
}

# The slope of the rotated fit of pivot `k` in ilr1 at the reference
# composition, by emmeans::emtrends() on the fit's own data, in the shape of
# tc_iv_ilr1_stats(). Where emtrends() gives none (a per-category fit has no
# single linear predictor; lmrob has no emmeans support) the row is NA and
# `no_slope` says why.
tc_iv_ref_slope <- function(object, k) {
  model <- object$models[[k]]
  meta <- tc_model_meta(object$model_type)
  ref <- tc_pivot_coords(as.data.frame(object$ref_comp), pivot = k)
  at <- c(as.list(as.data.frame(ref)),
          if (!is.null(object$total_ref)) list(t_total = object$total_ref))
  # The default step is a forward difference that misses a quadratic's slope by
  # the step times the curvature (1-5 percent here); a small step is exact to
  # about 1e-6.
  # The registry's `emm_mode` and `emm_df`, as `tc_ref_grid()` applies them;
  # `emtrends()` takes the model, not a grid, so it cannot call it.
  args <- list(model, specs = ~ 1, var = "ilr1", at = at, delta.var = 1e-6)
  if (!is.null(meta$emm_mode)) args$mode <- meta$emm_mode
  # Averaging over an interacting covariate's levels is this row's documented
  # behaviour, so emmeans' note about it is off.
  tc_quiet_emmeans()
  s <- if (isTRUE(meta$per_category)) {
    "a per-category model has no single linear predictor to differentiate."
  } else if (meta$emm_via == "qdrg") {
    # No native emmeans support (lmrob): the row is NA, with emmeans' reason.
    tryCatch(summary(do.call(emmeans::emtrends, args), infer = c(TRUE, TRUE)),
             error = function(e) trimws(strsplit(conditionMessage(e), "\n")[[1]][1]))
  } else {
    trend <- do.call(emmeans::emtrends, args)
    # The engine's forced df (rlmer: asymptotic), as on its prediction grids.
    if (!is.null(meta$emm_df)) trend <- stats::update(trend, df = meta$emm_df)
    summary(trend, infer = c(TRUE, TRUE))
  }
  if (!is.character(s) && !is.finite(s$SE[1])) {
    s <- "emmeans::emtrends() gave no finite standard error."
  }
  no_slope <- if (is.character(s)) s else NA_character_
  if (is.character(s)) {
    s <- list(ilr1.trend = NA_real_, SE = NA_real_, p.value = NA_real_, df = NA_real_)
  }
  level <- if (isTRUE(meta$per_category)) {
    unique(sub("~.*$", "", grep("~", names(stats::coef(model)), value = TRUE)))
  } else {
    NA_character_
  }
  tibble::tibble(level = level, est = s$ilr1.trend, se = s$SE, p = s$p.value,
                 df = s$df, no_slope = no_slope)
}

# The ilr1 row(s) of a fitted model's coefficient table: estimate, standard
# error, p-value, and the degrees of freedom the table itself implies. One row
# per fit, except a per-category model (multinomial), which carries one ilr1
# coefficient per non-reference outcome level (`level`; NA otherwise).
#
# Coefficient naming comes from tc_ilr_coef_names(), so the multinomial
# "<level>~ilr1" convention lives in one place. Base-R matrix indexing is the
# matrix-algebra exception to the tidy-iteration convention.
tc_iv_ilr1_stats <- function(model, model_type) {
  per_cat <- isTRUE(tc_model_meta(model_type)$per_category)
  mat <- stats::coef(summary(model))
  cn <- colnames(mat)
  nms <- tc_ilr_coef_names(rownames(mat), model_type, "ilr1")
  if (length(nms) == 0) {
    cli::cli_abort(c(
      "Could not locate the {.val ilr1} coefficient in the fitted model.",
      "i" = "The per-part effect reads the linear {.code ilr1} term of each rotated fit."
    ))
  }
  idx <- match(nms, rownames(mat))

  est <- unname(mat[idx, if ("coef" %in% cn) "coef" else "Estimate"])
  se <- unname(mat[idx, if ("se(coef)" %in% cn) "se(coef)" else "Std. Error"])
  # The degrees of freedom come from whatever produced the p-value, so the
  # interval and the p printed beside it can never disagree:
  #   Pr(>|t|) -> t on the backend's residual df
  #   Pr(>|z|) -> asymptotic
  #   no Pr column (lmer/rlmer) -> asymptotic, because the p is derived from
  #   the normal just below; a linear mixed fit then takes the
  #   Kenward-Roger or Satterthwaite row when emmeans would, further down.
  # Keying this off the statistic column ("t value") would give lmer a t
  # interval beside a normal p, and would call df.residual() on an rlmerMod,
  # which errors rather than warns. qt(p, df = Inf) is exactly
  # qnorm(p), so one critical-value path still serves every model class.
  p_idx <- grep("^Pr", cn)
  p <- if (length(p_idx) >= 1) {
    unname(mat[idx, p_idx[1]])
  } else {
    2 * stats::pnorm(-abs(est / se))
  }
  df <- if (length(p_idx) >= 1 && grepl("t", cn[p_idx[1]], fixed = TRUE)) {
    rdf <- tryCatch(stats::df.residual(model), error = function(e) NULL)
    if (is.null(rdf) || length(rdf) != 1 || !is.finite(rdf)) Inf else rdf
  } else {
    Inf
  }
  # A linear mixed fit takes the standard error and df emmeans would use
  # instead, so this row and predict_realloc() rest on the same approximation.
  mixed <- tc_iv_mixed_df(model, nms)
  if (!is.null(mixed)) {
    se <- mixed$se
    df <- mixed$df
    p <- 2 * stats::pt(-abs(est / se), df = df)
  }
  tibble::tibble(
    level = if (per_cat) sub("~.*$", "", nms) else NA_character_,
    est = est, se = se, p = p, df = df
  )
}

# SE and df of one lmer coefficient by the rule emmeans applies to the same fit
# (Kenward-Roger, Satterthwaite or asymptotic); NULL for asymptotic.
tc_iv_mixed_df <- function(model, coef_name) {
  if (!inherits(model, "lmerMod")) {
    return(NULL)
  }
  opt <- emmeans::get_emm_option
  n <- lme4::getME(model, "N")
  pbkr_ok <- !isTRUE(opt("disable.pbkrtest")) && n <= opt("pbkrtest.limit") &&
    requireNamespace("pbkrtest", quietly = TRUE)
  lmert_ok <- !isTRUE(opt("disable.lmerTest")) && n <= opt("lmerTest.limit") &&
    requireNamespace("lmerTest", quietly = TRUE)
  mode <- match.arg(tolower(opt("lmer.df")),
                    c("satterthwaite", "kenward-roger", "asymptotic"))
  if (mode == "kenward-roger" && !pbkr_ok) {
    mode <- "satterthwaite"
  }
  if (mode == "satterthwaite" && !lmert_ok) {
    mode <- if (pbkr_ok) "kenward-roger" else "asymptotic"
  }
  if (mode == "satterthwaite") {
    row <- stats::coef(summary(lmerTest::as_lmerModLmerTest(model),
                               ddf = "Satterthwaite"))[coef_name, ]
    return(list(se = unname(row["Std. Error"]), df = unname(row["df"])))
  }
  if (mode == "asymptotic") {
    return(NULL)
  }
  l_vec <- matrix(as.numeric(names(lme4::fixef(model)) == coef_name), nrow = 1)
  v_adj <- pbkrtest::vcovAdj(model)
  list(
    se = sqrt(as.numeric(l_vec %*% as.matrix(v_adj) %*% t(l_vec))),
    df = pbkrtest::Lb_ddf(l_vec, as.matrix(stats::vcov(model)), v_adj)
  )
}

# Build summary_comp_iv's gt table: the global and per-part rows ("<label> vs
# remaining"), a response-level column for a per-category fit, a moderation
# p-value column for a moderated one, and source notes for slope rows, the
# curvature test and the total test. The compact heading style uses absolute
# pixel sizes, because gt scales percentages unreliably in LaTeX, where its
# default title renders at ~20pt. The styling is written out inline
# deliberately, so that no builder depends on a shared helper for it. Do not
# "de-duplicate" it.
tc_pub_ready_comp_iv <- function(x, digits = 2) {
  per_part <- x$per_part
  global <- x$global
  per_part$component <- paste(per_part$label, "vs remaining")
  global$component <- global$label
  has_cat <- "response_level" %in% names(per_part)
  has_mod <- "mod_p.value" %in% names(per_part)
  # bind_rows() fills the columns global lacks with NA on its own.
  disp <- dplyr::bind_rows(global, per_part)
  disp <- dplyr::transmute(
    disp,
    parameter = .data$component,
    response_level = if (has_cat) .data$response_level else NULL,
    estimate = ifelse(is.na(.data$estimate), "-",
                      tc_fmt(.data$estimate, digits)),
    ci = ifelse(is.na(.data$conf.low), "-",
                paste0("[", tc_fmt(.data$conf.low, digits), ", ",
                       tc_fmt(.data$conf.high, digits), "]")),
    p.value = ifelse(is.na(.data$p.value), "-",
                     format.pval(.data$p.value, digits = 3, eps = 0.001)),
    mod_p.value = if (has_mod) {
      ifelse(is.na(.data$mod_p.value), "-",
             format.pval(.data$mod_p.value, digits = 3, eps = 0.001))
    } else {
      NULL
    }
  )
  labs <- list(
    parameter = "Analysis Component",
    estimate = x$effect_measure,
    ci = paste0(format(100 * x$conf_level), "% CI"),
    p.value = "p-value"
  )
  if (has_cat) {
    disp$response_level[is.na(disp$response_level)] <- "-"
    labs$response_level <- "Outcome"
  }
  if (has_mod) {
    labs$mod_p.value <- glue::glue("Moderation p ({x$moderator})")
  }
  effects <- gt::gt(disp)
  effects <- do.call(gt::cols_label, c(list(effects), labs))
  align_cols <- c("parameter", if (has_cat) "response_level")
  effects <- gt::cols_align(effects, align = "left", columns = align_cols)
  effects <- gt::tab_header(
    effects,
    title = "Compositional effects",
    subtitle = "Each part versus the remaining composition"
  )
  if (isTRUE(x$per_part_slopes)) {
    effects <- gt::tab_source_note(
      effects,
      "Per-part effects are slopes at the reference composition (the model has non-linear composition terms)."
    )
  }
  if (!is.null(x$curvature_test)) {
    cv <- x$curvature_test
    effects <- gt::tab_source_note(
      effects,
      glue::glue(
        "Curvature (Wald): chi-squared({cv$df}) = {tc_fmt(cv$statistic, 3)}, p = {format.pval(cv$p.value, digits = 3)}"
      )
    )
  }
  if (!is.null(x$total_test)) {
    effects <- gt::tab_source_note(
      effects,
      glue::glue(
        "Total (log geometric mean, Wald): beta = {tc_fmt(x$total_beta, 3)}, p = {format.pval(x$total_test$p.value, digits = 3, eps = 0.001)}"
      )
    )
  }
  list(
    effects = effects |>
      gt::tab_options(
        heading.title.font.size    = gt::px(13),
        heading.subtitle.font.size = gt::px(11),
        heading.padding            = gt::px(2)
      ) |>
      gt::tab_style(
        style = gt::cell_text(weight = "bold"),
        locations = gt::cells_title(groups = "title")
      ) |>
      gt::tab_style(
        style = gt::cell_text(color = "#6b6b6b"),
        locations = gt::cells_title(groups = "subtitle")
      )
  )
}
