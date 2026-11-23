# Goodness-of-fit and residual diagnostics for comp_iv() fits.
# fit_check.comp_iv() picks a tc_check_*() backend from the registry's
# `check_backend` field; every backend returns list(tests, raw, note), which
# fit_check.comp_iv() turns into the returned object.

#' Check the goodness of fit of a composition model
#'
#' Runs residual and goodness-of-fit diagnostics on a fitted [comp_iv()] model,
#' with the diagnostic backend chosen from the model class. Returns a
#' `fit_check` list holding a tidy tibble of test statistics and the raw
#' backend object. With `plot = TRUE` it also draws the backend's diagnostic
#' plot and prints the report.
#'
#' @details The backend follows the model class:
#'   \describe{
#'     \item{`lm`, `glm`, `lmer`, `glmer`}{\pkg{DHARMa} simulation-based
#'       quantile residuals (KS uniformity, dispersion, and outlier tests).}
#'     \item{`cox`}{[survival::cox.zph()] (proportional-hazards test and
#'       Schoenfeld residual plots).}
#'     \item{ordinal `clm`}{\pkg{sure} surrogate residuals, drawn by
#'       \pkg{sure}'s own plot method as a Q-Q plot and a
#'       residuals-against-fitted plot, plus [ordinal::nominal_test()] (a
#'       likelihood-ratio test of proportional odds per term; Christensen,
#'       \pkg{ordinal}). The surrogate residuals are random draws made under a
#'       fixed seed (123) with the session's random state restored, so they
#'       repeat from call to call.}
#'     \item{robust `lmrob` / `rlmer`}{the robustness weights of the
#'       MM-estimate (Yohai 1987) or of the robust mixed fit (Koller 2016): the
#'       share of observations downweighted below 0.5 (half the full weight)
#'       and the smallest weight. A Shapiro-Wilk test of the residuals follows.}
#'     \item{mixed ordinal `clmm`, multinomial `mblogit`, mixed-effects Cox
#'       `coxme`}{fit statistics (AIC and log-likelihood) and a note only. No
#'       residual diagnostic is run for these classes, nor a proportional-odds
#'       test for `clmm` or a proportional-hazards test for `coxme`, and
#'       nothing is drawn. For `coxme` the log-likelihood is the penalized one
#'       that [stats::logLik()] returns, with the fit's effective degrees of
#'       freedom in the AIC. For an `mblogit` fit with a random effect, the
#'       convergence line shows \pkg{mclogit}'s own flag, which can read as
#'       converged although the inner iterations warned. The log-likelihood
#'       row is `NA`.}
#'   }
#'   The header of the printed report carries one line on the fit's convergence
#'   (`Convergence: converged` or `Convergence: NOT converged`) with the number
#'   of backend warnings captured at fit time. The line is left out when the
#'   backend reports no convergence state (`lm`, `cox` and `coxme` fits). The
#'   same state is returned in the `convergence` element and has no row in
#'   `tests`.
#'
#'   \pkg{DHARMa}, \pkg{sure} (with \pkg{tinyplot} for its plots) and
#'   \pkg{robustlmm} are optional (Suggests) and checked at call time.
#'
#' @param object A `comp_iv` object.
#' @param plot Logical; if `TRUE` (default), draw the backend's diagnostic
#'   plot(s) and print the report. The tidy results are returned either way.
#'   The `clmm`, `mblogit` and `coxme` fits draw nothing.
#' @param ... Passed to [DHARMa::simulateResiduals()] for the DHARMa backend
#'   (e.g., `n`, the number of simulations, or `seed`); ignored by the other
#'   backends. DHARMa simulates under its own seed (123 by default) and then
#'   restores the session's random state. The results are therefore the same
#'   in a fresh session and after any [set.seed()]; only the arguments passed
#'   here, such as `seed` or `n`, change them.
#'
#' @return Invisibly, a `fit_check` object: a list with
#'   \describe{
#'     \item{`tests`}{a tibble with one row per check and columns `check`,
#'       `statistic`, `df` and `p.value` (`NA` where a check has none).}
#'     \item{`raw`}{the backend object: the \pkg{DHARMa} simulation, the
#'       [survival::cox.zph()] result, a list of the surrogate residuals and
#'       the [ordinal::nominal_test()] table (`clm`), or the robustness
#'       weights (`lmrob`, `rlmer`). It is `NULL` for `clmm`, `mblogit` and
#'       `coxme` fits. The diagnostic plot is drawn as a side effect of
#'       `plot = TRUE` and is not returned; `plot(x$raw)` re-draws it for the
#'       \pkg{DHARMa} and `cox.zph` backends only.}
#'     \item{`note`}{an explanatory note, or `NULL`.}
#'     \item{`backend`}{the backend used: `"dharma"`, `"coxzph"`,
#'       `"ordinal"`, `"robust"`, `"robust_mixed"` or `"limited"`.}
#'     \item{`model_type`}{the model class.}
#'     \item{`convergence`}{the fit's convergence state as [comp_iv()] stores
#'       it: a one-row tibble with `converged` (`TRUE`, `FALSE`, or `NA` when
#'       the backend reports none), `singular` and `n_warnings` (the number of
#'       backend warnings captured at fit time).}
#'   }
#'
#' @references
#' Christensen, R. H. B. ordinal: Regression models for ordinal data. R
#'   package. <https://CRAN.R-project.org/package=ordinal>
#'   \doi{10.32614/CRAN.package.ordinal}
#'
#' Dunn, P. K., & Smyth, G. K. (1996). Randomized quantile residuals.
#'   *Journal of Computational and Graphical Statistics*, 5(3), 236-244.
#'   \doi{10.1080/10618600.1996.10474708}
#'
#' Grambsch, P. M., & Therneau, T. M. (1994). Proportional hazards tests and
#'   diagnostics based on weighted residuals. *Biometrika*, 81(3), 515-526.
#'   \doi{10.1093/biomet/81.3.515}
#'
#' Greenwell, B. M., McCarthy, A. J., Boehmke, B. C., & Liu, D. (2018).
#'   Residuals and diagnostics for binary and ordinal regression models: An
#'   introduction to the sure package. *The R Journal*, 10(1).
#'   \doi{10.32614/RJ-2018-004}
#'
#' Hartig, F. DHARMa: Residual diagnostics for hierarchical (multi-level /
#'   mixed) regression models. R package.
#'   <https://CRAN.R-project.org/package=DHARMa>
#'   \doi{10.32614/CRAN.package.DHARMa}
#'
#' Koller, M. (2016). robustlmm: An R package for robust estimation of linear
#'   mixed-effects models. *Journal of Statistical Software*, 75(6), 1-24.
#'   \doi{10.18637/jss.v075.i06}
#'
#' Liu, D., & Zhang, H. (2018). Residuals and diagnostics for ordinal
#'   regression models: A surrogate approach. *Journal of the American
#'   Statistical Association*, 113(522), 845-854.
#'   \doi{10.1080/01621459.2017.1292915}
#'
#' Yohai, V. J. (1987). High breakdown-point and high efficiency robust
#'   estimates for regression. *The Annals of Statistics*, 15(2), 642-656.
#'   \doi{10.1214/aos/1176350366}
#'
#' @examples
#' cd <- declare_composition(timeuse, c("MVPA", "LPA", "SB", "SPT"),
#'                           closure = 1440)
#' m <- comp_iv(cd, outcome = "obesity", model = "glm", family = "binomial")
#' if (requireNamespace("DHARMa", quietly = TRUE)) {
#'   fc <- fit_check(m, plot = FALSE)
#'   fc$tests
#' }
#'
#' @seealso [comp_iv()], the model object this diagnoses.
#' @export
fit_check <- function(object, plot = TRUE, ...) {
  if (!rlang::is_bool(plot)) {
    cli::cli_abort("{.arg plot} must be TRUE or FALSE, not {.obj_type_friendly {plot}}.")
  }
  UseMethod("fit_check")
}

#' @export
fit_check.default <- function(object, plot = TRUE, ...) {
  cli::cli_abort(
    "{.fn fit_check} works on {.cls comp_iv} objects, not {.cls {class(object)[1]}}."
  )
}

#' @rdname fit_check
#' @export
fit_check.comp_iv <- function(object, plot = TRUE, ...) {
  backend <- tc_model_meta(object$model_type)$check_backend
  res <- switch(backend,
    dharma = tc_check_dharma(object$model, plot = plot, ...),
    coxzph = tc_check_coxzph(object$model, plot = plot),
    ordinal = tc_check_ordinal(object$model, plot = plot),
    robust = tc_check_robust(object$model, plot = plot),
    robust_mixed = tc_check_robust_mixed(object$model, plot = plot),
    tc_check_limited(object)
  )
  # The report carries the fit's convergence tibble and is printed when the
  # diagnostic plot was requested.
  out <- structure(
    list(tests = res$tests, raw = res$raw, note = res$note,
         backend = backend, model_type = object$model_type,
         convergence = object$convergence),
    class = "fit_check"
  )
  if (isTRUE(plot)) print(out)
  invisible(out)
}

# Print ------------------------------------------------------------------------

#' @export
print.fit_check <- function(x, ...) {
  label <- tc_report_meta(x$model_type)$label
  # The header, backend and note lines each end in a space; `{' '}` keeps it
  # where the line would otherwise end in whitespace.
  print(glue::glue("
<fit_check> {label}{' '}
  Backend: {x$backend}{' '}
"))
  # Only where the backend reports a state (lm, cox and coxme report none).
  status <- x$convergence
  if (!is.na(status$converged)) {
    state <- if (isTRUE(status$converged)) "converged" else "NOT converged"
    n <- status$n_warnings
    s <- if (n == 1) "" else "s"
    print(glue::glue("  Convergence: {state} ({n} backend warning{s})"))
  }
  if (!is.null(x$note)) print(glue::glue("  Note: {x$note} "))
  if (nrow(x$tests) > 0) {
    cat("\n")
    print(x$tests)
  }
  invisible(x)
}

# Evaluate `expr` under a fixed seed and put the session's random state back,
# so a randomized diagnostic gives the same result in a fresh session and after
# any set.seed(). The state is created first when the session has none.
tc_with_seed <- function(expr, seed = 123) {
  if (!exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
    stats::runif(1)
  }
  old_seed <- get(".Random.seed", envir = globalenv())
  on.exit(assign(".Random.seed", old_seed, envir = globalenv()), add = TRUE)
  set.seed(seed)
  expr
}

# A tidy one-check-per-row tibble with the package's fixed columns.
tc_check_tibble <- function(check, statistic = NA_real_, df = NA_real_,
                            p.value = NA_real_) {
  tibble::tibble(
    check = check,
    statistic = as.numeric(statistic),
    df = as.numeric(df),
    p.value = as.numeric(p.value)
  )
}

# Backends ---------------------------------------------------------------------

# DHARMa simulation-based quantile residuals: works for lm/glm/lmer/glmer.
tc_check_dharma <- function(model, plot = TRUE, n = 250, ...) {
  rlang::check_installed("DHARMa", reason = "for residual diagnostics of this model.")
  # One seeded block holds the simulation and the three tests: DHARMa seeds
  # its own simulation (123 unless `seed` is passed) and puts the state back,
  # so the bootstrapped outlier test would otherwise draw from whatever the
  # session's stream was. tc_with_seed() also makes sure a state exists first,
  # since DHARMa draws once after seeding when the session has none.
  res <- tc_with_seed({
    sr <- DHARMa::simulateResiduals(model, n = n, plot = FALSE, ...)
    list(
      sr = sr,
      u = DHARMa::testUniformity(sr, plot = FALSE),
      d = DHARMa::testDispersion(sr, plot = FALSE),
      o = DHARMa::testOutliers(sr, plot = FALSE)
    )
  })
  if (isTRUE(plot)) graphics::plot(res$sr)
  tests <- tc_check_tibble(
    check = c("KS uniformity", "Dispersion", "Outliers"),
    statistic = c(res$u$statistic, res$d$statistic, res$o$statistic),
    p.value = c(res$u$p.value, res$d$p.value, res$o$p.value)
  )
  list(tests = tests, raw = res$sr, note = NULL)
}

# survival::cox.zph proportional-hazards test + Schoenfeld residual plots.
tc_check_coxzph <- function(model, plot = TRUE) {
  z <- survival::cox.zph(model)
  if (isTRUE(plot)) graphics::plot(z)
  tab <- z$table
  tests <- tc_check_tibble(
    check = paste0("Proportional hazards: ", rownames(tab)),
    statistic = tab[, "chisq"], df = tab[, "df"], p.value = tab[, "p"]
  )
  list(tests = tests, raw = z, note = NULL)
}

# Ordinal clm: sure surrogate-residual plot + ordinal::nominal_test (the
# proportional-odds assumption). nominal_test's first row is the baseline model
# (no LR test), so only the per-term rows carry a p-value.
tc_check_ordinal <- function(model, plot = TRUE) {
  rlang::check_installed("sure", version = "0.3.0",
                         reason = "for ordinal surrogate residuals.")
  # Surrogate residuals are random draws: seeded, so a call is reproducible.
  sres <- tc_with_seed(sure::resids(model))
  if (isTRUE(plot)) {
    # sure's own plot() method (a Q-Q plot against the distribution the link
    # implies, and residuals against the fitted mean) draws with tinyplot.
    rlang::check_installed("tinyplot",
                           reason = "to draw the surrogate-residual plots.")
    graphics::plot(sres, what = c("qq", "fitted"), fit = model)
  }
  # Proportional-odds LR test, refitting once per term through the stored call.
  nt <- tryCatch(ordinal::nominal_test(model), error = function(e) NULL)
  po_rows <- tc_check_tibble(character(0))
  if (!is.null(nt)) {
    p_col <- grep("^Pr", names(nt), value = TRUE)[1]
    keep <- which(!is.na(nt[[p_col]]))
    if (length(keep) > 0) {
      po_rows <- tc_check_tibble(
        check = paste0("Proportional odds: ", rownames(nt)[keep]),
        statistic = nt[["LRT"]][keep], df = nt[["Df"]][keep],
        p.value = nt[[p_col]][keep]
      )
    }
  }
  note <- if (is.null(nt)) {
    "The proportional-odds test could not be run for this fit; re-run with `plot = TRUE` to inspect the surrogate-residual plot."
  } else if (nrow(po_rows) == 0) {
    "The proportional-odds test returned no estimable term; re-run with `plot = TRUE` to inspect the surrogate-residual plot."
  } else {
    "Proportional-odds rows: likelihood-ratio test of each term against a model letting its effect vary across the outcome thresholds; a small p-value is evidence that the proportional-odds assumption fails for that term."
  }
  fit_rows <- tc_check_tibble(
    check = c("AIC", "log-likelihood"),
    statistic = c(stats::AIC(model), as.numeric(stats::logLik(model)))
  )
  tests <- dplyr::bind_rows(po_rows, fit_rows)
  list(tests = tests, raw = list(surrogate = sres, nominal_test = nt),
       note = note)
}

# Robust lmrob: summarize the robustness weights (which observations the robust
# fit downweighted) plus a residual-normality check.
tc_check_robust <- function(model, plot = TRUE) {
  w <- model$rweights
  if (isTRUE(plot)) {
    tryCatch(graphics::plot(model, which = c(2, 3), ask = FALSE),
             error = function(e) NULL)
  }
  tc_check_robust_weights(w, stats::residuals(model))
}

# Robust rlmer: residual robustness weights (w_e) from robustlmm.
tc_check_robust_mixed <- function(model, plot = TRUE) {
  rlang::check_installed("robustlmm", reason = "for robust mixed-model diagnostics.")
  w <- robustlmm::getME(model, "w_e")
  if (isTRUE(plot)) {
    tryCatch(print(graphics::plot(model)), error = function(e) NULL)
  }
  tc_check_robust_weights(w, stats::residuals(model))
}

# Weight summary and residual Shapiro-Wilk test shared by both robust backends;
# the Shapiro-Wilk row is omitted when shapiro.test() refuses the residuals
# (fewer than 3 or more than 5000 of them, or all identical).
tc_check_robust_weights <- function(w, resid) {
  sw <- tryCatch(stats::shapiro.test(resid), error = function(e) NULL)
  tests <- tc_check_tibble(
    check = c("Downweighted obs (weight < 0.5)", "Minimum robustness weight"),
    statistic = c(mean(w < 0.5), min(w))
  )
  if (!is.null(sw)) {
    tests <- dplyr::bind_rows(
      tests,
      tc_check_tibble("Residual normality (Shapiro-Wilk)",
                      statistic = sw$statistic, p.value = sw$p.value)
    )
  }
  list(tests = tests, raw = w, note = NULL)
}

# Any registry row whose `check_backend` is "limited" (the default): AIC and
# log-likelihood rows plus the row's `limited_note`, or a generic one.
tc_check_limited <- function(object) {
  model <- object$model
  tests <- tc_check_tibble(
    check = c("AIC", "log-likelihood"),
    statistic = c(stats::AIC(model), as.numeric(stats::logLik(model)))
  )
  note <- tc_model_meta(object$model_type)$limited_note
  if (is.null(note)) {
    note <- "Residual diagnostics are not available for this model type; only fit statistics are reported."
  }
  list(tests = tests, raw = NULL, note = note)
}
