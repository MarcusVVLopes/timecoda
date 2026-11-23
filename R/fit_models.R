# Internal model machinery behind comp_iv(): the print method for its fits, the
# tc_iv_fit() builder, and the helpers that read a fitted model. Not exported.

#' Print a timecoda model fit
#'
#' Prints a compact summary of a `timecoda_fit` object, followed by the
#' underlying fitted model object. The summary gives the model class and (for
#' GLM-family engines) its family, the composition parts and total (with its
#' unit label, when one was declared) and the reference composition (the
#' sample's closed geometric mean, from which reallocations start).
#'
#' @param x A `timecoda_fit` object (or subclass, e.g., `comp_iv`), as returned
#'   by [comp_iv()].
#' @param ... Further arguments (unused; present for S3 consistency).
#'
#' @return `x`, invisibly.
#'
#' @examples
#' cd <- declare_composition(timeuse, c("MVPA", "LPA", "SB", "SPT"),
#'                           closure = 1440)
#' m <- comp_iv(cd, outcome = "obesity", model = "glm", family = "binomial")
#' print(m)
#'
#' @seealso [comp_iv()] to fit the model; [summary.comp_iv()] for the full
#'   effect table.
#' @export
print.timecoda_fit <- function(x, ...) {
  ref <- round(unlist(x$ref_comp), 1)
  print(glue::glue("
<timecoda_fit> compositional {tc_report_meta(x$model_type, x$family)$label}
  Parts: {glue::glue_collapse(x$parts, sep = ', ')} (total = {tc_with_units(x$total, x$comp$units)})
  Reference composition: {glue::glue_collapse(glue::glue('{names(ref)} = {ref}'), sep = ', ')}
"))
  n_warn <- nrow(x$fit_warnings)
  if (n_warn > 0) {
    print(glue::glue(
      "  Convergence: {n_warn} captured backend message{if (n_warn == 1) '' else 's'} (see summary())"
    ))
  }
  cat("\n")
  # The backend call holds the fitting frame itself; name it, so the printed
  # call stays one line instead of the whole data set.
  shown <- x$model
  if (isS4(shown)) {
    shown@call$data <- quote(data)
  } else {
    shown$call$data <- quote(data)
  }
  print(shown)
  invisible(x)
}

# The model builder ----------------------------------------------------------

# Fit one comp_iv() model at one pivot: build the formula, gather the backend's
# arguments, fit while capturing backend conditions, and bundle the result as a
# `timecoda_fit`. How each engine is called is read from its registry row
# (`fitter`, `fixed_args`, `response`, `random_style`), so nothing here is
# per-engine.
#
# comp_iv() calls this once per pivot, so everything pivot-invariant (the
# argument rules, the response checks, the reference composition) is done
# there, once, before the first call. `data` already carries this pivot's ilr
# columns; `source` labels the captured conditions ("canonical", or the rotated
# part's name). A `spec$model_formula` (formula/moderator interface) is used as
# given. `family` is the resolved family object; `family_expr` is the family
# as the user typed it, which the stored call records.
tc_iv_fit <- function(data, model, spec, parts, total, ref_comp, family = NULL,
                      family_expr = family, pivot = 1L, model_args = list(),
                      source = "canonical", env = rlang::caller_env(),
                      call = rlang::caller_env()) {
  meta <- tc_model_meta(model)
  formula <- spec$model_formula
  # An engine that takes its grouping through its own `random =` argument
  # (multinomial) must not also see the bar term in the formula.
  if (!is.null(formula) && identical(meta$random_style, "arg")) {
    frm_env <- environment(formula)
    formula <- reformulas::nobars(formula)
    environment(formula) <- frm_env
  }

  # Formula: outcome ~ ilr1 + ... + ilr{D-1} + covariates, with random-effect
  # bar terms folded in for the engines that carry them in the formula.
  if (is.null(formula)) {
    response <- if (identical(meta$response, "surv")) {
      paste0("survival::Surv(", spec$time, ", ", spec$event, ")")
    } else {
      spec$outcome
    }
    extra <- spec$covariates
    if (identical(meta$random_style, "bars") && length(spec$random) > 0) {
      extra <- c(extra, tc_random_terms(spec$random, data, call = call))
    }
    # The formula gets a small environment of its own, parented on the package
    # namespace, rather than the caller's frame: the fitted object keeps the
    # formula's environment, so a fit made inside a function would otherwise
    # carry (and save) everything in that function.
    formula <- stats::reformulate(
      c(tc_coord_names(length(parts)), extra),
      response = str2lang(response),
      env = rlang::new_environment(parent = topenv(environment()))
    )
  }
  # model.frame() looks up `weights`, `offset` and `subset` in the data and then
  # the formula's environment, so a `model_args` value naming an object local
  # to the caller (`weights = w`) is copied into an environment timecoda made
  # (the discrete path's, or the `moderator` formula's, marked `.tc_synth`),
  # never into a user's own formula's; data columns are left to the model frame.
  formula_env <- environment(formula)
  if (is.null(spec$model_formula) || isTRUE(formula_env$.tc_synth)) {
    free <- setdiff(unlist(purrr::map(model_args, all.vars)), names(data))
    purrr::walk(free, function(nm) {
      if (exists(nm, envir = env)) {
        assign(nm, get(nm, envir = env), envir = formula_env)
      }
    })
  }

  # Backend arguments the engine always passes, plus the two it may: the GLM
  # family comp_iv() resolved, and the multinomial engine's own
  # `random = ~ <lhs> | g`, from a plain grouping name (`~ 1 | g`) or a bar
  # term such as "(ilr1 | school)", whose left-hand side is kept.
  args <- meta$fixed_args
  if (isTRUE(meta$allows_family)) {
    args$family <- family
  }
  if (identical(meta$random_style, "arg") && length(spec$random) > 0) {
    # comp_iv() has already required exactly one `|` term for this style; a
    # plain grouping name is checked against the data and comes back as
    # "(1 | name)".
    term <- str2lang(tc_random_terms(spec$random, data, call = call))
    while (is.call(term) && identical(term[[1]], quote(`(`))) term <- term[[2]]
    if (!(is.call(term) && identical(term[[1]], quote(`|`)))) {
      cli::cli_abort(
        "{.arg random} for {.code model = \"{model}\"} must be a grouping name or one bar term; got {.val {spec$random}}.",
        call = call
      )
    }
    args$random <- stats::as.formula(call("~", term))
  }

  # `fitter` is spliced in quoted, so the fitted object's own `call` still names
  # the backend, and the formula as itself. The call is evaluated in a child of
  # the user's environment (`env`, the frame comp_iv() was called from), so a
  # `model_args` value named there (`control = ctrl`) resolves as in a direct
  # backend call; the frame is reached there as `prep$data`.
  fit_call <- rlang::expr(
    (!!meta$fitter)(!!formula, data = prep$data, !!!args, !!!model_args)
  )
  prep <- list(data = data)
  fit_env <- rlang::env(env, prep = prep)
  # Fit while capturing (but not muffling) the backend's warnings and messages,
  # so convergence trouble can be surfaced later; the user still sees them.
  captured <- list()
  record <- function(type, msg) {
    captured[[length(captured) + 1L]] <<- list(
      source = source, type = type, message = msg
    )
  }
  # A row with `quiet` discards what its backend prints to the console during
  # the fit; capture.output() lets the value, warnings, messages and errors
  # pass as they would without it.
  fitted <- withCallingHandlers(
    if (meta$quiet) {
      utils::capture.output(value <- eval(fit_call, envir = fit_env))
      value
    } else {
      eval(fit_call, envir = fit_env)
    },
    warning = function(w) record("warning", conditionMessage(w)),
    message = function(m) {
      msg <- sub("\n$", "", conditionMessage(m))
      # Skip trivial package-startup chatter; keep everything substantive.
      if (!grepl("^(Loading|Attaching|Registered S3|The following)", msg)) {
        record("message", msg)
      }
    }
  )
  fit_warnings <- tibble::tibble(
    source = purrr::map_chr(captured, "source"),
    type = purrr::map_chr(captured, "type"),
    message = purrr::map_chr(captured, "message")
  )

  # Convergence, by the registry's `convergence` style: "lme4" reads
  # @optinfo's messages and isSingular(), "flag" the fit's $converged,
  # "ordinal" its $convergence code (clm) or $optRes$convergence (clmm, 0 =
  # converged), "optinfo" robustlmm's optimizer code; "none" (lm, cox, coxme)
  # is NA.
  style <- meta$convergence
  converged <- NA
  singular <- NA
  if (style == "lme4") {
    msgs <- tryCatch(fitted@optinfo$conv$lme4$messages, error = function(e) NULL)
    converged <- length(msgs) == 0
    singular <- tryCatch(lme4::isSingular(fitted), error = function(e) NA)
  } else if (style == "flag") {
    conv <- tryCatch(fitted$converged, error = function(e) NULL)
    if (!is.null(conv)) converged <- isTRUE(conv)
  } else if (style == "ordinal") {
    conv <- tryCatch(fitted$convergence, error = function(e) NULL)
    if (is.null(conv)) conv <- fitted$optRes$convergence
    if (is.list(conv) && !is.null(conv$code)) {
      converged <- isTRUE(conv$code == 0)
    } else if (is.numeric(conv)) {
      converged <- isTRUE(conv == 0)
    }
  } else if (style == "optinfo") {
    code <- tryCatch(fitted@optinfo$conv$opt, error = function(e) NULL)
    if (!is.null(code)) converged <- isTRUE(code == 0)
  }

  # The stored call then holds the fitting frame itself in place of
  # `prep$data`, which exists only in `fit_env`: with the formula and the data
  # both held as values, the call replays anywhere, so update(), drop1(),
  # anova() between fits, survival::cox.zph() on a plain Cox fit (current
  # survival does not support it on a coxme fit) and emmeans' Satterthwaite
  # degrees of freedom (which re-fit through it) all work on the returned
  # backend object, and emmeans called on it directly recovers the frame.
  # emmeans reads the formula's `survival::Surv()` response as a
  # transformation; predict_realloc() removes it. The family, which the backend
  # received as the resolved object, is recorded as the user typed it, as a
  # direct call would.
  call_fields <- list(data = data)
  if (isTRUE(meta$allows_family) && !is.null(stats::getCall(fitted)$family)) {
    call_fields$family <- family_expr
  }
  stored_call <- rlang::call_modify(stats::getCall(fitted), !!!call_fields)
  if (isS4(fitted)) {
    fitted@call <- stored_call
  } else {
    fitted$call <- stored_call
  }
  covariates <- spec$covariates
  structure(
    list(
      model = fitted,
      model_type = model,
      parts = parts,
      total = total,
      pivot = pivot,
      ref_comp = ref_comp,
      covariates = covariates,
      # Composition-by-covariate interactions make reallocation contrasts
      # covariate-conditional; products of coordinates alone do not.
      has_ilr_interaction = length(
        tc_composition_terms(covariates, tc_coord_names(length(parts)))$crossed
      ) > 0,
      family = if (isTRUE(meta$allows_family)) family else NULL,
      data = data,
      fit_warnings = fit_warnings,
      convergence = tibble::tibble(
        converged = converged,
        singular = singular,
        n_warnings = as.integer(nrow(fit_warnings))
      )
    ),
    class = "timecoda_fit"
  )
}

# A `random` spec (grouping column names and/or full bar terms) as "(1 | name)"
# term strings. An element already containing "|" is used verbatim; every plain
# name must be a column of the data when one is given, and the missing ones are
# named together. Shared by tc_iv_fit() and the moderator sugar in
# tc_iv_resolve_spec(), which builds its formula as text before any model frame
# exists.
tc_random_terms <- function(random, data = NULL, call = rlang::caller_env()) {
  if (!is.character(random) || length(random) < 1 || anyNA(random)) {
    cli::cli_abort(
      "{.arg random} must be a character vector of grouping column names and/or random-effects terms.",
      call = call
    )
  }
  is_bar <- grepl("|", random, fixed = TRUE)
  missing <- if (is.null(data)) character(0) else setdiff(random[!is_bar], names(data))
  if (length(missing) > 0) {
    cli::cli_abort(
      "Grouping column{?s} {.val {missing}} not found in the {.cls compD} data.",
      call = call
    )
  }
  ifelse(is_bar, random, paste0("(1 | ", random, ")"))
}

# A GLM family, given as an object, a function or a name, as a family object.
# `env` is where a name is looked up: comp_iv() passes its own caller's
# environment, so a locally defined family function resolves.
tc_match_family <- function(family, env = rlang::caller_env(),
                            call = rlang::caller_env()) {
  if (is.null(family)) {
    cli::cli_abort(
      "{.arg family} is required (e.g. {.val binomial}).",
      call = call
    )
  }
  # Already a family object, e.g. `gaussian()`.
  if (inherits(family, "family")) {
    return(family)
  }
  if (is.character(family)) {
    family <- tryCatch(
      get(family, mode = "function", envir = env)(),
      error = function(e) NULL
    )
  }
  if (is.function(family)) {
    family <- family()
  }
  if (!inherits(family, "family")) {
    cli::cli_abort(
      "{.arg family} must be a GLM family object or name (e.g. {.val binomial}).",
      call = call
    )
  }
  family
}

# Reading a fitted model ------------------------------------------------------
#
# Everything below reads the fitted object tc_iv_fit() produced. It lives beside
# the builder rather than in model_backends.R (which holds only the registry)
# because every caller is part of fitting a comp_iv() model: the omnibus test in
# comp_iv(), the T-space total test, and the two moderation tests.

# Fixed-effect estimates and their covariance, extracted uniformly across model
# classes. Mixed fits (lmer/glmer/rlmer) use fixef(); everything else coef().
tc_model_coefs <- function(model, model_type) {
  b <- if (isTRUE(tc_model_meta(model_type)$fixef)) {
    lme4::fixef(model)
  } else {
    stats::coef(model)
  }
  list(b = b, V = as.matrix(stats::vcov(model)))
}

# The term labels that involve the composition: those whose variables are
# coordinates `only`, and those that cross a coordinate with another variable.
tc_composition_terms <- function(labels, coord_names) {
  vars <- purrr::map(labels, function(label) all.vars(str2lang(label)))
  has_coord <- purrr::map_lgl(vars, function(v) any(v %in% coord_names))
  has_other <- purrr::map_lgl(vars, function(v) !all(v %in% coord_names))
  list(only = labels[has_coord & !has_other],
       crossed = labels[has_coord & has_other])
}

# The data columns named in a fit's terms that cross a coordinate with another variable.
tc_interacting_covariates <- function(object) {
  coords <- tc_coord_names(length(object$parts))
  purrr::map(tc_composition_terms(object$covariates, coords)$crossed,
             function(term) setdiff(all.vars(str2lang(term)), coords)) |>
    purrr::list_c() |>
    unique() |>
    intersect(names(object$data))
}

# The coefficient names the term `labels` generate: each label itself, a call
# component (poly(), ns()) optionally followed by its basis column ("1", "1.0"),
# and, with `levels`, any non-coordinate component by its factor level ("sex"
# -> "sexmale", "factor(sex)" -> "factor(sex)male"); a coordinate matches
# exactly. A multinomial fit's "<level>~" prefix is matched past, so every
# category's coefficient is returned.
tc_ilr_coef_names <- function(b_names, model_type, labels, levels = FALSE) {
  predictor <- if (isTRUE(tc_model_meta(model_type)$per_category)) {
    sub("^[^~]*~", "", b_names)
  } else {
    b_names
  }
  # Split on `:` but not on the `::` of a namespaced call.
  components <- strsplit(labels, "(?<!:):(?!:)", perl = TRUE)
  patterns <- purrr::map_chr(components, function(cmp) {
    suffix <- ifelse(levels & !grepl("^ilr[0-9]+$", cmp), "[^:]*",
                     ifelse(endsWith(cmp, ")"), "[0-9.]*", ""))
    paste0("^", paste0("\\Q", cmp, "\\E", suffix, collapse = ":"), "$")
  })
  hit <- purrr::map(patterns, function(p) grepl(p, predictor, perl = TRUE))
  b_names[purrr::reduce(hit, `|`, .init = rep(FALSE, length(b_names)))]
}

# Joint (multi-df) Wald test that a named set of model coefficients are jointly
# zero. Used for the composition omnibus (every composition term), the
# curvature test, the T-space total test, and both moderation tests, which
# differ only in which coefficients they name. Works uniformly across every
# engine because it reads the fixed effects and their covariance through
# tc_model_coefs().
tc_joint_ilr_test <- function(model, model_type, coef_names,
                              term = "composition") {
  cf <- tc_model_coefs(model, model_type)
  b <- cf$b
  V <- cf$V
  if (!all(coef_names %in% names(b)) ||
      !all(coef_names %in% rownames(V))) {
    cli::cli_abort("Could not locate the {term} coefficients in the model.")
  }
  # An aliased (NA) coefficient repeats a column already tested and is left
  # out; with none left there is nothing to test.
  coef_names <- coef_names[!is.na(b[coef_names])]
  if (length(coef_names) == 0) {
    return(tibble::tibble(term = term, statistic = NA_real_, df = 0L,
                          p.value = NA_real_))
  }
  bb <- b[coef_names]
  VV <- V[coef_names, coef_names, drop = FALSE]
  stat <- as.numeric(t(bb) %*% solve(VV) %*% bb)
  df <- length(coef_names)
  tibble::tibble(
    term = term,
    statistic = stat,
    df = df,
    p.value = stats::pchisq(stat, df = df, lower.tail = FALSE)
  )
}
