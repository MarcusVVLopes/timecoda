# Module: formula interface for comp_iv().
#
# Advanced users can specify the model as a `formula` instead of discrete
# arguments; it is fitted as written, so interactions, non-linear terms and
# random terms pass through verbatim to the underlying model function. The
# composition is referenced by the marker `ilr()`, which the parser below
# expands into the scalar pivot coordinates `ilr1 ... ilr{D-1}` (what emmeans
# `at = list(ilr1 = ...)` and the per-part rotations need): inside the formula
# operators it becomes the sum `(ilr1 + ... + ilr{D-1})`, and inside any other
# call it is distributed, one copy of the call per coordinate.

#' Composition marker for `comp_iv()` formulas
#'
#' A marker, valid only inside a `formula` passed to [comp_iv()], that stands
#' for the whole composition: the pivot isometric log-ratio (ILR) coordinates
#' `ilr1, ..., ilr{D-1}`. The formula parser expands it before the model is
#' fitted, and calling `ilr()` directly is an error.
#'
#' @details A [comp_iv()] formula is fitted as written once `ilr()` is
#'   expanded, by one of two rules that depend on where the marker sits. No
#'   composition term is refused; only `ilr()` with arguments is an error.
#'
#'   Inside the formula operators `+ - * : / ^ %in% | ||` and parentheses,
#'   `ilr()` is replaced by the sum `(ilr1 + ... + ilr{D-1})`, and R's formula
#'   algebra applies to it. `outcome ~ ilr() * sex` fits the composition, `sex`
#'   and every composition-by-`sex` interaction; `(ilr() | school)` is an
#'   `lme4` random slope on every coordinate. `ilr()^2` means what
#'   `(ilr1 + ilr2 + ilr3)^2` means in R: the main effects and every pairwise
#'   product of coordinates, with no squares.
#'
#'   Inside any other call, `ilr()` is distributed: the call is repeated once
#'   per coordinate and the copies are added. `I(ilr()^2)` becomes
#'   `(I(ilr1^2) + ... + I(ilr{D-1}^2))` and `poly(ilr(), 2)` becomes
#'   `(poly(ilr1, 2) + ... + poly(ilr{D-1}, 2))`; `splines::ns()`, `log()` and
#'   `scale()` are treated alike. Every `ilr()` within one such call takes the
#'   same coordinate, so `I(ilr() * ilr())` gives the squares; the pairwise
#'   products come from `ilr()^2`.
#'
#'   The full quadratic `ilr()^2 + I(ilr()^2)` (main effects, squares and
#'   pairwise products) is the same model in every ILR basis, so its fit and its
#'   reallocations do not depend on the order of the parts. The fit of a
#'   squares-only or `poly(ilr(), 2)` specification depends on the part order
#'   given to [declare_composition()]. On the `timeuse` data, the odds ratio of
#'   obesity for 10 minutes moved from SB to MVPA under `ilr() + I(ilr()^2)`
#'   ranged from 0.828 to 0.834 across three part orders. [comp_iv()] fits the
#'   full quadratic without a formula through `curvature = TRUE`.
#'
#'   [comp_iv()] warns when the formula leaves out the composition's main
#'   effect, that is, when a coordinate appears neither as a term of its own
#'   nor as the first argument of `poly()`, `ns()` or `bs()`.
#'
#' @return None; `ilr()` signals an error when called outside a formula.
#'
#' @examples
#' cd <- declare_composition(timeuse, c("MVPA", "LPA", "SB", "SPT"),
#'                           closure = 1440)
#' # The composition enters via ilr(); sex and its interaction are
#' # ordinary terms
#' comp_iv(cd, obesity ~ ilr() * sex, model = "glm", family = "binomial")
#' # A full quadratic in the composition, which does not depend on the
#' # part order
#' comp_iv(cd, zbmi ~ ilr()^2 + I(ilr()^2) + sex, model = "lm")
#'
#' @seealso [comp_iv()], whose `formula` argument accepts this marker.
#' @export
ilr <- function() {
  cli::cli_abort(c(
    "{.fn ilr} is only valid inside a {.fn comp_iv} formula.",
    "i" = "It marks where the composition enters the model and is expanded to {.code ilr1 + ... + ilr(D-1)} internally."
  ))
}

# Does an unevaluated expression contain a bare `ilr()` call anywhere?
tc_expr_has_ilr <- function(expr) {
  if (is.call(expr)) {
    if (is.name(expr[[1]]) && identical(expr[[1]], as.name("ilr"))) {
      return(TRUE)
    }
    return(any(purrr::map_lgl(as.list(expr)[-1], tc_expr_has_ilr)))
  }
  FALSE
}

# Expand each `ilr()` in a formula right-hand side: within the formula operators
# it becomes `coord_sum`; any other call holding it is copied once per
# coordinate, each copy's `ilr()` taking that coordinate, and the copies summed.
tc_expand_ilr_expr <- function(expr, coord_sum, call = rlang::caller_env()) {
  if (!is.call(expr)) {
    return(expr)
  }
  if (identical(expr[[1]], as.name("ilr"))) {
    if (length(expr) != 1L) {
      cli::cli_abort("{.fn ilr} takes no arguments inside a formula.", call = call)
    }
    return(coord_sum)
  }
  op <- if (is.name(expr[[1]])) as.character(expr[[1]]) else ""
  is_op <- op %in% c("+", "-", "*", ":", "/", "^", "%in%", "|", "||", "(")
  if (!is_op && !tc_expr_has_ilr(expr)) {
    return(expr)
  }
  subs <- if (is_op) list(coord_sum) else purrr::map(all.vars(coord_sum), as.name)
  copies <- rep(list(expr), length(subs))
  # Loops, not purrr::map(), so a nested abort keeps its own message.
  for (j in seq_along(subs)) {
    for (k in seq_along(expr)[-1]) {
      copies[[j]][[k]] <- tc_expand_ilr_expr(expr[[k]], subs[[j]], call)
    }
  }
  if (length(copies) == 1L) {
    return(copies[[1]])
  }
  call("(", purrr::reduce(copies, function(a, b) call("+", a, b)))
}


# The comp_iv() formula path and its moderation test --------------------------

# Parse a user (or `moderator`-synthesized) formula into comp_iv()'s internal
# spec: `ilr()` expanded, the left-hand side kept verbatim, covariate and
# random terms split out, and every column the formula names checked.
tc_iv_parse_formula <- function(formula, comp, model, call = rlang::caller_env()) {
  if (!rlang::is_formula(formula) || is.null(rlang::f_lhs(formula))) {
    cli::cli_abort("{.arg formula} must be a two-sided formula ({.code outcome ~ ...}).",
                   call = call)
  }
  coord_names <- tc_coord_names(length(comp$parts))
  lhs <- rlang::f_lhs(formula)
  rhs <- rlang::f_rhs(formula)
  # `ilr()` is expanded; without the marker anywhere, the composition is added
  # as a main effect, even beside coordinate names typed by hand.
  coord_sum <- str2lang(paste0("(", paste(coord_names, collapse = " + "), ")"))
  expanded <- tc_expand_ilr_expr(rhs, coord_sum, call)
  if (!tc_expr_has_ilr(rhs)) {
    expanded <- call("+", coord_sum, expanded)
  }
  # A name that is not a data column but a scalar in the formula's environment
  # (a knots value, `T`) is replaced by its value, so the stored formula, its
  # replays and the reallocation grid need nothing outside the data. A scalar
  # that shares a visible function's name, such as `c`, keeps its name here,
  # since substitute() would replace the function too.
  env <- rlang::f_env(formula)
  scalars <- tc_caller_scalars(
    setdiff(all.vars(call("~", lhs, expanded)),
            c(coord_names, names(comp$data_raw))),
    env
  )
  values <- scalars$values[setdiff(scalars$present, scalars$functions)]
  lhs <- do.call(substitute, list(lhs, values))
  expanded <- do.call(substitute, list(expanded, values))
  model_formula <- rlang::new_formula(lhs, expanded, env = env)

  # The `(1 | g)` random terms, and the fixed-effect terms without them.
  random <- purrr::map_chr(reformulas::findbars(model_formula),
                           function(b) paste0("(", deparse1(b), ")"))
  term_labels <- attr(
    stats::terms(reformulas::nobars(model_formula)), "term.labels"
  )

  # A coordinate's linear term is a term of its own or the first argument of a
  # poly(), ns() or bs() basis, which spans it; a missing one is fitted, with a
  # warning.
  linear <- purrr::map_chr(term_labels, function(label) {
    e <- str2lang(label)
    is_basis <- rlang::is_call(e, c("poly", "ns", "bs")) && length(e) > 1L
    if (is_basis) deparse1(e[[2]]) else label
  })
  if (!all(coord_names %in% linear)) {
    cli::cli_warn(
      c("The formula does not include the composition's main effect ({.code ilr()}).",
        "!" = "Estimates may be biased.",
        "i" = "Add {.code ilr()} to the formula unless this is intended."),
      call = call
    )
  }

  outcome <- NULL
  if (identical(tc_model_meta(model)$response, "surv")) {
    if (!is.call(lhs) || !grepl("Surv", deparse(lhs)[1])) {
      cli::cli_abort(
        "For {.code model = \"{model}\"} the formula left-hand side must be a {.fn survival::Surv} call.",
        call = call
      )
    }
  } else {
    if (!is.name(lhs)) {
      cli::cli_abort(
        "The formula left-hand side must be a single outcome column name.",
        call = call
      )
    }
    outcome <- as.character(lhs)
    if (!outcome %in% names(comp$data_raw)) {
      cli::cli_abort("Outcome {.val {outcome}} is not a column of the {.cls compD} data.",
                     call = call)
    }
  }

  # Every other name must be a data column, carried into the model frame;
  # anything else, a vector of the caller's included, is refused, so every
  # variable stays row-aligned with the compD.
  vars <- setdiff(all.vars(model_formula), coord_names)
  need <- intersect(vars, names(comp$data_raw))
  absent <- setdiff(vars, c(need, scalars$present))
  if (length(absent) > 0) {
    cli::cli_abort("Column{?s} {.val {absent}} not found in the {.cls compD} data.",
                   call = call)
  }
  list(outcome = outcome, time = NULL, event = NULL,
       covariates = setdiff(term_labels, coord_names), random = random,
       model_formula = model_formula, moderator = NULL, need = need)
}

# Of `nms`, those bound in `env` to a length-1 non-function (`present`), their
# values (`values`), and those of them also naming a visible function, such as
# `c` (`functions`).
tc_caller_scalars <- function(nms, env) {
  nms <- as.character(nms)
  objs <- purrr::map(rlang::set_names(nms), function(nm) get0(nm, envir = env))
  present <- nms[purrr::map_lgl(objs, function(x) {
    !is.function(x) && length(x) == 1L
  })]
  functions <- purrr::keep(present, function(nm) {
    exists(nm, envir = env, mode = "function")
  })
  list(present = present, values = objs[present], functions = functions)
}

# Joint Wald test of the composition-by-moderator coefficients (H0: no
# moderation): every `ilr{k}:moderator` term of the canonical fit, or, with
# `per_part`, the `ilr1:moderator` term(s) of one part's rotated fit.
tc_moderation_test <- function(model, model_type, moderator, parts,
                               per_part = FALSE) {
  coords <- if (per_part) "ilr1" else tc_coord_names(length(parts))
  what <- if (per_part) "part" else "composition"
  # The moderator's own coefficient labels, as model.matrix() builds them from
  # the fitting data the call carries: the name itself for a numeric moderator,
  # else the name followed by each contrast column.
  col <- stats::getCall(model)$data[[moderator]]
  labels <- if (is.null(col) || is.numeric(col)) {
    moderator
  } else {
    f <- droplevels(as.factor(col))
    # droplevels() drops a coding set on the column itself; restore it.
    if (!is.null(attr(col, "contrasts"))) stats::contrasts(f) <- attr(col, "contrasts")
    cm <- stats::contrasts(f)
    suffix <- colnames(cm)
    if (is.null(suffix)) suffix <- seq_len(ncol(cm))
    paste0(moderator, suffix)
  }
  pairs <- c(outer(coords, labels, paste, sep = ":"),
             outer(labels, coords, paste, sep = ":"))
  # tc_ilr_coef_names() also matches a multinomial fit's "<level>~" names; the
  # coefficients keep the model's own order.
  b_names <- names(tc_model_coefs(model, model_type)$b)
  int_names <- intersect(b_names, tc_ilr_coef_names(b_names, model_type, pairs))
  if (length(int_names) == 0) {
    cli::cli_abort(
      "Could not locate any {what}-by-{.val {moderator}} interaction coefficients."
    )
  }
  # intToUtf8(215) is the multiplication sign; kept out of source as ASCII.
  tc_joint_ilr_test(model, model_type, int_names,
                    term = paste(what, intToUtf8(215L), moderator))
}
