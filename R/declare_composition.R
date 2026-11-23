# Module 0: declaring a compositional variable.
#
# `declare_composition()` is the entry point of the composition-first workflow.
# It validates a set of parts, summarizes and imputes compositional zeroes,
# optionally flags compositional outliers, and returns a `compD` object
# carrying the closed composition, the zero-free parts on their original scale
# and a zero-summary table. `summary.compD()` and `plot.compD()` describe and
# draw this object.

#' Declare a compositional variable
#'
#' Builds a `compD` object from a set of composition parts. It validates the
#' parts (zeroes permitted), summarizes zero occurrence, imputes any zeroes,
#' optionally flags compositional outliers and stores the closed composition.
#'
#' Compositional data carry only relative information, and a zero has no
#' log-ratio. When zeroes are present, [comp_impute()] imputes them on the
#' original (pre-closure) scale before [compositions::acomp()] closes the
#' result. Every downstream description, plot and coordinate is computed from
#' this imputed, closed composition. The zero-free values on the original scale
#' keep each row's absolute total and are stored as `parts_imputed` for the
#' `include_total` T-space covariate of [comp_iv()]. Without zeroes, no
#' imputation is performed and `parts_imputed` equals the input parts.
#'
#' @section Outlier flagging:
#' With `detect_outliers = TRUE`, [comp_outliers()] runs on the zero-free parts
#' after imputation, and its per-observation columns are appended to `data_raw`
#' as `.outlier_distance`, `.outlier_limit` and `.outlier`. Rows are flagged and
#' never dropped. Every field of the object keeps one entry per input
#' observation, and `data_raw`, `composition` and `parts_imputed` stay
#' row-aligned with the input and with each other.
#'
#' Filtering is left to the caller. Dropping the flagged rows means declaring
#' the composition again, which recomputes it (and `closure = "mean"`) on the
#' remaining rows. The `.outlier` columns must be dropped in the same step.
#' `declare_composition()` writes them itself and stops with an error on any
#' input that already carries them, whether or not detection is requested.
#'
#' ```r
#' kept <- cd$data_raw |>
#'   dplyr::filter(!.outlier) |>
#'   dplyr::select(-dplyr::starts_with(".outlier"))
#'
#' cd_kept <- declare_composition(kept, parts, closure = 1440)
#' ```
#'
#' @param data A data frame containing the composition parts.
#' @param parts Character vector of column names making up the composition
#'   (at least two parts). Required; there is no default.
#' @param closure Scale the composition is reported on. Either a positive number
#'   (e.g., `1440` for minutes-per-day, `1` for proportions) or one of the
#'   strings `"percent"` (equivalent to `closure = 100`) or `"mean"` (the
#'   average of the observed row totals of the original input parts). The
#'   resolved numeric value is stored as `compD$closure`, and every downstream
#'   summary, plot and prediction uses it. Defaults to `1` (unit proportions).
#' @param labels Optional display labels for the parts, reported in summary and
#'   plot output in place of the raw part names. Either a named character vector
#'   mapping part names to labels (e.g., `c(MVPA = "Exercise (min/day)")`; parts
#'   not named keep their raw name) or an unnamed character vector of length
#'   `length(parts)` giving a label per part in order. Labels must be unique
#'   across parts (a label may not repeat another part's label or raw name).
#'   Stored on the object as a complete named vector keyed by every part; `NULL`
#'   (default) means the raw part names are used throughout.
#' @param units Optional unit label for the closure scale, a single string such
#'   as `"min/day"` for a 1440-minute day. The print methods show it verbatim
#'   after the closure and after values reported on that scale. It also appears
#'   in the [summary.compD()] and [summary.realloc()] publication tables
#'   ([pub_ready()]). `NULL` (default) shows no unit.
#' @param impute_method Zero-imputation method, passed to [comp_impute()]
#'   (ignored when the composition has no zeroes). Defaults to `"lrEM"`
#'   (log-ratio EM); see [comp_impute()] for the alternatives.
#' @param detect_outliers Single logical; flag compositional outliers
#'   ([comp_outliers()]) on the zero-free parts. Default `FALSE`.
#' @param outlier_args A named list of further arguments for [comp_outliers()],
#'   most usefully `quantile` (the chi-squared cutoff probability) or `alpha`
#'   (MCD's subset proportion). With `detect_outliers = FALSE` it is only
#'   checked to be a named list. With `TRUE`, a name that [comp_outliers()]
#'   does not accept is an error.
#' @param ... Additional arguments forwarded to [comp_impute()], and from there
#'   to the chosen method. The most useful is `dl`, the per-part detection limit
#'   below which an imputed value must fall. By default `declare_composition()`
#'   derives it from each part's smallest observed positive value on the
#'   original scale. Override it with `dl = <vector or matrix of limits>`, or
#'   pass other method-specific controls such as `lrEM()`'s `rob`, `tolerance`
#'   or `max.iter`.
#'
#' @return An object of class `compD`: a list with elements
#'   \describe{
#'     \item{`data_raw`}{The original input data (a tibble), unchanged, so it
#'       may still contain zeroes. With `detect_outliers = TRUE` it also carries
#'       three appended columns: `.outlier_distance` (the Mahalanobis distance
#'       in ILR coordinates), `.outlier_limit` (the cutoff, constant down the
#'       column) and `.outlier` (logical).}
#'     \item{`parts`}{The character vector of part names.}
#'     \item{`closure`}{The resolved numeric closure (always a positive number,
#'       even when the user supplied `"percent"` or `"mean"`).}
#'     \item{`composition`}{The imputed, closed composition, a
#'       [compositions::acomp()] object, from which [summary.compD()] and
#'       [plot.compD()] compute everything they show.}
#'     \item{`parts_imputed`}{The zero-free parts on the original
#'       (pre-closure) scale, a tibble: the input values with any zeroes
#'       replaced by [comp_impute()]. It keeps each row's absolute total,
#'       which the `include_total` (Coenders T-space) covariate of [comp_iv()]
#'       requires.}
#'     \item{`imput_table`}{A tibble describing zero occurrence per part:
#'       `part`, `n`, `n_zero`, `pct_zero`, `min_nonzero`.}
#'     \item{`labels`}{Either `NULL` or a complete named character vector keyed
#'       by every part, giving the display label used in summaries and plots.}
#'     \item{`units`}{The `units` argument as given, or `NULL`.}
#'     \item{`impute_method`}{The `impute_method` argument, or `NA_character_`
#'       when no part had a zero (the method was requested but never run).
#'       [summary.compD()] reports it as the method used.}
#'   }
#'
#' @references
#' Palarea-Albaladejo, J., & Martín-Fernández, J.A. (2015). zCompositions -- R
#' package for multivariate imputation of left-censored data under a
#' compositional approach. *Chemometrics and Intelligent Laboratory Systems*,
#' 143, 85-96. \doi{10.1016/j.chemolab.2015.02.019}
#'
#' @examples
#' cd <- declare_composition(
#'   timeuse,
#'   parts = c("MVPA", "LPA", "SB", "SPT"),
#'   closure = 1440,
#'   units = "min/day"
#' )
#' cd$imput_table
#'
#' set.seed(1)
#' flagged <- declare_composition(
#'   timeuse,
#'   parts = c("MVPA", "LPA", "SB", "SPT"),
#'   closure = 1440,
#'   detect_outliers = TRUE
#' )
#' sum(flagged$data_raw$.outlier)
#'
#' @seealso [summary.compD()] and [plot.compD()] to describe and visualize the
#'   object; [comp_impute()] and [comp_outliers()] for the zero-imputation and
#'   outlier-detection steps on their own; [comp_iv()] to model it as a
#'   predictor of an outcome.
#' @family declare and describe composition
#' @export
declare_composition <- function(data, parts, closure = 1, labels = NULL,
                                units = NULL, impute_method = "lrEM",
                                detect_outliers = FALSE,
                                outlier_args = list(), ...) {
  labels <- tc_normalize_labels(labels, parts)
  if (!is.null(units) &&
      !(is.character(units) && length(units) == 1 && !is.na(units) &&
        nzchar(trimws(units)))) {
    cli::cli_abort(
      "{.arg units} must be a single non-empty string (e.g. {.val min/day}) or {.code NULL}."
    )
  }

  # Zeroes (if any) are imputed on the ORIGINAL (pre-closure) scale, so the
  # zero-free raw values keep each row's absolute total, which the Coenders
  # T-space total covariate in comp_iv(include_total =) needs; only then is
  # the composition closed to `closure`. Closure is a per-row rescaling, so
  # the closed composition is unaffected by the scale imputation runs on.
  #
  # comp_impute() validates the parts, so this function does not repeat the
  # pass; `closure` is resolved afterwards for the same reason, since
  # `"mean"` reads `data[parts]` and must not do so before validation. The
  # error is identical either way, because check_composition() reports under
  # its own name from both call paths.
  imp <- comp_impute(data, parts, method = impute_method, ...)
  # Resolve `closure` to a positive numeric total.
  if (is.character(closure)) {
    closure <- rlang::arg_match(closure, c("percent", "mean"))
    closure <- if (closure == "percent") {
      100
    } else {
      mean(rowSums(as.matrix(data[parts])), na.rm = TRUE)
    }
  } else {
    tc_check_positive_scalar(
      closure, "closure",
      hint = "It may also be {.val percent} or {.val mean}."
    )
  }

  # Outlier flagging is a log-ratio distance, so it needs the zero-free parts:
  # it follows imputation. The distance is closure-invariant, so running it
  # before closure costs nothing. Rows are flagged and never removed here:
  # removing observations changes every estimate the composition feeds.
  data <- tc_declare_outliers(
    data, imp$parts_imputed, parts, detect_outliers, outlier_args
  )

  composition <- compositions::acomp(
    as.matrix(imp$parts_imputed[parts]), total = closure
  )

  structure(
    list(
      data_raw = tibble::as_tibble(data),
      parts = parts,
      closure = closure,
      composition = composition,
      parts_imputed = imp$parts_imputed,
      imput_table = imp$imput_table,
      labels = labels,
      units = units,
      impute_method = imp$impute_method
    ),
    class = "compD"
  )
}

#' Print a compositional variable
#'
#' Prints the parts, the closure (with its unit label, when one was declared)
#' and the number of observations of a [declare_composition()] object, with a
#' line each for imputed zeroes and flagged outliers when there are any.
#' [summary.compD()] gives the descriptive tables and [plot.compD()] the
#' figures.
#'
#' @param x A `compD` object from [declare_composition()].
#' @param ... Unused; present for S3 consistency.
#'
#' @return `x`, invisibly.
#'
#' @examples
#' declare_composition(timeuse, c("MVPA", "LPA", "SB", "SPT"), closure = 1440)
#'
#' @seealso [summary.compD()] for the descriptive tables and the fuller
#'   provenance report; [plot.compD()] for the figures.
#' @export
print.compD <- function(x, ...) {
  n <- nrow(x$data_raw)
  print(glue::glue("
<compD> compositional variable
  Parts: {glue::glue_collapse(tc_part_labels(x$parts, x$labels), sep = ', ')}
  Closure: {tc_with_units(x$closure, x$units)}
  Observations: {n}
"))

  # Both lines report what was found, so each is silent when there was
  # nothing, as in print.summary_compD().
  zeroed <- tc_zeroed_parts(x$imput_table)
  if (nrow(zeroed) > 0) {
    print(glue::glue(
      "  Zeroes imputed ({x$impute_method}): {sum(zeroed$n_zero)} across {glue::glue_collapse(zeroed$part, sep = ', ')}"
    ))
  }
  if (tc_has_outliers(x)) {
    n_flagged <- sum(x$data_raw$.outlier)
    pct <- tc_fmt(100 * n_flagged / n)
    print(glue::glue("  Outliers flagged: {n_flagged} ({pct}%)"))
  }
  invisible(x)
}

# Run declare-time outlier flagging and return `data` with comp_outliers()'s
# columns appended, or unchanged when flagging was not asked for.
tc_declare_outliers <- function(data, parts_imputed, parts, detect_outliers,
                                args, call = rlang::caller_env()) {
  # Checked rather than coerced with isTRUE(): silently treating a stray
  # `detect_outliers = 1` as "no" would skip flagging the caller asked for.
  if (!is.logical(detect_outliers) || length(detect_outliers) != 1 ||
      is.na(detect_outliers)) {
    cli::cli_abort(
      "{.arg detect_outliers} must be {.code TRUE} or {.code FALSE}.",
      call = call
    )
  }
  if (!is.list(args) || (length(args) > 0 && !rlang::is_named(args))) {
    cli::cli_abort(
      "{.arg outlier_args} must be a named list of {.fn comp_outliers} arguments.",
      call = call
    )
  }
  # Dot-prefixed after broom::augment()'s convention, and refused rather than
  # overwritten, as ilr_transform() does for its own coordinates. Checked
  # before the `detect_outliers` branch because these columns are the only
  # record that flagging ran, so a set left over from an earlier declaration
  # would otherwise be read downstream as this run's result, and reported with
  # a cutoff by a computation that never happened.
  cols <- c(".outlier_distance", ".outlier_limit", ".outlier")
  clash <- intersect(cols, names(data))
  if (length(clash) > 0) {
    cli::cli_abort(
      c(
        "{.arg data} already contains column{?s} {.val {clash}}.",
        "i" = "{cli::qty(length(clash))}{.fn declare_composition} writes {?this column/these columns} itself; drop {?it/them} first (e.g. with {.code dplyr::select(-dplyr::starts_with(\".outlier\"))})."
      ),
      call = call
    )
  }
  if (!detect_outliers) {
    return(data)
  }

  tbl <- do.call(comp_outliers, c(list(parts_imputed, parts), args))
  data[cols] <- tbl[c("distance", "limit", "is_outlier")]
  data
}

# Validate and normalize `labels` into a complete named vector keyed by every
# part; see declare_composition()'s @param labels for the accepted input shapes.
tc_normalize_labels <- function(labels, parts, call = rlang::caller_env()) {
  if (is.null(labels)) {
    return(NULL)
  }
  if (!is.character(labels) || anyNA(labels)) {
    cli::cli_abort(
      "{.arg labels} must be a character vector without missing values.",
      call = call
    )
  }
  out <- stats::setNames(parts, parts)
  nm <- names(labels)
  if (is.null(nm)) {
    if (length(labels) != length(parts)) {
      cli::cli_abort(
        c(
          "An unnamed {.arg labels} must have one entry per part.",
          "i" = "{length(parts)} part{?s} but {length(labels)} label{?s} supplied."
        ),
        call = call
      )
    }
    out[] <- labels
  } else {
    bad <- setdiff(nm, parts)
    if (length(bad) > 0) {
      cli::cli_abort(
        c(
          "{.arg labels} names {.val {bad}}, which {?is/are} not a part of the composition.",
          "i" = "Composition parts are {.val {parts}}."
        ),
        call = call
      )
    }
    out[nm] <- unname(labels)
  }
  # Plots key facets and legends by label, so two parts sharing one label
  # would silently merge into a single panel.
  dup <- unique(out[duplicated(out)])
  if (length(dup) > 0) {
    cli::cli_abort(
      c("{.arg labels} must be unique across parts.",
        "i" = "{.val {dup}} {?is/are} used for more than one part (a raw part name counts as that part's label)."),
      call = call
    )
  }
  out
}

# Look up display labels for a vector of part names, falling back to the raw
# name for any part absent from `labels` (or when `labels` is NULL).
tc_part_labels <- function(parts, labels) {
  if (is.null(labels)) {
    return(parts)
  }
  out <- unname(labels[parts])
  out[is.na(out)] <- parts[is.na(out)]
  out
}

# Format a value on the closure scale and append the compD's unit label
# ("1440 min/day"), or return it without a label when none was declared. A
# number is formatted here so every print line renders the closure alike;
# text (values the caller already rounded) is left as given.
tc_with_units <- function(value, units) {
  if (is.numeric(value)) value <- format(value)
  if (is.null(units)) value else paste(value, units)
}

# Did declare_composition() flag outliers? The appended columns are the only
# record, so this predicate is the one place that knows which column marks it.
tc_has_outliers <- function(object) {
  ".outlier" %in% names(object$data_raw)
}

# Fixed-decimal formatting for the text reports and the gt builders.
tc_fmt <- function(x, digits = 1) {
  formatC(x, format = "f", digits = digits)
}

# The rows of an `imput_table` for parts that actually had zeroes; both print
# methods start from this, then report it at their own level of detail.
tc_zeroed_parts <- function(imput_table) {
  imput_table[imput_table$n_zero > 0, , drop = FALSE]
}

# Reconstruct the parts data frame from the imputed composition vector, the
# single source of truth for every description (never `data_raw`).
tc_comp_df <- function(object) {
  df <- as.data.frame(unclass(object$composition))
  names(df) <- object$parts
  tibble::as_tibble(df)
}
