# Module 0: compositional outlier detection.
#
# comp_outliers() wraps robCompositions::outCoDa(), the Filzmoser & Hron
# (2008) robust Mahalanobis method, behind one interface. Detection runs in
# pivot ILR coordinates, so it needs a zero-free composition:
# declare_composition() calls this after comp_impute().

#' Detect compositional outliers
#'
#' Flags observations whose *ratios between parts* are atypical, using the
#' robust Mahalanobis distance of Filzmoser and Hron (2008) computed in pivot
#' isometric log-ratio coordinates. An observation is an outlier when its
#' distance exceeds `sqrt(qchisq(quantile, df = D - 1))` for a `D`-part
#' composition. The distance is invariant to the closure of the composition
#' and to the choice of ILR pivot.
#'
#' @details
#' The centre and covariance come from the minimum covariance determinant
#' (MCD) estimator, which resists the outliers it is meant to find. Its
#' classical counterpart, the arithmetic mean and sample covariance, is itself
#' distorted by them (the masking effect) and is not offered.
#'
#' MCD estimates the covariance from a randomized subset search, so distances
#' vary slightly between calls. A call to [set.seed()] beforehand makes the
#' flags reproducible.
#'
#' Detection requires a zero-free composition. Zeroes can be imputed first with
#' [comp_impute()], or the composition built with [declare_composition()], which
#' imputes before it flags.
#'
#' @param data A data frame containing the composition parts, free of zeroes.
#' @param parts Character vector of column names making up the composition
#'   (at least two parts).
#' @param quantile Chi-squared probability setting the outlier cutoff.
#'   Defaults to `0.975`.
#' @param ... Additional arguments forwarded to [robCompositions::outCoDa()].
#'   The useful one is `alpha`, the proportion of observations MCD retains in
#'   its subset search (default `0.5`, the most robust; higher values trade
#'   robustness for efficiency). The `coda` argument of `outCoDa()` cannot be
#'   overridden, and passing it is an error.
#'
#' @return A tibble with one row per observation and the columns `.row`
#'   (position in `data`), `distance` (robust Mahalanobis distance in ILR
#'   coordinates), `limit` (the cutoff, constant down the column), and
#'   `is_outlier` (logical, `distance > limit`). Rows are never removed;
#'   acting on the flags is left to the caller.
#'
#' @references
#' Filzmoser, P., & Hron, K. (2008). Outlier detection for compositional data
#' using robust methods. *Mathematical Geosciences*, 40, 233-248.
#' \doi{10.1007/s11004-007-9141-5}
#'
#' Rousseeuw, P.J., & Van Driessen, K. (1999). A fast algorithm for the minimum
#' covariance determinant estimator. *Technometrics*, 41, 212-223.
#' \doi{10.1080/00401706.1999.10485670}
#'
#' @examples
#' set.seed(1)
#' out <- comp_outliers(timeuse, c("MVPA", "LPA", "SB", "SPT"))
#' sum(out$is_outlier)
#'
#' @seealso [declare_composition()], which calls this when
#'   `detect_outliers = TRUE`; [comp_impute()] to remove zeroes first.
#' @family declare and describe composition
#' @export
comp_outliers <- function(data, parts, quantile = 0.975, ...) {
  check_composition(data, parts, allow_zero = TRUE)
  # Checked here rather than through check_composition()'s strict default so
  # the message can name the fix, as this is the one likely way to reach it.
  zeroed <- parts[purrr::map_lgl(data[parts], ~ any(.x == 0))]
  if (length(zeroed) > 0) {
    cli::cli_abort(
      c(
        "Composition part{?s} {.val {zeroed}} contain{?s/} zeroes.",
        "i" = "A log-ratio of zero is undefined; impute them with {.fn comp_impute}, or build the composition with {.fn declare_composition}, first."
      )
    )
  }

  n <- nrow(data)
  n_coord <- length(parts) - 1L
  # MCD needs a subset of about (n + p + 1) / 2 observations to be
  # non-singular in p coordinates; below 2p it returns NaN distances behind a
  # warning rather than failing, which would be a silent wrong answer.
  n_min <- 2L * n_coord
  if (n < n_min) {
    cli::cli_abort(
      c(
        "Too few observations for robust outlier detection.",
        "i" = "{n} row{?s} for {length(parts)} part{?s} ({n_coord} coordinate{?s}); {n_min} {?is/are} the minimum."
      )
    )
  }

  dots <- list(...)
  # `coda` is part of this function's contract, not a tunable: on raw parts
  # the distance is not a compositional statistic, and for a closed
  # composition the covariance is singular and outCoDa() fails outright.
  # Passing it is therefore an error, never silently ignored.
  if ("coda" %in% names(dots)) {
    cli::cli_abort(
      c(
        "{.arg coda} cannot be overridden.",
        "i" = "{.fn comp_outliers} measures distance in pivot ILR coordinates by definition."
      )
    )
  }
  args <- utils::modifyList(list(quantile = quantile, coda = TRUE), dots)
  res <- do.call(
    robCompositions::outCoDa,
    c(list(as.data.frame(data[parts])), args)
  )
  # outCoDa() computes both a robust and a classical estimate and returns both
  # (its own `method` argument is ignored upstream in v2.6.0: `covEst(x,
  # "robust")` is hardcoded and `method` only relabels a plot axis). The robust
  # fields are read here; the classical ones are not offered.
  if (anyNA(res$mahalDist)) {
    cli::cli_abort(
      c(
        "{.fn outCoDa} returned undefined distances.",
        "i" = "The covariance estimate is singular; check for collinear or duplicated parts."
      )
    )
  }

  tibble::tibble(
    .row = seq_len(n),
    distance = as.numeric(res$mahalDist),
    # A constant column, not an attribute: it keeps the cutoff readable from
    # the table alone, as robCompositions' own plot frame does.
    limit = res$limit,
    is_outlier = as.logical(res$outlierIndex)
  )
}
