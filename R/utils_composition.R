# Module 5: helpers (utils_*.R)
#
# ILR + composition helpers, generalized to a D-part composition.
# Single ILR backend: `robCompositions::pivotCoord()` (do not mix in
# compositions/epicoda transforms; `compositions` is for descriptive helpers
# only).

#' Validate a compositional dataset
#'
#' Checks that the composition parts are present in the data, numeric,
#' strictly positive, free of missing values, and (optionally) that every row
#' sums to the expected total.
#'
#' @param data A data frame containing the composition.
#' @param parts Character vector of column names making up the composition
#'   (at least two parts).
#' @param total Expected row total for the parts (e.g., `1440` for minutes in
#'   a 24-hour day). If `NULL` (default), row totals are not checked.
#' @param tol Numeric tolerance for the row-total check, as a proportion of
#'   `total` (default `0.01`, i.e., 1%). The margin allows for the imprecision
#'   of accelerometer- or diary-derived time-use minutes, which rarely sum to
#'   exactly 1440.
#' @param allow_zero Logical. If `FALSE` (default), parts must be strictly
#'   positive. If `TRUE`, zeroes are permitted and only negative values are
#'   rejected, for data whose zeroes are imputed before any log-ratio
#'   transformation, as in [declare_composition()].
#'
#' @return `data`, invisibly, if all checks pass; otherwise an informative
#'   error.
#'
#' @examples
#' check_composition(timeuse, c("MVPA", "LPA", "SB", "SPT"), total = 1440)
#'
#' @family composition helpers
#' @export
check_composition <- function(
  data,
  parts,
  total = NULL,
  tol = 0.01,
  allow_zero = FALSE
) {
  if (!is.data.frame(data)) {
    cli::cli_abort(
      "{.arg data} must be a data frame, not {.cls {class(data)[1]}}."
    )
  }
  if (!is.character(parts) || length(parts) < 2) {
    cli::cli_abort(
      "{.arg parts} must be a character vector of at least two column names."
    )
  }
  missing_parts <- setdiff(parts, names(data))
  if (length(missing_parts) > 0) {
    cli::cli_abort(
      "Composition part{?s} {.val {missing_parts}} not found in {.arg data}."
    )
  }
  pdata <- data[parts]
  non_numeric <- parts[!purrr::map_lgl(pdata, is.numeric)]
  if (length(non_numeric) > 0) {
    cli::cli_abort(
      "Composition part{?s} {.val {non_numeric}} must be numeric."
    )
  }
  has_na <- parts[purrr::map_lgl(pdata, anyNA)]
  if (length(has_na) > 0) {
    cli::cli_abort(
      "Composition part{?s} {.val {has_na}} contain{?s/} missing values; compositional data must be complete."
    )
  }
  cmp <- if (allow_zero) `<` else `<=`
  offending <- parts[purrr::map_lgl(pdata, ~ any(cmp(.x, 0)))]
  if (length(offending) > 0) {
    if (allow_zero) {
      cli::cli_abort(
        c(
          "Composition part{?s} {.val {offending}} contain{?s/} negative values.",
          "i" = "All parts of a composition must be non-negative."
        )
      )
    } else {
      cli::cli_abort(
        c(
          "Composition part{?s} {.val {offending}} contain{?s/} zero or negative values.",
          "i" = "All parts of a composition must be strictly positive."
        )
      )
    }
  }
  if (!is.null(total)) {
    tc_check_positive_scalar(total, "total")
    row_sums <- rowSums(pdata)
    off <- which(abs(row_sums - total) > tol * total)
    if (length(off) > 0) {
      cli::cli_abort(
        c(
          "{length(off)} row{?s} of {.arg data} do{?es/} not sum to {.val {total}} across {.val {parts}}.",
          "i" = "First offending row: {off[1]} (sum = {round(row_sums[off[1]], 6)})."
        )
      )
    }
  }
  invisible(data)
}

#' Isometric log-ratio (pivot) coordinates of a composition
#'
#' Transforms a `D`-part composition into `D - 1` pivot ILR coordinates via
#' [robCompositions::pivotCoord()], appended to `data` as plain numeric
#' columns `ilr1 ... ilr{D-1}`. The first coordinate expresses the pivot part
#' relative to the geometric mean of the remaining parts. Reallocation
#' predictions are invariant to the choice of pivot when the model is linear
#' in the coordinates or holds their full quadratic ([ilr()]).
#'
#' @inheritParams check_composition
#' @param pivot The part forming the numerator of the first coordinate:
#'   either a column name from `parts` or its position. Defaults to the first
#'   part.
#'
#' @return `data` as a tibble with `length(parts) - 1` additional numeric
#'   columns named `ilr1`, `ilr2`, ....
#'
#' @references
#' Filzmoser, P., Hron, K., & Templ, M. (2018). *Applied Compositional Data
#'   Analysis: With Worked Examples in R.* Springer.
#'   \doi{10.1007/978-3-319-96422-5}
#'
#' @examples
#' ilr_transform(timeuse, c("MVPA", "LPA", "SB", "SPT"), total = 1440)
#'
#' @family composition helpers
#' @export
ilr_transform <- function(data, parts, total = NULL, pivot = 1) {
  check_composition(data, parts, total = total)
  pivot <- tc_match_pivot(pivot, parts)
  existing <- intersect(tc_coord_names(length(parts)), names(data))
  if (length(existing) > 0) {
    cli::cli_abort(
      "{.arg data} already contains column{?s} {.val {existing}}; rename or drop {?it/them} before calling {.fn ilr_transform}."
    )
  }
  coords <- tc_pivot_coords(data[parts], pivot)
  out <- tibble::as_tibble(data)
  out[colnames(coords)] <- coords
  out
}

#' Back-transform ILR coordinates to a composition
#'
#' Inverts [ilr_transform()]: converts `D - 1` pivot ILR coordinates back to
#' the `D`-part composition, closed to `total`.
#'
#' @param coords A data frame or matrix of ILR coordinates with `D - 1`
#'   columns, in the column order produced by [ilr_transform()].
#' @param parts Character vector of the `D` composition part names, in the
#'   same order passed to [ilr_transform()].
#' @param total Row total the composition is closed to (e.g., `1440`).
#'   Defaults to `1`, returning proportions.
#' @param pivot The pivot used in [ilr_transform()] (name or position).
#'
#' @return A tibble with one column per part, rows summing to `total`.
#'
#' @references
#' Filzmoser, P., Hron, K., & Templ, M. (2018). *Applied Compositional Data
#'   Analysis: With Worked Examples in R.* Springer.
#'   \doi{10.1007/978-3-319-96422-5}
#'
#' @examples
#' ilr <- ilr_transform(timeuse, c("MVPA", "LPA", "SB", "SPT"), total = 1440)
#' comp <- ilr_backtransform(
#'   ilr[c("ilr1", "ilr2", "ilr3")],
#'   parts = c("MVPA", "LPA", "SB", "SPT"),
#'   total = 1440
#' )
#' all.equal(comp$MVPA, timeuse$MVPA)
#'
#' @family composition helpers
#' @export
ilr_backtransform <- function(coords, parts, total = 1, pivot = 1) {
  if (is.matrix(coords)) {
    coords <- as.data.frame(coords)
  }
  if (!is.data.frame(coords)) {
    cli::cli_abort(
      "{.arg coords} must be a data frame or matrix of ILR coordinates."
    )
  }
  if (ncol(coords) != length(parts) - 1) {
    cli::cli_abort(
      "{.arg coords} has {ncol(coords)} column{?s} but {length(parts)} parts imply {length(parts) - 1} ILR coordinates."
    )
  }
  tc_check_positive_scalar(total, "total")
  pivot <- tc_match_pivot(pivot, parts)
  # pivotCoordInv() returns parts reordered pivot-first, so the columns are
  # renamed to reflect that order, then reindexed back to `parts`.
  comp <- as.matrix(robCompositions::pivotCoordInv(as.matrix(coords)))
  colnames(comp) <- c(parts[pivot], parts[-pivot])
  comp <- comp[, parts, drop = FALSE]
  comp <- comp / rowSums(comp) * total
  tibble::as_tibble(comp)
}

#' Geometric mean composition
#'
#' Computes the closed geometric mean of a composition on the simplex. The
#' sample's closed geometric mean, stored in a [comp_iv()] fit, is the default
#' reference composition of [predict_realloc()].
#'
#' @inheritParams check_composition
#' @param total Row total to close the result to (e.g., `1440` for minutes per
#'   day). Defaults to `1`, returning proportions.
#' @param robust Passed to [compositions::mean.acomp()]. Defaults to `FALSE`
#'   (the classical geometric mean); pass `TRUE` for the MCD-based robust
#'   centre.
#' @param ... Passed on to [compositions::mean.acomp()].
#'
#' @return A one-row tibble with one column per part, summing to `total`.
#'
#' @references
#' Aitchison, J. (1986). *The Statistical Analysis of Compositional Data.*
#'   Chapman & Hall.
#'
#' @examples
#' geom_mean_comp(timeuse, c("MVPA", "LPA", "SB", "SPT"), total = 1440)
#'
#' @family composition helpers
#' @export
geom_mean_comp <- function(data, parts, total = 1, robust = FALSE, ...) {
  check_composition(data, parts)
  tc_check_positive_scalar(total, "total")
  # Centre via the compositions API: the closed geometric mean of the acomp
  # (mean() dispatches to compositions' mean.acomp), then re-closed to
  # `total` with clo().
  centre <- mean(compositions::acomp(data[parts]), robust = robust, ...)
  g <- compositions::clo(centre, total = total)
  g <- stats::setNames(as.numeric(g), parts)
  tibble::as_tibble(as.data.frame(t(g)))
}

#' Compositional variation matrix
#'
#' Computes the variation matrix `T`, with entries
#' `T[i, j] = var(log(x_i / x_j))`. Small values indicate parts that are
#' highly proportional.
#'
#' @inheritParams check_composition
#' @param robust Passed to [compositions::variation.acomp()]. Defaults to
#'   `FALSE` (the classical variation matrix); pass `TRUE` for the MCD-based
#'   robust variant.
#' @param ... Passed on to [compositions::variation.acomp()].
#'
#' @return A symmetric numeric matrix with `parts` as row and column names.
#'
#' @references
#' Aitchison, J. (1986). *The Statistical Analysis of Compositional Data.*
#'   Chapman & Hall.
#'
#' @examples
#' variation_matrix(timeuse, c("MVPA", "LPA", "SB", "SPT"))
#'
#' @family composition helpers
#' @export
variation_matrix <- function(data, parts, robust = FALSE, ...) {
  check_composition(data, parts)
  v <- compositions::variation(compositions::acomp(data[parts]), robust = robust, ...)
  v <- unclass(v)
  dimnames(v) <- list(parts, parts)
  v
}

# Internal helpers ------------------------------------------------------------

# A single positive scalar (`total`/`closure`-style arguments). `hint` adds an
# info bullet for arguments that also accept something other than a number.
tc_check_positive_scalar <- function(x, arg_name, hint = NULL,
                                     call = rlang::caller_env()) {
  if (!is.numeric(x) || length(x) != 1 || is.na(x) || x <= 0) {
    cli::cli_abort(
      c(
        "{.arg {arg_name}} must be a single positive number.",
        if (!is.null(hint)) c("i" = hint)
      ),
      call = call
    )
  }
}

# Names of the ILR coordinate columns for a D-part composition.
tc_coord_names <- function(n_parts) {
  paste0("ilr", seq_len(n_parts - 1))
}

# Resolve `pivot` (name or position) to a position within `parts`.
tc_match_pivot <- function(pivot, parts, call = rlang::caller_env()) {
  if (is.character(pivot)) {
    idx <- match(pivot, parts)
    if (length(pivot) != 1 || is.na(idx)) {
      cli::cli_abort(
        "{.arg pivot} must be one of {.val {parts}}, not {.val {pivot}}.",
        call = call
      )
    }
    return(idx)
  }
  if (
    !is.numeric(pivot) ||
      length(pivot) != 1 ||
      is.na(pivot) ||
      pivot != as.integer(pivot) ||
      pivot < 1 ||
      pivot > length(parts)
  ) {
    cli::cli_abort(
      "{.arg pivot} must be a part name or a single integer between 1 and {length(parts)}.",
      call = call
    )
  }
  as.integer(pivot)
}

# Pivot ILR coordinates of a positive composition matrix/data frame, as a
# plain numeric matrix with columns ilr1 ... ilr{D-1}. The as.data.frame()
# conversion is required, not defensive: robCompositions::pivotCoord()
# mishandles a tibble's `[, i]` subsetting (1-column tibble instead of a
# vector), and callers may pass a tibble-backed `data[parts]`.
tc_pivot_coords <- function(x, pivot = 1L) {
  coords <- as.matrix(robCompositions::pivotCoord(
    as.data.frame(x),
    pivotvar = pivot
  ))
  colnames(coords) <- tc_coord_names(ncol(x))
  coords
}
