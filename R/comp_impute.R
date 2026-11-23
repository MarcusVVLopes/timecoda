# Module 0: compositional zero imputation.
#
# comp_impute() wraps zCompositions' family of continuous compositional
# imputation methods behind one interface. It derives the per-part detection
# limit from the data (the smallest observed positive raw value) unless the
# caller overrides it via `dl`, and works on the raw (pre-closure) scale:
# declare_composition() calls it before closing the composition, so absolute
# row totals survive imputation.

#' Impute compositional zeroes
#'
#' Replaces zero values in a set of composition parts using one of the
#' continuous imputation methods of `zCompositions`. [declare_composition()]
#' calls this function on the raw (pre-closure) scale, before it closes the
#' composition.
#'
#' @param data A data frame containing the composition parts.
#' @param parts Character vector of column names making up the composition.
#' @param method Imputation method: `"lrEM"` (default, log-ratio EM,
#'   [zCompositions::lrEM()]), `"multRepl"` (simple multiplicative replacement,
#'   [zCompositions::multRepl()]), `"lrDA"` (log-ratio data augmentation,
#'   [zCompositions::lrDA()]), `"multLN"` (multiplicative lognormal,
#'   [zCompositions::multLN()]) or `"lrSVD"` (SVD-based,
#'   [zCompositions::lrSVD()]). `"multLN"` is called with `rob = TRUE` (the
#'   ROS-based alternative) unless `rob` is supplied through `...`; its
#'   classical-MLE default (`rob = FALSE`) fails inside `comp_impute()` with
#'   "could not find function \"survreg\"". `zCompositions::multKM()`
#'   (Kaplan-Meier based) fails in the same way and is not offered.
#' @param ... Additional arguments forwarded to the chosen method (most
#'   usefully `dl`, to override the derived per-part detection limit, or a
#'   method-specific control such as `lrEM`'s `rob`/`tolerance`/`max.iter`).
#'   The internal defaults (`label = 0`, the derived `dl`, and
#'   `suppress.print = TRUE` where the method accepts it) fill in only the
#'   arguments not supplied.
#'
#' @return A list with elements `parts_imputed` (a tibble: the zero-free parts
#'   on the original scale), `imput_table` (a tibble summarizing zero
#'   occurrence per part: `part`, `n`, `n_zero`, `pct_zero`, `min_nonzero`), and
#'   `impute_method` (the `method` used, or `NA_character_` when no
#'   part had a zero). When no part has a zero, `parts_imputed` equals the
#'   input `parts` columns unchanged and no imputation method is called.
#'
#' @references
#' Palarea-Albaladejo, J., & Martín-Fernández, J.A. (2015). zCompositions -- R
#' package for multivariate imputation of left-censored data under a
#' compositional approach. *Chemometrics and Intelligent Laboratory Systems*,
#' 143, 85-96. \doi{10.1016/j.chemolab.2015.02.019}
#'
#' @examples
#' cd_data <- timeuse
#' cd_data$MVPA[1:3] <- 0
#' imp <- comp_impute(cd_data, c("MVPA", "LPA", "SB", "SPT"))
#' imp$imput_table
#'
#' @seealso [declare_composition()], which calls this internally.
#' @family declare and describe composition
#' @export
comp_impute <- function(data, parts,
                        method = c("lrEM", "multRepl", "lrDA", "multLN",
                                   "lrSVD"),
                        ...) {
  method <- rlang::arg_match(method)
  check_composition(data, parts, allow_zero = TRUE)
  part_mat <- as.matrix(data[parts])
  n <- nrow(part_mat)

  n_zero <- colSums(part_mat == 0)
  min_nonzero <- purrr::map_dbl(parts, function(p) {
    col <- part_mat[, p]
    nz <- col[col > 0]
    if (length(nz) == 0) NA_real_ else min(nz)
  })
  imput_table <- tibble::tibble(
    part = parts,
    n = n,
    n_zero = as.integer(n_zero),
    pct_zero = 100 * n_zero / n,
    min_nonzero = min_nonzero
  )

  imputed_mat <- part_mat
  needs_imputation <- any(part_mat == 0)
  if (needs_imputation) {
    all_zero <- parts[is.na(min_nonzero)]
    if (length(all_zero) > 0) {
      cli::cli_abort(
        c(
          "Cannot impute: part{?s} {.val {all_zero}} {?is/are} entirely zero.",
          "i" = "Every part needs at least one positive value to set a detection limit."
        )
      )
    }
    # Per-part detection limit on the original scale: the smallest observed
    # positive value. Imputed zeroes fall below it. Users may override `dl`
    # (or any other method-specific control) through `...`.
    dl <- matrix(
      rep(min_nonzero, each = n), nrow = n,
      dimnames = list(NULL, parts)
    )
    # Dispatch to the chosen zCompositions function. multKM() and multLN()'s
    # classical rob = FALSE path both call NADA::cenmle(), and NADA declares
    # survival under Depends, not Imports, so cenmle() finds survreg() only
    # when survival is attached with library(), which a package must not do
    # on a caller's behalf. multKM() is therefore not offered, and multLN()
    # defaults to rob = TRUE, the ROS-based path that does not use survreg().
    defaults <- switch(method,
      lrEM = list(label = 0, dl = dl, suppress.print = TRUE),
      multLN = list(label = 0, dl = dl, rob = TRUE),
      list(label = 0, dl = dl)
    )
    fn <- switch(method,
      lrEM = zCompositions::lrEM,
      multRepl = zCompositions::multRepl,
      lrDA = zCompositions::lrDA,
      multLN = zCompositions::multLN,
      lrSVD = zCompositions::lrSVD
    )
    imputed_mat <- as.matrix(do.call(
      fn, c(list(part_mat), utils::modifyList(defaults, list(...)))
    ))
    # zCompositions deletes observations (and parts) that are almost entirely
    # zero, behind a warning (see its `z.delete` argument). Returning fewer
    # rows than it was given would leave `composition` and `parts_imputed` out
    # of step with `data_raw`, which every downstream consumer binds against
    # column-wise, so a changed shape is an error here.
    if (!identical(dim(imputed_mat), dim(part_mat))) {
      empty <- unname(which(rowSums(part_mat > 0) == 0))
      # Worded here rather than with cli's `{?}` markers: the quantity they
      # read carries across bullets, so the counts in the first bullet would
      # drive the plurals in the third.
      detail <- if (length(empty) > 0) {
        one <- length(empty) == 1
        c("i" = paste0(
          if (one) "Row " else "Rows ", paste(empty, collapse = ", "),
          if (one) " is" else " are", " entirely zero; drop ",
          if (one) "it" else "them", " first."
        ))
      }
      cli::cli_abort(c(
        "{.fn {method}} returned {nrow(imputed_mat)} of {n} rows and {ncol(imputed_mat)} of {length(parts)} parts.",
        "i" = "{.pkg zCompositions} drops observations that are almost entirely zero, which would leave the composition out of step with {.arg data}.",
        detail
      ))
    }
    dimnames(imputed_mat) <- list(NULL, parts)
    if (any(imputed_mat <= 0)) {
      cli::cli_abort(
        c(
          "{.fn {method}} returned non-positive values.",
          "i" = "Compositional parts must stay strictly positive after imputation."
        )
      )
    }
  }

  list(
    parts_imputed = tibble::as_tibble(as.data.frame(imputed_mat)),
    imput_table = imput_table,
    # NA, not `method`, when nothing needed imputing: the requested method
    # was never run, so echoing it back would be misleading.
    impute_method = if (needs_imputation) method else NA_character_
  )
}
