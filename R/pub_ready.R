# Publication-ready gt tables from a timecoda summary.
#
# `pub_ready()` is the single entry point that turns a `timecoda` summary
# object into a publication-style `gt` table, with the numbers rearranged into
# a manuscript table's layout. `summary()` methods return tidy data and print
# plain text (like `summary.lm()`); building the table is a separate call.
# Dispatch is a `switch()` on the summary object's class, mirroring
# the package's `tc_model_meta()` registry pattern, one case per summary type,
# each built in its own type's file. `pub_ready()` itself is gated on
# `rlang::check_installed("gt")` and aborts if `gt` is missing, since calling
# it is an explicit request for a `gt` table.

#' Build publication-ready tables from a timecoda summary
#'
#' Turns the object returned by a `timecoda` `summary()` method into a named
#' list of [gt::gt()] tables, with the numbers arranged in the layout of a
#' manuscript or report table.
#'
#' @param x A summary object from a `timecoda` `summary()` method: a
#'   `summary_compD` (from [summary.compD()]), a `summary_comp_iv` (from
#'   [summary.comp_iv()]) or a `summary_realloc` (from [summary.realloc()]).
#' @param digits Number of decimal places in the rendered table: a single
#'   non-negative whole number. Default `2`.
#' @param p_flag Significance threshold below which the `summary_realloc`
#'   table flags an effect with `*` (relative output only). A single number
#'   from 0 to 1; default `0.05`. At `0` no effect is flagged. The other
#'   summary types ignore it.
#'
#' @return A named list of [gt::gt()] tables, one per table the summary carries,
#'   with no entry for a table it lacks.
#'   \describe{
#'     \item{`summary_compD`}{`desc` (each part's arithmetic Mean \eqn{\pm} SD
#'       of the original data, before zero imputation and closure, alongside
#'       its geometric mean on the closed composition and, when the composition
#'       is not already closed to proportions or percentages, its share of the
#'       total; each column header names its scale) and `var_matrix`.}
#'     \item{`summary_comp_iv`}{`effects` (the omnibus composition row followed
#'       by one "part vs remaining" row per part, with the effect, its
#'       confidence interval and p-value; a per-category model adds an Outcome
#'       column and a moderated fit a per-part moderation p-value).}
#'     \item{`summary_realloc`}{`effects` (one row per reallocation pair, with
#'       the effect, flagged `*` below `p_flag`, and its confidence interval
#'       under a column spanner for each `delta`; a `by` grouping or a
#'       multinomial outcome adds a column).}
#'   }
#'
#' @examplesIf rlang::is_installed("gt")
#' cd <- declare_composition(timeuse, c("MVPA", "LPA", "SB", "SPT"),
#'                           closure = 1440, units = "min/day")
#' s <- summary(cd)
#' pub_ready(s)$desc
#' pub_ready(s)$var_matrix
#'
#' m <- comp_iv(cd, outcome = "obesity", model = "glm", family = "binomial")
#' pub_ready(summary(m))$effects
#'
#' r <- predict_realloc(m, delta = c(10, 20, 30), on_invalid = "na")
#' pub_ready(summary(r))$effects
#'
#' @seealso [summary.compD()], [summary.comp_iv()] and [summary.realloc()] for
#'   the underlying summaries.
#' @export
pub_ready <- function(x, digits = 2, p_flag = 0.05) {
  if (!(rlang::is_integerish(digits, n = 1, finite = TRUE) && digits >= 0)) {
    cli::cli_abort("{.arg digits} must be a single non-negative whole number.")
  }
  if (!(is.numeric(p_flag) && length(p_flag) == 1 && !is.na(p_flag) &&
        p_flag >= 0 && p_flag <= 1)) {
    cli::cli_abort("{.arg p_flag} must be a single number from 0 to 1.")
  }
  rlang::check_installed("gt", reason = "to build a publication-ready table.")
  switch(class(x)[1],
    summary_compD = tc_pub_ready_compD(x, digits),
    summary_comp_iv = tc_pub_ready_comp_iv(x, digits),
    summary_realloc = tc_pub_ready_realloc(x, digits, p_flag),
    cli::cli_abort("No {.fn pub_ready} method yet for {.cls {class(x)[1]}}.")
  )
}
