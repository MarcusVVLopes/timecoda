# Module 0: describing a compD object.
#
# `summary.compD()` packages the descriptive steps of the
# "composition-variable" vignette (centre, per-part spread, variation matrix)
# into one dispatched method, returning a `summary_compD` object: tidy tibbles
# plus the metadata needed to format them. Printing it (automatic at the
# console, like `summary.lm()`) shows a text report; `pub_ready()` builds the
# formatted `gt` table on request. The geometric mean and variation matrix
# derive from the imputed, closed composition; the per-part spread statistics
# and row totals describe `data_raw`, the data as recorded.

#' Describe a compositional variable
#'
#' Describes a [declare_composition()] object by a per-part descriptive table
#' and the compositional [variation_matrix()]. The table leads with the closed
#' geometric mean of the imputed composition (the compositional centre, from
#' [geom_mean_comp()]); its mean, SD, quartiles and range describe the original
#' data, before zero imputation and closure. Printing the result shows a text
#' report, as with [summary.lm()], and [pub_ready()] builds a publication `gt`
#' table from it.
#'
#' @param object A `compD` object from [declare_composition()].
#' @param desc Logical; include the per-part descriptive table. Default `TRUE`.
#' @param var_matrix Logical; include the variation matrix. Default `TRUE`.
#' @param ... Unused; present for S3 consistency.
#'
#' @return An object of class `summary_compD`, printed as a text report by
#'   [print.summary_compD()] and passed to [pub_ready()] for a `gt` table. A
#'   list with elements:
#'   \describe{
#'     \item{`desc`}{One row per part: `part`, `label` (display label), then
#'       `geom_mean` (the closed geometric mean of the imputed composition), and
#'       `mean`, `sd`, `median`, `p25`, `p75`, `min`, `max` of the original
#'       data, before zero imputation and closure. `NULL` when `desc = FALSE`.}
#'     \item{`var_matrix`}{The variation matrix of the imputed composition,
#'       with leading `part` and `label` columns. `NULL` when
#'       `var_matrix = FALSE`.}
#'     \item{`parts`, `labels`, `closure`, `units`}{The part names, resolved
#'       display labels and closure of the `compD`, carried through for
#'       [pub_ready()], and its unit label (`NULL` when none was declared),
#'       shown by [print.summary_compD()].}
#'     \item{`row_totals`}{A one-row tibble (`mean`, `sd`, `min`, `max`) of the
#'       raw part sums per observation, before imputation or closure.
#'       [print.summary_compD()] reports it to show how much row totals
#'       varied in the original data.}
#'     \item{`n_obs`}{The number of observations in the composition.}
#'     \item{`imput_table`, `impute_method`}{Carried through from the `compD`
#'       object (see [declare_composition()]), for [print.summary_compD()]'s
#'       zero-detection report.}
#'     \item{`outliers`}{A one-row tibble (`n_outlier`, `pct`, `limit`)
#'       summarizing the `.outlier` columns [declare_composition()] appended
#'       to `data_raw`, for [print.summary_compD()]'s outlier line. `NULL`
#'       when outliers were not flagged.}
#'   }
#'
#' @examples
#' cd <- declare_composition(timeuse, c("MVPA", "LPA", "SB", "SPT"),
#'                           closure = 1440)
#' s <- summary(cd)
#' s$desc
#'
#' @seealso [declare_composition()] to build the `compD` object;
#'   [plot.compD()] for the figures; [pub_ready()] for a
#'   formatted `gt` table built from this summary.
#' @family declare and describe composition
#' @export
summary.compD <- function(object, desc = TRUE, var_matrix = TRUE, ...) {
  if (!inherits(object, "compD")) {
    cli::cli_abort("{.arg object} must be a {.cls compD} object.")
  }
  comp_df <- tc_comp_df(object)
  parts <- object$parts
  # Both tables are built in `parts` order (variation_matrix() names its
  # dimensions with `parts`), so one label lookup serves every use below.
  labs <- tc_part_labels(parts, object$labels)

  # The parts as recorded (zeroes included, before closure): the spread
  # statistics and the row totals describe these; only the geometric mean and
  # the variation matrix are on the imputed, closed composition.
  raw <- as.matrix(object$data_raw[parts])

  desc_tbl <- NULL
  if (isTRUE(desc)) {
    gm <- geom_mean_comp(comp_df, parts, total = object$closure)
    desc_tbl <- purrr::map2(parts, labs, function(p, lab) {
      x <- raw[, p]
      q <- stats::quantile(x, c(0.25, 0.75), names = FALSE)
      tibble::tibble(
        part = p,
        label = lab,
        # The compositional centre leads the per-part spread statistics.
        geom_mean = as.numeric(gm[[p]]),
        mean = mean(x),
        sd = stats::sd(x),
        median = stats::median(x),
        p25 = q[1],
        p75 = q[2],
        min = min(x),
        max = max(x)
      )
    }) |>
      purrr::list_rbind()
  }

  var_tbl <- NULL
  if (isTRUE(var_matrix)) {
    vm <- variation_matrix(comp_df, parts)
    var_tbl <- tibble::as_tibble(as.data.frame(vm), rownames = "part")
    var_tbl <- tibble::add_column(var_tbl, label = labs, .after = "part")
  }

  # How much the original observations' totals varied, independent of zero
  # handling.
  raw_totals <- rowSums(raw)

  # Condense declare_composition()'s `.outlier` columns for the text report;
  # NULL, and so silent, when outliers were not flagged.
  outliers <- NULL
  if (tc_has_outliers(object)) {
    flag <- object$data_raw$.outlier
    outliers <- tibble::tibble(
      n_outlier = sum(flag),
      pct = 100 * mean(flag),
      limit = object$data_raw$.outlier_limit[1]
    )
  }

  structure(
    list(
      desc = desc_tbl,
      var_matrix = var_tbl,
      parts = parts,
      labels = labs,
      closure = object$closure,
      units = object$units,
      row_totals = tibble::tibble(
        mean = mean(raw_totals), sd = stats::sd(raw_totals),
        min = min(raw_totals), max = max(raw_totals)
      ),
      n_obs = nrow(object$data_raw),
      imput_table = object$imput_table,
      impute_method = object$impute_method,
      outliers = outliers
    ),
    class = "summary_compD"
  )
}

#' Print a compositional variable summary
#'
#' Prints a [summary.compD()] object as a text report, automatically at the
#' console, as base R's `print.summary.lm()` does. The report describes the
#' composition by its parts, closure (with its unit label, when one was
#' declared), number of observations and how much the raw row totals varied
#' before imputation or closure. It states which parts had zeroes and, if any
#' were imputed, the method used, and how many observations were flagged as
#' compositional outliers when flagging was run. The requested tables follow as
#' plain tibble output.
#'
#' @param x A `summary_compD` object, as returned by [summary.compD()].
#' @param ... Unused; present for S3 consistency.
#'
#' @return `x`, invisibly.
#'
#' @examples
#' cd <- declare_composition(timeuse, c("MVPA", "LPA", "SB", "SPT"),
#'                           closure = 1440)
#' print(summary(cd))
#'
#' @seealso [summary.compD()] to build this object; [pub_ready()] for a
#'   formatted `gt` table.
#' @export
print.summary_compD <- function(x, ...) {
  rt <- x$row_totals
  print(glue::glue("
<compD> compositional variable summary
  Parts: {paste(x$labels, collapse = ', ')}
  Closure: {tc_with_units(x$closure, x$units)}
  Observations: {x$n_obs}

  Row totals (before closure): mean {tc_fmt(rt$mean)} (SD {tc_fmt(rt$sd)}), range {tc_fmt(rt$min)}-{tc_fmt(rt$max)}
"))

  zeroed <- tc_zeroed_parts(x$imput_table)
  if (nrow(zeroed) == 0) {
    print(glue::glue("  No zeroes detected."))
  } else {
    per_part <- glue::glue_collapse(
      glue::glue("{zeroed$part} ({zeroed$n_zero} of {zeroed$n}, {tc_fmt(zeroed$pct_zero)}%)"),
      sep = ", "
    )
    print(glue::glue("  Zeroes detected: {per_part}"))
    print(glue::glue("  Zeroes imputed via {x$impute_method}"))
  }

  # Silent when flagging never ran, as the zero block is when there are no
  # zeroes: the report states only what was found.
  if (!is.null(x$outliers)) {
    o <- x$outliers
    # The cutoff alone is uninterpretable, and `quantile` is not stored; it is
    # recovered exactly, since limit = sqrt(qchisq(quantile, df = D - 1)).
    quantile <- stats::pchisq(o$limit^2, df = length(x$parts) - 1)
    print(glue::glue(
      "  Outliers flagged: {o$n_outlier} of {x$n_obs} ({tc_fmt(o$pct)}%) beyond the robust Mahalanobis cutoff {tc_fmt(o$limit, 2)} (chi-squared {tc_fmt(quantile, 3)})"
    ))
  }

  if (!is.null(x$desc)) {
    cat("\n")
    print(glue::glue("
Descriptive statistics
  geom_mean: imputed, closed composition; other columns: original data
"))
    print(x$desc)
  }
  if (!is.null(x$var_matrix)) {
    cat("\n")
    print(glue::glue("Variation matrix"))
    print(x$var_matrix)
  }
  invisible(x)
}

# Build summary_compD's gt tables (per-part descriptives, variation matrix) and
# give both the compact heading style: absolute pixel sizes, because gt scales
# percentages unreliably in LaTeX, where its default title renders at ~20pt.
# The styling is written out here deliberately, not in a shared helper.
tc_pub_ready_compD <- function(x, digits = 2) {
  # Per-part descriptives: arithmetic "Mean +/- SD" of the original data beside
  # the geometric (compositional) mean, plus its share of the total.
  desc <- if (!is.null(x$desc)) {
    desc_disp <- tibble::tibble(
      part = x$desc$label,
      mean_sd = paste0(
        tc_fmt(x$desc$mean, digits), " \u00b1 ", tc_fmt(x$desc$sd, digits)
      ),
      geom_mean = x$desc$geom_mean
    )
    # Closed to proportions or percentages, the share column would repeat the
    # geometric mean; object$closure is the single source of truth for scale.
    show_pct <- !(isTRUE(all.equal(x$closure, 1)) ||
                    isTRUE(all.equal(x$closure, 100)))
    if (show_pct) {
      desc_disp$pct_total <- x$desc$geom_mean / x$closure * 100
    }
    # Each column header names its scale: the spread statistics describe the
    # original data, the geometric mean the composition closed to `closure`
    # (the declared unit label belongs to the closure scale).
    gt_tbl <- gt::gt(desc_disp) |>
      gt::tab_header(
        title = "Descriptive summary",
        subtitle = "Spread of the original data and the compositional centre per part"
      ) |>
      gt::fmt_number(columns = "geom_mean", decimals = digits) |>
      gt::cols_label(
        part = "Part",
        mean_sd = gt::md("Mean \u00b1 SD (original data)"),
        geom_mean = paste0(
          "Geometric mean (closed to ", tc_with_units(x$closure, x$units), ")"
        )
      )
    if (show_pct) {
      gt_tbl <- gt_tbl |>
        gt::fmt_number(columns = "pct_total", decimals = digits) |>
        gt::cols_label(pct_total = "% of total")
    }
    gt_tbl
  }

  # Variation matrix, display labels replacing the raw part names in the row
  # column and the per-part value columns.
  var_matrix <- if (!is.null(x$var_matrix)) {
    var_disp <- x$var_matrix
    var_disp$part <- var_disp$label
    var_disp$label <- NULL
    col_labels <- stats::setNames(
      as.list(c("Part", x$labels)), c("part", x$parts)
    )
    gt::gt(var_disp) |>
      gt::tab_header(
        title = "Variation matrix",
        subtitle = "Pairwise log-ratio variance, var(log(x_i / x_j))"
      ) |>
      gt::fmt_number(columns = dplyr::all_of(x$parts), decimals = digits) |>
      gt::cols_label(!!!col_labels)
  }

  purrr::compact(list(desc = desc, var_matrix = var_matrix)) |>
    purrr::map(\(g) {
      g |>
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
    })
}
