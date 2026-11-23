# Module 2: summarizing reallocations.
#
# `summary.realloc()` returns a `summary_realloc` object: the reallocations as a
# tidy tibble carrying the parts' display labels, plus the metadata the
# publication table needs. Printing it shows a text report; `pub_ready()`
# reshapes it into the wide table used in time-use papers: one row per
# reallocation pair, the effect (+ significance flag) and confidence interval
# grouped under a column spanner for each `delta`.

#' Summarize time-reallocation effects
#'
#' Tabulates the reallocations of a [predict_realloc()] result as a tidy
#' tibble at full numeric precision: one row per receiving part, donor part
#' and reallocation magnitude (`delta`), with the parts' display labels.
#' Printing the result shows a text report, as with [summary.lm()].
#' [pub_ready()] builds the wide publication table from it, with the effect
#' estimate and its confidence interval under a column spanner for each
#' `delta`.
#'
#' @param object A `realloc` object from [predict_realloc()]. An infeasible
#'   transfer (a row with `NA` estimates, the [predict_realloc()] default) is
#'   kept, and shows as `-` in the [pub_ready()] table.
#' @param labels Optional named character vector mapping raw part names to
#'   display labels (e.g., `c(MVPA = "Exercise", SB = "Sedentary")`). Parts
#'   not named keep their raw names, and the implicit "all remaining parts"
#'   side of a one-to-all reallocation is shown as "All other parts". Defaults
#'   to the labels of the [declare_composition()] object behind the fit (if
#'   any); pass a vector here to override them.
#' @param ... Unused; present for S3 consistency. Passing `display`, `digits`,
#'   `p_flag` or `conf_level` is an error: rounding and the significance flag
#'   are set in [pub_ready()], and the confidence level in [predict_realloc()].
#'
#' @return An object of class `summary_realloc`, printed as a text report by
#'   [print.summary_realloc()] and passed to [pub_ready()] for a `gt` table. A
#'   list with elements:
#'   \describe{
#'     \item{`table`}{A tibble with one row per reallocation: `to`, `from`, the
#'       `by` grouping column and `response_level` when present, the display
#'       labels `to_label` and `from_label`, then `delta`, `estimate`,
#'       `conf.low`, `conf.high` and, for relative output, `p.value`. The
#'       reference row of an absolute prediction is dropped; a zero reallocation
#'       (`delta = 0`) is kept.}
#'     \item{`effect_label`, `type`, `conf_level`, `parts`, `by`}{
#'       Kept for [print.summary_realloc()] and [pub_ready()].}
#'     \item{`units`}{The composition's unit label, shown beside each `delta`
#'       in the [pub_ready()] table; `NULL` when none was declared.}
#'   }
#'
#' @examples
#' cd <- declare_composition(timeuse, c("MVPA", "LPA", "SB", "SPT"),
#'                           closure = 1440)
#' fit <- comp_iv(cd, outcome = "obesity", model = "glm", family = "binomial")
#' r <- predict_realloc(fit, delta = c(10, 20, 30), on_invalid = "na")
#' s <- summary(r, labels = c(MVPA = "Exercise"))
#' s$table
#'
#' @seealso [predict_realloc()] to build the `realloc` object;
#'   [plot_realloc()] for the continuous sweep figure; [pub_ready()] for the
#'   publication `gt` table built from this summary.
#' @family reallocation and prediction
#' @export
summary.realloc <- function(object, labels = attr(object, "labels"), ...) {
  if (!inherits(object, "realloc")) {
    cli::cli_abort("{.arg object} must be a {.cls realloc} object from {.fn predict_realloc}.")
  }
  removed <- intersect(c("display", "digits", "p_flag", "conf_level"),
                       names(list(...)))
  if (length(removed) > 0) {
    cli::cli_abort(c(
      "{.fn summary} of a {.cls realloc} object does not take {.arg {removed}}.",
      "i" = if (!identical(removed, "conf_level")) "Build the formatted table with {.code pub_ready(summary(x), digits = 2)}.",
      "i" = if ("conf_level" %in% removed) "Set the confidence level in {.code predict_realloc(conf_level = )}."
    ))
  }
  parts <- attr(object, "parts")
  if (!is.null(labels)) {
    if (!is.character(labels) || is.null(names(labels))) {
      cli::cli_abort("{.arg labels} must be a named character vector (part = label).")
    }
    bad <- setdiff(names(labels), parts)
    if (length(bad) > 0) {
      cli::cli_abort(c(
        "{.arg labels} names {.val {bad}}, which {?is/are} not a part of the composition.",
        "i" = "Composition parts are {.val {parts}}."
      ))
    }
  }
  by <- attr(object, "by")
  group_keys <- c(as.character(by),
                  if ("response_level" %in% names(object)) "response_level")

  # Display labels. Plain part names map straight through the `labels` lookup;
  # the compact multi-part strings of a custom reallocation ("SPT +30",
  # "SB -15, LPA -15") have each embedded part name substituted in place.
  relabel <- function(x) {
    out <- ifelse(x == tc_rest_key(), tc_rest_label(), x)
    if (!is.null(labels)) {
      hit <- x %in% names(labels)
      out[hit] <- unname(labels[x[hit]])
      compound <- !hit & grepl("[+-]", x)
      out[compound] <- purrr::map_chr(x[compound], tc_relabel_compound,
                                      labels = labels)
    }
    out
  }
  # Receiver then donor (parts in composition order, the pooled rest last),
  # then group, then magnitude.
  lvl <- c(parts, tc_rest_key())
  table <- tibble::as_tibble(as.data.frame(object)) |>
    # The reference row that an absolute prediction prepends (from/to NA).
    dplyr::filter(!is.na(.data$to), !is.na(.data$from)) |>
    dplyr::mutate(to_label = relabel(.data$to),
                  from_label = relabel(.data$from)) |>
    dplyr::arrange(factor(.data$to, levels = lvl),
                   factor(.data$from, levels = lvl),
                   dplyr::pick(dplyr::all_of(group_keys)), .data$delta) |>
    dplyr::select(dplyr::any_of(c(
      "to", "from", group_keys, "to_label", "from_label", "delta",
      "estimate", "conf.low", "conf.high", "p.value"
    )))

  structure(
    list(
      table = table,
      effect_label = attr(object, "effect_label"),
      type = attr(object, "realloc_type"),
      conf_level = attr(object, "conf_level"),
      parts = parts,
      by = by,
      units = attr(object, "units")
    ),
    class = "summary_realloc"
  )
}

#' Print a time-reallocation summary
#'
#' Prints a text report of a [summary.realloc()] object: a short header on the
#' effect and its reference, then the summary table as plain tibble output. The
#' report prints automatically at the console, as base R's `print.summary.lm()`
#' does.
#'
#' @param x A `summary_realloc` object, as returned by [summary.realloc()].
#' @param ... Unused; present for S3 consistency.
#'
#' @return `x`, invisibly.
#'
#' @examples
#' cd <- declare_composition(timeuse, c("MVPA", "LPA", "SB", "SPT"),
#'                           closure = 1440)
#' fit <- comp_iv(cd, outcome = "obesity", model = "glm", family = "binomial")
#' r <- predict_realloc(fit, delta = c(10, 20), from = "SB", to = "MVPA")
#' print(summary(r))
#'
#' @seealso [summary.realloc()] to build this object; [pub_ready()] for the
#'   publication `gt` table.
#' @export
print.summary_realloc <- function(x, ...) {
  tbl <- x$table
  anchor <- if (identical(x$type, "absolute")) {
    "at each reallocated composition"
  } else {
    "versus the reference composition"
  }
  deltas <- glue::glue_collapse(format(sort(unique(tbl$delta)), trim = TRUE),
                                sep = ", ")
  noun <- if (nrow(tbl) == 1) "reallocation" else "reallocations"
  print(glue::glue("
Time-reallocation effects: {x$effect_label} {anchor}
  {nrow(tbl)} {noun} of {tc_with_units(deltas, x$units)}; {format(x$conf_level * 100)}% confidence intervals
"))
  if (!is.null(x$by)) {
    print(glue::glue("  Estimated per level of {x$by}"))
  }
  n_na <- sum(is.na(tbl$estimate))
  if (n_na > 0) {
    print(glue::glue("  No estimate (a part driven non-positive): {n_na}"))
  }
  cat("\n")
  print(tbl)
  invisible(x)
}

# Substitute display labels into one compact multi-part reallocation string
# ("SPT +30" / "SB -15, LPA -15"): split it into its "part delta" moves and look
# each part name up exactly, so no part name can match inside another or inside
# a label already substituted, and regex metacharacters in a name are inert.
tc_relabel_compound <- function(s, labels) {
  moves <- strsplit(s, ", ", fixed = TRUE)[[1]]
  part <- sub(" [+-][0-9.]+$", "", moves)
  hit <- part %in% names(labels)
  moves[hit] <- paste0(labels[part[hit]],
                       substring(moves[hit], nchar(part[hit]) + 1))
  paste(moves, collapse = ", ")
}

# Build summary_realloc's publication gt table: per `delta`, a column spanner
# over the effect formatted to `digits` (with `*` when p < `p_flag`) and a
# "low, high" CI string; a `by` group and a multinomial `response_level` are
# extra row keys. The heading uses absolute pixel sizes, because gt scales
# percentages unreliably in LaTeX, where its default title renders at ~20pt.
# The styling is written out inline deliberately, so that no builder depends
# on a shared helper for it. Do not "de-duplicate" it.
tc_pub_ready_realloc <- function(x, digits = 2, p_flag = 0.05) {
  df <- as.data.frame(x$table)
  est_fmt <- paste0("%.", digits, "f")
  ci_fmt <- paste0("%.", digits, "f, %.", digits, "f")
  flag <- if ("p.value" %in% names(df)) {
    ifelse(!is.na(df$p.value) & df$p.value < p_flag, "*", "")
  } else {
    rep("", nrow(df))
  }
  group_keys <- c(x$by, if ("response_level" %in% names(df)) "response_level")
  key <- c("to", "from", group_keys, "to_label", "from_label")
  long <- df[key]
  long$delta <- df$delta
  long$est <- ifelse(is.na(df$estimate), NA_character_,
                     paste0(sprintf(est_fmt, df$estimate), flag))
  long$ci <- ifelse(is.na(df$conf.low), NA_character_,
                    sprintf(ci_fmt, df$conf.low, df$conf.high))

  deltas <- sort(unique(long$delta))
  wide <- purrr::map(deltas, function(d) {
    sub <- long[long$delta == d, c(key, "est", "ci")]
    names(sub) <- c(key, paste0("est_", d), paste0("ci_", d))
    sub
  }) |>
    purrr::reduce(function(a, b) merge(a, b, by = key, all = TRUE))

  # Order rows by receiver then donor (parts first, the pooled rest last),
  # then group.
  lvl <- c(x$parts, tc_rest_key())
  ord_keys <- c(list(factor(wide$to, levels = lvl),
                     factor(wide$from, levels = lvl)),
                unname(as.list(wide[group_keys])))
  wide <- wide[do.call(order, ord_keys), , drop = FALSE]
  per_delta <- as.vector(rbind(paste0("est_", deltas), paste0("ci_", deltas)))
  col_order <- c(group_keys, "to_label", "from_label", per_delta)
  g <- gt::gt(tibble::as_tibble(wide[, col_order, drop = FALSE]))

  labels <- list(
    to_label = "Allocated to (+)",
    from_label = "Time from (-)"
  )
  for (k in group_keys) {
    labels[[k]] <- if (k == "response_level") "Outcome" else tools::toTitleCase(k)
  }
  ci_label <- paste0(format(x$conf_level * 100), "% CI")
  for (d in deltas) {
    labels[[paste0("est_", d)]] <- x$effect_label
    labels[[paste0("ci_", d)]] <- ci_label
    g <- gt::tab_spanner(
      g, label = tc_with_units(d, x$units),
      columns = c(paste0("est_", d), paste0("ci_", d))
    )
  }
  g <- do.call(gt::cols_label, c(list(g), labels))
  g <- gt::cols_align(g, align = "left",
                      columns = c(group_keys, "to_label", "from_label"))
  g <- gt::tab_header(
    g,
    title = "Time-reallocation effects",
    subtitle = paste(x$effect_label, if (identical(x$type, "absolute")) {
      "at each reallocated composition"
    } else {
      "versus the reference composition"
    })
  )
  g <- gt::sub_missing(g, missing_text = "-")
  if ("p.value" %in% names(x$table)) {
    g <- gt::tab_source_note(g, source_note = paste0("* p < ", p_flag))
  }
  list(
    effects = g |>
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
