# Module 4: the reallocation transfer plot.
#
# `plot_realloc()` draws a `predict_realloc()` sweep over a signed-delta
# sequence, assembling one panel per part pair (or per part, for
# "one-to-all") with `patchwork`. The layout is derived from the number of
# parts, never from a hardcoded layout_matrix or slice(). Which magnitudes each panel
# can sweep comes from the same scheme grid predict_realloc() builds.
# A `by` grouping or a multinomial outcome's categories get one such grid
# each, stacked as rows.

#' Plot a time-reallocation sweep
#'
#' Sweeps [predict_realloc()] over a sequence of reallocation magnitudes in
#' both directions and plots each estimate with its confidence interval as a
#' ribbon, one panel per part pair or per part. The figure shows exactly the
#' estimates [predict_realloc()] returns for those magnitudes.
#'
#' @details
#' The panel layout follows `scheme`:
#'
#' * `"one-to-one"` (default): a lower-triangular matrix, one panel per
#'   *unordered* part pair (no mirrored duplicate). The panels are assembled
#'   with [patchwork::wrap_plots()], and the unused upper-triangle cells are
#'   left out of the figure.
#' * `"one-to-all"`: a compact wrapped grid, one panel per part, showing the
#'   proportional reallocation between that focal part and all the others
#'   (the implicit side that [predict_realloc()] labels `"(rest)"`).
#'
#' In both schemes, each panel sweeps time in either direction, up to
#' `delta_max` units of the composition. The sweep stops where a reallocation
#' would drive a part to zero or below at the reference composition, so
#' `on_invalid` has no visible effect. The x-axis shows the *magnitude* of the
#' reallocation, labelled as a positive amount on both sides of zero. Colour
#' gives the direction, since each half of the sweep is coloured by the part
#' gaining time, from a shared `palette` (viridis by default). The pooled
#' `"(rest)"` side of a one-to-all panel is drawn in neutral grey. For relative
#' output, the y-axis shows odds, hazard or rate ratios (log scale, reference
#' line at 1) or a mean difference (reference line at 0), according to the
#' model; for absolute output, it shows the predicted response. All panels share
#' one y-axis range and one legend, and their titles use the part labels of the
#' [declare_composition()] object. The figure carries no caption.
#'
#' Both halves of a panel meet at the reference value (an effect of exactly 1,
#' or 0 on the identity scale, at zero reallocated). An effect is significant
#' at the level `1 - conf_level` exactly where its confidence ribbon excludes
#' the dashed reference line.
#'
#' The reallocations can be grouped by a `by` covariate, by the outcome
#' categories of a multinomial fit, or by both. The `by` covariate is passed
#' through `...`, or is the factor `moderator` of the fit, which
#' [predict_realloc()] uses by default. The figure then stacks one full panel
#' grid per group as rows, each titled with its group (e.g., `sex: female`,
#' `outcome: B`, or both joined by a comma). Every grid shares one y-axis
#' range, and the figure has one legend.
#'
#' @param object A `comp_iv` / `timecoda_fit` from [comp_iv()].
#' @param delta_max Maximum amount reallocated in each direction, in the
#'   composition's units. Each panel sweeps from `-delta_max` to `delta_max`,
#'   capped per direction where a part would be driven non-positive at the
#'   reference composition.
#' @param n_points Number of points in the sweep across both directions, a
#'   single positive whole number; each direction gets `ceiling(n_points / 2)`
#'   of them.
#' @param type `"relative"` (default) or `"absolute"`; passed to
#'   [predict_realloc()].
#' @param scheme How time is reallocated, passed to [predict_realloc()].
#'   `"one-to-one"` (default) draws the pairwise transfer triangle;
#'   `"one-to-all"` draws one panel per part, reallocating proportionally
#'   between that part and all the others.
#' @param exclude Optional character vector of part names (not display
#'   labels) whose panels are dropped from the figure. Use it, for example,
#'   for a part whose balance is not associated with the outcome, so its
#'   panels show only noise. Excluded parts remain in the composition and
#'   still absorb their proportional share of a `"one-to-all"` reallocation,
#'   so the filter changes the display only.
#' @param palette Optional named character vector of colours keyed by part
#'   name (one colour per part). Defaults to the viridis palette
#'   ([scales::viridis_pal()]) over the fitted parts. Each half of a panel is
#'   coloured by its receiving part.
#' @param legend Draw the shared legend(s). `FALSE` hides all guides.
#' @param ... Further arguments passed to [predict_realloc()]: `reference`
#'   (which also sets the composition at which the feasible sweep is capped),
#'   `at`, `by` (one grid per level, see Details), and `conf_level`. Unused
#'   by the print method.
#'
#' @return Invisibly, a `realloc_plots` object: a list with `panels`, the
#'   individual `ggplot` panels, and `figure`, the assembled `patchwork`
#'   (for grouped results, one titled row per group). `panels` is named by the
#'   part names of each panel, in their declared order, joined by a dot
#'   (`"MVPA.SB"`; `"SB"` under `scheme = "one-to-all"`). A grouped result
#'   prefixes the group level (`"female.MVPA.SB"`). Part names that contain a
#'   dot can give two panels the same name (e.g., the pairs `a` with `b.c` and
#'   `a.b` with `c`); such panels are picked by position. The figure is drawn
#'   as a side effect when the function is called, and printing the object
#'   draws it again. A panel can be reused on its own or re-assembled with
#'   [patchwork::wrap_plots()].
#'
#' @examples
#' cd <- declare_composition(timeuse, c("MVPA", "LPA", "SB", "SPT"),
#'                           closure = 1440)
#' fit <- comp_iv(cd, outcome = "obesity", model = "glm", family = "binomial")
#' plot_realloc(fit, delta_max = 60, n_points = 9)
#'
#' # Proportional (one-to-all) reallocation, one panel per part
#' plot_realloc(fit, delta_max = 60, n_points = 9, scheme = "one-to-all")
#'
#' # The pairwise triangle without the three MVPA panels
#' plot_realloc(fit, delta_max = 60, n_points = 9, exclude = "MVPA")
#'
#' # One panel on its own, and two panels assembled side by side
#' p <- plot_realloc(fit, delta_max = 60, n_points = 9)
#' names(p$panels)
#' p$panels[["MVPA.SB"]]
#' patchwork::wrap_plots(p$panels[c("MVPA.SB", "MVPA.LPA")], ncol = 2)
#'
#' @seealso [comp_iv()] to fit the model; [predict_realloc()], which this
#'   function sweeps over; [summary.realloc()] for the tabular counterpart.
#' @family reallocation and prediction
#' @export
plot_realloc <- function(object, delta_max = 60, n_points = 41,
                         type = c("relative", "absolute"),
                         scheme = c("one-to-one", "one-to-all"),
                         exclude = NULL,
                         palette = NULL, legend = TRUE,
                         ...) {
  tc_check_fit(object)
  type <- match.arg(type)
  scheme <- match.arg(scheme)
  if (!is.numeric(delta_max) || length(delta_max) != 1 || delta_max <= 0) {
    cli::cli_abort("{.arg delta_max} must be a single positive number.")
  }
  if (!rlang::is_integerish(n_points, n = 1, finite = TRUE) || n_points < 1) {
    cli::cli_abort("{.arg n_points} must be a single positive whole number.")
  }
  # `exclude` prunes panels only: the fit, the reference composition, and every
  # estimate stay whole, so a filtered figure still matches predict_realloc().
  parts <- object$parts
  if (!is.null(exclude)) {
    if (!is.character(exclude)) {
      cli::cli_abort("{.arg exclude} must be a character vector of part names.")
    }
    bad <- setdiff(exclude, parts)
    if (length(bad) > 0) {
      cli::cli_abort(c(
        "{.arg exclude} names {.val {bad}}, which {?is/are} not in the fitted composition.",
        "i" = "Composition parts are {.val {parts}}.",
        "i" = "{.arg exclude} takes part names, not the {.fn declare_composition} display labels."
      ))
    }
    parts <- setdiff(parts, exclude)
    min_parts <- if (scheme == "one-to-all") 1L else 2L
    if (length(parts) < min_parts) {
      cli::cli_abort(c(
        "{.arg exclude} leaves {length(parts)} part{?s}, too few to plot a {.val {scheme}} reallocation.",
        "i" = "{.val {scheme}} needs at least {min_parts} part{?s} to remain."
      ))
    }
  }
  labs <- tc_part_labels(parts, object$comp$labels)
  # The receiving-part colours, keyed by part. The default viridis palette
  # spans every fitted part, so `exclude` never recolours a remaining one. A
  # one-to-all figure's pooled "(rest)" side gets a neutral grey (so it reads
  # as "everything else"), unless a user palette names it.
  if (is.null(palette)) {
    palette <- stats::setNames(scales::viridis_pal()(length(object$parts)),
                               object$parts)[parts]
  } else {
    if (!is.character(palette) || is.null(names(palette))) {
      cli::cli_abort(
        "{.arg palette} must be a named character vector of colours keyed by part."
      )
    }
    unmatched <- setdiff(parts, names(palette))
    if (length(unmatched) > 0) {
      cli::cli_abort(
        "{.arg palette} is missing a colour for part{?s} {.val {unmatched}}."
      )
    }
  }
  if (scheme == "one-to-all" && !tc_rest_key() %in% names(palette)) {
    palette[[tc_rest_key()]] <- "grey55"
  }
  dots <- rlang::list2(...)
  # The anchor predict_realloc() will use, so the sweep is capped where its
  # estimates stop being feasible.
  reference <- tc_check_reference(dots[["reference"]], object)

  cells <- tc_realloc_cells(object, parts, labs, delta_max, n_points, type,
                            scheme = scheme, ..., .reference = reference)
  if (length(cells) == 0) {
    cli::cli_abort(
      "No reallocation is feasible for any part at this {.arg delta_max}."
    )
  }

  all_df <- purrr::map(cells, "data") |> dplyr::bind_rows()
  y_lab <- tc_effect_label(object, type)
  # A ratio effect (odds, hazard, rate) is drawn on a log axis around 1; a
  # difference around 0. The registry, not the model name, says which.
  use_log <- type == "relative" &&
    tc_report_meta(object$model_type, object$family)$exponentiate
  ref_line <- if (type == "relative") (if (use_log) 1 else 0) else NULL
  # A shared y-range (computed once, across every pair and every group) keeps
  # panels on a comparable scale, the way facet_grid's default fixed scales did.
  ylim <- range(c(all_df$conf.low, all_df$conf.high, ref_line), na.rm = TRUE)

  # Legend keys are identical in every panel (so patchwork collects one guide).
  # One-to-all panels add the pooled "(rest)" side as a further key.
  legend_limits <- if (scheme == "one-to-all") c(parts, tc_rest_key()) else parts
  legend_labels <- if (scheme == "one-to-all") c(labs, tc_rest_label()) else labs

  # One symmetric break vector, shared by every panel, so all panels carry the
  # same x-axis ticks over the full [-delta_max, delta_max] span (the sweeps
  # themselves are capped per cell, but the axis is not).
  x_breaks <- scales::extended_breaks()(c(-delta_max, delta_max))
  # The x-axis title, with the composition's unit label when one was declared.
  # Bracketed, so not `tc_with_units()`.
  x_title <- if (is.null(object$comp$units)) {
    "Time reallocated"
  } else {
    as.character(glue::glue("Time reallocated ({object$comp$units})"))
  }

  # One grid of panels for a set of cells, its rows moved down by `offset`.
  grid_parts <- function(cells, offset = 0) {
    # In relative mode a zero reallocation is the reference against
    # itself: exactly `ref_line` with a zero-width interval. The sweep itself
    # holds only non-zero deltas, so anchor each receiving side's line at
    # (0, ref_line) so both halves visibly meet the reference; the anchor rows'
    # other columns are NA. (Absolute mode already carries the reference row
    # that predict_realloc() returns.)
    if (type == "relative") {
      cells <- purrr::map(cells, function(cell) {
        anchor <- tibble::tibble(
          signed_delta = 0, delta = 0, receiving = unique(cell$data$receiving),
          estimate = ref_line, conf.low = ref_line, conf.high = ref_line
        )
        df <- dplyr::bind_rows(cell$data, anchor)
        cell$data <- df[order(df$signed_delta), , drop = FALSE]
        cell
      })
    }
    rows <- purrr::map_dbl(cells, "row")
    cols <- purrr::map_dbl(cells, "col")
    panels <- purrr::map(cells, function(cell) {
      # Axis titles go on the outer edge of whatever shape the cells form: the
      # leftmost occupied column for y, and the bottom-most cell of each column
      # for x. Derived from the cells themselves, so it holds for the triangle
      # and the wrapped one-to-all grid alike.
      tc_realloc_panel(
        cell$data,
        y_lab, ylim, use_log, ref_line, palette, legend_limits, legend_labels,
        delta_max = delta_max, x_breaks = x_breaks, x_title = x_title,
        show_y_title = cell$col == min(cols),
        show_x_title = cell$row == max(rows[cols == cell$col]),
        row_label = cell$row_label, col_label = cell$col_label
      )
    })
    design <- purrr::map(cells, function(cell) {
      patchwork::area(t = cell$row + offset, l = cell$col,
                      b = cell$row + offset, r = cell$col)
    }) |>
      purrr::reduce(c)
    list(panels = panels, design = design, n_row = max(rows), n_col = max(cols))
  }

  # A `by` grouping (explicit, or a moderated fit's default, resolved as
  # predict_realloc() resolves it) and the categories of a multinomial outcome
  # each get their own full grid, stacked as rows in the order
  # predict_realloc() returns them.
  group_keys <- intersect(
    c(tc_resolve_by(object, dots[["by"]], dots[["at"]]), "response_level"),
    names(all_df)
  )
  if (length(group_keys) == 0) {
    g <- grid_parts(cells)
    panels <- g$panels
    figure <- patchwork::wrap_plots(g$panels, design = g$design)
  } else {
    key_of <- function(df) {
      do.call(paste, c(unname(as.list(df[group_keys])), sep = "\r"))
    }
    groups <- dplyr::distinct(all_df[group_keys])
    key_names <- ifelse(group_keys == "response_level", "outcome", group_keys)
    grids <- purrr::map(seq_len(nrow(groups)), function(i) {
      key <- key_of(groups[i, ])
      g <- grid_parts(purrr::map(cells, function(cell) {
        cell$data <- cell$data[key_of(cell$data) == key, , drop = FALSE]
        cell
      }), offset = 1)
      values <- unlist(groups[i, ])
      title <- paste(paste0(key_names, ": ", values), collapse = ", ")
      # The grid's title is a panel-less plot spanning its top row; it is not
      # a panel, so it is not in the returned list.
      header <- ggplot2::ggplot() +
        ggplot2::labs(title = title) +
        ggplot2::theme_void(base_size = 9) +
        ggplot2::theme(plot.title = ggplot2::element_text(face = "bold"))
      list(
        panels = rlang::set_names(
          g$panels, paste(paste(values, collapse = "."), names(g$panels), sep = ".")
        ),
        grid = patchwork::wrap_plots(
          c(list(header), g$panels),
          design = c(patchwork::area(t = 1, l = 1, b = 1, r = g$n_col), g$design),
          heights = c(0.001, rep(1, g$n_row))
        )
      )
    })
    panels <- purrr::list_flatten(purrr::map(grids, "panels"))
    figure <- patchwork::wrap_plots(purrr::map(grids, "grid"), ncol = 1)
  }
  figure <- figure +
    patchwork::plot_layout(guides = if (legend) "collect" else "keep")
  if (!legend) {
    figure <- figure & ggplot2::theme(legend.position = "none")
    panels <- purrr::map(panels, function(p) {
      p + ggplot2::theme(legend.position = "none")
    })
  }
  out <- structure(list(panels = panels, figure = figure),
                   class = "realloc_plots")
  # Draw as a side effect (through the shared, ggtern-pollution-resilient
  # renderer) so the figure displays at the console, and return the object
  # invisibly so it stays editable without auto-printing a second time.
  tc_draw_plots(figure)
  invisible(out)
}

#' @rdname plot_realloc
#' @param x A `realloc_plots` object returned by `plot_realloc()`.
#' @export
print.realloc_plots <- function(x, ...) {
  tc_draw_plots(x$figure)
  invisible(x)
}

# The panels of each directed `scheme`, keyed as `tc_realloc_schemes`: a
# function of the parts and their labels that returns one entry per panel.
# An entry holds `name`, the panel's key in the returned list (the part names
# joined by "."), `sides`, the two directions it sweeps (`from`/`to` as
# passed to predict_realloc(), and the sign the x-axis gives that direction),
# its grid `row`/`col`, and the `row_label`/`col_label` of its title.
tc_realloc_layouts <- list(
  # One panel per unordered pair, in the lower triangle of a (D-1) x (D-1)
  # grid: rows are the later part (parts[-1]), columns the earlier one
  # (parts[-D]), so excluding a part closes the triangle up. A positive
  # delta moves time to the row part.
  "one-to-one" = function(parts, labs) {
    purrr::map(utils::combn(parts, 2, simplify = FALSE), function(p) {
      list(
        name = paste(p, collapse = "."),
        sides = list(list(from = p[1], to = p[2], sign = 1),
                     list(from = p[2], to = p[1], sign = -1)),
        row = match(p[2], parts) - 1L, col = match(p[1], parts),
        row_label = labs[match(p[2], parts)],
        col_label = labs[match(p[1], parts)]
      )
    })
  },
  # One panel per part, in a wrapped near-square grid. A positive delta: the
  # part gains from the rest; a negative one: it gives to the rest.
  "one-to-all" = function(parts, labs) {
    ncol <- as.integer(ceiling(sqrt(length(parts))))
    purrr::map(seq_along(parts), function(i) {
      list(
        name = parts[[i]],
        sides = list(list(from = NULL, to = parts[[i]], sign = 1),
                     list(from = parts[[i]], to = NULL, sign = -1)),
        row = ((i - 1L) %/% ncol) + 1L,
        col = ((i - 1L) %% ncol) + 1L,
        row_label = labs[[i]], col_label = tc_rest_label()
      )
    })
  }
)

# One cell per panel of `scheme`'s layout. Each side of a panel's two-way
# sweep is predicted at the magnitudes that the scheme's own grid finds
# feasible at `.reference`; a panel with no feasible magnitude either way is
# dropped.
# Every cell has the same shape (`data`, `row`, `col`, `row_label`,
# `col_label`) and the list is named by the layout's `name`, so everything
# downstream is scheme-agnostic. `.reference`
# follows `...` so a `reference` inside `...` still reaches predict_realloc().
# lapply(), not purrr::map(): map() would re-raise predict_realloc()'s own
# aborts wrapped in an "In index: i" message.
tc_realloc_cells <- function(object, parts, labs, delta_max, n_points, type,
                             scheme = "one-to-one", ...,
                             .reference = tc_check_reference(NULL, object)) {
  half <- seq(0, delta_max, length.out = ceiling(n_points / 2) + 1)[-1]
  moves_of <- tc_realloc_schemes[[scheme]]$moves
  layout <- tc_realloc_layouts[[scheme]](parts, labs)
  cells <- lapply(layout, function(cell) {
    sides <- lapply(cell$sides, function(side) {
      moves <- moves_of(object$parts, side$from, side$to)
      feasible <- tc_realloc_grid(.reference, moves, half, scheme,
                                  on_invalid = "na")$valid
      if (!any(feasible)) {
        return(NULL)
      }
      res <- predict_realloc(object, delta = half[feasible], from = side$from,
                             to = side$to, scheme = scheme, type = type, ...)
      res$signed_delta <- side$sign * res$delta
      # The part gaining time is the move's receiver: a part, or the pooled
      # rest of a one-to-all reallocation.
      res$receiving <- moves$to
      res
    })
    df <- dplyr::bind_rows(sides)
    if (nrow(df) == 0) {
      return(NULL)
    }
    c(list(data = df[order(df$signed_delta), , drop = FALSE]),
      cell[c("row", "col", "row_label", "col_label")])
  })
  names(cells) <- purrr::map_chr(layout, "name")
  purrr::compact(cells)
}

# Build one panel (no facet) for a single cell's reallocation sweep. The cell
# is self-identified by its title, which uses a bidirectional arrow
# (row <-> column) because the sweep runs both ways; the title is centred over
# the panel, and long titles break onto two lines so adjacent panels cannot
# overlap. Axis titles are only drawn on the outer edge of the grid to avoid
# repeating them on every panel. Line and ribbon are coloured by the receiving
# part (the behaviour gaining time); significance carries no extra encoding
# and is read off the ribbon against the reference line. The colour and fill
# scales are built identically in every panel (with `limits`/`drop = FALSE`,
# and `labels` so the shared legend shows the compD display labels, not raw
# part names) so patchwork collects them into one shared legend. Every panel
# spans the full [-delta_max, delta_max] with identical breaks (a shared
# x-scale + a coord_cartesian xlim zoom), so capped sweeps simply stop short of
# the axis edge rather than rescaling the panel. The x labels are the absolute
# amounts moved: both halves read as a positive magnitude, and direction is
# carried by colour alone.
tc_realloc_panel <- function(df, y_lab, ylim, use_log, ref_line, palette,
                             legend_limits, legend_labels, delta_max, x_breaks,
                             x_title, show_y_title, show_x_title,
                             row_label, col_label) {
  arrow <- intToUtf8(8596L)
  title <- if (nchar(row_label) + nchar(col_label) > 24) {
    paste0(row_label, " ", arrow, "\n", col_label)
  } else {
    paste(row_label, arrow, col_label)
  }
  p <- ggplot2::ggplot(
    df, ggplot2::aes(x = .data$signed_delta, y = .data$estimate)
  ) +
    ggplot2::geom_ribbon(
      ggplot2::aes(
        ymin = .data$conf.low, ymax = .data$conf.high,
        fill = .data$receiving,
        group = .data$receiving
      ),
      alpha = 0.2,
      # Draw a key for every part, beyond the panel's own pair, so all
      # panels carry an identical guide and patchwork merges them into one.
      show.legend = TRUE
    ) +
    ggplot2::geom_line(
      ggplot2::aes(colour = .data$receiving, group = .data$receiving),
      show.legend = TRUE
    ) +
    ggplot2::scale_colour_manual(
      name = "Behaviour gaining time",
      values = palette, limits = legend_limits, labels = legend_labels,
      drop = FALSE
    ) +
    ggplot2::scale_fill_manual(
      name = "Behaviour gaining time",
      values = palette, limits = legend_limits, labels = legend_labels,
      drop = FALSE
    ) +
    ggplot2::scale_x_continuous(breaks = x_breaks, labels = abs) +
    ggplot2::labs(
      title = title,
      x = if (isTRUE(show_x_title)) x_title else NULL,
      y = if (isTRUE(show_y_title)) y_lab else NULL
    ) +
    ggplot2::theme_minimal(base_size = 9) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 8, hjust = 0.5)
    )

  if (!is.null(ref_line)) {
    p <- p + ggplot2::geom_hline(yintercept = ref_line, linetype = "dashed")
  }
  if (isTRUE(use_log)) {
    p <- p + ggplot2::scale_y_continuous(trans = "log10")
  }
  p + ggplot2::coord_cartesian(xlim = c(-delta_max, delta_max), ylim = ylim)
}
