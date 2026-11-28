# Module 4: ternary plots and the variation heatmap.
#
# Every ternary plot in the package goes through `ggtern`. Panels are derived
# from the parts, one per triplet via the shared `tc_triplet_plots()` in
# `utils_plot.R`, never from a hardcoded layout.

#' Ternary plot(s) of a composition
#'
#' Draws ternary (three-part) plots of a time-use composition via `ggtern`,
#' one panel per triplet of parts (`choose(D, 3)` panels for a `D`-part
#' composition). Group colouring, group-mean points, crosshairs at the group
#' means and confidence regions are optional.
#'
#' @param x A `compD` object from [declare_composition()]. The ternary
#'   coordinates are always the object's closed (imputed) composition. A long
#'   part label is drawn in full where the figure has room beside the triangle.
#'   In a small figure it can run past the edge. A shorter label or a larger
#'   figure keeps it inside.
#' @param group Optional name of a column of the data originally given to
#'   [declare_composition()], used to colour the points and to compute
#'   per-group means.
#' @param mean_points Add the group (or overall) closed geometric mean as
#'   points.
#' @param crosshairs Add crosshairs at the mean composition(s).
#' @param conf_level One or more confidence levels, each drawn as a
#'   `geom_confidence_tern()` region; set `NULL` to omit them.
#' @param theme Optional `ggtern` theme added to every panel (e.g.,
#'   `ggtern::theme_bw()`). Default `NULL` keeps the `ggtern` default. Complete
#'   themes include `theme_bw()`, `theme_light()`, `theme_dark()`,
#'   `theme_minimal()`, `theme_classic()`, `theme_void()`, `theme_rgbw()`,
#'   `theme_rgbg()`, `theme_bluelight()` and `theme_tropical()`. Add-on themes
#'   such as `theme_hidegrid()`, `theme_showarrows()`, `theme_noticks()` or
#'   `theme_latex()` can be combined with them by `+` (e.g.,
#'   `ggtern::theme_bw() + ggtern::theme_hidegrid()`). All are documented in
#'   the `ggtern` package; browse them with `help(package = "ggtern")`
#'   (functions whose names start with `theme_`) or read `?ggtern::theme_bw`.
#'
#' @return Invisibly, a single `ggtern`/`ggplot` object when there is one
#'   triplet, otherwise a named list of them (names like `"MVPA.LPA.SB"`). The
#'   call draws the plot(s) as a side effect.
#'
#' @examples
#' \donttest{
#' cd <- declare_composition(timeuse, c("MVPA", "LPA", "SB", "SPT"),
#'                           closure = 1440)
#' plot_ternary(cd, group = "sex")
#' plot_ternary(cd, group = "sex", theme = ggtern::theme_bw())
#' }
#'
#' @seealso [declare_composition()] for the composition object;
#'   [plot.compD()], which delegates to this function.
#' @family visualization
#' @export
plot_ternary <- function(x, group = NULL, mean_points = TRUE,
                         crosshairs = TRUE, conf_level = 0.95, theme = NULL) {
  if (!inherits(x, "compD")) {
    cli::cli_abort("{.arg x} must be a {.cls compD} object.")
  }
  plots <- tc_build_ternary(x, group, mean_points, crosshairs, conf_level,
                            theme)
  tc_draw_plots(plots)
}

#' Heatmap of the compositional variation matrix
#'
#' Draws the compositional variation matrix as a tile heatmap, shading each
#' cell by the pairwise log-ratio variance `var(log(x_i / x_j))`, zero by
#' construction on the diagonal. Small values mark highly proportional parts
#' (parts that move together); with the default colours they are drawn darker
#' and large values lighter. The axes carry the part labels of the
#' [declare_composition()] object.
#'
#' @param x A `compD` object from [declare_composition()].
#' @param digits Decimal places for the value labels drawn in each cell.
#'   Default `2`.
#' @param show_values Logical; print the rounded variance in each cell.
#'   Default `TRUE`.
#' @param low,high Colours for the low (dark, highly proportional) and high
#'   (light) ends of the gradient. Defaults give a dark-to-light blue scale.
#' @param theme Optional `ggplot2` theme added over the default heatmap theme
#'   (e.g., `ggplot2::theme_bw()` or
#'   `ggplot2::theme(legend.position = "bottom")`). Default `NULL` adds none.
#'   Complete themes are listed in the `ggplot2` reference under "Themes"
#'   (`?ggplot2::theme_bw`). `plot_ternary()` takes a separate `ggtern` theme.
#'
#' @return Invisibly, a `ggplot` object. The call draws the heatmap as a side
#'   effect.
#'
#' @examples
#' cd <- declare_composition(timeuse, c("MVPA", "LPA", "SB", "SPT"),
#'                           closure = 1440)
#' plot_variation(cd)
#' plot_variation(cd, theme = ggplot2::theme_bw())
#'
#' @family visualization
#' @seealso [summary.compD()] for the variation matrix as a table.
#' @export
plot_variation <- function(x, digits = 2, show_values = TRUE,
                           low = "#132B43", high = "#EAF3FB", theme = NULL) {
  if (!inherits(x, "compD")) {
    cli::cli_abort("{.arg x} must be a {.cls compD} object.")
  }
  p <- tc_variation_heatmap(
    tc_comp_df(x), x$parts, tc_part_labels(x$parts, x$labels),
    digits, show_values, low, high, theme
  )
  tc_draw_plots(p)
}

# Internal helpers ------------------------------------------------------------

# Build the variation-matrix heatmap. Lower variances (more proportional parts)
# get darker fills; cells optionally carry the rounded value with a contrasting
# text colour.
tc_variation_heatmap <- function(comp_df, parts, labels, digits = 2,
                                 show_values = TRUE, low = "#132B43",
                                 high = "#EAF3FB", theme = NULL) {
  tc_check_theme(theme, "ggplot2")
  vm <- variation_matrix(comp_df, parts)
  long <- expand.grid(row = parts, col = parts, KEEP.OUT.ATTRS = FALSE,
                      stringsAsFactors = FALSE)
  long$value <- vm[cbind(long$row, long$col)]
  # Reverse the rows so the zero diagonal runs top-left to bottom-right, as in a
  # printed matrix; display labels drive the tick text.
  long$row <- factor(long$row, levels = rev(parts), labels = rev(labels))
  long$col <- factor(long$col, levels = parts, labels = labels)
  # Text colour follows the luminance of each cell's own fill, computed with
  # the same Lab interpolation ggplot2's gradient uses, so the values stay
  # readable whichever way round `low` and `high` are given.
  fill_hex <- scales::seq_gradient_pal(low, high, "Lab")(
    scales::rescale(long$value)
  )
  chan <- grDevices::col2rgb(fill_hex) / 255
  luminance <- 0.299 * chan[1, ] + 0.587 * chan[2, ] + 0.114 * chan[3, ]
  long$txt_col <- ifelse(luminance < 0.5, "white", "grey10")

  p <- ggplot2::ggplot(
    long,
    ggplot2::aes(x = .data$col, y = .data$row, fill = .data$value)
  ) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.5) +
    ggplot2::scale_fill_gradient(low = low, high = high, name = "Variation") +
    ggplot2::coord_fixed() +
    ggplot2::labs(
      x = NULL, y = NULL,
      title = "Compositional variation matrix",
      subtitle = "Pairwise log-ratio variance, var(log(x_i / x_j))"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      panel.grid = ggplot2::element_blank()
    )

  if (isTRUE(show_values)) {
    p <- p +
      ggplot2::geom_text(
        ggplot2::aes(
          label = formatC(.data$value, format = "f", digits = digits),
          colour = .data$txt_col
        ),
        size = 3, show.legend = FALSE
      ) +
      ggplot2::scale_colour_identity()
  }
  if (is.null(theme)) p else p + theme
}

# Build the ternary plot(s) without drawing them, so internal callers
# (plot.compD) can embed the result and control display themselves.
tc_build_ternary <- function(x, group = NULL, mean_points = TRUE,
                             crosshairs = TRUE, conf_level = 0.95,
                             theme = NULL) {
  tc_check_theme(theme, "ggtern")
  if (!is.null(conf_level) &&
      (!is.numeric(conf_level) || length(conf_level) == 0 ||
       anyNA(conf_level) || any(conf_level <= 0 | conf_level >= 1))) {
    cli::cli_abort(
      "{.arg conf_level} must be one or more numbers in (0, 1), or {.code NULL}.",
      call = rlang::caller_env()
    )
  }
  data <- tc_comp_df(x)
  parts <- x$parts
  # Axis titles carry the display labels, as every other plot builder in the
  # package does; keyed by part so a triplet can be looked up directly.
  axis_labels <- stats::setNames(tc_part_labels(parts, x$labels), parts)
  if (!is.null(group)) {
    if (!is.character(group) || length(group) != 1 ||
        !group %in% names(x$data_raw)) {
      cli::cli_abort(
        "{.arg group} must name a single column of the data given to {.fn declare_composition}."
      )
    }
    data[[group]] <- as.factor(x$data_raw[[group]])
  }
  tc_triplet_plots(parts, function(tri) {
    p <- tc_one_ternary(data, tri, axis_labels, group, mean_points, crosshairs,
                        conf_level)
    if (is.null(theme)) p else p + theme
  })
}

# The `theme` argument of the plot functions must be a theme object from `pkg`
# (ggtern for ternary plots, ggplot2 for the heatmap); NULL keeps the default.
tc_check_theme <- function(theme, pkg, call = rlang::caller_env()) {
  if (!is.null(theme) && !inherits(theme, "theme")) {
    cli::cli_abort(
      "{.arg theme} must be a {.pkg {pkg}} theme object such as {.code {pkg}::theme_bw()}.",
      call = call
    )
  }
}

# One ternary panel for a triplet of parts: the point cloud, an optional
# confidence region, and the group (or overall) geometric-mean composition.
# `axis_labels` is keyed by part name; `data` still uses the raw names.
tc_one_ternary <- function(data, tri, axis_labels, group, mean_points,
                           crosshairs, conf_level) {
  ax <- unname(axis_labels[tri])
  aes_args <- list(
    x = rlang::sym(tri[1]),
    y = rlang::sym(tri[2]),
    z = rlang::sym(tri[3])
  )
  if (!is.null(group)) {
    aes_args$colour <- rlang::sym(group)
    aes_args$group <- rlang::sym(group)
  }
  p <- ggtern::ggtern(data, do.call(ggplot2::aes, aes_args)) +
    ggplot2::geom_point(alpha = 0.4)

  if (!is.null(conf_level)) {
    p <- p + ggtern::geom_confidence_tern(
      breaks = conf_level, linewidth = 0.5, show.legend = FALSE
    )
  }
  if (!mean_points && !crosshairs) {
    return(p + ggplot2::labs(x = ax[1], y = ax[2], z = ax[3]))
  }

  # Geometric-mean composition per group level (or overall), on this triplet.
  if (is.null(group)) {
    means <- geom_mean_comp(data, tri)
  } else {
    by_level <- split(data, droplevels(data[[group]]), drop = TRUE)
    means <- purrr::imap(by_level, function(d, g) {
      m <- geom_mean_comp(d, tri)
      m[[group]] <- factor(g, levels = names(by_level))
      m
    }) |>
      purrr::list_rbind()
  }

  if (crosshairs) {
    # geom_crosshair_tern mishandles a grouped (multi-row, colour-mapped)
    # layer, so each mean row gets its own layer with the colour as a static
    # parameter. Colours are looked up by group level, not row position:
    # ggplot2's discrete scale drops unused factor levels, so sizing the
    # palette by nlevels() would shift every hue after the missing one.
    cols <- if (is.null(group)) {
      rep("black", nrow(means))
    } else {
      levs <- levels(means[[group]])
      pal <- stats::setNames(scales::hue_pal()(length(levs)), levs)
      unname(pal[as.character(means[[group]])])
    }
    cross_aes <- ggplot2::aes(
      x = .data[[tri[1]]], y = .data[[tri[2]]], z = .data[[tri[3]]]
    )
    p <- p + purrr::map(seq_len(nrow(means)), function(i) {
      ggtern::geom_crosshair_tern(
        data = means[i, , drop = FALSE], mapping = cross_aes,
        inherit.aes = FALSE, linetype = 2, linewidth = 0.5,
        colour = cols[i]
      )
    })
  }
  if (mean_points) {
    p <- p + ggplot2::geom_point(data = means, size = 3)
  }
  p + ggplot2::labs(x = ax[1], y = ax[2], z = ax[3])
}
