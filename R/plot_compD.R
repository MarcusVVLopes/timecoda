# `plot.compD()` packages the exploratory figures of the "composition-variable"
# vignette (per-part distributions and the ternary plot) into one dispatched
# method. All figures derive from the imputed composition vector stored in the
# object, so they are consistent with `summary.compD()` and free of zeroes.

#' Visualize a compositional variable
#'
#' Draws per-part density curves and boxplots of the untransformed parts of a
#' [declare_composition()] object, the ternary plots of the composition and the
#' heatmap of its variation matrix, all from the object's imputed composition.
#' When the object carries flagged compositional outliers
#' ([declare_composition()]'s `detect_outliers`), two further panels show the
#' outlier distances and the flagged observations in the simplex.
#'
#' @param x A `compD` object from [declare_composition()].
#' @param histogram Logical; include per-part density curves
#'   ([ggplot2::geom_density()], faceted by part). Default `TRUE`.
#' @param boxplot Logical; include per-part boxplots
#'   ([ggplot2::geom_boxplot()]). Default `TRUE`.
#' @param ternary Logical; include the ternary plot(s) ([plot_ternary()]).
#'   Default `TRUE`. Set `FALSE` for a two-part composition, which has no
#'   ternary representation; this also drops the outlier ternary.
#' @param variation Logical; include the variation-matrix heatmap
#'   ([plot_variation()]). Default `TRUE`.
#' @param outliers Logical; include the outlier panels. Default `TRUE`. The
#'   panels appear only when `x` carries flagged outliers, so a composition
#'   declared without `detect_outliers` is unaffected.
#' @param ... Further arguments passed to [plot_ternary()] (e.g., `group`).
#'   They do not reach the outlier panels, which always describe the whole
#'   sample (see Value).
#'
#' @return Invisibly, a named list of the requested plots, with elements
#'   `histogram`, `boxplot`, `ternary` (a single `ggtern` plot or, for more than
#'   three parts, a named list of `choose(D, 3)` panels) and `variation` (the
#'   variation-matrix heatmap). When `x` carries flagged outliers, two further
#'   elements follow: `outliers`, the Mahalanobis distance of each observation
#'   against the cutoff, and, unless `ternary = FALSE`, `outlier_ternary`, the
#'   ternary panel(s) with the flagged observations highlighted. Both describe
#'   the whole sample and are never split by `group`. Calling `plot()` draws the
#'   plots as a side effect.
#'
#' @examples
#' \donttest{
#' cd <- declare_composition(timeuse, c("MVPA", "LPA", "SB", "SPT"),
#'                           closure = 1440)
#' plots <- plot(cd)
#' plots$boxplot
#' }
#'
#' @seealso [declare_composition()] to build the `compD` object;
#'   [summary.compD()] for the descriptive tables; [plot_ternary()] and
#'   [plot_variation()] for the ternary and variation-matrix figures alone.
#' @family declare and describe composition
#' @export
plot.compD <- function(x, histogram = TRUE, boxplot = TRUE, ternary = TRUE,
                       variation = TRUE, outliers = TRUE, ...) {
  if (!inherits(x, "compD")) {
    cli::cli_abort("{.arg x} must be a {.cls compD} object.")
  }
  comp_df <- tc_comp_df(x)
  parts <- x$parts
  labs <- tc_part_labels(parts, x$labels)

  out <- list()

  if (isTRUE(histogram) || isTRUE(boxplot)) {
    long <- purrr::map(parts, function(p) {
      tibble::tibble(part = p, value = comp_df[[p]])
    }) |>
      purrr::list_rbind()
    # Display labels drive the facet/axis text; raw parts set the order.
    long$part <- factor(long$part, levels = parts, labels = labs)

    if (isTRUE(histogram)) {
      out$histogram <- ggplot2::ggplot(
        long,
        ggplot2::aes(x = .data$value, colour = .data$part, fill = .data$part)
      ) +
        ggplot2::geom_density(alpha = 0.3) +
        ggplot2::facet_wrap(ggplot2::vars(.data$part)) +
        ggplot2::labs(x = NULL, y = "Density") +
        ggplot2::theme_classic() +
        ggplot2::theme(legend.position = "none")
    }
    if (isTRUE(boxplot)) {
      out$boxplot <- ggplot2::ggplot(
        long,
        ggplot2::aes(x = .data$part, y = .data$value, colour = .data$part)
      ) +
        ggplot2::geom_boxplot() +
        ggplot2::labs(x = NULL, y = NULL) +
        ggplot2::theme_classic() +
        ggplot2::theme(legend.position = "none")
    }
  }

  if (isTRUE(ternary)) {
    out$ternary <- tc_build_ternary(x, ...)
  }

  if (isTRUE(variation)) {
    out$variation <- tc_variation_heatmap(comp_df, parts, labs)
  }

  # The outlier panels never take `...` and so are never split by `group`:
  # compositional outlyingness is measured against the covariance of the
  # entire sample, so a per-group display would show structure the detection
  # never used.
  if (isTRUE(outliers) && tc_has_outliers(x)) {
    out <- c(out, tc_outlier_panels(comp_df, parts, labs, x$data_raw,
                                    isTRUE(ternary)))
  }

  # Draw every panel so plot(cd) displays like any base plot, then return the
  # objects (invisibly) so they remain editable.
  tc_draw_plots(out)
}

# The outlier distance panel and, when `ternary`, the outlier ternary panel(s).
# Never split by `group`: outlyingness is measured against the whole sample.
tc_outlier_panels <- function(comp_df, parts, labels, data_raw, ternary) {
  # "Typical" first so it takes the base colour; the same scales in both
  # panels, so an observation reads the same way in each.
  status <- factor(
    ifelse(data_raw$.outlier, "Outlier", "Typical"),
    levels = c("Typical", "Outlier")
  )
  scales <- list(
    ggplot2::scale_colour_manual(
      values = c(Typical = "grey60", Outlier = "#D55E00")
    ),
    ggplot2::scale_shape_manual(values = c(Typical = 16, Outlier = 17)),
    ggplot2::labs(colour = NULL, shape = NULL)
  )

  # Mahalanobis distance against observation index, with the cutoff as a
  # dashed rule: a rebuild of robCompositions' own (base-graphics) diagnostic.
  panels <- list(
    outliers = ggplot2::ggplot(
      tibble::tibble(
        row = seq_len(nrow(data_raw)),
        distance = data_raw$.outlier_distance,
        status = status
      ),
      ggplot2::aes(
        x = .data$row, y = .data$distance,
        colour = .data$status, shape = .data$status
      )
    ) +
      ggplot2::geom_point(alpha = 0.7) +
      ggplot2::geom_hline(
        yintercept = data_raw$.outlier_limit[1], linetype = 2, colour = "grey30"
      ) +
      ggplot2::labs(x = "Observation", y = "Robust Mahalanobis distance") +
      scales +
      ggplot2::theme_classic()
  )

  # Gated on `ternary`, which switches off every ternary panel: it is how a
  # two-part composition is plotted at all. The distance panel has no limit.
  if (ternary) {
    data <- comp_df
    data$status <- status
    panels$outlier_ternary <- tc_triplet_plots(parts, function(tri) {
      ax <- unname(stats::setNames(labels, parts)[tri])
      ggtern::ggtern(
        data,
        ggplot2::aes(
          x = .data[[tri[1]]], y = .data[[tri[2]]], z = .data[[tri[3]]],
          colour = .data$status, shape = .data$status
        )
      ) +
        ggplot2::geom_point(alpha = 0.6) +
        ggplot2::labs(x = ax[1], y = ax[2], z = ax[3]) +
        scales
    })
  }
  panels
}
