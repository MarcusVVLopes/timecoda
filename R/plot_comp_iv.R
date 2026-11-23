# plot.comp_iv(): the predicted-outcome surface of a comp_iv() fit, one
# tiled (or smoothed) ternary per part triplet.

#' Visualize a composition-as-predictor model
#'
#' Draws ternary plots of the fitted composition effect of a [comp_iv()] fit,
#' one per part triplet (`choose(D, 3)` panels for a `D`-part composition), as
#' tiles of the computed values (default) or as a smoothed heatmap. The parts
#' outside the triplet are held at the reference composition, the sample's
#' closed geometric mean.
#'
#' Each panel evaluates the fit on a grid of compositions over its triplet
#' through the same `emmeans` computation as [predict_realloc()], averaged over
#' the covariates:
#' \describe{
#'   \item{`scale = "relative"` (default)}{Each grid composition is contrasted
#'     against the reference composition, so the surface is the reallocation
#'     effect measure: an odds ratio (logistic), hazard ratio (survival), rate
#'     ratio (Poisson) or mean difference (identity link). The reference
#'     composition itself maps to `1` (ratio scales) or `0` (difference
#'     scale).}
#'   \item{`scale = "absolute"`}{The predicted response (probability, mean or
#'     rate) at each grid composition.}
#' }
#' The survival engines (`"cox"`, `"coxme"`) estimate no baseline hazard, and
#' the ordinal and multinomial engines (`"clm"`, `"clmm"`, `"mblogit"`) predict
#' per outcome category. These engines support `scale = "relative"` only, as in
#' [predict_realloc()], and `scale = "absolute"` is an error for them. An
#' ordinal surface is the cumulative odds ratio versus the reference
#' composition. A multinomial surface is the relative risk ratio of the outcome
#' level named by `category` against the reference level.
#'
#' By default (`smooth = FALSE`), each tile is one computed composition filled
#' with its own predicted value, without interpolation, and [predict_realloc()]
#' gives the same values as a table.
#'
#' With `smooth = TRUE`, each panel shows a quadratic surface fitted through the
#' computed cells in the triplet's two isometric log-ratio coordinates. The fit
#' is made on the model's link scale (the log for a ratio) and mapped back. Its
#' values are drawn as finer tiles inside the same band, without contour lines.
#' The lattice refines the computed grid by the largest whole factor, up to 4,
#' that keeps it at or under 80 steps. At a `resolution` of 40 or below the
#' lattice is finer than the grid, and from 41 up it is the grid itself. The
#' computed cells stay in each panel's data under the smooth, and
#' [predict_realloc()] gives them exactly. The legend spans the computed values,
#' and a smoothed tile beyond them takes the end colour. A panel keeps its tiles
#' and warns when the quadratic cannot be fitted, which happens with fewer than
#' 12 usable computed cells, a constant surface or a degenerate arrangement of
#' compositions.
#'
#' In both pictures, the border strip of the triangle, where a part would be
#' near zero, is left blank, and a black point marks the reference composition.
#'
#' The surface is drawn only where the sample has data. A tile is kept when each
#' of its three triplet parts lies within that part's `quantiles` band over the
#' observed compositions (by default its 5th to 95th percentile). Tiles outside
#' are left blank, the axes keep the whole triangle, and the held parts, at the
#' reference composition (`ref_comp` of the fit), are not filtered. A panel
#' carries no caption; on a fit with `include_total = TRUE`, the total is held
#' at `total_ref`.
#'
#' On a fit with a composition-by-covariate interaction (a `moderator`), the
#' surface is marginal over the interacting covariate: factor levels are
#' averaged with equal weight, and numeric covariates are held at their mean. A
#' message says so once and names the covariate, and [predict_realloc()] with
#' `by` gives the effect within each level.
#'
#' @param x A `comp_iv` object from [comp_iv()].
#' @param scale Effect scale of the surface: `"relative"` (default) for the
#'   effect versus the reference composition, or `"absolute"` for the predicted
#'   response. `"absolute"` is unavailable for the survival, ordinal and
#'   multinomial engines.
#' @param resolution Number of steps along each edge of a triplet's simplex, a
#'   whole number of at least 3 (default 40). The grid holds the
#'   `choose(resolution - 1, 2)` compositions whose three shares are positive
#'   multiples of `1 / resolution` of the triplet's total (741 at the default),
#'   one tile each. Only those inside the `quantiles` band are evaluated, on a
#'   few `emmeans` reference grids per panel. A grid of fewer than about 10
#'   steps gives a very coarse mosaic and may leave a panel empty under the
#'   default `quantiles`.
#' @param quantiles The band of each part's observed distribution within which
#'   tiles are drawn: two probabilities in increasing order, by default
#'   `c(0.05, 0.95)`, the 5th to 95th percentiles. `NULL` draws the whole
#'   triangle, evaluating every grid composition.
#' @param smooth `FALSE` (default) draws one tile per computed composition;
#'   `TRUE` draws a quadratic surface fitted through the computed compositions
#'   as a smoothed heatmap.
#' @param category For a multinomial (`"mblogit"`) fit, the non-reference
#'   outcome level whose relative risk ratio against the reference level is
#'   drawn: required for such a fit, and refused for any other. `NULL`
#'   (default).
#' @param ... Present for S3 consistency. `plot()` takes no further arguments,
#'   and one given here is an error.
#'
#' @return Invisibly, a `ggplot` object built with `ggtern` for a three-part
#'   composition, otherwise a list of them named by triplet (`"MVPA.LPA.SB"`).
#'   The plots are drawn as a side effect.
#'
#' @examples
#' \donttest{
#' cd <- declare_composition(timeuse, c("MVPA", "LPA", "SB", "SPT"),
#'                           closure = 1440)
#' m <- comp_iv(cd, outcome = "obesity", model = "glm", family = "binomial")
#' plot(m)                     # odds ratio versus the reference composition
#' plot(m, scale = "absolute") # predicted probability
#' plot(m, smooth = TRUE)      # a quadratic surface through the tiles
#' if (requireNamespace("ordinal", quietly = TRUE)) {
#'   ord <- comp_iv(cd, outcome = "grade", model = "clm")
#'   plot(ord)                 # cumulative odds ratio
#' }
#' }
#'
#' @seealso [comp_iv()] to fit the model; [summary.comp_iv()] for the effect
#'   table; [predict_realloc()] for the reallocation table the surface draws.
#' @family composition as predictor
#' @export
plot.comp_iv <- function(x, scale = c("relative", "absolute"),
                         resolution = 40L, quantiles = c(0.05, 0.95),
                         smooth = FALSE, category = NULL, ...) {
  if (!inherits(x, "comp_iv")) {
    cli::cli_abort("{.arg x} must be a {.cls comp_iv} object.")
  }
  tc_check_plot_dots("plot.comp_iv", ...)
  scale <- match.arg(scale)
  meta <- tc_model_meta(x$model_type)
  if (scale == "absolute" && !meta$absolute_ok) {
    cli::cli_abort(c(
      "{.val {x$model_type}} models support {.code scale = \"relative\"} only.",
      "i" = "Their {.pkg emmeans} predictions have no single absolute-response scale ({meta$absolute_note})."
    ))
  }
  # A multinomial fit is drawn for one outcome level, named by `category`; the
  # messages name such a fit by its outcome, as optimal_composition() does.
  name <- if (meta$per_category) {
    deparse1(rlang::f_lhs(stats::formula(x$model)))
  } else {
    x$model_type
  }
  tc_check_category(rlang::set_names(list(x), name), category, "draw")
  if (!rlang::is_integerish(resolution, n = 1, finite = TRUE) ||
      resolution < 3) {
    cli::cli_abort(c(
      "{.arg resolution} must be a single whole number of at least 3.",
      "i" = "You supplied {.val {resolution}}."
    ))
  }
  if (!is.null(quantiles) &&
      (!is.numeric(quantiles) || length(quantiles) != 2 || anyNA(quantiles) ||
         quantiles[1] < 0 || quantiles[1] >= quantiles[2] || quantiles[2] > 1)) {
    cli::cli_abort(c(
      "{.arg quantiles} must be {.code NULL} or two probabilities in increasing order between 0 and 1.",
      "i" = "You supplied {.val {quantiles}}."
    ))
  }
  if (!rlang::is_bool(smooth)) {
    cli::cli_abort("{.arg smooth} must be {.code TRUE} or {.code FALSE}.")
  }
  bounds <- if (!is.null(quantiles)) {
    structure(
      purrr::map(tc_comp_df(x$comp), stats::quantile, probs = quantiles,
                 names = FALSE),
      quantiles = quantiles
    )
  }
  reference <- tc_check_reference(NULL, x)[x$parts]
  report <- if (scale == "relative") tc_report_meta(x$model_type, x$family)
  fill_name <- if (scale == "absolute") {
    "Predicted response"
  } else if (!is.null(category)) {
    paste(report$measure, "of", category)
  } else {
    report$measure
  }
  # The smooth is fitted where the model is linear in the coordinates: on the
  # log scale for a ratio, the family's link scale for an absolute prediction.
  fam <- if (scale == "absolute" && inherits(x$family, "family")) x$family
  link <- inv <- identity
  if (isTRUE(report$exponentiate)) {
    link <- log
    inv <- exp
  } else if (!is.null(fam)) {
    link <- fam$linkfun
    inv <- fam$linkinv
  }
  if (isTRUE(x$has_ilr_interaction)) {
    covs <- tc_interacting_covariates(x)
    cli::cli_inform(c(
      "The surface is marginal over the interacting covariate{?s} {.val {covs}}: factor levels are averaged with equal weight and numeric covariates are held at their mean.",
      "i" = "Use {.fn predict_realloc} with {.arg by} for the effect within each level."
    ))
  }
  out <- tc_triplet_plots(x$parts, function(tri) {
    tc_surface_panel(x, tri, reference, resolution, scale, fill_name, bounds,
                     smooth, link, inv, category = category)
  })
  tc_draw_plots(out)
}

# One panel: the emmeans surface on the grid over triplet `tri` (the other
# parts at `reference`), one tile per grid composition, with no caption.
# `bounds` keeps those inside the band: a c(lo, hi) per part (plot.comp_iv()'s
# `quantiles`), or a function of the full composition matrix giving the rows
# to keep (a search's region). `smooth` draws a quadratic through them on a
# finer lattice instead. `control`, when given, replaces `reference` as the
# composition a relative surface is contrasted with; `category` names a
# multinomial fit's level. `kept_by` names what set `bounds` in the warnings:
# "quantiles", "region" or "candidates" (a search's per-part range); `NULL`
# reads a function as a region and anything else as the `quantiles` band.
tc_surface_panel <- function(object, tri, reference, resolution, scale,
                             fill_name, bounds = NULL, smooth = FALSE, link,
                             inv, at = NULL, control = NULL, category = NULL,
                             kept_by = NULL) {
  if (is.null(kept_by)) {
    kept_by <- if (is.function(bounds)) "region" else "quantiles"
  }
  reference <- reference[object$parts]
  held <- setdiff(object$parts, tri)
  total3 <- object$total - sum(reference[held])
  # The triplet's shares on the i/r lattice, those inside the band.
  lattice <- function(r) {
    ii <- seq_len(r - 1)
    grid <- expand.grid(i = ii, j = ii)
    grid <- grid[grid$i + grid$j < r, , drop = FALSE]
    a <- grid$i / r
    b <- grid$j / r
    shares <- cbind(a, b, 1 - a - b) * total3
    if (is.null(bounds)) return(shares)
    if (is.function(bounds)) {
      # A region's own test, on the full compositions: the triplet's shares
      # with the other parts at `reference`.
      full <- matrix(reference, nrow = nrow(shares), ncol = length(reference),
                     byrow = TRUE, dimnames = list(NULL, names(reference)))
      full[, tri] <- shares
      return(shares[bounds(full), , drop = FALSE])
    }
    lo <- purrr::map_dbl(bounds[tri], 1)
    hi <- purrr::map_dbl(bounds[tri], 2)
    inside <- sweep(shares, 2, lo, ">=") & sweep(shares, 2, hi, "<=")
    shares[rowSums(inside) == ncol(shares), , drop = FALSE]
  }
  shares <- lattice(resolution)
  pred <- numeric()
  # The two warnings below name what kept the cells, a search's region or
  # candidates or the `quantiles` band, and the remedy the calling plot offers.
  remedy <- if (kept_by == "quantiles") {
    "Raise {.arg resolution} or widen {.arg quantiles}."
  } else {
    "Raise {.arg resolution}."
  }
  if (nrow(shares) == 0) {
    cli::cli_warn(c(
      switch(kept_by,
        quantiles = "No grid composition of {.val {tri}} lies within the {.arg quantiles} band; the panel is blank.",
        region = "No grid composition of {.val {tri}} lies within the search's region; the panel is blank.",
        candidates = "No grid composition of {.val {tri}} lies within the range of the search's candidates; the panel is blank."
      ),
      "i" = remedy
    ))
  } else {
    # Every cell on a few emmeans grids, the way the Goldilocks search predicts
    # its candidates; a relative surface is contrasted with `control`, else
    # the reference composition.
    comps <- matrix(reference, nrow = nrow(shares), ncol = length(reference),
                    byrow = TRUE, dimnames = list(NULL, names(reference)))
    comps[, tri] <- shares
    from <- if (!is.null(control)) control[object$parts] else reference
    cand <- list(z = tc_pivot_coords(comps, pivot = object$pivot))
    pred <- tc_optimal_predict(object, cand, at, scale, from,
                               category)$response
  }
  d <- tibble::tibble(
    triplet = paste(tri, collapse = "."),
    xpart = tri[1], ypart = tri[2], zpart = tri[3],
    x = shares[, 1], y = shares[, 2], z = shares[, 3], .pred = as.numeric(pred)
  )
  # The triplet's own two ILR coordinates, as columns u and v.
  uv <- function(shares) {
    df <- stats::setNames(as.data.frame(shares), c("x", "y", "z"))
    stats::setNames(ilr_transform(df, c("x", "y", "z"))[c("ilr1", "ilr2")],
                    c("u", "v"))
  }
  # A quadratic on the link scale (`link`, `inv`, required when smoothing)
  # through the finite values: it needs twice its six parameters, varying
  # values and a full-rank fit.
  fit <- NULL
  if (smooth) {
    value <- link(d$.pred)
    ok <- is.finite(value)
    if (sum(ok) >= 12 && diff(range(value[ok])) > 0) {
      fit <- stats::lm(value ~ stats::poly(u, v, degree = 2),
                       data = cbind(uv(d[ok, c("x", "y", "z")]),
                                    value = value[ok]))
    }
  }
  if (smooth && nrow(d) > 0 && (is.null(fit) || fit$rank != 6L)) {
    inside <- switch(kept_by, quantiles = "inside the band",
                     region = "inside the region",
                     candidates = "inside the candidates' range")
    cli::cli_warn(c(
      "The quadratic surface cannot be fitted for {.val {tri}} (fewer than 12 usable grid compositions {inside}, a constant surface, or a degenerate arrangement of compositions); drawing tiles.",
      "i" = remedy
    ))
  }
  smooth <- !is.null(fit) && fit$rank == 6L
  tile_at <- d
  tile_res <- resolution
  if (smooth) {
    # A whole multiple of the computed grid, so its lattice contains every
    # computed cell: up to four times finer while at or under 80 steps.
    tile_res <- resolution * max(1L, min(4L, 80L %/% resolution))
    fine <- lattice(tile_res)
    tile_at <- tibble::tibble(
      x = fine[, 1], y = fine[, 2], z = fine[, 3],
      .pred = inv(stats::predict(fit, uv(fine)))
    )
  }
  # Each tile is the hexagon whose corners are the centroids of the six
  # lattice triangles around its composition, listed in angular order.
  hex <- rbind(c(2, -1, -1), c(1, 1, -2), c(-1, 2, -1), c(-2, 1, 1),
               c(-1, -1, 2), c(1, -2, 1)) * total3 / (3 * tile_res)
  k <- rep(seq_len(nrow(tile_at)), each = 6)
  tiles <- tibble::tibble(
    cell = k, .pred = tile_at$.pred[k],
    x = tile_at$x[k] + hex[, 1], y = tile_at$y[k] + hex[, 2],
    z = tile_at$z[k] + hex[, 3]
  )
  # The corners are titled by the parts' display labels.
  corners <- tc_part_labels(tri, object$comp$labels)
  # A smooth's legend spans the computed values; its overshoot is clamped.
  fill_scale <- if (smooth) {
    ggplot2::scale_fill_viridis_c(name = fill_name, limits = range(d$.pred[ok]),
                                  oob = scales::squish)
  } else {
    ggplot2::scale_fill_viridis_c(name = fill_name)
  }
  ggtern::ggtern(d, ggplot2::aes(x = .data$x, y = .data$y, z = .data$z)) +
    ggplot2::geom_polygon(
      ggplot2::aes(group = .data$cell, fill = .data$.pred),
      data = tiles, colour = NA
    ) +
    ggplot2::geom_point(
      data = tibble::tibble(x = reference[[tri[1]]], y = reference[[tri[2]]],
                            z = reference[[tri[3]]]),
      colour = "black", size = 1.8
    ) +
    fill_scale +
    ggplot2::labs(x = corners[1], y = corners[2], z = corners[3])
}
