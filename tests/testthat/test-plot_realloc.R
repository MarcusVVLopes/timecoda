# plot_realloc(): panel layout, palette, legend, schemes and the exclude filter.

parts <- c("MVPA", "LPA", "SB", "SPT")

test_that("plot_realloc returns one panel per unordered pair (true triangle)", {
  cd <- declare_composition(timeuse, parts, closure = 1440)
  fit <- comp_iv(cd, "obesity", model = "glm", family = "binomial")
  p <- plot_realloc(fit, delta_max = 20)
  expect_s3_class(p, "realloc_plots")
  expect_s3_class(p$figure, "patchwork")
  expect_s3_class(p$figure, "ggplot")
  # One panel per unordered pair -- no mirrored duplicate, no empty notch.
  expect_equal(length(p$figure), choose(length(parts), 2))
})

test_that("plot_realloc panel count tracks the number of parts", {
  set.seed(42)
  for (d in c(3, 5)) {
    pn <- paste0("p", seq_len(d))
    sim <- as.data.frame(matrix(runif(200 * d, 5, 50), ncol = d))
    names(sim) <- pn
    sim$y <- rnorm(200)
    cd <- declare_composition(sim, pn)
    fit <- comp_iv(cd, "y", model = "lm")
    p <- plot_realloc(fit, delta_max = 3, type = "relative")
    expect_equal(length(p$figure), choose(d, 2))
  }
})

test_that("plot_realloc aborts when no reallocation is feasible", {
  cd <- declare_composition(timeuse, parts, closure = 1440)
  fit <- comp_iv(cd, "obesity", model = "glm", family = "binomial")
  # n_points = 3 -> smallest sweep point is delta_max / 2. With delta_max well
  # above the 1440-minute total, that already exceeds every part's reference
  # share, so no pair has a feasible delta in either direction.
  expect_error(
    plot_realloc(fit, delta_max = 3000, n_points = 3),
    "No reallocation is feasible"
  )
})

test_that("plot_realloc labels carry through to panel titles", {
  cd <- declare_composition(timeuse, parts, closure = 1440,
                            labels = c(MVPA = "Exercise"))
  fit <- comp_iv(cd, "obesity", model = "glm", family = "binomial")
  p <- plot_realloc(fit, delta_max = 20)$figure
  titles <- vapply(p$patches$plots, function(pl) pl$labels$title, character(1))
  titles <- c(titles, p$labels$title)
  expect_true(any(grepl("Exercise", titles)))
})

# Distinct line types drawn in a patchwork's constituent panels (the GeomLine
# sweep line only -- not the ribbon or the dashed reference hline).
panel_line_linetypes <- function(p) {
  panels <- p$figure$patches$plots
  lt <- unlist(lapply(panels, function(pl) {
    b <- ggplot2::ggplot_build(pl)
    is_line <- vapply(
      pl$layers, function(l) inherits(l$geom, "GeomLine"), logical(1)
    )
    b$data[[which(is_line)[1]]]$linetype
  }))
  unique(lt)
}

test_that("plot_realloc draws one line type (no significance encoding)", {
  cd <- declare_composition(timeuse, parts, closure = 1440)
  fit <- comp_iv(cd, "obesity", model = "glm", family = "binomial")
  p <- plot_realloc(fit, delta_max = 90, n_points = 25)
  expect_s3_class(p, "realloc_plots")
  expect_s3_class(p$figure, "patchwork")
  expect_equal(length(panel_line_linetypes(p)), 1)
})

test_that("plot_realloc sweep lines cross the reference at zero minutes", {
  cd <- declare_composition(timeuse, parts, closure = 1440)
  fit <- comp_iv(cd, "obesity", model = "glm", family = "binomial")
  p <- plot_realloc(fit, delta_max = 60, n_points = 9)
  pl <- p$panels[["MVPA.LPA"]]
  b <- ggplot2::ggplot_build(pl)
  is_line <- vapply(
    pl$layers, function(l) inherits(l$geom, "GeomLine"), logical(1)
  )
  ld <- b$data[[which(is_line)[1]]]
  # Both receiving sides carry an anchor at signed_delta = 0 whose y is the
  # reference (OR = 1, i.e. 0 after the log10 scale transform).
  anchors <- ld[ld$x == 0, ]
  expect_equal(nrow(anchors), 2)
  expect_equal(unique(anchors$y), 0, tolerance = 1e-10)
})

test_that("plot_realloc default palette is viridis over every part", {
  cd <- declare_composition(timeuse, parts, closure = 1440)
  fit <- comp_iv(cd, "obesity", model = "glm", family = "binomial")
  colours <- function(p, keys) {
    unname(p$figure$patches$plots[[1]]$scales$get_scales("colour")$map(keys))
  }
  full <- plot_realloc(fit, delta_max = 20, n_points = 3)
  expect_equal(colours(full, parts), scales::viridis_pal()(length(parts)))
  # Excluding a part leaves the remaining parts' colours unchanged.
  kept <- setdiff(parts, "MVPA")
  cut <- plot_realloc(fit, delta_max = 20, n_points = 3, exclude = "MVPA")
  expect_equal(colours(cut, kept), colours(full, kept))
})

test_that("plot_realloc rejects a palette missing a part", {
  cd <- declare_composition(timeuse, parts, closure = 1440)
  fit <- comp_iv(cd, "obesity", model = "glm", family = "binomial")
  bad <- c(MVPA = "red", LPA = "blue", SB = "green")
  expect_error(
    plot_realloc(fit, delta_max = 20, palette = bad),
    "missing a colour"
  )
})

test_that("plot_realloc legend = FALSE hides guides", {
  cd <- declare_composition(timeuse, parts, closure = 1440)
  fit <- comp_iv(cd, "obesity", model = "glm", family = "binomial")
  p <- plot_realloc(fit, delta_max = 20, legend = FALSE)
  # The whole-patchwork `& theme()` sets legend.position on every panel.
  positions <- vapply(
    p$figure$patches$plots,
    function(pl) {
      lp <- pl$theme$legend.position
      if (is.null(lp)) "" else as.character(lp)[1]
    },
    character(1)
  )
  expect_true(all(positions == "none"))
})

test_that("plot_realloc colours the line by the receiving part", {
  cd <- declare_composition(timeuse, parts, closure = 1440)
  fit <- comp_iv(cd, "obesity", model = "glm", family = "binomial")
  pal <- c(MVPA = "#111111", LPA = "#222222", SB = "#333333", SPT = "#444444")
  p <- plot_realloc(fit, delta_max = 20, palette = pal)
  cols <- unlist(lapply(p$figure$patches$plots, function(pl) {
    b <- ggplot2::ggplot_build(pl)
    is_line <- vapply(
      pl$layers, function(l) inherits(l$geom, "GeomLine"), logical(1)
    )
    b$data[[which(is_line)[1]]]$colour
  }))
  # Every drawn colour comes from the supplied palette (receiving-part colours).
  expect_true(all(tolower(unique(cols)) %in% tolower(pal)))
})

test_that("plot_realloc panels share the full [-delta_max, delta_max] x-axis", {
  cd <- declare_composition(timeuse, parts, closure = 1440)
  fit <- comp_iv(cd, "obesity", model = "glm", family = "binomial")
  delta_max <- 60
  p <- plot_realloc(fit, delta_max = delta_max, n_points = 9)
  # Every panel's rendered coordinate range spans the requested delta window
  # (via coord_cartesian's xlim zoom), regardless of how far its own (capped)
  # sweep actually reaches. The built x.range carries the coord's expansion.
  ranges <- lapply(p$figure$patches$plots, function(pl) {
    ggplot2::ggplot_build(pl)$layout$panel_params[[1]]$x.range
  })
  for (rg in ranges) {
    expect_lte(rg[1], -delta_max)
    expect_gte(rg[2], delta_max)
  }
})

test_that("plot_realloc legend uses the compD display labels", {
  cd <- declare_composition(
    timeuse, parts, closure = 1440,
    labels = c(MVPA = "MVPA (min/day)", LPA = "LPA (min/day)",
               SB = "SB (min/day)", SPT = "Sleep (min/day)")
  )
  fit <- comp_iv(cd, "obesity", model = "glm", family = "binomial")
  p <- plot_realloc(fit, delta_max = 20, n_points = 9)
  pl <- p$panels[["MVPA.LPA"]]
  colour_scale <- pl$scales$get_scales("colour")
  expect_equal(
    colour_scale$labels,
    tc_part_labels(parts, cd$labels)
  )
  expect_true("MVPA (min/day)" %in% colour_scale$labels)
})

# --- Scenario (b): Cox rejects type = "absolute" in plot_realloc -------------

test_that("plot_realloc rejects type = 'absolute' for a Cox fit with an informative error", {
  small_parts <- c("MVPA", "SB", "SPT")
  cd_cox <- declare_composition(timeuse, small_parts, closure = 1440)
  cox_fit <- comp_iv(cd_cox, time = "time", event = "event", model = "cox")
  expect_error(
    plot_realloc(cox_fit, delta_max = 20, type = "absolute"),
    "absolute|relative"
  )
})

test_that("plot_realloc scheme = 'one-to-all' gives one panel per part", {
  cd <- declare_composition(timeuse, parts, closure = 1440)
  fit <- comp_iv(cd, "obesity", model = "glm", family = "binomial")
  p <- plot_realloc(fit, delta_max = 20, n_points = 9, scheme = "one-to-all")
  expect_s3_class(p, "realloc_plots")
  expect_s3_class(p$figure, "patchwork")
  expect_equal(length(p$figure), length(parts))
})

test_that("plot_realloc one-to-all panel count tracks the number of parts", {
  set.seed(42)
  for (d in c(3, 5)) {
    pn <- paste0("p", seq_len(d))
    sim <- as.data.frame(matrix(runif(200 * d, 5, 50), ncol = d))
    names(sim) <- pn
    sim$y <- rnorm(200)
    cd <- declare_composition(sim, pn)
    fit <- comp_iv(cd, "y", model = "lm")
    p <- plot_realloc(fit, delta_max = 3, n_points = 9, scheme = "one-to-all")
    expect_equal(length(p$figure), d)
  }
})

test_that("plot_realloc one-to-all adds a neutral '(rest)' legend key", {
  cd <- declare_composition(timeuse, parts, closure = 1440)
  fit <- comp_iv(cd, "obesity", model = "glm", family = "binomial")
  p <- plot_realloc(fit, delta_max = 20, n_points = 9, scheme = "one-to-all")
  sc <- p$panels[["MVPA"]]$scales$get_scales("colour")
  expect_equal(sc$limits, c(parts, "(rest)"))
  expect_true("All other parts" %in% sc$labels)
  # The pooled side is drawn in neutral grey, not a viridis behaviour colour.
  cols <- unlist(lapply(p$figure$patches$plots, function(pl) {
    b <- ggplot2::ggplot_build(pl)
    is_line <- vapply(pl$layers, function(l) inherits(l$geom, "GeomLine"),
                      logical(1))
    b$data[[which(is_line)[1]]]$colour
  }))
  expect_true("grey55" %in% cols)
})

test_that("plot_realloc one-to-all matches predict_realloc estimates", {
  cd <- declare_composition(timeuse, parts, closure = 1440)
  fit <- comp_iv(cd, "obesity", model = "glm", family = "binomial")
  labs <- tc_part_labels(parts, cd$labels)
  # n_points = 3 -> the sweep visits delta = 30 and 60 in each direction.
  cells <- tc_realloc_cells(fit, parts, labs, delta_max = 60, n_points = 3,
                            type = "relative", scheme = "one-to-all")
  sb <- cells[["SB"]]$data
  gain <- sb[sb$signed_delta == 30, ]
  ref <- predict_realloc(fit, delta = 30, scheme = "one-to-all", to = "SB")
  expect_equal(gain$estimate, ref$estimate)
  expect_equal(gain$conf.low, ref$conf.low)
  expect_equal(gain$receiving, "SB")
  # The give direction pools the remaining parts under "(rest)".
  give <- sb[sb$signed_delta == -30, ]
  expect_equal(give$receiving, "(rest)")
  expect_equal(give$estimate,
               predict_realloc(fit, delta = 30, scheme = "one-to-all",
                               from = "SB")$estimate)
})

test_that("plot_realloc exclude drops that part's panels", {
  cd <- declare_composition(timeuse, parts, closure = 1440,
                            labels = c(MVPA = "Exercise"))
  fit <- comp_iv(cd, "obesity", model = "glm", family = "binomial")
  p <- plot_realloc(fit, delta_max = 20, n_points = 9, exclude = "MVPA")$figure
  expect_equal(length(p), choose(length(parts) - 1, 2))
  titles <- c(
    vapply(p$patches$plots, function(pl) pl$labels$title, character(1)),
    p$labels$title
  )
  expect_false(any(grepl("Exercise", titles)))
  # The excluded part also leaves the shared legend.
  sc <- p$patches$plots[[1]]$scales$get_scales("colour")
  expect_false("MVPA" %in% sc$limits)

  # One-to-all drops exactly one panel.
  pa <- plot_realloc(fit, delta_max = 20, n_points = 9,
                     scheme = "one-to-all", exclude = "MVPA")
  expect_equal(length(pa$figure), length(parts) - 1)
})

test_that("plot_realloc exclude is a display filter, not a refit", {
  cd <- declare_composition(timeuse, parts, closure = 1440)
  fit <- comp_iv(cd, "obesity", model = "glm", family = "binomial")
  labs <- tc_part_labels(parts, cd$labels)
  # Excluding MVPA prunes panels; SB still spreads its reallocation over every
  # other part, so the estimates must be identical to the unfiltered ones.
  kept <- setdiff(parts, "MVPA")
  cells <- tc_realloc_cells(fit, kept, tc_part_labels(kept, cd$labels),
                            delta_max = 60, n_points = 3, type = "relative",
                            scheme = "one-to-all")
  full <- tc_realloc_cells(fit, parts, labs, delta_max = 60, n_points = 3,
                           type = "relative", scheme = "one-to-all")
  sb_cut <- cells[["SB"]]$data
  sb_all <- full[["SB"]]$data
  expect_identical(sb_cut, sb_all)
})

test_that("plot_realloc rejects an unknown or over-eager exclude", {
  cd <- declare_composition(timeuse, parts, closure = 1440)
  fit <- comp_iv(cd, "obesity", model = "glm", family = "binomial")
  expect_error(plot_realloc(fit, exclude = "walking"), "not in the fitted")
  expect_error(plot_realloc(fit, exclude = 1), "character vector")
  expect_error(
    plot_realloc(fit, exclude = c("MVPA", "LPA", "SB")),
    "too few to plot"
  )
})

test_that("plot_realloc labels both x directions as positive minutes", {
  cd <- declare_composition(timeuse, parts, closure = 1440)
  fit <- comp_iv(cd, "obesity", model = "glm", family = "binomial")
  p <- plot_realloc(fit, delta_max = 60, n_points = 9)
  sc <- p$panels[["MVPA.LPA"]]$scales$get_scales("x")
  # The axis still runs negative -> positive; only the labels are magnitudes,
  # so a reader reads direction off the colour instead of the sign.
  expect_equal(sc$labels(c(-30, 0, 30)), c(30, 0, 30))
})

test_that("plot_realloc panel titles are bidirectional and centred", {
  cd <- declare_composition(timeuse, parts, closure = 1440)
  fit <- comp_iv(cd, "obesity", model = "glm", family = "binomial")
  p <- plot_realloc(fit, delta_max = 20, n_points = 9)$figure
  titles <- c(
    vapply(p$patches$plots, function(pl) pl$labels$title, character(1)),
    p$labels$title
  )
  both_ways <- intToUtf8(8596L)
  expect_true(all(grepl(both_ways, titles, fixed = TRUE)))
  expect_false(any(grepl(intToUtf8(8594L), titles, fixed = TRUE)))
  expect_equal(p$patches$plots[[1]]$theme$plot.title$hjust, 0.5)
})

test_that("realloc cells carry integer grid coordinates in both schemes", {
  # ceiling(sqrt(D)) is a double, so the one-to-all wrapped grid once produced
  # double row/col indices and aborted the layout's vapply(). Both builders
  # must agree on the coordinate type.
  cd <- declare_composition(timeuse, parts, closure = 1440)
  fit <- comp_iv(cd, "obesity", model = "glm", family = "binomial")
  labs <- tc_part_labels(parts, cd$labels)
  for (scheme in c("one-to-one", "one-to-all")) {
    cells <- tc_realloc_cells(fit, parts, labs, delta_max = 60, n_points = 3,
                              type = "relative", scheme = scheme)
    expect_true(all(vapply(cells, function(x) is.integer(x$row), logical(1))))
    expect_true(all(vapply(cells, function(x) is.integer(x$col), logical(1))))
  }
})

test_that("the reference line follows the effect scale, not the model name", {
  cd <- declare_composition(timeuse, parts, closure = 1440)
  # The zero-minute anchor rows sit exactly on the reference value.
  ref_of <- function(fit) {
    p <- plot_realloc(fit, delta_max = 60, n_points = 9)
    d <- p$panels[["MVPA.LPA"]]$data
    unique(d$estimate[d$signed_delta == 0])
  }
  expect_equal(ref_of(comp_iv(cd, "zbmi", model = "lm")), 0)
  expect_equal(ref_of(comp_iv(cd, "zbmi", model = "glm", family = "gaussian")), 0)
  expect_equal(ref_of(comp_iv(cd, "obesity", model = "glm", family = "binomial")), 1)
})

# --- The sweep follows the reference predict_realloc() uses -------------------

test_that("a user reference sets the sweep's feasible range", {
  cd <- declare_composition(timeuse, parts, closure = 1440)
  fit <- comp_iv(cd, "zbmi", model = "lm")
  # 15 min of MVPA can give 10 but not 20 or 30.
  ref <- c(MVPA = 15, LPA = 255, SB = 660, SPT = 510)
  labs <- tc_part_labels(parts, cd$labels)
  cells <- tc_realloc_cells(fit, parts, labs, delta_max = 30, n_points = 5,
                            type = "relative", reference = ref,
                            .reference = ref)
  mvpa_lpa <- cells[["MVPA.LPA"]]$data
  gives <- mvpa_lpa[mvpa_lpa$receiving == "LPA", ]
  expect_equal(gives$delta, 10)
  gains <- mvpa_lpa[mvpa_lpa$receiving == "MVPA", ]
  expect_equal(sort(gains$delta), c(10, 20, 30))
  # And the drawn estimates are predict_realloc()'s at that reference.
  expect_equal(
    gains$estimate[order(gains$delta)],
    predict_realloc(fit, delta = c(10, 20, 30), from = "LPA", to = "MVPA",
                    reference = ref)$estimate
  )
  # plot_realloc() passes the reference through without aborting.
  p <- plot_realloc(fit, delta_max = 30, n_points = 5, reference = ref)
  expect_s3_class(p, "realloc_plots")
  expect_s3_class(p$figure, "patchwork")
  pa <- plot_realloc(fit, delta_max = 30, n_points = 5, reference = ref,
                     scheme = "one-to-all")
  expect_equal(length(pa$figure), length(parts))
})

test_that("a user reference with more room extends the sweep", {
  cd <- declare_composition(timeuse, parts, closure = 1440)
  fit <- comp_iv(cd, "zbmi", model = "lm")
  # The fit's own reference has under 30 min of MVPA; this one has 60.
  expect_lt(fit$ref_comp$MVPA, 30)
  ref <- c(MVPA = 60, LPA = 240, SB = 630, SPT = 510)
  labs <- tc_part_labels(parts, cd$labels)
  cells <- tc_realloc_cells(fit, parts, labs, delta_max = 30, n_points = 5,
                            type = "relative", scheme = "one-to-all",
                            reference = ref, .reference = ref)
  gives <- cells[["MVPA"]]$data
  expect_equal(sort(gives$delta[gives$receiving == tc_rest_key()]),
               c(10, 20, 30))
})

test_that("delta_max must be a single positive number", {
  cd <- declare_composition(timeuse, parts, closure = 1440)
  fit <- comp_iv(cd, "zbmi", model = "lm")
  expect_error(plot_realloc(fit, delta_max = 0), "single positive number")
  expect_error(plot_realloc(fit, delta_max = -5), "single positive number")
  expect_error(plot_realloc(fit, delta_max = c(10, 20)),
               "single positive number")
})

test_that("n_points must be a single positive whole number", {
  cd <- declare_composition(timeuse, parts, closure = 1440)
  fit <- comp_iv(cd, "zbmi", model = "lm")
  for (bad in list("a", 0, -3, 2.5, c(5, 9), NA)) {
    expect_error(plot_realloc(fit, delta_max = 20, n_points = bad),
                 "`n_points` must be a single positive whole number.",
                 fixed = TRUE)
  }
})

test_that("n_points splits across the two directions", {
  cd <- declare_composition(timeuse, parts, closure = 1440)
  fit <- comp_iv(cd, "zbmi", model = "lm")
  labs <- tc_part_labels(parts, cd$labels)
  # ceiling(n_points / 2) points per direction; n_points = 1 sweeps delta_max.
  one <- tc_realloc_cells(fit, parts, labs, delta_max = 20, n_points = 1,
                          type = "relative")
  expect_equal(sort(unique(abs(one[["MVPA.LPA"]]$data$signed_delta))), 20)
  nine <- tc_realloc_cells(fit, parts, labs, delta_max = 20, n_points = 9,
                           type = "relative")
  swept <- nine[["MVPA.LPA"]]$data$signed_delta
  expect_equal(sort(unique(swept[swept > 0])), c(4, 8, 12, 16, 20))
})

test_that("the one-to-all pooled side uses the prediction's rest key", {
  cd <- declare_composition(timeuse, parts, closure = 1440)
  fit <- comp_iv(cd, "zbmi", model = "lm")
  r <- predict_realloc(fit, delta = 10, scheme = "one-to-all", from = "SB")
  expect_equal(r$to, tc_rest_key())
  labs <- tc_part_labels(parts, cd$labels)
  cells <- tc_realloc_cells(fit, parts, labs, delta_max = 20, n_points = 3,
                            type = "relative", scheme = "one-to-all")
  expect_true(tc_rest_key() %in% cells[["MVPA"]]$data$receiving)
  p <- plot_realloc(fit, delta_max = 20, n_points = 3, scheme = "one-to-all")
  sc <- p$panels[["MVPA"]]$scales$get_scales("colour")
  expect_equal(unname(sc$map(tc_rest_key())), "grey55")
})

# --- Grouped results: one grid per group -------------------------------------

# The plots of a patchwork, in layout order.
patch_plots <- function(p) lapply(seq_len(length(p)), function(i) p[[i]])

test_that("a by grouping stacks one titled grid per level on one y-range", {
  cd <- declare_composition(timeuse, parts, closure = 1440)
  fit <- comp_iv(cd, "zbmi", model = "lm", covariates = c("sex", "ilr1:sex"))
  levs <- levels(factor(timeuse$sex))
  p <- plot_realloc(fit, delta_max = 20, n_points = 5, by = "sex")
  expect_s3_class(p, "realloc_plots")
  expect_s3_class(p$figure, "patchwork")
  grids <- patch_plots(p$figure)
  expect_length(grids, length(levs))
  expect_equal(p$figure$patches$layout$ncol, 1)
  titles <- vapply(grids, function(g) g[[1]]$labels$title, character(1))
  expect_identical(titles, paste0("sex: ", levs))
  # Each grid holds the full panel triangle, drawing only its own level.
  for (i in seq_along(grids)) {
    panels <- patch_plots(grids[[i]])[-1]
    expect_length(panels, choose(length(parts), 2))
    for (pl in panels) {
      expect_true(all(stats::na.omit(pl$data$sex) == levs[i]))
    }
  }
  # One y-range across every panel of every grid.
  y_ranges <- unlist(lapply(grids, function(g) {
    lapply(patch_plots(g)[-1], function(pl) {
      list(ggplot2::ggplot_build(pl)$layout$panel_params[[1]]$y.range)
    })
  }), recursive = FALSE)
  expect_true(all(vapply(y_ranges, identical, logical(1), y_ranges[[1]])))
  # The drawn estimates are predict_realloc()'s, level by level (the first
  # panel's positive minutes move time from MVPA to LPA).
  r <- predict_realloc(fit, delta = 20, from = "MVPA", to = "LPA", by = "sex")
  first <- patch_plots(grids[[2]])[[2]]$data
  expect_equal(first$estimate[first$signed_delta == 20],
               r$estimate[r$sex == levs[2]])
})

test_that("without a grouping the figure is the single panel grid", {
  cd <- declare_composition(timeuse, parts, closure = 1440)
  fit <- comp_iv(cd, "zbmi", model = "lm", covariates = "sex")
  p <- plot_realloc(fit, delta_max = 20, n_points = 5)$figure
  expect_equal(length(p), choose(length(parts), 2))
  expect_false(any(vapply(patch_plots(p), inherits, logical(1), "patchwork")))
})

test_that("a multinomial fit draws one grid per outcome category", {
  skip_if_not_installed("mclogit")
  cd <- declare_composition(timeuse, parts, closure = 1440)
  cd$data_raw$grade_u <- factor(timeuse$grade, ordered = FALSE)
  fit <- suppressWarnings(
    comp_iv(cd, outcome = "grade_u", model = "mblogit", covariates = "sex")
  )
  cats <- unique(predict_realloc(fit, delta = 10, from = "SB",
                                 to = "MVPA")$response_level)
  p <- plot_realloc(fit, delta_max = 20, n_points = 3)
  titles <- vapply(patch_plots(p$figure), function(g) g[[1]]$labels$title,
                   character(1))
  expect_identical(titles, paste0("outcome: ", cats))
  levs <- levels(factor(timeuse$sex))
  pb <- plot_realloc(fit, delta_max = 20, n_points = 3, by = "sex")
  titles_b <- vapply(patch_plots(pb$figure), function(g) g[[1]]$labels$title,
                     character(1))
  expect_identical(titles_b, paste0("sex: ", rep(levs, each = length(cats)),
                                    ", outcome: ",
                                    rep(cats, times = length(levs))))
})

test_that("plot_realloc's axis title carries the declared unit label, with no caption", {
  x_titles <- function(p) {
    p <- p$figure
    unlist(c(purrr::map(p$patches$plots, function(pl) pl$labels$x),
             list(p$labels$x)))
  }
  cdu <- declare_composition(timeuse, parts, closure = 1440, units = "min/day")
  fit_u <- comp_iv(cdu, "obesity", model = "glm", family = "binomial")
  p_u <- plot_realloc(fit_u, delta_max = 20, n_points = 5)
  expect_true(all(x_titles(p_u) == "Time reallocated (min/day)"))
  expect_false(is.character(p_u$figure$patches$annotation$caption))
  expect_true(all(purrr::map_lgl(p_u$panels, function(panel) {
    is.null(panel$labels$caption)
  })))
  cd <- declare_composition(timeuse, parts, closure = 1440)
  fit <- comp_iv(cd, "obesity", model = "glm", family = "binomial")
  p <- plot_realloc(fit, delta_max = 20, n_points = 5)
  expect_true(all(x_titles(p) == "Time reallocated"))
  expect_length(x_titles(p), length(parts) - 1)
})

# --- The returned panels -----------------------------------------------------

pair_names <- function(parts) {
  purrr::map_chr(utils::combn(parts, 2, simplify = FALSE), paste,
                 collapse = ".")
}

test_that("one-to-one panels are named by their pair of parts", {
  cd <- declare_composition(timeuse, parts, closure = 1440)
  fit <- comp_iv(cd, "obesity", model = "glm", family = "binomial")
  p <- plot_realloc(fit, delta_max = 20, n_points = 5)
  expect_identical(names(p$panels), pair_names(parts))
  expect_length(p$panels, 6)
  expect_true(all(purrr::map_lgl(p$panels, inherits, "ggplot")))
  expect_false(any(purrr::map_lgl(p$panels, inherits, "patchwork")))
})

test_that("one-to-all panels are named by their part", {
  cd <- declare_composition(timeuse, parts, closure = 1440)
  fit <- comp_iv(cd, "obesity", model = "glm", family = "binomial")
  p <- plot_realloc(fit, delta_max = 20, n_points = 5, scheme = "one-to-all")
  expect_identical(names(p$panels), parts)
})

test_that("an excluded part leaves no panel name", {
  cd <- declare_composition(timeuse, parts, closure = 1440)
  fit <- comp_iv(cd, "obesity", model = "glm", family = "binomial")
  p <- plot_realloc(fit, delta_max = 20, n_points = 5, exclude = "MVPA")
  expect_identical(names(p$panels), c("LPA.SB", "LPA.SPT", "SB.SPT"))
})

test_that("grouped panels carry the group level as a name prefix", {
  cd <- declare_composition(timeuse, parts, closure = 1440)
  fit <- comp_iv(cd, "zbmi", model = "lm", covariates = c("sex", "ilr1:sex"))
  levs <- levels(factor(timeuse$sex))
  p <- plot_realloc(fit, delta_max = 20, n_points = 5, by = "sex")
  expect_length(p$panels, length(levs) * 6)
  expect_identical(
    names(p$panels),
    paste(rep(levs, each = 6), rep(pair_names(parts), times = length(levs)),
          sep = ".")
  )
  # The group headers are not panels: every panel draws its own data.
  expect_true(all(purrr::map_lgl(p$panels, function(pl) {
    is.data.frame(pl$data) && nrow(pl$data) > 0
  })))
})

test_that("a multinomial fit prefixes panel names with the outcome category", {
  skip_if_not_installed("mclogit")
  cd <- declare_composition(timeuse, parts, closure = 1440)
  cd$data_raw$grade_u <- factor(timeuse$grade, ordered = FALSE)
  fit <- suppressWarnings(
    comp_iv(cd, outcome = "grade_u", model = "mblogit", covariates = "sex")
  )
  cats <- unique(predict_realloc(fit, delta = 10, from = "SB",
                                 to = "MVPA")$response_level)
  p <- plot_realloc(fit, delta_max = 20, n_points = 3)
  expect_identical(
    names(p$panels),
    paste(rep(cats, each = 6), rep(pair_names(parts), times = length(cats)),
          sep = ".")
  )
})

test_that("printing a realloc_plots object draws it and returns it invisibly", {
  cd <- declare_composition(timeuse, parts, closure = 1440)
  fit <- comp_iv(cd, "obesity", model = "glm", family = "binomial")
  p <- plot_realloc(fit, delta_max = 20, n_points = 3)
  out <- expect_invisible(print(p))
  expect_identical(out, p)
})

test_that("legend = FALSE hides the legend on every returned panel", {
  cd <- declare_composition(timeuse, parts, closure = 1440)
  fit <- comp_iv(cd, "obesity", model = "glm", family = "binomial")
  p <- plot_realloc(fit, delta_max = 20, n_points = 3, legend = FALSE)
  positions <- purrr::map_chr(p$panels, function(pl) {
    as.character(pl$theme$legend.position)[1]
  })
  expect_true(all(positions == "none"))
})

test_that("realloc cells are named as the panels are", {
  cd <- declare_composition(timeuse, parts, closure = 1440)
  fit <- comp_iv(cd, "obesity", model = "glm", family = "binomial")
  labs <- tc_part_labels(parts, cd$labels)
  cells <- tc_realloc_cells(fit, parts, labs, delta_max = 20, n_points = 3,
                            type = "relative")
  p <- plot_realloc(fit, delta_max = 20, n_points = 3)
  expect_identical(names(cells), names(p$panels))
})

test_that("plot_realloc() refuses a conf_level outside (0, 1) through predict_realloc()", {
  cd <- declare_composition(timeuse, parts, closure = 1440)
  fit <- comp_iv(cd, "zbmi", model = "lm")
  expect_error(plot_realloc(fit, delta_max = 20, n_points = 3, conf_level = 95),
               "conf_level.*single number in \\(0, 1\\)")
})
