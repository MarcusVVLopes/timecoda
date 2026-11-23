parts <- c("MVPA", "LPA", "SB", "SPT")
cd <- declare_composition(timeuse, parts, closure = 1440)

# Name of the fill scale on a ggtern/ggplot object (the surface's legend title).
fill_scale_name <- function(p) {
  for (s in p$scales$scales) {
    if ("fill" %in% s$aesthetics) {
      return(s$name)
    }
  }
  NA_character_
}

test_that("plot.comp_iv returns one tiled ternary per part triplet", {
  m <- comp_iv(cd, "obesity", model = "glm", family = "binomial")
  p <- plot(m, resolution = 10, smooth = FALSE)
  expect_type(p, "list")
  expect_length(p, choose(length(parts), 3))
  expect_s3_class(p[[1]], "ggplot")
})

test_that("relative surface is named for the effect measure", {
  m <- comp_iv(cd, "obesity", model = "glm", family = "binomial")
  p <- plot(m, resolution = 10, smooth = FALSE)
  expect_equal(fill_scale_name(p[[1]]), "Odds ratio")
})

# The relative surface over one triplet, from the panel's data (not drawn).
panel_surface <- function(m, tri = m$parts[1:3]) {
  tc_surface_panel(m, tri, unlist(m$ref_comp)[m$parts], 6, "relative",
                   "value")$data
}

# The drawn fine tiles equal a quadratic refitted here on the link scale to
# the panel's computed cells, evaluated at each tile's centre (the mean of its
# six corners) and mapped back.
expect_drawn_on_link_scale <- function(panel, link, inv) {
  d <- panel$data
  uv <- function(s) {
    ilr_transform(as.data.frame(s[c("x", "y", "z")]), c("x", "y", "z"))
  }
  cells <- uv(d)
  fd <- data.frame(u = cells$ilr1, v = cells$ilr2, value = link(d$.pred))
  fit <- stats::lm(value ~ stats::poly(u, v, degree = 2), data = fd)
  tiles <- panel$layers[[1]]$data
  centres <- stats::aggregate(cbind(x, y, z, .pred) ~ cell, data = tiles,
                              FUN = mean)
  at <- uv(centres)
  refit <- inv(stats::predict(fit, data.frame(u = at$ilr1, v = at$ilr2)))
  expect_gt(nrow(centres), nrow(d))
  expect_equal(centres$.pred, unname(refit), tolerance = 1e-6)
}

# Linear-predictor difference between the first grid point and the reference
# composition, by hand from the coordinate coefficients of a covariate-free fit.
hand_contrast <- function(m, d) {
  ref <- unlist(m$ref_comp)[m$parts]
  comp <- ref
  comp[m$parts[1:3]] <- c(d$x[1], d$y[1], d$z[1])
  coords <- tc_coord_names(length(m$parts))
  z <- function(v) unlist(ilr_transform(as.data.frame(t(v)), m$parts)[coords])
  sum((z(comp) - z(ref)) * stats::coef(m$model)[coords])
}

test_that("relative surface is the odds ratio versus the reference (ratio scale)", {
  m <- comp_iv(cd, "obesity", model = "glm", family = "binomial")
  d <- panel_surface(m)
  expect_equal(d$.pred[1], exp(hand_contrast(m, d)), tolerance = 1e-6)
})

test_that("relative surface is the mean difference versus the reference", {
  m <- comp_iv(cd, "zbmi", model = "lm")
  d <- panel_surface(m)
  expect_equal(d$.pred[1], hand_contrast(m, d), tolerance = 1e-6)
})

test_that("absolute surface renders for a glm fit", {
  m <- comp_iv(cd, "obesity", model = "glm", family = "binomial")
  p <- plot(m, scale = "absolute", resolution = 10, smooth = FALSE)
  expect_type(p, "list")
  expect_length(p, choose(length(parts), 3))
  expect_s3_class(p[[1]], "ggplot")
  expect_equal(fill_scale_name(p[[1]]), "Predicted response")
})

test_that("absolute surface renders for an lm fit", {
  m <- comp_iv(cd, "zbmi", model = "lm")
  p <- plot(m, scale = "absolute", resolution = 10, smooth = FALSE)
  expect_s3_class(p[[1]], "ggplot")
})

test_that("cox rejects the absolute surface cleanly", {
  m <- comp_iv(cd, time = "time", event = "event", model = "cox")
  expect_error(plot(m, scale = "absolute"), "relative")
  # The message is compared with its line breaks read as spaces.
  err <- expect_error(plot(m, scale = "absolute"))
  expect_match(gsub("\\s+", " ", conditionMessage(err)),
               "(a Cox model leaves the baseline hazard unestimated)",
               fixed = TRUE)
  # But the relative (hazard-ratio) surface still renders.
  p <- plot(m, resolution = 10, smooth = FALSE)
  expect_s3_class(p[[1]], "ggplot")
})

test_that("relative surface is invariant to the ILR basis (part order)", {
  m1 <- comp_iv(cd, "obesity", model = "glm", family = "binomial")
  cd2 <- declare_composition(timeuse, rev(parts), closure = 1440)
  m2 <- comp_iv(cd2, "obesity", model = "glm", family = "binomial")
  # Evaluate both fits (which use different ILR bases) on the SAME grid of
  # physical compositions: one triplet, named in the same order for both.
  expect_equal(panel_surface(m1, parts[1:3])$.pred,
               panel_surface(m2, parts[1:3])$.pred, tolerance = 1e-4)
})

test_that("plot.comp_iv panel count tracks the number of parts", {
  for (d in c(3, 5)) {
    pn <- paste0("p", seq_len(d))
    m <- comp_iv(declare_composition(sim_parts(d), pn), "y", model = "lm")
    pl <- plot(m, resolution = 8, quantiles = NULL)
    if (d == 3) {
      expect_s3_class(pl, "ggplot")
    } else {
      expect_length(pl, choose(d, 3))
    }
  }
})

test_that("plot.comp_iv rejects non-comp_iv input", {
  expect_error(plot.comp_iv(list()), "comp_iv")
})

test_that("plot.comp_iv draws without warnings", {
  m <- comp_iv(cd, "obesity", model = "glm", family = "binomial")
  expect_no_warning(quiet_draw(plot(m)))
})

test_that("no panel carries a caption", {
  no_caption <- function(p) {
    expect_true(all(purrr::map_lgl(p, function(panel) is.null(panel$labels$caption))))
  }
  for (d in c(4, 5)) {
    pn <- paste0("p", seq_len(d))
    m <- comp_iv(declare_composition(sim_parts(d), pn, closure = 1440,
                                     units = "min/day"),
                 "y", model = "lm")
    no_caption(plot(m, resolution = 8, quantiles = NULL))
  }
  m <- comp_iv(cd, "obesity", model = "glm", family = "binomial")
  no_caption(plot(m, resolution = 10, smooth = FALSE))
  no_caption(plot(m, resolution = 8, quantiles = NULL, smooth = TRUE))
  tu <- as.data.frame(timeuse)
  set.seed(1)
  tu[parts] <- tu[parts] * stats::runif(nrow(tu), 0.8, 1.2)
  m_total <- comp_iv(declare_composition(tu, parts, closure = "mean"), "zbmi",
                     model = "lm", include_total = TRUE)
  no_caption(plot(m_total, resolution = 8, quantiles = NULL))
})

test_that("an argument plot() does not take is an error that names it", {
  m <- comp_iv(cd, "zbmi", model = "lm")
  flat <- function(err) gsub("\\s+", " ", conditionMessage(err))
  msg <- flat(expect_error(plot(m, coverage = 0.9), class = "rlang_error"))
  expect_match(
    msg, "`plot.comp_iv()` was given an argument it does not take: `coverage`.",
    fixed = TRUE
  )
  msg <- flat(expect_error(plot(m, sd = 2, coverage = 0.9),
                           class = "rlang_error"))
  expect_match(msg, "was given arguments it does not take: `sd` and `coverage`.",
               fixed = TRUE)
  # An unnamed one, after every named argument, by its position in `...`.
  msg <- flat(expect_error(
    plot(m, "relative", 10, NULL, FALSE, NULL, 0.5, sd = 2),
    class = "rlang_error"
  ))
  expect_match(
    msg,
    "an unnamed argument in position 1 of `...` and `sd`.",
    fixed = TRUE
  )
})

test_that("plot.comp_iv refuses a bad resolution or quantile band", {
  m <- comp_iv(cd, "zbmi", model = "lm")
  expect_error(plot(m, resolution = 13.5), "single whole number")
  expect_error(plot(m, resolution = 2), "at least 3")
  expect_error(plot(m, quantiles = c(0.95, 0.05)), "increasing order")
})

test_that("a panel with no grid composition inside the band is drawn blank", {
  m <- comp_iv(cd, "zbmi", model = "lm")
  w <- character()
  p <- withCallingHandlers(
    plot(m, resolution = 3, smooth = FALSE),
    warning = function(x) {
      w <<- c(w, conditionMessage(x))
      invokeRestart("muffleWarning")
    }
  )
  expect_length(p, choose(length(parts), 3))
  expect_length(w, choose(length(parts), 3))
  expect_match(w[1], "the panel is blank")
  expect_equal(nrow(p[[1]]$data), 0)
})

test_that("a blank panel of a search's region names the region, not `quantiles`", {
  m <- comp_iv(cd, "zbmi", model = "lm")
  reference <- tc_check_reference(NULL, m)[m$parts]
  w <- character()
  p <- withCallingHandlers(
    tc_surface_panel(m, parts[1:3], reference, 8, "relative", "Difference",
                     bounds = function(comps) rep(FALSE, nrow(comps)),
                     link = identity, inv = identity),
    warning = function(x) {
      w <<- c(w, conditionMessage(x))
      invokeRestart("muffleWarning")
    }
  )
  expect_length(w, 1)
  expect_match(w, "region; the panel is blank")
  expect_match(w, "Raise")
  expect_false(grepl("quantiles", w))
  expect_equal(nrow(p$data), 0)
})

test_that("smooth = TRUE draws a heatmap over the same computed cells", {
  m <- comp_iv(cd, "obesity", model = "glm", family = "binomial")
  tiled <- plot(m, resolution = 8, quantiles = NULL, smooth = FALSE)
  expect_no_warning(
    p <- plot(m, resolution = 8, quantiles = NULL, smooth = TRUE)
  )
  expect_length(p, choose(length(parts), 3))
  expect_s3_class(p[[1]], "ggplot")
  # Smoothed: more drawn tiles than computed cells.
  expect_gt(nrow(p[[1]]$layers[[1]]$data), 6 * nrow(p[[1]]$data))
  expect_identical(p[[1]]$data, tiled[[1]]$data)
  expect_error(plot(m, smooth = "yes"), "TRUE")
})

test_that("the default draws the computed values as tiles", {
  m <- comp_iv(cd, "obesity", model = "glm", family = "binomial")
  tiled <- plot(m, resolution = 8, quantiles = NULL, smooth = FALSE)
  default <- plot(m, resolution = 8, quantiles = NULL)
  expect_identical(default[[1]]$data, tiled[[1]]$data)
  expect_identical(default[[1]]$layers[[1]]$data, tiled[[1]]$layers[[1]]$data)
  expect_equal(nrow(default[[1]]$layers[[1]]$data), 6 * nrow(default[[1]]$data))
})

test_that("a moderated fit names once the covariate the surface averages", {
  m <- comp_iv(cd, "obesity", model = "glm", family = "binomial",
               moderator = "sex")
  old <- getOption("emmeans")
  msgs <- list()
  withCallingHandlers(plot(m, resolution = 4, quantiles = NULL, smooth = FALSE),
                      message = function(x) {
                        msgs[[length(msgs) + 1]] <<- x
                        invokeRestart("muffleMessage")
                      })
  expect_length(msgs, 1)
  expect_s3_class(msgs[[1]], "rlang_message")
  expect_match(conditionMessage(msgs[[1]]),
               "marginal over the interacting covariate", fixed = TRUE)
  expect_identical(getOption("emmeans"), old)
  expect_error(plot(m, resolution = 2), "at least 3")
  expect_identical(getOption("emmeans"), old)
})

test_that("a constant surface keeps its tiles when smoothed", {
  m <- comp_iv(cd, "zbmi", model = "lm")
  local_mocked_bindings(tc_optimal_predict = function(fit, cand, ...) {
    list(response = rep(1, nrow(cand$z)))
  })
  # 21 cells over the whole triangle at 8: only the constant guard can fire.
  expect_warning(
    p <- tc_surface_panel(m, parts[1:3], unlist(m$ref_comp)[parts], 8,
                          "relative", "value", smooth = TRUE,
                          link = identity, inv = identity),
    "cannot be fitted"
  )
  expect_equal(nrow(p$data), 21)
  expect_equal(nrow(p$layers[[1]]$data), 6 * nrow(p$data))
})

test_that("each panel marks the reference composition with one point", {
  m <- comp_iv(cd, "zbmi", model = "lm")
  p <- plot(m, resolution = 10, smooth = FALSE)[[1]]
  is_point <- purrr::map_lgl(p$layers, function(l) inherits(l$geom, "GeomPoint"))
  expect_equal(sum(is_point), 1)
  pt <- p$layers[[which(is_point)]]$data
  ref <- unlist(m$ref_comp)[parts[1:3]]
  expect_equal(unname(unlist(pt[c("x", "y", "z")])), unname(ref),
               tolerance = 1e-10)
})

test_that("the smooth stays close to the computed values", {
  m <- comp_iv(cd, "obesity", model = "glm", family = "binomial")
  # At 25 the default band keeps at least 12 cells in every panel.
  p <- plot(m, resolution = 25, smooth = TRUE)
  for (panel in p) {
    d <- panel$data
    uv <- ilr_transform(as.data.frame(d[c("x", "y", "z")]), c("x", "y", "z"))
    fit <- stats::lm(log(d$.pred) ~ stats::poly(uv$ilr1, uv$ilr2, degree = 2))
    expect_lte(max(abs(exp(stats::fitted(fit)) - d$.pred)),
               0.05 * diff(range(d$.pred)))
    fill <- panel$scales$get_scales("fill")
    expect_equal(fill$limits, range(d$.pred))
    beyond <- range(d$.pred) + c(-1, 1)
    expect_equal(fill$oob(beyond, fill$limits), range(d$.pred))
  }
  expect_drawn_on_link_scale(p[[1]], log, exp)
})

test_that("the absolute smooth is fitted on the link scale", {
  m <- comp_iv(cd, "obesity", model = "glm", family = "binomial")
  p <- plot(m, scale = "absolute", resolution = 25, smooth = TRUE)
  for (panel in p) {
    d <- panel$data
    uv <- ilr_transform(as.data.frame(d[c("x", "y", "z")]), c("x", "y", "z"))
    fit <- stats::lm(stats::qlogis(d$.pred) ~
                       stats::poly(uv$ilr1, uv$ilr2, degree = 2))
    expect_lte(max(abs(stats::plogis(stats::fitted(fit)) - d$.pred)),
               0.05 * diff(range(d$.pred)))
    tiles <- panel$layers[[1]]$data$.pred
    expect_true(all(tiles >= 0 & tiles <= 1))
  }
  expect_drawn_on_link_scale(p[[1]], stats::qlogis, stats::plogis)
})

# A multinomial fit of timeuse's grade, as an unordered factor.
mblogit_fit <- function() {
  d <- timeuse
  d$grade_unordered <- factor(as.character(d$grade), levels = levels(d$grade))
  suppressWarnings(comp_iv(declare_composition(d, parts, closure = 1440),
                           "grade_unordered", model = "mblogit",
                           covariates = "sex"))
}

test_that("the ordinal and multinomial engines draw a relative surface", {
  skip_if_not_installed("ordinal")
  m <- comp_iv(cd, "grade", model = "clm")
  p <- plot(m, resolution = 10)
  expect_length(p, choose(length(parts), 3))
  expect_equal(fill_scale_name(p[[1]]), "Cumulative odds ratio")
  expect_error(plot(m, scale = "absolute"),
               "support `scale = \"relative\"` only")
  expect_error(plot(m, category = "above"),
               "applies to multinomial fits only")
  skip_if_not_installed("mclogit")
  mb <- mblogit_fit()
  levs <- unique(sub("~.*$", "", names(stats::coef(mb$model))))
  p <- plot(mb, resolution = 10, category = levs[2])
  expect_length(p, choose(length(parts), 3))
  expect_equal(fill_scale_name(p[[1]]),
               paste("Relative risk ratio of", levs[2]))
  expect_error(plot(mb, scale = "absolute", category = levs[1]),
               "support `scale = \"relative\"` only")
  # `category` is required, the message naming the levels, and must be one;
  # the messages are read with their line breaks as spaces.
  flat <- function(err) gsub("\\s+", " ", conditionMessage(err))
  msg <- flat(expect_error(plot(mb)))
  expect_match(msg, "must name the outcome level to draw", fixed = TRUE)
  expect_match(msg, levs[1], fixed = TRUE)
  expect_match(msg, levs[2], fixed = TRUE)
  expect_match(flat(expect_error(plot(mb, category = "none"))),
               "non-reference outcome levels", fixed = TRUE)
})

test_that("category is refused for a fit that is not multinomial", {
  m <- comp_iv(cd, "zbmi", model = "lm")
  expect_error(plot(m, category = "yes"), "applies to multinomial fits only")
})

# The panel's values at sampled cells, against one emmeans grid per
# composition and the contrast with the reference row.
expect_per_point_surface <- function(m, scale = "relative", category = NULL) {
  ref <- unlist(m$ref_comp)[m$parts]
  tri <- m$parts[1:3]
  d <- tc_surface_panel(m, tri, ref, 10, scale, "value",
                        category = category)$data
  set.seed(3)
  i <- sort(sample(nrow(d), 8))
  comps <- matrix(ref, nrow = length(i) + 1, ncol = length(ref), byrow = TRUE,
                  dimnames = list(NULL, m$parts))
  comps[-1, tri] <- as.matrix(d[i, c("x", "y", "z")])
  emm <- tc_emm_at_points(m, comps, at = NULL, conf_level = 0.95,
                          category = category)
  hand <- if (scale == "relative") {
    tc_estimate_col(summary(emmeans::contrast(emm, method = "trt.vs.ctrl",
                                              ref = 1),
                            infer = c(FALSE, FALSE), type = "response",
                            adjust = "none"))
  } else {
    tc_estimate_col(summary(emm, infer = c(FALSE, FALSE),
                            type = "response"))[-1]
  }
  expect_equal(d$.pred[i], hand, tolerance = 1e-10)
}

test_that("the surface equals the per-point predictions at sampled cells", {
  lm_fit <- comp_iv(cd, "zbmi", model = "lm", covariates = "sex")
  glm_fit <- comp_iv(cd, "obesity", model = "glm", family = "binomial")
  for (scale in c("relative", "absolute")) {
    expect_per_point_surface(lm_fit, scale)
    expect_per_point_surface(glm_fit, scale)
  }
  expect_per_point_surface(comp_iv(cd, time = "time", event = "event",
                                   model = "cox"))
  skip_if_not_installed("lme4")
  lmer_fit <- comp_iv(cd, "zbmi", model = "lmer", random = "school")
  for (scale in c("relative", "absolute")) {
    expect_per_point_surface(lmer_fit, scale)
  }
  skip_if_not_installed("ordinal")
  expect_per_point_surface(comp_iv(cd, "grade", model = "clm"))
  skip_if_not_installed("mclogit")
  mb <- mblogit_fit()
  levs <- unique(sub("~.*$", "", names(stats::coef(mb$model))))
  expect_per_point_surface(mb, category = levs[1])
})

test_that("the surface plots title the corners by the display labels", {
  labs <- c(MVPA = "Moderate-to-vigorous PA", SB = "Sedentary behaviour")
  cdl <- declare_composition(timeuse, parts, closure = 1440, labels = labs)
  m <- comp_iv(cdl, "zbmi", model = "lm")
  p <- quiet_draw(plot(m, resolution = 10))[[1]]
  expect_identical(c(p$labels$x, p$labels$y, p$labels$z),
                   c("Moderate-to-vigorous PA", "LPA", "Sedentary behaviour"))
})
