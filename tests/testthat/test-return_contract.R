# Display-object return contract.
#
# Every formatted output (gt table, ggplot/patchwork) must be RETURNED to the
# caller for post-editing, not just printed as a side effect. This file audits
# that contract across the summary/plot methods and fit_check(), independent
# of the statistical-correctness tests owned by the per-function test files.
# Fixtures are kept small/fast and reuse the shipped `timeuse` data; plotting
# side effects are routed to a null device by setup.R.

parts <- c("MVPA", "LPA", "SB", "SPT")
cd <- declare_composition(timeuse, parts, closure = 1440)

# summary.compD / pub_ready() --------------------------------------------------

test_that("pub_ready(summary(cd)) returns gt elements inheriting gt_tbl", {
  skip_if_not_installed("gt")
  s <- summary(cd)
  g <- pub_ready(s)
  expect_named(g, c("desc", "var_matrix"))
  expect_s3_class(g$desc, "gt_tbl")
  expect_s3_class(g$var_matrix, "gt_tbl")
})

test_that("summary(cd) is a visible summary_compD object, not printed as a side effect of computation", {
  # summary() itself only builds the object; printing happens through the
  # separate print.summary_compD() dispatch, exercised in test-summary_compD.R.
  s <- summary(cd)
  expect_s3_class(s, "summary_compD")
})

# summary.comp_iv ---------------------------------------------------------------

test_that("pub_ready() on a summary.comp_iv returns a gt_tbl element", {
  skip_if_not_installed("gt")
  m <- comp_iv(cd, "obesity", model = "glm", family = "binomial")
  s <- summary(m)
  expect_s3_class(pub_ready(s)$effects, "gt_tbl")
})

test_that("summary.comp_iv returns a visible summary_comp_iv object", {
  m <- comp_iv(cd, "obesity", model = "glm", family = "binomial")
  expect_s3_class(summary(m), "summary_comp_iv")
  expect_silent(summary(m))
})

# summary.realloc ---------------------------------------------------------------

test_that("pub_ready(summary.realloc()) returns a gt_tbl element", {
  skip_if_not_installed("gt")
  m <- comp_iv(cd, "obesity", model = "glm", family = "binomial")
  r <- predict_realloc(m, delta = 20, from = "SB", to = "MVPA")
  s <- summary(r)
  expect_s3_class(pub_ready(s)$effects, "gt_tbl")
})

test_that("summary.realloc is silent when assigned", {
  m <- comp_iv(cd, "obesity", model = "glm", family = "binomial")
  r <- predict_realloc(m, delta = 20, from = "SB", to = "MVPA")
  expect_silent(s <- summary(r))
})

# plot.compD --------------------------------------------------------------------

test_that("plot.compD returns a named list of gg/ggplot objects", {
  p <- quiet_draw(plot(cd))
  expect_type(p, "list")
  expect_named(p, c("histogram", "boxplot", "ternary", "variation"))
  expect_s3_class(p$histogram, "gg")
  expect_s3_class(p$boxplot, "gg")
  expect_s3_class(p$variation, "ggplot")
  # 4 parts -> choose(4, 3) = 4 ternary triplets, returned as a named list.
  expect_type(p$ternary, "list")
  expect_length(p$ternary, choose(length(parts), 3))
  expect_s3_class(p$ternary[[1]], "ggplot")
})

# plot.comp_iv --------------------------------------------------------------------

test_that("plot.comp_iv returns ggplot panel(s)", {
  m <- comp_iv(cd, "obesity", model = "glm", family = "binomial")
  p <- quiet_draw(plot(m, resolution = 10, smooth = FALSE))
  expect_type(p, "list")
  expect_length(p, choose(length(parts), 3))
  expect_s3_class(p[[1]], "ggplot")
})

# plot_ternary / plot_variation (R/plot_ternary.R), plot_realloc (R/plot_realloc.R)

test_that("plot_ternary returns a named list of ggplot objects", {
  cd <- declare_composition(timeuse, parts, closure = 1440)
  p <- quiet_draw(plot_ternary(cd, group = "sex"))
  expect_type(p, "list")
  expect_length(p, choose(length(parts), 3))
  expect_s3_class(p[[1]], "ggplot")
})

test_that("plot_variation returns a single ggplot object", {
  p <- plot_variation(cd)
  expect_s3_class(p, "ggplot")
})

test_that("plot_realloc returns its panels and the assembled patchwork", {
  m <- comp_iv(cd, "obesity", model = "glm", family = "binomial")
  p <- plot_realloc(m, delta_max = 30, n_points = 5)
  expect_s3_class(p, "realloc_plots")
  expect_s3_class(p$figure, "patchwork")
  expect_s3_class(p$figure, "gg")
})

# fit_check ---------------------------------------------------------------------

test_that("fit_check() returns its object invisibly with documented elements", {
  skip_if_not_installed("DHARMa")
  m <- comp_iv(cd, "obesity", model = "glm", family = "binomial")
  expect_false(withVisible(fit_check(m, plot = FALSE))$visible)
  fc <- fit_check(m, plot = FALSE)
  expect_s3_class(fc, "fit_check")
  expect_named(fc, c("tests", "raw", "note", "backend", "model_type",
                     "convergence"))
  expect_s3_class(fc$tests, "tbl_df")
  # The diagnostic plot is drawn, not returned, so there is no `plot` element;
  # callers re-plot via plot(fc$raw).
  expect_false("plot" %in% names(fc))
  expect_false(is.null(fc$raw))
})
