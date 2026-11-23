parts <- c("MVPA", "LPA", "SB", "SPT")

test_that("plot_variation returns a ggplot heatmap of the variation matrix", {
  cd <- declare_composition(timeuse, parts, closure = 1440)
  p <- plot_variation(cd)
  expect_s3_class(p, "ggplot")
  # one tile per ordered pair of parts (D x D, diagonal included)
  expect_equal(nrow(ggplot2::ggplot_build(p)$data[[1]]), length(parts)^2)
})

test_that("the variation matrix has a zero diagonal (its minimum)", {
  cd <- declare_composition(timeuse, parts, closure = 1440)
  vm <- variation_matrix(tc_comp_df(cd), parts)
  expect_equal(max(diag(vm)), 0)
  expect_equal(min(vm), 0)
})

test_that("plot_variation honours part labels and rejects non-compD", {
  cd <- declare_composition(timeuse, parts, closure = 1440,
                            labels = c(MVPA = "Exercise"))
  p <- plot_variation(cd)
  x_labels <- levels(ggplot2::ggplot_build(p)$plot$data$col)
  expect_true("Exercise" %in% x_labels)
  expect_error(plot_variation(list()), "compD")
})

test_that("plot_variation applies a ggplot2 theme and rejects a non-theme", {
  cd <- declare_composition(timeuse, c("MVPA", "LPA", "SB", "SPT"), closure = 1440)
  p <- plot_variation(cd, theme = ggplot2::theme(legend.position = "bottom"))
  expect_s3_class(p, "ggplot")
  expect_identical(p$theme$legend.position, "bottom")
  expect_error(plot_variation(cd, theme = "theme_bw"), "theme")
})
