# plot_ternary(), plot_variation() and their builders in R/plot_ternary.R.

parts <- c("MVPA", "LPA", "SB", "SPT")

test_that("plot_ternary returns one panel per part triplet", {
  cd3 <- declare_composition(timeuse, parts[1:3], closure = 1440)
  single <- quiet_draw(plot_ternary(cd3, conf_level = NULL))
  expect_s3_class(single, "ggplot")

  cd <- declare_composition(timeuse, parts, closure = 1440)
  many <- quiet_draw(plot_ternary(cd, group = "sex", conf_level = NULL))
  expect_type(many, "list")
  expect_length(many, choose(length(parts), 3))
  expect_s3_class(many[[1]], "ggplot")
})

test_that("plot_ternary applies a ggtern theme to every panel and rejects a non-theme", {
  # Built without drawing: rendering under ggtern's theme elements emits
  # hundreds of the known ggtern/ggplot2 version-mismatch warnings.
  cd <- declare_composition(timeuse, parts[1:3], closure = 1440)
  p <- tc_build_ternary(cd, conf_level = NULL, theme = ggtern::theme_bw())
  expect_s3_class(p, "ggplot")
  expect_error(quiet_draw(plot_ternary(cd, theme = "theme_bw")), "theme")
  many <- declare_composition(timeuse, parts, closure = 1440)
  themed <- tc_build_ternary(many, conf_level = NULL, theme = ggtern::theme_dark())
  expect_length(themed, choose(length(parts), 3))
  expect_s3_class(themed[[1]], "ggplot")
})

test_that("plot_ternary accepts only a compD and validates group against its data", {
  cd <- declare_composition(timeuse, parts, closure = 1440)
  expect_error(plot_ternary(timeuse), "compD")
  expect_error(quiet_draw(plot_ternary(cd, group = "not_a_column")), "group")
  expect_error(plot_ternary(cd, group = c("sex", "sex")), "group")
})

test_that("ternary plots give a specific error for fewer than three parts", {
  cd2 <- declare_composition(timeuse, parts[1:2], closure = 1440)
  expect_error(quiet_draw(plot_ternary(cd2)), "at least three composition parts")
  expect_error(quiet_draw(plot(cd2, histogram = FALSE, boxplot = FALSE, variation = FALSE)),
               "at least three composition parts")
})

test_that("tc_check_theme() accepts NULL or a theme object and names the package", {
  expect_null(tc_check_theme(NULL, "ggplot2"))
  expect_null(tc_check_theme(ggplot2::theme_bw(), "ggplot2"))
  expect_error(tc_check_theme("theme_bw", "ggtern"), "ggtern")
})

# Labels and colours shared across the plot builders --------------------------

test_that("ternary axis titles use the compD display labels", {
  # Found by the combined read: every other plot builder resolved labels, so a
  # single quiet_draw(plot(cd)) showed "Moderate-to-vigorous PA" on the heatmap and "MVPA"
  # on the ternary.
  labs <- c(MVPA = "Moderate-to-vigorous PA", SB = "Sedentary behaviour")
  cd <- declare_composition(timeuse, parts, closure = 1440, labels = labs)
  p <- tc_build_ternary(cd, conf_level = NULL)[[1]]
  expect_identical(p$labels$x, "Moderate-to-vigorous PA")
  expect_identical(p$labels$z, "Sedentary behaviour")
  # A part with no supplied label keeps its raw name.
  expect_identical(p$labels$y, "LPA")
})

test_that("ternary axis titles fall back to part names without labels", {
  cd <- declare_composition(timeuse, parts, closure = 1440)
  p <- tc_build_ternary(cd, conf_level = NULL)[[1]]
  expect_identical(
    c(p$labels$x, p$labels$y, p$labels$z), c("MVPA", "LPA", "SB")
  )
})

test_that("heatmap value text stays readable when low/high are inverted", {
  # The text colour used to be keyed to the value rather than to the rendered
  # fill, so an inverted gradient put white text on a white tile.
  cd <- declare_composition(timeuse, parts, closure = 1440)
  comp_df <- tc_comp_df(cd)

  dark_low <- tc_variation_heatmap(comp_df, parts, parts,
                                   low = "#132B43", high = "#EAF3FB")
  light_low <- tc_variation_heatmap(comp_df, parts, parts,
                                    low = "#FFFFFF", high = "#132B43")
  small <- function(p) {
    v <- p$data$value
    unique(p$data$txt_col[v < mean(range(v))])
  }
  # Smallest cells are dark under the default gradient and light under the
  # inverted one, so the text colour flips with them.
  expect_identical(small(dark_low), "white")
  expect_identical(small(light_low), "grey10")
})

test_that("plot_variation draws as a side effect and returns its ggplot invisibly", {
  cd <- declare_composition(timeuse, parts, closure = 1440)
  out <- withVisible(plot_variation(cd))
  expect_false(out$visible)
  expect_s3_class(out$value, "ggplot")
})

test_that("ternary plots take one or more levels in (0, 1) and refuse others", {
  cd <- declare_composition(timeuse, parts[1:3], closure = 1440)
  refused <- "conf_level.*one or more numbers in \\(0, 1\\), or .*NULL"
  expect_error(plot_ternary(cd, conf_level = 95), refused)
  expect_error(plot_ternary(cd, conf_level = c(0.5, 1.5)), refused)
  expect_error(plot_ternary(cd, conf_level = c(0.5, NA)), refused)
  expect_error(plot_ternary(cd, conf_level = "0.95"), refused)
  expect_error(plot(cd, histogram = FALSE, boxplot = FALSE, variation = FALSE,
                    conf_level = 0),
               refused)
  # Several levels draw one region each.
  p <- quiet_draw(plot_ternary(cd, conf_level = c(0.5, 0.95)))
  expect_s3_class(p, "ggplot")
  region <- purrr::keep(p$layers, function(l) inherits(l$geom, "GeomConfidenceTern"))
  expect_length(region, 1)
  expect_identical(region[[1]]$stat_params$breaks, c(0.5, 0.95))
})
