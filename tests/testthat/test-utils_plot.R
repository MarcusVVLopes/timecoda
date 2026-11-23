# The shared plot helpers in R/utils_plot.R: tc_triplet_plots() and the
# tc_draw_plots() renderer's recovery from ggtern theme-element pollution.

parts <- c("MVPA", "LPA", "SB", "SPT")

test_that("plot_realloc renders despite ggtern theme-element pollution", {
  skip_if_not_installed("withr")
  # Simulate ggtern having globally registered a ternary theme element whose
  # declared type (rel) the installed ggplot2 rejects when building any plot.
  ggplot2::register_theme_elements(
    tern.axis.ticks.length.major = grid::unit(1, "pt"),
    element_tree = list(
      tern.axis.ticks.length.major = ggplot2::el_def("rel")
    )
  )
  withr::defer(ggplot2::reset_theme_settings())

  cd <- declare_composition(timeuse, parts, closure = 1440)
  fit <- comp_iv(cd, "obesity", model = "glm", family = "binomial")

  # The shared renderer should clear the bad registration and draw successfully
  # (onto setup.R's throwaway device), so no "Could not render" warning fires.
  expect_no_warning(p <- plot_realloc(fit, delta_max = 20))
  expect_s3_class(p, "realloc_plots")
  expect_s3_class(p$figure, "patchwork")
})

# --- Scenario (d): ggtern theme pollution across real plot calls -------------

test_that("plot_realloc and plot_variation still render after a real plot_ternary() call", {
  withr::defer(ggplot2::reset_theme_settings())

  small_parts <- c("MVPA", "SB", "SPT")
  cd <- declare_composition(timeuse, small_parts, closure = 1440)

  # Draw a real ternary plot first -- this is what registers ggtern's theme
  # elements in the running session (the actual pollution source, as opposed
  # to the synthetic registration used in the test above).
  invisible(quiet_draw(plot_ternary(cd, conf_level = NULL)))

  fit <- comp_iv(cd, "obesity", model = "glm", family = "binomial")

  # The shared tc_draw_plots() renderer should clear any bad
  # ggtern theme-element registration and draw successfully.
  p <- plot_realloc(fit, delta_max = 20, n_points = 11)
  expect_s3_class(p, "realloc_plots")
  expect_s3_class(p$figure, "patchwork")
  expect_s3_class(ggplot2::ggplot_build(p$panels[["MVPA.SB"]]), "ggplot_built")

  v <- plot_variation(cd)
  expect_s3_class(v, "ggplot")
  expect_s3_class(ggplot2::ggplot_build(v), "ggplot_built")
})

test_that("tc_triplet_plots() returns a bare plot for one triplet and a named list otherwise", {
  panel <- function(tri) paste(tri, collapse = "+")
  expect_identical(tc_triplet_plots(parts[1:3], panel), "MVPA+LPA+SB")
  many <- tc_triplet_plots(parts, panel)
  expect_length(many, choose(length(parts), 3))
  expect_identical(names(many)[1], "MVPA.LPA.SB")
})

test_that("every ternary panel draws its corner titles unclipped", {
  cd <- declare_composition(timeuse, parts, closure = 1440)
  clip_of <- function(panels) purrr::map_chr(panels, function(p) p$coordinates$clip)
  expect_identical(unname(clip_of(tc_build_ternary(cd, conf_level = NULL))),
                   rep("off", 4))
  fit <- comp_iv(cd, "zbmi", model = "lm")
  surface <- quiet_draw(plot(fit, resolution = 10, quantiles = NULL))
  expect_identical(unname(clip_of(surface)), rep("off", 4))
  # Each panel has its own coord, so a panel built alone keeps the default.
  alone <- ggtern::ggtern(timeuse, ggplot2::aes(x = MVPA, y = SB, z = SPT))
  expect_identical(alone$coordinates$clip, "on")
})
