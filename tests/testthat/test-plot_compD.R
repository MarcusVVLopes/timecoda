parts <- c("MVPA", "LPA", "SB", "SPT")

test_that("plot.compD returns the requested ggplot components", {
  cd <- declare_composition(timeuse, parts, closure = 1440)
  plots <- quiet_draw(plot(cd))
  expect_named(plots, c("histogram", "boxplot", "ternary", "variation"))
  expect_s3_class(plots$histogram, "ggplot")
  expect_s3_class(plots$boxplot, "ggplot")
  expect_s3_class(plots$variation, "ggplot")
})

test_that("the ternary element has choose(D, 3) panels", {
  for (d in c(3, 4)) {
    p <- paste0("p", seq_len(d))
    sim <- as.data.frame(matrix(runif(60 * d, 1, 10), ncol = d))
    names(sim) <- p
    cd <- declare_composition(sim, p)
    tern <- quiet_draw(plot(cd))$ternary
    n_expected <- choose(d, 3)
    if (n_expected == 1) {
      expect_s3_class(tern, "ggplot")
    } else {
      expect_length(tern, n_expected)
    }
  }
})

test_that("toggling flags drops components", {
  cd <- declare_composition(timeuse, parts, closure = 1440)
  expect_named(quiet_draw(plot(cd, histogram = FALSE)), c("boxplot", "ternary", "variation"))
  expect_named(quiet_draw(plot(cd, ternary = FALSE)), c("histogram", "boxplot", "variation"))
  expect_named(quiet_draw(plot(cd, variation = FALSE)), c("histogram", "boxplot", "ternary"))
})

test_that("the outlier panels appear only when outliers were flagged", {
  cd <- declare_composition(timeuse, parts, closure = 1440)
  expect_named(quiet_draw(plot(cd)), c("histogram", "boxplot", "ternary", "variation"))

  set.seed(1)
  flagged <- declare_composition(timeuse, parts, closure = 1440,
                                 detect_outliers = TRUE)
  p <- quiet_draw(plot(flagged))
  expect_named(p, c("histogram", "boxplot", "ternary", "variation",
                    "outliers", "outlier_ternary"))
  expect_s3_class(p$outliers, "ggplot")
  expect_length(p$outlier_ternary, choose(length(parts), 3))
})

test_that("outliers = FALSE drops both outlier panels", {
  set.seed(1)
  cd <- declare_composition(timeuse, parts, closure = 1440,
                            detect_outliers = TRUE)
  expect_named(quiet_draw(plot(cd, outliers = FALSE)),
               c("histogram", "boxplot", "ternary", "variation"))
})

test_that("the outlier ternary has choose(D, 3) panels for any D", {
  for (d in c(3, 4)) {
    pn <- paste0("p", seq_len(d))
    sim <- as.data.frame(matrix(runif(60 * d, 1, 10), ncol = d))
    names(sim) <- pn
    set.seed(1)
    cd <- declare_composition(sim, pn, detect_outliers = TRUE)
    tern <- quiet_draw(plot(cd))$outlier_ternary
    if (choose(d, 3) == 1) {
      expect_s3_class(tern, "ggplot")
    } else {
      expect_length(tern, choose(d, 3))
    }
  }
})

test_that("ternary = FALSE drops the outlier ternary but keeps the distance panel", {
  set.seed(1)
  cd <- declare_composition(timeuse, parts, closure = 1440,
                            detect_outliers = TRUE)
  expect_named(quiet_draw(plot(cd, ternary = FALSE)),
               c("histogram", "boxplot", "variation", "outliers"))
})

test_that("a two-part flagged composition plots with ternary = FALSE", {
  # A ternary needs three parts, so `ternary = FALSE` is how a two-part
  # composition is plotted at all; the outlier panels must not reintroduce
  # the abort that gating fixes.
  set.seed(1)
  two <- data.frame(a = runif(40, 1, 10), b = runif(40, 1, 10))
  cd <- declare_composition(two, c("a", "b"), detect_outliers = TRUE)
  p <- quiet_draw(plot(cd, ternary = FALSE))
  expect_named(p, c("histogram", "boxplot", "variation", "outliers"))
  expect_s3_class(p$outliers, "ggplot")
  # And it still aborts when a ternary really is asked for.
  expect_error(quiet_draw(plot(cd)), "at least three composition parts")
})

test_that("the outlier ternary titles its axes with the compD display labels", {
  set.seed(1)
  cd <- declare_composition(
    timeuse, parts, closure = 1440, detect_outliers = TRUE,
    labels = c(MVPA = "Moderate-to-vigorous PA", SB = "Sedentary behaviour")
  )
  p <- quiet_draw(plot(cd, histogram = FALSE, boxplot = FALSE, variation = FALSE))
  first <- p$outlier_ternary[[1]]  # the MVPA.LPA.SB triplet
  expect_identical(first$labels$x, "Moderate-to-vigorous PA")
  expect_identical(first$labels$y, "LPA")
  expect_identical(first$labels$z, "Sedentary behaviour")
})
