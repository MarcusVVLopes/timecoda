parts <- c("MVPA", "LPA", "SB", "SPT")

test_that("summary.compD returns a visible summary_compD object with tidy tibbles", {
  cd <- declare_composition(timeuse, parts, closure = 1440)
  s <- summary(cd)
  expect_s3_class(s, "summary_compD")
  expect_named(
    s,
    c("desc", "var_matrix", "parts", "labels", "closure", "units", "row_totals",
      "n_obs", "imput_table", "impute_method", "outliers")
  )
  expect_s3_class(s$desc, "tbl_df")
  expect_setequal(
    names(s$desc),
    c("part", "label", "geom_mean", "mean", "sd", "median",
      "p25", "p75", "min", "max")
  )
  # part name, display label, then the geometric mean (compositional centre)
  expect_equal(names(s$desc)[1:3], c("part", "label", "geom_mean"))
  expect_equal(nrow(s$desc), length(parts))
  expect_s3_class(s$var_matrix, "tbl_df")
  expect_equal(nrow(s$var_matrix), length(parts))
})

test_that("summary.compD honours the desc/var_matrix content flags", {
  cd <- declare_composition(timeuse, parts, closure = 1440)
  expect_null(summary(cd, var_matrix = FALSE)$var_matrix)
  expect_null(summary(cd, desc = FALSE)$desc)
})

test_that("print.summary_compD prints a text report and returns invisibly", {
  cd <- declare_composition(timeuse, parts, closure = 1440)
  s <- summary(cd)
  expect_output(print(s), "compositional variable summary")
  expect_output(print(s), "Row totals \\(before closure\\)")
  expect_output(print(s), "No zeroes detected")
  expect_output(print(s), "Descriptive statistics")
  expect_output(print(s), "Variation matrix")
  expect_invisible(print(s))
  expect_match(testthat::capture_output(print(s)), "Closure: 1440\n")
})

test_that("print.summary_compD shows a declared unit label after the closure", {
  cd <- declare_composition(timeuse, parts, closure = 1440, units = "min/day")
  s <- summary(cd)
  expect_identical(s$units, "min/day")
  expect_output(print(s), "Closure: 1440 min/day")
})

test_that("print.summary_compD reports which parts had zeroes and the imputation method", {
  z <- timeuse
  z$MVPA[1:5] <- 0
  cd <- declare_composition(z, parts, closure = 1440, impute_method = "multRepl")
  s <- summary(cd)
  expect_output(print(s), "Zeroes detected.*MVPA")
  expect_output(print(s), "Zeroes imputed via multRepl")
  expect_no_match(
    testthat::capture_output(print(s)), "No zeroes detected"
  )
})

test_that("summary(cd)'s row_totals reflect the raw pre-closure scale", {
  cd <- declare_composition(timeuse, parts, closure = 1440)
  s <- summary(cd)
  raw <- rowSums(timeuse[parts])
  expect_equal(s$row_totals$mean, mean(raw))
  expect_equal(s$row_totals$sd, stats::sd(raw))
  expect_equal(s$row_totals$min, min(raw))
  expect_equal(s$row_totals$max, max(raw))
})

test_that("summary desc spread statistics describe the original data, zeroes included", {
  z <- timeuse
  z$MVPA[1:5] <- 0
  cd <- declare_composition(z, parts, closure = 1440, impute_method = "multRepl")
  d <- summary(cd)$desc
  expect_equal(d$mean, unname(colMeans(z[parts])))
  expect_equal(d$min[d$part == "MVPA"], 0)
  expect_equal(d$max, unname(purrr::map_dbl(z[parts], max)))
  q <- purrr::map(z[parts], \(v) stats::quantile(v, c(0.25, 0.75), names = FALSE))
  expect_equal(d$p25, unname(purrr::map_dbl(q, 1)))
  expect_equal(d$p75, unname(purrr::map_dbl(q, 2)))
  # the geometric mean stays on the closed, imputed composition
  expect_equal(sum(d$geom_mean), 1440, tolerance = 1e-6)
})

test_that("summary desc spread statistics stay on the recorded scale under a percent closure", {
  cd <- declare_composition(timeuse, parts, closure = "percent")
  d <- summary(cd)$desc
  expect_equal(d$mean, unname(colMeans(timeuse[parts])))
  expect_equal(d$median, unname(purrr::map_dbl(timeuse[parts], stats::median)))
  expect_equal(sum(d$geom_mean), 100, tolerance = 1e-6)
})

test_that("print.summary_compD states the scale of the descriptive columns", {
  cd <- declare_composition(timeuse, parts, closure = 1440)
  out <- testthat::capture_output(print(summary(cd)))
  expect_match(out, "Descriptive statistics\n  geom_mean: imputed, closed composition; other columns: original data\n")
})

test_that("labels are reported in the desc and var_matrix tibbles", {
  cd <- declare_composition(timeuse, parts, closure = 1440,
                            labels = c(MVPA = "Exercise"))
  s <- summary(cd)
  expect_equal(s$desc$label[s$desc$part == "MVPA"], "Exercise")
  # parts without a supplied label fall back to their raw name
  expect_equal(s$desc$label[s$desc$part == "LPA"], "LPA")
  expect_true("label" %in% names(s$var_matrix))
})

test_that("summary desc sits on the object's own closure", {
  # Numeric closure stored at declare time
  cd_1440 <- declare_composition(timeuse, parts, closure = 1440)
  s_1440 <- summary(cd_1440)$desc
  expect_equal(sum(s_1440$geom_mean), 1440, tolerance = 1e-6)

  # String closure resolved at declare time propagates through to summary
  cd_prop <- declare_composition(timeuse, parts, closure = "percent")
  s_prop <- summary(cd_prop)$desc
  expect_equal(sum(s_prop$geom_mean), 100, tolerance = 1e-6)

  cd_mean <- declare_composition(timeuse, parts, closure = "mean")
  s_mean <- summary(cd_mean)$desc
  expect_equal(sum(s_mean$geom_mean),
               mean(rowSums(timeuse[parts])), tolerance = 1e-6)
})

# pub_ready()'s own tests live in test-pub_ready.R (R/pub_ready.R).

test_that("print.summary_compD says nothing about outliers when none were sought", {
  cd <- declare_composition(timeuse, parts, closure = 1440)
  s <- summary(cd)
  expect_null(s$outliers)
  expect_no_match(
    testthat::capture_output(print(s)), "Outlier", ignore.case = TRUE
  )
  expect_output(print(s), "Observations: 300")
})

test_that("print.summary_compD reports flagged outliers with the cutoff", {
  set.seed(1)
  cd <- declare_composition(timeuse, parts, closure = 1440,
                            detect_outliers = TRUE)
  s <- summary(cd)
  expect_equal(s$n_obs, nrow(timeuse))
  expect_named(s$outliers, c("n_outlier", "pct", "limit"))
  expect_equal(s$outliers$n_outlier, sum(cd$data_raw$.outlier))
  expect_equal(s$outliers$limit, cd$data_raw$.outlier_limit[1])
  expect_output(print(s), "Outliers flagged: \\d+ of 300")
  expect_output(print(s), "robust Mahalanobis cutoff")
  # The chi-squared probability is recovered from the cutoff, not stored.
  expect_output(print(s), "chi-squared 0\\.975")
})
