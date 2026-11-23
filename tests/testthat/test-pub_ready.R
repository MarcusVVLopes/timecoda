# R/pub_ready.R: pub_ready() (the summary()+pub_ready() split), and the
# compD builder's compact heading style it applies to every table.

parts <- c("MVPA", "LPA", "SB", "SPT")

# Checked before gt is needed, so these run without it.
test_that("pub_ready() refuses a digits that is not a single non-negative whole number", {
  s <- summary(declare_composition(timeuse, parts, closure = 1440))
  for (bad in list("two", -1, 2.5, c(1, 2), NA)) {
    expect_error(pub_ready(s, digits = bad),
                 "must be a single non-negative whole number", fixed = TRUE)
  }
})

test_that("pub_ready() refuses a p_flag that is not a single number from 0 to 1", {
  s <- summary(declare_composition(timeuse, parts, closure = 1440))
  for (bad in list("a", 2, -0.1, c(0.05, 0.01), NA_real_)) {
    expect_error(pub_ready(s, p_flag = bad),
                 "must be a single number from 0 to 1", fixed = TRUE)
  }
})

# gt is a Suggests dependency (display only).
skip_if_not_installed("gt")

test_that("pub_ready() builds the desc (Mean +/- SD) and var_matrix gt tables from the summary", {
  cd <- declare_composition(timeuse, parts, closure = 1440)
  s <- summary(cd)
  g <- pub_ready(s)
  expect_s3_class(g$desc, "gt_tbl")
  expect_s3_class(g$var_matrix, "gt_tbl")
  # desc is the compact Mean +/- SD + geometric-mean table, not the full
  # median/p25/p75/min/max spread -- that full spread still lives on the
  # tidy s$desc tibble, just not as a separate gt rendering.
  expect_setequal(names(g$desc$`_data`), c("part", "mean_sd", "geom_mean", "pct_total"))
})

test_that("pub_ready() omits a table the summary didn't compute", {
  cd <- declare_composition(timeuse, parts, closure = 1440)
  s <- summary(cd, var_matrix = FALSE)
  g <- pub_ready(s)
  expect_s3_class(g$desc, "gt_tbl")
  expect_null(g$var_matrix)

  s2 <- summary(cd, desc = FALSE)
  expect_null(pub_ready(s2)$desc)
})

test_that("pub_ready()'s desc table's % column follows the summary's closure", {
  cd_min <- declare_composition(timeuse, parts, closure = 1440)
  g_min <- pub_ready(summary(cd_min))
  expect_true("pct_total" %in% names(g_min$desc$`_data`))

  cd_prop <- declare_composition(timeuse, parts, closure = 1)
  g_prop <- pub_ready(summary(cd_prop))
  expect_false("pct_total" %in% names(g_prop$desc$`_data`))
  # At closure 100 the share would repeat the geometric mean, so it is dropped.
  cd_pct <- declare_composition(timeuse, parts, closure = "percent")
  expect_false("pct_total" %in% names(pub_ready(summary(cd_pct))$desc$`_data`))
})

test_that("pub_ready() errors for a class it has no case for", {
  expect_error(pub_ready(structure(list(), class = "not_a_timecoda_summary")))
})

test_that("pub_ready() applies the compact heading style to every table", {
  cd <- declare_composition(timeuse, parts, closure = 1440)
  g <- pub_ready(summary(cd))
  purrr::walk(g, \(tbl) {
    opts <- tbl[["_options"]]
    expect_identical(
      opts$value[[which(opts$parameter == "heading_title_font_size")]], "13px"
    )
    expect_identical(
      opts$value[[which(opts$parameter == "heading_subtitle_font_size")]], "11px"
    )
  })
})
