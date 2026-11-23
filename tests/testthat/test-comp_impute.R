parts <- c("MVPA", "LPA", "SB", "SPT")

test_that("comp_impute returns parts_imputed + imput_table with the expected shape", {
  imp <- comp_impute(timeuse, parts)
  expect_type(imp, "list")
  expect_named(imp, c("parts_imputed", "imput_table", "impute_method"))
  expect_setequal(
    names(imp$imput_table),
    c("part", "n", "n_zero", "pct_zero", "min_nonzero")
  )
  expect_equal(names(imp$parts_imputed), parts)
})

test_that("with no zeroes, parts_imputed equals the input parts unchanged and impute_method is NA", {
  imp <- comp_impute(timeuse, parts)
  expect_true(all(imp$imput_table$n_zero == 0))
  expect_equal(as.data.frame(imp$parts_imputed), as.data.frame(timeuse[parts]),
               ignore_attr = TRUE)
  # No method was actually invoked, so echoing the requested one would mislead.
  expect_identical(imp$impute_method, NA_character_)
})

test_that("default method (lrEM) matches zCompositions::lrEM directly", {
  z <- timeuse
  z$MVPA[1:5] <- 0
  imp <- comp_impute(z, parts)

  part_mat <- as.matrix(z[parts])
  min_nonzero <- vapply(parts, function(p) {
    col <- part_mat[, p]
    min(col[col > 0])
  }, numeric(1))
  dl <- matrix(rep(min_nonzero, each = nrow(z)), nrow = nrow(z),
               dimnames = list(NULL, parts))
  ref <- zCompositions::lrEM(part_mat, label = 0, dl = dl,
                             suppress.print = TRUE)

  expect_equal(as.matrix(imp$parts_imputed), unname(as.matrix(ref)),
               ignore_attr = TRUE)
})

test_that("every supported method runs, returns strictly positive values, and reports itself as used", {
  z <- timeuse
  z$MVPA[1:5] <- 0
  for (m in c("lrEM", "multRepl", "lrDA", "multLN", "lrSVD")) {
    imp <- comp_impute(z, parts, method = m)
    expect_true(all(as.matrix(imp$parts_imputed) > 0), info = m)
    expect_equal(names(imp$parts_imputed), parts, info = m)
    expect_equal(nrow(imp$parts_imputed), nrow(z), info = m)
    expect_identical(imp$impute_method, m, info = m)
  }
})

test_that("an unsupported method is rejected with an informative error", {
  expect_error(comp_impute(timeuse, parts, method = "not_a_method"))
})

test_that("an entirely-zero part is rejected before any imputation runs", {
  z <- timeuse
  z$MVPA <- 0
  expect_error(comp_impute(z, parts), "entirely zero")
})

test_that("a user-supplied dl overrides the derived default", {
  z <- timeuse
  z$MVPA[1:5] <- 0
  custom_dl <- c(MVPA = 1, LPA = 1, SB = 1, SPT = 1)
  imp_custom <- comp_impute(z, parts, dl = custom_dl)
  imp_auto <- comp_impute(z, parts)
  expect_false(isTRUE(all.equal(
    as.matrix(imp_custom$parts_imputed), as.matrix(imp_auto$parts_imputed)
  )))
})

test_that("comp_impute() rejects missing and negative values with a specific error", {
  d <- as.data.frame(timeuse[c("MVPA", "LPA", "SB", "SPT")])
  d_na <- d
  d_na$MVPA[1] <- NA
  expect_error(comp_impute(d_na, names(d)), "missing values")
  d_neg <- d
  d_neg$LPA[2] <- -1
  expect_error(comp_impute(d_neg, names(d)), "negative")
})

test_that("a row dropped by zCompositions aborts rather than misaligning", {
  # z.delete = TRUE silently deletes rows that are almost entirely zero,
  # behind a warning, which would leave composition/parts_imputed shorter
  # than data_raw and break every downstream column-bind.
  one <- timeuse
  one[5, parts] <- 0
  expect_error(suppressWarnings(comp_impute(one, parts)), "of 300 rows")
  expect_error(suppressWarnings(comp_impute(one, parts)),
               "Row 5 is entirely zero")

  two <- one
  two[9, parts] <- 0
  expect_error(suppressWarnings(comp_impute(two, parts)),
               "Rows 5, 9 are entirely zero")

  # And through the declare pipeline, where the misalignment would land.
  expect_error(suppressWarnings(declare_composition(one, parts, closure = 1440)),
               "of 300 rows")
})
