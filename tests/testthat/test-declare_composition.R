parts <- c("MVPA", "LPA", "SB", "SPT")

test_that("declare_composition returns a compD object with the expected elements", {
  cd <- declare_composition(timeuse, parts, closure = 1440)
  expect_s3_class(cd, "compD")
  expect_named(
    cd,
    c("data_raw", "parts", "closure", "composition", "parts_imputed",
      "imput_table", "labels", "units", "impute_method")
  )
  expect_null(cd$units)
  # timeuse has no zeroes, so no method was actually invoked.
  expect_identical(cd$impute_method, NA_character_)
  # With no zeroes, `parts_imputed` is the input parts unchanged (original scale).
  expect_equal(
    as.data.frame(cd$parts_imputed), as.data.frame(timeuse[parts]),
    ignore_attr = TRUE
  )
  expect_identical(cd$parts, parts)
  expect_equal(cd$closure, 1440)
  expect_s3_class(cd$composition, "acomp")
  expect_setequal(
    names(cd$imput_table),
    c("part", "n", "n_zero", "pct_zero", "min_nonzero")
  )
})

test_that("the composition has D parts for any D", {
  for (d in 3:5) {
    p <- paste0("p", seq_len(d))
    sim <- as.data.frame(matrix(runif(40 * d, 1, 10), ncol = d))
    names(sim) <- p
    cd <- declare_composition(sim, p)
    expect_equal(ncol(unclass(cd$composition)), d)
  }
})

test_that("with no zeroes the composition matches acomp and no imputation occurs", {
  cd <- declare_composition(timeuse, parts, closure = 1440)
  expect_true(all(cd$imput_table$n_zero == 0))
  ref <- compositions::acomp(as.matrix(timeuse[parts]), total = 1440)
  expect_equal(as.numeric(cd$composition), as.numeric(ref), tolerance = 1e-8)
})

test_that("zeroes are summarized and imputed to strictly positive values", {
  z <- timeuse
  z$MVPA[1:5] <- 0
  cd <- declare_composition(z, parts, closure = 1440)
  expect_equal(cd$imput_table$n_zero[cd$imput_table$part == "MVPA"], 5L)
  comp_mat <- matrix(as.numeric(cd$composition), nrow = nrow(z))
  expect_true(all(comp_mat > 0))
  # `parts_imputed` is zero-free on the original scale: untouched rows keep
  # their exact raw values, imputed zeroes become small positive values.
  imp <- as.data.frame(cd$parts_imputed)
  expect_true(all(as.matrix(imp) > 0))
  expect_equal(imp[-(1:5), ], as.data.frame(z[parts])[-(1:5), ],
               ignore_attr = TRUE)
  expect_true(all(imp$MVPA[1:5] > 0 & imp$MVPA[1:5] < min(z$MVPA[z$MVPA > 0])))
})

test_that("impute_method is carried through when zeroes make it actually run", {
  z <- timeuse
  z$MVPA[1:5] <- 0
  cd <- declare_composition(z, parts, closure = 1440, impute_method = "multRepl")
  expect_identical(cd$impute_method, "multRepl")
})

test_that("extra arguments are forwarded to lrEM (e.g. a custom dl)", {
  z <- timeuse
  z$MVPA[1:5] <- 0
  # A user-supplied detection limit (one per part) is honoured by lrEM.
  custom_dl <- c(MVPA = 1, LPA = 1, SB = 1, SPT = 1)
  cd <- declare_composition(z, parts, closure = 1440, dl = custom_dl)
  expect_s3_class(cd, "compD")
  comp_mat <- matrix(as.numeric(cd$composition), nrow = nrow(z))
  expect_true(all(comp_mat > 0))
  # The custom dl changes the imputed values vs the auto-derived default.
  cd_auto <- declare_composition(z, parts, closure = 1440)
  expect_false(isTRUE(all.equal(
    as.numeric(cd$composition), as.numeric(cd_auto$composition)
  )))
})

test_that("labels are normalized to a complete named vector keyed by parts", {
  expect_null(declare_composition(timeuse, parts, closure = 1440)$labels)

  named <- declare_composition(timeuse, parts, closure = 1440,
                               labels = c(MVPA = "Exercise"))$labels
  expect_equal(named[["MVPA"]], "Exercise")
  expect_equal(named[["LPA"]], "LPA")           # unsupplied parts keep raw name
  expect_named(named, parts)

  positional <- declare_composition(
    timeuse, parts, closure = 1440,
    labels = c("a", "b", "c", "d")
  )$labels
  expect_equal(unname(positional), c("a", "b", "c", "d"))

  expect_error(
    declare_composition(timeuse, parts, labels = c(NOPE = "x")),
    "not a part"
  )
  expect_error(
    declare_composition(timeuse, parts, labels = c("only", "three", "here")),
    "one entry per part"
  )
})

test_that("declare_composition errors on missing parts, NA, and negative values", {
  expect_error(declare_composition(timeuse, c(parts, "nope")), "not found")
  na_d <- timeuse
  na_d$MVPA[1] <- NA
  expect_error(declare_composition(na_d, parts), "missing")
  neg <- timeuse
  neg$MVPA[1] <- -1
  expect_error(declare_composition(neg, parts), "negative")
  expect_error(declare_composition(timeuse, parts, closure = -1), "positive")
})

test_that("units must be a single non-empty string", {
  bad <- list(1, c("min", "day"), NA_character_, "", " ")
  for (u in bad) {
    expect_error(
      declare_composition(timeuse, parts, closure = 1440, units = u),
      "single non-empty string"
    )
  }
  expect_error(
    declare_composition(timeuse, parts, closure = 1440, units = ""),
    "single non-empty string (e.g. \"min/day\") or `NULL`",
    fixed = TRUE
  )
})

test_that("string closure is resolved and stored as a numeric", {
  # "percent" -> 100
  cd_prop <- declare_composition(timeuse, parts, closure = "percent")
  expect_equal(cd_prop$closure, 100)
  # Composition is closed to 100
  row_totals <- rowSums(matrix(as.numeric(cd_prop$composition),
                               nrow = nrow(timeuse)))
  expect_equal(row_totals, rep(100, nrow(timeuse)), tolerance = 1e-6)

  # "mean" -> mean of observed row totals
  expected_mean <- mean(rowSums(as.matrix(timeuse[parts])), na.rm = TRUE)
  cd_mean <- declare_composition(timeuse, parts, closure = "mean")
  expect_equal(cd_mean$closure, expected_mean, tolerance = 1e-8)
  row_totals_m <- rowSums(matrix(as.numeric(cd_mean$composition),
                                 nrow = nrow(timeuse)))
  expect_equal(row_totals_m, rep(expected_mean, nrow(timeuse)), tolerance = 1e-6)

  # "nonsense" errors
  expect_error(
    declare_composition(timeuse, parts, closure = "nonsense"),
    "percent.*mean|mean.*percent"
  )
})

test_that("without detect_outliers nothing outlier-related is added", {
  cd <- declare_composition(timeuse, parts, closure = 1440)
  # The object keeps its usual shape: flagging only ever adds data_raw columns.
  expect_named(cd$data_raw, names(timeuse))
})

test_that("detect_outliers appends the three columns to data_raw", {
  set.seed(1)
  cd <- declare_composition(timeuse, parts, closure = 1440,
                            detect_outliers = TRUE)
  expect_named(
    cd$data_raw,
    c(names(timeuse), ".outlier_distance", ".outlier_limit", ".outlier")
  )
  expect_type(cd$data_raw$.outlier, "logical")
  expect_true(any(cd$data_raw$.outlier))
  # The flag is exactly the distance against the cutoff.
  expect_identical(
    cd$data_raw$.outlier,
    cd$data_raw$.outlier_distance > cd$data_raw$.outlier_limit
  )
  # And equals what the detector returns when called directly.
  set.seed(1)
  direct <- comp_outliers(timeuse, parts)
  expect_equal(cd$data_raw$.outlier_distance, direct$distance)
  expect_identical(cd$data_raw$.outlier, direct$is_outlier)
})

test_that("flagging never drops rows, so every field stays row-aligned", {
  set.seed(1)
  cd <- declare_composition(timeuse, parts, closure = 1440,
                            detect_outliers = TRUE)
  expect_equal(nrow(cd$data_raw), nrow(timeuse))
  expect_equal(nrow(cd$parts_imputed), nrow(timeuse))
  expect_equal(nrow(cd$composition), nrow(timeuse))
  # The parts themselves are untouched by flagging.
  expect_equal(as.data.frame(cd$data_raw[parts]),
               as.data.frame(timeuse[parts]), ignore_attr = TRUE)
})

test_that("flagging does not change the composition or its closure", {
  plain <- declare_composition(timeuse, parts, closure = "mean")
  set.seed(1)
  flagged <- declare_composition(timeuse, parts, closure = "mean",
                                 detect_outliers = TRUE)
  expect_equal(flagged$closure, plain$closure)
  expect_equal(unclass(flagged$composition), unclass(plain$composition))
})

test_that("a column name clash is rejected rather than overwritten", {
  clashing <- timeuse
  clashing$.outlier <- FALSE
  expect_error(
    declare_composition(clashing, parts, detect_outliers = TRUE),
    "already contains column"
  )
})

test_that("stale .outlier columns are refused even without detection", {
  # Re-declaring on rows filtered out of an earlier flagged composition, but
  # forgetting to drop its columns, used to carry them through silently:
  # tc_has_outliers() then read TRUE and summary() reported a cutoff for a
  # computation that never ran.
  set.seed(1)
  cd <- declare_composition(timeuse, parts, closure = 1440,
                            detect_outliers = TRUE)
  kept <- cd$data_raw[!cd$data_raw$.outlier, , drop = FALSE]

  expect_error(declare_composition(kept, parts, closure = 1440),
               "already contains column")
  expect_error(declare_composition(kept, parts, closure = 1440),
               "starts_with")

  # Dropping them, as the documented recipe does, works.
  dropped <- kept[, !startsWith(names(kept), ".outlier"), drop = FALSE]
  cd2 <- declare_composition(dropped, parts, closure = 1440)
  expect_false(tc_has_outliers(cd2))
  expect_null(summary(cd2)$outliers)
})

test_that("a bad outlier_args is rejected by its user-facing name", {
  expect_error(
    declare_composition(timeuse, parts, detect_outliers = TRUE,
                        outlier_args = list(1)),
    "`outlier_args`"
  )
})

test_that("outlier_args reaches comp_outliers", {
  set.seed(1)
  cd <- declare_composition(timeuse, parts, closure = 1440,
                            detect_outliers = TRUE,
                            outlier_args = list(quantile = 0.99))
  expect_equal(cd$data_raw$.outlier_limit[1],
               sqrt(stats::qchisq(0.99, df = length(parts) - 1)))
})

test_that("outliers are flagged after imputation, so zeroes do not abort", {
  zeroed <- timeuse
  zeroed$MVPA[1:3] <- 0
  set.seed(1)
  cd <- declare_composition(zeroed, parts, closure = 1440,
                            detect_outliers = TRUE)
  expect_equal(nrow(cd$data_raw), nrow(zeroed))
  expect_false(anyNA(cd$data_raw$.outlier_distance))
})

test_that("print.compD reports the headline facts and returns invisibly", {
  cd <- declare_composition(timeuse, parts, closure = 1440)
  expect_output(print(cd), "<compD> compositional variable")
  expect_output(print(cd), "Parts: MVPA, LPA, SB, SPT")
  expect_output(print(cd), "Closure: 1440")
  expect_output(print(cd), "Observations: 300")
  expect_invisible(print(cd))
  # Silent about zeroes and outliers when there are none.
  out <- testthat::capture_output(print(cd))
  expect_no_match(out, "Zeroes")
  expect_no_match(out, "Outliers")
  expect_match(out, "Closure: 1440\n")
})

test_that("a declared unit label is stored and printed after the closure", {
  cd <- declare_composition(timeuse, parts, closure = 1440, units = "min/day")
  expect_identical(cd$units, "min/day")
  expect_output(print(cd), "Closure: 1440 min/day")
})

test_that("print.compD uses display labels when they were supplied", {
  cd <- declare_composition(timeuse, parts, closure = 1440,
                            labels = c(MVPA = "Exercise"))
  expect_output(print(cd), "Parts: Exercise, LPA, SB, SPT")
})

test_that("print.compD reports zeroes and flagged outliers when present", {
  z <- timeuse
  z$MVPA[1:5] <- 0
  set.seed(1)
  cd <- declare_composition(z, parts, closure = 1440,
                            impute_method = "multRepl",
                            detect_outliers = TRUE)
  expect_output(print(cd), "Zeroes imputed \\(multRepl\\): 5 across MVPA")
  expect_output(print(cd), "Outliers flagged: \\d+ \\(")
})

test_that("a bad closure still names the string alternatives", {
  # The positive-number check carries a hint naming "percent"/"mean".
  expect_error(declare_composition(timeuse, parts, closure = -1),
               "single positive number")
  expect_error(declare_composition(timeuse, parts, closure = -1),
               "percent.*mean|mean.*percent")
})

test_that("detect_outliers must be a single non-missing logical", {
  # isTRUE() alone treated every one of these as "no", silently skipping the
  # flagging the caller asked for.
  for (bad in list(1, c(TRUE, TRUE), "yes", NA)) {
    expect_error(
      declare_composition(timeuse, parts, closure = 1440,
                          detect_outliers = bad),
      "`detect_outliers`"
    )
  }
})

test_that("outlier_args is validated even when detection is off", {
  expect_error(
    declare_composition(timeuse, parts, closure = 1440,
                        detect_outliers = FALSE, outlier_args = list(1)),
    "`outlier_args`"
  )
})

test_that("duplicate display labels are refused", {
  expect_error(
    declare_composition(timeuse, parts, closure = 1440,
                        labels = c(MVPA = "Activity", LPA = "Activity")),
    "unique"
  )
  # A label equal to another part's raw name would merge two parts too.
  expect_error(
    declare_composition(timeuse, parts, closure = 1440, labels = c(MVPA = "SB")),
    "unique"
  )
})
