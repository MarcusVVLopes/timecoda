# Custom multi-part reallocation vectors ("best-day" style) on predict_realloc().

# A shared glm fit on the shipped timeuse composition.
custom_fit <- function(covariates = NULL) {
  cd <- declare_composition(timeuse, c("MVPA", "LPA", "SB", "SPT"),
                            closure = 1440)
  comp_iv(cd, outcome = "obesity", model = "glm", family = "binomial",
          covariates = covariates)
}

test_that("a single named vector returns a tidy realloc with a compact label", {
  fit <- custom_fit()
  r <- predict_realloc(fit, realloc = c(SPT = 30, SB = -15, LPA = -15))
  expect_s3_class(r, "realloc")
  expect_true(all(c("from", "to", "delta", "estimate", "conf.low",
                    "conf.high", "p.value") %in% names(r)))
  expect_equal(nrow(r), 1L)
  # to = the positive side, from = the negative side; delta = minutes moved.
  expect_equal(r$to, "SPT +30")
  expect_match(r$from, "SB -15")
  expect_match(r$from, "LPA -15")
  expect_equal(r$delta, 30)
  expect_true(is.finite(r$estimate))
  expect_identical(attr(r, "scheme"), "custom")
  # Each amount is formatted on its own, not padded to a common width.
  r <- predict_realloc(fit, realloc = c(SPT = 17.5, SB = -15, LPA = -2.5))
  expect_equal(r$from, "SB -15, LPA -2.5")
})

test_that("a list of scenarios stacks into one row per scenario", {
  fit <- custom_fit()
  r <- predict_realloc(
    fit,
    realloc = list(c(SPT = 30, SB = -30), c(MVPA = 20, LPA = -20))
  )
  expect_equal(nrow(r), 2L)
  expect_setequal(r$to, c("SPT +30", "MVPA +20"))
})

test_that("a data frame input gives one scenario per row", {
  fit <- custom_fit()
  df <- data.frame(MVPA = c(10, 20), SB = c(-10, -20))
  r <- predict_realloc(fit, realloc = df)
  expect_equal(nrow(r), 2L)
  expect_setequal(r$to, c("MVPA +10", "MVPA +20"))
  expect_equal(sort(r$delta), c(10, 20))
})

test_that("an imbalanced scenario errors and names the imbalance", {
  fit <- custom_fit()
  expect_error(
    predict_realloc(fit, realloc = c(SPT = 30, SB = -20)),
    "sum to 0"
  )
  expect_error(
    predict_realloc(fit, realloc = c(SPT = 30, SB = -20)),
    "10"
  )
})

test_that("an unknown part errors and names the offender", {
  fit <- custom_fit()
  expect_error(
    predict_realloc(fit, realloc = c(SPT = 30, NOTAPART = -30)),
    "NOTAPART"
  )
})

test_that("a scenario without both a gain and a loss errors", {
  fit <- custom_fit()
  # An all-zero scenario is balanced (sums to 0) but moves nothing, so it is
  # rejected for lacking both a gain and a loss.
  expect_error(
    predict_realloc(fit, realloc = c(SPT = 0, SB = 0)),
    "positive and one negative"
  )
})

test_that("custom mode is mutually exclusive with from/to/scheme", {
  fit <- custom_fit()
  expect_error(
    predict_realloc(fit, realloc = c(SPT = 30, SB = -30),
                    from = "SB", to = "SPT"),
    "mutually exclusive"
  )
  expect_error(
    predict_realloc(fit, realloc = c(SPT = 30, SB = -30),
                    scheme = "one-to-all"),
    "mutually exclusive"
  )
})

test_that("custom realloc reproduces the one-to-one estimate to 1e-8", {
  fit <- custom_fit()
  cust <- predict_realloc(fit, realloc = c(MVPA = 10, SB = -10))
  o2o  <- predict_realloc(fit, from = "SB", to = "MVPA", delta = 10)
  expect_equal(cust$estimate, o2o$estimate, tolerance = 1e-8)
  expect_equal(cust$conf.low, o2o$conf.low, tolerance = 1e-8)
  expect_equal(cust$conf.high, o2o$conf.high, tolerance = 1e-8)
  expect_equal(cust$p.value, o2o$p.value, tolerance = 1e-8)
})

test_that("summary() substitutes part labels into scenarios", {
  fit <- custom_fit()
  r <- predict_realloc(fit, realloc = c(SPT = 30, SB = -15, LPA = -15))
  s <- summary(r, labels = c(SPT = "Sleep", SB = "Sedentary"))
  expect_true("to_label" %in% names(s$table))
  # Part labels reach the compact custom strings.
  expect_equal(s$table$to_label, "Sleep +30")
  expect_match(s$table$from_label, "Sedentary -15")
  expect_match(s$table$from_label, "LPA -15")
})

test_that("pub_ready() renders the custom scenarios", {
  skip_if_not_installed("gt")
  fit <- custom_fit()
  r <- predict_realloc(fit, realloc = c(SPT = 30, SB = -15, LPA = -15))
  g <- pub_ready(summary(r, labels = c(SPT = "Sleep")))$effects
  expect_s3_class(g, "gt_tbl")
  expect_equal(g[["_data"]]$to_label, "Sleep +30")
})

test_that("a scenario's error names its position among several", {
  fit <- custom_fit()
  expect_error(
    predict_realloc(fit, realloc = list(c(SPT = 30, SB = -30), c(SB = 5))),
    "(scenario 2)", fixed = TRUE
  )
  expect_error(
    predict_realloc(fit, realloc = list(c(SPT = 30, SB = -30),
                                        c(FOO = 5, SB = -5))),
    '`realloc` (scenario 2) names "FOO", which is not', fixed = TRUE
  )
})

test_that("a part named twice is reported, singly or in the plural", {
  fit <- custom_fit()
  expect_error(
    predict_realloc(fit, realloc = c(SB = 10, SB = -10, MVPA = 0)),
    'names part "SB" more than once', fixed = TRUE
  )
  expect_error(
    predict_realloc(fit, realloc = c(SB = 10, SB = -10, MVPA = 5, MVPA = -5)),
    'names parts "SB" and "MVPA" more than once', fixed = TRUE
  )
})

test_that("by = works group-wise in custom mode", {
  fit <- custom_fit(covariates = "sex")
  r <- predict_realloc(fit, realloc = c(SPT = 30, SB = -30), by = "sex")
  expect_true("sex" %in% names(r))
  expect_setequal(as.character(r$sex), levels(timeuse$sex))
  expect_equal(nrow(r), nlevels(timeuse$sex))
})

test_that("multinomial fits keep their per-category response_level in custom mode", {
  cd <- declare_composition(timeuse, c("MVPA", "LPA", "SB", "SPT"),
                            closure = 1440)
  cd$data_raw$grade_u <- factor(timeuse$grade, ordered = FALSE)
  fit <- comp_iv(cd, outcome = "grade_u", model = "mblogit")
  r <- predict_realloc(fit, realloc = c(SPT = 30, SB = -30))
  expect_true("response_level" %in% names(r))
  expect_gt(nrow(r), 1L)
  expect_true(all(r$to == "SPT +30"))
})

test_that("on_invalid handles a scenario that empties a part", {
  fit <- custom_fit()
  # SB has ~ a few hundred minutes at the geometric mean; moving far more than
  # it holds drives it non-positive.
  big <- c(SPT = 2000, SB = -2000)
  expect_error(predict_realloc(fit, realloc = big, on_invalid = "error"),
               "non-positive")
  na_out <- predict_realloc(fit, realloc = big)
  expect_true(is.na(na_out$estimate))
})
