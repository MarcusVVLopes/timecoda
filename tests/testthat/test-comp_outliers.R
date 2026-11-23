parts <- c("MVPA", "LPA", "SB", "SPT")

test_that("comp_outliers returns the documented table shape", {
  set.seed(1)
  tbl <- comp_outliers(timeuse, parts)
  expect_s3_class(tbl, "tbl_df")
  expect_named(tbl, c(".row", "distance", "limit", "is_outlier"))
  expect_equal(nrow(tbl), nrow(timeuse))
  expect_identical(tbl$.row, seq_len(nrow(timeuse)))
  expect_type(tbl$is_outlier, "logical")
  expect_false(anyNA(tbl$distance))
})

test_that("is_outlier is exactly distance > limit", {
  set.seed(1)
  tbl <- comp_outliers(timeuse, parts)
  expect_identical(tbl$is_outlier, tbl$distance > tbl$limit)
})

test_that("the cutoff is sqrt(qchisq(quantile, df = D - 1)) for any D", {
  for (d in 3:5) {
    p <- paste0("p", seq_len(d))
    sim <- as.data.frame(matrix(runif(60 * d, 1, 10), ncol = d))
    names(sim) <- p
    for (q in c(0.95, 0.975)) {
      set.seed(1)
      tbl <- comp_outliers(sim, p, quantile = q)
      expect_equal(tbl$limit[1], sqrt(stats::qchisq(q, df = d - 1)))
      # A constant column, not an attribute.
      expect_length(unique(tbl$limit), 1)
    }
  }
})

test_that("the distances are the robust ones, not outCoDa's classical pair", {
  # robCompositions 2.6.0 returns both estimates; comp_outliers() reads the
  # robust fields, which must differ from the classical ones it discards.
  set.seed(1)
  ours <- comp_outliers(timeuse, parts)
  set.seed(1)
  raw <- robCompositions::outCoDa(as.data.frame(timeuse[parts]))
  expect_equal(ours$distance, as.numeric(raw$mahalDist))
  expect_identical(ours$is_outlier, as.logical(raw$outlierIndex))
  expect_false(isTRUE(all.equal(ours$distance, as.numeric(raw$m2))))
})

test_that("an extreme-ratio observation is flagged", {
  planted <- timeuse
  # Inside every marginal range, but an atypical MVPA:SPT ratio.
  planted[3, parts] <- as.list(c(MVPA = 230, LPA = 60, SB = 900, SPT = 250))
  set.seed(1)
  expect_true(comp_outliers(planted, parts)$is_outlier[3])
})

test_that("zeroes are rejected with a message naming comp_impute()", {
  zeroed <- timeuse
  zeroed$MVPA[1] <- 0
  expect_error(comp_outliers(zeroed, parts), "comp_impute")
  expect_error(comp_outliers(zeroed, parts), "MVPA")
})

test_that("too few observations abort instead of returning NA distances", {
  # outCoDa() returns NaN distances behind a mere warning below 2 * (D - 1).
  small <- as.data.frame(matrix(runif(6 * 5, 1, 100), nrow = 6))
  names(small) <- paste0("p", 1:5)
  expect_error(comp_outliers(small, names(small)), "Too few observations")
  # And it works from 2 * (D - 1) = 8 rows up.
  ok <- as.data.frame(matrix(runif(8 * 5, 1, 100), nrow = 8))
  names(ok) <- paste0("p", 1:5)
  set.seed(1)
  expect_false(anyNA(comp_outliers(ok, names(ok))$distance))
})

test_that("... reaches outCoDa", {
  set.seed(1)
  strict <- comp_outliers(timeuse, parts, alpha = 0.75)
  set.seed(1)
  default <- comp_outliers(timeuse, parts)
  # alpha changes MCD's subset size, so the distances must move.
  expect_false(isTRUE(all.equal(strict$distance, default$distance)))
})

test_that("distances are invariant to closure and to the ILR pivot", {
  set.seed(1)
  base <- comp_outliers(timeuse, parts)$distance

  closed <- timeuse
  closed[parts] <- timeuse[parts] / rowSums(timeuse[parts])
  set.seed(1)
  expect_equal(comp_outliers(closed, parts)$distance, base)

  # Reordering the parts changes the pivot; the distance is a rotation
  # invariant, so it must not move.
  for (i in seq_along(parts)) {
    rotated <- c(parts[i], setdiff(parts, parts[i]))
    set.seed(1)
    expect_equal(comp_outliers(timeuse, rotated)$distance, base)
  }
})

test_that("coda cannot be overridden", {
  # On a closed composition outCoDa(coda = FALSE) hits a singular covariance
  # and fails inside robustbase; on unclosed data it would silently change the
  # cutoff's degrees of freedom from D-1 to D. Refused rather than ignored.
  expect_error(comp_outliers(timeuse, parts, coda = FALSE), "cannot be overridden")
  expect_error(comp_outliers(timeuse, parts, coda = TRUE), "cannot be overridden")
})
