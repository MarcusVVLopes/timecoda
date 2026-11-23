test_that("timeuse fixture has the expected structure", {
  data(timeuse, package = "timecoda", envir = environment())

  expect_s3_class(timeuse, "data.frame")
  expect_identical(nrow(timeuse), 300L)
  expect_named(
    timeuse,
    c("id", "sex", "school", "region", "MVPA", "LPA", "SB", "SPT",
      "zbmi", "obesity", "time", "event", "grade")
  )
})

test_that("the four parts sum to 1440 minutes for every row", {
  data(timeuse, package = "timecoda", envir = environment())

  total <- timeuse$MVPA + timeuse$LPA + timeuse$SB + timeuse$SPT
  expect_true(all(abs(total - 1440) < 1e-6))
})

test_that("compositional parts are strictly positive", {
  data(timeuse, package = "timecoda", envir = environment())

  parts <- timeuse[, c("MVPA", "LPA", "SB", "SPT")]
  expect_true(all(parts > 0))
})

test_that("grouping and clustering factors have the right levels", {
  data(timeuse, package = "timecoda", envir = environment())

  expect_s3_class(timeuse$sex, "factor")
  expect_setequal(levels(timeuse$sex), c("female", "male"))

  expect_s3_class(timeuse$school, "factor")
  expect_identical(nlevels(timeuse$school), 10L)

  expect_s3_class(timeuse$region, "factor")
  expect_identical(levels(timeuse$region), c("A", "B", "C"))
})

test_that("region deterministically nests school (multi-cluster demo)", {
  data(timeuse, package = "timecoda", envir = environment())

  # Every school maps to exactly one region -> region is a coarser level above
  # school, so random = c("region", "school") is an honest nesting.
  tab <- table(timeuse$school, timeuse$region)
  expect_true(all(rowSums(tab > 0) == 1L))

  # The documented, deterministic mapping.
  idx <- as.integer(sub("school", "", as.character(timeuse$school)))
  expected <- ifelse(idx <= 3L, "A", ifelse(idx <= 6L, "B", "C"))
  expect_identical(as.character(timeuse$region), expected)
})

test_that("grade is an ordered achievement band with real signal", {
  data(timeuse, package = "timecoda", envir = environment())

  expect_s3_class(timeuse$grade, "ordered")
  expect_identical(levels(timeuse$grade), c("below", "expected", "above"))
  # Tertile cut -> every band populated (balanced ~100 each).
  expect_true(all(table(timeuse$grade) > 0))
  expect_identical(nlevels(timeuse$grade), 3L)
})

test_that("outcomes cover each model class", {
  data(timeuse, package = "timecoda", envir = environment())

  # binary
  expect_s3_class(timeuse$obesity, "factor")
  expect_setequal(levels(timeuse$obesity), c("no", "yes"))
  expect_true(all(table(timeuse$obesity) > 0))

  # continuous
  expect_type(timeuse$zbmi, "double")

  # survival
  expect_true(all(timeuse$time > 0))
  expect_true(all(timeuse$event %in% c(0L, 1L)))
  expect_true(any(timeuse$event == 1L) && any(timeuse$event == 0L))

  # ordinal (new)
  expect_s3_class(timeuse$grade, "ordered")
})

test_that("existing columns are unchanged after the grade/region extension", {
  # Adding `grade` and `region` must not perturb any previously shipped column
  # (RNG discipline: the grade noise draw is appended after every existing
  # draw, and region is deterministic). These are exact-value regression
  # guards on the pre-extension columns.
  data(timeuse, package = "timecoda", envir = environment())

  expect_identical(timeuse$id[c(1L, 300L)], c("p001", "p300"))

  # Composition parts (rounded to 0.1 min in the generator).
  expect_equal(sum(timeuse$MVPA), 13777.9, tolerance = 1e-8)
  expect_equal(sum(timeuse$LPA), 72766.9, tolerance = 1e-8)
  expect_equal(sum(timeuse$SB), 198209.9, tolerance = 1e-8)
  expect_equal(sum(timeuse$SPT), 147245.3, tolerance = 1e-8)

  # Continuous / survival outcomes.
  expect_equal(mean(timeuse$zbmi), -0.02285, tolerance = 1e-6)
  expect_equal(sum(timeuse$time), 4974.86, tolerance = 1e-6)
  expect_identical(sum(timeuse$event), 250L)

  # Binary outcome composition.
  expect_identical(sum(timeuse$obesity == "yes"), 84L)

  # Exact first-row values across every existing column (byte-level guard).
  expect_identical(as.character(timeuse$sex[1L]), "male")
  expect_identical(as.character(timeuse$school[1L]), "school02")
  expect_equal(
    unname(unlist(timeuse[1L, c("MVPA", "LPA", "SB", "SPT")])),
    c(36.1, 320, 541.1, 542.8),
    tolerance = 1e-8
  )
  expect_equal(timeuse$zbmi[1L], -1.31, tolerance = 1e-8)
})
