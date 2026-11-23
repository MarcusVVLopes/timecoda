parts <- c("MVPA", "LPA", "SB", "SPT")

test_that("check_composition passes valid data and returns it invisibly", {
  expect_identical(check_composition(timeuse, parts, total = 1440), timeuse)
})

test_that("check_composition flags missing, non-positive, and wrong-total inputs", {
  expect_error(check_composition(timeuse, c(parts, "nope")), "not found")
  bad <- timeuse
  bad$MVPA[1] <- 0
  expect_error(check_composition(bad, parts), "positive")
  expect_error(check_composition(timeuse, parts, total = 1000), "sum")
  expect_error(check_composition(timeuse, "MVPA"), "at least two")
})

test_that("check_composition allow_zero permits zeroes but still rejects negatives", {
  z <- timeuse
  z$MVPA[1] <- 0
  expect_identical(check_composition(z, parts, allow_zero = TRUE), z)
  neg <- timeuse
  neg$MVPA[1] <- -1
  expect_error(check_composition(neg, parts, allow_zero = TRUE), "negative")
})

test_that("ilr_transform appends D-1 numeric coordinate columns", {
  for (d in 3:5) {
    p <- paste0("p", seq_len(d))
    sim <- as.data.frame(matrix(runif(20 * d, 1, 10), ncol = d))
    names(sim) <- p
    out <- ilr_transform(sim, p)
    coords <- grep("^ilr", names(out), value = TRUE)
    expect_length(coords, d - 1)
    expect_true(all(vapply(out[coords], is.numeric, logical(1))))
  }
})

test_that("ilr_backtransform inverts ilr_transform for any pivot", {
  for (piv in seq_along(parts)) {
    ilr <- ilr_transform(timeuse, parts, total = 1440, pivot = piv)
    comp <- ilr_backtransform(
      ilr[paste0("ilr", seq_len(length(parts) - 1))],
      parts = parts, total = 1440, pivot = piv
    )
    for (p in parts) {
      expect_equal(comp[[p]], timeuse[[p]], tolerance = 1e-8)
    }
  }
})

test_that("a pivot named by its part equals the same pivot by position", {
  coords <- paste0("ilr", seq_len(length(parts) - 1))
  for (k in seq_along(parts)) {
    by_pos <- ilr_transform(timeuse, parts, total = 1440, pivot = k)
    expect_identical(ilr_transform(timeuse, parts, total = 1440,
                                   pivot = parts[k]), by_pos)
    expect_identical(
      ilr_backtransform(by_pos[coords], parts, total = 1440, pivot = parts[k]),
      ilr_backtransform(by_pos[coords], parts, total = 1440, pivot = k)
    )
  }
  expect_error(ilr_transform(timeuse, parts, pivot = "sleep"),
               "must be one of")
  expect_error(ilr_transform(timeuse, parts, pivot = c("SB", "SPT")),
               "must be one of")
})

test_that("geom_mean_comp closes to the requested total", {
  gm <- geom_mean_comp(timeuse, parts, total = 1440)
  expect_equal(sum(gm), 1440, tolerance = 1e-8)
  expect_named(gm, parts)
})

test_that("variation_matrix is symmetric with part dimnames", {
  v <- variation_matrix(timeuse, parts)
  expect_equal(dim(v), c(4, 4))
  expect_equal(rownames(v), parts)
  expect_equal(v, t(v), tolerance = 1e-8)
})
