# Placeholder so the testthat harness has at least one test file while the
# package is still code-less. Per-function test files (one per R/ file)
# replace the need for this as functions are implemented.

test_that("package namespace loads", {
  expect_true("timecoda" %in% loadedNamespaces())
})
