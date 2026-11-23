# Numerical benchmark of the ILR and composition helpers against
# robCompositions and compositions. Their calls are namespace-qualified rather
# than attached, because attaching compositions alongside robCompositions
# corrupts later acomp()/ilr() calls; `x` is coerced to a data.frame because
# pivotCoord() mishandles tibble subsetting.

skip_if_not_installed("robCompositions")
skip_if_not_installed("compositions")
skip_on_cran()

parts <- c("MVPA", "LPA", "SB", "SPT")
x <- as.data.frame(timeuse[parts])
total <- 1440
coord_cols <- paste0("ilr", seq_len(length(parts) - 1))

test_that("ilr_transform() matches robCompositions::pivotCoord() directly", {
  ours <- ilr_transform(timeuse, parts, total = total)
  ref <- as.matrix(robCompositions::pivotCoord(x, pivotvar = 1))
  expect_equal(as.matrix(ours[coord_cols]), ref, tolerance = 1e-6, ignore_attr = TRUE)
})

test_that("ilr_transform() agrees with compositions::ilr() up to an orthonormal basis change", {
  ours <- ilr_transform(timeuse, parts, total = total)
  comp_ilr <- compositions::ilr(compositions::acomp(x))
  dist_ours <- as.matrix(stats::dist(as.matrix(ours[coord_cols])))
  dist_compositions <- as.matrix(stats::dist(unclass(comp_ilr)))
  expect_equal(dist_ours, dist_compositions, tolerance = 1e-6)
})

test_that("ilr_backtransform() matches robCompositions::pivotCoordInv() directly", {
  ours <- ilr_transform(timeuse, parts, total = total)
  back_ours <- ilr_backtransform(ours[coord_cols], parts = parts, total = total, pivot = 1)
  back_ref <- as.matrix(robCompositions::pivotCoordInv(as.matrix(ours[coord_cols])))
  back_ref <- back_ref / rowSums(back_ref) * total
  colnames(back_ref) <- parts
  expect_equal(as.matrix(back_ours[parts]), back_ref[, parts], tolerance = 1e-6, ignore_attr = TRUE)
})

test_that("geom_mean_comp() matches a hand-computed closed geometric mean", {
  gm_ours <- geom_mean_comp(timeuse, parts, total = total)
  gm_hand <- exp(colMeans(log(as.matrix(x))))
  gm_hand <- gm_hand / sum(gm_hand) * total
  expect_equal(as.numeric(gm_ours[parts]), unname(gm_hand[parts]), tolerance = 1e-6)
})

test_that("variation_matrix() matches a hand-computed var(log(x_i / x_j)) matrix", {
  v_ours <- variation_matrix(timeuse, parts)
  D <- length(parts)
  v_hand <- matrix(0, D, D, dimnames = list(parts, parts))
  for (i in seq_len(D)) {
    for (j in seq_len(D)) {
      v_hand[i, j] <- stats::var(log(x[[parts[i]]] / x[[parts[j]]]))
    }
  }
  expect_equal(v_ours, v_hand, tolerance = 1e-6, ignore_attr = TRUE)
})
