# Numerical benchmark of declare_composition(), comp_impute(), comp_outliers()
# and summary.compD() against compositions, zCompositions and robCompositions.
# compositions calls are namespace-qualified and pass `robust = FALSE`
# explicitly (with robCompositions loaded, their default is corrupted), and
# mean.acomp() is re-closed to the object's scale because it returns proportions.

skip_if_not_installed("compositions")
skip_if_not_installed("zCompositions")
skip_on_cran()

parts <- c("MVPA", "LPA", "SB", "SPT")
closure <- 1440
x <- as.data.frame(timeuse[parts])
cd <- declare_composition(timeuse, parts, closure = closure)
s <- summary(cd)

test_that("closure matches compositions::acomp(total =)", {
  ref <- unclass(compositions::acomp(x, total = closure))
  expect_equal(as.matrix(tc_comp_df(cd)), ref, tolerance = 1e-6,
               ignore_attr = TRUE)
})

test_that("closure strings resolve to the documented numbers", {
  expect_equal(declare_composition(timeuse, parts, closure = "mean")$closure,
               mean(rowSums(x)), tolerance = 1e-6)
  expect_equal(declare_composition(timeuse, parts, closure = "percent")$closure,
               100)
})

test_that("summary geometric mean matches compositions::mean.acomp() and a hand calculation", {
  ref <- as.numeric(unclass(compositions::mean.acomp(
    compositions::acomp(x, total = closure), robust = FALSE
  )))
  ref <- ref / sum(ref) * closure
  expect_equal(s$desc$geom_mean, ref, tolerance = 1e-6)
  hand <- exp(colMeans(log(as.matrix(tc_comp_df(cd)))))
  expect_equal(s$desc$geom_mean, unname(hand / sum(hand) * closure), tolerance = 1e-6)
})

test_that("summary variation matrix matches compositions::variation()", {
  ref <- unclass(compositions::variation(compositions::acomp(x), robust = FALSE))
  expect_equal(as.matrix(s$var_matrix[parts]), ref[parts, parts],
               tolerance = 1e-6, ignore_attr = TRUE)
})

test_that("summary descriptives and raw row totals match hand computation", {
  raw <- as.matrix(x)
  expect_equal(s$desc$mean, unname(colMeans(raw)), tolerance = 1e-6)
  expect_equal(s$desc$sd, unname(apply(raw, 2, stats::sd)), tolerance = 1e-6)
  expect_equal(s$desc$median, unname(apply(raw, 2, stats::median)),
               tolerance = 1e-6)
  rt <- rowSums(x)
  expect_equal(s$row_totals$mean, mean(rt), tolerance = 1e-6)
  expect_equal(s$row_totals$sd, stats::sd(rt), tolerance = 1e-6)
})

test_that("comp_impute() matches zCompositions called directly with the same dl", {
  z <- timeuse
  z$MVPA[1:8] <- 0
  z$SB[c(3, 20)] <- 0
  zx <- as.matrix(z[parts])
  min_pos <- apply(zx, 2, function(v) min(v[v > 0]))
  dl <- matrix(rep(min_pos, each = nrow(zx)), nrow = nrow(zx),
               dimnames = list(NULL, parts))

  lrem <- comp_impute(z, parts, method = "lrEM")
  ref <- as.matrix(zCompositions::lrEM(zx, label = 0, dl = dl,
                                        suppress.print = TRUE))
  expect_equal(as.matrix(lrem$parts_imputed), ref, tolerance = 1e-3,
               ignore_attr = TRUE)

  mr <- comp_impute(z, parts, method = "multRepl")
  ref_mr <- as.matrix(zCompositions::multRepl(zx, label = 0, dl = dl))
  expect_equal(as.matrix(mr$parts_imputed), ref_mr, tolerance = 1e-3,
               ignore_attr = TRUE)

  nz <- zx > 0
  expect_equal(as.matrix(lrem$parts_imputed)[nz], zx[nz], tolerance = 1e-6)
  expect_true(all(as.matrix(lrem$parts_imputed)[!nz] < dl[!nz]))
  expect_equal(lrem$imput_table$n_zero, as.integer(colSums(zx == 0)))
})

test_that("plot_variation() draws the compositions variation matrix", {
  ref <- unclass(compositions::variation(compositions::acomp(x), robust = FALSE))
  p <- plot_variation(cd)
  expect_equal(sort(p$data$value), sort(as.vector(ref)), tolerance = 1e-6)
})

# Outlier flagging (added 2026-09-18 with comp_outliers()). MCD re-randomizes
# its subset search per call, so every covMcd-backed call below is preceded by
# its own set.seed(); without that these checks are flaky, with it the
# invariances reproduce at machine precision.

test_that("comp_outliers() is a pass-through to robCompositions::outCoDa()", {
  set.seed(42)
  ours <- comp_outliers(timeuse, parts)
  set.seed(42)
  ref <- robCompositions::outCoDa(x)
  expect_equal(ours$distance, as.numeric(ref$mahalDist), tolerance = 1e-6)
  expect_identical(ours$is_outlier, as.logical(ref$outlierIndex))
  expect_equal(ours$limit[1], ref$limit, tolerance = 1e-6)
})

test_that("the distances reconstruct Filzmoser & Hron (2008) independently", {
  skip_if_not_installed("robustbase")
  # Robust Mahalanobis distance of the pivot ILR coordinates, built from
  # robustbase + stats directly rather than from anything robCompositions does.
  coords <- as.matrix(robCompositions::pivotCoord(x))
  set.seed(42)
  mcd <- robustbase::covMcd(coords, alpha = 0.5)
  hand <- sqrt(stats::mahalanobis(coords, center = mcd$center, cov = mcd$cov))
  set.seed(42)
  expect_equal(comp_outliers(timeuse, parts)$distance, hand, tolerance = 1e-6)
})

test_that("the cutoff is sqrt(qchisq(quantile, df = D - 1))", {
  set.seed(42)
  expect_equal(comp_outliers(timeuse, parts)$limit[1],
               sqrt(stats::qchisq(0.975, df = length(parts) - 1)),
               tolerance = 1e-6)
})

test_that("outlier distances are invariant to closure and to the pivot", {
  # MCD is affine-equivariant, so rescaling each row cannot move the distance;
  # a pivot change is an orthogonal rotation, which it is also invariant to
  #.
  set.seed(42)
  base <- comp_outliers(timeuse, parts)$distance

  x_closed <- as.data.frame(sweep(as.matrix(x), 1, rowSums(x), "/") * closure)
  set.seed(42)
  expect_equal(comp_outliers(x_closed, parts)$distance, base, tolerance = 1e-6)

  for (i in seq_along(parts)) {
    rotated <- c(parts[i], setdiff(parts, parts[i]))
    set.seed(42)
    expect_equal(comp_outliers(timeuse, rotated)$distance, base,
                 tolerance = 1e-6)
  }
})

test_that("declare_composition() and summary.compD() carry the numbers through", {
  set.seed(42)
  ref <- comp_outliers(timeuse, parts)
  set.seed(42)
  cd_flag <- declare_composition(timeuse, parts, closure = closure,
                                 detect_outliers = TRUE)
  expect_equal(cd_flag$data_raw$.outlier_distance, ref$distance,
               tolerance = 1e-6)
  expect_identical(cd_flag$data_raw$.outlier, ref$is_outlier)

  s_flag <- summary(cd_flag)
  expect_equal(s_flag$outliers$n_outlier, sum(ref$is_outlier))
  expect_equal(s_flag$outliers$limit, ref$limit[1], tolerance = 1e-6)
})
