# Numerical benchmark of optimal_composition(), the Goldilocks-zone search.
# Every reference is computed outside timecoda: pivot coordinates from a
# hand-built orthonormal basis, the minute grid enumerated over the whole
# simplex in integer arithmetic, the box limits from colMeans(), sd() and
# quantile(), the ellipsoid from colMeans(), cov(), mahalanobis() and the
# chi-square quantile, predictions from stats::predict() on the backend fit
# (factor covariates averaged with equal weight, as emmeans does), percentile
# ranks from rank(), zone centres as closed geometric means, and the per-point
# predictions of predict_realloc(). The box search follows Dumuid et al.
# (2021), the percentile box Brakenridge et al. (2024) and the ellipsoid
# Mellow et al. (2026). The quadratic outcome is simulated with its maximum
# inside the region; a linear fit's best candidate lies in the ellipsoid's
# outer shell.

parts <- c("MVPA", "LPA", "SB", "SPT")
D <- length(parts)
coords <- paste0("ilr", seq_len(D - 1))
total <- 1440

# Column j: part j against the geometric mean of parts j+1..D.
pivot_basis <- function(D) {
  V <- matrix(0, D, D - 1)
  for (j in seq_len(D - 1)) {
    V[j, j] <- sqrt((D - j) / (D - j + 1))
    V[(j + 1):D, j] <- -1 / sqrt((D - j) * (D - j + 1))
  }
  V
}
V <- pivot_basis(D)
to_ilr <- function(X, basis = V) {
  Z <- log(as.matrix(X)) %*% basis
  colnames(Z) <- paste0("ilr", seq_len(ncol(basis)))
  Z
}
from_ilr <- function(Z) {
  E <- exp(as.matrix(Z) %*% t(V))
  X <- E / rowSums(E) * total
  colnames(X) <- parts
  X
}
hand_gm <- function(X) {
  g <- exp(colMeans(log(as.matrix(X))))
  g / sum(g) * total
}
pct_rank <- function(x) rank(x, ties.method = "average") / length(x)
max_diff <- function(a, b) {
  a <- as.numeric(unlist(a))
  b <- as.numeric(unlist(b))
  if (length(a) != length(b) || !identical(is.na(a), is.na(b))) return(Inf)
  ok <- !is.na(a)
  max(c(0, abs(a[ok] - b[ok])))
}
n_xor <- function(a, b) {
  if (length(a) != length(b) || anyNA(a) || anyNA(b)) return(Inf)
  sum(xor(as.logical(a), as.logical(b)))
}
centre_of <- function(oc, k = "overall") {
  unlist(oc$zones[oc$zones$outcome == k, paste0("centre_", parts)])
}

# Every m-tuple of positive integers summing to K.
int_comps <- function(K, m) {
  if (m == 2) return(cbind(seq_len(K - 1), rev(seq_len(K - 1))))
  do.call(rbind, lapply(seq_len(K - m + 1), function(a) {
    cbind(a, int_comps(K - a, m - 1))
  }))
}
# Every composition of `total` into `d` parts that are positive multiples of
# `spacing`, the whole simplex, then filtered by `keep`.
simplex_grid <- function(spacing, d = D, keep = NULL) {
  M <- int_comps(round(total / spacing), d) * spacing
  dimnames(M) <- NULL
  if (!is.null(keep)) M <- M[keep(M), , drop = FALSE]
  M
}
# Each part within mean +/- sdk SD, or within the q quantiles, never beyond
# its observed range.
box_limits <- function(X, sdk = NULL, q = NULL) {
  if (is.null(q)) {
    lo <- colMeans(X) - sdk * apply(X, 2, stats::sd)
    hi <- colMeans(X) + sdk * apply(X, 2, stats::sd)
  } else {
    lo <- apply(X, 2, stats::quantile, probs = q[1], names = FALSE)
    hi <- apply(X, 2, stats::quantile, probs = q[2], names = FALSE)
  }
  lo <- pmax(lo, apply(X, 2, min))
  hi <- pmin(hi, apply(X, 2, max))
  function(M) rowSums(sweep(M, 2, lo, ">=") & sweep(M, 2, hi, "<=")) == ncol(M)
}
# Candidate sets compared as sets of rows (0 when equal), and the package's
# rows reordered to a hand matrix.
row_key <- function(M, scale = 1) {
  do.call(paste, c(as.data.frame(round(as.matrix(M) * scale, 6)), sep = "|"))
}
set_diff <- function(M, cand, scale = 1) {
  a <- row_key(M)
  b <- row_key(cand, scale)
  sum(!a %in% b) + sum(!b %in% a) + sum(duplicated(b))
}
align <- function(oc, M, cols = parts, scale = 1) {
  oc$candidates[match(row_key(M), row_key(oc$candidates[cols], scale)), ]
}
seeded_pick <- function(n) {
  set.seed(20261112)
  sort(sample(n, min(20, n)))
}

X0 <- as.matrix(timeuse[parts])
X0 <- X0 / rowSums(X0) * total
Z0 <- to_ilr(X0)
mu0 <- colMeans(Z0)
S0 <- stats::cov(Z0)
r2_0 <- stats::qchisq(0.80, df = D - 1)

# y = 10 - (z - z*)' A (z - z*) + sex effect + noise, z* inside the region.
zstar <- mu0 + 0.3 * sqrt(diag(S0))
rot <- function(a, i, j) {
  R <- diag(D - 1)
  R[i, i] <- cos(a)
  R[j, j] <- cos(a)
  R[i, j] <- -sin(a)
  R[j, i] <- sin(a)
  R
}
Qr <- rot(0.7, 1, 2) %*% rot(0.4, 2, 3) %*% rot(-0.5, 1, 3)
A <- Qr %*% diag(c(2, 1.5, 1)) %*% t(Qr)
set.seed(15)
dz <- sweep(Z0, 2, zstar)
tu <- timeuse
tu$y <- 10 - rowSums((dz %*% A) * dz) + ifelse(tu$sex == "male", 0.4, 0) +
  stats::rnorm(nrow(tu), 0, 0.5)
cd <- declare_composition(tu, parts, closure = total)

# Link-scale prediction at ILR nodes, factor covariates averaged with equal
# weight unless pinned, numeric ones at their mean.
hand_link <- function(fit, Z, pins = list()) {
  d <- fit$data
  covs <- intersect(setdiff(all.vars(stats::formula(fit$model)[[3]]), colnames(Z)), names(d))
  grid <- list()
  for (v in covs) {
    grid[[v]] <- if (!is.null(pins[[v]])) pins[[v]] else
      if (is.factor(d[[v]])) levels(d[[v]]) else mean(d[[v]])
  }
  combos <- if (length(grid)) expand.grid(grid, stringsAsFactors = FALSE) else data.frame(.k = 1)
  preds <- lapply(seq_len(nrow(combos)), function(i) {
    nd <- as.data.frame(Z)
    for (v in names(grid)) {
      nd[[v]] <- if (is.factor(d[[v]])) factor(combos[[v]][i], levels = levels(d[[v]])) else combos[[v]][i]
    }
    if (inherits(fit$model, "glm")) stats::predict(fit$model, newdata = nd, type = "link")
    else stats::predict(fit$model, newdata = nd)
  })
  Reduce(`+`, preds) / length(preds)
}
# The ellipsoid's candidates: every composition on multiples of `spacing`
# (the last part the remainder) inside the chi-square fence.
hand_grid <- function(spacing) {
  vals <- seq(spacing, total - (D - 1) * spacing, by = spacing)
  gr <- as.matrix(expand.grid(rep(list(vals), D - 1)))
  gr <- cbind(gr, total - rowSums(gr))
  gr <- gr[gr[, D] >= spacing, , drop = FALSE]
  nodes <- to_ilr(gr)
  d2 <- stats::mahalanobis(nodes, mu0, S0)
  list(Z = nodes[d2 <= r2_0, , drop = FALSE], d2 = d2[d2 <= r2_0])
}
# Candidates reordered to match the rows of a hand node matrix.
matched <- function(oc, M, cols = coords) {
  key <- function(X) do.call(order, as.data.frame(round(X, 8)))
  idx <- integer(nrow(M))
  idx[key(M)] <- key(as.matrix(oc$candidates[cols]))
  oc$candidates[idx, ]
}
count_messages <- function(expr, pattern) {
  n <- 0
  val <- withCallingHandlers(expr, message = function(m) {
    if (grepl(pattern, conditionMessage(m))) n <<- n + 1
    invokeRestart("muffleMessage")
  })
  list(value = val, n = n)
}
# A search against a hand candidate set: the same candidates, the relative
# value of each against the fit's reference day, the best 5% by rank, the
# closed geometric mean and the per-part range of the zone.
expect_hand_search <- function(oc, fit, M, direction) {
  k <- oc$outcomes[1]
  expect_equal(set_diff(M, oc$candidates[parts]), 0)
  cc <- align(oc, M)
  e <- hand_link(fit, to_ilr(M))
  e_from <- hand_link(fit, to_ilr(t(unlist(fit$ref_comp)[parts])))
  expect_lt(max_diff(cc[[paste0("estimate_", k)]], e - e_from), 1e-8)
  in_zone <- pct_rank(if (direction == "max") e else -e) >= 0.95
  expect_equal(n_xor(cc$in_zone, in_zone), 0)
  Xz <- M[in_zone, , drop = FALSE]
  ov <- oc$zones[oc$zones$outcome == "overall", ]
  expect_equal(ov$n, sum(in_zone))
  expect_lt(max_diff(centre_of(oc), hand_gm(Xz)), 1e-8)
  expect_equal(max_diff(ov[paste0("min_", parts)], apply(Xz, 2, min)), 0)
  expect_equal(max_diff(ov[paste0("max_", parts)], apply(Xz, 2, max)), 0)
}
quad_formula <- function(lhs) {
  stats::as.formula(paste(lhs, "~ ilr()^2 + I(ilr()^2) + sex"))
}

test_that("the ellipsoid grid, predictions, zone and centre equal the hand computation", {
  skip_on_cran()
  f_quad <- comp_iv(cd, formula = y ~ ilr()^2 + I(ilr()^2) + sex, model = "lm")
  oc <- suppressMessages(optimal_composition(f_quad, region = "ellipsoid",
                                             spacing = 30, type = "absolute"))
  par <- oc$region$parameters[[1]]
  expect_lt(max_diff(Z0, f_quad$data[coords]), 1e-12)
  expect_lt(max_diff(par$mu, mu0), 1e-8)
  expect_lt(max_diff(par$S, S0), 1e-8)
  expect_lt(abs(par$r2 - r2_0), 1e-8)

  lat <- hand_grid(30)
  expect_equal(nrow(oc$candidates), nrow(lat$Z))
  cc <- matched(oc, lat$Z)
  expect_lt(max_diff(cc[coords], lat$Z), 1e-8)
  expect_lt(max_diff(cc[parts], from_ilr(lat$Z)), 1e-8)
  expect_lt(max_diff(cc$mahalanobis, lat$d2), 1e-8)

  e <- hand_link(f_quad, lat$Z)
  expect_lt(max_diff(cc$estimate_y, e), 1e-8)
  node <- stats::setNames(as.numeric(cc[1, parts]), parts)
  pr <- predict_realloc(f_quad, delta = 1, from = "SB", to = "SPT",
                        type = "absolute", reference = node)
  expect_lt(abs(pr$estimate[1] - e[1]), 1e-8)

  rk <- pct_rank(e)
  expect_lt(max_diff(cc$rank_y, rk), 1e-12)
  in_zone <- rk >= 0.95
  expect_equal(n_xor(cc$in_zone, in_zone), 0)
  Xz <- from_ilr(lat$Z[in_zone, , drop = FALSE])
  cen <- hand_gm(Xz)
  expect_lt(max_diff(cen, from_ilr(t(colMeans(lat$Z[in_zone, , drop = FALSE])))), 1e-8)
  ov <- oc$zones[oc$zones$outcome == "overall", ]
  expect_equal(ov$n, sum(in_zone))
  expect_lt(max_diff(centre_of(oc), cen), 1e-8)
  expect_lt(max_diff(ov[paste0("min_", parts)], apply(Xz, 2, min)), 1e-8)
  expect_lt(max_diff(ov[paste0("max_", parts)], apply(Xz, 2, max)), 1e-8)
  expect_lt(max_diff(oc$optimum[1, parts], cen), 1e-8)
  expect_lt(abs(oc$optimum$estimate[1] - hand_link(f_quad, to_ilr(t(cen)))), 1e-8)

  # The fitted stationary point -H^-1 g lies inside the region.
  b <- stats::coef(f_quad$model)
  H <- diag(2 * b[paste0("I(", coords, "^2)")], D - 1)
  for (i in 1:(D - 2)) for (j in (i + 1):(D - 1)) {
    H[i, j] <- H[j, i] <- b[[paste0(coords[i], ":", coords[j])]]
  }
  z_hat <- -solve(H, b[coords])
  expect_lt(stats::mahalanobis(z_hat, mu0, S0), r2_0)
})

test_that("a linear fit's best candidate lies in the ellipsoid's outer shell", {
  skip_on_cran()
  f_lin <- comp_iv(cd, "zbmi", model = "lm", covariates = "sex")
  for (direction in c("max", "min")) {
    run <- count_messages(optimal_composition(f_lin, direction = direction,
                                              region = "ellipsoid", spacing = 10,
                                              type = "absolute"),
                          "linear in the composition")
    oc <- run$value
    expect_equal(run$n, 1)
    lat <- hand_grid(10)
    cc <- matched(oc, lat$Z)
    s <- if (direction == "max") 1 else -1
    e <- cc$estimate_zbmi
    z_best <- lat$Z[which.max(s * e), ]
    expect_gt(stats::mahalanobis(z_best, mu0, S0) / r2_0, 0.9)
    in_zone <- pct_rank(s * e) >= 0.95
    share <- mean(lat$d2[in_zone] > 0.9 * r2_0)
    ob <- oc$zones$on_boundary[oc$zones$outcome == "overall"]
    expect_lt(abs(ob - share), 1e-12)
  }
})

test_that("the link-scale zone of a binomial glm equals a ranking on the response scale", {
  skip_on_cran()
  f_glmp <- comp_iv(cd, formula = obesity ~ poly(ilr(), 2) + sex, model = "glm",
                    family = "binomial")
  oc <- optimal_composition(f_glmp, region = "ellipsoid", spacing = 30)
  lat <- hand_grid(30)
  cc <- matched(oc, lat$Z)
  rk <- pct_rank(stats::plogis(hand_link(f_glmp, lat$Z)))
  expect_lt(max_diff(cc$rank_obesity, rk), 1e-12)
  expect_equal(n_xor(cc$in_zone, rk >= 0.95), 0)

  oc_f <- optimal_composition(f_glmp, region = "ellipsoid", spacing = 30,
                              at = list(sex = "female"))
  cc_f <- matched(oc_f, lat$Z)
  nd <- as.data.frame(lat$Z)
  nd$sex <- factor("female", levels = levels(f_glmp$data$sex))
  rk_f <- pct_rank(stats::predict(f_glmp$model, newdata = nd, type = "response"))
  expect_lt(max_diff(cc_f$rank_obesity, rk_f), 1e-12)
  expect_equal(n_xor(cc_f$in_zone, rk_f >= 0.95), 0)
})

test_that("two outcomes intersect, q_overlap is the best common rank, min flips the ranks", {
  skip_on_cran()
  f_quad <- comp_iv(cd, formula = y ~ ilr()^2 + I(ilr()^2) + sex, model = "lm")
  f_glmp <- comp_iv(cd, formula = obesity ~ poly(ilr(), 2) + sex, model = "glm",
                    family = "binomial")
  fits <- list(y = f_quad, obesity = f_glmp)
  oc <- suppressMessages(optimal_composition(fits, region = "ellipsoid",
                                             spacing = 30))
  lat <- hand_grid(30)
  cc <- matched(oc, lat$Z)
  ry <- pct_rank(hand_link(f_quad, lat$Z))
  ro <- pct_rank(hand_link(f_glmp, lat$Z))
  expect_lt(max_diff(cc$rank_y, ry), 1e-12)
  expect_lt(max_diff(cc$rank_obesity, ro), 1e-12)
  expect_equal(n_xor(cc$in_zone, ry >= 0.95 & ro >= 0.95), 0)

  expect_false(any(ry >= 0.99 & ro >= 0.99))
  oc_e <- suppressMessages(optimal_composition(fits, zone = 0.99,
                                               region = "ellipsoid", spacing = 30))
  ov <- oc_e$zones[oc_e$zones$outcome == "overall", ]
  expect_equal(ov$n, 0)
  expect_false(any(oc_e$candidates$in_zone))
  expect_lt(abs(ov$q_overlap - max(pmin(ry, ro))), 1e-12)

  oc_d <- suppressMessages(optimal_composition(fits, direction = c("max", "min"),
                                               region = "ellipsoid", spacing = 30))
  cc_d <- matched(oc_d, lat$Z)
  ro_min <- pct_rank(-hand_link(f_glmp, lat$Z))
  expect_lt(max_diff(cc_d$rank_obesity, ro_min), 1e-12)
  expect_lt(max_diff(cc_d$rank_y, ry), 1e-12)
  expect_equal(n_xor(cc_d$in_zone, ry >= 0.95 & ro_min >= 0.95), 0)
})

test_that("a threshold band rules the zone by hand band membership", {
  skip_on_cran()
  f_quad <- comp_iv(cd, formula = y ~ ilr()^2 + I(ilr()^2) + sex, model = "lm")
  lat <- hand_grid(30)
  e <- hand_link(f_quad, lat$Z)
  band <- unname(stats::quantile(e, c(0.60, 0.90)))
  oc <- optimal_composition(f_quad, region = "ellipsoid", spacing = 30,
                            type = "absolute", threshold = list(y = band))
  cc <- matched(oc, lat$Z)
  in_band <- e >= band[1] & e <= band[2]
  expect_equal(n_xor(cc$in_zone, in_band), 0)
  expect_lt(max_diff(centre_of(oc), hand_gm(from_ilr(lat$Z[in_band, ]))), 1e-8)
})

test_that("the default search equals the box search of Dumuid et al. (2021) by hand", {
  skip_on_cran()
  in_box <- box_limits(X0, sdk = 3)
  quad <- comp_iv(cd, formula = quad_formula("zbmi"), model = "lm")
  fits <- list(
    comp_iv(cd, "zbmi", model = "lm", covariates = "sex"),
    quad,
    comp_iv(cd, "obesity", model = "glm", family = "binomial", covariates = "sex")
  )
  M30 <- simplex_grid(30, keep = in_box)
  for (fit in fits) {
    oc <- suppressMessages(optimal_composition(fit, direction = "min", spacing = 30))
    expect_hand_search(oc, fit, M30, "min")
  }
  # The default call: the 10-minute grid of the box of the mean +/- 3 SD.
  oc <- optimal_composition(quad, direction = "min")
  expect_identical(c(oc$total, oc$spacing, oc$sd), c(1440, 10, 3))
  expect_hand_search(oc, quad, simplex_grid(10, keep = in_box), "min")
})

test_that("the ellipsoid equals the chi-square fence of Mellow et al. (2026) in either basis", {
  skip_on_cran()
  inside <- function(M) stats::mahalanobis(to_ilr(M), mu0, S0) <= r2_0
  M30 <- simplex_grid(30, keep = inside)
  quad <- comp_iv(cd, formula = quad_formula("zbmi"), model = "lm")
  oc <- optimal_composition(quad, direction = "min", region = "ellipsoid",
                            spacing = 30)
  expect_identical(oc$region$parameters[[1]]$r2, r2_0)
  expect_hand_search(oc, quad, M30, "min")
  # The sequential-binary-partition basis of the authors' app gives the same
  # fence: the Mahalanobis distance does not depend on the ILR basis.
  sbp <- matrix(0, D, D - 1)
  for (j in seq_len(D - 1)) {
    sbp[j, j] <- 1
    sbp[(j + 1):D, j] <- -1
  }
  psi <- compositions::gsi.buildilrBase(sbp)
  Zs <- to_ilr(X0, psi)
  all30 <- simplex_grid(30)
  d_pivot <- stats::mahalanobis(to_ilr(all30), mu0, S0)
  d_sbp <- stats::mahalanobis(to_ilr(all30, psi), colMeans(Zs), stats::cov(Zs))
  expect_lt(max(abs(d_pivot - d_sbp)), 1e-8)
  expect_equal(n_xor(d_pivot <= r2_0, d_sbp <= r2_0), 0)
})

test_that("the percentile box of Brakenridge et al. (2024) equals a hand enumeration", {
  skip_on_cran()
  quad <- comp_iv(cd, formula = quad_formula("zbmi"), model = "lm")
  oc <- optimal_composition(quad, direction = "min", coverage = c(0.01, 0.99))
  expect_null(oc$sd)
  in_box <- box_limits(X0, q = c(0.01, 0.99))
  expect_hand_search(oc, quad, simplex_grid(10, keep = in_box), "min")
})

test_that("the grid, the zone and the reference day follow the declared closure", {
  skip_on_cran()
  search <- function(closure) {
    cdx <- declare_composition(timeuse, parts, closure = closure)
    fit <- comp_iv(cdx, formula = quad_formula("zbmi"), model = "lm")
    list(fit = fit, oc = optimal_composition(fit, direction = "min"))
  }
  base <- search(1440)$oc
  M <- as.matrix(base$candidates[parts])
  for (closure in c(24, 1)) {
    run <- search(closure)
    oc <- run$oc
    scale <- 1440 / closure
    expect_identical(oc$total, closure)
    expect_lt(abs(oc$spacing - closure / 144), 1e-15 * closure)
    expect_equal(nrow(oc$candidates), nrow(M))
    expect_equal(set_diff(M, oc$candidates[parts], scale), 0)
    cc <- align(oc, M, scale = scale)
    expect_lt(max_diff(cc[parts], M / scale), 1e-9 * closure)
    expect_lt(max_diff(cc$estimate_zbmi, base$candidates$estimate_zbmi), 1e-8)
    expect_equal(n_xor(cc$in_zone, base$candidates$in_zone), 0)
    expect_lt(max_diff(centre_of(oc), centre_of(base) / scale), 1e-8 * closure)
    from <- unlist(oc$optimum[1, paste0("from_", parts)])
    expect_lt(max_diff(from, unlist(run$fit$ref_comp)[parts]), 1e-12 * closure)
    gm <- geom_mean_comp(as.data.frame(X0 / scale), parts, total = closure)
    expect_lt(max_diff(from, unlist(gm)[parts]), 1e-9 * closure)
  }
})

test_that("the ellipsoid's part limits keep every candidate of the whole simplex", {
  skip_on_cran()
  tu3 <- timeuse
  tu3$PA <- tu3$MVPA + tu3$LPA
  p3 <- c("PA", "SB", "SPT")
  X3 <- as.matrix(tu3[p3])
  cases <- list(
    list(parts = parts, X = X0,
         fit = comp_iv(cd, "zbmi", model = "lm", covariates = "sex")),
    list(parts = p3, X = X3 / rowSums(X3) * total,
         fit = comp_iv(declare_composition(tu3, p3, closure = total), "zbmi",
                       model = "lm", covariates = "sex"))
  )
  for (case in cases) {
    basis <- pivot_basis(length(case$parts))
    Z <- to_ilr(case$X, basis)
    r2 <- stats::qchisq(0.80, df = length(case$parts) - 1)
    M <- simplex_grid(30, length(case$parts), keep = function(M) {
      stats::mahalanobis(to_ilr(M, basis), colMeans(Z), stats::cov(Z)) <= r2
    })
    oc <- suppressMessages(optimal_composition(case$fit, region = "ellipsoid",
                                               spacing = 30))
    expect_equal(set_diff(M, oc$candidates[case$parts]), 0)
  }
})

# On the default relative scale a candidate's estimate is its link-scale
# contrast with `from`, the fit's reference composition. The per-point route
# for the same contrast is the custom reallocation from `from` to the
# candidate, anchored at `from`; for an absolute search it is the reference
# row of a one-minute move anchored at the candidate.
per_point_relative <- function(fit, comps, category = NULL, at = NULL) {
  from <- unlist(fit$ref_comp)[parts]
  change <- as.data.frame(sweep(as.matrix(comps[parts]), 2, from))
  pr <- predict_realloc(fit, realloc = change, reference = from,
                        type = "relative", at = at)
  if (!is.null(category)) pr <- pr[pr$response_level == category, ]
  pr
}
per_point_absolute <- function(fit, comps, at = NULL) {
  linkfun <- if (inherits(fit$family, "family")) fit$family$linkfun else identity
  vapply(seq_len(nrow(comps)), function(i) {
    node <- stats::setNames(as.numeric(comps[i, parts]), parts)
    r <- predict_realloc(fit, delta = 1, from = "SB", to = "SPT", type = "absolute",
                         reference = node, at = at)
    linkfun(r$estimate[1])
  }, numeric(1))
}
# Twenty seeded candidates of a search, on both scales where the engine has
# an absolute response; the engines without one refuse it. The 30-minute grid
# keeps the file fast.
expect_one_grid <- function(fit, category = NULL, spacing = 30, at = NULL) {
  oc <- suppressMessages(optimal_composition(fit, spacing = spacing,
                                             category = category, at = at))
  k <- oc$outcomes[1]
  pick <- seeded_pick(nrow(oc$candidates))
  est <- oc$candidates[[paste0("estimate_", k)]][pick]
  if (oc$measures[[k]] != "Mean difference") est <- exp(est)
  pr <- per_point_relative(fit, oc$candidates[pick, ], category, at)
  expect_lt(max_diff(est, pr$estimate), 1e-6)
  if (fit$model_type %in% c("cox", "coxme", "clm", "clmm", "mblogit")) {
    expect_error(optimal_composition(fit, spacing = spacing, category = category,
                                     type = "absolute"), "relative")
    return(invisible(list(oc = oc)))
  }
  oca <- suppressMessages(optimal_composition(fit, spacing = spacing, at = at,
                                              type = "absolute"))
  expect_identical(oca$candidates[parts], oc$candidates[parts])
  pp <- per_point_absolute(fit, oca$candidates[pick, ], at)
  expect_lt(max_diff(oca$candidates[[paste0("estimate_", k)]][pick], pp), 1e-10)
  invisible(list(oc = oc, oca = oca, pp = pp, pick = pick))
}

test_that("candidate values equal the per-point route on the base engines", {
  skip_on_cran()
  set.seed(2026)
  tv <- tu
  tv[parts] <- tv[parts] * stats::runif(nrow(tv), 0.8, 1.2)
  cd_tv <- declare_composition(tv, parts, closure = "mean")
  expect_one_grid(comp_iv(cd, formula = quad_formula("zbmi"), model = "lm"))
  expect_one_grid(comp_iv(cd, formula = quad_formula("obesity"), model = "glm",
                          family = "binomial"))
  expect_one_grid(comp_iv(cd, "zbmi", model = "lmer", random = "school",
                          covariates = "sex"))
  expect_one_grid(comp_iv(cd, "obesity", model = "glmer", family = "binomial",
                          random = "school", covariates = "sex"))
  expect_one_grid(comp_iv(cd, formula = quad_formula("survival::Surv(time, event)"),
                          model = "cox"))
  expect_one_grid(comp_iv(cd_tv, "zbmi", model = "lm", include_total = TRUE))
  expect_one_grid(comp_iv(cd, formula = zbmi ~ splines::ns(ilr(), 3) + sex,
                          model = "lm"))
})

test_that("candidate values of a robust linear fit equal the per-point route", {
  skip_on_cran()
  skip_if_not_installed("robustbase")
  set.seed(1)
  expect_one_grid(comp_iv(cd, formula = quad_formula("zbmi"), model = "lmrob"))
})

test_that("candidate values of a robust mixed fit equal the per-point route", {
  skip_on_cran()
  skip_if_not_installed("robustlmm")
  fit <- comp_iv(cd, "zbmi", model = "rlmer", random = "school", covariates = "sex")
  run <- expect_one_grid(fit)
  # The interval at the centre is asymptotic, as on the per-point route.
  optimum <- run$oca$optimum
  centre <- stats::setNames(as.numeric(optimum[1, parts]), parts)
  r <- predict_realloc(fit, delta = 1, from = "SB", to = "SPT",
                       type = "absolute", reference = centre)
  expect_equal(r$df[1], Inf)
  half <- stats::qnorm(0.975) * r$std.error[1]
  expect_lt(abs(optimum$estimate[1] - r$estimate[1]), 1e-8)
  expect_lt(abs(optimum$conf.low[1] - (r$estimate[1] - half)), 1e-8)
  expect_lt(abs(optimum$conf.high[1] - (r$estimate[1] + half)), 1e-8)
})

test_that("candidate values of a mixed-effects Cox fit equal the per-point route", {
  skip_on_cran()
  skip_if_not_installed("coxme")
  expect_one_grid(comp_iv(cd, time = "time", event = "event", model = "coxme",
                          random = "school", covariates = "sex"))
})

test_that("candidate values of ordinal fits equal the per-point route", {
  skip_on_cran()
  skip_if_not_installed("ordinal")
  expect_one_grid(comp_iv(cd, formula = quad_formula("grade"), model = "clm"))
  expect_one_grid(comp_iv(cd, "grade", model = "clmm", random = "school",
                          covariates = "sex"))
})

test_that("candidate values of a multinomial fit equal the per-point route", {
  skip_on_cran()
  skip_if_not_installed("mclogit")
  d <- timeuse
  d$grade_unordered <- factor(as.character(d$grade), levels = levels(d$grade))
  cd_m <- declare_composition(d, parts, closure = total)
  fit <- suppressWarnings(comp_iv(cd_m, "grade_unordered", model = "mblogit",
                                  covariates = "sex"))
  levs <- unique(sub("~.*$", "", names(stats::coef(fit$model))))
  expect_one_grid(fit, category = levs[1])
})

test_that("a factor interaction and a covariate held at two values stay on the one grid", {
  skip_on_cran()
  crossed <- comp_iv(cd, formula = zbmi ~ ilr()^2 + I(ilr()^2) + sex * region,
                     model = "lm")
  expect_one_grid(crossed)
  pins <- list(sex = c("female", "male"))
  quad <- comp_iv(cd, formula = quad_formula("zbmi"), model = "lm")
  run <- expect_one_grid(quad, at = pins)
  given <- optimal_composition(quad, region = as.data.frame(run$oca$candidates[run$pick, parts]),
                               zone = 0.5, type = "absolute", at = pins)
  expect_lt(max_diff(given$candidates$estimate_zbmi, run$pp), 1e-10)
})

test_that("ranks, zones and centres are the same on both scales", {
  skip_on_cran()
  fits <- list(
    zbmi = comp_iv(cd, formula = quad_formula("zbmi"), model = "lm"),
    obesity = comp_iv(cd, formula = quad_formula("obesity"), model = "glm",
                      family = "binomial")
  )
  same <- function(rel, abso) {
    cr <- rel$candidates
    ca <- abso$candidates
    cols <- grep("^(rank_|in_zone)", names(cr), value = TRUE)
    expect_identical(cr[c(parts, coords, cols)], ca[c(parts, coords, cols)])
    cols <- grep("^(centre|min|max|best)_|^(n|on_boundary|q_overlap)$",
                 names(rel$zones), value = TRUE)
    expect_identical(rel$zones[cols], abso$zones[cols])
    gain <- c("gain", "gain.low", "gain.high", "p.value")
    expect_identical(rel$optimum[c(parts, gain)], abso$optimum[c(parts, gain)])
    expect_identical(unname(as.list(rel$optimum[c("estimate", "conf.low", "conf.high")])),
                     unname(as.list(rel$optimum[c("gain", "gain.low", "gain.high")])))
  }
  for (fit in fits) {
    same(suppressMessages(optimal_composition(fit, spacing = 30)),
         suppressMessages(optimal_composition(fit, spacing = 30, type = "absolute")))
  }
  same(suppressMessages(optimal_composition(fits, spacing = 30, direction = "min")),
       suppressMessages(optimal_composition(fits, spacing = 30, direction = "min",
                                            type = "absolute")))
})

test_that("the relative optimum equals the per-point gain from the reference day", {
  skip_on_cran()
  fits <- list(
    comp_iv(cd, formula = quad_formula("zbmi"), model = "lm"),
    comp_iv(cd, formula = quad_formula("obesity"), model = "glm", family = "binomial"),
    comp_iv(cd, formula = quad_formula("survival::Surv(time, event)"), model = "cox")
  )
  for (fit in fits) {
    op <- suppressMessages(optimal_composition(fit, spacing = 60))$optimum
    from <- stats::setNames(unlist(op[1, paste0("from_", parts)]), parts)
    pr <- predict_realloc(fit, realloc = unlist(op[1, parts]) - from,
                          reference = from, type = "relative")
    expect_lt(max_diff(op[1, c("estimate", "conf.low", "conf.high")],
                       pr[1, c("estimate", "conf.low", "conf.high")]), 1e-6)
    expect_lt(abs(op$p.value[1] - pr$p.value[1]), 1e-6)
  }
})

test_that("a Cox fit's relative values are hazard ratios from its coefficients", {
  skip_on_cran()
  fit <- comp_iv(cd, formula = quad_formula("survival::Surv(time, event)"),
                 model = "cox")
  oc <- suppressMessages(optimal_composition(fit, spacing = 60))
  beta <- stats::coef(fit$model)
  pick <- unique(round(seq(1, nrow(oc$candidates), length.out = 5)))
  from <- unlist(fit$ref_comp)[parts]
  X <- rbind(as.data.frame(oc$candidates[pick, parts]), as.data.frame(t(from)))
  z <- as.data.frame(ilr_transform(X, parts, pivot = fit$pivot))[coords]
  # Each coefficient's term evaluated from its own name; the sex term is the
  # same in both rows and cancels, and a Cox model has no intercept.
  terms <- lapply(gsub("(?<!:):(?!:)", "*", names(beta), perl = TRUE), str2lang)
  eta <- vapply(seq_len(nrow(z)), function(i) {
    x <- vapply(terms, function(term) {
      if (!any(all.vars(term) %in% coords)) return(0)
      as.numeric(eval(term, as.list(z[i, ])))
    }, numeric(1))
    sum(x * beta)
  }, numeric(1))
  hand <- exp(eta[seq_along(pick)] - eta[length(eta)])
  est <- oc$candidates[[paste0("estimate_", oc$outcomes)]][pick]
  expect_lt(max_diff(exp(est), hand), 1e-8)
})

test_that("a relative threshold is a ratio versus the reference day", {
  skip_on_cran()
  fit <- comp_iv(cd, formula = quad_formula("obesity"), model = "glm",
                 family = "binomial")
  oc <- optimal_composition(fit, threshold = list(obesity = 0.8),
                            direction = "min")
  member <- exp(oc$candidates$estimate_obesity) <= 0.8
  expect_gt(sum(member), 0)
  expect_identical(oc$candidates$in_zone_obesity, member)
  expect_identical(oc$candidates$in_zone, member)
  expect_identical(oc$zones$threshold_high[1], 0.8)
})

# The surface of plot() at 20 seeded cells of its first panel equals
# predict_realloc() for the same moves from the same day.
expect_plot_surface <- function(opt, category = NULL) {
  p <- suppressMessages(quiet_draw(plot(opt, resolution = 30)))
  expect_length(p, choose(D, 3))
  cells <- p[[1]]$data
  tri <- c(cells$xpart[1], cells$ypart[1], cells$zpart[1])
  idx <- seeded_pick(nrow(cells))
  comps <- matrix(centre_of(opt), nrow = length(idx), ncol = D, byrow = TRUE,
                  dimnames = list(NULL, parts))
  comps[, tri] <- as.matrix(cells[idx, c("x", "y", "z")])
  from <- stats::setNames(unlist(opt$optimum[1, paste0("from_", parts)]), parts)
  pr <- predict_realloc(opt$fits[[1]], realloc = as.data.frame(sweep(comps, 2, from)),
                        reference = from, type = "relative")
  if (!is.null(category)) pr <- pr[pr$response_level == category, ]
  expect_lt(max_diff(cells$.pred[idx], pr$estimate), 1e-8)
}

test_that("plot() of an ordinal search draws the per-point surface", {
  skip_on_cran()
  skip_if_not_installed("ordinal")
  fit <- comp_iv(cd, formula = quad_formula("grade"), model = "clm")
  expect_plot_surface(suppressMessages(optimal_composition(fit, spacing = 30)))
})

test_that("plot() of a multinomial search draws the per-point surface", {
  skip_on_cran()
  skip_if_not_installed("mclogit")
  d <- timeuse
  d$grade_unordered <- factor(as.character(d$grade), levels = levels(d$grade))
  cd_m <- declare_composition(d, parts, closure = total)
  fit <- suppressWarnings(comp_iv(cd_m, "grade_unordered", model = "mblogit",
                                  covariates = "sex"))
  levs <- unique(sub("~.*$", "", names(stats::coef(fit$model))))
  opt <- suppressMessages(optimal_composition(fit, spacing = 30, category = levs[1]))
  expect_plot_surface(opt, category = levs[1])
})

# Three outcomes on one 30-minute grid: y_opp, simulated lowest where z-BMI is
# highest (little MVPA, much SB), conflicts with z-BMI, so their best 5% share
# no candidate; obesity is ruled by a threshold, an odds ratio of at most 0.8
# versus the reference day. Dumuid et al. (2021) and Brakenridge et al. (2024)
# widened such zones by hand in 5-point steps until they met. The searches are
# run once, on first use, and shared by the tests below.
widen_searches <- local({
  runs <- NULL
  function() {
    if (is.null(runs)) {
      tw <- timeuse
      lr <- log(tw$SB / tw$MVPA)
      set.seed(20261118)
      tw$y_opp <- -0.6 * lr + 0.08 * lr^2 + stats::rnorm(nrow(tw), sd = 0.3)
      cdw <- declare_composition(tw, parts, closure = total)
      fits <- list(
        zbmi = comp_iv(cdw, formula = quad_formula("zbmi"), model = "lm"),
        y_opp = comp_iv(cdw, formula = quad_formula("y_opp"), model = "lm"),
        obesity = comp_iv(cdw, formula = quad_formula("obesity"), model = "glm",
                          family = "binomial")
      )
      search <- function(...) {
        optimal_composition(fits, direction = "min", spacing = 30,
                            threshold = list(obesity = 0.8), ...)
      }
      runs <<- list(
        plain = suppressMessages(search()),
        widened = count_messages(search(widen = TRUE), "were widened to"),
        direct = suppressMessages(search(zone = 0.6)),
        step_above = suppressMessages(search(zone = 0.65))
      )
    }
    runs
  }
})
overall_zone <- function(oc) oc$zones[oc$zones$outcome == "overall", ]

test_that("a widened search equals a search run directly at the level it reached", {
  skip_on_cran()
  runs <- widen_searches()
  w <- runs$widened$value
  d <- runs$direct
  expect_equal(runs$widened$n, 1)
  expect_identical(w$zone_used, c(all = 0.6))
  expect_identical(w$candidates, d$candidates)
  expect_identical(w$zones, d$zones)
  expect_identical(w$optimum, d$optimum)
  settings <- setdiff(names(w), c("candidates", "zones", "optimum", "zone", "widen"))
  expect_identical(w[settings], d[settings])
})

test_that("widening stops at the first 5-point step at which the zones share a candidate", {
  skip_on_cran()
  runs <- widen_searches()
  expect_equal(overall_zone(runs$plain)$n, 0)
  expect_equal(overall_zone(runs$step_above)$n, 0)
  expect_gt(overall_zone(runs$direct)$n, 0)
  expect_equal(overall_zone(runs$widened$value)$n, overall_zone(runs$direct)$n)
})

test_that("the widened zone equals a hand computation from the ranks", {
  skip_on_cran()
  runs <- widen_searches()
  cand <- runs$plain$candidates
  w <- runs$widened$value
  feasible <- exp(cand$estimate_obesity) <= 0.8
  common <- pmin(cand$rank_zbmi, cand$rank_y_opp)
  q <- max(common[feasible])
  expect_lt(abs(q - overall_zone(runs$plain)$q_overlap), 1e-9)
  steps <- seq(95, 50, by = -5) / 100
  level <- steps[steps <= q][1]
  expect_identical(level, w$zone_used[["all"]])
  member <- feasible & common >= level
  expect_equal(overall_zone(w)$n, sum(member))
  expect_equal(n_xor(w$candidates$in_zone, member), 0)
  Xz <- as.matrix(cand[member, parts])
  expect_lt(max_diff(centre_of(w), hand_gm(Xz)), 1e-9)
  expect_lt(max_diff(overall_zone(w)[paste0("min_", parts)], apply(Xz, 2, min)), 1e-9)
  expect_lt(max_diff(overall_zone(w)[paste0("max_", parts)], apply(Xz, 2, max)), 1e-9)
})

test_that("a threshold-ruled outcome keeps its zone when the others are widened", {
  skip_on_cran()
  runs <- widen_searches()
  p <- runs$plain
  w <- runs$widened$value
  expect_identical(w$candidates$in_zone_obesity, p$candidates$in_zone_obesity)
  expect_identical(w$candidates$in_zone_obesity, exp(p$candidates$estimate_obesity) <= 0.8)
  expect_identical(w$zones[w$zones$outcome == "obesity", ],
                   p$zones[p$zones$outcome == "obesity", ])
  expect_gt(w$zones$n[w$zones$outcome == "zbmi"], p$zones$n[p$zones$outcome == "zbmi"])
})
