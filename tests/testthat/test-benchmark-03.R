# Numerical benchmark of comp_iv() for model = "lm" and "glm", its omnibus
# Wald test and summary.comp_iv()'s per-part "vs remaining" effects, against
# stats::lm()/stats::glm() fitted by hand on robCompositions::pivotCoord()
# coordinates and against codaredistlm, whose internal lm uses the same pivot
# basis. Also checks that the omnibus and per-part effects do not depend on
# the order the parts are given in (D = 3, 4, 5). Reference calls are
# namespace-qualified so each number is visibly computed outside timecoda.

skip_if_not_installed("robCompositions")
skip_on_cran()

parts <- c("MVPA", "LPA", "SB", "SPT")
closure <- 1440
tu <- as.data.frame(timeuse)
cd <- declare_composition(timeuse, parts, closure = closure)
ilr_names <- paste0("ilr", 1:3)

hand_frame <- function(k = 1, data = tu, p = parts) {
  z <- robCompositions::pivotCoord(as.data.frame(data[p]), pivotvar = k)
  colnames(z) <- paste0("ilr", seq_len(length(p) - 1))
  cbind(data, z)
}
hand_fit <- function(kind, k = 1, cov = "sex") {
  f <- stats::as.formula(paste(
    if (kind == "lm") "zbmi" else "obesity", "~",
    paste(c(ilr_names, cov), collapse = " + ")
  ))
  if (kind == "lm") {
    stats::lm(f, data = hand_frame(k))
  } else {
    stats::glm(f, data = hand_frame(k), family = stats::binomial())
  }
}
hand_wald <- function(m) {
  b <- stats::coef(m)[ilr_names]
  as.numeric(t(b) %*% solve(stats::vcov(m)[ilr_names, ilr_names]) %*% b)
}

fit_lm <- comp_iv(cd, "zbmi", model = "lm", covariates = "sex")
fit_lm0 <- comp_iv(cd, "zbmi", model = "lm")
fit_glm <- comp_iv(cd, "obesity", model = "glm", covariates = "sex",
                   family = "binomial")

test_that("lm and glm fits match hand fits on pivot coordinates", {
  ref <- hand_fit("lm")
  expect_equal(coef(fit_lm$model), coef(ref), tolerance = 1e-6)
  expect_equal(vcov(fit_lm$model), vcov(ref), tolerance = 1e-6)
  ref0 <- hand_fit("lm", cov = NULL)
  expect_equal(coef(fit_lm0$model), coef(ref0), tolerance = 1e-6)
  expect_equal(vcov(fit_lm0$model), vcov(ref0), tolerance = 1e-6)
  refg <- hand_fit("glm")
  expect_equal(coef(fit_glm$model), coef(refg), tolerance = 1e-3)
  expect_equal(vcov(fit_glm$model), vcov(refg), tolerance = 1e-3)
})

test_that("the omnibus is the chi-squared Wald test of the ILR coefficients", {
  fits <- list(fit_lm, fit_lm0, fit_glm)
  refs <- list(hand_fit("lm"), hand_fit("lm", cov = NULL), hand_fit("glm"))
  tols <- c(1e-6, 1e-6, 1e-3)
  for (i in seq_along(fits)) {
    st <- hand_wald(refs[[i]])
    expect_equal(fits[[i]]$omnibus$statistic, st, tolerance = tols[i])
    expect_equal(fits[[i]]$omnibus$df, 3)
    expect_equal(fits[[i]]$omnibus$p.value,
                 stats::pchisq(st, 3, lower.tail = FALSE), tolerance = tols[i])
  }
})

test_that("per-part effects match hand refits at each pivot", {
  for (kind in c("lm", "glm")) {
    fit <- if (kind == "lm") fit_lm else fit_glm
    tol <- if (kind == "lm") 1e-6 else 1e-3
    for (lev in c(0.95, 0.90)) {
      pp <- summary(fit, conf_level = lev)$per_part
      for (k in seq_along(parts)) {
        m <- hand_fit(kind, k)
        ct <- stats::coef(summary(m))["ilr1", ]
        crit <- if (kind == "lm") {
          stats::qt(1 - (1 - lev) / 2, stats::df.residual(m))
        } else {
          stats::qnorm(1 - (1 - lev) / 2)
        }
        ref <- c(ct[[1]], ct[[1]] - crit * ct[[2]], ct[[1]] + crit * ct[[2]])
        if (kind == "glm") ref <- exp(ref)
        expect_equal(
          c(pp$estimate[k], pp$conf.low[k], pp$conf.high[k]), ref,
          tolerance = tol
        )
        expect_equal(pp$p.value[k], ct[[4]], tolerance = tol)
      }
    }
  }
})

test_that("the lm matches codaredistlm's internal fit (same pivot basis)", {
  skip_if_not_installed("codaredistlm")
  tu_num <- tu
  tu_num$male <- as.numeric(tu$sex == "male")
  fit <- comp_iv(declare_composition(tu_num, parts, closure = closure),
                 "zbmi", model = "lm", covariates = "male")
  utils::capture.output(suppressMessages(
    obj <- codaredistlm::predict_delta_comps(
      dataf = tu_num, y = "zbmi", comps = parts, covars = "male",
      deltas = c(0, 10) / 1440, comparisons = "one-v-one"
    )
  ))
  ref <- attr(obj, "lm")
  cf <- coef(fit$model)
  expect_equal(cf, coef(ref)[names(cf)], tolerance = 1e-6)
  expect_equal(vcov(fit$model), vcov(ref)[names(cf), names(cf)],
               tolerance = 1e-6)
  expect_equal(unname(fitted(fit$model)), unname(fitted(ref)), tolerance = 1e-6)
  # Prediction at the package's reference composition (the geometric mean).
  nd <- as.data.frame(fit$ref_comp)
  nd <- cbind(nd, stats::setNames(
    as.data.frame(robCompositions::pivotCoord(nd[parts])), ilr_names
  ), male = mean(tu_num$male))
  expect_equal(unname(stats::predict(fit$model, newdata = nd)),
               unname(attr(obj, "mean_pred")$fit), tolerance = 1e-6)
})

test_that("omnibus and per-part effects do not depend on the part order", {
  set.seed(2026)
  tu5 <- tu
  u <- stats::runif(nrow(tu5), 0.3, 0.7)
  tu5$SB1 <- tu5$SB * u
  tu5$SB2 <- tu5$SB * (1 - u)
  sets <- list(c("MVPA", "LPA", "SB"), parts,
               c("MVPA", "LPA", "SB1", "SB2", "SPT"))
  for (p in sets) {
    dat <- if (length(p) == 5) tu5 else tu
    fits <- purrr::map(list(p, rev(p), c(p[-1], p[1])), function(pp) {
      cdx <- declare_composition(dat, pp, closure = closure)
      list(lm = comp_iv(cdx, "zbmi", model = "lm", covariates = "sex"),
           glm = comp_iv(cdx, "obesity", model = "glm", covariates = "sex",
                         family = "binomial"))
    })
    for (kind in c("lm", "glm")) {
      est <- purrr::map(fits, function(f) {
        s <- summary(f[[kind]])$per_part
        stats::setNames(s$estimate, s$parameter)[paste(p, "vs remaining")]
      })
      for (j in 2:3) {
        expect_equal(fits[[j]][[kind]]$omnibus$statistic,
                     fits[[1]][[kind]]$omnibus$statistic, tolerance = 1e-8)
        expect_equal(est[[j]], est[[1]], tolerance = 1e-8)
      }
    }
  }
})

test_that("the glm per-part estimate is exp() of each rotated fit's ilr1", {
  s <- summary(fit_glm)
  ref <- purrr::map_dbl(parts, function(p) {
    exp(stats::coef(fit_glm$models[[p]])[["ilr1"]])
  })
  expect_equal(s$per_part$estimate, ref, tolerance = 1e-12)
})
