# Numerical benchmark of the curvature test that comp_iv() stores when the model
# holds the full quadratic in the ILR coordinates: the joint Wald test that
# every square and every pairwise product is zero. Compared with a hand
# b' V^-1 b on the fitted coefficients, car::linearHypothesis() with the
# chi-squared test (lm, glm, lmer), independent refits on coordinates built
# with robCompositions::pivotCoord() (lm, glm, lmer, Cox).
# The statistic is also checked against the order of the parts and the two ways
# of writing the quadratic.

skip_on_cran()

parts <- c("MVPA", "LPA", "SB", "SPT")
closure <- 1440
cd <- declare_composition(timeuse, parts, closure = closure)

fixed <- function(m) {
  if (inherits(m, "merMod")) lme4::fixef(m) else stats::coef(m)
}
quad_names <- function(nms) {
  nms[grepl("^I\\(ilr[0-9]+\\^2\\)$", nms) | grepl("^ilr[0-9]+:ilr[0-9]+$", nms)]
}
hand_wald <- function(m) {
  b <- fixed(m)
  V <- as.matrix(stats::vcov(m))
  nm <- quad_names(names(b))
  stat <- as.numeric(t(b[nm]) %*% solve(V[nm, nm]) %*% b[nm])
  list(stat = stat, df = length(nm),
       p = stats::pchisq(stat, length(nm), lower.tail = FALSE), names = nm)
}
car_wald <- function(m, nm) {
  b <- fixed(m)
  L <- matrix(0, length(nm), length(b), dimnames = list(nm, names(b)))
  L[cbind(nm, nm)] <- 1
  lh <- car::linearHypothesis(m, hypothesis.matrix = L,
                              rhs = rep(0, length(nm)), test = "Chisq")
  list(stat = lh[["Chisq"]][2], df = lh[["Df"]][2], p = lh[["Pr(>Chisq)"]][2])
}
expect_curvature <- function(ct, ref, tolerance) {
  expect_equal(ct$statistic, ref$stat, tolerance = tolerance)
  expect_equal(as.numeric(ct$df), as.numeric(ref$df), tolerance = tolerance)
  expect_equal(ct$p.value, ref$p, tolerance = tolerance)
}

# Coordinates built outside timecoda, the quadratic written term by term.
ilr_frame <- function(extra) {
  z <- robCompositions::pivotCoord(as.data.frame(tc_comp_df(cd)), pivotvar = 1)
  z <- as.data.frame(as.matrix(z))
  names(z) <- paste0("ilr", seq_len(ncol(z)))
  cbind(z, timeuse[extra])
}
quad_rhs <- local({
  co <- paste0("ilr", seq_len(length(parts) - 1))
  paste(c(co, paste0("I(", co, "^2)"),
          utils::combn(co, 2, paste, collapse = ":")), collapse = " + ")
})

fit_lm <- comp_iv(cd, "zbmi", model = "lm", curvature = TRUE)
fit_glm <- comp_iv(cd, "obesity", model = "glm", family = "binomial",
                   curvature = TRUE)
fit_cox <- comp_iv(cd, time = "time", event = "event", model = "cox",
                   curvature = TRUE)

test_that("the curvature test tests every square and product, six for four parts", {
  for (fit in list(fit_lm, fit_glm, fit_cox)) {
    expect_equal(fit$curvature_test$df, 6)
    expect_equal(length(quad_names(names(fixed(fit$model)))), 6L)
  }
})

test_that("the curvature test equals a hand Wald statistic on an lm", {
  expect_curvature(fit_lm$curvature_test, hand_wald(fit_lm$model), 1e-6)
})

test_that("the curvature test equals car::linearHypothesis() on an lm", {
  skip_if_not_installed("car")
  ref <- car_wald(fit_lm$model, hand_wald(fit_lm$model)$names)
  expect_curvature(fit_lm$curvature_test, ref, 1e-6)
})

test_that("the curvature test equals a refit of the lm on pivotCoord() coordinates", {
  m <- stats::lm(stats::as.formula(paste("zbmi ~", quad_rhs)),
                 data = ilr_frame("zbmi"))
  expect_curvature(fit_lm$curvature_test, hand_wald(m), 1e-6)
})

test_that("summary() carries the curvature test of the fit", {
  cols <- c("statistic", "df", "p.value")
  expect_equal(unlist(summary(fit_lm)$curvature_test[cols]),
               unlist(fit_lm$curvature_test[cols]), tolerance = 1e-6)
})

test_that("the curvature test equals a hand Wald statistic on a binomial glm", {
  expect_curvature(fit_glm$curvature_test, hand_wald(fit_glm$model), 1e-6)
})

test_that("the curvature test equals car::linearHypothesis() on a binomial glm", {
  skip_if_not_installed("car")
  ref <- car_wald(fit_glm$model, hand_wald(fit_glm$model)$names)
  expect_curvature(fit_glm$curvature_test, ref, 1e-6)
})

test_that("the curvature test equals a refit of the glm on pivotCoord() coordinates", {
  m <- stats::glm(stats::as.formula(paste("obesity ~", quad_rhs)),
                  data = ilr_frame("obesity"), family = stats::binomial())
  expect_curvature(fit_glm$curvature_test, hand_wald(m), 1e-3)
})

test_that("the curvature test equals car::linearHypothesis() and a hand Wald statistic on an lmer", {
  skip_if_not_installed("lme4")
  skip_if_not_installed("car")
  fit <- comp_iv(cd, "zbmi", model = "lmer", random = "school",
                 curvature = TRUE)
  h <- hand_wald(fit$model)
  expect_curvature(fit$curvature_test, h, 1e-6)
  expect_curvature(fit$curvature_test, car_wald(fit$model, h$names), 1e-6)
  m <- lme4::lmer(
    stats::as.formula(paste("zbmi ~", quad_rhs, "+ (1 | school)")),
    data = ilr_frame(c("zbmi", "school"))
  )
  expect_curvature(fit$curvature_test, hand_wald(m), 1e-3)
})

test_that("the curvature test equals a hand Wald statistic on a Cox fit and its refit", {
  expect_curvature(fit_cox$curvature_test, hand_wald(fit_cox$model), 1e-6)
  m <- survival::coxph(
    stats::as.formula(paste("survival::Surv(time, event) ~", quad_rhs)),
    data = ilr_frame(c("time", "event"))
  )
  expect_curvature(fit_cox$curvature_test, hand_wald(m), 1e-3)
})

test_that("the curvature test is the same for curvature = TRUE and the formula spelling", {
  cols <- c("statistic", "df", "p.value")
  same <- function(a, b) {
    expect_equal(unlist(a$curvature_test[cols]), unlist(b$curvature_test[cols]),
                 tolerance = 1e-6)
  }
  same(fit_lm, comp_iv(cd, formula = zbmi ~ ilr()^2 + I(ilr()^2), model = "lm"))
  same(fit_glm, comp_iv(cd, formula = obesity ~ ilr()^2 + I(ilr()^2),
                        model = "glm", family = "binomial"))
  same(fit_cox, comp_iv(cd, formula = survival::Surv(time, event) ~
                          ilr()^2 + I(ilr()^2), model = "cox"))
})

test_that("the curvature test does not depend on the order of the parts", {
  cols <- c("statistic", "df", "p.value")
  same <- function(a, b) {
    expect_equal(unlist(a$curvature_test[cols]), unlist(b$curvature_test[cols]),
                 tolerance = 1e-6)
  }
  for (ord in list(rev(parts), parts[c(3, 1, 4, 2)])) {
    cd_ord <- declare_composition(timeuse, ord, closure = closure)
    same(fit_lm, comp_iv(cd_ord, "zbmi", model = "lm", curvature = TRUE))
    same(fit_glm, comp_iv(cd_ord, "obesity", model = "glm", family = "binomial",
                          curvature = TRUE))
    same(fit_cox, comp_iv(cd_ord, time = "time", event = "event", model = "cox",
                          curvature = TRUE))
  }
})
