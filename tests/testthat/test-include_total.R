# Tests for the Coenders total-as-covariate ("T-space") model:
# comp_iv(include_total = TRUE). Compositions with a variable row total carry
# both relative (balance) and absolute (total) information; include_total adds
# the normalized log geometric mean of the zero-free PRE-CLOSURE parts,
# t = sqrt(D) * log((x1...xD)^(1/D)) (Coenders et al. 2017, eq. 3.2), as a
# `t_total` covariate.

# A variable-total 4-part composition: a per-row lognormal scaling factor makes
# the daily total vary, so the total genuinely carries information. Two outcomes
# are simulated: one that depends on both a balance and the log total (signal),
# and one that depends on the balance only (no total signal).
make_vartotal_data <- function(seed = 42, n = 300) {
  set.seed(seed)
  base <- cbind(
    MVPA = stats::rlnorm(n, log(60), 0.4),
    LPA  = stats::rlnorm(n, log(200), 0.4),
    SB   = stats::rlnorm(n, log(500), 0.4),
    SPT  = stats::rlnorm(n, log(480), 0.3)
  )
  sc <- stats::rlnorm(n, 0, 0.3)     # variable total per row
  raw <- base * sc
  df <- as.data.frame(raw)
  balance <- as.numeric(scale(log(raw[, "MVPA"]) - log(raw[, "SB"])))
  logtot <- as.numeric(scale(log(rowSums(raw))))
  df$y_signal <- 0.6 * balance + 1.2 * logtot + stats::rnorm(n)
  df$y_nosig <- 0.8 * balance + stats::rnorm(n)
  df$sex <- factor(sample(c("F", "M"), n, replace = TRUE))
  df
}

vt_parts <- c("MVPA", "LPA", "SB", "SPT")

test_that("include_total adds a t_total coefficient and a total_test", {
  df <- make_vartotal_data()
  cd <- declare_composition(df, vt_parts, closure = "mean")
  fit <- comp_iv(cd, outcome = "y_signal", model = "lm", include_total = TRUE)

  expect_true("t_total" %in% names(stats::coef(fit$model)))
  expect_s3_class(fit$total_test, "tbl_df")
  expect_identical(fit$total_test$term, "total")
  expect_identical(fit$total_test$df, 1L)
  expect_false(is.null(fit$total_ref))
  # t_total is sqrt(D) * the log geometric mean of the RAW (pre-closure)
  # parts (Coenders et al. 2017, eq. 3.2); total_ref stores its mean.
  expect_equal(
    fit$total_ref,
    mean(sqrt(length(vt_parts)) * rowMeans(log(df[vt_parts])))
  )
  # That mean is t_total at the geometric-mean composition of the raw parts.
  geo <- purrr::map_dbl(df[vt_parts], function(x) exp(mean(log(x))))
  expect_equal(fit$total_ref, sqrt(length(vt_parts)) * mean(log(geo)))
})

test_that("include_total works when the raw data contain zeroes (lrEM path)", {
  df <- make_vartotal_data(seed = 7)
  df$MVPA[1:8] <- 0
  cd <- declare_composition(df, vt_parts, closure = "mean")
  fit <- comp_iv(cd, outcome = "y_signal", model = "lm", include_total = TRUE)
  expect_true("t_total" %in% names(stats::coef(fit$model)))
  expect_true(is.finite(fit$total_ref))
  # t_total comes from the compD's zero-free `imputed` element (original
  # scale), matching the prototype's TSpace construction up to scaling.
  expect_equal(
    fit$total_ref,
    mean(sqrt(length(vt_parts)) * rowMeans(log(as.data.frame(cd$parts_imputed))))
  )
})

test_that("rescaling t_total leaves the total test invariant", {
  # Any positive rescaling of t (the paper's sqrt(D) vs a sqrt(D-1)
  # convention) must give the same p-value, with proportional coefficients.
  df <- make_vartotal_data(seed = 11)
  d <- length(vt_parts)
  df$t_alt <- sqrt(d - 1) * rowMeans(log(df[vt_parts]))
  cd <- declare_composition(df, vt_parts, closure = "mean")
  fit <- comp_iv(cd, outcome = "y_signal", model = "lm", include_total = TRUE)
  fit2 <- comp_iv(cd, outcome = "y_signal", model = "lm", covariates = "t_alt")
  b1 <- stats::coef(fit$model)[["t_total"]]
  b2 <- stats::coef(fit2$model)[["t_alt"]]
  expect_equal(b1 * sqrt(d), b2 * sqrt(d - 1), tolerance = 1e-8)
  expect_equal(
    fit$total_test$p.value,
    summary(fit2$model)$coefficients["t_alt", "Pr(>|t|)"],
    tolerance = 1e-8
  )
})

test_that("total_test rejects H0 when the total carries signal", {
  df <- make_vartotal_data()
  cd <- declare_composition(df, vt_parts, closure = "mean")
  fit <- comp_iv(cd, outcome = "y_signal", model = "lm", include_total = TRUE)
  expect_lt(fit$total_test$p.value, 0.05)
})

test_that("total_test does not reject H0 when the total carries no signal", {
  df <- make_vartotal_data()
  cd <- declare_composition(df, vt_parts, closure = "mean")
  fit <- comp_iv(cd, outcome = "y_nosig", model = "lm", include_total = TRUE)
  expect_gt(fit$total_test$p.value, 0.05)
})

test_that("summary() surfaces total_test and prints a header line", {
  df <- make_vartotal_data()
  cd <- declare_composition(df, vt_parts, closure = "mean")
  fit <- comp_iv(cd, outcome = "y_signal", model = "lm", include_total = TRUE)
  s <- summary(fit)
  expect_identical(s$total_test, fit$total_test)
  expect_false("total_verdict" %in% names(s))
  report <- paste(capture.output(print(s)), collapse = "\n")
  expect_match(report, "Total (log geometric mean, Wald): beta", fixed = TRUE)
  expect_no_match(report, "matters", fixed = TRUE)
  expect_no_match(report, "suffices", fixed = TRUE)
})

test_that("shipped constant-total timeuse errors under include_total", {
  cd <- declare_composition(timeuse, vt_parts, closure = 1440)
  expect_error(
    comp_iv(cd, outcome = "obesity", model = "glm", family = "binomial",
            include_total = TRUE),
    "variable"
  )
})

test_that("include_total = FALSE is unchanged (no t_total, NULL total_test)", {
  df <- make_vartotal_data()
  cd <- declare_composition(df, vt_parts, closure = "mean")
  fit0 <- comp_iv(cd, outcome = "y_signal", model = "lm")
  fitF <- comp_iv(cd, outcome = "y_signal", model = "lm", include_total = FALSE)

  expect_false("t_total" %in% names(stats::coef(fit0$model)))
  expect_null(fit0$total_test)
  expect_null(fit0$total_ref)
  # Explicit FALSE reproduces the default fit exactly (same coefficients).
  expect_equal(stats::coef(fit0$model), stats::coef(fitF$model))
  expect_identical(names(stats::coef(fit0$model)),
                   c("(Intercept)", "ilr1", "ilr2", "ilr3"))
})

test_that("reallocation on an include_total fit runs and is finite", {
  df <- make_vartotal_data()
  cd <- declare_composition(df, vt_parts, closure = "mean")
  fit <- comp_iv(cd, outcome = "y_signal", model = "lm", include_total = TRUE)
  pr <- predict_realloc(fit, delta = 30, from = "SB", to = "MVPA")
  expect_s3_class(pr, "realloc")
  expect_true(all(is.finite(pr$estimate)))
  expect_true(all(is.finite(pr$conf.low)))
  expect_true(all(is.finite(pr$conf.high)))
})

test_that("reallocation pins t_total: constant-total contrast matches no-total fit when total carries no signal", {
  # When the total carries no signal, holding it constant (include_total pins it
  # at the reference value) must leave the reallocation estimate essentially
  # equal to the plain closed-composition fit's estimate.
  df <- make_vartotal_data()
  cd <- declare_composition(df, vt_parts, closure = "mean")
  fit_it <- comp_iv(cd, outcome = "y_nosig", model = "lm", include_total = TRUE)
  fit_no <- comp_iv(cd, outcome = "y_nosig", model = "lm", include_total = FALSE)
  pr_it <- predict_realloc(fit_it, delta = 30, from = "SB", to = "MVPA")
  pr_no <- predict_realloc(fit_no, delta = 30, from = "SB", to = "MVPA")
  expect_equal(pr_it$estimate, pr_no$estimate, tolerance = 0.05)
})

test_that("reallocation grid pins t_total at the reference value", {
  # Build the emmeans reference grid the way tc_emm_at_points() does, and check
  # that t_total is pinned at the fit's stored reference (mean log geo-mean of
  # the raw parts) rather than averaged over the data. Inspecting the ref_grid
  # (before emmeans() collapses the nuisance coordinates) is what exposes the
  # pinned covariate as a grid level.
  df <- make_vartotal_data()
  cd <- declare_composition(df, vt_parts, closure = "mean")
  fit <- comp_iv(cd, outcome = "y_signal", model = "lm", include_total = TRUE)

  ref <- unlist(fit$ref_comp)
  ilr <- tc_pivot_coords(as.data.frame(rbind(ref)), pivot = fit$pivot)
  rg <- emmeans::ref_grid(
    fit$model,
    at = c(as.list(ilr[1, ]), list(t_total = fit$total_ref)),
    data = fit$data
  )
  grid <- as.data.frame(rg@grid)
  expect_true("t_total" %in% names(grid))
  expect_true(all(grid$t_total == fit$total_ref))
})

test_that("rotated per-part fits also carry t_total", {
  df <- make_vartotal_data()
  cd <- declare_composition(df, vt_parts, closure = "mean")
  fit <- comp_iv(cd, outcome = "y_signal", model = "lm", include_total = TRUE)
  for (p in vt_parts) {
    expect_true("t_total" %in% names(stats::coef(fit$models[[p]])),
                info = paste("part", p))
  }
})

test_that("include_total works through the formula interface", {
  df <- make_vartotal_data()
  cd <- declare_composition(df, vt_parts, closure = "mean")
  fit <- comp_iv(cd, y_signal ~ ilr() + sex, include_total = TRUE)
  expect_true("t_total" %in% names(stats::coef(fit$model)))
  expect_lt(fit$total_test$p.value, 0.05)
})

test_that("include_total works with a GLM (binomial) engine", {
  df <- make_vartotal_data()
  df$obesity <- as.integer(df$y_signal > stats::median(df$y_signal))
  cd <- declare_composition(df, vt_parts, closure = "mean")
  fit <- comp_iv(cd, outcome = "obesity", model = "glm", family = "binomial",
                 include_total = TRUE)
  expect_true("t_total" %in% names(stats::coef(fit$model)))
  pr <- predict_realloc(fit, delta = 20, from = "SB", to = "MVPA")
  expect_true(all(is.finite(pr$estimate)))
})

test_that("t_total stays out of the omnibus of a curvature formula", {
  df <- make_vartotal_data()
  cd <- declare_composition(df, vt_parts, closure = "mean")
  fit <- comp_iv(cd, y_signal ~ ilr() + I(ilr()^2), include_total = TRUE)
  expect_identical(sum(names(stats::coef(fit$model)) == "t_total"), 1L)
  expect_identical(fit$omnibus$df, 2L * (length(vt_parts) - 1L))
})

test_that("a multinomial fit tests one t_total coefficient per category", {
  skip_if_not_installed("mclogit")
  df <- make_vartotal_data()
  df$band <- cut(df$y_signal, 3, labels = c("low", "mid", "high"))
  cd <- declare_composition(df, vt_parts, closure = "mean")
  fit <- comp_iv(cd, outcome = "band", model = "mblogit", include_total = TRUE)
  expect_identical(fit$total_test$df, 2L)
})

test_that("include_total refuses a model that names t_total itself", {
  df <- make_vartotal_data()
  df$t_total <- 1
  cd <- declare_composition(df, vt_parts, closure = "mean")
  expect_error(
    comp_iv(cd, outcome = "y_signal", covariates = "t_total",
            include_total = TRUE),
    "cannot also name"
  )
  expect_error(
    comp_iv(cd, y_signal ~ ilr() + t_total, include_total = TRUE),
    "cannot also name"
  )
  expect_error(
    comp_iv(cd, outcome = "y_signal", model = "lm", include_total = TRUE,
            model_args = list(weights = t_total)),
    "cannot also name"
  )
  # The refusal reads the name itself, so a scalar of that name in the calling
  # environment, which the formula parser would substitute, is refused too.
  t_total <- 2
  cd0 <- declare_composition(make_vartotal_data(), vt_parts, closure = "mean")
  expect_error(
    comp_iv(cd0, y_signal ~ ilr() + t_total, include_total = TRUE),
    "cannot also name"
  )
})

test_that("a one-row composition is refused for a reason that names it", {
  df <- make_vartotal_data()[1, ]
  cd <- declare_composition(df, vt_parts, closure = "mean")
  expect_error(
    comp_iv(cd, outcome = "y_signal", model = "lm", include_total = TRUE),
    "fewer than two observations"
  )
})
