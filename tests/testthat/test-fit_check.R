# fit_check() goodness-of-fit diagnostics. Diagnostic backends are optional
# (Suggests), so each backend test skips when its package is absent. Plotting is
# a side effect; setup.R opens a null device so plots don't pollute the run.

parts <- c("MVPA", "LPA", "SB", "SPT")

make_fixture <- function(seed = 1) {
  set.seed(seed)
  d <- timeuse
  p <- tc_comp_df(declare_composition(d, parts, closure = 1440))
  lp <- 0.9 * log(p$MVPA / p$SB) - 0.4 * log(p$LPA / p$SPT) +
    as.numeric(d$school) * 0.1
  d$grade <- ordered(cut(lp + stats::rnorm(nrow(d)), 4, labels = FALSE))
  d$catf <- factor(ifelse(lp + stats::rnorm(nrow(d)) < -0.3, "lo",
                          ifelse(lp + stats::rnorm(nrow(d)) > 0.3, "hi", "mid")))
  d$rz <- lp + stats::rnorm(nrow(d))
  declare_composition(d, parts, closure = 1440)
}
cd <- make_fixture()

expect_fit_check <- function(fc, backend) {
  expect_s3_class(fc, "fit_check")
  expect_identical(fc$backend, backend)
  expect_s3_class(fc$tests, "tbl_df")
  expect_setequal(names(fc$tests), c("check", "statistic", "df", "p.value"))
}

test_that("registry assigns a check backend to every model type", {
  expect_identical(tc_model_meta("glm")$check_backend, "dharma")
  expect_identical(tc_model_meta("cox")$check_backend, "coxzph")
  expect_identical(tc_model_meta("clm")$check_backend, "ordinal")
  expect_identical(tc_model_meta("clmm")$check_backend, "limited")
  expect_identical(tc_model_meta("lmrob")$check_backend, "robust")
  expect_identical(tc_model_meta("rlmer")$check_backend, "robust_mixed")
  expect_identical(tc_model_meta("mblogit")$check_backend, "limited")
})

test_that("fit_check rejects non-model input", {
  expect_error(fit_check(1:5), "comp_iv")
})

test_that("DHARMa backend: glm / lm / lmer", {
  skip_if_not_installed("DHARMa")
  m <- comp_iv(cd, "obesity", model = "glm", family = "binomial")
  fc <- fit_check(m, plot = FALSE)
  expect_fit_check(fc, "dharma")
  # the table holds the backend's own checks; convergence is its own element
  expect_identical(fc$tests$check, c("KS uniformity", "Dispersion", "Outliers"))
  expect_identical(fc$convergence, m$convergence)
  expect_s3_class(fc$raw, "DHARMa")
  # plotting path runs without error (null device from setup.R)
  expect_no_error(quiet_smooth(fit_check(m, plot = TRUE)))

  ml <- comp_iv(cd, "zbmi", model = "lmer", random = "school")
  expect_fit_check(fit_check(ml, plot = FALSE), "dharma")
})

test_that("cox backend: proportional-hazards test with a GLOBAL row", {
  m <- comp_iv(cd, time = "time", event = "event", model = "cox")
  fc <- fit_check(m, plot = FALSE)
  expect_fit_check(fc, "coxzph")
  expect_true(any(grepl("GLOBAL", fc$tests$check)))
  # one row per ilr term plus GLOBAL
  expect_equal(nrow(fc$tests), length(parts) - 1 + 1)
})

test_that("ordinal backend: surrogate residuals + fit stats", {
  skip_if_not_installed("sure")
  skip_if_not_installed("ordinal")
  m <- comp_iv(cd, "grade", model = "clm")
  fc <- fit_check(m, plot = FALSE)
  expect_fit_check(fc, "ordinal")
  expect_true("AIC" %in% fc$tests$check)
  expect_true(is.numeric(fc$raw$surrogate) &&
                length(fc$raw$surrogate) == nrow(cd$data_raw))
  expect_true(is.character(fc$note) && grepl("Proportional-odds", fc$note))
  expect_match(fc$note, "is evidence that", fixed = TRUE)
  expect_no_match(fc$note, "means", fixed = TRUE)
  # The surrogate draws repeat and leave the session's random stream alone.
  set.seed(42)
  before <- .Random.seed
  again <- fit_check(m, plot = FALSE)
  expect_identical(.Random.seed, before)
  expect_identical(again$raw$surrogate, fc$raw$surrogate)
})

test_that("ordinal backend: the proportional-odds test has a row per ILR coordinate", {
  skip_if_not_installed("sure")
  skip_if_not_installed("ordinal")
  m <- comp_iv(cd, "grade", model = "clm")
  fc <- fit_check(m, plot = FALSE)
  po <- fc$tests[grepl("^Proportional odds: ", fc$tests$check), ]
  expect_identical(po$check, paste0("Proportional odds: ilr", 1:3))
  expect_true(all(is.finite(po$p.value)))
  # The same test run directly on the pivot coordinates. The data go into the
  # call itself, so the refits inside nominal_test() can find them.
  direct <- ordinal::nominal_test(do.call(
    ordinal::clm,
    list(formula = grade ~ ilr1 + ilr2 + ilr3, data = m$data, Hess = TRUE)
  ))
  expect_equal(po$p.value, direct[["Pr(>Chi)"]][-1], tolerance = 1e-6)
})

test_that("robust backends: weight summaries for lmrob and rlmer", {
  skip_if_not_installed("robustbase")
  lr <- comp_iv(cd, "rz", model = "lmrob")
  fc <- fit_check(lr, plot = FALSE)
  expect_fit_check(fc, "robust")
  expect_true(any(grepl("robustness weight", fc$tests$check, ignore.case = TRUE)))

  skip_if_not_installed("robustlmm")
  rl <- comp_iv(cd, "rz", model = "rlmer", random = "school")
  expect_fit_check(fit_check(rl, plot = FALSE), "robust_mixed")
})

test_that("limited backends (clmm, mblogit) return fit stats + a note", {
  skip_if_not_installed("ordinal")
  cm <- comp_iv(cd, "grade", model = "clmm", random = "school")
  fc <- fit_check(cm, plot = FALSE)
  expect_fit_check(fc, "limited")
  # fit stats only; the convergence state is its own element
  expect_identical(fc$tests$check, c("AIC", "log-likelihood"))
  expect_identical(fc$convergence, cm$convergence)
  expect_true(is.character(fc$note) && grepl("clmm|mixed ordinal", fc$note))

  skip_if_not_installed("mclogit")
  mb <- suppressWarnings(comp_iv(cd, "catf", model = "mblogit"))
  fcm <- fit_check(mb, plot = FALSE)
  expect_fit_check(fcm, "limited")
  expect_true(grepl("multinomial", fcm$note))
})

test_that("the convergence line is omitted for lm and reads NOT converged for a failed fit", {
  skip_if_not_installed("DHARMa")
  m <- comp_iv(cd, "zbmi", model = "lm")
  fc <- fit_check(m, plot = FALSE)
  expect_identical(fc$convergence, m$convergence)
  expect_true(is.na(fc$convergence$converged))
  expect_false(any(grepl("Convergence", fc$tests$check)))
  expect_false(any(grepl("Convergence", utils::capture.output(print(fc)))))

  ok <- fit_check(comp_iv(cd, "obesity", model = "glm", family = "binomial"),
                  plot = FALSE)
  expect_identical(utils::capture.output(print(ok))[3],
                   "  Convergence: converged (0 backend warnings)")

  starved <- suppressWarnings(
    comp_iv(cd, "obesity", model = "glm", family = "binomial",
            model_args = list(control = stats::glm.control(maxit = 1)))
  )
  fcs <- suppressWarnings(fit_check(starved, plot = FALSE))
  expect_false(any(grepl("Convergence", fcs$tests$check)))
  expect_identical(fcs$convergence, starved$convergence)
  n <- starved$convergence$n_warnings
  expect_identical(
    utils::capture.output(print(fcs))[3],
    paste0("  Convergence: NOT converged (", n, " backend warning",
           if (n == 1) "" else "s", ")")
  )
  # one warning reads in the singular
  one <- fcs
  one$convergence$n_warnings <- 1L
  expect_identical(utils::capture.output(print(one))[3],
                   "  Convergence: NOT converged (1 backend warning)")
})

test_that("n sets the number of DHARMa simulations", {
  skip_if_not_installed("DHARMa")
  m <- comp_iv(cd, "zbmi", model = "lm")
  expect_equal(fit_check(m, plot = FALSE, n = 50)$raw$nSim, 50)
  expect_equal(fit_check(m, plot = FALSE)$raw$nSim, 250)
})

test_that("print() reports the model, backend and note", {
  skip_if_not_installed("DHARMa")
  fc <- fit_check(comp_iv(cd, "zbmi", model = "lm"), plot = FALSE)
  out <- utils::capture.output(print(fc))
  expect_identical(out[1:3], c(
    "<fit_check> linear model ",
    "  Backend: dharma ",
    ""
  ))
  expect_false(any(grepl("Plot:", out)))
  expect_false(any(grepl("Note:", out)))
  expect_true(any(grepl("KS uniformity", out)))

  skip_if_not_installed("ordinal")
  fcl <- fit_check(comp_iv(cd, "grade", model = "clmm", random = "school"), plot = FALSE)
  lim_out <- utils::capture.output(print(fcl))
  expect_match(lim_out[3], "^  Convergence: ")
  expect_identical(lim_out[4], paste0("  Note: ", fcl$note, " "))
})

test_that("plot = TRUE prints the report; plot = FALSE prints nothing", {
  skip_if_not_installed("DHARMa")
  m <- comp_iv(cd, "obesity", model = "glm", family = "binomial")
  expect_identical(utils::capture.output(fit_check(m, plot = FALSE)), character(0))
  shown <- utils::capture.output(quiet_smooth(fit_check(m, plot = TRUE)))
  expect_identical(shown[1], "<fit_check> generalized linear model ")
})

test_that("plot = TRUE runs on the cox, ordinal, robust and limited backends", {
  m_cox <- comp_iv(cd, time = "time", event = "event", model = "cox")
  expect_no_error(utils::capture.output(fit_check(m_cox, plot = TRUE)))

  skip_if_not_installed("ordinal")
  cm <- comp_iv(cd, "grade", model = "clmm", random = "school")
  utils::capture.output(fc_plot <- fit_check(cm, plot = TRUE))
  expect_identical(fc_plot$tests, fit_check(cm, plot = FALSE)$tests)

  skip_if_not_installed("sure")
  m_clm <- comp_iv(cd, "grade", model = "clm")
  expect_no_error(utils::capture.output(fit_check(m_clm, plot = TRUE)))

  skip_if_not_installed("robustbase")
  lr <- comp_iv(cd, "rz", model = "lmrob")
  expect_no_error(utils::capture.output(fit_check(lr, plot = TRUE)))

  skip_if_not_installed("robustlmm")
  rl <- comp_iv(cd, "rz", model = "rlmer", random = "school")
  expect_no_error(utils::capture.output(fit_check(rl, plot = TRUE)))
})

test_that("plot must be a single TRUE or FALSE", {
  m <- comp_iv(cd, "zbmi", model = "lm")
  bad <- list("yes", 1, c(TRUE, FALSE))
  purrr::walk(bad, function(v) {
    expect_error(fit_check(m, plot = v), "plot. must be TRUE or FALSE, not")
  })
})

test_that("seed reaches DHARMa, whose own default seed is 123", {
  skip_if_not_installed("DHARMa")
  m <- comp_iv(cd, "zbmi", model = "lm")
  s1 <- fit_check(m, plot = FALSE, seed = 1)$tests
  s2 <- fit_check(m, plot = FALSE, seed = 2)$tests
  expect_false(identical(s1, s2))
  expect_identical(fit_check(m, plot = FALSE, seed = 123)$tests,
                   fit_check(m, plot = FALSE)$tests)
})

test_that("DHARMa results do not depend on whether the session has drawn a random number", {
  skip_if_not_installed("DHARMa")
  m <- comp_iv(cd, "obesity", model = "glm", family = "binomial",
               covariates = c("sex", "zbmi"))
  set.seed(1)
  used <- fit_check(m, plot = FALSE)$tests
  # A fresh session has no random state; DHARMa used to draw once more there.
  rm(".Random.seed", envir = globalenv())
  fresh <- fit_check(m, plot = FALSE)$tests
  expect_identical(fresh, used)
  set.seed(99)
  expect_identical(fit_check(m, plot = FALSE)$tests, used)
})

test_that("the bootstrapped outlier test runs under the fixed seed as well", {
  skip_if_not_installed("DHARMa")
  d <- timeuse
  d$count <- as.integer(round(exp(pmin(pmax(d$zbmi, -3), 3))))
  fit <- comp_iv(declare_composition(d, parts, closure = 1440), "count",
                 model = "glmer", family = "poisson", random = "school")
  set.seed(1)
  t1 <- fit_check(fit, plot = FALSE)$tests
  set.seed(2)
  t2 <- fit_check(fit, plot = FALSE)$tests
  # Only a fit with outliers exercises the bootstrap.
  expect_gt(t1$statistic[t1$check == "Outliers"], 0)
  expect_identical(t1, t2)
})
