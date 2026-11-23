# Backend convergence capture + optimizer control.
#
# comp_iv() wraps every engine call so backend warnings/messages are recorded
# (but re-emitted) on the fit, probes per-class convergence status, and accepts
# an engine control object passed through `model_args`.
# These tests exercise the capture tibble, the convergence status, the summary
# and print surfacing, the model_args pass-through, and the fit_check row.

parts <- c("MVPA", "LPA", "SB", "SPT")
cd <- declare_composition(timeuse, parts, closure = 1440)

test_that("glm with a starved iteration budget captures a non-convergence warning", {
  fit <- suppressWarnings(
    comp_iv(cd, "obesity", model = "glm", family = "binomial",
            model_args = list(control = stats::glm.control(maxit = 1)))
  )
  expect_s3_class(fit$fit_warnings, "tbl_df")
  expect_setequal(names(fit$fit_warnings), c("source", "type", "message"))
  expect_gt(nrow(fit$fit_warnings), 0)
  # The canonical fit's condition is labelled "canonical".
  expect_true("canonical" %in% fit$fit_warnings$source)
  expect_false(isTRUE(fit$convergence$converged))
  expect_identical(fit$convergence$n_warnings, nrow(fit$fit_warnings))

  s <- summary(fit)
  expect_true("convergence" %in% names(s))
  expect_s3_class(s$convergence$status, "tbl_df")
  expect_gt(nrow(s$convergence$warnings), 0)

  out <- utils::capture.output(print(fit))
  expect_true(any(grepl("Convergence:", out)))
})

test_that("clean glm has zero-row fit_warnings and converged TRUE", {
  fit <- comp_iv(cd, "obesity", model = "glm", family = "binomial")
  expect_equal(nrow(fit$fit_warnings), 0)
  expect_true(isTRUE(fit$convergence$converged))
  expect_identical(fit$convergence$n_warnings, 0L)
  # No convergence line in print() when the fit is clean.
  out <- utils::capture.output(print(fit))
  expect_false(any(grepl("Convergence:", out)))
})

test_that("lm fit reports NA convergence and rejects unknown backend args", {
  fit <- comp_iv(cd, "zbmi", model = "lm")
  expect_true(is.na(fit$convergence$converged))
  expect_equal(nrow(fit$fit_warnings), 0)
  # `control` is no longer a comp_iv() argument: engine settings go through
  # model_args, so a top-level `control =` is an unused argument.
  expect_error(
    comp_iv(cd, "zbmi", model = "lm", control = list(maxit = 1)),
    "unused argument"
  )
})

test_that("clean glm summary prints 'Convergence: OK'", {
  fit <- comp_iv(cd, "obesity", model = "glm", family = "binomial")
  out <- utils::capture.output(print(summary(fit)))
  expect_true(any(grepl("Convergence: OK", out)))
})

test_that("lmer control threads the optimizer through and a clean fit is OK", {
  skip_if_not_installed("lme4")
  fit <- comp_iv(cd, "zbmi", model = "lmer", random = "school",
                 model_args = list(
                   control = lme4::lmerControl(optimizer = "bobyqa")
                 ))
  expect_identical(fit$model@optinfo$optimizer, "bobyqa")
  expect_equal(nrow(fit$fit_warnings), 0)
  # 'Convergence: OK' is emitted via cat() for a clean iterative fit.
  out <- utils::capture.output(suppressMessages(print(summary(fit))))
  expect_true(any(grepl("Convergence: OK", out)))
})

test_that("glmer with a starved optimizer budget captures a convergence warning", {
  skip_if_not_installed("lme4")
  fit <- suppressWarnings(
    comp_iv(cd, "obesity", model = "glmer", family = "binomial",
            random = "school",
            model_args = list(
              control = lme4::glmerControl(optCtrl = list(maxfun = 5))
            ))
  )
  # Either a captured warning or a flagged non-convergence/singularity.
  captured <- nrow(fit$fit_warnings) > 0
  flagged <- isFALSE(fit$convergence$converged) || isTRUE(fit$convergence$singular)
  expect_true(captured || flagged)
})

test_that("fit_check on the failed glm carries the fit's convergence state", {
  skip_if_not_installed("DHARMa")
  fit <- suppressWarnings(
    comp_iv(cd, "obesity", model = "glm", family = "binomial",
            model_args = list(control = stats::glm.control(maxit = 1)))
  )
  fc <- suppressWarnings(fit_check(fit, plot = FALSE))
  expect_false(any(grepl("^Convergence", fc$tests$check)))
  expect_identical(fc$convergence, fit$convergence)
  expect_equal(fc$convergence$n_warnings, nrow(fit$fit_warnings))
  out <- utils::capture.output(print(fc))
  expect_true(any(grepl("^  Convergence: NOT converged \\(", out)))
})
