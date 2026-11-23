# Extra model families (ordinal / multinomial / robust / additive) and the
# tc_model_meta() registry that drives them. Backends are optional (Suggests),
# so every engine test skips when its package is absent.

parts <- c("MVPA", "LPA", "SB", "SPT")

# A fixture with signal-bearing outcomes derived from the composition, so the
# fits are well conditioned (a degenerate ordinal fit fails to compute vcov).
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

test_that("tc_model_meta returns filled defaults and rejects unknown types", {
  m <- tc_model_meta("clm")
  expect_true(all(c("label", "scale", "emm_mode", "per_category") %in%
                    names(m)))
  expect_identical(m$emm_mode, "latent")
  expect_identical(tc_model_meta("lm")$emm_via, "ref_grid")
  expect_false(tc_model_meta("lmrob")$per_category)
  expect_error(tc_model_meta("nope"), "Unknown model type")
})

test_that("effect scale and measure match each family", {
  expect_false(tc_report_meta("lmrob", NULL)$exponentiate)
  expect_identical(tc_report_meta("lmrob", NULL)$measure, "Mean difference")
  expect_true(tc_report_meta("clm", NULL)$exponentiate)
  expect_identical(tc_report_meta("clm", NULL)$measure, "Cumulative odds ratio")
  expect_identical(tc_report_meta("mblogit", NULL)$measure, "Relative risk ratio")
  # existing types unchanged
  expect_identical(tc_report_meta("glm", "binomial")$measure, "Odds ratio")
  expect_identical(tc_report_meta("cox", NULL)$measure, "Hazard ratio")
})

test_that("ordinal: clm and clmm fit, summarize, and reallocate", {
  skip_if_not_installed("ordinal")
  clm <- comp_iv(cd, outcome = "grade", model = "clm")
  expect_identical(clm$model_type, "clm")
  expect_equal(clm$omnibus$df, length(parts) - 1)
  s <- summary(clm)
  expect_equal(nrow(s$per_part), length(parts))
  expect_true(all(s$per_part$estimate > 0))          # odds ratios
  rr <- predict_realloc(clm, delta = 30, from = "SB", to = "MVPA")
  expect_true(is.finite(rr$estimate))
  # reference vs reference is OR ~ 1
  expect_equal(predict_realloc(clm, delta = 1e-6, from = "SB",
                               to = "MVPA")$estimate, 1, tolerance = 1e-4)

  clmm <- comp_iv(cd, outcome = "grade", model = "clmm", random = "school")
  expect_identical(clmm$model_type, "clmm")
  expect_true(is.finite(predict_realloc(clmm, delta = 20, from = "SB",
                                        to = "MVPA")$estimate))
})

test_that("ordinal per-part OR matches a hand-rotated clm refit", {
  skip_if_not_installed("ordinal")
  m <- comp_iv(cd, outcome = "grade", model = "clm")
  s <- summary(m)
  df <- dplyr::bind_cols(tc_comp_df(cd), cd$data_raw["grade"])
  d2 <- ilr_transform(df, parts, total = 1440, pivot = 2)
  f2 <- ordinal::clm(grade ~ ilr1 + ilr2 + ilr3, data = d2)
  expect_equal(s$per_part$estimate[2],
               unname(exp(stats::coef(f2)["ilr1"])), tolerance = 1e-6)
})

test_that("ordinal requires an ordered factor and rejects absolute predictions", {
  skip_if_not_installed("ordinal")
  expect_error(comp_iv(cd, outcome = "rz", model = "clm"), "ordered factor")
  m <- comp_iv(cd, outcome = "grade", model = "clm")
  expect_error(
    predict_realloc(m, delta = 20, from = "SB", to = "MVPA", type = "absolute"),
    "relative"
  )
})

test_that("robust linear: lmrob and rlmer fit and reallocate", {
  skip_if_not_installed("robustbase")
  lr <- comp_iv(cd, outcome = "rz", model = "lmrob")
  expect_identical(lr$model_type, "lmrob")
  expect_equal(lr$omnibus$df, length(parts) - 1)
  s <- summary(lr)
  expect_equal(nrow(s$per_part), length(parts))
  # lmrob has no native emmeans method: the qdrg path must still predict
  expect_true(is.finite(predict_realloc(lr, delta = 30, from = "SB",
                                        to = "MVPA")$estimate))

  skip_if_not_installed("robustlmm")
  rl <- comp_iv(cd, outcome = "rz", model = "rlmer", random = "school")
  expect_identical(rl$model_type, "rlmer")
  expect_true(is.finite(predict_realloc(rl, delta = 20, from = "SB",
                                        to = "MVPA")$estimate))
})

test_that("multinomial carries response_level through the whole stack", {
  skip_if_not_installed("mclogit")
  m <- suppressWarnings(comp_iv(cd, outcome = "catf", model = "mblogit"))
  expect_identical(m$model_type, "mblogit")
  k <- nlevels(cd$data_raw$catf)
  # omnibus df = (K - 1)(D - 1)
  expect_equal(m$omnibus$df, (k - 1) * (length(parts) - 1))

  s <- summary(m)
  expect_true("response_level" %in% names(s$per_part))
  expect_equal(nrow(s$per_part), length(parts) * (k - 1))
  expect_true(all(s$per_part$estimate > 0))          # relative-risk ratios

  rr <- suppressWarnings(
    predict_realloc(m, delta = c(15, 30), from = "SB", to = "MVPA",
                    on_invalid = "na")
  )
  expect_true("response_level" %in% names(rr))
  expect_equal(nrow(rr), 2 * (k - 1))
  # reference vs reference RRR ~ 1 for every category
  r0 <- suppressWarnings(
    predict_realloc(m, delta = 1e-6, from = "SB", to = "MVPA")
  )
  expect_equal(r0$estimate, rep(1, k - 1), tolerance = 1e-4)

  wide <- summary(rr)$table
  expect_true("response_level" %in% names(wide))
})

test_that("the registry drives the family / random argument rules", {
  expect_error(comp_iv(cd, outcome = "rz", model = "lmrob", family = "gaussian"),
               "does not apply")
  # a clustered engine requires `random`, an unclustered one refuses it
  expect_error(comp_iv(cd, outcome = "zbmi", model = "lmer"),
               "needs a grouping term")
  expect_error(comp_iv(cd, outcome = "zbmi", model = "lm", random = "school"),
               "does not accept")
  # the formula interface satisfies the requirement too
  expect_s3_class(comp_iv(cd, zbmi ~ (1 | school), model = "lmer"), "comp_iv")
})

test_that("clmm refuses Hess in model_args before fitting, and clm accepts it", {
  skip_if_not_installed("ordinal")
  expect_error(
    comp_iv(cd, outcome = "grade", model = "clmm", random = "school",
            model_args = list(Hess = FALSE)),
    "Hess.*clmm.*sets itself"
  )
  # ordinal::clm() computes the Hessian whatever `Hess` says, so the fit is
  # the default one.
  plain <- comp_iv(cd, outcome = "grade", model = "clm")
  no_hess <- comp_iv(cd, outcome = "grade", model = "clm",
                     model_args = list(Hess = FALSE))
  expect_equal(stats::coef(no_hess$model), stats::coef(plain$model),
               tolerance = 0)
  expect_null(stats::getCall(plain$model)$Hess)
})


# The core families' effect scale and the composition omnibus test (moved here
# here with their helpers' definitions). The registry
# test above covers the extra families; these cover lm/glm/cox/lmer/glmer.

test_that("tc_report_meta: lm reports a mean difference", {
  expect_false(tc_report_meta("lm", NULL)$exponentiate)
  expect_identical(tc_report_meta("lm", NULL)$measure, "Mean difference")
})

test_that("tc_report_meta: glm follows the family rule", {
  expect_true(tc_report_meta("glm", "binomial")$exponentiate)
  expect_identical(tc_report_meta("glm", "binomial")$measure, "Odds ratio")
  expect_true(tc_report_meta("glm", "poisson")$exponentiate)
  expect_identical(tc_report_meta("glm", "poisson")$measure, "Rate ratio")
  expect_false(tc_report_meta("glm", "gaussian")$exponentiate)
  expect_identical(tc_report_meta("glm", "gaussian")$measure, "Mean difference")
})

test_that("tc_report_meta: a GLM effect follows the link, not the family name", {
  log_binom <- tc_report_meta("glm", stats::binomial(link = "log"))
  expect_true(log_binom$exponentiate)
  expect_identical(log_binom$measure, "Risk ratio")
  expect_identical(tc_report_meta("glm", stats::quasipoisson())$measure,
                   "Rate ratio")
  log_gauss <- tc_report_meta("glm", stats::gaussian(link = "log"))
  expect_true(log_gauss$exponentiate)
  expect_identical(log_gauss$measure, "Ratio of means")
  expect_identical(tc_report_meta("glm", stats::Gamma(link = "log"))$measure,
                   "Ratio of means")
  expect_false(tc_report_meta("glm", stats::quasibinomial(link = "identity"))$exponentiate)
  # A family object and its name give the same answer.
  expect_identical(tc_report_meta("glm", stats::poisson()),
                   tc_report_meta("glm", "poisson"))
  expect_error(tc_report_meta("glm", stats::binomial(link = "probit")),
               "Unsupported link")
})

test_that("tc_report_meta: the label names a link only when it is not the default", {
  expect_identical(tc_report_meta("glm", stats::gaussian(link = "log"))$label,
                   "generalized linear model (gaussian, log link)")
  expect_identical(tc_report_meta("glm", stats::binomial())$label,
                   "generalized linear model (binomial)")
  expect_identical(tc_report_meta("glm", "binomial")$label,
                   "generalized linear model (binomial)")
  expect_identical(tc_report_meta("glmer", stats::poisson(link = "identity"))$label,
                   "generalized linear mixed model (poisson, identity link)")
})

test_that("tc_report_meta: a GLM engine without a family answers only the label", {
  r <- tc_report_meta("glm")
  expect_identical(r$label, "generalized linear model")
  expect_identical(r$exponentiate, NA)
  expect_identical(r$measure, NA_character_)
})

test_that("tc_report_meta: cox reports a hazard ratio", {
  expect_true(tc_report_meta("cox", NULL)$exponentiate)
  expect_identical(tc_report_meta("cox", NULL)$measure, "Hazard ratio")
})

test_that("tc_report_meta: lmer/glmer mirror lm/glm", {
  expect_false(tc_report_meta("lmer", NULL)$exponentiate)
  expect_identical(tc_report_meta("lmer", NULL)$measure, "Mean difference")
  expect_true(tc_report_meta("glmer", "binomial")$exponentiate)
  expect_identical(tc_report_meta("glmer", "binomial")$measure, "Odds ratio")
})

# The omnibus is assembled in comp_iv() from tc_ilr_coef_names() and
# tc_joint_ilr_test(), so it is checked on the stored result.
test_that("comp_iv stores the composition omnibus as a joint Wald test", {
  m <- comp_iv(cd, outcome = "zbmi", model = "lm")
  expect_identical(m$omnibus$term, "composition")
  expect_equal(m$omnibus$df, length(parts) - 1)
})

# comp_iv() spells its `model` choices out literally so that ?comp_iv and
# argument completion show them. This keeps that list from drifting away from
# the registry: a key added to one but not the other fails here.
test_that("comp_iv()'s model choices are exactly the registry keys, in order", {
  expect_identical(eval(formals(comp_iv)$model), names(tc_model_registry))
})
