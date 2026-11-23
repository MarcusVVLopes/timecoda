# The optional engines -- ordinal (clm/clmm), multinomial (mblogit), robust
# (lmrob/rlmer) and mixed-effects Cox (coxme) -- fitted through comp_iv() and
# the shared tc_iv_fit() builder (R/fit_models.R), from their tc_model_registry
# rows. Checks the resulting model class and that the composition coefficients
# ilr1..ilr{D-1} are present. Each backend is optional (Suggests), so every
# test skips when its package is unavailable. test-model_backends.R covers the
# registry itself and end-to-end summarize/reallocate behaviour; these focus
# narrowly on the engine fit, except the coxme block at the end, which also
# covers that engine's summary, reallocation and diagnostics against hand fits.

parts <- c("MVPA", "LPA", "SB", "SPT")
coords <- paste0("ilr", seq_len(length(parts) - 1))

# timeuse's `grade` is an ordered factor (composition-driven latent score),
# and `region` deterministically nests `school` -- used here as the grouping
# factor for the clustered engines.
cd <- declare_composition(timeuse, parts, closure = 1440)

test_that("ordinal: clm fits with the composition coefficients on ilr1..ilr{D-1}", {
  skip_if_not_installed("ordinal")
  m <- comp_iv(cd, outcome = "grade", model = "clm")
  expect_identical(m$model_type, "clm")
  expect_s3_class(m$model, "clm")
  expect_true(all(coords %in% names(stats::coef(m$model))))
})

test_that("ordinal: clmm fits clustered via random = \"school\"", {
  skip_if_not_installed("ordinal")
  m <- comp_iv(cd, outcome = "grade", model = "clmm", random = "school")
  expect_identical(m$model_type, "clmm")
  expect_s3_class(m$model, "clmm")
  expect_true(all(coords %in% names(stats::coef(m$model))))
})

test_that("multinom: mblogit fits with one coefficient block per non-reference category", {
  skip_if_not_installed("mclogit")
  # multinom takes an unordered factor outcome; grade coerced accordingly.
  d <- timeuse
  d$grade_unordered <- factor(as.character(d$grade), levels = levels(d$grade))
  cd2 <- declare_composition(d, parts, closure = 1440)
  m <- suppressWarnings(comp_iv(cd2, outcome = "grade_unordered", model = "mblogit"))
  expect_identical(m$model_type, "mblogit")
  b_names <- names(stats::coef(m$model))
  predictor <- sub("^[^~]*~", "", b_names)
  expect_true(all(coords %in% predictor))
})

test_that("rlm: lmrob fits with the composition coefficients on ilr1..ilr{D-1}", {
  skip_if_not_installed("robustbase")
  m <- comp_iv(cd, outcome = "zbmi", model = "lmrob")
  expect_identical(m$model_type, "lmrob")
  expect_s3_class(m$model, "lmrob")
  expect_true(all(coords %in% names(stats::coef(m$model))))
})

test_that("clustered mblogit takes its grouping variable from either spec form", {
  skip_if_not_installed("mclogit")
  d <- timeuse
  d$grade_unordered <- factor(as.character(d$grade), levels = levels(d$grade))
  cd2 <- declare_composition(d, parts, closure = 1440)
  by_name <- suppressWarnings(
    comp_iv(cd2, outcome = "grade_unordered", model = "mblogit", random = "school")
  )
  by_term <- suppressWarnings(
    comp_iv(cd2, outcome = "grade_unordered", model = "mblogit",
            random = "(1 | school)")
  )
  expect_identical(deparse(stats::getCall(by_name$model)$random), "~1 | school")
  expect_equal(stats::coef(by_name$model), stats::coef(by_term$model))
})

test_that("mblogit refuses more than one grouping term", {
  skip_if_not_installed("mclogit")
  d <- timeuse
  d$grade_unordered <- factor(as.character(d$grade), levels = levels(d$grade))
  cd2 <- declare_composition(d, parts, closure = 1440)
  expect_error(
    comp_iv(cd2, outcome = "grade_unordered", model = "mblogit",
            random = c("school", "region")),
    "mblogit.*accepts one grouping term"
  )
})

test_that("a clustered mblogit given a formula with a bar term fits like the discrete form", {
  skip_if_not_installed("mclogit")
  d <- timeuse
  d$grade_unordered <- factor(as.character(d$grade), levels = levels(d$grade))
  cd2 <- declare_composition(d, parts, closure = 1440)
  by_formula <- suppressWarnings(
    comp_iv(cd2, grade_unordered ~ ilr() + (1 | school), model = "mblogit")
  )
  by_arg <- suppressWarnings(
    comp_iv(cd2, outcome = "grade_unordered", model = "mblogit", random = "school")
  )
  expect_equal(stats::coef(by_formula$model), stats::coef(by_arg$model))
  expect_identical(deparse(stats::getCall(by_formula$model)$random), "~1 | school")
})

test_that("a clustered mblogit keeps a bar term's random slope", {
  skip_if_not_installed("mclogit")
  d <- timeuse
  d$grade_unordered <- factor(as.character(d$grade), levels = levels(d$grade))
  cd2 <- declare_composition(d, parts, closure = 1440)
  # Every clustered mblogit fit of these data warns that mclogit's inner
  # iterations did not converge, hence suppressWarnings() on each fit.
  slope <- suppressWarnings(
    comp_iv(cd2, outcome = "grade_unordered", model = "mblogit",
            random = "(ilr1 | school)")
  )
  expect_identical(deparse(stats::getCall(slope$model)$random), "~ilr1 | school")
  expect_true(slope$convergence$converged)
  hand <- suppressWarnings(mclogit::mblogit(
    grade_unordered ~ ilr1 + ilr2 + ilr3, data = slope$data,
    random = ~ ilr1 | school, trace = FALSE
  ))
  expect_equal(stats::coef(slope$model), stats::coef(hand), tolerance = 1e-10)
  # The formula form fits the same model as the discrete one.
  by_formula <- suppressWarnings(
    comp_iv(cd2, grade_unordered ~ ilr() + (ilr1 | school), model = "mblogit")
  )
  expect_equal(stats::coef(by_formula$model), stats::coef(slope$model),
               tolerance = 1e-8)
  # The slope is fitted: one covariance matrix over both categories'
  # intercepts and ilr1 slopes, and fixed effects that differ from the
  # intercept-only fit's.
  intercept <- suppressWarnings(
    comp_iv(cd2, outcome = "grade_unordered", model = "mblogit",
            random = "school")
  )
  expect_false(isTRUE(all.equal(stats::coef(slope$model),
                                stats::coef(intercept$model))))
  vc <- slope$model$VarCov
  expect_identical(names(vc), "school")
  expect_identical(rownames(vc$school),
                   c("expected~1", "above~1", "expected~ilr1", "above~ilr1"))
  expect_identical(dim(intercept$model$VarCov$school), c(2L, 2L))
  all_slopes <- suppressWarnings(
    comp_iv(cd2, grade_unordered ~ ilr() + (ilr() | school), model = "mblogit")
  )
  expect_match(deparse(stats::getCall(all_slopes$model)$random),
               "ilr1 + ilr2 + ilr3", fixed = TRUE)
})

test_that("a clustered mblogit refuses a misspelt group, a || term and two terms", {
  skip_if_not_installed("mclogit")
  d <- timeuse
  d$grade_unordered <- factor(as.character(d$grade), levels = levels(d$grade))
  cd2 <- declare_composition(d, parts, closure = 1440)
  expect_error(
    comp_iv(cd2, outcome = "grade_unordered", model = "mblogit", random = "schol"),
    "not found in the"
  )
  expect_error(
    comp_iv(cd2, outcome = "grade_unordered", model = "mblogit",
            random = "(ilr1 || school)"),
    "correlated random effects only"
  )
  expect_error(
    comp_iv(cd2, grade_unordered ~ ilr() + (ilr() || school), model = "mblogit"),
    "correlated random effects only"
  )
  expect_error(
    comp_iv(cd2, outcome = "grade_unordered", model = "mblogit",
            random = "(1 | school) + (1 | region)"),
    "one grouping term"
  )
})

test_that("mblogit refuses a numeric outcome", {
  skip_if_not_installed("mclogit")
  expect_error(comp_iv(cd, "zbmi", model = "mblogit"), "is numeric")
  expect_error(comp_iv(cd, zbmi ~ ilr(), model = "mblogit"), "is numeric")
})

test_that("model_args may not name an argument the engine sets itself", {
  skip_if_not_installed("mclogit")
  d <- timeuse
  d$grade_unordered <- factor(as.character(d$grade), levels = levels(d$grade))
  cd2 <- declare_composition(d, parts, closure = 1440)
  expect_error(
    comp_iv(cd2, "grade_unordered", model = "mblogit",
            model_args = list(trace = TRUE)),
    "sets itself"
  )
  expect_error(
    comp_iv(cd, time = "time", event = "event", model = "cox",
            model_args = list(model = FALSE)),
    "sets itself"
  )
})

test_that("the ordinal engines take the logit link only", {
  skip_if_not_installed("ordinal")
  expect_error(
    comp_iv(cd, "grade", model = "clm", model_args = list(link = "probit")),
    "logit link"
  )
  m <- comp_iv(cd, "grade", model = "clm", model_args = list(link = "logit"))
  expect_identical(m$model$link, "logit")
})

test_that("clmm and rlmer fits report their convergence", {
  skip_if_not_installed("ordinal")
  cm <- comp_iv(cd, "grade", model = "clmm", random = "school")
  expect_true(cm$convergence$converged)
  expect_output(print(summary(cm)), "Convergence: OK")
  skip_if_not_installed("robustlmm")
  rl <- comp_iv(cd, "zbmi", model = "rlmer", random = "school")
  expect_true(rl$convergence$converged)
})

test_that("an lmrob fit that did not converge is reported as such", {
  skip_if_not_installed("robustbase")
  expect_error(
    suppressWarnings(comp_iv(
      cd, "zbmi", model = "lmrob",
      model_args = list(control = robustbase::lmrob.control(max.it = 1))
    )),
    "did not converge"
  )
})

test_that("the clm, clmm, mblogit, lmrob and rlmer registry rows keep their fields", {
  keys <- c("clm", "clmm", "mblogit", "lmrob", "rlmer")
  field <- function(f) {
    purrr::map(rlang::set_names(keys), function(k) tc_model_meta(k)[[f]])
  }
  expect_identical(
    purrr::map_chr(field("fitter"), deparse),
    c(clm = "ordinal::clm", clmm = "ordinal::clmm",
      mblogit = "mclogit::mblogit", lmrob = "robustbase::lmrob",
      rlmer = "robustlmm::rlmer")
  )
  expect_identical(field("fixed_args"), list(
    clm = list(), clmm = list(Hess = TRUE), mblogit = list(trace = FALSE),
    lmrob = list(), rlmer = list()
  ))
  expect_identical(field("outcome_check"), list(
    clm = "ordered", clmm = "ordered", mblogit = "unordered",
    lmrob = NULL, rlmer = NULL
  ))
  expect_identical(field("emm_mode"), list(
    clm = "latent", clmm = "latent", mblogit = NULL, lmrob = NULL, rlmer = NULL
  ))
  expect_identical(field("link_tran"), list(
    clm = "log", clmm = "log", mblogit = "log", lmrob = NULL, rlmer = NULL
  ))
  chr <- function(f) unlist(field(f))
  expect_identical(chr("scale"), c(clm = "odds", clmm = "odds", mblogit = "rrr",
                                   lmrob = "difference", rlmer = "difference"))
  expect_identical(chr("random_style"), c(clm = "none", clmm = "bars",
                                          mblogit = "arg", lmrob = "none",
                                          rlmer = "bars"))
  expect_identical(chr("convergence"), c(clm = "ordinal", clmm = "ordinal",
                                         mblogit = "flag", lmrob = "flag",
                                         rlmer = "optinfo"))
  expect_identical(chr("check_backend"), c(clm = "ordinal", clmm = "limited",
                                           mblogit = "limited",
                                           lmrob = "robust",
                                           rlmer = "robust_mixed"))
  expect_identical(chr("emm_via"), c(clm = "ref_grid", clmm = "ref_grid",
                                     mblogit = "ref_grid", lmrob = "qdrg",
                                     rlmer = "ref_grid"))
  expect_identical(chr("icc"), c(clm = "none", clmm = "clmm", mblogit = "none",
                                 lmrob = "none", rlmer = "lmm"))
  expect_identical(field("emm_df"), list(
    clm = NULL, clmm = NULL, mblogit = NULL, lmrob = NULL, rlmer = Inf
  ))
  expect_identical(chr("per_category"), c(clm = FALSE, clmm = FALSE,
                                          mblogit = TRUE, lmrob = FALSE,
                                          rlmer = FALSE))
  expect_identical(chr("absolute_ok"), c(clm = FALSE, clmm = FALSE,
                                         mblogit = FALSE, lmrob = TRUE,
                                         rlmer = TRUE))
})

# Mixed-effects Cox (coxme) ----------------------------------------------------
# The clustered counterpart of `cox`. No other test file fits it, so this
# block covers its consumers too, each against a hand coxme() fit on the pivot
# coordinates. The fit is built once, on first use, after the skip.

coxme_fit <- local({
  fit <- NULL
  function() {
    if (is.null(fit)) {
      fit <<- comp_iv(cd, time = "time", event = "event", model = "coxme",
                      random = "school")
    }
    fit
  }
})

# A hand coxme() fit of pivot `k`'s coordinates, outside comp_iv().
coxme_by_hand <- function(k = 1L) {
  d <- dplyr::bind_cols(tc_comp_df(cd), timeuse[c("time", "event", "school")])
  coxme::coxme(
    survival::Surv(time, event) ~ ilr1 + ilr2 + ilr3 + (1 | school),
    data = ilr_transform(d, parts, pivot = k)
  )
}

test_that("the coxme registry row describes a relative-only clustered survival engine", {
  m <- tc_model_meta("coxme")
  expect_identical(m$fitter, quote(coxme::coxme))
  expect_identical(m$response, "surv")
  expect_identical(m$random_style, "bars")
  expect_true(m$requires_random)
  expect_true(m$fixef)
  expect_false(m$absolute_ok)
  expect_identical(m$pkg, "coxme")
  expect_identical(m$check_backend, "limited")
  expect_identical(tc_report_meta("coxme")$measure, "Hazard ratio")
})

test_that("coxme fits Surv(time, event) ~ ilr1 + ... + (1 | school) and matches a hand fit", {
  skip_if_not_installed("coxme")
  m <- coxme_fit()
  expect_identical(m$model_type, "coxme")
  expect_s3_class(m$model, "coxme")
  expect_identical(deparse(m$model$formulaList$fixed),
                   "survival::Surv(time, event) ~ ilr1 + ilr2 + ilr3")
  expect_identical(deparse(m$model$formulaList$random[[1]]), "(1 | school)")
  # The stored call holds the formula and the fitting frame as values.
  expect_s3_class(stats::getCall(m$model)$formula, "formula")
  expect_identical(deparse(stats::getCall(m$model)$formula),
                   "survival::Surv(time, event) ~ ilr1 + ilr2 + ilr3 + (1 | school)")
  expect_identical(stats::getCall(m$model)$data, m$data)

  hand <- coxme_by_hand()
  cf <- tc_model_coefs(m$model, "coxme")
  expect_equal(cf$b, coxme::fixef(hand), tolerance = 1e-8)
  expect_equal(cf$V, as.matrix(stats::vcov(hand)), tolerance = 1e-8)
  # The omnibus is the Wald test of the three fixed ILR effects.
  b <- coxme::fixef(hand)
  wald <- as.numeric(t(b) %*% solve(as.matrix(stats::vcov(hand))) %*% b)
  expect_equal(m$omnibus$statistic, wald, tolerance = 1e-8)
  expect_identical(m$omnibus$df, 3L)
})

test_that("coxme summaries report hazard ratios equal to exp(coef) of each rotated fit", {
  skip_if_not_installed("coxme")
  m <- coxme_fit()
  s <- summary(m)
  expect_identical(s$effect_measure, "Hazard ratio")
  expect_output(print(s), "Hazard ratio")
  expect_output(print(m), "mixed-effects Cox proportional hazards model")
  # One rotated coxme per part, each on its own pivot's coordinates.
  expect_identical(names(m$models), parts)
  expect_true(all(purrr::map_lgl(m$models, inherits, "coxme")))
  hand_ilr1 <- purrr::map_dbl(seq_along(parts), function(k) {
    coxme::fixef(coxme_by_hand(k))[["ilr1"]]
  })
  expect_equal(unname(purrr::map_dbl(m$models, function(x) coxme::fixef(x)[["ilr1"]])),
               hand_ilr1, tolerance = 1e-8)
  expect_equal(s$per_part$estimate, exp(hand_ilr1), tolerance = 1e-8)
  skip_if_not_installed("gt")
  expect_s3_class(pub_ready(s)$effects, "gt_tbl")
})

test_that("coxme relative reallocations equal exp(d'beta), from any calling frame", {
  skip_if_not_installed("coxme")
  m <- coxme_fit()
  top <- predict_realloc(m, delta = 15, from = "SB", to = "MVPA")
  # Fitted and predicted from fresh function frames: emmeans recovers the data
  # through the fit, not through the stored call.
  fit_here <- function() {
    comp_iv(declare_composition(timeuse, parts, closure = 1440),
            time = "time", event = "event", model = "coxme", random = "school")
  }
  predict_here <- function(f) predict_realloc(f, delta = 15, from = "SB", to = "MVPA")
  expect_equal(as.data.frame(predict_here(fit_here())), as.data.frame(top),
               tolerance = 0)

  ref <- unlist(m$ref_comp)
  moved <- ref + c(MVPA = 15, LPA = 0, SB = -15, SPT = 0)[parts]
  z <- function(x) unlist(ilr_transform(as.data.frame(t(x)), parts)[coords])
  dz <- z(moved) - z(ref)
  hand <- coxme_by_hand()
  expect_equal(top$estimate, exp(sum(dz * coxme::fixef(hand))), tolerance = 1e-8)
  expect_identical(top$df, Inf)

  # The predicted-outcome surface is the same contrast at each grid point.
  pg <- tc_surface_panel(m, c("MVPA", "LPA", "SB"), ref[parts], 6,
                         "relative", "Hazard ratio")$data
  grid1 <- c(MVPA = pg$x[1], LPA = pg$y[1], SB = pg$z[1], SPT = ref[["SPT"]])
  expect_equal(pg$.pred[1],
               exp(sum((z(grid1) - z(ref)) * coxme::fixef(hand))),
               tolerance = 1e-8)
})

test_that("coxme refuses absolute predictions, surfaces and searches", {
  skip_if_not_installed("coxme")
  m <- coxme_fit()
  expect_error(
    predict_realloc(m, delta = 15, from = "SB", to = "MVPA", type = "absolute"),
    "support `type = \"relative\"` only"
  )
  expect_error(plot(m, scale = "absolute"), "support `scale = \"relative\"` only")
  expect_error(optimal_composition(m, type = "absolute"),
               "support .type = \"relative\". only")
})

test_that("fit_check on a coxme fit reports fit statistics and a note", {
  skip_if_not_installed("coxme")
  m <- coxme_fit()
  fc <- fit_check(m, plot = FALSE)
  expect_identical(fc$backend, "limited")
  expect_identical(fc$tests$check, c("AIC", "log-likelihood"))
  expect_true(is.na(fc$convergence$converged))
  expect_equal(fc$tests$statistic[2], as.numeric(stats::logLik(m$model)))
  expect_identical(fc$note, tc_model_meta("coxme")$limited_note)
  expect_null(fc$raw)
})

test_that("a coxme formula with a bar term fits like the discrete form", {
  skip_if_not_installed("coxme")
  m <- coxme_fit()
  f <- comp_iv(cd, survival::Surv(time, event) ~ ilr() + (1 | school),
               model = "coxme")
  expect_equal(coxme::fixef(f$model), coxme::fixef(m$model), tolerance = 0)
  expect_equal(f$omnibus, m$omnibus, tolerance = 0)
  expect_equal(
    as.data.frame(predict_realloc(f, delta = 15, from = "SB", to = "MVPA")),
    as.data.frame(predict_realloc(m, delta = 15, from = "SB", to = "MVPA")),
    tolerance = 0
  )
})

test_that("coxme needs a grouping term", {
  skip_if_not_installed("coxme")
  expect_error(
    comp_iv(cd, time = "time", event = "event", model = "coxme"),
    "clustered engine and needs a grouping term"
  )
  expect_error(
    comp_iv(cd, survival::Surv(time, event) ~ ilr(), model = "coxme"),
    "clustered engine and needs a grouping term"
  )
})

test_that("coxme refuses an outcome beside time and event", {
  skip_if_not_installed("coxme")
  expect_error(
    comp_iv(cd, outcome = "zbmi", time = "time", event = "event",
            model = "coxme", random = "school"),
    "outcome.* does not apply to .*model = \"coxme\""
  )
})
