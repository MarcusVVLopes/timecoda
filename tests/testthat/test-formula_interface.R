parts <- c("MVPA", "LPA", "SB", "SPT")
cd <- declare_composition(timeuse, parts, closure = 1440)

test_that("ilr() errors if called directly", {
  expect_error(ilr(), "only valid inside")
})

test_that("the parser expands ilr() to the coordinate sum", {
  spec <- tc_iv_parse_formula(zbmi ~ ilr() + sex, cd, "lm")
  txt <- paste(deparse(spec$model_formula), collapse = " ")
  expect_match(txt, "ilr1 \\+ ilr2 \\+ ilr3")
  expect_match(txt, "sex")
})

test_that("the parser auto-adds the composition when ilr() is absent", {
  spec <- tc_iv_parse_formula(zbmi ~ sex + region, cd, "lm")
  txt <- paste(deparse(spec$model_formula), collapse = " ")
  expect_match(txt, "ilr1 \\+ ilr2 \\+ ilr3")
  expect_identical(spec$covariates, c("sex", "region"))
})

test_that("the parser distributes ilr() across interactions", {
  spec <- tc_iv_parse_formula(zbmi ~ ilr() * sex, cd, "lm")
  expect_true(all(c("ilr1:sex", "ilr3:sex") %in% spec$covariates))
})

test_that("ilr() is a sum inside formula operators and distributed inside other calls", {
  for (k in c(3L, 5L)) {
    coords <- tc_coord_names(k + 1L)
    cs <- str2lang(paste0("(", paste(coords, collapse = " + "), ")"))
    each <- function(fmt) paste0("(", paste(sprintf(fmt, coords), collapse = " + "), ")")
    expand <- function(e) deparse1(tc_expand_ilr_expr(e, cs))
    sum_txt <- deparse1(cs)
    expect_identical(expand(quote(ilr()^2)), paste0(sum_txt, "^2"))
    expect_identical(expand(quote(ilr() * sex)), paste(sum_txt, "* sex"))
    expect_identical(expand(quote(ilr() %in% g)), paste(sum_txt, "%in% g"))
    expect_identical(expand(quote((ilr() | g))), paste0("(", sum_txt, " | g)"))
    expect_identical(expand(quote(I(ilr()^2))), each("I(%s^2)"))
    expect_identical(expand(quote(poly(ilr(), 2))), each("poly(%s, 2)"))
    expect_identical(expand(quote(splines::ns(ilr(), 3))), each("splines::ns(%s, 3)"))
    # Every ilr() in one call takes the same coordinate: squares, not products.
    expect_identical(expand(quote(I(ilr() * ilr()))), each("I(%1$s * %1$s)"))
    expect_error(tc_expand_ilr_expr(quote(I(ilr(x)^2)), cs), "takes no arguments")
  }
})

test_that("an lme4 random slope on the composition is split out as a random term", {
  spec <- tc_iv_parse_formula(zbmi ~ ilr() + (ilr() | school), cd, "lmer")
  expect_identical(spec$random, "((ilr1 + ilr2 + ilr3) | school)")
  expect_error(comp_iv(cd, zbmi ~ ilr() + (ilr() | school), model = "lm"),
               "does not accept random-effect terms")
})

test_that("a formula without the composition's main effect warns once and fits", {
  n_warn <- function(f) {
    n <- 0L
    withCallingHandlers(comp_iv(cd, f, model = "lm"), warning = function(w) {
      n <<- n + 1L
      invokeRestart("muffleWarning")
    })
    n
  }
  expect_identical(n_warn(zbmi ~ I(ilr()^2)), 1L)
  expect_warning(comp_iv(cd, zbmi ~ I(ilr()^2), model = "lm"),
                 "does not include the composition's main effect")
  expect_identical(n_warn(zbmi ~ poly(ilr(), 2)), 0L)
  expect_identical(n_warn(zbmi ~ ilr() + I(ilr()^2)), 0L)
  expect_identical(n_warn(zbmi ~ sex), 0L)
  # With no term in the coordinates alone, the omnibus tests every coefficient
  # crossing a coordinate with a covariate.
  fit <- suppressWarnings(comp_iv(cd, zbmi ~ ilr():sex, model = "lm"))
  nm <- names(stats::coef(fit$model))[-1]
  b <- stats::coef(fit$model)[nm]
  hand <- as.numeric(t(b) %*% solve(stats::vcov(fit$model)[nm, nm]) %*% b)
  expect_equal(c(fit$omnibus$statistic, fit$omnibus$df), c(hand, length(nm)),
               tolerance = 1e-10)
  region <- suppressWarnings(comp_iv(cd, zbmi ~ ilr():region, model = "lm"))
  expect_identical(region$omnibus$df,
                   (length(parts) - 1L) * nlevels(factor(timeuse$region)))
  # A call on the covariate side names its levels the same way.
  called <- suppressWarnings(comp_iv(cd, zbmi ~ ilr():factor(sex), model = "lm"))
  expect_equal(called$omnibus, fit$omnibus, tolerance = 1e-10)
})

test_that("non-linear composition terms fit, and the omnibus tests all of them", {
  d1 <- length(parts) - 1L
  coords <- tc_coord_names(length(parts))
  squares <- sprintf("I(%s^2)", coords)
  products <- combn(coords, 2, paste, collapse = ":")
  q <- comp_iv(cd, zbmi ~ ilr()^2 + I(ilr()^2), model = "lm")
  expect_setequal(names(stats::coef(q$model)), c("(Intercept)", coords, squares, products))
  expect_equal(q$omnibus$df, d1 * (d1 + 3L) / 2L)
  expect_false(q$has_ilr_interaction)
  s <- comp_iv(cd, zbmi ~ ilr() + I(ilr()^2), model = "lm")
  expect_equal(s$omnibus$df, 2L * d1)
  # A repeated linear term is aliased and left out: the same column space as s.
  aliased <- comp_iv(cd, zbmi ~ ilr() + poly(ilr(), 2), model = "lm")
  expect_equal(aliased$omnibus, s$omnibus, tolerance = 1e-8)
  p <- comp_iv(cd, obesity ~ poly(ilr(), 2), model = "glm", family = "binomial")
  expect_identical(names(stats::coef(p$model))[-1],
                   c(t(outer(sprintf("poly(%s, 2)", coords), 1:2, paste0))))
  expect_equal(p$omnibus$df, 2L * d1)
  expect_output(print(summary(q)), "slopes at the reference composition")
  # A reallocation needs no `at`, and equals the hand contrast.
  ref <- tc_check_reference(NULL, q)[parts]
  moved <- ref + c(10, 0, -10, 0)
  nd <- as.data.frame(tc_pivot_coords(as.data.frame(rbind(ref, moved)), 1L))
  hand <- diff(stats::predict(q$model, newdata = nd))
  est <- predict_realloc(q, delta = 10, from = "SB", to = "MVPA")$estimate
  expect_equal(est, unname(hand), tolerance = 1e-10)
  # The full quadratic does not depend on the part order; squares alone do.
  cd2 <- declare_composition(timeuse, c("SB", "SPT", "MVPA", "LPA"), closure = 1440)
  realloc <- function(f, comp) {
    predict_realloc(comp_iv(comp, f, model = "lm"), delta = 10, from = "SB",
                    to = "MVPA")$estimate
  }
  expect_equal(realloc(zbmi ~ ilr()^2 + I(ilr()^2), cd2), est, tolerance = 1e-10)
  expect_false(isTRUE(all.equal(realloc(zbmi ~ ilr() + I(ilr()^2), cd2),
                                realloc(zbmi ~ ilr() + I(ilr()^2), cd))))
  # Every coordinate is pinned on the surface grid: no emmeans interaction note.
  panel <- function(fit, res, scale) {
    tc_surface_panel(fit, parts[1:3], ref, res, scale, "value")
  }
  expect_no_message(panel(q, 4, "relative"))
  expect_no_message(panel(q, 4, "absolute"))
  # A composition-by-covariate interaction is averaged over on the surface;
  # emmeans' note stays off there too (plot.comp_iv() says so once itself).
  mod <- comp_iv(cd, "zbmi", moderator = "sex", model = "lm")
  expect_no_message(panel(mod, 3, "relative"))
})

test_that("only a term crossing a coordinate with a covariate is an interaction", {
  crossed <- comp_iv(cd, zbmi ~ ilr() + I(ilr1^2):sex, model = "lm")
  expect_true(crossed$has_ilr_interaction)
  expect_false(comp_iv(cd, zbmi ~ ilr() + ilr1:ilr2, model = "lm")$has_ilr_interaction)
  ns_fit <- comp_iv(cd, zbmi ~ splines::ns(ilr(), 3), model = "lm")
  expect_false(ns_fit$has_ilr_interaction)
  expect_s3_class(predict_realloc(ns_fit, delta = 10, from = "SB", to = "MVPA"),
                  "realloc")
})

test_that("tc_ilr_coef_names returns the coefficients a set of term labels generates", {
  cn <- c("ilr1", "ilr2", "ilr3")
  b <- c("(Intercept)", "ilr1", "ilr2", "ilr3", "sexmale", "ilr1:sexmale",
         "I(ilr1^2)", "ilr1:ilr2", "poly(ilr2, 2)1", "poly(ilr2, 2)2", "ilr10")
  expect_identical(tc_ilr_coef_names(b, "lm", cn), cn)
  expect_identical(tc_ilr_coef_names(b, "lm", c(cn, "I(ilr1^2)", "ilr1:ilr2")),
                   c(cn, "I(ilr1^2)", "ilr1:ilr2"))
  expect_identical(tc_ilr_coef_names(b, "lm", "poly(ilr2, 2)"),
                   c("poly(ilr2, 2)1", "poly(ilr2, 2)2"))
  expect_identical(tc_ilr_coef_names(b, "lm", "ilr1:sexmale"), "ilr1:sexmale")
  mb <- c("B~(Intercept)", "B~ilr1", "B~ilr2", "B~ilr3", "C~ilr1", "C~ilr2",
          "C~ilr3", "B~I(ilr1^2)")
  expect_identical(tc_ilr_coef_names(mb, "mblogit", c(cn, "I(ilr1^2)")), mb[-1])
  expect_identical(tc_ilr_coef_names(mb, "mblogit", "ilr1"), c("B~ilr1", "C~ilr1"))
})

test_that("every name a formula or moderator uses must be a column or a scalar", {
  expect_error(comp_iv(cd, zbmi ~ ilr() + nope, model = "lm"),
               "Column \"nope\" not found in the .*compD.* data")
  expect_error(comp_iv(cd, zbmi ~ ilr() + nope + (1 | nope2), model = "lmer"),
               "Columns \"nope\" and \"nope2\" not found")
  expect_error(comp_iv(cd, outcome = "zbmi", model = "lmer", moderator = "sex",
                       random = "nope"),
               "Column \"nope\" not found")
  expect_error(comp_iv(cd, zbmi ~ ilr() + c, model = "lm"),
               "Column \"c\" not found")
  # `T` and a caller's knot value are scalars, not columns.
  kn <- stats::median(timeuse$time)
  expect_s3_class(
    comp_iv(cd, zbmi ~ ilr() + poly(time, 2, raw = T) + I(time > kn), model = "lm"),
    "comp_iv"
  )
  # A vector of the caller's must be a column instead.
  env_vec <- stats::rnorm(nrow(timeuse))
  expect_error(comp_iv(cd, zbmi ~ ilr() + env_vec, model = "lm"),
               "Column \"env_vec\" not found")
})

test_that("an environment vector covariate is refused at fit time", {
  grp <- factor(rep(c("a", "b", "c"), length.out = nrow(timeuse)))
  expect_error(comp_iv(cd, zbmi ~ ilr() * grp, model = "lm"),
               "Column \"grp\" not found")
})

test_that("a scalar from the formula's environment is stored as its value", {
  kn <- 16
  fit <- comp_iv(cd, zbmi ~ ilr() + I(time > kn), model = "lm")
  expect_match(deparse1(stats::formula(fit$model)), "I(time > 16)", fixed = TRUE)
  expect_s3_class(predict_realloc(fit, delta = 15, from = "SB", to = "MVPA"),
                  "realloc")
  # A scalar named like a function stays a name, so the function call survives.
  c <- 2
  expect_s3_class(
    comp_iv(cd, zbmi ~ ilr() + poly(time, c) + I(c(1) * as.numeric(region)),
            model = "lm"),
    "comp_iv"
  )
})

test_that("a saved moderator fit is about the size of its formula twin", {
  mod <- comp_iv(cd, outcome = "zbmi", model = "lm", moderator = "sex")
  # The twin's formula is built in the global environment, which serializes as
  # a reference; this test's own environment would carry `mod` and `cd`.
  twin <- comp_iv(cd, stats::as.formula("zbmi ~ ilr() * sex", env = globalenv()),
                  model = "lm")
  ratio <- length(serialize(mod, NULL)) / length(serialize(twin, NULL))
  expect_lt(abs(ratio - 1), 0.1)
})

test_that("tc_expr_has_ilr detects the marker", {
  expect_true(tc_expr_has_ilr(quote(ilr() + x)))
  expect_true(tc_expr_has_ilr(quote(a * ilr())))
  expect_false(tc_expr_has_ilr(quote(a + b)))
})

test_that("tc_moderation_test finds the composition-by-moderator interaction", {
  m <- comp_iv(cd, outcome = "zbmi", moderator = "sex", model = "lm")
  out <- tc_moderation_test(m$model, m$model_type, "sex", parts)
  expect_s3_class(out, "tbl_df")
  expect_identical(names(out), c("term", "statistic", "df", "p.value"))
  # df = (L - 1) * (D - 1) for a 2-level factor moderator.
  expect_equal(out$df, (nlevels(cd$data_raw$sex) - 1) * (length(parts) - 1))
  # comp_iv() stores the identical test on the fit.
  expect_equal(out$statistic, m$moderation$statistic, tolerance = 1e-8)
  # A coding set on the column itself is read too, and the test does not
  # depend on it.
  d <- timeuse
  stats::contrasts(d$sex) <- stats::contr.sum(2)
  m_sum <- comp_iv(declare_composition(d, parts, closure = 1440),
                   outcome = "zbmi", moderator = "sex", model = "lm")
  expect_equal(m_sum$moderation$statistic, m$moderation$statistic)
})

test_that("tc_moderation_test errors when no interaction terms are present", {
  m <- comp_iv(cd, outcome = "zbmi", model = "lm")
  expect_error(
    tc_moderation_test(m$model, m$model_type, "sex", parts),
    "Could not locate"
  )
  # A prefix of the moderator's name is not the moderator.
  m <- comp_iv(cd, outcome = "zbmi", model = "lm", moderator = "sex")
  expect_error(tc_moderation_test(m$model, m$model_type, "se", parts),
               "Could not locate")
})

test_that("a multinomial fit's moderation test reads its per-category names", {
  skip_if_not_installed("mclogit")
  m <- comp_iv(cd, outcome = "region", model = "mblogit", moderator = "sex")
  ints <- grep(":sexmale$", names(stats::coef(m$model)), value = TRUE)
  # (K - 1) * (L - 1) * (D - 1) for 3 outcome categories and a 2-level factor.
  expect_equal(m$moderation$df, 2 * 1 * 3)
  expect_equal(m$moderation$statistic,
               tc_joint_ilr_test(m$model, "mblogit", ints)$statistic,
               tolerance = 0)
})

test_that("the moderator path checks outcome, time and event as the discrete path does", {
  expect_error(
    comp_iv(cd, outcome = "zbmi", time = "time", event = "event",
            model = "cox", moderator = "sex"),
    "outcome.* does not apply to .*model = \"cox\""
  )
  expect_error(comp_iv(cd, model = "lm", moderator = "sex"),
               "outcome.* is required")
  expect_error(comp_iv(cd, model = "cox", moderator = "sex"),
               "time.* and .*event.* must each name a single column")
})

test_that("a survival engine's formula needs a Surv() left-hand side, named for that engine", {
  expect_error(comp_iv(cd, time ~ ilr(), model = "cox"),
               "model = \"cox\".*left-hand side must be a .*survival::Surv")
  skip_if_not_installed("coxme")
  expect_error(comp_iv(cd, time ~ ilr() + (1 | school), model = "coxme"),
               "model = \"coxme\".*left-hand side must be a .*survival::Surv")
})

test_that("the formula path requires a random term for every clustered engine", {
  skip_if_not_installed("ordinal")
  expect_error(comp_iv(cd, grade ~ ilr(), model = "clmm"), "random")
  skip_if_not_installed("robustlmm")
  expect_error(comp_iv(cd, zbmi ~ ilr(), model = "rlmer"), "random")
})

test_that("a formula with a random-effect term fits without a warning", {
  expect_no_warning(comp_iv(cd, zbmi ~ ilr() + (1 | school), model = "lmer"))
})

test_that("a random term longer than one deparse line parses as one term", {
  pn <- paste0("p", 1:10)
  sim <- sim_parts(10)
  sim$g <- rep(letters[1:8], length.out = nrow(sim))
  spec <- tc_iv_parse_formula(y ~ ilr() + (ilr() | g),
                              declare_composition(sim, pn), "lmer")
  expect_length(spec$random, 1L)
  expect_match(spec$random, "ilr9) | g)", fixed = TRUE)
})
