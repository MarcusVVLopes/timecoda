parts <- c("MVPA", "LPA", "SB", "SPT")
cd <- declare_composition(timeuse, parts, closure = 1440)

test_that("comp_iv returns a comp_iv/timecoda_fit with the right engine", {
  m <- comp_iv(cd, "obesity", model = "glm", family = "binomial")
  expect_s3_class(m, "comp_iv")
  expect_s3_class(m, "timecoda_fit")
  expect_identical(m$model_type, "glm")
  expect_s3_class(m$family, "family")
  expect_identical(m$family$family, "binomial")
  expect_identical(m$family$link, "logit")
  expect_true(all(paste0("ilr", 1:3) %in% names(m$data)))
  expect_s3_class(m$comp, "compD")
})

test_that("model argument selects each engine", {
  expect_identical(comp_iv(cd, "zbmi", model = "lm")$model_type, "lm")
  expect_identical(
    comp_iv(cd, time = "time", event = "event", model = "cox")$model_type,
    "cox"
  )
  lme <- comp_iv(cd, "zbmi", model = "lmer", random = "school")
  expect_identical(lme$model_type, "lmer")
  expect_s4_class(lme$model, "lmerMod")
})

test_that("omnibus is attached", {
  m <- comp_iv(cd, "obesity", model = "glm", family = "binomial")
  expect_s3_class(m$omnibus, "tbl_df")
  expect_setequal(names(m$omnibus), c("term", "statistic", "df", "p.value"))
  expect_equal(m$omnibus$df, length(parts) - 1)
  # The surface is not precomputed at fit time; plot.comp_iv() builds it.
  expect_null(m$pred_grid)
})

test_that("number of ILR coordinates tracks the number of parts", {
  for (d in 3:5) {
    p <- paste0("p", seq_len(d))
    sim <- as.data.frame(matrix(runif(60 * d, 1, 10), ncol = d))
    names(sim) <- p
    sim$y <- rnorm(60)
    m <- comp_iv(declare_composition(sim, p), "y", model = "lm")
    expect_true(all(paste0("ilr", seq_len(d - 1)) %in% names(m$data)))
    expect_equal(m$omnibus$df, d - 1)
    expect_length(m$models, d)
  }
})

test_that("comp_iv stores one rotated model per part", {
  m <- comp_iv(cd, "obesity", model = "glm", family = "binomial")
  expect_type(m$models, "list")
  expect_length(m$models, length(parts))
  expect_named(m$models, parts)
  expect_true(all(vapply(m$models, inherits, logical(1), "glm")))
  # The canonical fit is reused as the first part's rotated model.
  expect_identical(m$models[[parts[1]]], m$model)
})

test_that("the canonical model is unchanged by storing rotated models", {
  m <- comp_iv(cd, "obesity", model = "glm", family = "binomial")
  m2 <- comp_iv(cd, "obesity", model = "glm", family = "binomial")
  expect_equal(stats::coef(m$model), stats::coef(m2$model))
  expect_identical(m$pivot, 1L)
})

test_that("each rotated model's ilr1 equals a hand-rotated refit", {
  m <- comp_iv(cd, "obesity", model = "glm", family = "binomial")
  df <- dplyr::bind_cols(tc_comp_df(cd), timeuse["obesity"])
  for (k in seq_along(parts)) {
    dk <- ilr_transform(df, parts, pivot = k)
    fk <- stats::glm(obesity ~ ilr1 + ilr2 + ilr3, data = dk,
                     family = stats::binomial())
    expect_equal(
      unname(stats::coef(m$models[[parts[k]]])["ilr1"]),
      unname(stats::coef(fk)["ilr1"]),
      tolerance = 1e-8
    )
  }
})

test_that("rotated fits are exact when the first rows repeat one composition", {
  d <- timeuse
  d[1:6, parts] <- d[rep(1, 6), parts]
  cd_rep <- declare_composition(d, parts, closure = 1440)
  comps <- tc_comp_df(cd_rep)[parts]
  # The first D + 2 rows span fewer than D - 1 coordinate directions.
  expect_lt(qr(tc_pivot_coords(comps[1:6, ], pivot = 1L))$rank,
            length(parts) - 1)
  fit <- comp_iv(cd_rep, "zbmi", model = "lm")
  df <- dplyr::bind_cols(comps, d["zbmi"])
  for (k in seq_along(parts)) {
    hand <- stats::lm(zbmi ~ ilr1 + ilr2 + ilr3,
                      data = ilr_transform(df, parts, pivot = k))
    expect_equal(stats::coef(fit$models[[parts[k]]]), stats::coef(hand),
                 tolerance = 1e-10)
  }
})

test_that("the omnibus test is invariant to the pivot basis", {
  m <- comp_iv(cd, "obesity", model = "glm", family = "binomial")
  coords <- paste0("ilr", seq_len(length(parts) - 1))
  for (k in seq_along(parts)) {
    om_k <- tc_joint_ilr_test(m$models[[parts[k]]], "glm", coords)
    expect_equal(om_k$statistic, m$omnibus$statistic, tolerance = 1e-6)
    expect_equal(om_k$p.value, m$omnibus$p.value, tolerance = 1e-6)
  }
})

test_that("comp_iv accepts multiple random effects", {
  tu2 <- timeuse
  set.seed(1)
  tu2$region <- factor(sample(paste0("r", 1:5), nrow(tu2), replace = TRUE))
  cd2 <- declare_composition(tu2, parts, closure = 1440)

  m <- comp_iv(cd2, "zbmi", model = "lmer",
               random = c("school", "region"))
  expect_identical(m$model_type, "lmer")
  expect_setequal(names(lme4::ranef(m$model)), c("school", "region"))
  expect_length(m$models, length(parts))

  # A verbatim lme4 term works too.
  m2 <- comp_iv(cd2, "zbmi", model = "lmer",
                random = "(1 | school)")
  expect_setequal(names(lme4::ranef(m2$model)), "school")
})

test_that("predict_realloc still works on a comp_iv object", {
  m <- comp_iv(cd, "obesity", model = "glm", family = "binomial")
  out <- predict_realloc(m, delta = 30, from = "SB", to = "MVPA")
  expect_s3_class(out, "tbl_df")
  expect_true(is.finite(out$estimate))
})

test_that("comp_iv rejects bad input informatively", {
  expect_error(comp_iv(timeuse, "obesity", model = "glm"), "compD")
  expect_error(comp_iv(cd, model = "lm"), "outcome")
  expect_error(comp_iv(cd, "zbmi", model = "glmer", random = "school"), "family")
})

test_that("a survival engine refuses an outcome beside time and event", {
  expect_error(
    comp_iv(cd, outcome = "zbmi", time = "time", event = "event", model = "cox"),
    "outcome.* does not apply to .*model = \"cox\""
  )
  # Positionally, as the second argument.
  expect_error(
    comp_iv(cd, "zbmi", time = "time", event = "event", model = "cox"),
    "outcome.* does not apply"
  )
})

test_that("missing survival columns are named together, without an index wrapper", {
  one <- tryCatch(comp_iv(cd, time = "nope", event = "event", model = "cox"),
                  error = identity)
  expect_match(conditionMessage(one),
               "`time` must name a single column of the .*compD.* data")
  two <- tryCatch(comp_iv(cd, time = "nope", event = "nope2", model = "cox"),
                  error = identity)
  expect_match(conditionMessage(two),
               "`time` and `event` must each name a single column of the .*compD.* data")
  expect_no_match(conditionMessage(one), "In index")
  expect_no_match(conditionMessage(two), "In index")
  expect_error(comp_iv(cd, time = "time", event = c("event", "event"), model = "cox"),
               "`event` must name a single column")
})

test_that("comp_iv refuses a GLM link whose effect has no ratio or difference reading", {
  expect_error(
    comp_iv(cd, "obesity", model = "glm", family = stats::binomial(link = "probit")),
    "probit.*not supported"
  )
  expect_error(
    comp_iv(cd, "obesity", model = "glmer", random = "school",
            family = stats::binomial(link = "cloglog")),
    "cloglog.*not supported"
  )
  # Gamma's default link is the inverse.
  d <- timeuse
  d$y_pos <- exp(d$zbmi)
  cd_pos <- declare_composition(d, parts, closure = 1440)
  expect_error(comp_iv(cd_pos, "y_pos", model = "glm", family = "Gamma"),
               "inverse.*not supported")
})

test_that("formula interface matches the discrete-argument interface", {
  fa <- comp_iv(cd, outcome = "obesity", covariates = "sex", model = "glm",
                family = "binomial")
  ff <- comp_iv(cd, obesity ~ sex, model = "glm", family = "binomial")
  expect_equal(ff$omnibus$statistic, fa$omnibus$statistic, tolerance = 1e-10)
  expect_equal(summary(ff)$per_part$estimate,
               summary(fa)$per_part$estimate, tolerance = 1e-8)
  expect_equal(
    predict_realloc(ff, from = "SB", to = "MVPA", delta = c(5, 15))$estimate,
    predict_realloc(fa, from = "SB", to = "MVPA", delta = c(5, 15))$estimate,
    tolerance = 1e-8
  )
})

test_that("omitting ilr() adds the composition as a main effect", {
  f_imp <- comp_iv(cd, obesity ~ sex, model = "glm", family = "binomial")
  f_exp <- comp_iv(cd, obesity ~ ilr() + sex, model = "glm", family = "binomial")
  expect_equal(f_imp$omnibus$statistic, f_exp$omnibus$statistic, tolerance = 1e-10)
})

test_that("cox and lmer formulas match their argument forms", {
  ca <- comp_iv(cd, time = "time", event = "event", model = "cox")
  cf <- comp_iv(cd, survival::Surv(time, event) ~ 1, model = "cox")
  expect_equal(cf$omnibus$statistic, ca$omnibus$statistic, tolerance = 1e-8)

  la <- comp_iv(cd, outcome = "zbmi", model = "lmer", random = "school")
  lf <- comp_iv(cd, zbmi ~ (1 | school), model = "lmer")
  expect_equal(lf$omnibus$statistic, la$omnibus$statistic, tolerance = 1e-8)
})

test_that("advanced interaction and non-linear terms fit", {
  m <- comp_iv(cd, obesity ~ ilr() * sex + poly(zbmi, 2), model = "glm",
               family = "binomial")
  expect_s3_class(m, "comp_iv")
  expect_true(m$has_ilr_interaction)
})

test_that("moderator equals the explicit interaction formula and adds a Wald test", {
  m1 <- comp_iv(cd, outcome = "obesity", moderator = "sex", model = "glm",
                family = "binomial")
  m2 <- comp_iv(cd, obesity ~ ilr() * sex, model = "glm", family = "binomial")
  expect_equal(m1$omnibus$statistic, m2$omnibus$statistic, tolerance = 1e-10)
  expect_equal(m1$moderator, "sex")
  expect_s3_class(m1$moderation, "tbl_df")
  expect_equal(nrow(m1$moderation), 1L)
  # df = (L-1)*(D-1) = 1 * 3 for a 2-level moderator, 4-part composition
  expect_equal(m1$moderation$df, 3L)
  expect_true(m1$has_ilr_interaction)
})

test_that("formula/moderator argument conflicts are rejected", {
  expect_error(
    comp_iv(cd, formula = obesity ~ sex, outcome = "obesity", model = "glm",
            family = "binomial"),
    "formula"
  )
  expect_error(
    comp_iv(cd, formula = obesity ~ sex, moderator = "sex", model = "glm",
            family = "binomial"),
    "formula"
  )
  expect_error(
    comp_iv(cd, outcome = "obesity", moderator = "sex", covariates = "sex",
            model = "glm", family = "binomial"),
    "covariates"
  )
  expect_error(
    comp_iv(cd, zbmi ~ sex + (1 | school), model = "glm"),
    "lmer"
  )
})

# --- curvature: the full quadratic written for the discrete arguments ---------

fixed_coefs <- function(m) {
  if (inherits(m$model, "merMod")) lme4::fixef(m$model) else stats::coef(m$model)
}

test_that("curvature = TRUE equals the explicit full-quadratic formula", {
  pairs <- list(
    lm = list(
      comp_iv(cd, "zbmi", model = "lm", curvature = TRUE),
      comp_iv(cd, formula = zbmi ~ ilr()^2 + I(ilr()^2), model = "lm")
    ),
    lm_cov = list(
      comp_iv(cd, "zbmi", model = "lm", covariates = "sex", curvature = TRUE),
      comp_iv(cd, formula = zbmi ~ ilr()^2 + I(ilr()^2) + sex, model = "lm")
    ),
    glm = list(
      comp_iv(cd, "obesity", model = "glm", family = "binomial",
              curvature = TRUE),
      comp_iv(cd, formula = obesity ~ ilr()^2 + I(ilr()^2), model = "glm",
              family = "binomial")
    ),
    lmer = list(
      comp_iv(cd, "zbmi", model = "lmer", random = "school", curvature = TRUE),
      comp_iv(cd, formula = zbmi ~ ilr()^2 + I(ilr()^2) + (1 | school),
              model = "lmer")
    )
  )
  for (nm in names(pairs)) {
    sugar <- pairs[[nm]][[1]]
    explicit <- pairs[[nm]][[2]]
    # 3 main effects, 3 products and 3 squares for a 4-part composition
    expect_length(grep("ilr", names(fixed_coefs(sugar))), 9L)
    expect_equal(fixed_coefs(sugar), fixed_coefs(explicit), tolerance = 1e-12,
                 label = nm)
    expect_equal(sugar$omnibus, explicit$omnibus, tolerance = 1e-12, label = nm)
    expect_equal(sugar$omnibus$df, 9L, label = nm)
    expect_equal(
      predict_realloc(sugar, from = "SB", to = "MVPA", delta = c(10, 30))$estimate,
      predict_realloc(explicit, from = "SB", to = "MVPA", delta = c(10, 30))$estimate,
      tolerance = 1e-12, label = nm
    )
    expect_null(sugar$moderator)
  }
})

test_that("curvature with a moderator keeps the moderation on the main effects", {
  sugar <- comp_iv(cd, "zbmi", model = "lm", moderator = "sex",
                   curvature = TRUE)
  explicit <- comp_iv(cd, formula = zbmi ~ ilr()^2 + I(ilr()^2) + ilr() * sex,
                      model = "lm")
  expect_equal(stats::coef(sugar$model), stats::coef(explicit$model),
               tolerance = 1e-12)
  expect_equal(sugar$moderator, "sex")
  expect_s3_class(sugar$moderation, "tbl_df")
  # (L-1)(D-1) = 1 * 3: the coordinate-by-sex terms only
  expect_equal(sugar$moderation$df, 3L)
  auto <- predict_realloc(sugar, from = "SB", to = "MVPA", delta = 15)
  expect_true("sex" %in% names(auto))
  expect_equal(nrow(auto), nlevels(factor(timeuse$sex)))
})

test_that("curvature is refused where it cannot apply", {
  expect_error(comp_iv(cd, "zbmi", model = "lm", curvature = "yes"),
               "curvature.*must be")
  expect_error(comp_iv(cd, "zbmi", model = "lm", curvature = NA),
               "curvature.*must be")
  expect_error(
    comp_iv(cd, formula = zbmi ~ ilr() + sex, model = "lm", curvature = TRUE),
    "curvature.*discrete arguments only"
  )
  expect_error(
    comp_iv(cd, zbmi ~ ilr() + sex, model = "lm", curvature = TRUE),
    "curvature.*discrete arguments only"
  )
  expect_error(
    comp_iv(cd, "zbmi", model = "lm", covariates = "ilr1:sex",
            curvature = TRUE),
    "cannot also name a coordinate.*ilr1"
  )
  expect_error(
    comp_iv(cd, "zbmi", model = "lm", covariates = "ilr():sex",
            curvature = TRUE),
    "marker.*ilr\\(\\)"
  )
})

test_that("curvature combines with include_total", {
  set.seed(7)
  d <- timeuse
  d[parts] <- d[parts] * stats::rlnorm(nrow(d), 0, 0.2)
  cd_vt <- declare_composition(d, parts, closure = "mean")
  sugar <- comp_iv(cd_vt, "zbmi", model = "lm", curvature = TRUE,
                   include_total = TRUE)
  explicit <- comp_iv(cd_vt, formula = zbmi ~ ilr()^2 + I(ilr()^2),
                      model = "lm", include_total = TRUE)
  expect_true("t_total" %in% names(stats::coef(sugar$model)))
  expect_equal(stats::coef(sugar$model), stats::coef(explicit$model),
               tolerance = 1e-12)
  expect_equal(sugar$total_test, explicit$total_test, tolerance = 1e-12)
  expect_error(
    comp_iv(cd_vt, "zbmi", model = "lm", covariates = "t_total",
            curvature = TRUE, include_total = TRUE),
    "t_total"
  )
})

test_that("a full quadratic stores the joint test of its squares and products", {
  sugar <- comp_iv(cd, "zbmi", model = "lm", covariates = "sex",
                   curvature = TRUE)
  explicit <- comp_iv(cd, formula = zbmi ~ ilr()^2 + I(ilr()^2) + sex,
                      model = "lm")
  expect_s3_class(sugar$curvature_test, "tbl_df")
  expect_identical(sugar$curvature_test$term, "curvature")
  # 3 squares and 3 products for a 4-part composition
  expect_equal(sugar$curvature_test$df, 6L)
  expect_equal(sugar$curvature_test, explicit$curvature_test, tolerance = 1e-12)
  # The hand Wald statistic b' V^-1 b on the six quadratic coefficients.
  b <- stats::coef(sugar$model)
  quad <- grep("\\^2\\)$|^ilr[0-9]:ilr[0-9]$", names(b), value = TRUE)
  expect_length(quad, 6L)
  v <- stats::vcov(sugar$model)[quad, quad]
  hand <- as.numeric(t(b[quad]) %*% solve(v) %*% b[quad])
  expect_equal(sugar$curvature_test$statistic, hand, tolerance = 1e-8)
  expect_equal(sugar$curvature_test$p.value,
               stats::pchisq(hand, df = 6, lower.tail = FALSE),
               tolerance = 1e-8)
  # The same model in every ILR basis: the parts given in another order.
  cd_rev <- declare_composition(timeuse, rev(parts), closure = 1440)
  reordered <- comp_iv(cd_rev, "zbmi", model = "lm", covariates = "sex",
                       curvature = TRUE)
  expect_equal(reordered$curvature_test$statistic,
               sugar$curvature_test$statistic, tolerance = 1e-6)
  # With a moderator, the test covers the terms in the coordinates alone.
  moderated <- comp_iv(cd, "zbmi", model = "lm", moderator = "sex",
                       curvature = TRUE)
  expect_equal(moderated$curvature_test$df, 6L)
  moderated_formula <- comp_iv(
    cd, formula = zbmi ~ ilr()^2 + I(ilr()^2) + ilr() * sex, model = "lm"
  )
  expect_equal(moderated$curvature_test, moderated_formula$curvature_test,
               tolerance = 1e-12)
})

test_that("only a full quadratic carries a curvature test", {
  expect_null(comp_iv(cd, "zbmi", model = "lm")$curvature_test)
  # Squares only, products only, and a polynomial basis.
  expect_null(
    comp_iv(cd, formula = zbmi ~ ilr() + I(ilr()^2), model = "lm")$curvature_test
  )
  expect_null(comp_iv(cd, formula = zbmi ~ ilr()^2, model = "lm")$curvature_test)
  expect_null(
    comp_iv(cd, formula = zbmi ~ poly(ilr(), 2), model = "lm")$curvature_test
  )
  # A cubic term is a term of another degree.
  expect_null(
    comp_iv(cd, formula = zbmi ~ ilr()^2 + I(ilr()^2) + I(ilr1^3),
            model = "lm")$curvature_test
  )
})

test_that("a two-part composition tests its one square", {
  sim <- sim_parts(2)
  m <- comp_iv(declare_composition(sim, c("p1", "p2")), "y", model = "lm",
               curvature = TRUE)
  expect_equal(m$curvature_test$df, 1L)
  st <- summary(m$model)$coefficients["I(ilr1^2)", "t value"]
  expect_equal(m$curvature_test$statistic, st^2, tolerance = 1e-8)
})

test_that("a multinomial full quadratic is tested jointly over its categories", {
  skip_if_not_installed("mclogit")
  d <- timeuse
  d$grade_unordered <- factor(as.character(d$grade), levels = levels(d$grade))
  m <- suppressWarnings(comp_iv(declare_composition(d, parts, closure = 1440),
                                "grade_unordered", model = "mblogit",
                                curvature = TRUE))
  k <- nlevels(d$grade_unordered)
  expect_equal(m$curvature_test$df, 6L * (k - 1L))
})

test_that("a malformed covariate gets one message with or without include_total", {
  set.seed(7)
  d <- timeuse
  d[parts] <- d[parts] * stats::rlnorm(nrow(d), 0, 0.2)
  cd_vt <- declare_composition(d, parts, closure = "mean")
  plain <- expect_error(comp_iv(cd_vt, "zbmi", covariates = "sex +"),
                        "not found")
  total <- expect_error(comp_iv(cd_vt, "zbmi", covariates = "sex +",
                                include_total = TRUE),
                        "not found")
  expect_identical(conditionMessage(total), conditionMessage(plain))
})

test_that("the messages that give an example put it in parentheses", {
  flat <- function(err) gsub("\\s+", " ", conditionMessage(err))
  expect_match(
    flat(expect_error(comp_iv(cd, "zbmi", model = "lm", model_args = list(1)))),
    "Name them as the backend does (e.g. `model_args = list(weights = w)`).",
    fixed = TRUE
  )
  expect_match(
    flat(expect_error(comp_iv(cd, "zbmi", model = "lmer"))),
    "Supply `random` (e.g. `random = \"school\"`) or a `(1 | group)` term",
    fixed = TRUE
  )
  skip_if_not_installed("ordinal")
  expect_match(
    flat(expect_error(comp_iv(cd, "zbmi", model = "clm"))),
    "(e.g. `d$zbmi <- ordered(d$zbmi)`).",
    fixed = TRUE
  )
})

# --- model_args: backend arguments, captured unevaluated ---------------------

test_that("model_args reaches the backend and is captured unevaluated", {
  set.seed(1)
  d <- timeuse
  d$wt <- runif(nrow(d), 0.5, 1.5)
  cdw <- declare_composition(d, parts, closure = 1440)

  # A BARE column name must work, exactly as in a direct lm() call: the
  # expression travels unevaluated and the column is carried into the model
  # frame for the backend to resolve.
  mw <- comp_iv(cdw, "zbmi", model = "lm", model_args = list(weights = wt))
  m0 <- comp_iv(cdw, "zbmi", model = "lm")
  expect_false(isTRUE(all.equal(stats::coef(mw$model), stats::coef(m0$model))))
  # every rotated refit is weighted too, not just the canonical fit
  expect_true(all(purrr::map_lgl(mw$models, ~ !is.null(.x$weights))))

  # capturing unevaluated keeps the stored call's argument as typed, and the
  # call replays with the same weights
  expect_identical(mw$model$call$weights, quote(wt))
  expect_equal(stats::coef(stats::update(mw$model)), stats::coef(mw$model),
               tolerance = 0)
})

test_that("model_args passes scalar backend settings through", {
  m <- comp_iv(cd, "obesity", model = "glm", family = "binomial",
               model_args = list(epsilon = 1e-10))
  expect_s3_class(m, "comp_iv")
  expect_equal(m$model$control$epsilon, 1e-10)
})

test_that("model_args rejects malformed input with its own message", {
  expect_error(comp_iv(cd, "zbmi", model = "lm", model_args = list(1, 2)),
               "must be named")
  expect_error(comp_iv(cd, "zbmi", model = "lm", model_args = 42),
               "must be a list")
})

test_that("model_args accepts a pre-built list", {
  ma <- list(epsilon = 1e-10)
  m <- comp_iv(cd, "obesity", model = "glm", family = "binomial",
               model_args = ma)
  expect_equal(m$model$control$epsilon, 1e-10)
})

test_that("a formula-path fit without random terms stores random = NULL", {
  cd <- declare_composition(timeuse, parts, closure = 1440)
  m <- comp_iv(cd, zbmi ~ ilr() + sex, model = "lm")
  expect_null(m$random)
})

test_that("model_args resolves an object local to the function that calls comp_iv()", {
  cd <- declare_composition(timeuse, parts, closure = 1440)
  fit_with_local_control <- function() {
    ctrl <- stats::glm.control(maxit = 50)
    comp_iv(cd, "obesity", model = "glm", family = "binomial",
            model_args = list(control = ctrl))
  }
  m <- fit_with_local_control()
  expect_s3_class(m, "comp_iv")
  plain <- comp_iv(cd, "obesity", model = "glm", family = "binomial")
  expect_equal(stats::coef(m$model), stats::coef(plain$model))
})

test_that("model_args resolves weights held by the function that calls comp_iv()", {
  cd <- declare_composition(timeuse, parts, closure = 1440)
  fit_with_local_weights <- function() {
    w <- rep(2, nrow(timeuse))
    comp_iv(cd, "zbmi", model = "lm", model_args = list(weights = w))
  }
  m <- fit_with_local_weights()
  expect_s3_class(m, "comp_iv")
  expect_equal(unname(stats::weights(m$model)), rep(2, nrow(timeuse)))
  plain <- comp_iv(cd, "zbmi", model = "lm")
  expect_equal(stats::coef(m$model), stats::coef(plain$model))
})

test_that("a discrete-path covariate must be a column or a scalar of the caller's", {
  expect_error(comp_iv(cd, "zbmi", covariates = "nonexistent"),
               "Column .*nonexistent.* not found in the .*compD.* data")
  # A vector, local or global, is not row-aligned with the composition.
  w_local <- stats::rnorm(nrow(timeuse))
  expect_error(comp_iv(cd, "zbmi", covariates = "w_local"),
               "Column .*w_local.* not found")
  assign("tc_test_w_global", stats::rnorm(nrow(timeuse)), envir = globalenv())
  expect_error(comp_iv(cd, "zbmi", covariates = "tc_test_w_global"),
               "Column .*tc_test_w_global.* not found")
  rm("tc_test_w_global", envir = globalenv())
  # A scalar is replaced by its value, giving the literal's fit.
  kn <- 16
  scalar <- comp_iv(cd, "zbmi", covariates = "splines::ns(time, knots = kn)")
  literal <- comp_iv(cd, "zbmi", covariates = "splines::ns(time, knots = 16)")
  expect_identical(scalar$covariates, literal$covariates)
  expect_equal(stats::coef(scalar$model), stats::coef(literal$model))
})

test_that("a caller-local covariate scalar reaches moderator and curvature fits", {
  fit_local <- function(...) {
    k <- 3
    comp_iv(cd, "zbmi", model = "lm", covariates = "I(time * k)", ...)
  }
  moderated <- fit_local(moderator = "sex")
  curved <- fit_local(curvature = TRUE)
  expect_match(deparse1(stats::formula(moderated$model)), "I(time * 3)",
               fixed = TRUE)
  expect_match(deparse1(stats::formula(curved$model)), "I(time * 3)",
               fixed = TRUE)
  # A vector of the caller's, and an unknown name, are refused on both paths.
  vector_local <- function(...) {
    w_local <- stats::rnorm(nrow(timeuse))
    comp_iv(cd, "zbmi", model = "lm", covariates = "I(w_local * 2)", ...)
  }
  expect_error(vector_local(moderator = "sex"), "Column .*w_local.* not found")
  expect_error(vector_local(curvature = TRUE), "Column .*w_local.* not found")
  expect_error(comp_iv(cd, "zbmi", model = "lm", moderator = "sex",
                       covariates = "nope"),
               "Column .*nope.* not found")
  expect_error(comp_iv(cd, "zbmi", model = "lm", curvature = TRUE,
                       covariates = "I(nope * 2)"),
               "Column .*nope.* not found")
})

test_that("a covariate scalar named like a function is replaced as a value", {
  c <- 3
  discrete <- comp_iv(cd, "zbmi", model = "lm", covariates = "I(time * c)")
  expect_identical(discrete$covariates, "I(time * 3)")
  moderated <- comp_iv(cd, "zbmi", model = "lm", moderator = "sex",
                       covariates = "I(time * c)")
  expect_match(deparse1(stats::formula(moderated$model)), "I(time * 3)",
               fixed = TRUE)
  # The call c(0) keeps its function name beside the replaced value.
  both <- comp_iv(cd, "zbmi", model = "lm", covariates = "I(time * c + c(0))")
  expect_identical(both$covariates, "I(time * 3 + c(0))")
  expect_equal(stats::coef(both$model)[["I(time * 3 + c(0))"]],
               stats::coef(discrete$model)[["I(time * 3)"]])
})

test_that("a caller-local model_args object reaches a moderator fit", {
  fit_local <- function() {
    w <- rep(2, nrow(timeuse))
    comp_iv(cd, "zbmi", model = "lm", moderator = "sex",
            model_args = list(weights = w))
  }
  d <- timeuse
  d$w2 <- 2
  cd_w <- declare_composition(d, parts, closure = 1440)
  column <- comp_iv(cd_w, "zbmi", model = "lm", moderator = "sex",
                    model_args = list(weights = w2))
  expect_equal(stats::coef(fit_local()$model), stats::coef(column$model))
})

test_that("comp_iv() refuses a composition part named as a model variable", {
  refused <- "must not name a composition part"
  expect_error(comp_iv(cd, "zbmi", covariates = "SB"), refused)
  expect_error(comp_iv(cd, "zbmi", covariates = "SB:sex"), refused)
  expect_error(comp_iv(cd, "zbmi", moderator = "SB"), refused)
  expect_error(comp_iv(cd, formula = zbmi ~ ilr() + SB), refused)
  expect_error(comp_iv(cd, formula = zbmi ~ ilr() * SB), refused)
  expect_error(comp_iv(cd, "zbmi", model_args = list(weights = SB)), refused)
  expect_error(comp_iv(cd, "zbmi", covariates = "SB"), "\"SB\" is a part")
  expect_error(comp_iv(cd, "zbmi", covariates = c("SB", "LPA")),
               "\"SB\" and \"LPA\" are parts")
})

test_that("include_total must be TRUE or FALSE", {
  for (bad in list(1, "TRUE", NA, c(TRUE, FALSE))) {
    expect_error(comp_iv(cd, "zbmi", model = "lm", include_total = bad),
                 "include_total.*must be")
  }
})

test_that("a malformed covariate gets one message with or without curvature", {
  plain <- expect_error(comp_iv(cd, "zbmi", covariates = "sex +"),
                        "not found")
  curved <- expect_error(comp_iv(cd, "zbmi", covariates = "sex +",
                                 curvature = TRUE),
                         "not found")
  expect_identical(conditionMessage(curved), conditionMessage(plain))
  # Only the term that does not parse is named.
  expect_error(comp_iv(cd, "zbmi", covariates = c("region", "sex +"),
                       curvature = TRUE),
               "Column \"sex \\+\" not found")
})

test_that("a malformed covariate gets one message with or without a moderator", {
  plain <- expect_error(comp_iv(cd, "zbmi", covariates = "region +"),
                        "not found")
  moderated <- expect_error(comp_iv(cd, "zbmi", covariates = "region +",
                                    moderator = "sex"),
                            "not found")
  expect_identical(conditionMessage(moderated), conditionMessage(plain))
})
