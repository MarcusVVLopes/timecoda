# Internal model machinery behind comp_iv(): the tc_iv_fit() builder (through
# the fits comp_iv() returns), tc_random_terms(), the intraclass
# correlation summary() reports,
# tc_match_family() and print.timecoda_fit(). comp_iv() itself is tested in
# test-comp_iv.R.

parts <- c("MVPA", "LPA", "SB", "SPT")
cd <- declare_composition(timeuse, parts, closure = 1440)

test_that("tc_iv_fit builds the documented timecoda_fit field contract", {
  m <- comp_iv(cd, outcome = "zbmi", model = "glm", family = "gaussian",
               covariates = "sex")
  expect_s3_class(m, "timecoda_fit")
  # The per-pivot fields tc_iv_fit() writes, in order, ahead of the
  # whole-fit fields comp_iv() adds.
  expect_identical(
    names(m)[1:12],
    c("model", "model_type", "parts", "total", "pivot", "ref_comp",
      "covariates", "has_ilr_interaction", "family", "data",
      "fit_warnings", "convergence")
  )
  expect_identical(m$model_type, "glm")
  expect_identical(m$parts, parts)
  expect_identical(m$total, 1440)
  expect_identical(m$pivot, 1L)
  expect_identical(m$ref_comp, geom_mean_comp(tc_comp_df(cd), parts, total = 1440))
  expect_identical(m$covariates, "sex")
  expect_s3_class(m$family, "family")
  expect_identical(m$family$family, "gaussian")
  expect_true(all(c(parts, "ilr1", "ilr2", "ilr3", "zbmi", "sex") %in% names(m$data)))
})

test_that("tc_iv_fit leaves family NULL off the GLM engines and detects ILR interactions", {
  plain <- comp_iv(cd, outcome = "zbmi")
  expect_null(plain$family)
  expect_false(plain$has_ilr_interaction)

  interacted <- comp_iv(cd, outcome = "zbmi", covariates = c("sex", "ilr1:sex"))
  expect_true(interacted$has_ilr_interaction)
})

test_that("comp_iv gives D-1 ILR columns and a reference composition on the closure", {
  for (d in c(3, 5)) {
    p <- paste0("p", seq_len(d))
    sim <- as.data.frame(matrix(stats::runif(40 * d, 1, 10), ncol = d))
    names(sim) <- p
    sim$y <- stats::rnorm(40)
    m <- comp_iv(declare_composition(sim, p, closure = 1), outcome = "y")
    coords <- paste0("ilr", seq_len(d - 1))
    expect_true(all(coords %in% names(m$data)))
    expect_false(paste0("ilr", d) %in% names(m$data))
    expect_named(m$ref_comp, p)
    expect_equal(sum(unlist(m$ref_comp)), 1, tolerance = 1e-6)
  }
})

test_that("each rotated fit is on its own pivot's coordinates", {
  m <- comp_iv(cd, outcome = "zbmi")
  ilr1 <- purrr::map(m$models, function(x) stats::model.frame(x)$ilr1)
  expect_false(isTRUE(all.equal(ilr1[[1]], ilr1[[2]])))
  # Pivot k's ilr1 is the balance of part k against the rest.
  z2 <- ilr_transform(tc_comp_df(cd), parts, pivot = 2)$ilr1
  expect_equal(ilr1[[2]], z2, tolerance = 1e-10)
})

# The formula is built inside tc_iv_fit() from the registry row, so it is
# checked on the fitted model rather than on a standalone helper.
test_that("tc_iv_fit builds outcome ~ ilr1 + ... + covariates", {
  cd <- declare_composition(timeuse, parts, closure = 1440)
  m <- comp_iv(cd, outcome = "zbmi", covariates = c("sex", "region"))
  expect_identical(deparse(stats::formula(m$model)),
                   "zbmi ~ ilr1 + ilr2 + ilr3 + sex + region")

  # D - 1 coordinates, derived from the parts rather than assumed.
  cd3 <- declare_composition(timeuse, parts[1:3], closure = "mean")
  m3 <- comp_iv(cd3, outcome = "zbmi")
  expect_identical(deparse(stats::formula(m3$model)), "zbmi ~ ilr1 + ilr2")

  # A survival response passes through as the LHS.
  mc <- comp_iv(cd, time = "time", event = "event", model = "cox")
  expect_identical(deparse(stats::formula(mc$model)),
                   "survival::Surv(time, event) ~ ilr1 + ilr2 + ilr3")
})

test_that("tc_random_terms expands grouping names and keeps bar terms verbatim", {
  expect_identical(tc_random_terms("school"), "(1 | school)")
  expect_identical(tc_random_terms("(1 | school)"), "(1 | school)")
  expect_identical(tc_random_terms(c("province", "school")),
                   c("(1 | province)", "(1 | school)"))
  expect_identical(tc_random_terms(c("school", "(1 | region)"), timeuse),
                   c("(1 | school)", "(1 | region)"))
  # Given data, a plain name must be one of its columns; a bar term is not
  # checked. Without data it only reformats, for the moderator sugar in
  # tc_iv_resolve_spec(), which runs before any model frame exists.
  expect_identical(tc_random_terms("(1 | not_a_column)", timeuse),
                   "(1 | not_a_column)")
  expect_identical(tc_random_terms("anything"), "(1 | anything)")
})

test_that("tc_random_terms names every missing grouping column in one message", {
  one <- tryCatch(tc_random_terms("not_a_column", timeuse), error = identity)
  expect_match(conditionMessage(one),
               "Grouping column .*not_a_column.* not found in the .*compD.* data")
  two <- tryCatch(tc_random_terms(c("nope1", "(1 | region)", "nope2"), timeuse),
                  error = identity)
  expect_match(conditionMessage(two),
               "Grouping columns .*nope1.* and .*nope2.* not found")
  # Straight from tc_random_terms(), not through an iteration wrapper.
  expect_no_match(conditionMessage(one), "In index")
  expect_no_match(conditionMessage(two), "In index")
  expect_error(comp_iv(cd, "zbmi", model = "lmer", random = "not_a_column"),
               "Grouping column .*not_a_column.* not found")
})

test_that("tc_random_terms refuses anything but a character vector without NA", {
  for (bad in list(NA_character_, character(0), 1, list("school"),
                   c("school", NA))) {
    expect_error(tc_random_terms(bad), "must be a character vector")
  }
})

test_that("summary()$icc reads the between-cluster share from the fit's own variance components", {
  fit <- comp_iv(cd, "zbmi", model = "lmer", random = "school")
  icc <- summary(fit)$icc
  expect_s3_class(icc, "tbl_df")
  expect_identical(names(icc), c("group", "variance", "measure", "value"))
  expect_identical(icc$group, "school")
  expect_identical(icc$measure, "ICC")
  v <- lme4::VarCorr(fit$model)$school[1, 1]
  expect_equal(icc$variance, v)
  expect_equal(icc$value, v / (v + stats::sigma(fit$model)^2))
  expect_equal(icc$value, 0.2638, tolerance = 1e-3)
  # One row per grouping term.
  two <- comp_iv(cd, "zbmi", model = "lmer", random = c("region", "school"))
  icc2 <- summary(two)$icc
  expect_setequal(icc2$group, c("region", "school"))
  expect_identical(icc2$measure, c("ICC", "ICC"))
  # Each row is that term's share of the total variance, so the rows cannot
  # sum past 1.
  vc2 <- lme4::VarCorr(two$model)
  total2 <- vc2$region[1, 1] + vc2$school[1, 1] + stats::sigma(two$model)^2
  expect_equal(icc2$value[icc2$group == "school"], vc2$school[1, 1] / total2)
  expect_equal(icc2$value[icc2$group == "region"], vc2$region[1, 1] / total2)
  expect_lte(sum(icc2$value), 1)
  # A binomial glmer: the latent-scale form.
  g <- comp_iv(cd, "obesity", model = "glmer", family = "binomial",
               random = "school")
  iccg <- summary(g)$icc
  expect_identical(iccg$measure, "latent ICC")
  vg <- lme4::VarCorr(g$model)$school[1, 1]
  expect_equal(iccg$value, vg / (vg + pi^2 / 3))
  expect_equal(iccg$value, 0.00139, tolerance = 1e-2)
  # print() reports one line per row.
  s <- summary(fit)
  expect_output(print(s), "Intraclass correlation \\(school\\): 0.264")
  expect_output(print(summary(g)), "Latent intraclass correlation \\(school\\): 0.001")
  out <- paste(utils::capture.output(print(summary(two))), collapse = "\n")
  expect_match(out, "Intraclass correlation \\(school\\): 0.264")
  expect_match(out, "Intraclass correlation \\(region\\): 0.000")
})

test_that("summary()$icc gives a number near zero for a singular fit, never NA", {
  set.seed(1)
  d <- timeuse
  d$y_indep <- stats::rnorm(nrow(d))
  fit <- suppressMessages(
    comp_iv(declare_composition(d, parts, closure = 1440), "y_indep",
            model = "lmer", random = "school")
  )
  expect_true(isTRUE(fit$convergence$singular))
  icc <- summary(fit)$icc
  expect_false(is.na(icc$value))
  expect_gte(icc$value, 0)
  expect_lt(icc$value, 1e-6)
})

test_that("summary()$icc is NULL without a random term and not available where no measure exists", {
  lm_fit <- comp_iv(cd, "zbmi", model = "lm")
  expect_null(summary(lm_fit)$icc)
  out <- paste(utils::capture.output(print(summary(lm_fit))), collapse = "\n")
  expect_no_match(out, "ntraclass|hazard ratio")
  d <- timeuse
  d$count <- as.integer(round(exp(pmin(pmax(d$zbmi, -3), 3))))
  pois <- comp_iv(declare_composition(d, parts, closure = 1440), "count",
                  model = "glmer", family = "poisson", random = "school")
  iccp <- summary(pois)$icc
  expect_identical(iccp$group, "school")
  expect_identical(iccp$measure, "not available")
  expect_true(is.na(iccp$value))
  expect_output(print(summary(pois)), "Intraclass correlation: not available for this model")
  # The mixed ordinal model gives the latent ICC from its own link, here from
  # a bar-term spec.
  skip_if_not_installed("ordinal")
  cm <- comp_iv(cd, "grade", model = "clmm", random = "(1 | school)")
  iccm <- summary(cm)$icc
  vm <- lme4::VarCorr(cm$model)$school[1, 1]
  expect_identical(iccm$group, "school")
  expect_identical(iccm$measure, "latent ICC")
  expect_equal(iccm$value, vm / (vm + pi^2 / 3), tolerance = 1e-12)
  expect_output(print(summary(cm)),
                "Latent intraclass correlation (school): 0.069", fixed = TRUE)
  # The robust mixed model reads its own variance components.
  skip_if_not_installed("robustlmm")
  r <- suppressWarnings(
    comp_iv(cd, "zbmi", model = "rlmer", random = "(1 | school)")
  )
  iccr <- summary(r)$icc
  v <- lme4::VarCorr(r$model)$school[1, 1]
  expect_identical(iccr$measure, "ICC")
  expect_equal(iccr$value, v / (v + stats::sigma(r$model)^2), tolerance = 1e-12)
})

test_that("summary()$icc takes the latent residual from the binomial link", {
  # comp_iv() accepts only the logit and log binomial links, so the probit and
  # cloglog cases are hand glmer fits of an outcome that does cluster.
  d <- ilr_transform(timeuse, parts, total = 1440)
  d$high <- as.integer(d$zbmi > 0)
  # Each is put in place of the model of a logit comp_iv() fit, whose summary()
  # reads the fitted model's own family.
  base <- comp_iv(cd, "obesity", model = "glmer", family = "binomial",
                  random = "school")
  hand <- function(link) {
    lme4::glmer(high ~ ilr1 + ilr2 + ilr3 + (1 | school), data = d,
                family = stats::binomial(link))
  }
  hand_icc <- function(model) {
    fit <- base
    fit$model <- model
    summary(fit)$icc
  }
  pr <- hand("probit")
  iccpr <- hand_icc(pr)
  expect_identical(iccpr$measure, "latent ICC")
  vpr <- lme4::VarCorr(pr)$school[1, 1]
  expect_gt(vpr, 0.1)
  expect_equal(iccpr$value, vpr / (vpr + 1))
  # No latent logistic or normal variable exists for any other link.
  cl <- hand("cloglog")
  iccc <- hand_icc(cl)
  expect_identical(iccc$group, "school")
  expect_identical(iccc$measure, "not available")
  expect_true(is.na(iccc$value))
  # The log link, which comp_iv() does accept, when it converges.
  lg <- tryCatch(
    suppressWarnings(suppressMessages(
      comp_iv(cd, "obesity", model = "glmer", family = stats::binomial("log"),
              random = "school")
    )),
    error = function(e) NULL
  )
  if (!is.null(lg)) {
    expect_identical(summary(lg)$icc$measure, "not available")
  }
})

test_that("summary()$icc reports the median hazard ratio for a coxme fit", {
  skip_if_not_installed("coxme")
  fit <- comp_iv(cd, time = "time", event = "event", model = "coxme",
                 random = "school")
  icc <- summary(fit)$icc
  expect_identical(icc$group, "school")
  expect_identical(icc$measure, "median hazard ratio")
  v <- unname(coxme::VarCorr(fit$model)$school[1])
  expect_equal(icc$variance, v)
  expect_equal(icc$value, exp(sqrt(2 * v) * stats::qnorm(0.75)))
  expect_equal(icc$value, 1.2311, tolerance = 1e-4)
  expect_output(print(summary(fit)), "Median hazard ratio \\(school\\): 1.231")
})

test_that("tc_match_family resolves strings, functions, and family objects", {
  expect_s3_class(tc_match_family("binomial"), "family")
  expect_identical(tc_match_family("binomial")$family, "binomial")
  expect_s3_class(tc_match_family(stats::gaussian), "family")
  expect_s3_class(tc_match_family(stats::poisson()), "family")
})

test_that("tc_match_family errors informatively on NULL or invalid input", {
  expect_error(tc_match_family(NULL), "required")
  expect_error(tc_match_family("not_a_family_fn"), "family")
  expect_error(tc_match_family(42), "family")
})

test_that("print.timecoda_fit prints the model class, parts, and reference composition", {
  out <- comp_iv(cd, outcome = "zbmi")
  expect_output(print(out), "timecoda_fit")
  expect_output(print(out), "linear model")
  expect_output(print(out), "MVPA")
  expect_output(print(out), "Reference composition")

  # print.timecoda_fit returns its argument invisibly.
  expect_identical(withr::with_output_sink(nullfile(), print(out)), out)
})

test_that("print.timecoda_fit appends the GLM family to the label", {
  out <- comp_iv(cd, outcome = "obesity", model = "glm", family = "binomial")
  expect_output(print(out), "generalized linear model \\(binomial\\)")
})

test_that("a declared unit label follows the total in print() and summary()", {
  cdu <- declare_composition(timeuse, parts, closure = 1440, units = "min/day")
  m <- comp_iv(cdu, outcome = "zbmi")
  expect_output(print(m), "(total = 1440 min/day)", fixed = TRUE)
  expect_identical(summary(m)$units, "min/day")
  expect_output(print(summary(m)), "(total = 1440 min/day)", fixed = TRUE)
  plain <- comp_iv(cd, outcome = "zbmi")
  expect_null(summary(plain)$units)
  expect_output(print(plain), "(total = 1440)", fixed = TRUE)
})

test_that("tc_joint_ilr_test recovers the omnibus composition Wald test", {
  m <- comp_iv(cd, outcome = "zbmi", model = "lm")
  b_names <- names(tc_model_coefs(m$model, "lm")$b)
  coords <- intersect(paste0("ilr", 1:3), b_names)
  out <- tc_joint_ilr_test(m$model, "lm", coords)
  expect_s3_class(out, "tbl_df")
  expect_identical(names(out), c("term", "statistic", "df", "p.value"))
  expect_identical(out$term, "composition")
  expect_equal(out$df, length(coords))
  expect_true(out$statistic >= 0)
  expect_true(out$p.value >= 0 && out$p.value <= 1)
  # Matches the omnibus stored by comp_iv() itself.
  expect_equal(out$statistic, m$omnibus$statistic, tolerance = 1e-8)
})

test_that("tc_joint_ilr_test errors when coefficients cannot be located", {
  m <- comp_iv(cd, outcome = "zbmi", model = "lm")
  expect_error(
    tc_joint_ilr_test(m$model, "lm", c("not_a_coef")),
    "Could not locate"
  )
})

test_that("the stored call carries the GLM family as it was typed", {
  by_name <- comp_iv(cd, outcome = "obesity", model = "glm", family = "binomial")
  expect_identical(stats::getCall(by_name$model)$family, "binomial")
  by_call <- comp_iv(cd, outcome = "obesity", model = "glm",
                     family = binomial(link = "logit"))
  expect_identical(deparse(stats::getCall(by_call$model)$family),
                   "binomial(link = \"logit\")")
  # The fit is the one the resolved family gives.
  direct <- stats::glm(obesity ~ ilr1 + ilr2 + ilr3, data = by_name$data,
                       family = stats::binomial())
  expect_equal(stats::coef(by_name$model), stats::coef(direct),
               tolerance = 1e-12)
  expect_identical(by_name$family$link, "logit")
  expect_equal(stats::coef(by_call$model), stats::coef(direct),
               tolerance = 1e-12)
})

test_that("a family typed as a local name resolves where comp_iv() was called", {
  fit_with_local_family <- function(model, ...) {
    my_family <- function() stats::binomial(link = "logit")
    comp_iv(cd, outcome = "obesity", model = model, family = "my_family", ...)
  }
  m <- fit_with_local_family("glm")
  ref <- comp_iv(cd, outcome = "obesity", model = "glm", family = "binomial")
  expect_equal(stats::coef(m$model), stats::coef(ref$model), tolerance = 0)
  expect_identical(stats::getCall(m$model)$family, "my_family")
  # lme4::glmer() looks a family name up two frames up, so it is handed the
  # resolved family object rather than the name.
  mm <- fit_with_local_family("glmer", random = "school")
  expect_identical(stats::getCall(mm$model)$family, "my_family")
})

test_that("the stored call carries the formula and the fitting frame as values", {
  m <- comp_iv(cd, outcome = "zbmi", covariates = "sex")
  f <- stats::getCall(m$model)$formula
  expect_s3_class(f, "formula")
  expect_match(paste(deparse(f), collapse = " "), "ilr1 \\+ ilr2 \\+ ilr3 \\+ sex")
  expect_identical(stats::getCall(m$model)$data, m$data)
  cox <- comp_iv(cd, time = "time", event = "event", model = "cox")
  expect_match(paste(deparse(stats::getCall(cox$model)$formula), collapse = " "),
               "survival::Surv(time, event) ~ ilr1 + ilr2 + ilr3", fixed = TRUE)
  expect_identical(stats::getCall(cox$model)$data, cox$data)
  # Each rotated fit's call carries its own pivot's frame.
  expect_equal(stats::getCall(m$models[[3]])$data$ilr1,
               ilr_transform(m$data[parts], parts, pivot = 3)$ilr1,
               tolerance = 1e-10)
  # The formula, moderator and formula-interface paths store the same shape.
  by_formula <- comp_iv(cd, formula = zbmi ~ ilr() + sex)
  expect_s3_class(stats::getCall(by_formula$model)$formula, "formula")
  expect_true(is.data.frame(stats::getCall(by_formula$model)$data))
  moderated <- comp_iv(cd, "zbmi", moderator = "sex")
  expect_s3_class(stats::getCall(moderated$model)$formula, "formula")
  expect_true(is.data.frame(stats::getCall(moderated$model)$data))
})

# The returned backend fit replays through R's own tools. One engine per
# backend family; coxme only when installed.
replay_fits <- local({
  fits <- NULL
  function() {
    if (is.null(fits)) {
      fits <<- list(
        lm = comp_iv(cd, "zbmi", model = "lm"),
        glm = comp_iv(cd, "obesity", model = "glm", family = "binomial"),
        lmer = comp_iv(cd, "zbmi", model = "lmer", random = "school"),
        cox = comp_iv(cd, time = "time", event = "event", model = "cox")
      )
      if (requireNamespace("coxme", quietly = TRUE)) {
        fits$coxme <<- comp_iv(cd, time = "time", event = "event",
                               model = "coxme", random = "school")
      }
    }
    fits
  }
})

fixed_coefs <- function(model) {
  if (inherits(model, "merMod")) {
    lme4::fixef(model)
  } else if (inherits(model, "coxme")) {
    coxme::fixef(model)
  } else {
    stats::coef(model)
  }
}

test_that("update() on the returned backend fit reproduces it", {
  for (nm in names(replay_fits())) {
    m <- replay_fits()[[nm]]$model
    call <- stats::getCall(m)
    expect_s3_class(call$formula, "formula")
    expect_true(is.data.frame(call$data))
    refit <- stats::update(m)
    expect_identical(class(refit), class(m), label = nm)
    expect_equal(fixed_coefs(refit), fixed_coefs(m), tolerance = 0, label = nm)
  }
})

test_that("a reduced fit from update() compares with anova()", {
  for (nm in names(replay_fits())) {
    m <- replay_fits()[[nm]]$model
    reduced <- suppressMessages(stats::update(m, . ~ . - ilr3))
    expect_false("ilr3" %in% names(fixed_coefs(reduced)), label = nm)
    expect_no_error(suppressMessages(stats::anova(m, reduced)))
  }
})

test_that("cox.zph() runs on the returned Cox fit and matches fit_check()", {
  fit <- replay_fits()$cox
  z <- survival::cox.zph(fit$model)$table
  expect_identical(rownames(z),
                   c(paste0("ilr", seq_len(length(parts) - 1)), "GLOBAL"))
  expect_true(all(z[, "p"] >= 0 & z[, "p"] <= 1))
  tests <- fit_check(fit, plot = FALSE)$tests
  ph <- tests[grepl("^Proportional hazards: ", tests$check), ]
  expect_identical(ph$check, paste0("Proportional hazards: ", rownames(z)))
  expect_equal(unname(ph$statistic), unname(z[, "chisq"]), tolerance = 0)
  expect_equal(unname(ph$df), unname(z[, "df"]), tolerance = 0)
  expect_equal(unname(ph$p.value), unname(z[, "p"]), tolerance = 0)
})

test_that("a clustered mblogit fit writes nothing to the console", {
  skip_if_not_installed("mclogit")
  d <- timeuse
  d$grade_unordered <- factor(as.character(d$grade), levels = levels(d$grade))
  cd2 <- declare_composition(d, parts, closure = 1440)
  # mclogit's inner iterations warn on these data; the warnings are still
  # captured on the fit.
  out <- utils::capture.output(
    fit <- suppressWarnings(
      comp_iv(cd2, "grade_unordered", model = "mblogit", random = "school")
    )
  )
  expect_identical(out, character(0))
  expect_gt(nrow(fit$fit_warnings), 0)
})

test_that("a composition the model cannot estimate leaves the omnibus test empty", {
  # Every day is the same composition, so each coordinate is constant and
  # aliased with the intercept.
  d <- timeuse
  d$MVPA <- 60
  d$LPA <- 300
  d$SB <- 480
  d$SPT <- 600
  expect_warning(
    fit <- comp_iv(declare_composition(d, parts, closure = 1440), "zbmi",
                   model = "lm"),
    "No coefficient of the composition can be estimated"
  )
  expect_true(all(is.na(stats::coef(fit$model)[c("ilr1", "ilr2", "ilr3")])))
  expect_identical(
    fit$omnibus,
    tibble::tibble(term = "composition", statistic = NA_real_, df = 0L,
                   p.value = NA_real_)
  )
})

test_that("a composition estimated in part gives no such warning", {
  # A covariate equal to the first coordinate aliases that coordinate only.
  d <- timeuse
  z <- tc_pivot_coords(as.data.frame(tc_comp_df(cd)[parts]), pivot = 1L)
  d$c1 <- z[, 1]
  expect_no_warning(
    fit <- comp_iv(declare_composition(d, parts, closure = 1440),
                   zbmi ~ c1 + ilr(), model = "lm")
  )
  expect_true(is.na(stats::coef(fit$model)[["ilr1"]]))
  expect_identical(fit$omnibus$df, 2L)
})

test_that("weights local to the caller replay through update()", {
  fit_with_local_weights <- function() {
    w <- rep(2, nrow(timeuse))
    comp_iv(cd, "zbmi", model = "lm", model_args = list(weights = w))
  }
  m <- fit_with_local_weights()
  refit <- stats::update(m$model)
  expect_equal(unname(stats::weights(refit)), rep(2, nrow(timeuse)))
  expect_equal(stats::coef(refit), stats::coef(m$model), tolerance = 0)
})

test_that("reallocations do not depend on the data held in the stored call", {
  for (nm in c("lm", "cox")) {
    fit <- replay_fits()[[nm]]
    bare <- fit
    bare$model$call$data <- quote(unreachable_frame)
    expect_identical(
      as.data.frame(predict_realloc(fit, delta = 10)),
      as.data.frame(predict_realloc(bare, delta = 10)),
      label = nm
    )
  }
})

test_that("print() of the fit names the data in the backend call", {
  out <- utils::capture.output(print(replay_fits()$lm))
  expect_true(any(grepl("data = data)", out, fixed = TRUE)))
  expect_lt(length(out), 20)
  out_lmer <- utils::capture.output(print(replay_fits()$lmer))
  expect_lt(length(out_lmer), 25)
  # The fit itself keeps the frame in its call.
  expect_true(is.data.frame(stats::getCall(replay_fits()$lm$model)$data))
})

test_that("print() counts the captured backend messages and points at summary()", {
  m <- comp_iv(cd, "zbmi", model = "lm")
  expect_no_match(paste(utils::capture.output(print(m)), collapse = "\n"),
                  "captured backend")
  m$fit_warnings <- tibble::tibble(source = "canonical", type = "warning",
                                   message = "a")
  expect_output(print(m), "1 captured backend message (see summary())",
                fixed = TRUE)
  m$fit_warnings <- m$fit_warnings[c(1, 1), ]
  expect_output(print(m), "2 captured backend messages (see summary())",
                fixed = TRUE)
})
