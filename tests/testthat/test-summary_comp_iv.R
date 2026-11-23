parts <- c("MVPA", "LPA", "SB", "SPT")
cd <- declare_composition(timeuse, parts, closure = 1440)

test_that("summary.comp_iv returns tidy global and per-part tibbles", {
  skip_if_not_installed("gt")
  m <- comp_iv(cd, "obesity", model = "glm", family = "binomial")
  s <- summary(m)
  expect_s3_class(s, "summary_comp_iv")
  expect_true(all(c("global", "per_part", "moderation", "convergence",
                    "total_test", "effect_measure", "conf_level") %in% names(s)))
  expect_null(s$moderation)
  expect_null(s$total_test)
  expect_s3_class(s$per_part, "tbl_df")
  expect_s3_class(s$global, "tbl_df")
  expect_s3_class(pub_ready(s)$effects, "gt_tbl")
  expect_setequal(
    names(s$per_part),
    c("parameter", "label", "estimate", "conf.low", "conf.high", "p.value")
  )
  expect_equal(nrow(s$per_part), length(parts))
  expect_true(all(grepl("vs remaining$", s$per_part$parameter)))
})

test_that("labels flow from the compD into the per-part table", {
  cd2 <- declare_composition(timeuse, parts, closure = 1440,
                             labels = c(MVPA = "Exercise"))
  m <- comp_iv(cd2, "obesity", model = "glm", family = "binomial")
  s <- summary(m)
  expect_equal(s$per_part$label[grepl("^MVPA", s$per_part$parameter)],
               "Exercise")
})

test_that("the global row carries the omnibus p-value", {
  m <- comp_iv(cd, "obesity", model = "glm", family = "binomial")
  s <- summary(m)
  expect_equal(nrow(s$global), 1L)
  expect_identical(s$global$parameter, "Global Composition")
  expect_equal(s$global$p.value, m$omnibus$p.value)
  expect_true(is.na(s$global$estimate))
})

test_that("logistic per-part effects are odds ratios (exponentiated, > 0)", {
  m <- comp_iv(cd, "obesity", model = "glm", family = "binomial")
  s <- summary(m)
  expect_true(all(s$per_part$estimate > 0))
  # Matches exp(ilr1) of a hand-rotated refit.
  df <- dplyr::bind_cols(tc_comp_df(cd), timeuse["obesity"])
  d2 <- ilr_transform(df, parts, pivot = 2)
  f2 <- stats::glm(obesity ~ ilr1 + ilr2 + ilr3, data = d2,
                   family = stats::binomial())
  expect_equal(s$per_part$estimate[2],
               unname(exp(stats::coef(f2)["ilr1"])), tolerance = 1e-8)
})

test_that("the effect measure follows the GLM link", {
  pois <- comp_iv(cd, "event", model = "glm", family = "poisson")
  expect_identical(summary(pois)$effect_measure, "Rate ratio")
  binom <- comp_iv(cd, "obesity", model = "glm", family = "binomial")
  expect_identical(summary(binom)$effect_measure, "Odds ratio")
})

test_that("a log-link gaussian reports exponentiated ratios of means", {
  d <- timeuse
  d$y_pos <- exp(d$zbmi / 2)
  cd_pos <- declare_composition(d, parts, closure = 1440)
  m <- comp_iv(cd_pos, "y_pos", model = "glm",
               family = stats::gaussian(link = "log"))
  s <- summary(m)
  expect_identical(s$effect_measure, "Ratio of means")
  expect_equal(
    s$per_part$estimate,
    unname(purrr::map_dbl(parts, function(p) exp(stats::coef(m$models[[p]])[["ilr1"]]))),
    tolerance = 1e-12
  )
  expect_true(all(s$per_part$conf.low > 0))
})

test_that("linear per-part effects are NOT exponentiated", {
  m <- comp_iv(cd, "zbmi", model = "lm")
  s <- summary(m)
  df <- dplyr::bind_cols(tc_comp_df(cd), timeuse["zbmi"])
  d2 <- ilr_transform(df, parts, pivot = 2)
  f2 <- stats::lm(zbmi ~ ilr1 + ilr2 + ilr3, data = d2)
  expect_equal(s$per_part$estimate[2],
               unname(stats::coef(f2)["ilr1"]), tolerance = 1e-8)
})

test_that("Cox per-part effects are hazard ratios", {
  m <- comp_iv(cd, time = "time", event = "event", model = "cox")
  s <- summary(m)
  expect_equal(nrow(s$per_part), length(parts))
  expect_true(all(s$per_part$estimate > 0))
})

test_that("summary.comp_iv returns a visible object and prints via print()", {
  m <- comp_iv(cd, "zbmi", model = "lm")
  s <- summary(m)
  expect_s3_class(s, "summary_comp_iv")
  # summary() itself computes and returns; the text report is print()'s job.
  expect_silent(summary(m))
  expect_output(print(s), "Omnibus composition test")
})

test_that("summary.comp_iv rejects non-comp_iv input", {
  expect_error(summary.comp_iv(list()), "comp_iv")
})

test_that("comp_iv() stores the resolved random spec on the fit", {
  m_fixed <- comp_iv(cd, "zbmi", model = "lm")
  expect_null(m_fixed$random)

  m_clustered <- comp_iv(cd, "zbmi", model = "lmer", random = "school")
  expect_identical(m_clustered$random, "school")
})

test_that("summary.comp_iv reports per-part moderation tests for a 2-level moderator", {
  skip_if_not_installed("gt")
  m <- comp_iv(cd, "obesity", model = "glm", family = "binomial",
              moderator = "sex")
  s <- summary(m)
  expect_true(all(c("mod_statistic", "mod_df", "mod_p.value") %in%
                 names(s$per_part)))
  expect_equal(nrow(s$per_part), length(parts))
  expect_true(all(is.finite(s$per_part$mod_p.value)))
  expect_true(all(s$per_part$mod_df == 1))
  expect_s3_class(pub_ready(s)$effects, "gt_tbl")
})

test_that("summary.comp_iv reports per-part moderation tests for a 3-level factor moderator", {
  m <- comp_iv(cd, "obesity", model = "glm", family = "binomial",
              moderator = "region")
  s <- summary(m)
  expect_true(all(s$per_part$mod_df == nlevels(cd$data_raw$region) - 1))
  expect_true(all(is.finite(s$per_part$mod_p.value)))
})

test_that("per-part moderation columns are absent for an unmoderated fit", {
  m <- comp_iv(cd, "obesity", model = "glm", family = "binomial")
  s <- summary(m)
  expect_false(any(c("mod_statistic", "mod_df", "mod_p.value") %in%
                   names(s$per_part)))
})

test_that("the canonical-pivot per-part moderation p matches a hand-computed Wald test", {
  m <- comp_iv(cd, "obesity", model = "glm", family = "binomial",
              moderator = "sex")
  s <- summary(m)

  # The canonical fit (m$model) is the pivot-1 model, i.e. exactly
  # object$models[[parts[1]]] -- so its own ilr1:sex coefficient should match
  # the first per-part row's moderation test computed directly from coef/vcov.
  b <- stats::coef(m$model)
  V <- stats::vcov(m$model)
  int_name <- grep("^ilr1:sex", names(b), value = TRUE)
  expect_length(int_name, 1)
  stat <- as.numeric(b[int_name]^2 / V[int_name, int_name])
  p <- stats::pchisq(stat, df = 1, lower.tail = FALSE)

  expect_equal(s$per_part$mod_p.value[1], p, tolerance = 1e-8)
  expect_equal(s$per_part$mod_statistic[1], stat, tolerance = 1e-8)
})

# --- the per-part interval and its p-value must share a distribution ---------
# Both directions of this were regressions caught by /code-review medium: keying
# the df off the statistic column ("t value") rather than the p-value column
# gave lmer a t interval beside a normal p, and called df.residual() on an
# rlmerMod, which errors rather than warns.

test_that("the per-part df comes from whatever produced the p-value", {
  expect_consistent <- function(fit, model_type) {
    st <- tc_iv_ilr1_stats(fit$models[[1]], model_type)
    mat <- switch(model_type, cox = summary(fit$models[[1]])$coefficients,
                  stats::coef(summary(fit$models[[1]])))
    pcol <- grep("^Pr", colnames(mat), value = TRUE)
    # a t-based p-value means a finite df; anything else must be asymptotic
    if (length(pcol) == 1 && grepl("t", pcol, fixed = TRUE)) {
      expect_true(is.finite(st$df[1]))
    } else {
      expect_identical(st$df[1], Inf)
    }
  }
  expect_consistent(comp_iv(cd, "zbmi", model = "lm"), "lm")
  expect_consistent(
    comp_iv(cd, "zbmi", model = "glm", family = "gaussian"), "glm")
  expect_consistent(
    comp_iv(cd, "obesity", model = "glm", family = "binomial"), "glm")
  # glmer reports Pr(>|z|), so its row is asymptotic, as its emmeans
  # predictions are; the lmer row is the Kenward-Roger case tested next.
  expect_consistent(
    comp_iv(cd, "obesity", model = "glmer", family = "binomial",
            random = "school"), "glmer")
})

test_that("lmer per-part rows carry the Kenward-Roger SE and df that emmeans uses", {
  skip_if_not_installed("pbkrtest")
  fit <- comp_iv(cd, "zbmi", model = "lmer", random = "school")
  s <- summary(fit)
  for (k in seq_along(parts)) {
    m <- fit$models[[k]]
    st <- tc_iv_ilr1_stats(m, "lmer")
    et <- summary(emmeans::emtrends(m, ~ 1, var = "ilr1"))
    expect_equal(st$est, et$ilr1.trend, tolerance = 1e-10)
    expect_equal(st$se, et$SE, tolerance = 1e-10)
    expect_equal(st$df, et$df, tolerance = 1e-10)
    expect_true(is.finite(st$df))
    expect_equal(c(s$per_part$conf.low[k], s$per_part$conf.high[k]),
                 c(et$lower.CL, et$upper.CL), tolerance = 1e-10)
    # The p-value is referred to the same t distribution as the interval.
    expect_equal(s$per_part$p.value[k],
                 2 * stats::pt(-abs(st$est / st$se), df = st$df))
  }
  # The same conditions under which emmeans falls back to asymptotic make the
  # row asymptotic: the plain covariance and df = Inf.
  m1 <- fit$models[[1]]
  se_plain <- sqrt(as.matrix(stats::vcov(m1))["ilr1", "ilr1"])
  opt_names <- c("lmer.df", "disable.pbkrtest", "pbkrtest.limit",
                 "disable.lmerTest", "lmerTest.limit")
  old <- stats::setNames(purrr::map(opt_names, emmeans::get_emm_option),
                         opt_names)
  withr::defer(do.call(emmeans::emm_options, old))
  emmeans::emm_options(lmer.df = "asymptotic")
  st <- tc_iv_ilr1_stats(m1, "lmer")
  expect_identical(st$df, Inf)
  expect_equal(st$se, se_plain)
  emmeans::emm_options(lmer.df = old$lmer.df, disable.pbkrtest = TRUE,
                       disable.lmerTest = TRUE)
  expect_identical(tc_iv_ilr1_stats(m1, "lmer")$df, Inf)
  emmeans::emm_options(disable.pbkrtest = FALSE, disable.lmerTest = FALSE,
                       pbkrtest.limit = 10, lmerTest.limit = 10)
  expect_identical(tc_iv_ilr1_stats(m1, "lmer")$df, Inf)
})

test_that("lmer per-part rows follow emmeans onto Satterthwaite df", {
  skip_if_not_installed("lmerTest")
  fit <- comp_iv(cd, "zbmi", model = "lmer", random = "school")
  opt_names <- c("lmer.df", "disable.pbkrtest")
  old <- stats::setNames(purrr::map(opt_names, emmeans::get_emm_option),
                         opt_names)
  withr::defer(do.call(emmeans::emm_options, old))
  settings <- list(list(disable.pbkrtest = TRUE),
                   list(disable.pbkrtest = FALSE, lmer.df = "satterthwaite"))
  for (setting in settings) {
    do.call(emmeans::emm_options, setting)
    s <- summary(fit)
    for (k in seq_along(parts)) {
      m <- fit$models[[k]]
      st <- tc_iv_ilr1_stats(m, "lmer")
      et <- summary(emmeans::emtrends(m, ~ 1, var = "ilr1"))
      expect_equal(st$se, et$SE, tolerance = 1e-10)
      expect_equal(st$df, et$df, tolerance = 1e-10)
      expect_true(is.finite(st$df))
      expect_equal(c(s$per_part$conf.low[k], s$per_part$conf.high[k]),
                   c(et$lower.CL, et$upper.CL), tolerance = 1e-10)
    }
    do.call(emmeans::emm_options, old)
  }
})

test_that("lmer per-part rows follow emmeans when only one df package is usable", {
  skip_if_not_installed("pbkrtest")
  fit <- comp_iv(cd, "zbmi", model = "lmer", random = "school")
  opt_names <- c("lmer.df", "disable.pbkrtest", "pbkrtest.limit",
                 "disable.lmerTest", "lmerTest.limit")
  old <- stats::setNames(purrr::map(opt_names, emmeans::get_emm_option),
                         opt_names)
  withr::defer(do.call(emmeans::emm_options, old))
  # Satterthwaite asked for without lmerTest goes to Kenward-Roger when
  # pbkrtest is usable; Kenward-Roger without either package is asymptotic.
  settings <- list(
    "kenward-roger" = list(lmer.df = "satterthwaite", disable.lmerTest = TRUE),
    asymptotic = list(lmer.df = "kenward-roger", disable.pbkrtest = TRUE,
                      disable.lmerTest = TRUE)
  )
  for (nm in names(settings)) {
    do.call(emmeans::emm_options, settings[[nm]])
    s <- summary(fit)
    for (k in seq_along(parts)) {
      m <- fit$models[[k]]
      st <- tc_iv_ilr1_stats(m, "lmer")
      et <- summary(emmeans::emtrends(m, ~ 1, var = "ilr1"))
      expect_match(attr(et, "mesg"), paste("method:", nm), all = FALSE)
      expect_equal(st$se, et$SE, tolerance = 1e-10)
      expect_equal(st$df, et$df, tolerance = 1e-10)
      expect_identical(is.finite(st$df), nm == "kenward-roger")
      ci <- if (nm == "asymptotic") {
        c(et$asymp.LCL, et$asymp.UCL)
      } else {
        c(et$lower.CL, et$upper.CL)
      }
      expect_equal(c(s$per_part$conf.low[k], s$per_part$conf.high[k]), ci,
                   tolerance = 1e-10)
    }
    do.call(emmeans::emm_options, old)
  }
})

test_that("summary() works for rlmer, whose df.residual() errors", {
  skip_if_not_installed("robustlmm")
  fit <- suppressWarnings(
    comp_iv(cd, "zbmi", model = "rlmer", random = "school")
  )
  s <- summary(fit)
  expect_s3_class(s, "summary_comp_iv")
  expect_equal(nrow(s$per_part), length(parts))
  expect_true(all(is.finite(s$per_part$conf.low)))
})

test_that("rlmer reallocations and surfaces are asymptotic, like its per-part rows", {
  skip_if_not_installed("robustlmm")
  rl <- comp_iv(cd, "zbmi", model = "rlmer", random = "school")
  r <- predict_realloc(rl, delta = 15, from = "SB", to = "MVPA")
  expect_identical(r$df, Inf)
  z <- stats::qnorm(0.975)
  expect_equal(c(r$conf.low, r$conf.high),
               r$estimate + c(-1, 1) * z * r$std.error, tolerance = 1e-10)
  expect_no_warning(p <- plot(rl, resolution = 6, quantiles = NULL, smooth = FALSE))
  expect_length(p, choose(length(parts), 3))
})

test_that("Convergence: OK is only claimed when the engine reported success", {
  g <- comp_iv(cd, "obesity", model = "glm", family = "binomial")
  expect_output(print(summary(g)), "Convergence: OK")
  # lm has no convergence concept (NA), so it says nothing either way
  l <- comp_iv(cd, "zbmi", model = "lm")
  expect_false(any(grepl("Convergence",
                         utils::capture.output(print(summary(l))))))
})

test_that("a captured backend message containing braces prints without error", {
  cd <- declare_composition(timeuse, parts, closure = 1440)
  m <- comp_iv(cd, "obesity", model = "glm", family = "binomial")
  m$fit_warnings <- tibble::tibble(
    source = "canonical", type = "warning",
    message = "optimizer said {something} in braces"
  )
  m$convergence$n_warnings <- 1L
  expect_message(print(summary(m)), "\\{something\\}")
})

test_that("with non-linear composition terms the per-part rows are slopes at the reference", {
  q <- comp_iv(cd, zbmi ~ ilr()^2 + I(ilr()^2), model = "lm")
  s <- summary(q)
  # The analytic derivative of the MVPA-first fit in ilr1 at the reference.
  z <- as.data.frame(tc_pivot_coords(as.data.frame(q$ref_comp), pivot = 1L))
  b <- stats::coef(q$models[["MVPA"]])
  exact <- b[["ilr1"]] + 2 * b[["I(ilr1^2)"]] * z$ilr1 +
    b[["ilr1:ilr2"]] * z$ilr2 + b[["ilr1:ilr3"]] * z$ilr3
  expect_equal(s$per_part$estimate[1], exact, tolerance = 1e-6)
  expect_lt(s$per_part$estimate[parts == "LPA"], 0)
  p <- summary(comp_iv(cd, obesity ~ poly(ilr(), 2), model = "glm",
                       family = "binomial"))
  expect_true(all(is.finite(unlist(
    p$per_part[c("estimate", "conf.low", "conf.high", "p.value")]
  ))))
  expect_output(print(p), "slopes at the reference composition")
  # A composition crossed with a moderator keeps the plain ilr1 coefficient.
  mod <- comp_iv(cd, "zbmi", moderator = "sex", model = "lm")
  expect_equal(summary(mod)$per_part$estimate[2],
               unname(stats::coef(mod$models[["LPA"]])["ilr1"]),
               tolerance = 1e-12)
})

test_that("the gt table notes per-part slopes only when the rows are slopes", {
  skip_if_not_installed("gt")
  notes <- function(f) {
    all <- unlist(pub_ready(summary(comp_iv(cd, f, model = "lm")))$effects[["_source_notes"]])
    sum(grepl("slopes at the reference composition", all, fixed = TRUE))
  }
  expect_identical(c(notes(zbmi ~ ilr()^2 + I(ilr()^2)), notes(zbmi ~ ilr())),
                   c(1L, 0L))
})

test_that("summary() reports the curvature test of a full quadratic only", {
  curved <- comp_iv(cd, "zbmi", model = "lm", curvature = TRUE)
  s <- summary(curved)
  expect_identical(s$curvature_test, curved$curvature_test)
  ct <- curved$curvature_test
  line <- glue::glue(
    "  Curvature (Wald): chi-squared(6) = {tc_fmt(ct$statistic, 3)}, p = {format.pval(ct$p.value, digits = 3)}"
  )
  out <- utils::capture.output(print(s))
  expect_true(line %in% out)
  # The line follows the moderation line on a moderated fit.
  mod <- utils::capture.output(print(summary(
    comp_iv(cd, "zbmi", model = "lm", moderator = "sex", curvature = TRUE)
  )))
  expect_match(mod[grep("^\\s*Moderation by sex", mod) + 1],
               "^\\s*Curvature \\(Wald\\): chi-squared\\(6\\)")
  s0 <- summary(comp_iv(cd, "zbmi", model = "lm"))
  expect_null(s0$curvature_test)
  expect_false(any(grepl("Curvature", utils::capture.output(print(s0)))))
  skip_if_not_installed("gt")
  notes <- unlist(pub_ready(s)$effects[["_source_notes"]])
  expect_true(as.character(trimws(line)) %in% notes)
  notes0 <- unlist(pub_ready(s0)$effects[["_source_notes"]])
  expect_false(any(grepl("Curvature", notes0)))
})

# The derivative in ilr1, at the reference composition, of the linear predictor
# of the full quadratic fitted with part k as the pivot.
quadratic_slope <- function(fit, k, b) {
  z <- as.data.frame(tc_pivot_coords(as.data.frame(fit$ref_comp), pivot = k))
  b[["ilr1"]] + 2 * b[["I(ilr1^2)"]] * z$ilr1 + b[["ilr1:ilr2"]] * z$ilr2 +
    b[["ilr1:ilr3"]] * z$ilr3
}

test_that("a curved ordinal fit reports each part's slope at the reference", {
  skip_if_not_installed("ordinal")
  fit <- comp_iv(cd, "grade", model = "clm", curvature = TRUE)
  s <- summary(fit)
  hand <- purrr::map_dbl(seq_along(parts), function(k) {
    quadratic_slope(fit, k, stats::coef(fit$models[[k]]))
  })
  # The summary's cumulative odds ratio is the exponentiated slope of the
  # latent linear predictor.
  expect_equal(log(s$per_part$estimate[1]), hand[1], tolerance = 1e-6)
  expect_equal(log(s$per_part$estimate), hand, tolerance = 1e-5)
  expect_true(all(s$per_part$conf.low < s$per_part$estimate &
                    s$per_part$estimate < s$per_part$conf.high))
  expect_true(s$per_part_slopes)
  expect_null(s$per_part_no_slope)
  expect_output(print(s), "slopes at the reference composition")
})

test_that("a curved robust mixed fit reports asymptotic slopes at the reference", {
  skip_if_not_installed("robustlmm")
  fit <- comp_iv(cd, zbmi ~ ilr()^2 + I(ilr()^2) + (1 | school),
                 model = "rlmer")
  s <- summary(fit)
  hand <- purrr::map_dbl(seq_along(parts), function(k) {
    quadratic_slope(fit, k, lme4::fixef(fit$models[[k]]))
  })
  expect_equal(s$per_part$estimate[1], hand[1], tolerance = 1e-6)
  expect_equal(s$per_part$estimate, hand, tolerance = 1e-5)
  # The interval and the p-value both use the normal distribution.
  z <- stats::qnorm(0.975)
  se <- (s$per_part$conf.high - s$per_part$conf.low) / (2 * z)
  expect_true(all(is.finite(se)))
  expect_equal(s$per_part$p.value,
               2 * stats::pnorm(-abs(s$per_part$estimate / se)),
               tolerance = 1e-8)
})

test_that("an engine with no slope at the reference gives NA rows and the reason", {
  na_rows <- function(s) {
    all(is.na(unlist(s$per_part[c("estimate", "conf.low", "conf.high",
                                  "p.value")])))
  }
  skip_if_not_installed("robustbase")
  rob <- summary(comp_iv(cd, zbmi ~ ilr()^2 + I(ilr()^2), model = "lmrob"))
  expect_true(na_rows(rob))
  expect_match(rob$per_part_no_slope, "lmrob", fixed = TRUE)
  expect_output(print(rob), paste("NA rows:", rob$per_part_no_slope),
                fixed = TRUE)
  skip_if_not_installed("mclogit")
  d <- timeuse
  d$grade_unordered <- factor(as.character(d$grade), levels = levels(d$grade))
  cd2 <- declare_composition(d, parts, closure = 1440)
  mb <- summary(suppressWarnings(
    comp_iv(cd2, "grade_unordered", model = "mblogit", curvature = TRUE)
  ))
  expect_true(na_rows(mb))
  expect_identical(nrow(mb$per_part), 2L * length(parts))
  expect_setequal(mb$per_part$response_level, c("expected", "above"))
  expect_identical(
    mb$per_part_no_slope,
    "a per-category model has no single linear predictor to differentiate."
  )
  expect_output(print(mb), "NA rows: a per-category model", fixed = TRUE)
})

test_that("a clustered engine with no intraclass correlation reports none", {
  skip_if_not_installed("mclogit")
  d <- timeuse
  d$grade_unordered <- factor(as.character(d$grade), levels = levels(d$grade))
  cd2 <- declare_composition(d, parts, closure = 1440)
  # mclogit's inner iterations warn on these data; the warnings stay on the fit.
  fit <- suppressWarnings(
    comp_iv(cd2, grade_unordered ~ ilr() + (1 | school), model = "mblogit")
  )
  s <- summary(fit)
  expect_identical(
    s$icc,
    tibble::tibble(group = "school", variance = NA_real_,
                   measure = "not available", value = NA_real_)
  )
  out <- utils::capture.output(suppressMessages(print(s)))
  expect_true("  Intraclass correlation: not available for this model" %in% out)
  expect_true(all(is.finite(s$per_part$estimate)))
  expect_null(s$per_part_no_slope)
})

test_that("pub_ready() carries the total test of an include_total fit as a source note", {
  skip_if_not_installed("gt")
  set.seed(7)
  vt <- timeuse
  vt[parts] <- vt[parts] * stats::rlnorm(nrow(vt), 0, 0.3)
  cdt <- declare_composition(vt, parts, closure = "mean")
  s <- summary(comp_iv(cdt, "zbmi", model = "lm", include_total = TRUE))
  notes <- unlist(pub_ready(s)$effects[["_source_notes"]])
  expect_length(notes, 1)
  expect_match(notes, "Total (log geometric mean, Wald): beta = ", fixed = TRUE)
  expect_true(any(grepl(notes, capture.output(print(s)), fixed = TRUE)))
})

test_that("summary() refuses a conf_level outside (0, 1)", {
  m <- comp_iv(cd, "zbmi", model = "lm")
  for (bad in list(95, 0, 1, NA_real_, c(0.9, 0.95))) {
    expect_error(summary(m, conf_level = bad),
                 "conf_level.*single number in \\(0, 1\\)")
  }
})

test_that("pub_ready() labels the interval with the confidence level as given", {
  skip_if_not_installed("gt")
  m <- comp_iv(cd, "zbmi", model = "lm")
  ci_label <- function(level) {
    tab <- pub_ready(summary(m, conf_level = level))$effects
    labs <- as.character(unlist(tab[["_boxhead"]]$column_label))
    labs[grepl("CI", labs, fixed = TRUE)]
  }
  expect_identical(ci_label(0.975), "97.5% CI")
  expect_identical(ci_label(0.95), "95% CI")
})
