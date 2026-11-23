parts <- c("MVPA", "LPA", "SB", "SPT")

cd <- declare_composition(timeuse, parts, closure = 1440)
glm_fit <- comp_iv(cd, "obesity", model = "glm", family = "binomial")
lm_fit <- comp_iv(cd, "zbmi", model = "lm")
cox_fit <- comp_iv(cd, time = "time", event = "event", model = "cox")

test_that("relative mode returns the documented columns", {
  out <- predict_realloc(glm_fit, delta = 30, from = "SB", to = "MVPA")
  expect_s3_class(out, "realloc")
  expect_s3_class(out, "tbl_df")
  expect_named(out, c("from", "to", "delta", "estimate", "std.error", "df",
                      "conf.low", "conf.high", "statistic", "p.value"))
  expect_true(all(out$conf.low <= out$estimate & out$estimate <= out$conf.high))
})

test_that("absolute mode returns a reference row and predictions", {
  out <- predict_realloc(lm_fit, delta = 30, from = "SB", to = "MVPA",
                         type = "absolute")
  expect_true(all(c("from", "to", "delta", "estimate") %in% names(out)))
  expect_true(is.na(out$from[1]) && out$delta[1] == 0)
})

test_that("reference-vs-reference contrast is ~1 on the response scale", {
  # A near-zero reallocation must give an OR indistinguishable from 1.
  out <- predict_realloc(glm_fit, delta = 1e-4, from = "SB", to = "MVPA")
  expect_equal(out$estimate, 1, tolerance = 1e-3)
})

test_that("reallocation results are invariant to the ILR basis", {
  refit <- function(piv) {
    d <- ilr_transform(timeuse, parts, total = 1440, pivot = piv)
    m <- glm(obesity ~ ilr1 + ilr2 + ilr3, data = d, family = binomial)
    structure(
      list(model = m, model_type = "glm", parts = parts, total = 1440,
           pivot = as.integer(piv),
           ref_comp = geom_mean_comp(timeuse, parts, total = 1440),
           covariates = NULL, has_ilr_interaction = FALSE,
           family = "binomial", data = d),
      class = "timecoda_fit"
    )
  }
  a <- predict_realloc(refit(1), delta = 30, from = "SB", to = "MVPA")
  b <- predict_realloc(refit(3), delta = 30, from = "SB", to = "MVPA")
  expect_equal(a$estimate, b$estimate, tolerance = 1e-8)
})

test_that("Cox supports relative output but rejects absolute mode", {
  hr <- predict_realloc(cox_fit, delta = 30, from = "SB", to = "MVPA")
  expect_true(is.finite(hr$estimate))
  expect_error(
    predict_realloc(cox_fit, delta = 30, from = "SB", to = "MVPA",
                    type = "absolute"),
    "relative"
  )
})

test_that("non-positive reallocations get NA estimates by default", {
  out <- predict_realloc(glm_fit, delta = c(30, 1e6), from = "SB", to = "MVPA")
  expect_equal(nrow(out), 2L)
  expect_true(is.na(out$estimate[out$delta == 1e6]))
  expect_false(is.na(out$estimate[out$delta == 30]))
  expect_equal(out, predict_realloc(glm_fit, delta = c(30, 1e6), from = "SB",
                                    to = "MVPA", on_invalid = "na"))
})

test_that("on_invalid = 'error' rejects a non-positive reallocation", {
  expect_error(
    predict_realloc(glm_fit, delta = 1e6, from = "MVPA", to = "SB",
                    on_invalid = "error"),
    "non-positive"
  )
})

test_that("on_invalid = 'omit' drops only the impossible cells", {
  out <- predict_realloc(glm_fit, delta = c(30, 1e6), from = "SB", to = "MVPA",
                         on_invalid = "omit")
  expect_equal(nrow(out), 1L)
  expect_equal(out$delta, 30)
  expect_true(is.finite(out$estimate))
})

test_that("on_invalid = 'na' keeps the row with NA results", {
  out <- predict_realloc(glm_fit, delta = c(30, 1e6), from = "SB", to = "MVPA",
                         on_invalid = "na")
  expect_equal(nrow(out), 2L)
  expect_true(is.na(out$estimate[out$delta == 1e6]))
  expect_true(is.na(out$p.value[out$delta == 1e6]))
  expect_false(is.na(out$estimate[out$delta == 30]))
})

test_that("one-to-all reallocation keeps the composition closed and positive", {
  ref <- unlist(geom_mean_comp(timeuse, parts, total = 1440))
  moves <- tc_realloc_focals(parts, from = NULL, to = "MVPA")
  g <- tc_realloc_grid(ref, moves, delta = c(10, 30), scheme = "one-to-all")
  expect_true(all(abs(rowSums(g$comps) - 1440) < 1e-8))
  expect_true(all(g$comps > 0))
  expect_true(all(g$to == "MVPA") && all(g$from == tc_rest_key()))
  # The focal part gains exactly delta; the rest give it up in proportion.
  expect_equal(unname(g$comps[, "MVPA"] - ref[["MVPA"]]), c(10, 30))
  rest <- setdiff(parts, "MVPA")
  expect_equal(unname(g$comps[1, rest] / ref[rest]),
               rep(1 - 10 / sum(ref[rest]), length(rest)))
})

test_that("a giving focal part at a user reference keeps the rest's ratios", {
  ref <- c(MVPA = 30, LPA = 240, SB = 660, SPT = 510)
  g <- tc_realloc_grid(ref, tc_realloc_focals(parts, from = "SB", to = NULL),
                       delta = 25, scheme = "one-to-all")
  rest <- setdiff(parts, "SB")
  expect_equal(sum(g$comps), 1440)
  expect_equal(unname(g$comps[1, "SB"]), 635)
  expect_equal(unname(g$comps[1, rest] / ref[rest]),
               rep(1 + 25 / sum(ref[rest]), length(rest)))
})

test_that("each scheme's grid reproduces its closed-form reallocation", {
  ref <- unlist(geom_mean_comp(timeuse, parts, total = 1440))
  pairs <- tc_realloc_pairs(parts, from = "SB", to = "MVPA")
  g <- tc_realloc_grid(ref, pairs, delta = c(-15, 15), scheme = "one-to-one")
  expect_equal(unname(g$comps[, "SB"]), ref[["SB"]] - c(-15, 15))
  expect_equal(unname(g$comps[, "MVPA"]), ref[["MVPA"]] + c(-15, 15))
  expect_equal(unname(g$comps[, c("LPA", "SPT")]),
               unname(rbind(ref[c("LPA", "SPT")], ref[c("LPA", "SPT")])))
  # Feasibility is flagged, not aborted, when asked: MVPA cannot give 1e3 min.
  bad <- tc_realloc_grid(ref, tc_realloc_pairs(parts, "MVPA", "SB"),
                         delta = c(10, 1e3), scheme = "one-to-one",
                         on_invalid = "na")
  expect_equal(bad$valid, c(TRUE, FALSE))
})

test_that("every directed scheme has a definition and a figure layout", {
  schemes <- eval(formals(predict_realloc)$scheme)
  expect_setequal(names(tc_realloc_schemes), schemes)
  expect_setequal(names(tc_realloc_layouts), schemes)
  expect_setequal(eval(formals(plot_realloc)$scheme), schemes)
})

test_that("one-to-all output labels the implicit side and is basis-invariant", {
  refit <- function(piv) {
    d <- ilr_transform(timeuse, parts, total = 1440, pivot = piv)
    m <- glm(obesity ~ ilr1 + ilr2 + ilr3, data = d, family = binomial)
    structure(
      list(model = m, model_type = "glm", parts = parts, total = 1440,
           pivot = as.integer(piv),
           ref_comp = geom_mean_comp(timeuse, parts, total = 1440),
           covariates = NULL, has_ilr_interaction = FALSE,
           family = "binomial", data = d),
      class = "timecoda_fit"
    )
  }
  a <- predict_realloc(refit(1), delta = 30, scheme = "one-to-all", to = "MVPA")
  b <- predict_realloc(refit(3), delta = 30, scheme = "one-to-all", to = "MVPA")
  expect_true(all(a$to == "MVPA") && all(a$from == "(rest)"))
  expect_equal(a$estimate, b$estimate, tolerance = 1e-8)
})

test_that("one-to-all defaults to every part as the focal receiver", {
  r <- predict_realloc(glm_fit, delta = 30, scheme = "one-to-all")
  expect_equal(nrow(r), length(parts))
  expect_setequal(r$to, parts)
  expect_true(all(r$from == "(rest)"))
})

test_that("one-to-all rejects supplying both from and to", {
  expect_error(
    predict_realloc(glm_fit, delta = 30, scheme = "one-to-all",
                    from = "SB", to = "MVPA"),
    "only one"
  )
})

test_that("predict_realloc(reference=) validates the reference composition", {
  ref <- c(MVPA = 30, LPA = 240, SB = 660, SPT = 510)
  out <- predict_realloc(lm_fit, delta = 30, reference = ref,
                         from = "SB", to = "MVPA")
  expect_s3_class(out, "tbl_df")
  expect_error(
    predict_realloc(lm_fit, delta = 30,
                    reference = c(MVPA = 30, LPA = 240, SB = 660, SPT = 1),
                    from = "SB", to = "MVPA"),
    "total"
  )
  expect_error(
    predict_realloc(lm_fit, delta = 30, reference = c(MVPA = 30, LPA = 240),
                    from = "SB", to = "MVPA"),
    "one value per part"
  )
  # A missing value, or a one-row data frame lacking a part, gets the same
  # message rather than a base R error.
  expect_error(
    predict_realloc(lm_fit, delta = 30, reference = replace(ref, "SPT", NA),
                    from = "SB", to = "MVPA"),
    "one value per part"
  )
  expect_error(
    predict_realloc(lm_fit, delta = 30, reference = data.frame(MVPA = 30),
                    from = "SB", to = "MVPA"),
    "one value per part"
  )
})

test_that("a supplied reference moves the anchor; NULL is the fit's own", {
  ref <- c(MVPA = 30, LPA = 240, SB = 660, SPT = 510)
  via_arg <- predict_realloc(lm_fit, delta = 30, reference = ref,
                             from = "SB", to = "MVPA")
  default_ref <- predict_realloc(lm_fit, delta = 30, from = "SB", to = "MVPA")
  expect_false(isTRUE(all.equal(default_ref$estimate, via_arg$estimate)))
  own <- predict_realloc(lm_fit, delta = 30, from = "SB", to = "MVPA",
                         reference = unlist(lm_fit$ref_comp))
  expect_equal(as.data.frame(own), as.data.frame(default_ref))
})

test_that("by = computes group-wise reallocations matching a manual at-pin", {
  mod <- comp_iv(cd, outcome = "obesity", moderator = "sex", model = "glm",
                 family = "binomial")
  # Auto-by defaults to the moderator.
  auto <- predict_realloc(mod, from = "SB", to = "MVPA", delta = 15)
  expect_true("sex" %in% names(auto))
  expect_equal(nrow(auto), nlevels(factor(timeuse$sex)))

  # Explicit by matches the auto default.
  expl <- predict_realloc(mod, from = "SB", to = "MVPA", delta = 15, by = "sex")
  expect_equal(expl$estimate, auto$estimate, tolerance = 1e-10)

  # And each group equals the manual single-level at-pin.
  lev <- levels(factor(timeuse$sex))[1]
  manual <- predict_realloc(mod, from = "SB", to = "MVPA", delta = 15,
                            at = list(sex = lev))
  grp <- auto[auto$sex == lev, ]
  expect_equal(grp$estimate, manual$estimate, tolerance = 1e-8)
})

test_that("summary.realloc carries the group column", {
  mod <- comp_iv(cd, outcome = "obesity", moderator = "sex", model = "glm",
                 family = "binomial")
  r <- predict_realloc(mod, from = "SB", to = "MVPA", delta = c(10, 20),
                       on_invalid = "na")
  s <- summary(r)
  expect_true("sex" %in% names(s$table))
})

test_that("by rejects a numeric covariate", {
  fit <- comp_iv(cd, obesity ~ ilr() * zbmi, model = "glm", family = "binomial")
  expect_error(
    predict_realloc(fit, from = "SB", to = "MVPA", delta = 15, by = "zbmi"),
    "factor or character"
  )
})

# --- Scenario (b): Cox rejects type = "absolute" -----------------------------

test_that("predict_realloc rejects type = 'absolute' for a Cox fit with an informative error", {
  expect_error(
    predict_realloc(cox_fit, delta = 30, from = "SB", to = "MVPA",
                    type = "absolute"),
    "absolute|relative"
  )
})

test_that("the relative-only message names the engine's own reason", {
  expect_error(
    predict_realloc(cox_fit, delta = 30, from = "SB", to = "MVPA",
                    type = "absolute"),
    "Cox\\s+model\\s+leaves\\s+the\\s+baseline\\s+hazard\\s+unestimated"
  )
  expect_error(plot(cox_fit, scale = "absolute"),
               "Cox\\s+model\\s+leaves\\s+the\\s+baseline\\s+hazard\\s+unestimated")
  skip_if_not_installed("ordinal")
  clm_fit <- comp_iv(cd, "grade", model = "clm")
  msg <- tryCatch(
    predict_realloc(clm_fit, delta = 30, from = "SB", to = "MVPA",
                    type = "absolute"),
    error = conditionMessage
  )
  expect_match(msg, "ordinal\\s+predictions\\s+are\\s+per\\s+outcome\\s+category")
  expect_no_match(msg, "hazard")
})

# --- Scenario (c): ILR basis invariance via two declare_composition() calls --

test_that("predict_realloc agrees across two compD declarations with parts in different orders", {
  small_parts <- c("MVPA", "SB", "SPT")
  reordered_parts <- c("SPT", "MVPA", "SB")

  cd_a <- declare_composition(timeuse, small_parts, closure = 1440)
  cd_b <- declare_composition(timeuse, reordered_parts, closure = 1440)

  fit_a <- comp_iv(cd_a, "obesity", model = "glm", family = "binomial")
  fit_b <- comp_iv(cd_b, "obesity", model = "glm", family = "binomial")

  out_a <- predict_realloc(fit_a, delta = c(10, 20), from = "SB", to = "MVPA",
                           on_invalid = "na")
  out_b <- predict_realloc(fit_b, delta = c(10, 20), from = "SB", to = "MVPA",
                           on_invalid = "na")

  # Match rows on from/to/delta rather than assuming identical row order.
  key_a <- paste(out_a$from, out_a$to, out_a$delta, sep = "|")
  key_b <- paste(out_b$from, out_b$to, out_b$delta, sep = "|")
  out_b <- out_b[match(key_a, key_b), ]

  expect_equal(out_a$estimate, out_b$estimate, tolerance = 1e-8)
  expect_equal(out_a$conf.low, out_b$conf.low, tolerance = 1e-8)
  expect_equal(out_a$conf.high, out_b$conf.high, tolerance = 1e-8)
})

# --- Scenario (a): multinomial end-to-end (predict_realloc leg) --------------

test_that("predict_realloc carries response_level through a multinom fit on shipped grade", {
  skip_if_not_installed("mclogit")
  small_parts <- c("MVPA", "SB", "SPT")
  d <- timeuse
  d$grade_cat <- factor(as.character(d$grade), levels = levels(d$grade))
  cd_multi <- declare_composition(d, small_parts, closure = 1440)
  fit <- suppressWarnings(
    comp_iv(cd_multi, "grade_cat", model = "mblogit")
  )
  expect_identical(fit$model_type, "mblogit")

  k <- nlevels(d$grade_cat)
  rr <- suppressWarnings(
    predict_realloc(fit, delta = 20, from = "SB", to = "MVPA")
  )
  expect_true("response_level" %in% names(rr))
  expect_equal(nrow(rr), k - 1)
  expect_true(all(rr$estimate > 0))

  # Reference-vs-reference contrast: RRR ~ 1 for every non-reference category.
  r0 <- suppressWarnings(
    predict_realloc(fit, delta = 1e-6, from = "SB", to = "MVPA")
  )
  expect_equal(r0$estimate, rep(1, k - 1), tolerance = 1e-3)
})

# --- Error messages ------------------------------------------------------------

test_that("an unknown from/to part is named, singly or in the plural", {
  expect_error(
    predict_realloc(lm_fit, delta = 15, from = "SB", to = "FOO"),
    '`from`/`to` names "FOO", which is not in the fitted composition'
  )
  expect_error(
    predict_realloc(lm_fit, delta = 15, from = c("SB", "BAR"), to = "FOO"),
    '"BAR" and "FOO", which are not in the fitted'
  )
  expect_error(
    predict_realloc(lm_fit, delta = 15, scheme = "one-to-all", to = "FOO"),
    '"FOO", which is not in the fitted'
  )
  expect_error(
    predict_realloc(lm_fit, delta = 15, scheme = "one-to-all",
                    to = c("FOO", "BAR")),
    '"FOO" and "BAR", which are not in the fitted'
  )
})

test_that("from and to must differ", {
  expect_error(predict_realloc(lm_fit, delta = 10, from = "SB", to = "SB"),
               "must differ")
})

# --- Covariates: at, by --------------------------------------------------------

test_that("by = a three-level factor gives one row set per level", {
  fit <- comp_iv(cd, "zbmi", model = "lm", covariates = c("sex", "region"))
  levs <- levels(factor(timeuse$region))
  expect_length(levs, 3L)
  r <- predict_realloc(fit, delta = c(15, 30), from = "SB", to = "MVPA",
                       by = "region")
  expect_equal(nrow(r), 2 * length(levs))
  expect_setequal(r$region, levs)
  expect_equal(names(r)[1:4], c("from", "to", "delta", "region"))
})

test_that("at and by together pin one covariate and group by another", {
  fit <- comp_iv(cd, "zbmi", model = "lm", covariates = c("sex", "region"))
  both <- predict_realloc(fit, delta = 15, from = "SB", to = "MVPA",
                          at = list(region = "A"), by = "sex")
  expect_equal(nrow(both), nlevels(factor(timeuse$sex)))
  lev <- levels(factor(timeuse$sex))[1]
  manual <- predict_realloc(fit, delta = 15, from = "SB", to = "MVPA",
                            at = list(region = "A", sex = lev))
  expect_equal(both$estimate[both$sex == lev], manual$estimate,
               tolerance = 1e-10)
})

test_that("at and by cannot name the same covariate", {
  fit <- comp_iv(cd, "zbmi", model = "lm", covariates = c("sex", "ilr1:sex"))
  expect_error(
    predict_realloc(fit, delta = 15, from = "SB", to = "MVPA",
                    at = list(sex = "female"), by = "sex"),
    "`at` and `by` both name \"sex\""
  )
})

test_that("the moderator default for by yields only to an at that pins it", {
  mod <- comp_iv(cd, outcome = "zbmi", model = "lm", moderator = "sex",
                 covariates = "region")
  other <- predict_realloc(mod, delta = 15, from = "SB", to = "MVPA",
                           at = list(region = "A"))
  expect_identical(attr(other, "by"), "sex")
  expect_setequal(other$sex, levels(factor(timeuse$sex)))
  manual <- predict_realloc(mod, delta = 15, from = "SB", to = "MVPA",
                            at = list(region = "A", sex = "female"))
  expect_equal(other$estimate[other$sex == "female"], manual$estimate,
               tolerance = 1e-10)
  pinned <- predict_realloc(mod, delta = 15, from = "SB", to = "MVPA",
                            at = list(sex = "female"))
  expect_null(attr(pinned, "by"))
  expect_false("sex" %in% names(pinned))
})

test_that("an abort inside a by group keeps its own message", {
  fit <- comp_iv(cd, "zbmi", model = "lm", covariates = "sex")
  err <- expect_error(
    predict_realloc(fit, delta = 5000, from = "SB", to = "MVPA", by = "sex",
                    on_invalid = "error"),
    "makes a part non-positive"
  )
  expect_false(inherits(err, "purrr_error_indexed"))
})

test_that("an interaction without at names the interacting covariates", {
  fit <- comp_iv(cd, zbmi ~ ilr() * sex + ilr():region, model = "lm")
  expect_error(predict_realloc(fit, delta = 15, from = "SB", to = "MVPA"),
               "at = list\\(sex = \"female\", region = \"A\"\\)")
})

test_that("every interacting covariate must be pinned by at or by", {
  fit <- comp_iv(cd, zbmi ~ ilr() * sex + ilr():region, model = "lm")
  expect_error(predict_realloc(fit, delta = 15, from = "SB", to = "MVPA",
                               by = "sex"),
               "a value for \"region\".*at = list\\(region = \"A\"\\)")
  out <- predict_realloc(fit, delta = 15, from = "SB", to = "MVPA",
                         at = list(sex = "female", region = "A"))
  expect_equal(nrow(out), 1L)
})

test_that("an unpinned interacting covariate is named with a value of its own type", {
  d <- timeuse
  d$active <- d$zbmi > 0
  d$k <- as.integer(d$region)
  cd1 <- declare_composition(d, parts, closure = 1440)
  # A logical covariate gives its first value, a numeric column wrapped in
  # factor() its smallest value, and any other numeric its mean to three
  # significant digits.
  cases <- list(
    list(fit = comp_iv(cd1, "zbmi", model = "lm", moderator = "active"),
         text = "at = list(active = FALSE)", at = list(active = FALSE)),
    list(fit = comp_iv(cd1, zbmi ~ ilr() * factor(k), model = "lm"),
         text = "at = list(k = 1)", at = list(k = 1)),
    list(fit = comp_iv(cd1, "obesity", model = "glm", family = "binomial",
                       moderator = "zbmi"),
         # The mean lies on a rounding edge (-0.0228 or -0.0229 by platform),
         # so the expected text is rounded here as the message rounds it.
         text = paste0("at = list(zbmi = ",
                       deparse1(signif(mean(timeuse$zbmi), 3)), ")"),
         at = list(zbmi = signif(mean(timeuse$zbmi), 3)))
  )
  for (case in cases) {
    err <- expect_error(
      predict_realloc(case$fit, delta = 15, from = "SB", to = "MVPA"),
      "covariate-conditional"
    )
    expect_match(conditionMessage(err), case$text, fixed = TRUE)
    # The suggested `at` is accepted and gives the estimate at that value.
    out <- predict_realloc(case$fit, delta = 15, from = "SB", to = "MVPA",
                           at = case$at)
    expect_true(is.finite(out$estimate))
  }
  # The value matters: the other level of the logical gives another estimate.
  other <- predict_realloc(cases[[1]]$fit, delta = 15, from = "SB", to = "MVPA",
                           at = list(active = TRUE))
  first <- predict_realloc(cases[[1]]$fit, delta = 15, from = "SB", to = "MVPA",
                           at = list(active = FALSE))
  expect_false(isTRUE(all.equal(other$estimate, first$estimate)))
})

test_that("a reallocation of a moderated fit prints no emmeans note", {
  mod <- comp_iv(cd, outcome = "zbmi", model = "lm", moderator = "sex")
  expect_no_message(predict_realloc(mod, delta = 15, from = "SB", to = "MVPA"))
})

# --- Compositions that share a pivot coordinate --------------------------------

test_that("a zero reallocation is exactly the null contrast", {
  d0 <- predict_realloc(lm_fit, delta = 0, from = "SB", to = "MVPA")
  expect_identical(d0$estimate, 0)
  expect_identical(d0$std.error, 0)
  # The Wald test of a zero-SE contrast is undefined: NA, not NaN.
  expect_true(is.na(d0$statistic) && !is.nan(d0$statistic))
  expect_true(is.na(d0$p.value) && !is.nan(d0$p.value))
  r0 <- predict_realloc(glm_fit, delta = c(0, 15), from = "SB", to = "MVPA")
  expect_identical(r0$estimate[1], 1)
  expect_identical(r0$std.error[1], 0)
  expect_true(is.na(r0$statistic[1]) && !is.nan(r0$statistic[1]))
  expect_true(is.na(r0$p.value[1]) && !is.nan(r0$p.value[1]))
  expect_true(is.finite(r0$statistic[2]) && is.finite(r0$p.value[2]))
  expect_equal(r0$estimate[2],
               predict_realloc(glm_fit, delta = 15, from = "SB",
                               to = "MVPA")$estimate)
})

test_that("a zero-minute row is the reference value with no test, on every engine", {
  lmer_fit <- comp_iv(cd, "zbmi", model = "lmer", random = "school")
  r <- predict_realloc(lmer_fit, delta = c(0, 15), from = "SB", to = "MVPA")
  expect_equal(nrow(r), 2L)
  expect_identical(r$estimate[1], 0)
  expect_identical(r$std.error[1], 0)
  expect_true(is.na(r$df[1]))
  expect_true(is.finite(r$estimate[2]) && r$estimate[2] != 0)
  g0 <- predict_realloc(glm_fit, delta = 0, from = "SB", to = "MVPA")
  expect_identical(g0$estimate, 1)
  expect_identical(g0$std.error, 0)
  expect_identical(c(g0$conf.low, g0$conf.high), c(1, 1))
  expect_true(is.na(g0$df))
})

test_that("a repeated delta and two equal reference parts are contrasted", {
  rep_d <- predict_realloc(glm_fit, delta = c(15, 15), from = "SB", to = "MVPA")
  expect_equal(nrow(rep_d), 2L)
  expect_identical(rep_d$estimate[1], rep_d$estimate[2])
  eq_ref <- c(MVPA = 60, LPA = 300, SB = 300, SPT = 780)
  all_pairs <- predict_realloc(lm_fit, delta = 15, reference = eq_ref)
  expect_equal(nrow(all_pairs), 12L)
  expect_false(anyNA(all_pairs$estimate))
  one <- predict_realloc(lm_fit, delta = 15, from = "LPA", to = "SB",
                         reference = eq_ref)
  expect_equal(all_pairs$estimate[all_pairs$from == "LPA" &
                                    all_pairs$to == "SB"],
               one$estimate, tolerance = 1e-12)
})

# --- Absolute output with an infeasible row ----------------------------------

test_that("absolute mode with on_invalid = 'na' keeps the reference row", {
  out <- predict_realloc(lm_fit, delta = c(20, 1e6), from = "SB", to = "MVPA",
                         type = "absolute", on_invalid = "na")
  expect_equal(nrow(out), 3L)
  expect_equal(out$delta, c(0, 20, 1e6))
  expect_false(is.na(out$estimate[1]))
  expect_false(is.na(out$estimate[2]))
  expect_true(is.na(out$estimate[3]) && is.na(out$conf.low[3]))
  valid_only <- predict_realloc(lm_fit, delta = 20, from = "SB", to = "MVPA",
                                type = "absolute")
  expect_equal(out$estimate[1:2], valid_only$estimate)
})

# --- Cox: a namespaced survival response ---------------------------------------

test_that("a namespaced survival::Surv() response is not back-transformed twice", {
  d <- cox_fit$data
  comps <- rbind(unlist(cox_fit$ref_comp),
                 unlist(cox_fit$ref_comp) + c(10, 0, -10, 0))
  relative <- function(model) {
    f <- cox_fit
    f$model <- model
    em <- tc_emm_at_points(f, comps, at = NULL, conf_level = 0.95)
    summary(emmeans::contrast(em, method = "trt.vs.ctrl", ref = 1),
            infer = c(TRUE, TRUE), type = "response", adjust = "none")
  }
  namespaced <- survival::coxph(
    survival::Surv(time, event) ~ ilr1 + ilr2 + ilr3, data = d
  )
  bare <- local({
    Surv <- survival::Surv
    survival::coxph(Surv(time, event) ~ ilr1 + ilr2 + ilr3, data = d)
  })
  expect_no_error(ns <- relative(namespaced))
  br <- relative(bare)
  expect_equal(ns$ratio, br$ratio, tolerance = 1e-12)
  expect_equal(ns$asymp.LCL, br$asymp.LCL, tolerance = 1e-12)
  expect_equal(ns$p.value, br$p.value, tolerance = 1e-12)
})

test_that("a declared but unused factor level is skipped by the by grouping", {
  d <- timeuse
  d$sex <- factor(d$sex, levels = c("female", "male", "other"))
  cd3 <- declare_composition(d, parts, closure = 1440)
  fit_by <- comp_iv(cd3, "obesity", model = "glm", family = "binomial",
                    covariates = "sex")
  r <- predict_realloc(fit_by, delta = 30, from = "SB", to = "MVPA", by = "sex")
  expect_setequal(r$sex, c("female", "male"))
  fit_mod <- comp_iv(cd3, "obesity", model = "glm", family = "binomial",
                     moderator = "sex")
  r_mod <- predict_realloc(fit_mod, delta = 30, from = "SB", to = "MVPA")
  expect_setequal(r_mod$sex, c("female", "male"))
})

test_that("by must name a covariate of the model, not any data column", {
  d <- timeuse
  d$obesity_f <- factor(d$obesity)
  cd3 <- declare_composition(d, parts, closure = 1440)
  fit <- comp_iv(cd3, "obesity_f", model = "glm", family = "binomial",
                 covariates = "sex")
  expect_error(
    predict_realloc(fit, delta = 30, from = "SB", to = "MVPA", by = "obesity_f"),
    "not a covariate in the fitted model"
  )
})

test_that("a missing delta is reported by predict_realloc() itself", {
  fit <- comp_iv(cd, "obesity", model = "glm", family = "binomial")
  expect_error(predict_realloc(fit, from = "SB", to = "MVPA"),
               "`delta` is missing")
})

test_that("by accepts a covariate that enters only inside an interaction", {
  fit <- comp_iv(cd, zbmi ~ ilr() * sex + ilr():region, model = "lm")
  r <- predict_realloc(fit, 15, "SB", "MVPA", at = list(sex = "female"),
                       by = "region")
  expect_setequal(unique(r$region), unique(as.character(timeuse$region)))
  expect_equal(nrow(r), length(unique(timeuse$region)))
})

test_that("predict_realloc() refuses a conf_level outside (0, 1)", {
  for (bad in list(95, 0, 1, NA_real_, c(0.9, 0.95), "0.95")) {
    expect_error(predict_realloc(lm_fit, 10, "SB", "MVPA", conf_level = bad),
                 "conf_level.*single number in \\(0, 1\\)")
  }
})

test_that("the emmeans column readers abort on a summary without the column", {
  s <- data.frame(contrast = "a", SE = 1)
  expect_error(tc_estimate_col(s), class = "rlang_error")
  expect_error(tc_ci_col(s, "lower"), class = "rlang_error")
  expect_error(tc_ci_col(s, "upper"), class = "rlang_error")
})
