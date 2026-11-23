parts <- c("MVPA", "LPA", "SB", "SPT")
cd <- declare_composition(timeuse, parts, closure = 1440)
glm_fit <- comp_iv(cd, "obesity", model = "glm", family = "binomial")

test_that("summary.realloc returns a visible summary_realloc with a tidy table", {
  r <- predict_realloc(glm_fit, delta = c(10, 20, 30), on_invalid = "na")
  expect_visible(summary(r))
  s <- summary(r)
  expect_s3_class(s, "summary_realloc")
  expect_named(s, c("table", "effect_label", "type", "conf_level", "parts",
                    "by", "units"))
  expect_s3_class(s$table, "tbl_df")
  expect_named(s$table, c("to", "from", "to_label", "from_label", "delta",
                          "estimate", "conf.low", "conf.high", "p.value"))
  # One row per reallocation: 12 ordered pairs x 3 deltas.
  expect_equal(nrow(s$table), 12 * 3)
  expect_equal(s$effect_label, "Odds ratio")
  expect_equal(s$type, "relative")
})

test_that("the summary table keeps full numeric precision", {
  r <- predict_realloc(glm_fit, delta = 30, from = "SB", to = "MVPA")
  s <- summary(r)
  expect_type(s$table$estimate, "double")
  expect_identical(s$table$estimate, r$estimate)
  expect_identical(s$table$conf.low, r$conf.low)
  expect_identical(s$table$p.value, r$p.value)
})

test_that("printing the summary shows a text report and the table", {
  r <- predict_realloc(glm_fit, delta = c(10, 30), from = "SB", to = "MVPA")
  s <- summary(r)
  out <- capture.output(print(s))
  expect_match(out[1], "Odds ratio versus the reference composition")
  expect_true(any(grepl("2 reallocations of 10, 30; 95% confidence", out)))
  capture.output(shown <- withVisible(print(s)))
  expect_false(shown$visible)
  expect_identical(shown$value, s)
  one <- capture.output(print(summary(predict_realloc(glm_fit, delta = 10,
                                                      from = "SB", to = "MVPA"))))
  expect_true(any(grepl("1 reallocation of 10;", one)))
})

test_that("the report counts infeasible reallocations and names the grouping", {
  r <- predict_realloc(glm_fit, delta = c(10, 1e6), from = "SB", to = "MVPA",
                       on_invalid = "na")
  out <- capture.output(print(summary(r)))
  expect_true(any(grepl("No estimate (a part driven non-positive): 1", out,
                        fixed = TRUE)))

  mod <- comp_iv(cd, outcome = "obesity", moderator = "sex", model = "glm",
                 family = "binomial")
  rb <- predict_realloc(mod, from = "SB", to = "MVPA", delta = 15)
  outb <- capture.output(print(summary(rb)))
  expect_true(any(grepl("Estimated per level of sex", outb, fixed = TRUE)))
})

test_that("an absolute summary drops the reference row", {
  lm_fit <- comp_iv(cd, "zbmi", model = "lm")
  r <- predict_realloc(lm_fit, delta = c(15, 30), from = "SB", to = "MVPA",
                       type = "absolute")
  s <- summary(r)
  expect_equal(nrow(s$table), 2L)
  expect_false("p.value" %in% names(s$table))
  expect_match(capture.output(print(s))[1], "at each reallocated composition")
})

test_that("labels relabel parts and the implicit rest", {
  r <- predict_realloc(glm_fit, delta = 30, scheme = "one-to-all", to = "MVPA",
                       on_invalid = "na")
  w <- summary(r, labels = c(MVPA = "Exercise"))$table
  expect_true(all(w$to_label == "Exercise"))
  expect_true(all(w$from_label == "All other parts"))
  expect_true(all(w$from == tc_rest_key()))
})

test_that("labels default to the compD labels carried by the fit", {
  cd2 <- declare_composition(timeuse, parts, closure = 1440,
                             labels = c(MVPA = "Exercise"))
  fit2 <- comp_iv(cd2, "obesity", model = "glm", family = "binomial")
  r <- predict_realloc(fit2, delta = 30, scheme = "one-to-all", to = "MVPA",
                       on_invalid = "na")
  w <- summary(r)$table
  expect_true(all(w$to_label == "Exercise"))
})

test_that("one-to-all summary yields one row per focal part and delta", {
  r <- predict_realloc(glm_fit, delta = c(10, 20, 30), scheme = "one-to-all",
                       on_invalid = "na")
  w <- summary(r)$table
  expect_equal(nrow(w), length(parts) * 3)
  expect_equal(unique(w$to), parts)
  expect_true(all(w$from_label == "All other parts"))
})

test_that("compound labels substitute whole part names only, once each", {
  labs <- c(MVPA = "Vigorous PA", PA = "Light activity")
  expect_identical(tc_relabel_compound("MVPA +10, PA -10", labs),
                   "Vigorous PA +10, Light activity -10")
  expect_identical(tc_relabel_compound("A.B +5, AXB -5", c("A.B" = "dotted")),
                   "dotted +5, AXB -5")
  expect_identical(tc_relabel_compound("SB -15, LPA -15", c(SPT = "Sleep")),
                   "SB -15, LPA -15")
})

test_that("summary.realloc rejects bad labels and non-realloc input", {
  r <- predict_realloc(glm_fit, delta = 30, from = "SB", to = "MVPA")
  expect_error(summary(r, labels = c(XYZ = "x")), "not a part")
  expect_error(summary(r, labels = c(XYZ = "x", ABC = "y")), "which are not")
  expect_error(summary.realloc(list()), "realloc")
})

test_that("summary.realloc refuses display, digits and p_flag, pointing to pub_ready()", {
  r <- predict_realloc(glm_fit, delta = 30, from = "SB", to = "MVPA")
  one <- expect_error(summary(r, display = FALSE), "pub_ready")
  expect_match(conditionMessage(one), "does not take `display`.", fixed = TRUE)
  flag <- expect_error(summary(r, p_flag = 0.1), "pub_ready")
  expect_match(conditionMessage(flag), "does not take `p_flag`.", fixed = TRUE)
  two <- expect_error(summary(r, display = FALSE, digits = 3), "pub_ready")
  expect_match(conditionMessage(two), "does not take `display` and `digits`.",
               fixed = TRUE)
  three <- expect_error(summary(r, p_flag = 0.1, digits = 3, display = FALSE),
                        "pub_ready")
  expect_match(conditionMessage(three),
               "does not take `display`, `digits`, and `p_flag`.",
               fixed = TRUE)
  # The confidence level is set where the estimates are made.
  expect_error(summary(r, conf_level = 0.9), "predict_realloc(conf_level = )",
               fixed = TRUE)
})

# --- pub_ready() ------------------------------------------------------------

test_that("pub_ready() builds the wide effects table from the summary", {
  skip_if_not_installed("gt")
  r <- predict_realloc(glm_fit, delta = c(10, 20, 30), on_invalid = "na")
  g <- pub_ready(summary(r))
  expect_named(g, "effects")
  expect_s3_class(g$effects, "gt_tbl")
  d <- g$effects[["_data"]]
  expect_true(all(c("to_label", "from_label", "est_10", "ci_10", "est_20",
                    "ci_20", "est_30", "ci_30") %in% names(d)))
  # One row per ordered pair; the across-delta "vs" columns are gone.
  expect_equal(nrow(d), 12L)
  expect_false(any(grepl("_vs_", names(d))))
})

test_that("pub_ready() columns and spanners derive from the deltas present", {
  skip_if_not_installed("gt")
  r <- predict_realloc(glm_fit, delta = c(15, 45), on_invalid = "na")
  d <- pub_ready(summary(r))$effects[["_data"]]
  expect_true(all(c("est_15", "ci_15", "est_45", "ci_45") %in% names(d)))
  expect_false(any(grepl("est_30", names(d))))
})

test_that("pub_ready() digits sets the formatted precision", {
  skip_if_not_installed("gt")
  r <- predict_realloc(glm_fit, delta = 30, from = "SB", to = "MVPA")
  d <- pub_ready(summary(r), digits = 4)$effects[["_data"]]
  # Four decimals in the formatted estimate (plus an optional trailing *).
  expect_match(d$est_30[1], "^[0-9]+\\.[0-9]{4}\\*?$")
  expect_match(d$ci_30[1], "^[0-9]+\\.[0-9]{4}, [0-9]+\\.[0-9]{4}$")
})

test_that("the pub_ready() significance flag follows its p_flag", {
  skip_if_not_installed("gt")
  r <- predict_realloc(glm_fit, delta = c(10, 20, 30), on_invalid = "na")
  s <- summary(r)
  flagged <- function(p_flag) {
    d <- pub_ready(s, p_flag = p_flag)$effects[["_data"]]
    est <- unlist(d[grep("^est_", names(d))])
    sum(grepl("\\*", est))
  }
  expect_gt(flagged(0.999), flagged(1e-12))
  expect_equal(flagged(0), 0)
  expect_equal(pub_ready(s)$effects[["_data"]],
               pub_ready(s, p_flag = 0.05)$effects[["_data"]])
  note <- pub_ready(s, p_flag = 0.1)$effects[["_source_notes"]][[1]]
  expect_match(as.character(note), "* p < 0.1", fixed = TRUE)
})

test_that("the pub_ready() subtitle matches the printed report", {
  skip_if_not_installed("gt")
  lm_fit <- comp_iv(cd, "zbmi", model = "lm")
  ab <- predict_realloc(lm_fit, delta = 15, from = "SB", to = "MVPA",
                        type = "absolute")
  expect_identical(pub_ready(summary(ab))$effects[["_heading"]]$subtitle,
                   "Predicted response at each reallocated composition")
  rel <- predict_realloc(glm_fit, delta = 15, from = "SB", to = "MVPA")
  expect_identical(pub_ready(summary(rel))$effects[["_heading"]]$subtitle,
                   "Odds ratio versus the reference composition")
})

test_that("pub_ready() keeps the group and response_level row keys", {
  skip_if_not_installed("gt")
  mod <- comp_iv(cd, outcome = "obesity", moderator = "sex", model = "glm",
                 family = "binomial")
  r <- predict_realloc(mod, from = "SB", to = "MVPA", delta = c(10, 20),
                       on_invalid = "na")
  s <- summary(r)
  expect_true("sex" %in% names(s$table))
  d <- pub_ready(s)$effects[["_data"]]
  expect_true("sex" %in% names(d))
  expect_equal(nrow(d), nlevels(timeuse$sex))
})

# --- Multinomial end-to-end (summary.comp_iv + summary.realloc) ---------------

test_that("multinomial end-to-end: summary(fit), predict_realloc, and summary(realloc) keep the response_level distinction", {
  skip_if_not_installed("gt")
  skip_if_not_installed("mclogit")
  small_parts <- c("MVPA", "SB", "SPT")
  d <- timeuse
  d$grade_cat <- factor(as.character(d$grade), levels = levels(d$grade))
  cd_multi <- declare_composition(d, small_parts, closure = 1440)
  fit <- suppressWarnings(
    comp_iv(cd_multi, "grade_cat", model = "mblogit")
  )
  k <- nlevels(d$grade_cat)

  # summary(fit)'s per-part table carries response_level.
  s_fit <- summary(fit)
  expect_true("response_level" %in% names(s_fit$per_part))
  expect_equal(nrow(s_fit$per_part), length(small_parts) * (k - 1))

  # predict_realloc() output carries response_level, one row-set per
  # non-reference category.
  rr <- suppressWarnings(
    predict_realloc(fit, delta = c(10, 20), from = "SB", to = "MVPA",
                    on_invalid = "na")
  )
  expect_true("response_level" %in% names(rr))
  expect_equal(nrow(rr), 2 * (k - 1))
  expect_setequal(rr$response_level, levels(d$grade_cat)[-1])

  # summary(realloc) and its publication table keep the category distinction.
  s_rr <- summary(rr)
  expect_true("response_level" %in% names(s_rr$table))
  expect_equal(nrow(s_rr$table), 2 * (k - 1))
  g <- pub_ready(s_rr)$effects
  expect_s3_class(g, "gt_tbl")
  expect_true("response_level" %in% names(g[["_data"]]))
  expect_equal(nrow(g[["_data"]]), k - 1)
})

test_that("the declared unit label travels to the summary and the pub_ready() spanners", {
  cdu <- declare_composition(timeuse, parts, closure = 1440, units = "min/day")
  fit_u <- comp_iv(cdu, "obesity", model = "glm", family = "binomial")
  r <- predict_realloc(fit_u, delta = c(10, 30), from = "SB", to = "MVPA")
  r0 <- predict_realloc(glm_fit, delta = c(10, 30), from = "SB", to = "MVPA")
  expect_identical(attr(r, "units"), "min/day")
  expect_identical(summary(r)$units, "min/day")
  expect_null(attr(r0, "units"))
  expect_match(capture.output(print(summary(r)))[2], "of 10, 30 min/day;", fixed = TRUE)
  expect_match(capture.output(print(summary(r0)))[2], "of 10, 30;", fixed = TRUE)
  skip_if_not_installed("gt")
  spanners <- function(x) unlist(pub_ready(summary(x))$effects[["_spanners"]]$spanner_label)
  expect_identical(spanners(r), c("10 min/day", "30 min/day"))
  expect_identical(spanners(r0), c("10", "30"))
})

test_that("summary keeps a zero reallocation the user asked for", {
  lm_fit <- comp_iv(cd, "zbmi", model = "lm")
  r <- predict_realloc(lm_fit, c(0, 10), "SB", "MVPA")
  s <- summary(r)
  expect_equal(s$table$delta, c(0, 10))
  expect_identical(s$table$estimate, r$estimate)
  # An absolute prediction drops its reference row and keeps the zero move.
  ra <- predict_realloc(lm_fit, c(0, 10), "SB", "MVPA", type = "absolute")
  expect_equal(summary(ra)$table$delta, c(0, 10))
})
