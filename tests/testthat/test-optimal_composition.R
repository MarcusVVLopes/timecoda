parts <- c("MVPA", "LPA", "SB", "SPT")
cd <- declare_composition(timeuse, parts, closure = 1440)
quad <- comp_iv(cd, formula = zbmi ~ ilr()^2 + I(ilr()^2) + sex, model = "lm")
risk <- comp_iv(cd, formula = obesity ~ ilr()^2 + I(ilr()^2), model = "glm",
                family = "binomial")
lin <- comp_iv(cd, "zbmi", model = "lm", covariates = "sex")

count_messages <- function(expr) {
  n <- 0L
  withCallingHandlers(expr, message = function(m) {
    n <<- n + 1L
    invokeRestart("muffleMessage")
  })
  n
}

# The printed report as one string, its lines joined and each run of blanks
# made one, so that a phrase is found wherever the console width wrapped it.
printed <- function(x) {
  gsub("\\s+", " ", paste(utils::capture.output(print(x)), collapse = " "))
}

# Every composition on multiples of `spacing` (the last part the remainder),
# over the whole simplex.
simplex_grid <- function(spacing, nm = parts, total = 1440) {
  d <- length(nm)
  v <- seq(spacing, total - (d - 1) * spacing, by = spacing)
  g <- as.matrix(expand.grid(rep(list(v), d - 1)))
  g <- cbind(g, total - rowSums(g))
  g <- g[g[, d] >= spacing, , drop = FALSE]
  colnames(g) <- nm
  g
}

# The search's estimates at 20 of its candidates against the per-point route,
# one emmeans grid per composition: the link, less the link at `from` when the
# search is relative.
expect_per_point <- function(opt, fit, at = NULL, category = NULL, n = 20) {
  cand <- opt$candidates
  idx <- unique(round(seq(1, nrow(cand), length.out = n)))
  link <- function(comps) {
    stats::predict(tc_emm_at_points(fit, comps, at = at, conf_level = 0.95,
                                    category = category))
  }
  ref <- link(as.matrix(cand[idx, opt$parts]))
  if (opt$type == "relative") {
    ref <- ref - link(t(unlist(fit$ref_comp)[opt$parts]))
  }
  expect_equal(cand[[paste0("estimate_", opt$outcomes[1])]][idx], ref,
               tolerance = 1e-10)
}

test_that("the result carries the documented elements", {
  opt <- optimal_composition(quad, spacing = 30)
  expect_s3_class(opt, "optimal_comp")
  expect_named(opt, c(
    "candidates", "zones", "optimum", "fits", "direction", "zone",
    "threshold", "widen", "type", "category", "region", "coverage", "sd",
    "spacing", "n_candidates", "zone_used", "parts", "labels", "units", "total",
    "outcomes", "linear", "measures", "by", "at", "conf_level"
  ))
  expect_false(opt$widen)
  expect_identical(opt$zone_used, c(all = 0.95))
  expect_true(all(c(parts, "ilr1", "ilr2", "ilr3", "mahalanobis",
                    "estimate_zbmi", "rank_zbmi", "in_zone_zbmi", "in_zone") %in%
                    names(opt$candidates)))
  expect_equal(rowSums(as.matrix(opt$candidates[parts])),
               rep(1440, nrow(opt$candidates)))
  expect_true(all(is.na(opt$candidates$mahalanobis)))
  expect_identical(opt$n_candidates[["all"]], nrow(opt$candidates))
  expect_identical(opt$zones$outcome, c("zbmi", "overall"))
  expect_identical(opt$zones$n[1], sum(opt$candidates$in_zone_zbmi))
  expect_identical(opt$optimum$centre_of, "overall")
  expect_identical(opt$measures, c(zbmi = "Mean difference"))
  expect_false(opt$linear[["zbmi"]])
})

test_that("the default call is the box of the mean +/- 3 SD on a grid of total / 144", {
  opt <- optimal_composition(quad, direction = "min")
  expect_identical(opt$region$name, "box")
  expect_identical(opt$sd, 3)
  expect_null(opt$coverage)
  expect_identical(opt$spacing, 1440 / 144)
  cand <- as.matrix(opt$candidates[parts])
  expect_true(all(abs(cand / 10 - round(cand / 10)) < 1e-9))
  expect_output(print(opt), "Region +per-part box of the mean \\+/- 3 SD\n  Grid +multiples of 10; ")
})

test_that("the total and the reference day come from the fit's composition", {
  for (closure in c(24, 1)) {
    cd_k <- declare_composition(timeuse, parts, closure = closure,
                                units = if (closure == 24) "h/day")
    fit <- comp_iv(cd_k, formula = zbmi ~ ilr()^2 + I(ilr()^2) + sex, model = "lm")
    opt <- optimal_composition(fit, direction = "min")
    expect_identical(opt$total, fit$comp$closure)
    expect_equal(opt$spacing, fit$comp$closure / 144)
    expect_equal(rowSums(as.matrix(opt$candidates[parts])),
                 rep(closure, nrow(opt$candidates)), tolerance = 1e-12)
    ref <- geom_mean_comp(tc_comp_df(fit$comp), parts, total = fit$comp$closure)
    expect_equal(unlist(opt$optimum[paste0("from_", parts)]), unlist(ref),
                 ignore_attr = TRUE)
    expect_equal(sum(unlist(opt$optimum[parts])), closure)
    overall <- opt$zones[opt$zones$outcome == "overall", ]
    expect_equal(sum(unlist(overall[paste0("centre_", parts)])), closure)
  }
  # A search split by `by` measures every level from the same day.
  mod <- comp_iv(cd, "zbmi", model = "lm", moderator = "sex", curvature = TRUE)
  opt <- suppressMessages(optimal_composition(mod, spacing = 40))
  from <- as.matrix(opt$optimum[paste0("from_", parts)])
  expect_identical(nrow(from), 2L)
  expect_equal(from[1, ], from[2, ])
  expect_equal(from[1, ], unlist(geom_mean_comp(tc_comp_df(cd), parts, total = 1440)),
               ignore_attr = TRUE)
})

test_that("the box search equals a hand enumeration of Dumuid et al. (2021)", {
  comp <- as.matrix(tc_comp_df(cd)[parts])
  lo <- pmax(colMeans(comp) - 3 * apply(comp, 2, stats::sd), apply(comp, 2, min))
  hi <- pmin(colMeans(comp) + 3 * apply(comp, 2, stats::sd), apply(comp, 2, max))
  grid <- simplex_grid(30)
  grid <- grid[rowSums(sweep(grid, 2, lo, ">=") & sweep(grid, 2, hi, "<=")) == 4, ]
  # The link by stats::predict(), sex averaged with equal weight.
  z <- as.data.frame(tc_pivot_coords(grid, pivot = quad$pivot))
  by_sex <- purrr::map(levels(timeuse$sex), function(s) {
    stats::predict(quad$model, newdata = cbind(z, sex = factor(s, levels(timeuse$sex))))
  })
  link <- Reduce(`+`, by_sex) / length(by_sex)
  in_zone <- rank(-link) / length(link) >= 0.95
  opt <- optimal_composition(quad, direction = "min", spacing = 30)
  key <- function(m) apply(round(m, 6), 1, paste, collapse = "|")
  expect_setequal(key(as.matrix(opt$candidates[parts])), key(grid))
  expect_setequal(key(as.matrix(opt$candidates[opt$candidates$in_zone, parts])),
                  key(grid[in_zone, ]))
  centre <- exp(colMeans(log(grid[in_zone, ])))
  expect_equal(unlist(opt$optimum[parts]), centre / sum(centre) * 1440,
               ignore_attr = TRUE, tolerance = 1e-8)
})

test_that("the ellipsoid is cut at the chi-square quantile", {
  opt <- optimal_composition(quad, region = "ellipsoid", spacing = 30)
  z <- as.matrix(ilr_transform(tc_comp_df(cd), parts)[c("ilr1", "ilr2", "ilr3")])
  par <- opt$region$parameters$all
  expect_equal(par$mu, colMeans(z), ignore_attr = TRUE, tolerance = 1e-12)
  expect_equal(par$S, stats::cov(z), ignore_attr = TRUE, tolerance = 1e-12)
  expect_identical(par$r2, stats::qchisq(0.80, df = 3))
  expect_equal(par$share,
               mean(stats::mahalanobis(z, colMeans(z), stats::cov(z)) <= par$r2))
  expect_identical(opt$coverage, 0.80)
  expect_null(opt$sd)
  expect_true(all(opt$candidates$mahalanobis <= par$r2))
  expect_match(printed(opt), sprintf("80%% chi-square level, holding %.1f%% of the sample",
                                     100 * par$share), fixed = TRUE)
  # The bounded enumeration equals the whole simplex filtered by the ellipsoid.
  grid <- simplex_grid(30)
  zg <- as.matrix(ilr_transform(as.data.frame(grid), parts)[c("ilr1", "ilr2", "ilr3")])
  inside <- stats::mahalanobis(zg, par$mu, par$S) <= par$r2
  key <- function(m) apply(round(m, 6), 1, paste, collapse = "|")
  expect_setequal(key(as.matrix(opt$candidates[parts])), key(grid[inside, ]))
})

test_that("the search works for three and five parts", {
  cd3 <- declare_composition(timeuse, c("MVPA", "SB", "SPT"), closure = 1440)
  fit3 <- comp_iv(cd3, formula = zbmi ~ ilr()^2 + I(ilr()^2) + sex, model = "lm")
  opt3 <- optimal_composition(fit3, spacing = 20)
  expect_true(all(c("MVPA", "SB", "SPT", "ilr1", "ilr2") %in% names(opt3$candidates)))
  expect_false("ilr3" %in% names(opt3$candidates))
  expect_equal(sum(unlist(opt3$optimum[c("MVPA", "SB", "SPT")])), 1440)
  expect_per_point(opt3, fit3)

  set.seed(15)
  t5 <- timeuse
  share <- pmin(pmax(0.6 + stats::rnorm(nrow(t5), 0, 0.05), 0.4), 0.8)
  t5$SB1 <- t5$SB * share
  t5$SB2 <- t5$SB * (1 - share)
  parts5 <- c("MVPA", "LPA", "SB1", "SB2", "SPT")
  cd5 <- declare_composition(t5, parts5, closure = 1440)
  opt5 <- suppressMessages(optimal_composition(comp_iv(cd5, "zbmi", model = "lm"),
                                               spacing = 60))
  expect_true(all(c(parts5, paste0("ilr", 1:4)) %in% names(opt5$candidates)))
  expect_length(grep("^centre_", names(opt5$zones)), 5)
})

test_that("the tied grids equal the per-point predictions", {
  expect_per_point(optimal_composition(quad, spacing = 30, type = "absolute"), quad)
  expect_per_point(optimal_composition(risk, spacing = 30, type = "absolute"), risk)
  cox <- comp_iv(cd, formula = survival::Surv(time, event) ~ ilr()^2 +
                   I(ilr()^2) + sex, model = "cox")
  expect_per_point(optimal_composition(cox, spacing = 30), cox)
  # The total covariate, held at its reference value.
  tv <- timeuse
  set.seed(2026)
  tv[parts] <- tv[parts] * stats::runif(nrow(tv), 0.8, 1.2)
  tot <- comp_iv(declare_composition(tv, parts, closure = "mean"), "zbmi",
                 model = "lm", covariates = "sex", include_total = TRUE)
  expect_per_point(suppressMessages(
    optimal_composition(tot, spacing = 30, type = "absolute")
  ), tot)
  # A pinned covariate, a numeric covariate held at its mean, and a moderator
  # split by `by`.
  expect_per_point(optimal_composition(quad, spacing = 30, at = list(sex = "male")),
                   quad, at = list(sex = "male"))
  num <- comp_iv(cd, formula = zbmi ~ ilr()^2 + I(ilr()^2) + time + sex, model = "lm")
  expect_per_point(optimal_composition(num, spacing = 30, type = "absolute"), num)
  mod <- comp_iv(cd, "zbmi", model = "lm", moderator = "sex", curvature = TRUE)
  by_sex <- optimal_composition(mod, by = "sex", spacing = 30)
  for (lev in c("female", "male")) {
    one <- by_sex
    one$candidates <- by_sex$candidates[by_sex$candidates$sex == lev, ]
    expect_per_point(one, mod, at = list(sex = lev))
  }
  mixed <- comp_iv(cd, "zbmi", model = "lmer", random = "school",
                   covariates = "sex", curvature = TRUE)
  expect_per_point(optimal_composition(mixed, spacing = 30, type = "absolute"),
                   mixed)
  skip_if_not_installed("robustbase")
  set.seed(1)
  rob <- comp_iv(cd, formula = zbmi ~ ilr()^2 + I(ilr()^2) + sex, model = "lmrob")
  expect_per_point(optimal_composition(rob, spacing = 30, type = "absolute"), rob)
  skip_if_not_installed("mclogit")
  d <- timeuse
  d$grade_unordered <- factor(as.character(d$grade), levels = levels(d$grade))
  mb <- suppressWarnings(comp_iv(declare_composition(d, parts, closure = 1440),
                                 "grade_unordered", model = "mblogit",
                                 covariates = "sex"))
  levs <- unique(sub("~.*$", "", names(stats::coef(mb$model))))
  expect_per_point(suppressMessages(
    optimal_composition(mb, spacing = 30, category = levs[1])
  ), mb, category = levs[1])
})

test_that("repeated days and days sharing a coordinate are predicted per day", {
  # The first two rows are one day; the last two share `ilr1` (the same MVPA,
  # the other three parts permuted).
  days <- data.frame(MVPA = c(40, 40, 30, 30, 90), LPA = c(240, 240, 240, 510, 200),
                     SB = c(660, 660, 660, 660, 600), SPT = c(500, 500, 510, 240, 550))
  z <- tc_pivot_coords(days)
  expect_equal(z[3, "ilr1"], z[4, "ilr1"])
  opt <- optimal_composition(quad, region = days, zone = 0.5, type = "absolute")
  expect_per_point(opt, quad)
  expect_equal(opt$candidates$estimate_zbmi[1], opt$candidates$estimate_zbmi[2],
               tolerance = 1e-12)
})

test_that("emmeans ties a coordinate to the grid through a cov.reduce offset formula", {
  set.seed(3)
  dat <- data.frame(ilr1 = stats::rnorm(80), ilr2 = stats::rnorm(80))
  dat$y <- 1 + dat$ilr1 - 2 * dat$ilr2 + dat$ilr2^2 + stats::rnorm(80, 0, 0.1)
  m <- stats::lm(y ~ ilr1 + ilr2 + I(ilr2^2), data = dat)
  pts <- cbind(ilr1 = c(-1, 0.2, 1.5), ilr2 = c(0.7, -0.4, 1.1))
  lookup <- function(v, j) {
    out <- pts[match(round(v, 12), round(pts[, 1], 12)), j]
    out[is.na(out)] <- 0
    out
  }
  tie <- stats::as.formula("ilr2 ~ 0 + offset(lookup(ilr1, 2))",
                           env = environment())
  rg <- emmeans::ref_grid(m, at = list(ilr1 = pts[, 1]), cov.reduce = list(tie))
  em <- emmeans::emmeans(rg, ~ ilr1)
  expect_equal(stats::predict(em), unname(stats::predict(m, as.data.frame(pts))),
               tolerance = 1e-12)
})

test_that("the zone centre is the inverse ILR of the zone's mean coordinates", {
  opt <- optimal_composition(quad, spacing = 30)
  zone <- opt$candidates[opt$candidates$in_zone, ]
  hand <- ilr_backtransform(
    matrix(colMeans(zone[c("ilr1", "ilr2", "ilr3")]), nrow = 1), parts,
    total = 1440
  )
  expect_equal(unlist(opt$zones[2, paste0("centre_", parts)]), unlist(hand),
               ignore_attr = TRUE, tolerance = 1e-8)
  expect_equal(unlist(opt$optimum[parts]), unlist(hand), ignore_attr = TRUE,
               tolerance = 1e-8)
  expect_equal(unlist(opt$zones[2, paste0("min_", parts)]),
               purrr::map_dbl(zone[parts], min), ignore_attr = TRUE)
})

test_that("the gain is the emmeans contrast from the reference composition", {
  opt <- optimal_composition(quad, spacing = 30, type = "absolute")
  centre <- unlist(opt$optimum[parts])
  emm <- tc_emm_at_points(quad, rbind(unlist(quad$ref_comp), centre), at = NULL,
                          conf_level = 0.95)
  g <- summary(emmeans::contrast(emm, method = "trt.vs.ctrl", ref = 1),
               infer = c(TRUE, TRUE), adjust = "none")
  expect_equal(opt$optimum$gain, g$estimate, tolerance = 1e-10)
  expect_equal(opt$optimum$p.value, g$p.value, tolerance = 1e-10)
  expect_equal(unlist(opt$optimum[paste0("from_", parts)]), unlist(quad$ref_comp),
               ignore_attr = TRUE)

  day <- c(MVPA = 30, LPA = 240, SB = 660, SPT = 500)
  opt_day <- optimal_composition(quad, spacing = 30, from = day,
                                 type = "absolute")
  expect_equal(unlist(opt_day$optimum[paste0("from_", parts)]),
               day / sum(day) * 1440, ignore_attr = TRUE)
  expect_false(isTRUE(all.equal(opt_day$optimum$gain, opt$optimum$gain)))

  at_centre <- optimal_composition(quad, spacing = 30, from = centre,
                                   type = "absolute")
  expect_identical(at_centre$optimum$gain, 0)
  expect_identical(c(at_centre$optimum$gain.low, at_centre$optimum$gain.high), c(0, 0))
  expect_identical(at_centre$optimum$p.value, NA_real_)
  expect_equal(at_centre$optimum$estimate, opt$optimum$estimate, tolerance = 1e-10)
})

test_that("each fit's gain is measured from its own reference", {
  cd60 <- declare_composition(timeuse[1:60, ], parts, closure = 1440)
  fit60 <- comp_iv(cd60, formula = zbmi ~ ilr()^2 + I(ilr()^2), model = "lm")
  opt <- suppressMessages(optimal_composition(list(full = quad, sub = fit60),
                                              direction = "min", spacing = 40))
  from_cols <- paste0("from_", parts)
  expect_equal(unlist(opt$optimum[1, from_cols]), unlist(quad$ref_comp),
               ignore_attr = TRUE)
  expect_equal(unlist(opt$optimum[2, from_cols]), unlist(fit60$ref_comp),
               ignore_attr = TRUE)
  expect_output(print(opt), "each fit's reference day")
})

test_that("direction = \"min\" keeps the lowest predictions", {
  opt <- optimal_composition(quad, spacing = 30, direction = "min",
                             type = "absolute")
  est <- opt$candidates$estimate_zbmi
  expect_lt(max(est[opt$candidates$in_zone]), min(est[!opt$candidates$in_zone]))
  expect_equal(mean(opt$candidates$in_zone), 0.05, tolerance = 0.1)
})

test_that("two outcomes intersect their zones, and an empty overlap reports q_overlap", {
  opt <- optimal_composition(list(zbmi = quad, obesity = risk), direction = "min",
                             spacing = 30)
  cand <- opt$candidates
  expect_identical(cand$in_zone, cand$in_zone_zbmi & cand$in_zone_obesity)
  expect_identical(opt$zones$outcome, c("zbmi", "obesity", "overall"))
  expect_gt(sum(cand$in_zone), 0)

  expect_message(
    empty <- optimal_composition(list(zbmi = quad, obesity = risk),
                                 direction = c("max", "min"), spacing = 30),
    "share no candidate"
  )
  overall <- empty$zones[empty$zones$outcome == "overall", ]
  expect_identical(overall$n, 0L)
  hand <- max(pmin(empty$candidates$rank_zbmi, empty$candidates$rank_obesity))
  expect_equal(overall$q_overlap, hand)
  expect_lt(overall$q_overlap, 0.95)
  expect_identical(empty$optimum$centre_of, c("outcome", "outcome"))
  expect_output(print(empty), "no overlap")
})

test_that("an unnamed list is named by each fit's outcome", {
  opt <- optimal_composition(list(quad, risk), direction = "min", spacing = 40)
  expect_identical(opt$outcomes, c("zbmi", "obesity"))
  expect_error(optimal_composition(list(quad, lin), spacing = 40),
               "share the outcome name")
})

test_that("a threshold rules an outcome on the response scale", {
  opt <- optimal_composition(quad, spacing = 30, type = "absolute",
                             threshold = list(zbmi = c(-0.5, 0)))
  est <- opt$candidates$estimate_zbmi
  expect_identical(opt$candidates$in_zone_zbmi, est >= -0.5 & est <= 0)
  expect_identical(opt$zones$rule[1], "threshold")
  expect_identical(c(opt$zones$threshold_low[1], opt$zones$threshold_high[1]),
                   c(-0.5, 0))

  cut <- optimal_composition(risk, spacing = 30, direction = "min",
                             type = "absolute",
                             threshold = list(obesity = 0.15))
  prob <- stats::plogis(cut$candidates$estimate_obesity)
  expect_identical(cut$candidates$in_zone_obesity, prob <= 0.15)
  expect_error(optimal_composition(quad, threshold = list(bmi = 1)),
               "named list")
})

test_that("the box evaluates the hand-built grid", {
  comp <- as.matrix(tc_comp_df(cd)[parts])
  lo <- pmax(colMeans(comp) - 3 * apply(comp, 2, stats::sd), apply(comp, 2, min))
  hi <- pmin(colMeans(comp) + 3 * apply(comp, 2, stats::sd), apply(comp, 2, max))
  grid <- simplex_grid(60)
  inside <- rowSums(sweep(grid, 2, lo, ">=") & sweep(grid, 2, hi, "<=")) == 4
  expect_identical(count_messages(
    opt <- optimal_composition(quad, region = "box", sd = 3, spacing = 60)
  ), 0L)
  expect_identical(nrow(opt$candidates), sum(inside))
  expect_equal(sort(opt$candidates$SB), sort(grid[inside, "SB"]))
  expect_true(all(is.na(opt$candidates$mahalanobis)))
})

test_that("the ellipsoid evaluates the hand-built grid", {
  opt <- optimal_composition(quad, region = "ellipsoid", spacing = 120)
  cand <- as.matrix(opt$candidates[parts])
  expect_true(all(abs(cand / 120 - round(cand / 120)) < 1e-9))
  expect_equal(rowSums(cand), rep(1440, nrow(cand)))
  par <- opt$region$parameters$all
  expect_true(all(opt$candidates$mahalanobis <= par$r2))
  grid <- simplex_grid(120)
  z <- as.matrix(ilr_transform(as.data.frame(grid), parts)[c("ilr1", "ilr2", "ilr3")])
  inside <- stats::mahalanobis(z, par$mu, par$S) <= par$r2
  expect_identical(nrow(opt$candidates), sum(inside))
  expect_equal(sort(opt$candidates$SB), sort(grid[inside, "SB"]))
})

test_that("the box takes percentile limits through coverage", {
  comp <- as.matrix(tc_comp_df(cd)[parts])
  opt <- optimal_composition(quad, coverage = c(0.05, 0.95), spacing = 30)
  expect_equal(opt$coverage, c(0.05, 0.95))
  expect_null(opt$sd)
  expect_equal(opt$region$parameters$all$lo,
               apply(comp, 2, stats::quantile, 0.05, names = FALSE),
               ignore_attr = TRUE)
  expect_true(all(sweep(as.matrix(opt$candidates[parts]), 2,
                        opt$region$parameters$all$lo, ">=")))
  expect_output(print(opt), "5th to the 95th percentile")
  central <- optimal_composition(quad, coverage = 0.8, spacing = 30)
  expect_equal(central$coverage, c(0.1, 0.9))

  wide <- optimal_composition(quad, sd = 3, spacing = 30)
  cand <- as.matrix(wide$candidates[parts])
  expect_true(all(sweep(cand, 2, apply(comp, 2, min), ">=")))
  expect_true(all(sweep(cand, 2, apply(comp, 2, max), "<=")))
})

test_that("the region arguments are checked", {
  expect_error(optimal_composition(quad, coverage = 0.9, sd = 3), "not both")
  expect_error(optimal_composition(quad, sd = NULL), "needs .*sd.* or .*coverage")
  expect_error(optimal_composition(quad, region = "ellipsoid", sd = 3),
               "applies to .*box")
  expect_error(optimal_composition(quad, region = "ellipsoid", coverage = c(0.1, 0.9)),
               "single share")
  expect_error(optimal_composition(quad, coverage = c(0.9, 0.1)), "increasing")
  expect_error(optimal_composition(quad, sd = -1), "positive")
  expect_error(optimal_composition(quad, spacing = 400), "at most 360")
  expect_error(optimal_composition(quad, sd = 0.01, spacing = 60),
               "lies inside the region")
  # The two memory guards: the grid before the region filters it, and the
  # candidates inside.
  err <- expect_error(optimal_composition(quad, region = "ellipsoid",
                                          spacing = 1))
  expect_match(gsub("\\s+", " ", conditionMessage(err)),
               "before the region filters it, above the limit of 5,000,000.",
               fixed = TRUE)
  expect_error(optimal_composition(quad, spacing = 3),
               "above the limit of 1,000,000")
})

test_that("the observed region and a data-frame region use their own rows", {
  cd60 <- declare_composition(timeuse[1:60, ], parts, closure = 1440)
  fit60 <- comp_iv(cd60, formula = zbmi ~ ilr()^2 + I(ilr()^2), model = "lm")
  obs <- optimal_composition(fit60, region = "observed")
  expect_identical(nrow(obs$candidates), 60L)
  expect_equal(as.matrix(obs$candidates[parts]), as.matrix(tc_comp_df(cd60)[parts]),
               ignore_attr = TRUE)
  expect_true(all(is.na(obs$candidates$mahalanobis)))
  expect_null(obs$spacing)
  expect_null(obs$sd)
  expect_null(obs$coverage)
  expect_message(optimal_composition(fit60, region = "observed", coverage = 0.9),
                 "ignored")

  days <- data.frame(MVPA = c(30, 60, 90), LPA = 240, SB = c(660, 630, 600),
                     SPT = 510)
  given <- optimal_composition(fit60, region = days, zone = 0.5)
  expect_identical(nrow(given$candidates), 3L)
  expect_equal(given$candidates$MVPA, c(30, 60, 90))
  expect_error(optimal_composition(fit60, region = days[1:2, 1:3]),
               "needs a column for every part")
  expect_error(optimal_composition(fit60, region = days[0, ]), "at least one row")
  off_total <- days
  off_total$SB <- off_total$SB + 60
  expect_error(optimal_composition(fit60, region = off_total), "do not sum")
  # One per cent off the total (14.4 of 1440) is refused as well.
  one_pct <- days
  one_pct$SB <- one_pct$SB + 14.4
  expect_error(optimal_composition(fit60, region = one_pct), "do not sum")
  expect_error(optimal_composition(fit60, region = days, spacing = 10),
               "applies to the")
  expect_error(optimal_composition(fit60, region = "observed", sd = 3),
               "applies to .*box")
})

test_that("by gives one region and one zone per level", {
  mod <- comp_iv(cd, "zbmi", model = "lm", moderator = "sex")
  by_default <- suppressMessages(optimal_composition(mod, spacing = 40))
  expect_identical(by_default$by, "sex")
  expect_identical(unique(by_default$zones$sex), c("female", "male"))
  expect_error(suppressMessages(optimal_composition(list(m = mod, q = quad), spacing = 40)),
               "resolve")
  num_mod <- comp_iv(cd, "obesity", model = "glm", family = "binomial",
                     moderator = "zbmi")
  expect_error(suppressMessages(optimal_composition(num_mod, spacing = 40)),
               "covariate-conditional")
  expect_no_error(suppressMessages(
    optimal_composition(mod, spacing = 40, at = list(sex = "female"))
  ))
  opt <- suppressMessages(optimal_composition(mod, by = "sex", region = "ellipsoid",
                                              spacing = 40))
  expect_identical(opt$zones$sex, rep(c("female", "male"), each = 2))
  expect_identical(sort(unique(opt$candidates$sex)), c("female", "male"))
  z_f <- tc_pivot_coords(tc_comp_df(cd)[timeuse$sex == "female", parts])
  expect_equal(opt$region$parameters$female$mu, colMeans(z_f),
               ignore_attr = TRUE, tolerance = 1e-12)
  expect_false(isTRUE(all.equal(opt$region$parameters$female$mu,
                                opt$region$parameters$male$mu)))
  expect_output(print(opt), "sex = male")
  expect_match(printed(opt), "of each level's sample", fixed = TRUE)
})

test_that("plot() draws a by level or the at values", {
  curved <- comp_iv(cd, formula = zbmi ~ ilr()^2 + I(ilr()^2) + ilr():sex + sex,
                    model = "lm")
  opt <- suppressMessages(optimal_composition(curved, by = "sex",
                                              region = "ellipsoid", spacing = 30))
  p <- suppressMessages(quiet_draw(plot(opt, level = "male", resolution = 20)))
  expect_length(p, choose(4, 3))
  # At most one outline per panel; a thin slice through the zone has none.
  n_paths <- purrr::map_int(p, function(panel) {
    sum(purrr::map_lgl(panel$layers, function(l) inherits(l$geom, "GeomPath")))
  })
  expect_true(all(n_paths <= 1L) && any(n_paths == 1L))
  expect_null(p[[1]]$labels$caption)
  # The level drawn sets the slice and the region, so the cells differ.
  p_f <- suppressMessages(quiet_draw(plot(opt, level = "female", resolution = 20)))
  expect_false(identical(p[[1]]$data, p_f[[1]]$data))
  expect_error(plot(opt, level = "other"), "level")

  mod <- comp_iv(cd, "zbmi", model = "lm", moderator = "sex")
  pinned <- suppressMessages(optimal_composition(mod, region = "ellipsoid",
                                                 spacing = 30,
                                                 at = list(sex = "female")))
  q <- suppressMessages(quiet_draw(plot(pinned, resolution = 20)))
  expect_length(q, choose(4, 3))
  expect_null(q[[1]]$labels$caption)
  expect_error(plot(pinned, level = "male"), "with .*by")

  # A fit linear in the composition is drawn as tiles, with no smoothing
  # attempted and so no warning that a quadratic cannot be fitted.
  flat <- suppressMessages(optimal_composition(mod, by = "sex", region = "ellipsoid",
                                               spacing = 30))
  expect_no_warning(
    r <- suppressMessages(quiet_draw(plot(flat, level = "male", resolution = 15)))
  )
  expect_length(r, choose(4, 3))
  # Tiles: one hexagon per computed cell.
  expect_equal(nrow(r[[1]]$layers[[1]]$data), 6 * nrow(r[[1]]$data))
})

test_that("a factor() term over a numeric column is averaged on the tied grids", {
  t_num <- timeuse
  t_num$region <- as.integer(t_num$region)
  cd_num <- declare_composition(t_num, parts, closure = 1440)
  fit <- comp_iv(cd_num, zbmi ~ ilr() + factor(region), model = "lm")
  opt <- suppressMessages(optimal_composition(fit, spacing = 30,
                                              type = "absolute"))
  expect_per_point(opt, fit)
})

test_that("candidates emmeans cannot estimate are left out of the zone, with a warning", {
  t_gap <- timeuse
  t_gap$g2 <- factor(ifelse(t_gap$sex == "male" & t_gap$region == "C", "B",
                            as.character(t_gap$region)))
  cd_gap <- declare_composition(t_gap, parts, closure = 1440)
  fit <- comp_iv(cd_gap, formula = zbmi ~ ilr()^2 + I(ilr()^2) + sex * g2,
                 model = "lm")
  expect_warning(
    opt <- suppressMessages(optimal_composition(
      fit, spacing = 60, at = list(sex = "male", g2 = "C")
    )),
    "could not be predicted .* are left out"
  )
  expect_true(all(is.na(opt$candidates$estimate_zbmi)))
  expect_true(all(is.na(opt$candidates$rank_zbmi)))
  expect_false(any(opt$candidates$in_zone))
  expect_identical(opt$zones$n, c(0L, 0L))
})

test_that("type = \"absolute\" is refused for engines without an absolute prediction", {
  cox <- comp_iv(cd, time = "time", event = "event", model = "cox")
  expect_error(optimal_composition(cox, type = "absolute"),
               "support .type = \"relative\". only")
  skip_if_not_installed("ordinal")
  ord <- comp_iv(cd, "grade", model = "clm")
  expect_error(optimal_composition(ord, type = "absolute"),
               "support .type = \"relative\". only")
})

test_that("a linear fit informs once and reports the boundary share", {
  expect_identical(count_messages(optimal_composition(lin, spacing = 30)), 1L)
  expect_message(opt <- optimal_composition(lin, spacing = 30),
                 "linear in the composition")
  expect_message(optimal_composition(lin, spacing = 30),
                 "refit with `curvature = TRUE`")
  expect_true(opt$linear[["zbmi"]])
  expect_true(all(opt$zones$on_boundary >= 0 & opt$zones$on_boundary <= 1))
  expect_identical(count_messages(optimal_composition(quad, spacing = 30)), 0L)
  # The linear line replaces the outer-shell line, whatever the share.
  text <- printed(opt)
  expect_match(text, "has no interior optimum: the zone is the edge of the region",
               fixed = TRUE)
  expect_match(text, "Refit with `curvature = TRUE`", fixed = TRUE)
  expect_false(grepl("outer shell", text, fixed = TRUE))
  box <- optimal_composition(quad, direction = "min", coverage = c(0.05, 0.95),
                             spacing = 20)
  expect_gt(box$zones$on_boundary[1], 0.5)
  expect_match(printed(box), "% of the zone lies near the edge of the region (its outer shell)",
               fixed = TRUE)
})

test_that("a curvature = TRUE fit is searched as a curved fit", {
  curved <- comp_iv(cd, "zbmi", model = "lm", covariates = "sex",
                    curvature = TRUE)
  expect_identical(count_messages(opt <- optimal_composition(curved, spacing = 30)), 0L)
  expect_false(opt$linear[["zbmi"]])
  text <- paste(utils::capture.output(print(opt)), collapse = "\n")
  expect_false(grepl("curvature = TRUE", text, fixed = TRUE))
})

test_that("arguments are checked", {
  cd3 <- declare_composition(timeuse, c("MVPA", "SB", "SPT"), closure = 1440)
  fit3 <- comp_iv(cd3, "obesity", model = "glm", family = "binomial")
  expect_error(optimal_composition(list(zbmi = quad, obesity = fit3)),
               "share the first fit's parts")
  expect_error(optimal_composition(quad, spacing = 0), "spacing")
  expect_error(optimal_composition(quad, zone = 1), "zone")
  expect_error(optimal_composition(quad, direction = "best"), "direction")
  expect_error(optimal_composition(quad, region = "hull"), "region")
  expect_error(optimal_composition(quad, from = c(MVPA = 30)), "from")
  expect_error(optimal_composition("fit"), "timecoda_fit")
})

test_that("print() reports the zone and plot() outlines it", {
  opt <- optimal_composition(quad, region = "ellipsoid", spacing = 30,
                             direction = "min")
  expect_output(print(opt), "Goldilocks-zone search")
  expect_output(print(opt), "Gain from the reference day")
  p <- suppressMessages(quiet_draw(plot(opt, resolution = 20, smooth = FALSE)))
  expect_type(p, "list")
  expect_length(p, choose(4, 3))
  paths <- purrr::map_lgl(p, function(panel) {
    any(purrr::map_lgl(panel$layers, function(l) inherits(l$geom, "GeomPath")))
  })
  expect_true(any(paths))
  expect_true(all(purrr::map_lgl(p, function(panel) is.null(panel$labels$caption))))
  expect_error(plot(opt, outcome = "bmi"), "outcome")
})

test_that("the relative and absolute searches rank and zone alike", {
  check_types <- function(object, ...) {
    rel <- suppressMessages(optimal_composition(object, type = "relative", ...))
    abs <- suppressMessages(optimal_composition(object, type = "absolute", ...))
    flags <- grep("^(in_zone|rank_)", names(rel$candidates), value = TRUE)
    expect_identical(rel$candidates[flags], abs$candidates[flags])
    zone_cols <- c("outcome", "n", "q_overlap",
                   grep("^(centre|min|max|best)_", names(rel$zones), value = TRUE))
    expect_identical(rel$zones[zone_cols], abs$zones[zone_cols])
    expect_equal(rel$optimum$gain, abs$optimum$gain)
    expect_identical(rel$optimum[c("estimate", "conf.low", "conf.high")],
                     rlang::set_names(rel$optimum[c("gain", "gain.low", "gain.high")],
                                      c("estimate", "conf.low", "conf.high")))
    list(rel = rel, abs = abs)
  }
  both <- check_types(quad, direction = "min", spacing = 30)
  check_types(risk, direction = "min", spacing = 30)
  check_types(list(zbmi = quad, obesity = risk), direction = "min",
              spacing = 30)
  expect_identical(both$rel$type, "relative")
  expect_output(print(both$rel), "Scale +Mean difference versus the reference day")
  expect_output(print(both$abs), "Scale +predicted response")
  expect_output(print(both$abs), "Predicted zbmi at the zone centre")
  rel_text <- paste(utils::capture.output(print(both$rel)), collapse = "\n")
  expect_false(grepl("Predicted zbmi", rel_text, fixed = TRUE))
})

test_that("a relative threshold is a ratio versus from", {
  cut <- optimal_composition(risk, spacing = 30, direction = "min",
                             threshold = list(obesity = 0.8))
  expect_identical(cut$candidates$in_zone_obesity,
                   exp(cut$candidates$estimate_obesity) <= 0.8)
  expect_identical(cut$zones$threshold_high[1], 0.8)
  expect_output(print(cut), "Odds ratio versus the reference day <= 0.8")
  expect_error(optimal_composition(risk, threshold = list(obesity = c(-1, 1))),
               "must be positive")
})

# One engine searched on a coarse grid in the ellipsoid: a non-empty zone,
# finite values on the effect scale named by `measure`, and a printed report.
search_engine <- function(fit, measure, ...) {
  opt <- suppressMessages(optimal_composition(fit, region = "ellipsoid",
                                              spacing = 40, ...))
  expect_gt(opt$zones$n[opt$zones$outcome == "overall"], 0)
  expect_true(all(is.finite(opt$candidates[[paste0("estimate_", opt$outcomes)]])))
  expect_true(all(is.finite(unlist(opt$optimum[c("estimate", "gain")]))))
  expect_identical(unname(opt$measures), measure)
  expect_output(print(opt), measure)
  opt
}

test_that("a Cox fit is searched on the hazard ratio and plotted", {
  cox <- comp_iv(cd, time = "time", event = "event", model = "cox")
  opt <- search_engine(cox, "Hazard ratio", direction = "min")
  p <- suppressMessages(quiet_draw(plot(opt, resolution = 10)))
  expect_length(p, choose(4, 3))
  expect_null(p[[1]]$labels$caption)
  expect_identical(p[[1]]$scales$get_scales("fill")$name,
                   "Hazard ratio versus the reference day")
})

test_that("a mixed-effects Cox fit is searched and plotted", {
  skip_if_not_installed("coxme")
  fit <- comp_iv(cd, time = "time", event = "event", model = "coxme",
                 random = "school")
  opt <- search_engine(fit, "Hazard ratio", direction = "min")
  expect_length(suppressMessages(quiet_draw(plot(opt, resolution = 10))),
                choose(4, 3))
})

# plot() of a relative search: one panel per triplet, the surface at 20 of the
# first panel's cells equal to predict_realloc() for the same moves from the
# same day, and the outlined cells within the zone's limits.
expect_surface <- function(opt, category = NULL) {
  p <- suppressMessages(quiet_draw(plot(opt, resolution = 30)))
  expect_length(p, choose(4, 3))
  cells <- p[[1]]$data
  tri <- c(cells$xpart[1], cells$ypart[1], cells$zpart[1])
  fit <- opt$fits[[1]]
  from <- stats::setNames(unlist(opt$optimum[1, paste0("from_", parts)]), parts)
  idx <- unique(round(seq(1, nrow(cells), length.out = 20)))
  comps <- matrix(unlist(opt$optimum[1, parts]), nrow = length(idx), ncol = 4,
                  byrow = TRUE, dimnames = list(NULL, parts))
  comps[, tri] <- as.matrix(cells[idx, c("x", "y", "z")])
  pr <- predict_realloc(fit, realloc = as.data.frame(sweep(comps, 2, from)),
                        reference = from, type = "relative")
  if (!is.null(category)) pr <- pr[pr$response_level == category, ]
  expect_equal(cells$.pred[idx], pr$estimate, tolerance = 1e-8)
  zone <- opt$zones[opt$zones$outcome == opt$outcomes[1], ]
  outlined <- unlist(purrr::map(p, function(panel) {
    path <- purrr::keep(panel$layers, function(l) inherits(l$geom, "GeomPath"))
    if (length(path) > 0) path[[1]]$data$.pred
  }))
  expect_gt(length(outlined), 0)
  expect_true(all(is.na(zone$threshold_low) | outlined >= zone$threshold_low))
  expect_true(all(is.na(zone$threshold_high) | outlined <= zone$threshold_high))
  p
}

test_that("ordinal fits are searched on the cumulative odds ratio and plotted", {
  skip_if_not_installed("ordinal")
  ord <- comp_iv(cd, "grade", model = "clm", curvature = TRUE)
  opt <- search_engine(ord, "Cumulative odds ratio")
  p <- expect_surface(opt)
  expect_identical(p[[1]]$scales$get_scales("fill")$name,
                   "Cumulative odds ratio versus the reference day")
  mixed <- comp_iv(cd, "grade", model = "clmm", random = "school")
  opt_mixed <- search_engine(mixed, "Cumulative odds ratio")
  expect_length(suppressMessages(quiet_draw(plot(opt_mixed, resolution = 10))),
                choose(4, 3))
})

test_that("the ordinal and multinomial surfaces take by, at and a threshold", {
  skip_if_not_installed("ordinal")
  ord <- comp_iv(cd, "grade", model = "clm", covariates = "sex", curvature = TRUE)
  split <- suppressMessages(optimal_composition(ord, by = "sex", region = "ellipsoid",
                                                spacing = 40,
                                                threshold = list(grade = 1)))
  p <- suppressMessages(quiet_draw(plot(split, level = "male", resolution = 20)))
  expect_length(p, choose(4, 3))
  expect_null(p[[1]]$labels$caption)
  zone <- split$zones[split$zones$outcome == "grade" & split$zones$sex == "male", ]
  outlined <- unlist(purrr::map(p, function(panel) {
    path <- purrr::keep(panel$layers, function(l) inherits(l$geom, "GeomPath"))
    if (length(path) > 0) path[[1]]$data$.pred
  }))
  expect_true(all(outlined >= zone$threshold_low))
  skip_if_not_installed("mclogit")
  d <- timeuse
  d$grade_unordered <- factor(as.character(d$grade), levels = levels(d$grade))
  mb <- suppressWarnings(comp_iv(declare_composition(d, parts, closure = 1440),
                                 "grade_unordered", model = "mblogit",
                                 covariates = "sex"))
  levs <- unique(sub("~.*$", "", names(stats::coef(mb$model))))
  pinned <- suppressMessages(optimal_composition(mb, category = levs[2],
                                                 region = "ellipsoid", spacing = 40,
                                                 at = list(sex = "female")))
  q <- suppressMessages(quiet_draw(plot(pinned, resolution = 10)))
  expect_length(q, choose(4, 3))
  expect_null(q[[1]]$labels$caption)
})

test_that("direction = \"max\" favours the higher categories of an ordinal fit", {
  skip_if_not_installed("ordinal")
  ord <- comp_iv(cd, "grade", model = "clm")
  mean_category <- function(direction) {
    opt <- suppressMessages(optimal_composition(ord, direction = direction,
                                                spacing = 60))
    centre <- unlist(opt$optimum[parts])
    z <- as.data.frame(tc_pivot_coords(t(centre), pivot = ord$pivot))
    prob <- stats::predict(ord$model, newdata = z, type = "prob")$fit
    sum(prob * seq_along(prob))
  }
  expect_gt(mean_category("max"), mean_category("min"))
})

test_that("a multinomial fit is searched for one category", {
  skip_if_not_installed("mclogit")
  d <- timeuse
  d$grade_unordered <- factor(as.character(d$grade), levels = levels(d$grade))
  cd2 <- declare_composition(d, parts, closure = 1440)
  mb <- suppressWarnings(comp_iv(cd2, "grade_unordered", model = "mblogit"))
  levs <- unique(sub("~.*$", "", names(stats::coef(mb$model))))
  opt <- search_engine(mb, "Relative risk ratio", category = levs[1])
  expect_identical(opt$category, levs[1])
  curved <- suppressWarnings(comp_iv(cd2, "grade_unordered", model = "mblogit",
                                     curvature = TRUE))
  p <- expect_surface(search_engine(curved, "Relative risk ratio",
                                    category = levs[1]),
                      category = levs[1])
  expect_identical(p[[1]]$scales$get_scales("fill")$name,
                   paste("Relative risk ratio of", levs[1], "versus the reference day"))
  expect_error(optimal_composition(mb, spacing = 60), "must name the outcome level")
  expect_error(optimal_composition(mb, spacing = 60), levs[length(levs)],
               fixed = TRUE)
  expect_error(optimal_composition(mb, spacing = 60, category = "none"),
               "must be one of")
  expect_error(optimal_composition(lin, category = levs[1]),
               "multinomial fits only")
})

test_that("a factor interacting with another covariate stays on the tied grids", {
  fit <- comp_iv(cd, formula = zbmi ~ ilr()^2 + I(ilr()^2) + sex * region,
                 model = "lm")
  opt <- suppressMessages(optimal_composition(fit, spacing = 40,
                                              type = "absolute"))
  expect_per_point(opt, fit)
})

test_that("a covariate held at several values is averaged on the tied grids", {
  pins <- list(sex = c("female", "male"))
  opt <- optimal_composition(quad, spacing = 40, type = "absolute", at = pins)
  expect_per_point(opt, quad, at = pins)
  idx <- round(seq(1, nrow(opt$candidates), length.out = 20))
  nodes <- as.data.frame(opt$candidates[idx, parts])
  given <- optimal_composition(quad, region = nodes, zone = 0.5,
                               type = "absolute", at = pins)
  expect_equal(given$candidates$estimate_zbmi, opt$candidates$estimate_zbmi[idx],
               tolerance = 1e-10)
})

test_that("the ellipsoid takes another chi-square level through coverage", {
  opt <- optimal_composition(quad, region = "ellipsoid", coverage = 0.95,
                             spacing = 30)
  par <- opt$region$parameters$all
  expect_identical(par$r2, stats::qchisq(0.95, 3))
  coords <- c("ilr1", "ilr2", "ilr3")
  z <- as.matrix(ilr_transform(tc_comp_df(quad$comp), parts,
                               pivot = quad$pivot)[coords])
  zc <- as.matrix(ilr_transform(as.data.frame(opt$candidates[parts]), parts,
                                pivot = quad$pivot)[coords])
  # A small allowance for floating-point noise at the fence.
  expect_true(all(stats::mahalanobis(zc, colMeans(z), stats::cov(z)) <=
                    par$r2 + 1e-8))
  narrow <- optimal_composition(quad, region = "ellipsoid", coverage = 0.80,
                                spacing = 30)
  expect_gt(nrow(opt$candidates), nrow(narrow$candidates))
  expect_identical(opt$coverage, 0.95)
  expect_match(printed(opt), "95% chi-square", fixed = TRUE)
})

test_that("conf_level sets the gain's interval and nothing else", {
  for (fit in list(quad, risk)) {
    wide <- optimal_composition(fit, spacing = 30)
    narrow <- optimal_composition(fit, spacing = 30, conf_level = 0.90)
    expect_identical(narrow$candidates, wide$candidates)
    expect_identical(narrow$zones, wide$zones)
    expect_identical(narrow$optimum[parts], wide$optimum[parts])
    expect_identical(narrow$optimum$gain, wide$optimum$gain)
    expect_gt(narrow$optimum$gain.low, wide$optimum$gain.low)
    expect_lt(narrow$optimum$gain.high, wide$optimum$gain.high)
    expect_identical(narrow$conf_level, 0.90)
    # The same move by predict_realloc(), at the same level.
    from <- stats::setNames(unlist(narrow$optimum[1, paste0("from_", parts)]), parts)
    pr <- predict_realloc(fit, realloc = unlist(narrow$optimum[1, parts]) - from,
                          reference = from, type = "relative", conf_level = 0.90)
    expect_equal(unlist(narrow$optimum[1, c("gain", "gain.low", "gain.high")]),
                 unlist(pr[1, c("estimate", "conf.low", "conf.high")]),
                 ignore_attr = TRUE, tolerance = 1e-8)
  }
})

test_that("the ellipsoid split by a factor reports each level's share of its sample", {
  mod <- comp_iv(cd, "zbmi", model = "lm", moderator = "sex", curvature = TRUE)
  opt <- suppressMessages(optimal_composition(mod, by = "sex", region = "ellipsoid",
                                              spacing = 30))
  levs <- c("female", "male")
  coords <- c("ilr1", "ilr2", "ilr3")
  hand <- purrr::map_dbl(levs, function(lev) {
    z <- as.matrix(ilr_transform(tc_comp_df(cd)[timeuse$sex == lev, ], parts,
                                 pivot = mod$pivot)[coords])
    mean(stats::mahalanobis(z, colMeans(z), stats::cov(z)) <=
           stats::qchisq(0.80, df = 3))
  })
  share <- purrr::map_dbl(opt$region$parameters[levs], "share")
  expect_true(all(share > 0 & share < 1))
  expect_equal(share, hand, ignore_attr = TRUE)
  # The two levels' shares differ after rounding, so the report gives their range.
  pct <- round(100 * hand, 1)
  expect_false(pct[1] == pct[2])
  expect_match(printed(opt), sprintf(
    "ILR ellipsoid at the 80%% chi-square level, holding %s%% to %s%% of each level's sample",
    format(min(pct), nsmall = 1), format(max(pct), nsmall = 1)
  ), fixed = TRUE)
})

test_that("the ellipsoid label gives one share when the levels' shares round alike", {
  label <- tc_optimal_regions$ellipsoid$label
  level <- function(share) list(share = share, coverage = 0.8)
  expect_identical(
    as.character(label(list(a = level(0.80012), b = level(0.80031)), NULL)),
    "ILR ellipsoid at the 80% chi-square level, holding 80.0% of each level's sample"
  )
  expect_identical(
    as.character(label(list(a = level(0.75), b = level(0.8)), NULL)),
    "ILR ellipsoid at the 80% chi-square level, holding 75.0% to 80.0% of each level's sample"
  )
  expect_identical(
    as.character(label(list(all = level(0.8)), NULL)),
    "ILR ellipsoid at the 80% chi-square level, holding 80.0% of the sample"
  )
})

# TRUE for each panel drawn smoothed: more drawn tiles than computed cells.
smoothed_panels <- function(p) {
  purrr::map_lgl(p, function(panel) {
    nrow(panel$layers[[1]]$data) > 6 * nrow(panel$data)
  })
}

test_that("plot() of a search takes smooth = NULL, TRUE or FALSE for any fit", {
  draw <- function(fit, ...) {
    opt <- suppressMessages(optimal_composition(fit, direction = "min",
                                                spacing = 30))
    suppressMessages(quiet_draw(plot(opt, ...)))
  }
  # NULL: tiles for the linear fit, a smooth for the curved one.
  expect_false(any(smoothed_panels(draw(lin))))
  expect_true(all(smoothed_panels(draw(quad))))
  expect_true(all(smoothed_panels(draw(lin, smooth = TRUE))))
  expect_false(any(smoothed_panels(draw(quad, smooth = FALSE))))
  expect_false(any(smoothed_panels(draw(lin, smooth = FALSE))))
  expect_true(all(smoothed_panels(draw(quad, smooth = TRUE))))
  opt <- suppressMessages(optimal_composition(quad, spacing = 30))
  expect_error(plot(opt, smooth = "yes"), "`NULL`, `TRUE` or `FALSE`")
})

test_that("plot() of a search names an argument it does not take", {
  opt <- suppressMessages(optimal_composition(quad, spacing = 30))
  flat <- function(err) gsub("\\s+", " ", conditionMessage(err))
  msg <- flat(expect_error(plot(opt, coverage = c(0.01, 0.99)),
                           class = "rlang_error"))
  expect_match(
    msg,
    "`plot.optimal_comp()` was given an argument it does not take: `coverage`.",
    fixed = TRUE
  )
  expect_match(msg, "set in `optimal_composition()`", fixed = TRUE)
  msg <- flat(expect_error(plot(opt, sd = 1, quantiles = c(0.1, 0.9)),
                           class = "rlang_error"))
  expect_match(msg, "was given arguments it does not take: `sd` and `quantiles`.",
               fixed = TRUE)
})

test_that("an ellipsoid search is drawn on the cells inside the ellipsoid", {
  opt <- suppressMessages(optimal_composition(quad, region = "ellipsoid",
                                              spacing = 30, direction = "min"))
  p <- suppressMessages(quiet_draw(plot(opt, resolution = 30, smooth = FALSE)))
  par <- opt$region$parameters$all
  centre <- stats::setNames(unlist(opt$zones[opt$zones$outcome == "overall",
                                             paste0("centre_", parts)]), parts)
  box <- purrr::map(rlang::set_names(parts), function(q) range(opt$candidates[[q]]))
  fewer <- purrr::map_lgl(p, function(panel) {
    cells <- panel$data
    tri <- c(cells$xpart[1], cells$ypart[1], cells$zpart[1])
    comps <- matrix(centre, nrow = nrow(cells), ncol = 4, byrow = TRUE,
                    dimnames = list(NULL, parts))
    comps[, tri] <- as.matrix(cells[c("x", "y", "z")])
    z <- tc_pivot_coords(comps, pivot = quad$pivot)
    expect_true(all(tc_optimal_regions$ellipsoid$inside(comps, z, par)))
    # The same panel within the candidates' per-part range holds more cells.
    in_box <- tc_surface_panel(quad, tri, centre, 30, "relative", "value",
                               box, link = identity, inv = identity,
                               control = unlist(quad$ref_comp)[parts])$data
    nrow(cells) < nrow(in_box)
  })
  expect_true(any(fewer))
})

test_that("a blank panel of a search over given or observed days names the candidates", {
  cd60 <- declare_composition(timeuse[1:60, ], parts, closure = 1440)
  fit60 <- comp_iv(cd60, formula = zbmi ~ ilr()^2 + I(ilr()^2), model = "lm")
  warnings_of <- function(opt, ...) {
    w <- character()
    withCallingHandlers(
      suppressMessages(quiet_draw(plot(opt, ...))),
      warning = function(x) {
        w <<- c(w, conditionMessage(x))
        invokeRestart("muffleWarning")
      }
    )
    w
  }
  days <- data.frame(MVPA = c(20, 40, 60), LPA = 240, SB = c(680, 660, 640),
                     SPT = 500)
  given <- suppressMessages(optimal_composition(fit60, region = days,
                                                zone = 0.5))
  w <- warnings_of(given)
  expect_length(w, choose(4, 3))
  expect_match(w, "within the range of the search's candidates; the panel is blank")
  expect_match(w, "Raise `resolution`")
  expect_false(any(grepl("quantiles", w)))
  # On the observed days, a three-step lattice leaves a panel blank, and too
  # few cells for the smooth in the others.
  observed <- suppressMessages(optimal_composition(fit60, region = "observed"))
  w <- warnings_of(observed, resolution = 3, smooth = TRUE)
  blank <- grepl("the panel is blank", w)
  expect_gt(sum(blank), 0)
  expect_gt(sum(!blank), 0)
  expect_match(w[blank], "within the range of the search's candidates")
  expect_match(w[!blank], "inside the candidates' range")
  expect_false(any(grepl("quantiles", w)))
})

test_that("a search over the observed days draws the lattice within the candidates' range", {
  opt <- suppressMessages(optimal_composition(risk, direction = "min",
                                              region = "observed"))
  p <- suppressMessages(quiet_draw(plot(opt, resolution = 20, smooth = FALSE)))
  centre <- unlist(opt$zones[opt$zones$outcome == "overall",
                             paste0("centre_", parts)])
  names(centre) <- parts
  ii <- seq_len(19)
  steps <- expand.grid(i = ii, j = ii)
  steps <- steps[steps$i + steps$j < 20, ]
  cut <- purrr::map_lgl(p, function(panel) {
    cells <- panel$data
    tri <- c(cells$xpart[1], cells$ypart[1], cells$zpart[1])
    # The triplet's lattice, the other part at the zone centre, kept where
    # each part lies within the range of the search's candidates.
    total3 <- 1440 - sum(centre[setdiff(parts, tri)])
    lattice <- cbind(steps$i / 20, steps$j / 20,
                     1 - steps$i / 20 - steps$j / 20) * total3
    inside <- purrr::map(1:3, function(m) {
      r <- range(opt$candidates[[tri[m]]])
      lattice[, m] >= r[1] & lattice[, m] <= r[2]
    })
    kept <- lattice[Reduce(`&`, inside), , drop = FALSE]
    expect_equal(unname(as.matrix(cells[c("x", "y", "z")])), unname(kept),
                 tolerance = 1e-10)
    nrow(kept) < nrow(lattice)
  })
  # The range leaves part of the triangle blank.
  expect_true(any(cut))
})

test_that("the smoothed response surface of a logistic search is quadratic on the logit scale", {
  opt <- suppressMessages(optimal_composition(risk, direction = "min",
                                              region = "observed",
                                              type = "absolute"))
  p <- suppressMessages(quiet_draw(plot(opt, resolution = 20, smooth = TRUE)))
  # A full quadratic in the triplet's own two ILR coordinates.
  quadratic <- function(cells, value, at) {
    uv <- function(m) {
      d <- stats::setNames(as.data.frame(m), c("x", "y", "z"))
      stats::setNames(ilr_transform(d, c("x", "y", "z"))[c("ilr1", "ilr2")],
                      c("u", "v"))
    }
    q <- stats::lm(value ~ u + v + I(u^2) + I(v^2) + I(u * v),
                   data = cbind(uv(cells), value = value))
    stats::predict(q, uv(at))
  }
  for (panel in p) {
    cells <- panel$data
    tiles <- panel$layers[[1]]$data
    # Each tile is a hexagon around its composition: the mean of its corners.
    centres <- dplyr::summarise(dplyr::group_by(tiles, .data$cell),
                                x = mean(.data$x), y = mean(.data$y),
                                z = mean(.data$z), .pred = .data$.pred[1])
    expect_gt(nrow(centres), nrow(cells))
    xyz <- cells[c("x", "y", "z")]
    at <- centres[c("x", "y", "z")]
    expect_equal(unname(centres$.pred),
                 unname(stats::plogis(quadratic(xyz, stats::qlogis(cells$.pred),
                                                at))),
                 tolerance = 1e-8)
    # A quadratic on the probability scale draws other values.
    expect_gt(max(abs(centres$.pred - quadratic(xyz, cells$.pred, at))), 1e-3)
  }
})

test_that("print() gives the region's limits and marks the box limits cut to the data", {
  comps <- as.matrix(tc_comp_df(cd)[parts])
  n_clipped <- function(k) {
    m <- colMeans(comps)
    s <- apply(comps, 2, stats::sd)
    sum(m - k * s < apply(comps, 2, min)) + sum(m + k * s > apply(comps, 2, max))
  }
  out <- function(...) {
    opt <- suppressMessages(optimal_composition(quad, spacing = 30, ...))
    list(opt = opt, lines = utils::capture.output(print(opt)))
  }
  # The lines of the table under `title`, down to the blank line that ends it.
  table_of <- function(lines, title) {
    start <- which(lines == title)
    expect_length(start, 1)
    n <- which(lines[-seq_len(start)] == "")[1] - 1
    lines[start + seq_len(n)]
  }
  # Each part's row of a table (under its header), split into its cells.
  cells <- function(tab) {
    purrr::map(tab[1 + seq_along(parts)], function(r) strsplit(trimws(r), " +")[[1]])
  }
  marked <- function(v, cut) paste0(sprintf("%.1f", v), ifelse(cut, "*", ""))
  footnote <- function(p) {
    kinds <- c(if (any(p$lo_clipped)) "minimum", if (any(p$hi_clipped)) "maximum")
    paste("  * the part's observed", paste(kinds, collapse = " or "))
  }
  sd3 <- out()
  par <- sd3$opt$region$parameters$all
  tab <- table_of(sd3$lines, "Region limits")
  rows <- cells(tab)
  expect_identical(purrr::map_chr(rows, 1), parts)
  expect_identical(purrr::map_chr(rows, 2), marked(par$lo, par$lo_clipped))
  expect_identical(purrr::map_chr(rows, 3), marked(par$hi, par$hi_clipped))
  expect_equal(sum(par$lo_clipped) + sum(par$hi_clipped), n_clipped(3))
  expect_identical(sum(endsWith(unlist(rows), "*")), n_clipped(3))
  expect_identical(tab[length(tab)], footnote(par))
  sd5 <- out(sd = 5)
  expect_identical(n_clipped(5), 8L)
  tab <- table_of(sd5$lines, "Region limits")
  expect_identical(sum(endsWith(unlist(cells(tab)), "*")), 8L)
  expect_identical(tab[length(tab)], "  * the part's observed minimum or maximum")
  # A percentile box is never cut, and another region prints its candidates'
  # range.
  pct <- out(coverage = c(0.01, 0.99))
  expect_false(any(unlist(pct$opt$region$parameters$all[c("lo_clipped", "hi_clipped")])))
  tab <- table_of(pct$lines, "Region limits")
  expect_length(tab, 1 + length(parts))
  expect_false(any(endsWith(unlist(cells(tab)), "*")))
  ell <- out(region = "ellipsoid")
  tab <- table_of(ell$lines, "Candidate range")
  expect_length(tab, 1 + length(parts))
  rows <- cells(tab)
  cand <- ell$opt$candidates
  expect_identical(purrr::map_chr(rows, 2),
                   sprintf("%.1f", purrr::map_dbl(parts, function(q) min(cand[[q]]))))
  expect_identical(purrr::map_chr(rows, 3),
                   sprintf("%.1f", purrr::map_dbl(parts, function(q) max(cand[[q]]))))
  # With `by`, one table per level, each marking its own cut limits.
  mod <- comp_iv(cd, "zbmi", model = "lm", moderator = "sex", curvature = TRUE)
  split <- suppressMessages(optimal_composition(mod, spacing = 30))
  lines <- utils::capture.output(print(split))
  n <- purrr::map_int(c("female", "male"), function(lev) {
    p <- split$region$parameters[[lev]]
    rows <- cells(table_of(lines, paste0("Region limits, sex = ", lev)))
    expect_identical(purrr::map_chr(rows, 2), marked(p$lo, p$lo_clipped))
    expect_identical(purrr::map_chr(rows, 3), marked(p$hi, p$hi_clipped))
    sum(p$lo_clipped) + sum(p$hi_clipped)
  })
  expect_gt(sum(n), 0)
})

test_that("print() keeps every line within the console width", {
  withr::local_options(width = 80)
  labelled <- declare_composition(
    timeuse, parts, closure = 1440, units = "min/day",
    labels = c(MVPA = "Moderate-to-vigorous PA", LPA = "Light PA",
               SB = "Sedentary behaviour", SPT = "Sleep")
  )
  bmi <- comp_iv(labelled, "zbmi", model = "lm", covariates = "sex",
                 curvature = TRUE)
  ob <- comp_iv(labelled, "obesity", model = "glm", family = "binomial",
                covariates = "sex", curvature = TRUE)
  mod <- comp_iv(labelled, "zbmi", model = "lm", moderator = "sex",
                 curvature = TRUE)
  two <- suppressMessages(optimal_composition(list(zbmi = bmi, obesity = ob),
                                              direction = "min", spacing = 30))
  split <- suppressMessages(optimal_composition(mod, direction = "min",
                                                spacing = 30))
  purrr::walk(list(two, split), function(opt) {
    expect_lte(max(nchar(utils::capture.output(print(opt)))), 80)
  })
  # Rules that differ are given per outcome, each on its own line.
  opposite <- suppressMessages(optimal_composition(
    list(zbmi = bmi, obesity = ob), direction = c("max", "min"), spacing = 30
  ))
  expect_output(print(opposite), paste0(
    "Zone rule +zbmi: best 5% of predictions, higher is better\n",
    " +obesity: best 5% of predictions, lower is better"
  ))
})

test_that("the default grid equals the per-point predictions for every engine", {
  skip_on_cran()
  # The default region and spacing: the box of the mean +/- 3 SD on multiples
  # of 10 minutes.
  for (type in c("absolute", "relative")) {
    expect_per_point(suppressMessages(optimal_composition(risk, type = type)), risk)
  }
  cox <- comp_iv(cd, formula = survival::Surv(time, event) ~ ilr()^2 +
                   I(ilr()^2) + sex, model = "cox")
  expect_per_point(suppressMessages(optimal_composition(cox)), cox)
  skip_if_not_installed("lme4")
  mixed <- comp_iv(cd, "zbmi", model = "lmer", random = "school",
                   covariates = "sex", curvature = TRUE)
  for (type in c("absolute", "relative")) {
    expect_per_point(suppressMessages(optimal_composition(mixed, type = type)),
                     mixed)
  }
  skip_if_not_installed("ordinal")
  ord <- comp_iv(cd, "grade", model = "clm", curvature = TRUE)
  expect_per_point(suppressMessages(optimal_composition(ord)), ord)
  skip_if_not_installed("mclogit")
  d <- timeuse
  d$grade_unordered <- factor(as.character(d$grade), levels = levels(d$grade))
  mb <- suppressWarnings(comp_iv(declare_composition(d, parts, closure = 1440),
                                 "grade_unordered", model = "mblogit",
                                 covariates = "sex"))
  levs <- unique(sub("~.*$", "", names(stats::coef(mb$model))))
  expect_per_point(suppressMessages(optimal_composition(mb, category = levs[1])),
                   mb, category = levs[1])
})

# Outcomes whose gradient opposes z-BMI's: with lr = log(SB / MVPA), `y_opp`
# is lowest at lr = 3.75 (little MVPA, much SB). `y_sex` follows it for girls
# and follows z-BMI for boys.
opposed <- local({
  d <- timeuse
  lr <- log(d$SB / d$MVPA)
  noise <- withr::with_seed(20261118, stats::rnorm(nrow(d), sd = 0.3))
  d$y_opp <- -0.6 * lr + 0.08 * lr^2 + noise
  d$y_sex <- ifelse(d$sex == "female", d$y_opp, d$zbmi + noise)
  cdo <- declare_composition(d, parts, closure = 1440)
  list(
    opp = comp_iv(cdo, formula = y_opp ~ ilr()^2 + I(ilr()^2) + sex,
                  model = "lm"),
    zbmi_sex = comp_iv(cdo, formula = zbmi ~ (ilr()^2 + I(ilr()^2)) * sex,
                       model = "lm"),
    y_sex = comp_iv(cdo, formula = y_sex ~ (ilr()^2 + I(ilr()^2)) * sex,
                    model = "lm")
  )
})
search_opposed <- function(...) {
  optimal_composition(list(zbmi = quad, y_opp = opposed$opp),
                      direction = "min", spacing = 60, ...)
}

test_that("widen = FALSE leaves zones that share no candidate empty", {
  msgs <- character()
  kept <- withCallingHandlers(search_opposed(), message = function(m) {
    msgs <<- c(msgs, conditionMessage(m))
    invokeRestart("muffleMessage")
  })
  overall <- kept$zones[kept$zones$outcome == "overall", ]
  expect_identical(overall$n, 0L)
  expect_identical(kept$zone_used, c(all = 0.95))
  expect_false(kept$widen)
  expect_identical(kept$optimum$centre_of, c("outcome", "outcome"))
  expect_length(msgs, 1)
  expect_match(msgs, "share no candidate; they would meet at")
  expect_no_match(msgs, "Widening|widened")
  expect_identical(search_opposed(widen = FALSE)[c("candidates", "zones", "optimum")],
                   kept[c("candidates", "zones", "optimum")])
})

test_that("the level at which the zones would meet is printed rounded down", {
  msgs <- character()
  kept <- withCallingHandlers(search_opposed(), message = function(m) {
    msgs <<- c(msgs, conditionMessage(m))
    invokeRestart("muffleMessage")
  })
  q <- kept$zones$q_overlap[kept$zones$outcome == "overall"]
  msg <- gsub("\\s+", " ", msgs)
  meet <- as.numeric(sub(".*would meet at `zone = ([0-9.]+)`.*", "\\1", msg))
  expect_identical(meet, floor(q * 1000) / 1000)
  plain <- paste(utils::capture.output(print(kept)), collapse = "\n")
  expect_match(plain,
               paste0("the zones meet at zone = ", format(meet), "(\n|$)"))
  # A search at the printed level overlaps; one a thousandth above does not.
  overall_n <- function(zone) {
    o <- suppressMessages(search_opposed(zone = zone))
    o$zones$n[o$zones$outcome == "overall"]
  }
  expect_gt(overall_n(meet), 0L)
  expect_identical(overall_n(meet + 0.001), 0L)
})

test_that("widen = TRUE takes the first 5-point step at which the zones overlap", {
  kept <- suppressMessages(search_opposed())
  q <- kept$zones$q_overlap[kept$zones$outcome == "overall"]
  steps <- c(0.9, 0.85, 0.8, 0.75, 0.7, 0.65, 0.6, 0.55, 0.5)
  msg <- gsub("\\s+", " ", conditionMessage(
    expect_message(widened <- search_opposed(widen = TRUE))
  ))
  expect_match(msg, "share no candidate at `zone = 0.95`; they were widened to `zone = 0.6`, the best 40% of each outcome's predictions.",
               fixed = TRUE)
  level <- widened$zone_used[["all"]]
  expect_identical(level, 0.6)
  expect_identical(level, steps[steps <= q][1])
  expect_true(widened$widen)
  expect_identical(widened$zone, 0.95)
  expect_gt(widened$zones$n[widened$zones$outcome == "overall"], 0L)
  expect_identical(widened$optimum$centre_of, c("overall", "overall"))
  expect_identical(widened$zones$q_overlap, kept$zones$q_overlap)
})

test_that("a widened search equals a search at the level reached", {
  widened <- suppressMessages(search_opposed(widen = TRUE))
  direct <- suppressMessages(search_opposed(zone = widened$zone_used[["all"]]))
  expect_identical(widened$candidates, direct$candidates)
  expect_identical(widened$zones, direct$zones)
  expect_identical(widened$optimum, direct$optimum)
  # plot() reads the widened zones as it reads a direct search's.
  p_w <- suppressMessages(quiet_draw(plot(widened, resolution = 10, smooth = FALSE)))
  p_d <- suppressMessages(quiet_draw(plot(direct, resolution = 10, smooth = FALSE)))
  expect_identical(purrr::map(p_w, function(p) purrr::map(p$layers, "data")),
                   purrr::map(p_d, function(p) purrr::map(p$layers, "data")))
})

test_that("widening never changes an outcome ruled by a threshold", {
  args <- list(object = list(zbmi = quad, y_opp = opposed$opp, obesity = risk),
               direction = "min", spacing = 60,
               threshold = list(obesity = 1.5))
  kept <- suppressMessages(do.call(optimal_composition, args))
  widened <- suppressMessages(do.call(optimal_composition, c(args, widen = TRUE)))
  expect_lt(widened$zone_used[["all"]], 0.95)
  expect_identical(widened$candidates$in_zone_obesity,
                   kept$candidates$in_zone_obesity)
  expect_identical(widened$zones[widened$zones$outcome == "obesity", ],
                   kept$zones[kept$zones$outcome == "obesity", ])
  expect_false(identical(widened$candidates$in_zone_zbmi,
                         kept$candidates$in_zone_zbmi))
})

test_that("with by, each level is widened on its own", {
  widened <- suppressMessages(optimal_composition(
    list(zbmi = opposed$zbmi_sex, y = opposed$y_sex), direction = "min",
    spacing = 60, by = "sex", widen = TRUE
  ))
  expect_identical(names(widened$zone_used), c("female", "male"))
  expect_lt(widened$zone_used[["female"]], 0.95)
  expect_identical(widened$zone_used[["male"]], 0.95)
  overall <- widened$zones[widened$zones$outcome == "overall", ]
  expect_true(all(overall$n > 0))
  text <- paste(utils::capture.output(print(widened)), collapse = "\n")
  expect_match(text, "\nsex = female: [0-9,]+ candidates\n  No candidate is in every zone")
  expect_match(text, "\nsex = male: [0-9,]+ candidates\n\n")
})

test_that("zones that do not meet by the best 50% keep the result at zone", {
  three <- list(low = quad, high = quad, cut = lin)
  dirs <- c("min", "max", "min")
  first <- suppressMessages(optimal_composition(three, direction = dirs,
                                                spacing = 60))
  cut <- unname(stats::quantile(first$candidates$estimate_cut, 0.1))
  kept <- suppressMessages(optimal_composition(
    three, direction = dirs, spacing = 60, threshold = list(cut = cut)
  ))
  msgs <- character()
  widened <- withCallingHandlers(
    optimal_composition(three, direction = dirs, spacing = 60,
                        threshold = list(cut = cut), widen = TRUE),
    message = function(m) {
      msgs <<- c(msgs, conditionMessage(m))
      invokeRestart("muffleMessage")
    }
  )
  expect_lt(widened$zones$q_overlap[widened$zones$outcome == "overall"], 0.5)
  expect_identical(widened$zone_used, c(all = 0.95))
  expect_identical(widened[c("candidates", "zones", "optimum")],
                   kept[c("candidates", "zones", "optimum")])
  msg <- gsub("\\s+", " ", msgs[grepl("share no candidate", msgs)])
  expect_length(msg, 1)
  expect_match(msg, "Widening in 5-point steps to `zone = 0.5`, the best 50%",
               fixed = TRUE)
  expect_identical(utils::capture.output(print(widened)),
                   utils::capture.output(print(kept)))
  # From zone = 0.5 there is no step to take, and nothing is said of one.
  msgs <- character()
  low <- withCallingHandlers(
    optimal_composition(three, direction = dirs, spacing = 60, zone = 0.5,
                        threshold = list(cut = cut), widen = TRUE),
    message = function(m) {
      msgs <<- c(msgs, conditionMessage(m))
      invokeRestart("muffleMessage")
    }
  )
  expect_identical(low$zone_used, c(all = 0.5))
  expect_identical(low$zones$n[low$zones$outcome == "overall"], 0L)
  expect_no_match(msgs, "Widening|widened")
})

test_that("widen changes nothing when there is nothing to widen", {
  same <- function(a, b) {
    expect_identical(a[c("candidates", "zones", "optimum", "zone_used")],
                     b[c("candidates", "zones", "optimum", "zone_used")])
  }
  expect_identical(count_messages(one <- optimal_composition(quad, spacing = 60,
                                                             widen = TRUE)), 0L)
  same(one, optimal_composition(quad, spacing = 60))
  meet <- list(zbmi = quad, obesity = risk)
  same(suppressMessages(optimal_composition(meet, direction = "min",
                                            spacing = 60, widen = TRUE)),
       suppressMessages(optimal_composition(meet, direction = "min",
                                            spacing = 60)))
  low <- suppressMessages(search_opposed(zone = 0.5, widen = TRUE))
  same(low, suppressMessages(search_opposed(zone = 0.5)))
  expect_identical(low$zone_used, c(all = 0.5))
  expect_error(search_opposed(widen = NA), "must be `TRUE` or `FALSE`")
  expect_error(search_opposed(widen = "yes"), "must be `TRUE` or `FALSE`")
  expect_error(search_opposed(widen = c(TRUE, FALSE)), "must be `TRUE` or `FALSE`")
})

test_that("print() states the share a widened search reached", {
  kept <- suppressMessages(search_opposed())
  widened <- suppressMessages(search_opposed(widen = TRUE))
  text <- paste(utils::capture.output(print(widened)), collapse = "\n")
  expect_match(text, "\n  Zone rule +best 5% of predictions, lower is better\n")
  expect_match(printed(widened), "No candidate is in every zone at zone = 0.95; the zones are widened to zone = 0.6, the best 40% of each outcome's predictions",
               fixed = TRUE)
  expect_match(text, "\nOverall zone, widened to the best 40%: ", fixed = TRUE)
  expect_match(text, "Gain from the reference day to the zone centre", fixed = TRUE)
  plain <- paste(utils::capture.output(print(kept)), collapse = "\n")
  expect_no_match(plain, "widened")
  expect_match(plain, "\nOverall zone: no overlap; the zones meet at zone = ",
               fixed = TRUE)
})
