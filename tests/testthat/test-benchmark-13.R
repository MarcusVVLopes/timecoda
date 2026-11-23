# Numerical benchmark of the total-as-covariate model,
# comp_iv(include_total = TRUE): the covariate is
# t_total = sqrt(D) * rowMeans(log(x)) on the zero-free original-scale parts
# (Coenders, Martín-Fernández and Ferrer-Rosell 2017, eq. 3.2) and completes
# the pivot basis orthonormally; every fit equals the model fitted by hand on
# robCompositions::pivotCoord() coordinates plus that column (stats::lm(),
# stats::glm(), survival::coxph(), lme4::lmer()); the total test is the one-df
# Wald chi-square, invariant to the covariate's scale, to a common rescaling
# of the parts and to the part order; the omnibus leaves the total out;
# reallocations and predictions hold the total at total_ref; and a constant
# total is refused. Reference calls are namespace-qualified so each number is
# visibly computed outside timecoda.

parts <- c("MVPA", "LPA", "SB", "SPT")
D <- length(parts)
coord_names <- function(p) paste0("ilr", seq_len(length(p) - 1))
coords <- coord_names(parts)

# A variable-total copy of timeuse: every row's parts scaled by one factor.
set.seed(20260923)
tu <- timeuse
tu$count <- as.integer(round(exp(pmin(pmax(tu$zbmi, -3), 3))))
tu$grade_unordered <- factor(as.character(tu$grade), levels = levels(tu$grade))
invisible(stats::runif(nrow(tu)))
invisible(stats::rnorm(nrow(tu)))
mult <- stats::runif(nrow(tu), 0.8, 1.2)
tu_total <- tu
tu_total[parts] <- tu_total[parts] * mult
cd_total <- declare_composition(tu_total, parts, closure = "mean")

zero_rows <- c(3, 17, 42, 101, 250)
tu_zero <- tu_total
tu_zero$MVPA[zero_rows] <- 0
cd_zero <- declare_composition(tu_zero, parts, closure = "mean")

set.seed(5)
n5 <- 300
parts5 <- paste0("P", 1:5)
raw5 <- cbind(
  P1 = stats::rlnorm(n5, log(60), 0.4), P2 = stats::rlnorm(n5, log(200), 0.4),
  P3 = stats::rlnorm(n5, log(400), 0.4), P4 = stats::rlnorm(n5, log(300), 0.4),
  P5 = stats::rlnorm(n5, log(480), 0.3)
) * stats::rlnorm(n5, 0, 0.3)
df5 <- as.data.frame(raw5)
df5$y <- 0.5 * as.numeric(scale(log(raw5[, "P1"] / raw5[, "P3"]))) +
  0.8 * as.numeric(scale(log(rowSums(raw5)))) + stats::rnorm(n5)
cd5 <- declare_composition(df5, parts5, closure = "mean")

# Max abs difference; Inf when the lengths or the NA / Inf patterns differ.
max_diff <- function(a, b) {
  a <- as.numeric(unlist(a))
  b <- as.numeric(unlist(b))
  if (length(a) != length(b) || !identical(is.na(a), is.na(b))) return(Inf)
  ok <- !is.na(a)
  if (!identical(a[ok & !is.finite(a)], b[ok & !is.finite(a)])) return(Inf)
  ok <- ok & is.finite(a)
  max(c(0, abs(a[ok] - b[ok])))
}
coefs <- function(m) if (inherits(m, "merMod")) lme4::fixef(m) else stats::coef(m)
vc <- function(m, nm = names(coefs(m))) as.matrix(stats::vcov(m))[nm, nm]
wald <- function(b, V, nm) {
  as.numeric(t(b[nm]) %*% solve(V[nm, nm, drop = FALSE]) %*% b[nm])
}
eq32 <- function(cd) {
  raw <- as.matrix(cd$parts_imputed[cd$parts])
  sqrt(ncol(raw)) * rowMeans(log(raw))
}
# The data columns, the pivot coordinates of the closed composition, and the
# hand covariate.
hand_frame <- function(cd) {
  p <- cd$parts
  z <- robCompositions::pivotCoord(as.data.frame(unclass(cd$composition)),
                                   pivotvar = 1)
  colnames(z) <- coord_names(p)
  out <- cbind(as.data.frame(cd$data_raw[setdiff(names(cd$data_raw), p)]), z)
  out$t_total <- eq32(cd)
  out
}
hf_total <- hand_frame(cd_total)
ilr_rhs <- paste(coords, collapse = " + ")
hand_formula <- function(lhs, extra = "") {
  stats::as.formula(paste(lhs, "~", ilr_rhs, extra, "+ t_total"))
}

fit_cache <- new.env()
fit_of <- function(name) {
  if (!is.null(fit_cache[[name]])) return(fit_cache[[name]])
  fit_cache[[name]] <- switch(name,
    l = comp_iv(cd_total, "zbmi", model = "lm", include_total = TRUE),
    l_cov = comp_iv(cd_total, "zbmi", model = "lm", covariates = "sex",
                    include_total = TRUE),
    g_bin = comp_iv(cd_total, "obesity", model = "glm", family = "binomial",
                    include_total = TRUE),
    g_pois = comp_iv(cd_total, "count", model = "glm", family = "poisson",
                     include_total = TRUE),
    cox = comp_iv(cd_total, model = "cox", time = "time", event = "event",
                  include_total = TRUE),
    lmer = comp_iv(cd_total, "zbmi", model = "lmer", random = "school",
                   include_total = TRUE),
    f_mod = comp_iv(cd_total, zbmi ~ ilr() * sex, include_total = TRUE),
    f_mod_arg = comp_iv(cd_total, "zbmi", model = "lm", moderator = "sex",
                        include_total = TRUE),
    f_sq = comp_iv(cd_total, zbmi ~ ilr() + I(ilr()^2), include_total = TRUE),
    l_zero = comp_iv(cd_zero, "zbmi", model = "lm", include_total = TRUE),
    d5 = comp_iv(cd5, "y", model = "lm", include_total = TRUE)
  )
}

test_that("the total covariate is eq. 3.2 on the zero-free original-scale parts", {
  skip_on_cran()
  skip_if_not_installed("survival")
  skip_if_not_installed("lme4")
  for (nm in c("l", "l_cov", "g_bin", "g_pois", "cox", "lmer", "f_mod", "f_sq",
               "l_zero", "d5")) {
    f <- fit_of(nm)
    raw <- as.matrix(f$comp$parts_imputed[f$parts])
    k <- ncol(raw)
    expect_lt(max_diff(f$data$t_total, sqrt(k) * rowMeans(log(raw))), 1e-10)
    expect_lt(max_diff(f$data$t_total, sqrt(k) * log(apply(raw, 1, prod)^(1 / k))),
              1e-10)
  }
  # The zero rows: finite, and not the same formula on the raw data (log 0).
  t_pkg <- fit_of("l_zero")$data$t_total
  t_raw <- sqrt(D) * rowMeans(log(as.matrix(cd_zero$data_raw[parts])))
  expect_true(all(is.finite(t_pkg[zero_rows])))
  expect_true(all(abs(t_pkg[zero_rows] - t_raw[zero_rows]) > 0))
  # No hardcoded D: the factor on five parts is sqrt(5).
  expect_lt(max_diff(fit_of("d5")$data$t_total /
                       rowMeans(log(as.matrix(cd5$parts_imputed[parts5]))),
                     rep(sqrt(5), n5)), 1e-10)
})

test_that("the total covariate completes the pivot basis orthonormally", {
  skip_on_cran()
  skip_if_not_installed("robCompositions")
  for (nm in c("l", "l_zero", "d5")) {
    f <- fit_of(nm)
    cd <- f$comp
    raw <- as.matrix(cd$parts_imputed[cd$parts])
    z_raw <- as.matrix(robCompositions::pivotCoord(as.data.frame(raw), pivotvar = 1))
    z_cl <- as.matrix(robCompositions::pivotCoord(
      as.data.frame(unclass(cd$composition)), pivotvar = 1))
    t <- f$data$t_total
    expect_lt(max_diff(rowSums(z_raw^2) + t^2, rowSums(log(raw)^2)), 1e-8)
    expect_lt(max_diff(rowSums(z_cl^2) + t^2, rowSums(log(raw)^2)), 1e-8)
    expect_lt(max_diff(z_cl, z_raw), 1e-8)
    expect_lt(max_diff(as.matrix(f$data[coord_names(cd$parts)]), z_cl), 1e-12)
  }
})

test_that("lm and glm fits equal hand refits and the total test is the 1-df Wald", {
  skip_on_cran()
  skip_if_not_installed("robCompositions")
  hand <- list(
    l = stats::lm(hand_formula("zbmi"), data = hf_total),
    l_cov = stats::lm(hand_formula("zbmi", "+ sex"), data = hf_total),
    g_bin = stats::glm(hand_formula("obesity"), family = stats::binomial(),
                       data = hf_total),
    g_pois = stats::glm(hand_formula("count"), family = stats::poisson(),
                        data = hf_total),
    l_zero = stats::lm(hand_formula("zbmi"), data = hand_frame(cd_zero)),
    d5 = stats::lm(stats::as.formula(paste(
      "y ~", paste(coord_names(parts5), collapse = " + "), "+ t_total")),
      data = hand_frame(cd5))
  )
  for (nm in names(hand)) {
    f <- fit_of(nm)
    h <- hand[[nm]]
    b <- coefs(f$model)
    expect_setequal(names(b), names(coefs(h)))
    expect_lt(max_diff(b, coefs(h)[names(b)]), 1e-8)
    expect_lt(max_diff(vc(f$model), vc(h, names(b))), 1e-8)
    w <- coefs(h)[["t_total"]]^2 / vc(h)["t_total", "t_total"]
    expect_lt(max_diff(f$total_test$statistic, w), 1e-10)
    expect_lt(max_diff(f$total_test$p.value,
                       stats::pchisq(w, 1, lower.tail = FALSE)), 1e-10)
  }
  # summary.lm() reports a t-test of the same coefficient: t^2 is the Wald.
  st <- summary(hand$l)$coefficients["t_total", ]
  expect_lt(max_diff(st[["t value"]]^2, fit_of("l")$total_test$statistic), 1e-8)
})

test_that("Cox and lmer fits equal hand refits on the same columns", {
  skip_on_cran()
  skip_if_not_installed("survival")
  skip_if_not_installed("lme4")
  skip_if_not_installed("robCompositions")
  h_cox <- survival::coxph(hand_formula("survival::Surv(time, event)"),
                           data = hf_total)
  f <- fit_of("cox")
  b <- coefs(f$model)
  expect_lt(max_diff(b, coefs(h_cox)[names(b)]), 1e-6)
  expect_lt(max_diff(vc(f$model), vc(h_cox, names(b))), 1e-6)
  w <- coefs(h_cox)[["t_total"]]^2 / vc(h_cox)["t_total", "t_total"]
  expect_lt(max_diff(f$total_test$statistic, w), 1e-8)

  f <- fit_of("lmer")
  reml <- f$model@call$REML
  reml <- if (is.null(reml)) TRUE else eval(reml)
  h_lmer <- lme4::lmer(hand_formula("zbmi", "+ (1 | school)"), data = hf_total,
                       REML = reml)
  b <- coefs(f$model)
  expect_lt(max_diff(b, coefs(h_lmer)[names(b)]), 1e-6)
  expect_lt(max_diff(vc(f$model), vc(h_lmer, names(b))), 1e-6)
})

test_that("the total test is invariant to its scale, a common rescaling and the part order", {
  skip_on_cran()
  l <- fit_of("l")
  b1 <- stats::coef(l$model)
  # The sqrt(D - 1) convention as a plain covariate.
  tu_alt <- tu_total
  tu_alt$t_alt <- sqrt(D - 1) * rowMeans(log(as.matrix(cd_total$parts_imputed[parts])))
  f_alt <- comp_iv(declare_composition(tu_alt, parts, closure = "mean"), "zbmi",
                   model = "lm", covariates = "t_alt")
  b2 <- stats::coef(f_alt$model)
  w_alt <- b2[["t_alt"]]^2 / stats::vcov(f_alt$model)["t_alt", "t_alt"]
  expect_lt(max_diff(b1[["t_total"]] * sqrt(D), b2[["t_alt"]] * sqrt(D - 1)), 1e-8)
  expect_lt(max_diff(b1[coords], b2[coords]), 1e-8)
  expect_lt(max_diff(w_alt, l$total_test$statistic), 1e-10)
  expect_lt(max_diff(stats::pchisq(w_alt, 1, lower.tail = FALSE),
                     l$total_test$p.value), 1e-10)
  # Every part doubled: t shifts by sqrt(D) log 2 and only the intercept moves.
  tu_c2 <- tu_total
  tu_c2[parts] <- tu_total[parts] * 2
  f_c2 <- comp_iv(declare_composition(tu_c2, parts, closure = "mean"), "zbmi",
                  model = "lm", include_total = TRUE)
  b_c2 <- stats::coef(f_c2$model)
  expect_lt(max_diff(f_c2$data$t_total - l$data$t_total,
                     rep(sqrt(D) * log(2), nrow(tu_total))), 1e-8)
  expect_lt(max_diff(b_c2[c(coords, "t_total")], b1[c(coords, "t_total")]), 1e-8)
  expect_lt(max_diff(b_c2[["(Intercept)"]] - b1[["(Intercept)"]],
                     -b1[["t_total"]] * sqrt(D) * log(2)), 1e-8)
  expect_lt(max_diff(f_c2$total_test[c("statistic", "df", "p.value")],
                     l$total_test[c("statistic", "df", "p.value")]), 1e-8)
  # Part order: the basis changes, the fit and the total test do not.
  for (ord in list(c("SB", "SPT", "MVPA", "LPA"), c("SPT", "LPA", "MVPA", "SB"))) {
    f_o <- comp_iv(declare_composition(tu_total, ord, closure = "mean"), "zbmi",
                   model = "lm", include_total = TRUE)
    expect_lt(max_diff(stats::fitted(f_o$model), stats::fitted(l$model)), 1e-8)
    expect_lt(max_diff(f_o$total_test[c("statistic", "df", "p.value")],
                       l$total_test[c("statistic", "df", "p.value")]), 1e-8)
    expect_lt(max_diff(stats::coef(f_o$model)[["t_total"]], b1[["t_total"]]), 1e-8)
  }
})

test_that("the omnibus test leaves the total out", {
  skip_on_cran()
  skip_if_not_installed("survival")
  skip_if_not_installed("lme4")
  coord_only <- function(nm) {
    nm[vapply(strsplit(nm, ":", fixed = TRUE),
              function(s) all(grepl("ilr[0-9]", s)), logical(1))]
  }
  for (nm in c("l", "l_cov", "g_bin", "cox", "lmer", "f_sq", "f_mod", "f_mod_arg")) {
    f <- fit_of(nm)
    b <- coefs(f$model)
    only <- coord_only(names(b))
    expect_false("t_total" %in% only)
    expect_identical(as.integer(f$omnibus$df),
                     if (nm == "f_sq") 2L * (D - 1L) else D - 1L)
    expect_identical(length(only), as.integer(f$omnibus$df))
    expect_lt(max_diff(f$omnibus$statistic, wald(b, vc(f$model), only)), 1e-8)
  }
  # The moderator argument stores the moderation test; the formula spelling of
  # the same model does not. t_total enters neither interaction.
  f_mod <- fit_of("f_mod")
  f_arg <- fit_of("f_mod_arg")
  expect_null(f_mod$moderation)
  expect_false(is.null(f_arg$moderation))
  b_fm <- stats::coef(f_mod$model)
  b_fa <- stats::coef(f_arg$model)
  expect_setequal(names(b_fa), names(b_fm))
  expect_lt(max_diff(b_fa, b_fm[names(b_fa)]), 1e-8)
  expect_identical(sum(grepl("t_total", names(b_fm), fixed = TRUE)), 1L)
  expect_identical(sum(grepl("t_total", names(b_fa), fixed = TRUE)), 1L)
  int_nm <- grep(":", names(b_fa), fixed = TRUE, value = TRUE)
  expect_false(any(grepl("t_total", int_nm)))
  expect_lt(max_diff(f_arg$moderation$statistic, wald(b_fa, vc(f_arg$model), int_nm)),
            1e-8)
})

test_that("reallocations and predictions hold the total at total_ref", {
  skip_on_cran()
  skip_if_not_installed("lme4")
  skip_if_not_installed("robCompositions")
  for (nm in c("l", "g_bin", "lmer", "l_zero", "d5")) {
    f <- fit_of(nm)
    raw <- as.matrix(f$comp$parts_imputed[f$parts])
    expect_lt(max_diff(f$total_ref, mean(f$data$t_total)), 1e-10)
    expect_lt(max_diff(f$total_ref,
                       sqrt(ncol(raw)) * mean(log(exp(colMeans(log(raw)))))), 1e-10)
  }
  # The reference composition and 15 minutes moved from SB to MVPA.
  comp_closed <- as.matrix(unclass(cd_total$composition))
  g <- exp(colMeans(log(comp_closed)))
  ref <- g / sum(g) * fit_of("l")$total
  new <- ref
  new[["SB"]] <- new[["SB"]] - 15
  new[["MVPA"]] <- new[["MVPA"]] + 15
  zz <- as.matrix(robCompositions::pivotCoord(as.data.frame(rbind(ref, new)),
                                              pivotvar = 1))
  colnames(zz) <- coords
  expect_lt(max_diff(unlist(fit_of("l")$ref_comp)[parts], ref[parts]), 1e-9)
  design <- function(f) {
    X <- cbind(1, zz, f$total_ref)
    colnames(X) <- c("(Intercept)", coords, "t_total")
    X
  }
  for (nm in c("l", "g_bin")) {
    f <- fit_of(nm)
    X <- design(f)
    b <- coefs(f$model)[colnames(X)]
    eta <- as.numeric(X %*% b)
    d <- X[2, ] - X[1, ]
    se_link <- sqrt(as.numeric(t(d) %*% vc(f$model, colnames(X)) %*% d))
    est <- if (nm == "g_bin") exp(eta[2] - eta[1]) else eta[2] - eta[1]
    se <- if (nm == "g_bin") est * se_link else se_link
    inv <- if (nm == "g_bin") stats::plogis else identity
    pr <- predict_realloc(f, delta = 15, from = "SB", to = "MVPA")
    expect_lt(max_diff(pr$estimate, est), 1e-8)
    expect_lt(max_diff(pr$std.error, se), 1e-8)
    pa <- predict_realloc(f, delta = 15, from = "SB", to = "MVPA",
                          type = "absolute", reference = f$ref_comp)
    expect_lt(max_diff(pa$estimate, inv(eta)), 1e-8)
  }
  f <- fit_of("lmer")
  X <- design(f)
  pa <- predict_realloc(f, delta = 15, from = "SB", to = "MVPA", type = "absolute")
  expect_lt(max_diff(pa$estimate, as.numeric(X %*% coefs(f$model)[colnames(X)])), 1e-6)
})

test_that("every single-predictor fit carries a one-df total test", {
  skip_on_cran()
  skip_if_not_installed("survival")
  skip_if_not_installed("lme4")
  for (nm in c("l", "l_cov", "g_bin", "g_pois", "cox", "lmer", "f_mod",
               "f_mod_arg", "f_sq", "l_zero", "d5")) {
    expect_identical(as.integer(fit_of(nm)$total_test$df), 1L)
  }
})

test_that("a per-category fit tests one total coefficient per non-baseline level", {
  skip_on_cran()
  skip_if_not_installed("mclogit")
  f <- comp_iv(cd_total, "grade_unordered", model = "mblogit", include_total = TRUE)
  k <- nlevels(tu_total$grade_unordered)
  b <- stats::coef(f$model)
  V <- as.matrix(stats::vcov(f$model))
  idx <- which(grepl("~t_total$", names(b)))
  expect_length(idx, k - 1)
  expect_identical(as.integer(f$total_test$df), k - 1L)
  expect_lt(max_diff(f$total_test$statistic,
                     as.numeric(t(b[idx]) %*% solve(V[idx, idx]) %*% b[idx])), 1e-8)
})

test_that("a constant or near-constant total is refused", {
  skip_on_cran()
  near <- function(eps) {
    x <- timeuse
    x[parts] <- x[parts] * (1 + eps * seq(-1, 1, length.out = nrow(x)))
    declare_composition(x, parts, closure = 1440)
  }
  fit_bin <- function(cd) {
    comp_iv(cd, "obesity", model = "glm", family = "binomial", include_total = TRUE)
  }
  expect_error(fit_bin(declare_composition(timeuse, parts, closure = 1440)), "variable")
  expect_error(fit_bin(near(1e-9)), "variable")
  expect_s3_class(fit_bin(near(1e-6)), "comp_iv")
})
