# Numerical checks of comp_iv() for model = "cox", "coxme", "lmer" and
# "glmer", of summary()'s per-part effects, predict_realloc() and fit_check()
# on those fits, against survival::coxph(), coxme::coxme(), lme4::lmer() and
# lme4::glmer() fitted by hand on robCompositions::pivotCoord() coordinates;
# the reallocation contrast exp(d'b) and its standard error computed by hand;
# the Kenward-Roger degrees of freedom from pbkrtest and the Satterthwaite
# ones from lmerTest; survival::cox.zph(),
# DHARMa and car::linearHypothesis() on the hand fits. Reference calls are
# namespace-qualified so each number is visibly computed outside timecoda.

skip_on_cran()
skip_if_not_installed("robCompositions")
skip_if_not_installed("lme4")

closure <- 1440
parts <- c("MVPA", "LPA", "SB", "SPT")
parts3 <- c("PA", "SB", "SPT")
tu <- as.data.frame(timeuse)
tu$count <- as.integer(round(exp(pmin(pmax(tu$zbmi, -3), 3))))
tu$PA <- tu$MVPA + tu$LPA
deltas <- c(10, 30)

coord_names <- function(p) paste0("ilr", seq_len(length(p) - 1))
max_diff <- function(a, b) {
  a <- as.numeric(unlist(a))
  b <- as.numeric(unlist(b))
  if (length(a) != length(b)) return(Inf)
  max(abs(a - b))
}
expect_close <- function(a, b, tol) expect_lte(max_diff(a, b), tol)

hand_frame <- function(k = 1, p = parts) {
  z <- robCompositions::pivotCoord(as.data.frame(tu[p]), pivotvar = k)
  colnames(z) <- coord_names(p)
  cbind(tu, z)
}
pivot_row <- function(x, p) {
  z <- robCompositions::pivotCoord(as.data.frame(t(x[p])))
  stats::setNames(as.numeric(as.matrix(z)), coord_names(p))
}
hand_ref <- function(p) {
  x <- as.matrix(tu[p])
  x <- x / rowSums(x)
  g <- exp(colMeans(log(x)))
  g / sum(g) * closure
}

specs <- list(
  cox      = list(engine = "cox", cov = NULL, p = parts),
  cox_sex  = list(engine = "cox", cov = "sex", p = parts),
  cox_d3   = list(engine = "cox", cov = NULL, p = parts3),
  lmer     = list(engine = "lmer", cov = NULL, p = parts),
  lmer_sex = list(engine = "lmer", cov = "sex", p = parts),
  lmer_d3  = list(engine = "lmer", cov = NULL, p = parts3),
  glmer_b  = list(engine = "glmer_b", cov = NULL, p = parts),
  glmer_p  = list(engine = "glmer_p", cov = NULL, p = parts),
  coxme    = list(engine = "coxme", cov = NULL, p = parts)
)
available <- function(nms) {
  if (!requireNamespace("coxme", quietly = TRUE)) nms <- setdiff(nms, "coxme")
  nms
}

pkg_fit <- function(s, reverse = FALSE) {
  p <- if (reverse) rev(s$p) else s$p
  comp <- declare_composition(tu, p, closure = closure)
  switch(s$engine,
    cox = comp_iv(comp, time = "time", event = "event", model = "cox",
                  covariates = s$cov),
    coxme = comp_iv(comp, time = "time", event = "event", model = "coxme",
                    random = "school", covariates = s$cov),
    lmer = comp_iv(comp, "zbmi", model = "lmer", random = "school",
                   covariates = s$cov),
    glmer_b = comp_iv(comp, "obesity", model = "glmer", family = "binomial",
                      random = "school", covariates = s$cov),
    glmer_p = comp_iv(comp, "count", model = "glmer", family = "poisson",
                      random = "school", covariates = s$cov)
  )
}
hand_fit <- function(s, k = 1, control = NULL, model_frame = TRUE,
                     data = hand_frame(k, s$p)) {
  h <- data
  rhs <- paste(c(coord_names(s$p), s$cov), collapse = " + ")
  f <- function(lhs, bars = TRUE) {
    stats::as.formula(paste(lhs, "~", rhs, if (bars) "+ (1 | school)"))
  }
  switch(s$engine,
    cox = survival::coxph(f("survival::Surv(time, event)", FALSE), data = h,
                          model = model_frame),
    coxme = coxme::coxme(f("survival::Surv(time, event)"), data = h),
    lmer = lme4::lmer(f("zbmi"), data = h),
    glmer_b = if (is.null(control)) {
      lme4::glmer(f("obesity"), data = h, family = stats::binomial())
    } else {
      lme4::glmer(f("obesity"), data = h, family = stats::binomial(),
                  control = control)
    },
    glmer_p = lme4::glmer(f("count"), data = h, family = stats::poisson())
  )
}
hand_b <- function(m) {
  if (inherits(m, "merMod")) {
    lme4::fixef(m)
  } else if (inherits(m, "coxme")) {
    coxme::fixef(m)
  } else {
    stats::coef(m)
  }
}
hand_V <- function(m) as.matrix(stats::vcov(m))

# Each fit and each hand fit is built once, on first use. The glmer hand fits
# use the package's own model frame at each pivot: its coordinates come from the
# closed composition and differ from pivotCoord() on the raw parts by about
# 1e-16, which moves glmer's Nelder-Mead optimum by about 1e-8 and its
# finite-difference Hessian covariance by about 1e-6. The coordinates themselves
# are checked against pivotCoord() below.
cache <- new.env()
fit_of <- function(nm) {
  key <- paste0("fit_", nm)
  if (is.null(cache[[key]])) cache[[key]] <- pkg_fit(specs[[nm]])
  cache[[key]]
}
hands_of <- function(nm) {
  key <- paste0("hand_", nm)
  if (is.null(cache[[key]])) {
    s <- specs[[nm]]
    same_frame <- s$engine %in% c("glmer_b", "glmer_p")
    cache[[key]] <- purrr::map(seq_along(s$p), function(k) {
      if (same_frame) {
        hand_fit(s, k, data = stats::model.frame(fit_of(nm)$models[[s$p[k]]]))
      } else {
        hand_fit(s, k)
      }
    })
  }
  cache[[key]]
}

moves_for <- function(p) {
  big <- if (length(p) == 3) "PA" else "MVPA"
  list(
    list(scheme = "one-to-one", from = "SB", to = big),
    list(scheme = "one-to-all", focal = big, role = "to"),
    list(scheme = "one-to-all", focal = "SB", role = "from")
  )
}
unit_of <- function(mv, ref) {
  u <- stats::setNames(numeric(length(ref)), names(ref))
  if (mv$scheme == "one-to-one") {
    u[[mv$from]] <- -1
    u[[mv$to]] <- 1
  } else {
    rest <- setdiff(names(ref), mv$focal)
    sgn <- if (mv$role == "to") 1 else -1
    u[[mv$focal]] <- sgn
    u[rest] <- -sgn * ref[rest] / sum(ref[rest])
  }
  u
}
pkg_realloc <- function(fit, mv, type = "relative") {
  if (mv$scheme == "one-to-one") {
    predict_realloc(fit, delta = deltas, from = mv$from, to = mv$to, type = type)
  } else if (mv$role == "to") {
    predict_realloc(fit, delta = deltas, scheme = "one-to-all", to = mv$focal,
                    type = type)
  } else {
    predict_realloc(fit, delta = deltas, scheme = "one-to-all", from = mv$focal,
                    type = type)
  }
}
# The contrast vector over the fixed effects for each delta.
contrast_rows <- function(m, s, mv) {
  ref <- hand_ref(s$p)
  zref <- pivot_row(ref, s$p)
  b <- hand_b(m)
  purrr::map(deltas, function(d) {
    L <- stats::setNames(numeric(length(b)), names(b))
    L[coord_names(s$p)] <- pivot_row(ref + d * unit_of(mv, ref), s$p) - zref
    L
  })
}
quad <- function(L, V) sqrt(as.numeric(t(L) %*% V[names(L), names(L)] %*% L))

test_that("the fitted pivot coordinates equal pivotCoord() at every pivot", {
  for (nm in c("cox", "cox_d3", "glmer_b")) {
    s <- specs[[nm]]
    fit <- fit_of(nm)
    expect_close(fit$data[coord_names(s$p)], hand_frame(1, s$p)[coord_names(s$p)], 1e-12)
    for (k in seq_along(s$p)[-1]) {
      mf <- stats::model.frame(fit$models[[s$p[k]]])
      expect_close(mf[coord_names(s$p)], hand_frame(k, s$p)[coord_names(s$p)], 1e-12)
    }
  }
  expect_close(fit_of("cox")$ref_comp[parts], hand_ref(parts), 1e-9)
  expect_close(fit_of("cox_d3")$ref_comp[parts3], hand_ref(parts3), 1e-9)
})

test_that("fixed effects and their covariance match hand fits on pivot coordinates", {
  for (nm in available(names(specs))) {
    s <- specs[[nm]]
    fit <- fit_of(nm)
    hands <- hands_of(nm)
    cf <- tc_model_coefs(fit$model, fit$model_type)
    expect_identical(names(cf$b), names(hand_b(hands[[1]])))
    expect_close(cf$b, hand_b(hands[[1]]), 1e-6)
    expect_close(cf$V, hand_V(hands[[1]])[names(cf$b), names(cf$b)], 1e-6)
    for (k in seq_along(s$p)) {
      b <- tc_model_coefs(fit$models[[s$p[k]]], fit$model_type)$b
      expect_close(b, hand_b(hands[[k]])[names(b)], 1e-6)
    }
  }
})

test_that("keeping the Cox model frame changes no estimate", {
  for (nm in c("cox", "cox_sex")) {
    a <- hands_of(nm)[[1]]
    b <- hand_fit(specs[[nm]], model_frame = FALSE)
    expect_identical(coef(a), coef(b))
    expect_identical(vcov(a), vcov(b))
    expect_identical(a$loglik, b$loglik)
  }
})

test_that("the omnibus is the chi-squared Wald test of the ILR block", {
  for (nm in available(names(specs))) {
    s <- specs[[nm]]
    m <- hands_of(nm)[[1]]
    cn <- coord_names(s$p)
    b <- hand_b(m)[cn]
    st <- as.numeric(t(b) %*% solve(hand_V(m)[cn, cn]) %*% b)
    om <- fit_of(nm)$omnibus
    expect_equal(om$df, length(cn))
    expect_close(om$statistic, st, 1e-6)
    expect_close(om$p.value, stats::pchisq(st, length(cn), lower.tail = FALSE), 1e-6)
  }
})

test_that("the omnibus equals car::linearHypothesis() on the hand fit", {
  skip_if_not_installed("car")
  for (nm in available(names(specs))) {
    m <- hands_of(nm)[[1]]
    lh <- car::linearHypothesis(m, paste(coord_names(specs[[nm]]$p), "= 0"))
    expect_close(fit_of(nm)$omnibus$statistic, lh$Chisq[2], 1e-6)
    expect_close(fit_of(nm)$omnibus$p.value, lh[["Pr(>Chisq)"]][2], 1e-6)
  }
})

lmer_specs <- c("lmer", "lmer_sex", "lmer_d3")

test_that("per-part effects are the asymptotic Wald effects of the hand fit at each pivot", {
  for (nm in available(setdiff(names(specs), lmer_specs))) {
    s <- specs[[nm]]
    for (lev in c(0.95, 0.90)) {
      pp <- summary(fit_of(nm), conf_level = lev)$per_part
      expect_identical(pp$parameter, paste(s$p, "vs remaining"))
      crit <- stats::qnorm(1 - (1 - lev) / 2)
      ref <- do.call(rbind, purrr::map(hands_of(nm), function(m) {
        est <- hand_b(m)[["ilr1"]]
        se <- sqrt(hand_V(m)["ilr1", "ilr1"])
        c(exp(c(est, est - crit * se, est + crit * se)),
          2 * stats::pnorm(-abs(est / se)))
      }))
      expect_close(pp$estimate, ref[, 1], 1e-6)
      expect_close(cbind(pp$conf.low, pp$conf.high), ref[, 2:3], 1e-6)
      expect_close(pp$p.value, ref[, 4], 1e-6)
    }
  }
})

test_that("lmer per-part effects carry the Kenward-Roger SE and df of the hand fit", {
  skip_if_not_installed("pbkrtest")
  for (nm in lmer_specs) {
    s <- specs[[nm]]
    for (lev in c(0.95, 0.90)) {
      pp <- summary(fit_of(nm), conf_level = lev)$per_part
      expect_identical(pp$parameter, paste(s$p, "vs remaining"))
      ref <- do.call(rbind, purrr::map(hands_of(nm), function(m) {
        b <- hand_b(m)
        L <- matrix(as.numeric(names(b) == "ilr1"), 1)
        est <- b[["ilr1"]]
        se <- sqrt(as.numeric(L %*% as.matrix(pbkrtest::vcovAdj(m)) %*% t(L)))
        df <- pbkrtest::get_Lb_ddf(m, L)
        crit <- stats::qt(1 - (1 - lev) / 2, df)
        c(est, est - crit * se, est + crit * se, 2 * stats::pt(-abs(est / se), df))
      }))
      expect_close(pp$estimate, ref[, 1], 1e-6)
      expect_close(cbind(pp$conf.low, pp$conf.high), ref[, 2:3], 1e-6)
      expect_close(pp$p.value, ref[, 4], 1e-6)
    }
    # The same SE and df emmeans reports for the ilr1 trend of each rotated fit.
    for (k in seq_along(s$p)) {
      m <- fit_of(nm)$models[[s$p[k]]]
      st <- tc_iv_ilr1_stats(m, "lmer")
      et <- summary(emmeans::emtrends(m, ~ 1, var = "ilr1"))
      expect_close(c(st$se, st$df), c(et$SE, et$df), 1e-6)
    }
  }
})

test_that("ratio-scale reallocations equal exp(d'b) with SE sqrt(d'Vd) on the link scale", {
  for (nm in available(c("cox", "cox_sex", "cox_d3", "glmer_b", "glmer_p", "coxme"))) {
    s <- specs[[nm]]
    m <- hands_of(nm)[[1]]
    b <- hand_b(m)
    V <- hand_V(m)
    for (mv in moves_for(s$p)) {
      pr <- pkg_realloc(fit_of(nm), mv)
      expect_identical(pr$delta, deltas)
      Ls <- contrast_rows(m, s, mv)
      est <- purrr::map_dbl(Ls, function(L) sum(L * b))
      se <- purrr::map_dbl(Ls, quad, V = V)
      z <- stats::qnorm(0.975)
      expect_close(pr$estimate, exp(est), 1e-6)
      expect_close(pr$std.error, exp(est) * se, 1e-6)
      expect_close(cbind(pr$conf.low, pr$conf.high),
                   cbind(exp(est - z * se), exp(est + z * se)), 1e-6)
      expect_close(pr$statistic, est / se, 1e-6)
      expect_close(pr$p.value, 2 * stats::pnorm(-abs(est / se)), 1e-6)
      expect_true(all(is.infinite(pr$df)))
    }
  }
})

test_that("lmer reallocations use Kenward-Roger df and the adjusted covariance", {
  skip_if_not_installed("pbkrtest")
  for (nm in c("lmer", "lmer_sex", "lmer_d3")) {
    s <- specs[[nm]]
    m <- hands_of(nm)[[1]]
    b <- hand_b(m)
    Va <- as.matrix(pbkrtest::vcovAdj(m))
    V0 <- hand_V(m)
    for (mv in moves_for(s$p)) {
      pr <- pkg_realloc(fit_of(nm), mv)
      Ls <- contrast_rows(m, s, mv)
      est <- purrr::map_dbl(Ls, function(L) sum(L * b))
      se <- purrr::map_dbl(Ls, quad, V = Va)
      df <- purrr::map_dbl(Ls, function(L) pbkrtest::get_Lb_ddf(m, matrix(L, 1)))
      crit <- stats::qt(0.975, df)
      expect_close(pr$estimate, est, 1e-6)
      expect_close(pr$std.error, se, 1e-6)
      expect_close(pr$df, df, 1e-6)
      expect_close(cbind(pr$conf.low, pr$conf.high),
                   cbind(est - crit * se, est + crit * se), 1e-6)
      expect_close(pr$p.value, 2 * stats::pt(-abs(est / se), df), 1e-6)
      # Asymptotic df gives the plain covariance and df = Inf.
      old <- emmeans::get_emm_option("lmer.df")
      emmeans::emm_options(lmer.df = "asymptotic")
      pr_as <- pkg_realloc(fit_of(nm), mv)
      emmeans::emm_options(lmer.df = old)
      expect_true(all(is.infinite(pr_as$df)))
      expect_close(pr_as$std.error, purrr::map_dbl(Ls, quad, V = V0), 1e-6)
      expect_close(pr_as$estimate, est, 1e-6)
    }
  }
})

test_that("lmer reallocations under Satterthwaite df match lmerTest", {
  skip_if_not_installed("lmerTest")
  old <- emmeans::get_emm_option("lmer.df")
  emmeans::emm_options(lmer.df = "satterthwaite")
  on.exit(emmeans::emm_options(lmer.df = old), add = TRUE)
  for (nm in c("lmer", "lmer_sex", "lmer_d3")) {
    s <- specs[[nm]]
    m <- hands_of(nm)[[1]]
    # A hand lmerTest fit of the same frame; its Satterthwaite SE is the plain
    # model covariance.
    mt <- lmerTest::lmer(stats::formula(m), data = m@frame)
    for (mv in moves_for(s$p)) {
      pr <- pkg_realloc(fit_of(nm), mv)
      ct <- purrr::map(contrast_rows(m, s, mv), function(L) {
        lmerTest::contest(mt, matrix(L, 1), joint = FALSE)
      })
      expect_close(pr$estimate, purrr::map_dbl(ct, "Estimate"), 1e-6)
      expect_close(pr$std.error, purrr::map_dbl(ct, "Std. Error"), 1e-6)
      expect_close(pr$df, purrr::map_dbl(ct, "df"), 1e-6)
    }
  }
})

# glmer is left out, and coxme's omnibus: refitting in another basis is a
# different iterative optimization, whose optimum (same log-likelihood to 1e-9)
# is reproduced only to the optimizer's tolerance.
test_that("reallocation estimates do not depend on the order of the parts", {
  for (nm in available(c("cox", "cox_sex", "cox_d3", "lmer", "lmer_sex", "lmer_d3", "coxme"))) {
    s <- specs[[nm]]
    fr <- pkg_fit(s, reverse = TRUE)
    if (nm != "coxme") {
      expect_close(fr$omnibus$statistic, fit_of(nm)$omnibus$statistic, 1e-8)
    }
    for (mv in moves_for(s$p)) {
      cols <- c("estimate", "std.error", "conf.low", "conf.high", "p.value")
      expect_close(pkg_realloc(fit_of(nm), mv)[cols], pkg_realloc(fr, mv)[cols], 1e-8)
    }
  }
})

test_that("absolute reallocations equal population-level predictions of the hand fit", {
  for (nm in c("lmer", "lmer_sex", "lmer_d3", "glmer_b", "glmer_p")) {
    s <- specs[[nm]]
    m <- hands_of(nm)[[1]]
    mv <- moves_for(s$p)[[1]]
    pa <- pkg_realloc(fit_of(nm), mv, type = "absolute")
    expect_true(is.na(pa$from[1]))
    expect_equal(pa$delta, c(0, deltas))
    ref <- hand_ref(s$p)
    comps <- c(list(ref), purrr::map(deltas, function(d) ref + d * unit_of(mv, ref)))
    nd <- as.data.frame(do.call(rbind, purrr::map(comps, pivot_row, p = s$p)))
    # Factor covariates are averaged with equal weights on the link scale.
    levs <- if (is.null(s$cov)) list(NULL) else as.list(levels(tu[[s$cov]]))
    eta <- rowMeans(do.call(cbind, purrr::map(levs, function(l) {
      if (!is.null(l)) nd[[s$cov]] <- factor(l, levels = levels(tu[[s$cov]]))
      stats::predict(m, newdata = nd, re.form = NA, type = "link")
    })))
    fam <- stats::family(m)
    expect_close(pa$estimate, fam$linkinv(eta), 1e-6)
    if (is.null(s$cov)) {
      expect_close(pa$estimate,
                   stats::predict(m, newdata = nd, re.form = NA, type = "response"),
                   1e-6)
    }
    if (s$engine != "lmer") {
      X <- cbind(1, as.matrix(nd))
      se <- sqrt(rowSums((X %*% hand_V(m)) * X))
      z <- stats::qnorm(0.975)
      expect_close(cbind(pa$conf.low, pa$conf.high),
                   cbind(fam$linkinv(eta - z * se), fam$linkinv(eta + z * se)), 1e-6)
      expect_close(pa$std.error, se * fam$mu.eta(eta), 1e-6)
      expect_true(all(is.infinite(pa$df)))
    }
  }
})

test_that("absolute lmer predictions use Kenward-Roger df per row", {
  skip_if_not_installed("pbkrtest")
  for (nm in c("lmer", "lmer_sex", "lmer_d3")) {
    s <- specs[[nm]]
    m <- hands_of(nm)[[1]]
    b <- hand_b(m)
    Va <- as.matrix(pbkrtest::vcovAdj(m))
    mv <- moves_for(s$p)[[1]]
    pa <- pkg_realloc(fit_of(nm), mv, type = "absolute")
    ref <- hand_ref(s$p)
    comps <- c(list(ref), purrr::map(deltas, function(d) ref + d * unit_of(mv, ref)))
    Ls <- purrr::map(comps, function(x) {
      L <- stats::setNames(numeric(length(b)), names(b))
      L[["(Intercept)"]] <- 1
      L[coord_names(s$p)] <- pivot_row(x, s$p)
      if (!is.null(s$cov)) L[["sexmale"]] <- 0.5
      L
    })
    est <- purrr::map_dbl(Ls, function(L) sum(L * b))
    se <- purrr::map_dbl(Ls, quad, V = Va)
    df <- purrr::map_dbl(Ls, function(L) pbkrtest::get_Lb_ddf(m, matrix(L, 1)))
    crit <- stats::qt(0.975, df)
    expect_close(pa$df, df, 1e-6)
    expect_close(pa$std.error, se, 1e-6)
    expect_close(cbind(pa$conf.low, pa$conf.high),
                 cbind(est - crit * se, est + crit * se), 1e-6)
  }
})

test_that("the Cox diagnostic rows are survival::cox.zph() of the hand fit", {
  for (nm in c("cox", "cox_sex", "cox_d3")) {
    tests <- fit_check(fit_of(nm), plot = FALSE)$tests
    z <- survival::cox.zph(hands_of(nm)[[1]])$table
    ph <- tests[grepl("^Proportional hazards", tests$check), ]
    expect_identical(ph$check, paste0("Proportional hazards: ", rownames(z)))
    expect_close(ph$statistic, z[, "chisq"], 1e-6)
    expect_close(ph$df, z[, "df"], 0)
    expect_close(ph$p.value, z[, "p"], 1e-6)
  }
})

test_that("the coxme fit statistics are AIC() and the penalized logLik() of the hand fit", {
  skip_if_not_installed("coxme")
  tests <- fit_check(fit_of("coxme"), plot = FALSE)$tests
  h <- hands_of("coxme")[[1]]
  expect_close(tests$statistic[tests$check == "AIC"], stats::AIC(h), 1e-6)
  expect_close(tests$statistic[tests$check == "log-likelihood"],
               as.numeric(stats::logLik(h, type = "penalized")), 1e-6)
})

test_that("the mixed-model DHARMa rows equal DHARMa on the hand fit with the same seed", {
  skip_if_not_installed("DHARMa")
  for (nm in c("lmer", "lmer_sex", "glmer_b", "glmer_p")) {
    tests <- fit_check(fit_of(nm), plot = FALSE)$tests
    # The package simulates and tests under set.seed(123); DHARMa reseeds its
    # simulation with 123 too, and the bootstrapped outlier test then draws
    # from that same stream.
    set.seed(123)
    sr <- DHARMa::simulateResiduals(hands_of(nm)[[1]], n = 250, seed = 123,
                                    plot = FALSE)
    u <- DHARMa::testUniformity(sr, plot = FALSE)
    d <- DHARMa::testDispersion(sr, plot = FALSE)
    o <- DHARMa::testOutliers(sr, plot = FALSE)
    got <- tests[tests$check %in% c("KS uniformity", "Dispersion", "Outliers"), ]
    expect_identical(got$check, c("KS uniformity", "Dispersion", "Outliers"))
    expect_close(cbind(got$statistic, got$p.value),
                 cbind(c(u$statistic, d$statistic, o$statistic),
                       c(u$p.value, d$p.value, o$p.value)), 1e-6)
  }
})

test_that("a starved glmer reports the hand fit's convergence state and messages", {
  ctrl <- lme4::glmerControl(optCtrl = list(maxfun = 5))
  fit <- suppressWarnings(
    comp_iv(declare_composition(tu, parts, closure = closure), "obesity",
            model = "glmer", family = "binomial", random = "school",
            model_args = list(control = ctrl))
  )
  captured <- list()
  hands <- purrr::map(seq_along(parts), function(k) {
    withCallingHandlers(
      hand_fit(specs$glmer_b, k, control = ctrl),
      warning = function(w) {
        captured[[length(captured) + 1]] <<- list(k = k, msg = conditionMessage(w))
        invokeRestart("muffleWarning")
      },
      message = function(m) {
        captured[[length(captured) + 1]] <<- list(
          k = k, msg = sub("\n$", "", conditionMessage(m))
        )
        invokeRestart("muffleMessage")
      }
    )
  })
  h1 <- hands[[1]]
  expect_identical(fit$convergence$converged,
                   length(h1@optinfo$conv$lme4$messages) == 0)
  expect_identical(fit$convergence$singular, lme4::isSingular(h1))
  expect_equal(fit$convergence$n_warnings, length(captured))
  canonical <- purrr::keep(captured, function(x) x$k == 1)
  expect_identical(fit$fit_warnings$message[fit$fit_warnings$source == "canonical"],
                   purrr::map_chr(canonical, "msg"))
  expect_close(lme4::fixef(fit$model), lme4::fixef(h1), 1e-6)
})
