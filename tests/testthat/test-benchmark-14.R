# Numerical benchmark of the optional comp_iv() engines: clm and clmm
# (ordinal), mblogit (multinomial baseline-logit), lmrob and rlmer (robust
# linear and robust linear mixed). Each fit equals its backend fitted by hand
# on robCompositions::pivotCoord() coordinates; MASS::polr() and
# nnet::multinom() agree as independent implementations; the omnibus is the
# joint Wald test of the coordinate coefficients; each per-part summary() row
# is the ilr1 coefficient of the hand refit on that part's pivot; a relative
# reallocation is exp(d'b) or d'b with link-scale standard error
# sqrt(d' V d); the result does not depend on the part order; fit_check()
# reproduces ordinal::nominal_test(), sure::resids() and the robustness
# weights. lmrob starts from random subsamples, so every lmrob comparison
# fits both sides under the same set.seed(). Reference calls are
# namespace-qualified so each number is visibly computed outside timecoda.

parts <- c("MVPA", "LPA", "SB", "SPT")
D <- length(parts)
coord_names <- function(p) paste0("ilr", seq_len(length(p) - 1))
coords <- coord_names(parts)
closure <- 1440
seed_lr <- 20261016

# timeuse with its ordered `grade`, the same bands unordered, and a numeric
# outcome whose first ten rows are shifted so the robust fits downweight them.
set.seed(1)
tu <- timeuse
lp <- 0.9 * log(tu$MVPA / tu$SB) - 0.4 * log(tu$LPA / tu$SPT) +
  as.numeric(tu$school) * 0.1
tu$rz <- lp + stats::rnorm(nrow(tu))
outlier_rows <- 1:10
tu$rz[outlier_rows] <- tu$rz[outlier_rows] + 8
tu$grade_unordered <- factor(as.character(tu$grade), levels = levels(tu$grade))
cd <- declare_composition(tu, parts, closure = closure)

comp_mat <- as.matrix(tu[parts])
comp_mat <- comp_mat / rowSums(comp_mat) * closure
piv <- function(x, k, p = colnames(x)) {
  z <- as.matrix(robCompositions::pivotCoord(as.data.frame(x), pivotvar = k))
  colnames(z) <- coord_names(p)
  z
}
hand_frames <- lapply(seq_len(D), function(k) cbind(tu, piv(comp_mat, k)))
d <- hand_frames[[1]]

# The closed geometric mean and a 15-minute SB -> MVPA move from it.
gm <- exp(colMeans(log(comp_mat)))
ref_hand <- gm / sum(gm) * closure
new_hand <- ref_hand
new_hand[["SB"]] <- new_hand[["SB"]] - 15
new_hand[["MVPA"]] <- new_hand[["MVPA"]] + 15
zz <- piv(rbind(ref_hand, new_hand), 1, parts)
dz <- stats::setNames(zz[2, ] - zz[1, ], coords)

max_diff <- function(a, b) {
  a <- as.numeric(unlist(a))
  b <- as.numeric(unlist(b))
  if (length(a) != length(b) || !identical(is.na(a), is.na(b))) return(Inf)
  ok <- !is.na(a)
  if (!identical(a[ok & !is.finite(a)], b[ok & !is.finite(a)])) return(Inf)
  ok <- ok & is.finite(a)
  max(c(0, abs(a[ok] - b[ok])))
}
coefs <- function(m) {
  if (inherits(m, c("rlmerMod", "merMod"))) lme4::fixef(m) else stats::coef(m)
}
vc <- function(m, nm = names(coefs(m))) as.matrix(stats::vcov(m))[nm, nm, drop = FALSE]
wald <- function(b, V, nm) {
  as.numeric(t(b[nm]) %*% solve(V[nm, nm, drop = FALSE]) %*% b[nm])
}
ilr_coefs <- function(nm) nm[grepl("(^|~)ilr[0-9]+$", nm)]
quiet <- function(expr) suppressWarnings(suppressMessages(expr))
# The backend arguments the fitted object's stored call carries beyond formula
# and data, so a hand refit uses the same settings.
call_extras <- function(m) {
  cl <- as.list(stats::getCall(m))[-1]
  lapply(cl[setdiff(names(cl), c("formula", "data"))], eval, envir = globalenv())
}
# The data travel as a value inside the call, so nominal_test() can refit.
hand_fit <- function(fun, formula, data, extras = list()) {
  do.call(fun, c(list(formula = formula, data = data), extras))
}
# The frame a fitted model was given, as stored in its call.
pkg_frame <- function(m) as.data.frame(eval(stats::getCall(m)$data))

f_ord <- grade ~ ilr1 + ilr2 + ilr3
f_ordm <- grade ~ ilr1 + ilr2 + ilr3 + (1 | school)
f_mb <- grade_unordered ~ ilr1 + ilr2 + ilr3
f_rz <- rz ~ ilr1 + ilr2 + ilr3
f_rzm <- rz ~ ilr1 + ilr2 + ilr3 + (1 | school)

fit_cache <- new.env()
fit_of <- function(name) {
  if (!is.null(fit_cache[[name]])) return(fit_cache[[name]])
  fit_cache[[name]] <- switch(name,
    o = comp_iv(cd, "grade", model = "clm"),
    o_cov = comp_iv(cd, "grade", model = "clm", covariates = "sex"),
    om = comp_iv(cd, "grade", model = "clmm", random = "school"),
    mb = quiet(comp_iv(cd, "grade_unordered", model = "mblogit")),
    mb_cl = quiet(comp_iv(cd, "grade_unordered", model = "mblogit",
                          random = "school")),
    lr = {
      set.seed(seed_lr)
      comp_iv(cd, "rz", model = "lmrob")
    },
    lr_cov = {
      set.seed(seed_lr)
      comp_iv(cd, "rz", model = "lmrob", covariates = "sex")
    },
    rl = comp_iv(cd, "rz", model = "rlmer", random = "school")
  )
}
# The hand refit of each fit on pivot k. clmm's covariance is a
# finite-difference Hessian that moves in the 7th digit when the coordinates
# move in the 16th, so its refits use the package's own frames, which equal
# the pivotCoord() columns to 1e-12 (checked below).
hand_cache <- new.env()
hand_of <- function(name, k = 1) {
  key <- paste0(name, "_", k)
  if (!is.null(hand_cache[[key]])) return(hand_cache[[key]])
  f <- fit_of(name)
  ex <- call_extras(f$model)
  fr <- hand_frames[[k]]
  hand_cache[[key]] <- switch(name,
    o = hand_fit(ordinal::clm, f_ord, fr, ex),
    o_cov = hand_fit(ordinal::clm, grade ~ ilr1 + ilr2 + ilr3 + sex, fr, ex),
    om = hand_fit(ordinal::clmm, f_ordm, pkg_frame(f$models[[k]]), ex),
    mb = quiet(hand_fit(mclogit::mblogit, f_mb, fr, ex)),
    mb_cl = quiet(hand_fit(mclogit::mblogit, f_mb, fr, ex)),
    lr = {
      # comp_iv() fits pivot 1, then pivots 2..D, on one stream.
      set.seed(seed_lr)
      fits <- lapply(hand_frames, function(x) hand_fit(robustbase::lmrob, f_rz, x, ex))
      for (j in seq_len(D)) hand_cache[[paste0(name, "_", j)]] <- fits[[j]]
      fits[[k]]
    },
    lr_cov = {
      set.seed(seed_lr)
      hand_fit(robustbase::lmrob, rz ~ ilr1 + ilr2 + ilr3 + sex, fr, ex)
    },
    rl = hand_fit(robustlmm::rlmer, f_rzm, fr, ex)
  )
}
engine_pkgs <- c(o = "ordinal", o_cov = "ordinal", om = "ordinal", mb = "mclogit",
                 mb_cl = "mclogit", lr = "robustbase", lr_cov = "robustbase",
                 rl = "robustlmm")
tol <- c(o = 1e-8, o_cov = 1e-8, om = 1e-6, mb = 1e-8, mb_cl = 1e-6, lr = 1e-8,
         lr_cov = 1e-8, rl = 1e-6)
ratio <- c(o = TRUE, om = TRUE, mb = TRUE, lr = FALSE, rl = FALSE)
skip_engine <- function(name) {
  skip_on_cran()
  skip_if_not_installed("robCompositions")
  skip_if_not_installed(engine_pkgs[[name]])
}

test_that("the fit coordinates are the pivotCoord() coordinates", {
  skip_on_cran()
  skip_if_not_installed("robCompositions")
  skip_if_not_installed("ordinal")
  for (k in seq_len(D)) {
    it <- ilr_transform(tu, parts, total = closure, pivot = k)
    expect_lt(max_diff(as.matrix(it[coords]), hand_frames[[k]][coords]), 1e-12)
    expect_lt(max_diff(as.matrix(pkg_frame(fit_of("om")$models[[k]])[coords]),
                       hand_frames[[k]][coords]), 1e-12)
  }
})

test_that("every engine's fit equals its backend fitted by hand", {
  for (nm in names(tol)) skip_engine(nm)
  for (nm in names(tol)) {
    f <- fit_of(nm)
    h <- hand_of(nm)
    b <- coefs(f$model)
    expect_setequal(names(b), names(coefs(h)))
    expect_lt(max_diff(b, coefs(h)[names(b)]), tol[[nm]])
    expect_lt(max_diff(vc(f$model), vc(h, names(b))), tol[[nm]])
  }
  expect_lt(max_diff(fit_of("o")$model$alpha, hand_of("o")$alpha), 1e-8)
  expect_true(all(grepl("^(expected|above)~", names(stats::coef(fit_of("mb")$model)))))
  expect_lt(max_diff(fit_of("lr")$model$rweights, hand_of("lr")$rweights), 1e-8)
  # rlmer is deterministic.
  h2 <- hand_fit(robustlmm::rlmer, f_rzm, d, call_extras(fit_of("rl")$model))
  expect_identical(lme4::fixef(h2), lme4::fixef(hand_of("rl")))
  # Its intraclass correlation is var_g / (var_g + sigma^2) of the hand fit.
  vg <- as.numeric(lme4::VarCorr(hand_of("rl"))$school[1, 1])
  expect_lt(max_diff(summary(fit_of("rl"))$icc$value,
                     vg / (vg + stats::sigma(hand_of("rl"))^2)), 1e-6)
  # A clustered mblogit passes a random slope through as its own `random`.
  sl <- quiet(comp_iv(cd, "grade_unordered", model = "mblogit",
                      random = "(ilr1 | school)"))
  h_sl <- quiet(hand_fit(mclogit::mblogit, f_mb, d,
                         list(random = ~ ilr1 | school, trace = FALSE)))
  expect_lt(max_diff(coefs(sl$model), coefs(h_sl)[names(coefs(sl$model))]), 1e-6)
  expect_lt(max_diff(vc(sl$model), vc(h_sl, names(coefs(sl$model)))), 1e-6)
})

test_that("MASS::polr() and nnet::multinom() agree with clm and mblogit", {
  skip_engine("o")
  skip_engine("mb")
  skip_if_not_installed("MASS")
  skip_if_not_installed("nnet")
  o <- fit_of("o")$model
  # Both parameterize logit P(Y <= j) = alpha_j - x'beta: same signs.
  pl <- MASS::polr(f_ord, data = d, Hess = TRUE)
  expect_lt(max_diff(stats::coef(pl)[coords], o$beta[coords]), 1e-3)
  expect_lt(max_diff(unname(pl$zeta), unname(o$alpha)), 1e-3)
  se_pl <- sqrt(diag(stats::vcov(pl)))
  se_o <- sqrt(diag(vc(o)))
  expect_lt(max_diff(se_pl[coords], se_o[coords]), 1e-3)
  expect_lt(max_diff(unname(se_pl[names(pl$zeta)]), unname(se_o[names(o$alpha)])), 1e-3)
  expect_lt(max_diff(stats::predict(pl, newdata = d, type = "probs"),
                     stats::predict(o, newdata = d[coords], type = "prob")$fit), 1e-3)
  # Both use the first level ("below") as the reference category.
  mn <- nnet::multinom(f_mb, data = d, trace = FALSE, Hess = TRUE, reltol = 1e-12,
                       maxit = 1000)
  cm <- stats::coef(mn)
  b <- stats::coef(fit_of("mb")$model)
  se <- sqrt(diag(vc(fit_of("mb")$model)))
  pick <- function(x) t(sapply(rownames(cm), function(l) x[paste0(l, "~", colnames(cm))]))
  expect_identical(mn$lev[1], levels(d$grade_unordered)[1])
  expect_lt(max_diff(pick(b), cm), 1e-3)
  expect_lt(max_diff(pick(se), summary(mn)$standard.errors), 1e-3)
})

test_that("the omnibus is the joint Wald test of the coordinate coefficients", {
  K <- nlevels(tu$grade_unordered)
  for (nm in names(tol)) skip_engine(nm)
  for (nm in names(tol)) {
    f <- fit_of(nm)
    h <- hand_of(nm)
    nm_ilr <- ilr_coefs(names(coefs(h)))
    want_df <- if (nm %in% c("mb", "mb_cl")) (K - 1) * (D - 1) else D - 1
    expect_equal(f$omnibus$df, want_df)
    expect_length(nm_ilr, want_df)
    w <- wald(coefs(h), vc(h), nm_ilr)
    expect_lt(max_diff(f$omnibus$statistic, w), tol[[nm]])
    expect_lt(max_diff(f$omnibus$p.value, stats::pchisq(w, want_df, lower.tail = FALSE)),
              tol[[nm]])
  }
})

test_that("each per-part summary row is the ilr1 coefficient of the pivot-k refit", {
  for (nm in names(ratio)) skip_engine(nm)
  for (nm in names(ratio)) {
    f <- fit_of(nm)
    pp <- summary(f)$per_part
    levs <- if (nm == "mb") c("expected", "above") else NA_character_
    for (k in seq_len(D)) {
      h <- hand_of(nm, k)
      expect_lt(max_diff(coefs(f$models[[k]]), coefs(h)[names(coefs(f$models[[k]]))]),
                tol[[nm]])
      for (lev in levs) {
        cn <- if (is.na(lev)) "ilr1" else paste0(lev, "~ilr1")
        row <- if (is.na(lev)) pp[k, ] else
          pp[pp$parameter == paste(parts[k], "vs remaining") & pp$response_level == lev, ]
        b1 <- coefs(h)[[cn]]
        se1 <- sqrt(vc(h)[cn, cn])
        expect_lt(max_diff(row$estimate, if (ratio[[nm]]) exp(b1) else b1), tol[[nm]])
        # clm, clmm and mblogit rows are asymptotic; lmrob's use its
        # residual df, as summary.lmrob() does.
        if (nm %in% c("o", "om", "mb")) {
          ci <- exp(b1 + c(-1, 1) * stats::qnorm(0.975) * se1)
          expect_lt(max_diff(c(row$conf.low, row$conf.high), ci), tol[[nm]])
        }
        if (nm == "lr") {
          ci <- b1 + c(-1, 1) * stats::qt(0.975, stats::df.residual(h)) * se1
          expect_lt(max_diff(c(row$conf.low, row$conf.high), ci), 1e-8)
          expect_lt(max_diff(row$p.value,
                             stats::coef(summary(h))["ilr1", "Pr(>|t|)"]), 1e-8)
        }
        # Every row's interval and p-value rest on one distribution: the
        # interval is centred on b, and the df it implies gives the p-value.
        lo <- if (ratio[[nm]]) log(row$conf.low) else row$conf.low
        hi <- if (ratio[[nm]]) log(row$conf.high) else row$conf.high
        expect_lt(abs((lo + hi) / 2 - b1), tol[[nm]])
        q <- (hi - lo) / (2 * se1)
        df_row <- if (abs(q - stats::qnorm(0.975)) < 1e-6) Inf else
          stats::uniroot(function(v) stats::qt(0.975, v) - q, c(0.5, 1e7),
                         tol = 1e-12)$root
        expect_lt(max_diff(row$p.value, 2 * stats::pt(-abs(b1 / se1), df_row)), 1e-6)
      }
    }
  }
})

test_that("a relative reallocation is the hand contrast on the link scale", {
  check_row <- function(row, b, V, is_ratio, tl) {
    est <- sum(dz * b[coords])
    se <- sqrt(as.numeric(t(dz) %*% V[coords, coords] %*% dz))
    ci <- est + c(-1, 1) * stats::qt(0.975, row$df) * se
    expect_lt(max_diff(row$estimate, if (is_ratio) exp(est) else est), tl)
    expect_lt(max_diff(if (is_ratio) row$std.error / row$estimate else row$std.error, se), tl)
    expect_lt(max_diff(c(row$conf.low, row$conf.high), if (is_ratio) exp(ci) else ci), tl)
    expect_lt(max_diff(row$p.value, 2 * stats::pt(-abs(est / se), row$df)), tl)
  }
  for (nm in names(ratio)) skip_engine(nm)
  for (nm in names(ratio)) {
    f <- fit_of(nm)
    h <- hand_of(nm)
    expect_lt(max_diff(unlist(f$ref_comp)[parts], ref_hand[parts]), 1e-9)
    pr <- predict_realloc(f, delta = 15, from = "SB", to = "MVPA")
    if (nm == "mb") {
      expect_setequal(pr$response_level, c("expected", "above"))
      for (lev in c("expected", "above")) {
        idx <- grep(paste0("^", lev, "~"), names(coefs(h)))
        b_c <- stats::setNames(coefs(h)[idx], sub("^[^~]*~", "", names(coefs(h))[idx]))
        V_c <- vc(h)[idx, idx]
        dimnames(V_c) <- list(names(b_c), names(b_c))
        check_row(pr[pr$response_level == lev, ], b_c, V_c, TRUE, 1e-8)
      }
    } else {
      check_row(pr, coefs(h), vc(h), ratio[[nm]], tol[[nm]])
    }
    # A zero move is the reference value, with no test.
    p0 <- predict_realloc(f, delta = 0, from = "SB", to = "MVPA")
    expect_true(all(p0$estimate == if (ratio[[nm]]) 1 else 0))
    expect_true(all(p0$std.error == 0))
    expect_true(all(is.na(p0$statistic) & is.na(p0$p.value) & is.na(p0$df)))
  }
})

test_that("absolute predictions: refused per category, x'b for the robust engines", {
  for (nm in c("o", "om", "mb")) skip_engine(nm)
  for (nm in c("o", "om", "mb")) {
    expect_error(predict_realloc(fit_of(nm), delta = 15, from = "SB", to = "MVPA",
                                 type = "absolute"), "per outcome category")
  }
  X <- cbind(1, zz)
  colnames(X) <- c("(Intercept)", coords)
  for (nm in c("lr", "rl")) skip_engine(nm)
  for (nm in c("lr", "rl")) {
    pa <- predict_realloc(fit_of(nm), delta = 15, from = "SB", to = "MVPA",
                          type = "absolute")
    h <- hand_of(nm)
    expect_lt(max_diff(pa$estimate, X %*% coefs(h)[colnames(X)]), tol[[nm]])
    expect_lt(max_diff(pa$std.error,
                       sqrt(rowSums((X %*% vc(h, colnames(X))) * X))), tol[[nm]])
  }
})

test_that("reallocations do not depend on the part order", {
  orders <- list(c("SB", "SPT", "MVPA", "LPA"), c("SPT", "LPA", "MVPA", "SB"))
  by_level <- function(x) {
    if ("response_level" %in% names(x)) x[order(x$response_level), ] else x
  }
  for (nm in c("o", "mb", "lr", "rl")) skip_engine(nm)
  for (nm in c("o", "mb", "lr", "rl")) {
    ref <- by_level(predict_realloc(fit_of(nm), delta = 15, from = "SB", to = "MVPA"))
    for (p in orders) {
      cd_o <- declare_composition(tu, p, closure = closure)
      f <- switch(nm,
        o = comp_iv(cd_o, "grade", model = "clm"),
        mb = quiet(comp_iv(cd_o, "grade_unordered", model = "mblogit")),
        lr = {
          set.seed(seed_lr)
          comp_iv(cd_o, "rz", model = "lmrob")
        },
        rl = comp_iv(cd_o, "rz", model = "rlmer", random = "school")
      )
      pr <- by_level(predict_realloc(f, delta = 15, from = "SB", to = "MVPA"))
      tl <- if (nm == "rl") 1e-6 else 1e-8
      expect_lt(max_diff(pr$estimate, ref$estimate), tl)
      expect_lt(max_diff(pr$std.error, ref$std.error), tl)
      # lmrob's iterations stop at rel.tol = 1e-7, so its omnibus (about 525)
      # is compared under a tight tolerance below.
      if (nm != "lr") {
        expect_lt(max_diff(f$omnibus$statistic, fit_of(nm)$omnibus$statistic), tl)
      }
    }
  }
  ctl <- robustbase::lmrob.control(rel.tol = 1e-12)
  tight <- lapply(c(list(parts), orders), function(p) {
    set.seed(seed_lr)
    comp_iv(declare_composition(tu, p, closure = closure), "rz", model = "lmrob",
            model_args = list(control = ctl))
  })
  for (i in 2:3) {
    expect_lt(max_diff(tight[[i]]$omnibus$statistic, tight[[1]]$omnibus$statistic), 1e-8)
    expect_lt(max_diff(stats::fitted(tight[[i]]$model), stats::fitted(tight[[1]]$model)), 1e-8)
  }
})

test_that("fit_check() reproduces nominal_test(), sure::resids() and the weights", {
  skip_engine("o")
  skip_if_not_installed("sure")
  fc <- fit_check(fit_of("o"), plot = FALSE)
  nt <- ordinal::nominal_test(hand_of("o"))
  keep <- which(!is.na(nt[["Pr(>Chi)"]]))
  po <- fc$tests[grepl("^Proportional odds: ", fc$tests$check), ]
  expect_identical(po$check, paste0("Proportional odds: ", rownames(nt)[keep]))
  expect_lt(max_diff(po$statistic, nt$LRT[keep]), 1e-8)
  expect_lt(max_diff(po$df, nt$Df[keep]), 1e-8)
  expect_lt(max_diff(po$p.value, nt[["Pr(>Chi)"]][keep]), 1e-8)
  # The surrogate draws run under seed 123 and leave the session's stream as
  # it was.
  set.seed(9)
  before <- .Random.seed
  fc2 <- fit_check(fit_of("o"), plot = FALSE)
  expect_identical(.Random.seed, before)
  expect_identical(fc2$raw$surrogate, fc$raw$surrogate)
  set.seed(123)
  expect_lt(max_diff(fc$raw$surrogate, sure::resids(fit_of("o")$model)), 1e-8)

  for (nm in c("lr", "rl")) skip_engine(nm)
  for (nm in c("lr", "rl")) {
    fc <- fit_check(fit_of(nm), plot = FALSE)
    h <- hand_of(nm)
    w <- if (nm == "lr") h$rweights else robustlmm::getME(h, "w_e")
    sw <- stats::shapiro.test(stats::residuals(h))
    st <- function(lbl) fc$tests[grepl(lbl, fc$tests$check, fixed = TRUE), ]
    expect_lt(max_diff(st("Downweighted obs")$statistic, mean(w < 0.5)), tol[[nm]])
    expect_lt(max_diff(st("Minimum robustness weight")$statistic, min(w)), tol[[nm]])
    expect_lt(max_diff(st("Shapiro-Wilk")$statistic, sw$statistic), tol[[nm]])
    expect_lt(max_diff(st("Shapiro-Wilk")$p.value, sw$p.value), tol[[nm]])
  }
  expect_identical(min(fit_of("lr")$model$rweights), 0)

  for (nm in c("om", "mb", "mb_cl")) skip_engine(nm)
  for (nm in c("om", "mb", "mb_cl")) {
    fc <- fit_check(fit_of(nm), plot = FALSE)
    h <- hand_of(nm)
    expect_lt(max_diff(fc$tests$statistic[fc$tests$check == "AIC"], stats::AIC(h)), 1e-6)
    expect_lt(max_diff(fc$tests$statistic[fc$tests$check == "log-likelihood"],
                       as.numeric(stats::logLik(h))), 1e-6)
    expect_null(fc$raw)
    expect_true(is.character(fc$note) && nzchar(fc$note))
  }
})

test_that("five parts: clm and mblogit fits, omnibus, per-part row and contrast", {
  skip_engine("o")
  skip_engine("mb")
  set.seed(5)
  n5 <- 400
  parts5 <- paste0("P", 1:5)
  coords5 <- coord_names(parts5)
  raw5 <- cbind(P1 = stats::rlnorm(n5, log(60), 0.5), P2 = stats::rlnorm(n5, log(250), 0.4),
                P3 = stats::rlnorm(n5, log(500), 0.3), P4 = stats::rlnorm(n5, log(150), 0.5),
                P5 = stats::rlnorm(n5, log(480), 0.2))
  raw5 <- raw5 / rowSums(raw5) * closure
  lp5 <- 0.8 * log(raw5[, "P1"] / raw5[, "P3"]) - 0.5 * log(raw5[, "P2"] / raw5[, "P5"]) +
    0.3 * log(raw5[, "P4"] / raw5[, "P1"])
  lat5 <- lp5 + stats::rlogis(n5)
  df5 <- as.data.frame(raw5)
  df5$y_ord <- cut(lat5, stats::quantile(lat5, c(0, 1 / 3, 2 / 3, 1)),
                   include.lowest = TRUE, labels = c("low", "mid", "high"),
                   ordered_result = TRUE)
  pr5 <- cbind(1, exp(0.6 * log(raw5[, "P1"] / raw5[, "P3"])),
               exp(-0.5 * log(raw5[, "P2"] / raw5[, "P5"]) + 0.3 * log(raw5[, "P4"] / raw5[, "P1"])))
  pr5 <- pr5 / rowSums(pr5)
  df5$y_cat <- factor(c("a", "b", "c")[apply(pr5, 1, function(p) sample(3, 1, prob = p))],
                      levels = c("a", "b", "c"))
  cd5 <- declare_composition(df5, parts5, closure = closure)
  d5 <- cbind(df5, piv(raw5, 1))
  g5 <- exp(colMeans(log(raw5)))
  ref5 <- g5 / sum(g5) * closure
  new5 <- ref5
  new5[["P3"]] <- new5[["P3"]] - 15
  new5[["P1"]] <- new5[["P1"]] + 15
  z5 <- piv(rbind(ref5, new5), 1, parts5)
  dz5 <- stats::setNames(z5[2, ] - z5[1, ], coords5)
  rhs <- paste(coords5, collapse = " + ")

  c5 <- comp_iv(cd5, "y_ord", model = "clm")
  m5 <- quiet(comp_iv(cd5, "y_cat", model = "mblogit"))
  h_c5 <- hand_fit(ordinal::clm, stats::as.formula(paste("y_ord ~", rhs)), d5,
                   call_extras(c5$model))
  h_m5 <- quiet(hand_fit(mclogit::mblogit, stats::as.formula(paste("y_cat ~", rhs)), d5,
                         call_extras(m5$model)))
  for (pair in list(list(c5, h_c5, 4), list(m5, h_m5, 8))) {
    f <- pair[[1]]
    h <- pair[[2]]
    b <- coefs(f$model)
    expect_lt(max_diff(b, coefs(h)[names(b)]), 1e-8)
    expect_lt(max_diff(vc(f$model), vc(h, names(b))), 1e-8)
    nm_ilr <- ilr_coefs(names(coefs(h)))
    expect_equal(f$omnibus$df, pair[[3]])
    expect_length(nm_ilr, pair[[3]])
    expect_lt(max_diff(f$omnibus$statistic, wald(coefs(h), vc(h), nm_ilr)), 1e-8)
  }
  expect_lt(max_diff(summary(c5)$per_part$estimate[1], exp(coefs(h_c5)[["ilr1"]])), 1e-8)
  pr <- predict_realloc(c5, delta = 15, from = "P3", to = "P1")
  expect_lt(max_diff(pr$estimate, exp(sum(dz5 * coefs(h_c5)[coords5]))), 1e-8)
  expect_lt(max_diff(pr$std.error / pr$estimate,
                     sqrt(as.numeric(t(dz5) %*% vc(h_c5, coords5) %*% dz5))), 1e-8)
  pp <- summary(m5)$per_part
  pr <- predict_realloc(m5, delta = 15, from = "P3", to = "P1")
  for (lev in c("b", "c")) {
    row1 <- pp[pp$parameter == "P1 vs remaining" & pp$response_level == lev, ]
    expect_lt(max_diff(row1$estimate, exp(coefs(h_m5)[[paste0(lev, "~ilr1")]])), 1e-8)
    nmc <- paste0(lev, "~", coords5)
    r <- pr[pr$response_level == lev, ]
    expect_lt(max_diff(r$estimate, exp(sum(dz5 * coefs(h_m5)[nmc]))), 1e-8)
    expect_lt(max_diff(r$std.error / r$estimate,
                       sqrt(as.numeric(t(dz5) %*% vc(h_m5, nmc) %*% dz5))), 1e-8)
  }
})

test_that("each engine refuses an outcome or a grouping it cannot fit", {
  skip_on_cran()
  skip_if_not_installed("ordinal")
  skip_if_not_installed("mclogit")
  skip_if_not_installed("robustbase")
  skip_if_not_installed("robustlmm")
  expect_error(comp_iv(cd, "grade_unordered", model = "clm"), "ordered factor")
  expect_error(comp_iv(cd, "rz", model = "clm"), "ordered factor")
  expect_error(comp_iv(cd, "grade", model = "mblogit"), "ordered factor")
  expect_error(comp_iv(cd, "rz", model = "lmrob", family = "gaussian"), "does not apply")
  expect_error(comp_iv(cd, "rz", model = "rlmer"), "grouping term")
  expect_error(comp_iv(cd, "grade", model = "clmm"), "grouping term")
  expect_error(comp_iv(cd, "grade_unordered", model = "mblogit",
                       random = c("school", "region")), "one grouping term")
})
