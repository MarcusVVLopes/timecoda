# Numerical benchmark of moderation, group-wise reallocation and the formula
# interface of comp_iv(): a formula fit equals its discrete-argument twin and
# moderator = equals the explicit `ilr() * moderator` formula; every fit equals
# the model fitted by hand on robCompositions::pivotCoord() coordinates; the
# joint and per-part moderation Wald tests equal b' V^-1 b on the hand fit and
# car::linearHypothesis(); predict_realloc(by =) equals the hand contrast at
# each level (Kenward-Roger covariance and df from pbkrtest for lmer); and the
# results are invariant to the part order. Reference calls are
# namespace-qualified so each number is visibly computed outside timecoda.

skip_on_cran()

closure <- 1440
parts <- c("MVPA", "LPA", "SB", "SPT")
parts3 <- c("A", "SB", "SPT")
deltas <- c(10, 30)
tu <- as.data.frame(timeuse)
tu$A <- tu$MVPA + tu$LPA
cd <- declare_composition(tu, parts, closure = closure)
cd3 <- declare_composition(tu, parts3, closure = closure)

# emmeans notes "involvement in interactions" for every reference grid of a
# moderated fit; only that note is muffled.
quiet_note <- function(expr) {
  withCallingHandlers(expr, message = function(m) {
    if (grepl("involvement in interactions", conditionMessage(m), fixed = TRUE)) {
      invokeRestart("muffleMessage")
    }
  })
}
# Max abs difference; Inf when the lengths or the NA patterns differ.
max_diff <- function(a, b) {
  a <- as.numeric(unlist(a))
  b <- as.numeric(unlist(b))
  if (length(a) != length(b) || !identical(is.na(a), is.na(b))) return(Inf)
  ok <- !is.na(a)
  if (!identical(a[ok & !is.finite(a)], b[ok & !is.finite(a)])) return(Inf)
  ok <- ok & is.finite(a)
  max(c(0, abs(a[ok] - b[ok])))
}
coord_names <- function(p) paste0("ilr", seq_len(length(p) - 1))
coefs <- function(m) if (inherits(m, "merMod")) lme4::fixef(m) else stats::coef(m)
vc <- function(m) {
  b <- coefs(m)
  as.matrix(stats::vcov(m))[names(b), names(b)]
}
hand_frame <- function(k, p) {
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
  g <- exp(colMeans(log(x / rowSums(x))))
  g / sum(g) * closure
}

specs <- list(
  lm_mod = list(engine = "lm", y = "zbmi", mod = "sex", p = parts),
  lm_mod_region = list(engine = "lm", y = "zbmi", mod = "region", p = parts),
  glm_mod = list(engine = "glm", y = "obesity", mod = "sex", p = parts),
  glm_mod_num = list(engine = "glm", y = "obesity", mod = "zbmi", p = parts),
  cox_mod = list(engine = "cox", y = NULL, mod = "sex", p = parts),
  lmer_mod = list(engine = "lmer", y = "zbmi", mod = "sex", p = parts),
  lm3_mod = list(engine = "lm", y = "zbmi", mod = "sex", p = parts3)
)
mod_fit <- function(s, comp) {
  switch(s$engine,
    lm = comp_iv(comp, s$y, model = "lm", moderator = s$mod),
    glm = comp_iv(comp, s$y, model = "glm", family = "binomial",
                  moderator = s$mod),
    cox = comp_iv(comp, time = "time", event = "event", model = "cox",
                  moderator = s$mod),
    lmer = comp_iv(comp, s$y, model = "lmer", random = "school",
                   moderator = s$mod)
  )
}
hand_fit <- function(s, k) {
  data <- hand_frame(k, s$p)
  rhs <- paste0("(", paste(coord_names(s$p), collapse = " + "), ") * ", s$mod)
  f <- function(lhs, bars = FALSE) {
    stats::as.formula(paste(lhs, "~", rhs, if (bars) "+ (1 | school)"))
  }
  switch(s$engine,
    lm = stats::lm(f(s$y), data = data),
    glm = stats::glm(f(s$y), data = data, family = stats::binomial()),
    cox = survival::coxph(f("survival::Surv(time, event)"), data = data),
    lmer = lme4::lmer(f(s$y, TRUE), data = data)
  )
}
fits <- lapply(specs, function(s) mod_fit(s, if (length(s$p) == 3) cd3 else cd))
hands <- lapply(specs, function(s) lapply(seq_along(s$p), function(k) hand_fit(s, k)))
tol_fit <- function(s) if (s$engine == "lmer") 1e-6 else 1e-10
int_names <- function(b, mod, coord = "ilr[0-9]+") {
  grep(paste0("^", coord, ":", mod), names(b), value = TRUE)
}
wald <- function(m, nm) {
  b <- coefs(m)[nm]
  as.numeric(t(b) %*% solve(vc(m)[nm, nm, drop = FALSE]) %*% b)
}
rel_cols <- c("estimate", "std.error", "df", "conf.low", "conf.high",
              "statistic", "p.value")

test_that("formula, discrete and moderator interfaces fit the same model", {
  skip_on_cran()
  same <- function(a, b, at = NULL, pp_cols = NULL) {
    if (!identical(names(coefs(a$model)), names(coefs(b$model)))) return(Inf)
    pa <- summary(a)$per_part
    pb <- summary(b)$per_part
    if (is.null(pp_cols)) {
      pp_cols <- intersect(names(pa)[vapply(pa, is.numeric, logical(1))], names(pb))
    }
    ra <- quiet_note(predict_realloc(a, delta = deltas, from = "SB",
                                     to = a$parts[1], at = at))
    rb <- quiet_note(predict_realloc(b, delta = deltas, from = "SB",
                                     to = b$parts[1], at = at))
    num <- names(ra)[vapply(ra, is.numeric, logical(1))]
    max(max_diff(coefs(a$model), coefs(b$model)), max_diff(vc(a$model), vc(b$model)),
        max_diff(a$omnibus[c("statistic", "df", "p.value")],
                 b$omnibus[c("statistic", "df", "p.value")]),
        if (identical(pa$parameter, pb$parameter)) max_diff(pa[pp_cols], pb[pp_cols]) else Inf,
        max_diff(unlist(a$ref_comp), unlist(b$ref_comp)),
        max_diff(ra[num], rb[num]))
  }
  pp <- c("estimate", "conf.low", "conf.high", "p.value")
  d_add <- comp_iv(cd, "zbmi", model = "lm", covariates = "sex")
  f_add <- comp_iv(cd, zbmi ~ ilr() + sex, model = "lm")
  worst <- max(
    same(f_add, d_add),
    same(comp_iv(cd, zbmi ~ sex, model = "lm"), d_add),
    same(fits$lm_mod, comp_iv(cd, zbmi ~ ilr() * sex, model = "lm"),
         at = list(sex = "female"), pp_cols = pp),
    same(comp_iv(cd, zbmi ~ ilr() + sex + (1 | school), model = "lmer"),
         comp_iv(cd, "zbmi", model = "lmer", random = "school", covariates = "sex")),
    same(comp_iv(cd, survival::Surv(time, event) ~ ilr() + sex, model = "cox"),
         comp_iv(cd, time = "time", event = "event", model = "cox", covariates = "sex")),
    same(comp_iv(cd3, zbmi ~ ilr() + sex, model = "lm"),
         comp_iv(cd3, "zbmi", model = "lm", covariates = "sex")),
    same(fits$lm3_mod, comp_iv(cd3, zbmi ~ ilr() * sex, model = "lm"),
         at = list(sex = "female"), pp_cols = pp),
    same(comp_iv(cd, zbmi ~ ilr() + region, model = "lm"),
         comp_iv(cd, "zbmi", model = "lm", covariates = "region")),
    same(fits$lm_mod_region, comp_iv(cd, zbmi ~ ilr() * region, model = "lm"),
         at = list(region = "A"), pp_cols = pp)
  )
  expect_identical(worst, 0)
})

test_that("every moderated fit equals the model fitted by hand", {
  skip_on_cran()
  # Coefficient names, coefficients and covariance of the canonical and every
  # rotated fit; scaled by the tolerance so one expectation covers lmer.
  worst <- max(vapply(names(specs), function(nm) {
    s <- specs[[nm]]
    d <- max(vapply(seq_along(s$p), function(k) {
      m <- fits[[nm]]$models[[s$p[k]]]
      h <- hands[[nm]][[k]]
      if (!identical(names(coefs(m)), names(coefs(h)))) return(Inf)
      max(max_diff(coefs(m), coefs(h)), max_diff(vc(m), vc(h)))
    }, numeric(1)))
    d / tol_fit(s)
  }, numeric(1)))
  f_poly <- comp_iv(cd, obesity ~ ilr() + poly(zbmi, 2), model = "glm",
                    family = "binomial")
  h_poly <- stats::glm(obesity ~ (ilr1 + ilr2 + ilr3) + poly(zbmi, 2),
                       data = hand_frame(1, parts), family = stats::binomial())
  poly_d <- if (identical(names(coef(f_poly$model)), names(coef(h_poly)))) {
    max(max_diff(coef(f_poly$model), coef(h_poly)),
        max_diff(vcov(f_poly$model), vcov(h_poly)))
  } else {
    Inf
  }
  expect_lte(max(worst, poly_d / 1e-10), 1)
})

test_that("the moderation tests are the joint Wald test on the hand fit", {
  skip_on_cran()
  # Overall: the ilr{k}:moderator block, (L - 1)(D - 1) df; per part: the
  # ilr1:moderator coefficient(s) of the fit with that part pivoted first,
  # (b / se)^2 for a single coefficient.
  out <- lapply(names(specs), function(nm) {
    s <- specs[[nm]]
    md <- fits[[nm]]$moderation
    h <- hands[[nm]][[1]]
    nms <- int_names(coefs(h), s$mod)
    st <- wald(h, nms)
    pp <- summary(fits[[nm]])$per_part
    ref <- vapply(seq_along(s$p), function(k) {
      hk <- hands[[nm]][[k]]
      nk <- int_names(coefs(hk), s$mod, coord = "ilr1")
      stk <- wald(hk, nk)
      z2 <- if (length(nk) == 1) (coefs(hk)[[nk]])^2 / vc(hk)[nk, nk] else stk
      c(stk, stats::pchisq(stk, length(nk), lower.tail = FALSE), z2, length(nk))
    }, numeric(4))
    n_lev <- if (s$mod == "region") 3 else 2
    list(
      df = c(md$df - length(nms), md$df - (n_lev - 1) * (length(s$p) - 1),
             pp$mod_df - ref[4, ], ref[4, ] - (n_lev - 1)),
      stat = max(abs(md$statistic - st),
                 abs(md$p.value - stats::pchisq(st, md$df, lower.tail = FALSE)),
                 max_diff(pp$mod_statistic, ref[1, ]),
                 max_diff(pp$mod_p.value, ref[2, ]),
                 max_diff(pp$mod_statistic, ref[3, ])) / tol_fit(s)
    )
  })
  expect_identical(max(abs(as.numeric(unlist(lapply(out, `[[`, "df"))))), 0)
  expect_lte(max(vapply(out, `[[`, numeric(1), "stat")), 1)
})

test_that("the moderation tests equal car::linearHypothesis()", {
  skip_on_cran()
  skip_if_not_installed("car")
  car_chisq <- function(m, nm) {
    lh <- car::linearHypothesis(m, hypothesis.matrix = nm, test = "Chisq")
    c(lh$Chisq[2], lh[["Pr(>Chisq)"]][2])
  }
  worst <- max(vapply(names(specs), function(nm) {
    s <- specs[[nm]]
    md <- fits[[nm]]$moderation
    h <- hands[[nm]][[1]]
    pp <- summary(fits[[nm]])$per_part
    per_part <- vapply(seq_along(s$p), function(k) {
      hk <- hands[[nm]][[k]]
      car_chisq(hk, int_names(coefs(hk), s$mod, coord = "ilr1"))
    }, numeric(2))
    max(max_diff(c(md$statistic, md$p.value), car_chisq(h, int_names(coefs(h), s$mod))),
        max_diff(rbind(pp$mod_statistic, pp$mod_p.value), per_part)) /
      max(tol_fit(s), 1e-8)
  }, numeric(1)))
  expect_lte(worst, 1)
})

test_that("predict_realloc(by =) equals the hand contrast at each level", {
  skip_on_cran()
  hand_level <- function(m, s, lev, kr) {
    b <- coefs(m)
    V <- vc(m)
    if (kr) {
      V <- as.matrix(pbkrtest::vcovAdj(m))
      dimnames(V) <- list(names(b), names(b))
    }
    ref <- hand_ref(s$p)
    zref <- pivot_row(ref, s$p)
    ratio <- s$engine %in% c("glm", "cox")
    t(vapply(deltas, function(d) {
      new <- ref
      new[["SB"]] <- new[["SB"]] - d
      new[["MVPA"]] <- new[["MVPA"]] + d
      dz <- pivot_row(new, s$p) - zref
      L <- stats::setNames(numeric(length(b)), names(b))
      L[names(dz)] <- dz
      int <- paste0(names(dz), ":", s$mod, lev)
      if (all(int %in% names(b))) L[int] <- dz
      est <- sum(L * b)
      se <- sqrt(as.numeric(t(L) %*% V %*% L))
      df <- if (kr) {
        pbkrtest::get_Lb_ddf(m, matrix(L, nrow = 1))
      } else if (s$engine == "lm") {
        stats::df.residual(m)
      } else {
        Inf
      }
      crit <- stats::qt(0.975, df)
      z <- est / se
      out <- c(est, se, df, est - crit * se, est + crit * se, z,
               2 * stats::pt(-abs(z), df))
      if (ratio) out[c(1, 4, 5)] <- exp(out[c(1, 4, 5)])
      if (ratio) out[2] <- exp(est) * se
      out
    }, numeric(7)))
  }
  d_fits <- c("lm_mod", "lm_mod_region", "glm_mod", "cox_mod")
  if (requireNamespace("pbkrtest", quietly = TRUE)) d_fits <- c(d_fits, "lmer_mod")
  res <- lapply(d_fits, function(nm) {
    s <- specs[[nm]]
    pr <- quiet_note(predict_realloc(fits[[nm]], delta = deltas, from = "SB",
                                     to = "MVPA"))
    levs <- levels(tu[[s$mod]])
    if (nrow(pr) != length(levs) * length(deltas)) return(c(Inf, Inf))
    d <- vapply(levs, function(lev) {
      got <- pr[pr[[s$mod]] == lev, rel_cols]
      at_l <- quiet_note(predict_realloc(fits[[nm]], delta = deltas, from = "SB",
                                         to = "MVPA",
                                         at = stats::setNames(list(lev), s$mod)))
      c(max_diff(got, at_l[rel_cols]),
        max_diff(got, hand_level(hands[[nm]][[1]], s, lev, s$engine == "lmer")) /
          tol_fit(s))
    }, numeric(2))
    apply(d, 1, max)
  })
  # Without the interaction the covariate cancels: one contrast at every level.
  cancel <- max(vapply(c("lm", "glm"), function(engine) {
    fit <- if (engine == "lm") {
      comp_iv(cd, "zbmi", model = "lm", covariates = "sex")
    } else {
      comp_iv(cd, "obesity", model = "glm", family = "binomial", covariates = "sex")
    }
    pr0 <- predict_realloc(fit, delta = deltas, from = "SB", to = "MVPA")
    prb <- predict_realloc(fit, delta = deltas, from = "SB", to = "MVPA", by = "sex")
    max(vapply(levels(tu$sex), function(lev) {
      max_diff(prb[prb$sex == lev, rel_cols], pr0[rel_cols])
    }, numeric(1)))
  }, numeric(1)))
  expect_identical(max(vapply(res, `[`, numeric(1), 1)), 0)
  expect_lte(max(vapply(res, `[`, numeric(1), 2), cancel / 1e-12), 1)
})

test_that("moderation and by-level reallocations do not depend on the part order", {
  skip_on_cran()
  orders <- list(rev(parts), c(parts[-1], parts[1]), c("SB", "MVPA", "SPT", "LPA"))
  canon <- function(pr, mod) pr[order(pr[[mod]], pr$delta), rel_cols]
  worst <- max(vapply(c("lm_mod", "glm_mod"), function(nm) {
    s <- specs[[nm]]
    base <- canon(quiet_note(predict_realloc(fits[[nm]], delta = deltas, from = "SB",
                                             to = "MVPA")), s$mod)
    pb <- summary(fits[[nm]])$per_part
    d <- max(vapply(orders, function(o) {
      rf <- mod_fit(s, declare_composition(tu, o, closure = closure))
      pa <- summary(rf)$per_part
      m <- match(pb$parameter, pa$parameter)
      max(max_diff(rf$moderation[c("statistic", "df", "p.value")],
                   fits[[nm]]$moderation[c("statistic", "df", "p.value")]),
          max_diff(canon(quiet_note(predict_realloc(rf, delta = deltas, from = "SB",
                                                    to = "MVPA")), s$mod), base),
          max_diff(pa[m, c("mod_statistic", "mod_p.value")],
                   pb[c("mod_statistic", "mod_p.value")]))
    }, numeric(1)))
    d / if (s$engine == "glm") 1e-3 else 1e-8
  }, numeric(1)))
  expect_lte(worst, 1)
})
