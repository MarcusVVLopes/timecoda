# Numerical benchmark of non-linear composition terms in comp_iv() formulas
# (ilr()^2, I(ilr()^2), poly(ilr(), 2), splines::ns(ilr(), 3), the full
# quadratic ilr()^2 + I(ilr()^2), a formula without the main effect): a marker
# formula equals its hand-typed twin; every fit equals the model fitted by hand
# on robCompositions::pivotCoord() coordinates at every pivot; the omnibus is
# the joint Wald test over every composition term (car::linearHypothesis());
# the full quadratic is invariant to the part order and squares alone are not;
# predict_realloc() equals the hand contrast; the per-part rows are the slopes
# at the reference composition (the analytic derivative of the hand fit); and a
# formula without the main effect warns once. Reference
# calls are namespace-qualified so each number is visibly computed outside
# timecoda.

skip_on_cran()

closure <- 1440
parts <- c("MVPA", "LPA", "SB", "SPT")
parts3 <- c("A", "SB", "SPT")
deltas <- c(10, 30)
tu <- as.data.frame(timeuse)
tu$A <- tu$MVPA + tu$LPA
cd <- declare_composition(tu, parts, closure = closure)
cd3 <- declare_composition(tu, parts3, closure = closure)

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
spread <- function(m) max(apply(as.matrix(m), 1, function(r) max(r) - min(r)))
coord_names <- function(p) paste0("ilr", seq_len(length(p) - 1))
coefs <- function(m) if (inherits(m, "merMod")) lme4::fixef(m) else stats::coef(m)
vc <- function(m) {
  b <- coefs(m)
  as.matrix(stats::vcov(m))[names(b), names(b)]
}
wald <- function(m, nm) {
  b <- coefs(m)[nm]
  as.numeric(t(b) %*% solve(vc(m)[nm, nm, drop = FALSE]) %*% b)
}
tol_fit <- function(engine) if (engine == "lmer") 1e-6 else 1e-10
hand_frame <- function(k = 1, p = parts) {
  z <- robCompositions::pivotCoord(as.data.frame(tu[p]), pivotvar = k)
  colnames(z) <- coord_names(p)
  cbind(tu, z)
}
pivot_row <- function(x, p, k = 1) {
  z <- robCompositions::pivotCoord(as.data.frame(t(x[p])), pivotvar = k)
  stats::setNames(as.numeric(as.matrix(z)), coord_names(p))
}
hand_ref <- function(p) {
  x <- as.matrix(tu[p])
  g <- exp(colMeans(log(x / rowSums(x))))
  g / sum(g) * closure
}
# Coordinates of the reference and of `deltas` minutes moved from SB to MVPA.
coord_rows <- function() {
  ref <- hand_ref(parts)
  comps <- rbind(ref, t(vapply(deltas, function(d) {
    new <- ref
    new[["SB"]] <- new[["SB"]] - d
    new[["MVPA"]] <- new[["MVPA"]] + d
    new
  }, numeric(length(parts)))))
  as.data.frame(t(apply(comps, 1, pivot_row, p = parts)))
}

# S squares, P poly(), X products without squares, Q full quadratic, ns
# splines, Sx squares crossed with sex, noL squares without the main effect.
lin <- function(cn) paste(cn, collapse = " + ")
sq <- function(cn) paste0("I(", cn, "^2)", collapse = " + ")
hand_rhs <- function(form, cn) {
  switch(form,
    L = lin(cn),
    S = paste(lin(cn), "+", sq(cn)),
    P = paste0("poly(", cn, ", 2)", collapse = " + "),
    X = paste0("(", lin(cn), ")^2"),
    Q = paste0("(", lin(cn), ")^2 + ", sq(cn)),
    ns = paste0("splines::ns(", cn, ", 3)", collapse = " + "),
    Sx = paste0("(", lin(cn), " + ", sq(cn), ") * sex"),
    noL = sq(cn)
  )
}
marker_rhs <- c(L = "ilr()", S = "ilr() + I(ilr()^2)", P = "poly(ilr(), 2)",
                X = "ilr()^2", Q = "ilr()^2 + I(ilr()^2)",
                ns = "splines::ns(ilr(), 3)", Sx = "(ilr() + I(ilr()^2)) * sex")
typed_rhs <- function(form, cn) {
  switch(form,
    S = paste("ilr() +", sq(cn)),
    Q = paste0("ilr() + (", lin(cn), ")^2 + ", sq(cn))
  )
}
lhs_of <- c(lm = "zbmi", glm = "obesity", cox = "survival::Surv(time, event)",
            lmer = "zbmi")
with_rand <- function(rhs, engine) {
  if (engine == "lmer") paste(rhs, "+ (1 | school)") else rhs
}
pkg_fit <- function(rhs, engine, comp) {
  f <- stats::as.formula(paste(lhs_of[[engine]], "~", with_rand(rhs, engine)))
  switch(engine,
    lm = comp_iv(comp, f, model = "lm"),
    glm = comp_iv(comp, f, model = "glm", family = "binomial"),
    cox = comp_iv(comp, f, model = "cox"),
    lmer = comp_iv(comp, f, model = "lmer")
  )
}
hand_cache <- new.env()
hand_fit <- function(form, engine, k = 1, p = parts) {
  key <- paste(form, engine, k, length(p))
  if (!is.null(hand_cache[[key]])) return(hand_cache[[key]])
  data <- hand_frame(k, p)
  f <- stats::as.formula(paste(lhs_of[[engine]], "~",
                               with_rand(hand_rhs(form, coord_names(p)), engine)))
  hand_cache[[key]] <- switch(engine,
    lm = stats::lm(f, data = data),
    glm = stats::glm(f, data = data, family = stats::binomial()),
    cox = survival::coxph(f, data = data),
    lmer = lme4::lmer(f, data = data)
  )
}

sp <- function(form, engine, how = "marker", p = parts) {
  list(form = form, engine = engine, how = how, p = p)
}
specs <- list(
  m_L = sp("L", "lm"), m_S = sp("S", "lm"), m_P = sp("P", "lm"),
  m_X = sp("X", "lm"), m_Q = sp("Q", "lm"), m_ns = sp("ns", "lm"),
  m_Sx = sp("Sx", "lm"), g_S = sp("S", "glm"), g_Q = sp("Q", "glm"),
  c_S = sp("S", "cox"), l_S = sp("S", "lmer"),
  m_S3 = sp("S", "lm", p = parts3), m_Q3 = sp("Q", "lm", p = parts3),
  h_S = sp("S", "lm", "typed"), h_Q = sp("Q", "lm", "typed"),
  hg_S = sp("S", "glm", "typed"), hg_Q = sp("Q", "glm", "typed"),
  h_S3 = sp("S", "lm", "typed", parts3), h_Q3 = sp("Q", "lm", "typed", parts3)
)
hand_of <- function(s, k = 1) hand_fit(s$form, s$engine, k, s$p)
fits <- lapply(specs, function(s) {
  rhs <- if (s$how == "marker") marker_rhs[[s$form]] else typed_rhs(s$form, coord_names(s$p))
  pkg_fit(rhs, s$engine, if (length(s$p) == 3) cd3 else cd)
})
rel_cols <- c("estimate", "std.error", "df", "conf.low", "conf.high",
              "statistic", "p.value")

test_that("a marker formula and its hand-typed twin fit the same model", {
  skip_on_cran()
  same <- function(a, b) {
    nm <- sort(names(coefs(a$model)))
    if (!identical(nm, sort(names(coefs(b$model))))) return(Inf)
    pa <- summary(a)$per_part
    pb <- summary(b)$per_part
    cols <- names(pa)[vapply(pa, is.numeric, logical(1))]
    ra <- predict_realloc(a, delta = deltas, from = "SB", to = a$parts[1])
    rb <- predict_realloc(b, delta = deltas, from = "SB", to = b$parts[1])
    num <- names(ra)[vapply(ra, is.numeric, logical(1))]
    max(max_diff(coefs(a$model)[nm], coefs(b$model)[nm]),
        max_diff(vc(a$model)[nm, nm], vc(b$model)[nm, nm]),
        max_diff(a$omnibus[c("statistic", "df", "p.value")],
                 b$omnibus[c("statistic", "df", "p.value")]),
        if (identical(pa$parameter, pb$parameter)) max_diff(pa[cols], pb[cols]) else Inf,
        max_diff(ra[num], rb[num]))
  }
  twins <- list(c("m_S", "h_S"), c("m_Q", "h_Q"), c("g_S", "hg_S"),
                c("g_Q", "hg_Q"), c("m_S3", "h_S3"), c("m_Q3", "h_Q3"))
  expect_identical(max(vapply(twins, function(tw) same(fits[[tw[1]]], fits[[tw[2]]]),
                              numeric(1))), 0)
})

test_that("non-linear fits equal the models fitted by hand at every pivot", {
  skip_on_cran()
  # Coefficients and covariance of every rotated fit, scaled by the tolerance
  # so one expectation covers lmer; poly() and ns() bases on fitted values,
  # logLik and coefficients.
  worst <- max(vapply(c("m_S", "m_Q", "m_X", "m_Sx", "g_S", "g_Q", "c_S", "l_S"),
                      function(nm) {
    s <- specs[[nm]]
    max(vapply(seq_along(s$p), function(k) {
      m <- fits[[nm]]$models[[s$p[k]]]
      h <- hand_of(s, k)
      b <- coefs(m)
      if (!identical(sort(names(b)), sort(names(coefs(h))))) return(Inf)
      max(max_diff(b, coefs(h)[names(b)]), max_diff(vc(m), vc(h)[names(b), names(b)]))
    }, numeric(1))) / tol_fit(s$engine)
  }, numeric(1)))
  bases <- max(vapply(c("m_P", "m_ns"), function(nm) {
    s <- specs[[nm]]
    max(vapply(seq_along(s$p), function(k) {
      m <- fits[[nm]]$models[[s$p[k]]]
      h <- hand_of(s, k)
      if (!identical(names(coefs(m)), names(coefs(h)))) return(Inf)
      max(max_diff(stats::fitted(m), stats::fitted(h)), max_diff(coefs(m), coefs(h)),
          abs(as.numeric(stats::logLik(m)) - as.numeric(stats::logLik(h))))
    }, numeric(1))) / 1e-10
  }, numeric(1)))
  expect_lte(max(worst, bases), 1)
})

# Every non-intercept coefficient whose term involves the coordinates only
# (the hand formulas' one covariate is sex), and the df each form implies.
comp_only <- function(b) {
  n <- setdiff(names(b), "(Intercept)")
  n[!grepl("sex", n, fixed = TRUE)]
}
want_df <- function(form, D) {
  switch(form,
    S = 2 * (D - 1), P = 2 * (D - 1), Q = (D - 1) * (D + 2) / 2,
    X = (D - 1) + (D - 1) * (D - 2) / 2, ns = 3 * (D - 1), Sx = 2 * (D - 1)
  )
}
c_keys <- c("m_S", "m_Q", "m_P", "m_X", "m_ns", "m_Sx", "g_S", "g_Q", "c_S",
            "l_S", "m_S3", "m_Q3")

test_that("the omnibus is the joint Wald test over every composition term", {
  skip_on_cran()
  out <- lapply(c_keys, function(nm) {
    s <- specs[[nm]]
    om <- fits[[nm]]$omnibus
    h <- hand_of(s)
    nms <- comp_only(coefs(h))
    wd <- want_df(s$form, length(s$p))
    st <- wald(h, nms)
    list(df = c(om$df - wd, length(nms) - wd),
         stat = max(abs(om$statistic - wald(fits[[nm]]$model, nms)) / 1e-10,
                    max(abs(om$statistic - st),
                        abs(om$p.value - stats::pchisq(st, wd, lower.tail = FALSE))) /
                      tol_fit(s$engine)))
  })
  # Squares and poly() span one column space: the same hypothesis.
  same_space <- max_diff(fits$m_P$omnibus[c("statistic", "df", "p.value")],
                         fits$m_S$omnibus[c("statistic", "df", "p.value")]) / 1e-10
  expect_identical(max(abs(unlist(lapply(out, `[[`, "df")))), 0)
  expect_lte(max(vapply(out, `[[`, numeric(1), "stat"), same_space), 1)
})

test_that("the omnibus equals car::linearHypothesis()", {
  skip_on_cran()
  skip_if_not_installed("car")
  car_chisq <- function(m, nm) {
    # A 0/1 selection matrix: car cannot parse names such as I(ilr1^2).
    L <- matrix(0, length(nm), length(coefs(m)), dimnames = list(NULL, names(coefs(m))))
    L[cbind(seq_along(nm), match(nm, names(coefs(m))))] <- 1
    lh <- car::linearHypothesis(m, hypothesis.matrix = L, test = "Chisq")
    c(lh$Chisq[2], lh[["Pr(>Chisq)"]][2], lh$Df[2])
  }
  worst <- max(vapply(c_keys, function(nm) {
    s <- specs[[nm]]
    om <- fits[[nm]]$omnibus
    h <- hand_of(s)
    max_diff(c(om$statistic, om$p.value, om$df), car_chisq(h, comp_only(coefs(h)))) /
      max(tol_fit(s$engine), 1e-8)
  }, numeric(1)))
  expect_lte(worst, 1)
})

test_that("the full quadratic does not depend on the part order; squares alone do", {
  skip_on_cran()
  orders <- list(list(c("SB", "SPT", "MVPA", "LPA"), c("SPT", "MVPA", "LPA", "SB")),
                 list(c("SB", "SPT", "A"), c("SPT", "A", "SB")))
  quantities <- function(fit, to) {
    rel <- predict_realloc(fit, delta = deltas, from = "SB", to = to)
    ab <- predict_realloc(fit, delta = deltas, from = "SB", to = to, type = "absolute")
    list(contrast = rel$estimate,
         rest = c(rel$std.error, rel$conf.low, rel$conf.high, ab$estimate,
                  fit$omnibus$statistic),
         logLik = as.numeric(stats::logLik(fit$model)))
  }
  spreads <- function(nm) {
    s <- specs[[nm]]
    three <- length(s$p) == 3
    to <- if (three) "A" else "MVPA"
    per <- c(list(quantities(fits[[nm]], to)),
             lapply(orders[[if (three) 2 else 1]], function(o) {
               quantities(pkg_fit(marker_rhs[[s$form]], s$engine,
                                  declare_composition(tu, o, closure = closure)), to)
             }))
    vapply(names(per[[1]]), function(q) spread(do.call(cbind, lapply(per, `[[`, q))),
           numeric(1))
  }
  invariant <- max(vapply(c("m_L", "m_Q", "g_Q", "m_Q3"), function(nm) {
    max(spreads(nm)) / if (specs[[nm]]$engine == "glm") 1e-8 else 1e-10
  }, numeric(1)))
  # Squares-only and poly() fits are a different model in each basis: the
  # reallocation contrast and the logLik move with the part order.
  order_dependent <- vapply(c("m_S", "m_P", "g_S"), function(nm) {
    sp <- spreads(nm)
    min(sp[["contrast"]], sp[["logLik"]])
  }, numeric(1))
  expect_lte(invariant, 1)
  expect_gt(min(order_dependent), 1e-6)
})

test_that("reallocations of non-linear fits equal the hand contrast", {
  skip_on_cran()
  notes <- 0
  count_notes <- function(expr) {
    withCallingHandlers(expr, message = function(m) {
      if (grepl("involvement in interactions", conditionMessage(m), fixed = TRUE)) {
        notes <<- notes + 1
      }
    })
  }
  nd <- coord_rows()
  # Design rows from the model's own terms, contrasted with the reference row.
  hand_rel <- function(fit, engine) {
    m <- fit$model
    tt <- stats::delete.response(stats::terms(m))
    X <- stats::model.matrix(tt, stats::model.frame(tt, nd))[, names(coefs(m)), drop = FALSE]
    L <- X[-1, , drop = FALSE] - X[rep(1, length(deltas)), , drop = FALSE]
    est <- as.numeric(L %*% coefs(m))
    se <- sqrt(rowSums((L %*% vc(m)) * L))
    df <- if (engine == "lm") stats::df.residual(m) else Inf
    crit <- stats::qt(0.975, df)
    z <- est / se
    out <- cbind(est, se, df, est - crit * se, est + crit * se, z,
                 2 * stats::pt(-abs(z), df))
    if (engine != "lm") {
      out[, c(1, 4, 5)] <- exp(out[, c(1, 4, 5)])
      out[, 2] <- exp(est) * se
    }
    out
  }
  relative <- max(vapply(c("m_S", "m_Q", "g_S", "g_Q", "c_S"), function(nm) {
    s <- specs[[nm]]
    pr <- count_notes(predict_realloc(fits[[nm]], delta = deltas, from = "SB", to = "MVPA"))
    utils::capture.output(count_notes(print(pr)))
    max_diff(pr[rel_cols], hand_rel(fits[[nm]], s$engine)) /
      if (s$engine == "glm") 1e-8 else 1e-10
  }, numeric(1)))
  absolute <- max(vapply(c("m_S", "g_S"), function(nm) {
    s <- specs[[nm]]
    ab <- count_notes(predict_realloc(fits[[nm]], delta = deltas, from = "SB",
                                      to = "MVPA", type = "absolute"))
    max_diff(ab$estimate, stats::predict(hand_of(s), newdata = nd, type = "response")) /
      if (s$engine == "glm") 1e-8 else 1e-10
  }, numeric(1)))
  expect_lte(max(relative, absolute), 1)
  # Coordinate products are not covariate interactions; squares crossed with
  # sex are. No emmeans interaction note on any of these calls.
  expect_identical(
    c(notes, fits$m_Q$has_ilr_interaction, fits$g_Q$has_ilr_interaction,
      fits$m_Sx$has_ilr_interaction),
    c(0, 0, 0, 1)
  )
})

test_that("per-part rows report each rotated fit's slope at the reference", {
  skip_on_cran()
  worst <- max(vapply(c("m_S", "m_Q"), function(nm) {
    pp <- summary(fits[[nm]])$per_part
    if (!identical(pp$parameter, paste(parts, "vs remaining"))) return(Inf)
    # The analytic derivative of the hand fit in ilr1 at the reference.
    max_diff(pp$estimate, vapply(seq_along(parts), function(k) {
      z <- pivot_row(hand_ref(parts), parts, k)
      b <- coefs(hand_of(specs[[nm]], k))
      cross <- paste0("ilr1:", names(z)[-1])
      b[["ilr1"]] + 2 * b[["I(ilr1^2)"]] * z[["ilr1"]] +
        sum(ifelse(cross %in% names(b), b[cross], 0) * z[-1])
    }, numeric(1)))
  }, numeric(1)))
  expect_lte(worst, 1e-6)
})

test_that("a formula without the composition's main effect warns once", {
  skip_on_cran()
  with_warnings <- function(expr) {
    w <- character(0)
    fit <- withCallingHandlers(expr, warning = function(cnd) {
      w <<- c(w, conditionMessage(cnd))
      invokeRestart("muffleWarning")
    })
    list(fit = fit, warnings = w)
  }
  no_main <- with_warnings(comp_iv(cd, zbmi ~ I(ilr()^2), model = "lm"))
  others <- vapply(list(zbmi ~ poly(ilr(), 2), zbmi ~ ilr() + I(ilr()^2),
                        zbmi ~ ilr(), zbmi ~ sex), function(f) {
    length(with_warnings(comp_iv(cd, f, model = "lm"))$warnings)
  }, integer(1))
  expect_identical(c(length(no_main$warnings), sum(grepl("main effect", no_main$warnings)),
                     others), c(1L, 1L, 0L, 0L, 0L, 0L))
  m <- no_main$fit$model
  h <- hand_fit("noL", "lm")
  diff <- if (identical(names(coef(m)), names(coef(h)))) {
    max(max_diff(coef(m), coef(h)), max_diff(vcov(m), vcov(h)))
  } else {
    Inf
  }
  expect_lte(diff, 1e-10)
})
