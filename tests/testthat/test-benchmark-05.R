# Numerical benchmark of predict_realloc() for model = "lm" and "glm"
# (binomial and Poisson), one-to-one and one-to-all, relative and absolute,
# against a hand computation from coef()/vcov() of the fitted model on
# robCompositions::pivotCoord() coordinates, and against
# codaredistlm::predict_delta_comps(), which uses the same pivot basis and the
# same geometric-mean reference. Also checks the reference-versus-reference
# contrast, invariance to the part order (ILR basis), and that
# summary.realloc(), pub_ready() and plot_realloc() carry exactly the
# predict_realloc() numbers. Reference calls are namespace-qualified so each
# number is visibly computed outside timecoda.

skip_if_not_installed("robCompositions")
skip_on_cran()

parts <- c("MVPA", "LPA", "SB", "SPT")
closure <- 1440
deltas <- c(15, 30)
ilr_names <- paste0("ilr", 1:3)
tu <- as.data.frame(timeuse)
tu$count <- round(exp(pmin(pmax(tu$zbmi, -3), 3)))
cd <- declare_composition(tu, parts, closure = closure)

fits <- list(
  lm0 = comp_iv(cd, "zbmi", model = "lm"),
  lmc = comp_iv(cd, "zbmi", model = "lm", covariates = "sex"),
  glmb = comp_iv(cd, "obesity", model = "glm", family = "binomial"),
  glmp = comp_iv(cd, "count", model = "glm", family = "poisson")
)

# Max-abs-difference expectation with matching NA patterns.
expect_close <- function(a, b, tol) {
  a <- as.numeric(unlist(a))
  b <- as.numeric(unlist(b))
  expect_identical(is.na(a), is.na(b))
  ok <- !is.na(a) & is.finite(a)
  expect_identical(a[!ok & !is.na(a)], b[!ok & !is.na(a)])
  expect_lte(max(c(0, abs(a[ok] - b[ok]))), tol)
}

comp_mat <- as.matrix(tu[parts])
comp_mat <- comp_mat / rowSums(comp_mat) * closure
gm <- exp(colMeans(log(comp_mat)))
ref_hand <- gm / sum(gm) * closure

ilr_of <- function(comps) {
  z <- as.matrix(robCompositions::pivotCoord(as.data.frame(comps),
                                             pivotvar = 1))
  colnames(z) <- ilr_names
  z
}

moves_1to1 <- expand.grid(from = parts, to = parts, stringsAsFactors = FALSE)
moves_1to1 <- moves_1to1[moves_1to1$from != moves_1to1$to, ]
moves_1toall <- rbind(
  data.frame(from = "(rest)", to = parts, focal = parts, role = "gain"),
  data.frame(from = parts, to = "(rest)", focal = parts, role = "give")
)

new_comp <- function(ref, move, delta, scheme) {
  x <- ref
  if (scheme == "one-to-one") {
    x[move$from] <- x[move$from] - delta
    x[move$to] <- x[move$to] + delta
  } else {
    f <- move$focal
    rest <- setdiff(names(ref), f)
    w <- ref[rest] / sum(ref[rest])
    s <- if (move$role == "gain") 1 else -1
    x[f] <- x[f] + s * delta
    x[rest] <- x[rest] - s * delta * w
  }
  x
}

# Relative: d'b (lm) or exp(d'b) (glm), SE sqrt(d'Vd), t CI on the residual
# df (lm) or normal CI (glm); emmeans reports a glm's std.error on the
# response scale, ratio * SE_link. Absolute: x'b back-transformed, with the
# covariates in `cov_x` (emmeans averages a factor over its levels with equal
# weights).
hand_realloc <- function(fit, scheme, type, cov_x = NULL) {
  moves <- if (scheme == "one-to-one") moves_1to1 else moves_1toall
  grid <- expand.grid(i = seq_len(nrow(moves)), delta = deltas)
  comps <- t(vapply(seq_len(nrow(grid)), function(r) {
    new_comp(ref_hand, moves[grid$i[r], ], grid$delta[r], scheme)
  }, numeric(length(ref_hand))))
  colnames(comps) <- names(ref_hand)
  valid <- rowSums(comps <= 0) == 0
  m <- fit$model
  b <- stats::coef(m)
  V <- stats::vcov(m)
  is_glm <- inherits(m, "glm")
  df <- if (is_glm) Inf else stats::df.residual(m)
  crit <- if (is_glm) stats::qnorm(0.975) else stats::qt(0.975, df)
  out <- data.frame(from = moves$from[grid$i], to = moves$to[grid$i],
                    delta = grid$delta, stringsAsFactors = FALSE)
  if (type == "relative") {
    est <- se <- rep(NA_real_, nrow(grid))
    d <- sweep(ilr_of(comps[valid, , drop = FALSE]), 2,
               ilr_of(rbind(ref_hand))[1, ])
    D <- matrix(0, nrow(d), length(b), dimnames = list(NULL, names(b)))
    D[, ilr_names] <- d
    est[valid] <- as.numeric(D %*% b)
    se[valid] <- sqrt(rowSums((D %*% V) * D))
    stat <- est / se
    tr <- if (is_glm) exp else identity
    out$estimate <- tr(est)
    out$std.error <- if (is_glm) exp(est) * se else se
    out$conf.low <- tr(est - crit * se)
    out$conf.high <- tr(est + crit * se)
    out$statistic <- stat
    out$p.value <- if (is_glm) 2 * stats::pnorm(-abs(stat)) else
      2 * stats::pt(-abs(stat), df)
    out$se_link <- se
  } else {
    out <- rbind(data.frame(from = NA_character_, to = NA_character_,
                            delta = 0), out)
    allc <- rbind(ref_hand, comps)
    valid <- c(TRUE, valid)
    X <- matrix(0, nrow(allc), length(b), dimnames = list(NULL, names(b)))
    X[, "(Intercept)"] <- 1
    X[valid, ilr_names] <- ilr_of(allc[valid, , drop = FALSE])
    for (nm in names(cov_x)) X[, nm] <- cov_x[[nm]]
    eta <- as.numeric(X %*% b)
    se <- sqrt(rowSums((X %*% V) * X))
    eta[!valid] <- NA
    se[!valid] <- NA
    fam <- if (is_glm) stats::family(m) else stats::gaussian()
    out$estimate <- fam$linkinv(eta)
    out$std.error <- fam$mu.eta(eta) * se
    out$conf.low <- fam$linkinv(eta - crit * se)
    out$conf.high <- fam$linkinv(eta + crit * se)
  }
  out$df <- ifelse(valid, df, NA_real_)
  out
}

pkg_realloc <- function(fit, scheme, type, ...) {
  r <- if (scheme == "one-to-one") {
    predict_realloc(fit, delta = deltas, type = type, on_invalid = "na", ...)
  } else {
    dplyr::bind_rows(
      predict_realloc(fit, delta = deltas, scheme = "one-to-all", to = parts,
                      type = type, on_invalid = "na", ...),
      predict_realloc(fit, delta = deltas, scheme = "one-to-all",
                      from = parts, type = type, on_invalid = "na", ...)
    )
  }
  as.data.frame(r)
}
key_of <- function(df) paste(df$from, df$to, df$delta, sep = "|")
align <- function(pkg, hand) {
  ref_rows <- which(is.na(pkg$from))
  if (length(ref_rows) > 1) pkg <- pkg[-ref_rows[-1], ]
  out <- pkg[match(key_of(hand), key_of(pkg)), ]
  expect_false(anyNA(out$delta))
  out
}

test_that("reallocations match the hand contrast d'b and SE sqrt(d'Vd)", {
  rel_cols <- c("estimate", "std.error", "conf.low", "conf.high",
                "statistic", "p.value", "df")
  abs_cols <- c("estimate", "std.error", "conf.low", "conf.high", "df")
  for (nm in names(fits)) {
    fit <- fits[[nm]]
    cov_x <- if (nm == "lmc") list(sexmale = 0.5) else NULL
    for (scheme in c("one-to-one", "one-to-all")) {
      for (type in c("relative", "absolute")) {
        hand <- hand_realloc(fit, scheme, type, cov_x = cov_x)
        pkg <- align(pkg_realloc(fit, scheme, type), hand)
        cols <- if (type == "relative") rel_cols else abs_cols
        for (col in cols) expect_close(pkg[[col]], hand[[col]], 1e-6)
        if (type == "relative" && inherits(fit$model, "glm")) {
          expect_close(pkg$std.error / pkg$estimate, hand$se_link, 1e-6)
        }
      }
    }
  }
})

test_that("covariates cancel in the relative contrast but not the absolute", {
  hand <- hand_realloc(fits$lmc, "one-to-one", "relative")
  for (at in list(NULL, list(sex = "female"), list(sex = "male"))) {
    pkg <- align(pkg_realloc(fits$lmc, "one-to-one", "relative", at = at), hand)
    expect_close(pkg[c("estimate", "std.error", "conf.low", "conf.high")],
                 hand[c("estimate", "std.error", "conf.low", "conf.high")],
                 1e-6)
  }
  for (lev in c("female", "male")) {
    hand_a <- hand_realloc(fits$lmc, "one-to-one", "absolute",
                           cov_x = list(sexmale = as.numeric(lev == "male")))
    pkg <- align(pkg_realloc(fits$lmc, "one-to-one", "absolute",
                             at = list(sex = lev)), hand_a)
    expect_close(pkg[c("estimate", "std.error", "conf.low", "conf.high")],
                 hand_a[c("estimate", "std.error", "conf.low", "conf.high")],
                 1e-6)
  }
})

test_that("the reference contrasted with itself is exactly 0 (link) or 1 (ratio)", {
  for (nm in c("lm0", "glmb")) {
    fit <- fits[[nm]]
    ra <- as.data.frame(predict_realloc(fit, delta = 0, from = "SB",
                                        to = "MVPA", type = "absolute"))
    expect_identical(unlist(ra[2, c("estimate", "std.error", "conf.low",
                                    "conf.high")]),
                     unlist(ra[1, c("estimate", "std.error", "conf.low",
                                    "conf.high")]))
    r0 <- tryCatch(
      suppressWarnings(predict_realloc(fit, delta = 0, from = "SB",
                                       to = "MVPA")),
      error = function(e) e
    )
    expect_s3_class(r0, "realloc")
    if (inherits(r0, "realloc")) {
      expect_identical(r0$estimate, if (nm == "glmb") 1 else 0)
      expect_identical(r0$std.error, 0)
    }
  }
})

test_that("reallocation estimates do not depend on the part order", {
  canon <- function(tab) {
    ref_rows <- which(is.na(tab$from))
    if (length(ref_rows) > 1) tab <- tab[-ref_rows[-1], ]
    tab[order(key_of(tab)), ]
  }
  cols <- c("estimate", "std.error", "conf.low", "conf.high")
  for (nm in c("lm0", "lmc")) {
    cov <- if (nm == "lmc") "sex" else NULL
    for (o in list(rev(parts), c("SB", "MVPA", "SPT", "LPA"))) {
      refit <- comp_iv(declare_composition(tu, o, closure = closure), "zbmi",
                       model = "lm", covariates = cov)
      expect_false(isTRUE(all.equal(refit$data$ilr1, fits[[nm]]$data$ilr1)))
      for (scheme in c("one-to-one", "one-to-all")) {
        for (type in c("relative", "absolute")) {
          a <- canon(pkg_realloc(fits[[nm]], scheme, type))
          b <- canon(pkg_realloc(refit, scheme, type))
          expect_identical(key_of(b), key_of(a))
          expect_close(b[cols], a[cols], 1e-8)
        }
      }
    }
  }
})

test_that("lm reallocations match codaredistlm::predict_delta_comps()", {
  skip_if_not_installed("codaredistlm")
  for (nm in c("lm0", "lmc")) {
    fit <- fits[[nm]]
    covars <- if (nm == "lmc") "sex" else NULL
    crit <- stats::qt(0.975, stats::df.residual(fit$model))
    for (cmp in c("one-v-one", "prop-realloc")) {
      utils::capture.output(suppressMessages(suppressWarnings(
        cr <- codaredistlm::predict_delta_comps(
          dataf = tu, y = "zbmi", comps = parts, covars = covars,
          deltas = deltas / closure, comparisons = cmp
        )
      )))
      cr <- as.data.frame(cr)
      d <- round(cr$delta * closure, 8)
      if (cmp == "one-v-one") {
        pk <- as.data.frame(predict_realloc(fit, delta = deltas,
                                            on_invalid = "na"))
        keys <- paste(cr$`comp-`, cr$`comp+`, d, sep = "|")
      } else {
        pk <- as.data.frame(predict_realloc(fit, delta = deltas,
                                            scheme = "one-to-all", to = parts,
                                            on_invalid = "na"))
        keys <- paste("(rest)", cr$`comp+`, d, sep = "|")
      }
      pk <- pk[match(keys, key_of(pk)), ]
      expect_false(anyNA(pk$delta))
      # Rows that drive a part negative are NA here; codaredistlm warns and
      # still returns a number for them, so they are not compared.
      ok <- !is.na(pk$estimate)
      expect_close(pk$estimate[ok], cr$delta_pred[ok], 1e-6)
      expect_close(pk$conf.low[ok], cr$ci_lo[ok], 1e-6)
      expect_close(pk$conf.high[ok], cr$ci_up[ok], 1e-6)
      expect_close(pk$std.error[ok], (cr$ci_up[ok] - cr$delta_pred[ok]) / crit,
                   1e-6)
    }
  }
})

test_that("summary.realloc() and pub_ready() carry the predict_realloc() numbers", {
  skip_if_not_installed("gt")
  cases <- list(
    predict_realloc(fits$lm0, delta = deltas, on_invalid = "na"),
    predict_realloc(fits$glmb, delta = deltas, on_invalid = "na"),
    predict_realloc(fits$glmp, delta = deltas, scheme = "one-to-all",
                    on_invalid = "na"),
    predict_realloc(fits$lm0, delta = deltas, type = "absolute",
                    on_invalid = "na"),
    predict_realloc(fits$lmc, delta = deltas, by = "sex", on_invalid = "na")
  )
  for (r in cases) {
    s <- summary(r)
    by <- attr(r, "by")
    rk <- function(df) do.call(paste, c(df[c("to", "from", by, "delta")],
                                        sep = "|"))
    rr <- as.data.frame(r)
    rr <- rr[!is.na(rr$from) & !is.na(rr$to) & rr$delta != 0, ]
    st <- as.data.frame(s$table)
    m <- match(rk(st), rk(rr))
    expect_identical(nrow(st), nrow(rr))
    expect_false(anyNA(m))
    expect_identical(anyDuplicated(m), 0L)
    cols <- intersect(c("estimate", "conf.low", "conf.high", "p.value"),
                      names(st))
    expect_identical(unname(as.list(st[cols])), unname(as.list(rr[m, cols])))
    for (dg in c(2, 10)) {
      gd <- as.data.frame(pub_ready(s, digits = dg)$effects[["_data"]])
      label_key <- function(df) do.call(paste, c(df[c("to_label", "from_label",
                                                      by)], sep = "|"))
      gi <- match(label_key(st), label_key(gd))
      expect_false(anyNA(gi))
      expect_identical(nrow(gd), length(unique(label_key(st))))
      fmt <- function(x) formatC(x, format = "f", digits = dg)
      flag <- if ("p.value" %in% names(st)) {
        ifelse(!is.na(st$p.value) & st$p.value < 0.05, "*", "")
      } else {
        ""
      }
      want_est <- ifelse(is.na(st$estimate), NA_character_,
                         paste0(fmt(st$estimate), flag))
      want_ci <- ifelse(is.na(st$conf.low), NA_character_,
                        paste0(fmt(st$conf.low), ", ", fmt(st$conf.high)))
      got_est <- purrr::map_chr(seq_len(nrow(st)), function(i) {
        gd[gi[i], paste0("est_", st$delta[i])]
      })
      got_ci <- purrr::map_chr(seq_len(nrow(st)), function(i) {
        gd[gi[i], paste0("ci_", st$delta[i])]
      })
      expect_identical(got_est, want_est)
      expect_identical(got_ci, want_ci)
    }
  }
})

test_that("plot_realloc() draws exactly the predict_realloc() estimates", {
  delta_max <- 30
  n_points <- 7
  half <- seq(0, delta_max, length.out = ceiling(n_points / 2) + 1)[-1]
  arrow <- intToUtf8(8596L)
  panels_of <- function(p) {
    p <- p$figure
    purrr::map(seq_along(p), function(i) {
      b <- ggplot2::ggplot_build(p[[i]])
      sides <- trimws(strsplit(p[[i]]$labels$title, arrow, fixed = TRUE)[[1]])
      list(row = sides[1], col = sides[2], ribbon = b$data[[1]],
           line = b$data[[2]])
    })
  }
  check <- function(fit, scheme, exclude = NULL, log_y = FALSE) {
    p <- plot_realloc(fit, delta_max = delta_max, n_points = n_points,
                      scheme = scheme, exclude = exclude)
    pd <- panels_of(p)
    kept <- setdiff(parts, exclude)
    expect_length(pd, if (scheme == "one-to-one") choose(length(kept), 2) else
      length(kept))
    inv <- if (log_y) function(y) 10^y else identity
    tol <- if (log_y) 1e-12 else 0
    ref_line <- if (log_y) 1 else 0
    for (panel in pd) {
      for (sign in c(1, -1)) {
        mv <- if (scheme == "one-to-one") {
          if (sign > 0) list(from = panel$col, to = panel$row) else
            list(from = panel$row, to = panel$col)
        } else {
          if (sign > 0) list(from = NULL, to = panel$row) else
            list(from = panel$row, to = NULL)
        }
        want <- as.data.frame(predict_realloc(fit, delta = half,
                                              from = mv$from, to = mv$to,
                                              scheme = scheme,
                                              on_invalid = "na"))
        want <- want[!is.na(want$estimate), ]
        ln <- panel$line[sign * panel$line$x > 0, ]
        rb <- panel$ribbon[sign * panel$ribbon$x > 0, ]
        expect_identical(nrow(ln), nrow(want))
        expect_identical(nrow(rb), nrow(want))
        ml <- match(sign * want$delta, ln$x)
        mr <- match(sign * want$delta, rb$x)
        expect_false(anyNA(c(ml, mr)))
        expect_close(inv(ln$y[ml]), want$estimate, tol)
        expect_close(inv(rb$ymin[mr]), want$conf.low, tol)
        expect_close(inv(rb$ymax[mr]), want$conf.high, tol)
      }
      zero <- c(panel$line$y[panel$line$x == 0],
                panel$ribbon$ymin[panel$ribbon$x == 0],
                panel$ribbon$ymax[panel$ribbon$x == 0])
      expect_close(inv(zero), rep(ref_line, length(zero)), tol)
    }
    pd
  }
  same_panels <- function(a, b) {
    key <- function(pd) purrr::map_chr(pd, function(x) paste(x$row, x$col))
    m <- match(key(b), key(a))
    expect_false(anyNA(m))
    for (i in seq_along(b)) {
      expect_identical(a[[m[i]]]$line[c("x", "y")], b[[i]]$line[c("x", "y")])
      expect_identical(a[[m[i]]]$ribbon[c("x", "ymin", "ymax")],
                       b[[i]]$ribbon[c("x", "ymin", "ymax")])
    }
  }
  same_panels(check(fits$lm0, "one-to-one"),
              check(fits$lm0, "one-to-one", exclude = "SPT"))
  # Excluded parts still absorb their proportional share of a one-to-all move.
  same_panels(check(fits$lm0, "one-to-all"),
              check(fits$lm0, "one-to-all", exclude = "SPT"))
  check(fits$glmb, "one-to-one", log_y = TRUE)
})
