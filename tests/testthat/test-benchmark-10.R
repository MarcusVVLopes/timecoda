# Numerical benchmark of the reallocation anchors and schemes: one-to-all
# (proportional) reallocation at a user-supplied reference composition,
# custom multi-part reallocations, plot_realloc() with a reference and
# `exclude`, and summary() / pub_ready() of those results. The reference is a
# hand computation from coef()/vcov() of the fitted model on
# robCompositions::pivotCoord() coordinates (Kenward-Roger covariance and df
# from pbkrtest for lmer), and codaredistlm::predict_delta_comps() with
# comparisons = "prop-realloc" for the linear models. Reference calls are
# namespace-qualified so each number is visibly computed outside timecoda.

skip_on_cran()
skip_if_not_installed("robCompositions")

parts <- c("MVPA", "LPA", "SB", "SPT")
parts3 <- c("A", "SB", "SPT")
closure <- 1440
ref_u <- c(MVPA = 30, LPA = 240, SB = 660, SPT = 510)
tu <- as.data.frame(timeuse)
tu$A <- tu$MVPA + tu$LPA
cd <- declare_composition(tu, parts, closure = closure)

scen4 <- list(c(SPT = 30, SB = -15, LPA = -15), c(MVPA = 10, SB = -10),
              c(LPA = 20, SPT = -20), c(SPT = 40, MVPA = -40))
spec <- function(fit, kind, p = parts, total = closure, reference = ref_u,
                 deltas = c(5, 15, 30, -15), scen = scen4, cov_x = NULL,
                 pair = c("SB", "MVPA")) {
  list(fit = fit, kind = kind, p = p, total = total, ref_u = reference,
       deltas = deltas, scen = scen, cov_x = cov_x, pair = pair)
}
specs <- list(
  lm = spec(comp_iv(cd, "zbmi", model = "lm"), "lm"),
  lm_cov = spec(comp_iv(cd, "zbmi", model = "lm", covariates = "sex"), "lm",
                cov_x = list(sexmale = 0.5)),
  glm_b = spec(comp_iv(cd, "obesity", model = "glm", family = "binomial"),
               "glm"),
  cox = spec(comp_iv(cd, time = "time", event = "event", model = "cox"), "cox"),
  lm_d3 = spec(comp_iv(declare_composition(tu, parts3, closure = closure),
                       "zbmi", model = "lm"), "lm", p = parts3,
               reference = c(A = 270, SB = 660, SPT = 510),
               scen = list(c(SPT = 30, SB = -15, A = -15), c(A = 10, SB = -10),
                           c(A = 20, SPT = -20)),
               pair = c("SB", "A")),
  lm_prop = spec(comp_iv(declare_composition(tu, parts, closure = 1), "zbmi",
                         model = "lm"), "lm", total = 1,
                 reference = ref_u / closure, deltas = c(0.005, 0.01, -0.005),
                 scen = list(c(SPT = 0.02, SB = -0.01, LPA = -0.01),
                             c(MVPA = 0.005, SB = -0.005),
                             c(SPT = 0.1, SB = -0.03, LPA = -0.07)))
)
# Kenward-Roger is emmeans' default df for lmer only when pbkrtest is present.
if (requireNamespace("pbkrtest", quietly = TRUE)) {
  specs$lmer <- spec(comp_iv(cd, "zbmi", model = "lmer", random = "school"),
                     "lmer")
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

ilr_of <- function(comps, p) {
  z <- as.matrix(robCompositions::pivotCoord(as.data.frame(rbind(comps)),
                                             pivotvar = 1))
  colnames(z) <- paste0("ilr", seq_len(length(p) - 1))
  z
}
unit_1a <- function(ref, focal) {
  rest <- setdiff(names(ref), focal)
  u <- stats::setNames(numeric(length(ref)), names(ref))
  u[focal] <- 1
  u[rest] <- -ref[rest] / sum(ref[rest])
  u
}
hand_ref <- function(p, total) {
  x <- as.matrix(tu[p])
  g <- exp(colMeans(log(x / rowSums(x))))
  g / sum(g) * total
}
model_q <- function(s) {
  m <- s$fit$model
  if (s$kind == "lmer") {
    b <- lme4::fixef(m)
    va <- pbkrtest::vcovAdj(m)
    v <- as.matrix(va)
    dimnames(v) <- list(names(b), names(b))
    ddf <- function(L) pbkrtest::get_Lb_ddf(m, matrix(L, nrow = 1))
  } else {
    b <- stats::coef(m)
    v <- as.matrix(stats::vcov(m))
    df <- if (s$kind == "lm") stats::df.residual(m) else Inf
    ddf <- function(L) df
  }
  list(b = b, v = v[names(b), names(b)], ddf = ddf,
       ratio = s$kind %in% c("glm", "cox"))
}
mq <- lapply(specs, model_q)
rel_cols <- c("estimate", "std.error", "df", "conf.low", "conf.high",
              "statistic", "p.value")
abs_cols <- c("estimate", "std.error", "df", "conf.low", "conf.high")
key_of <- function(df) paste(df$from, df$to, df$delta, sep = "|")

# Relative contrast L'b (or exp(L'b)), SE sqrt(L'VL), against `ref`; a row
# with a non-positive part is NA.
hand_contrast <- function(nm, comps, ref) {
  q <- mq[[nm]]
  p <- specs[[nm]]$p
  comps <- rbind(comps)
  out <- matrix(NA_real_, nrow(comps), 7, dimnames = list(NULL, rel_cols))
  zr <- ilr_of(ref, p)[1, ]
  for (i in which(rowSums(comps <= 0) == 0)) {
    L <- stats::setNames(numeric(length(q$b)), names(q$b))
    L[names(zr)] <- ilr_of(comps[i, ], p)[1, ] - zr
    est <- sum(L * q$b)
    se <- sqrt(as.numeric(t(L) %*% q$v %*% L))
    df <- q$ddf(L)
    crit <- stats::qt(0.975, df)
    tr <- if (q$ratio) exp else identity
    out[i, ] <- c(tr(est), if (q$ratio) exp(est) * se else se, df,
                  tr(est - crit * se), tr(est + crit * se), est / se,
                  2 * stats::pt(-abs(est / se), df))
  }
  as.data.frame(out)
}
# Every focal part as receiver ("to") and as donor ("from"), every delta.
hand_1a <- function(nm, ref) {
  s <- specs[[nm]]
  g <- expand.grid(focal = s$p, role = c("to", "from"), delta = s$deltas,
                   stringsAsFactors = FALSE)
  comps <- t(vapply(seq_len(nrow(g)), function(i) {
    ref + (if (g$role[i] == "to") 1 else -1) * g$delta[i] *
      unit_1a(ref, g$focal[i])
  }, numeric(length(ref))))
  colnames(comps) <- names(ref)
  keys <- data.frame(from = ifelse(g$role == "to", "(rest)", g$focal),
                     to = ifelse(g$role == "to", g$focal, "(rest)"),
                     delta = g$delta)
  list(keys = keys, comps = comps, grid = g)
}
pkg_1a <- function(nm, reference, type = "relative", fit = specs[[nm]]$fit) {
  s <- specs[[nm]]
  as.data.frame(dplyr::bind_rows(
    predict_realloc(fit, delta = s$deltas, scheme = "one-to-all", to = s$p,
                    reference = reference, type = type, on_invalid = "na"),
    predict_realloc(fit, delta = s$deltas, scheme = "one-to-all", from = s$p,
                    reference = reference, type = type, on_invalid = "na")
  ))
}
full_vec <- function(v, p) {
  x <- stats::setNames(numeric(length(p)), p)
  x[names(v)] <- v
  x
}

test_that("one-to-all reallocations at a user reference equal the hand contrast", {
  skip_on_cran()
  num <- comp <- sym <- 0
  for (nm in names(specs)) {
    s <- specs[[nm]]
    h <- hand_1a(nm, s$ref_u)
    hand <- cbind(h$keys, hand_contrast(nm, h$comps, s$ref_u))
    pk <- pkg_1a(nm, s$ref_u)
    pk <- pk[match(key_of(hand), key_of(pk)), ]
    num <- max(num, max_diff(pk[rel_cols], hand[rel_cols]))
    # The package builds the hand compositions, and the remaining parts keep
    # their mutual ratios.
    for (role in c("to", "from")) for (f in s$p) {
      moves <- tc_realloc_schemes[["one-to-all"]]$moves(
        s$p, if (role == "from") f, if (role == "to") f)
      pc <- tc_realloc_grid(s$ref_u, moves, s$deltas, "one-to-all",
                            on_invalid = "na")$comps
      rest <- setdiff(s$p, f)
      comp <- max(comp,
                  max_diff(pc, h$comps[h$grid$focal == f & h$grid$role == role, ]),
                  max_diff(pc[, rest] / pc[, rest[1]],
                           matrix(s$ref_u[rest] / s$ref_u[rest[1]], nrow(pc),
                                  length(rest), byrow = TRUE)))
    }
    # from = focal, delta = d is to = focal, delta = -d.
    a <- as.data.frame(predict_realloc(s$fit, delta = s$deltas,
                                       scheme = "one-to-all", from = s$p,
                                       reference = s$ref_u, on_invalid = "na"))
    b <- as.data.frame(predict_realloc(s$fit, delta = -s$deltas,
                                       scheme = "one-to-all", to = s$p,
                                       reference = s$ref_u, on_invalid = "na"))
    b <- b[match(paste(a$from, -a$delta), paste(b$to, b$delta)), ]
    sym <- max(sym, max_diff(a[rel_cols], b[rel_cols]))
  }
  expect_lte(num, 1e-6)
  expect_lte(comp, 1e-12)
  expect_lte(sym, 1e-12)
})

test_that("absolute one-to-all predictions at a user reference equal x'b", {
  skip_on_cran()
  worst <- 0
  for (nm in c("lm", "lm_cov", "glm_b", "lm_d3", "lm_prop")) {
    s <- specs[[nm]]
    q <- mq[[nm]]
    h <- hand_1a(nm, s$ref_u)
    comps <- rbind(s$ref_u, h$comps)
    ok <- rowSums(comps <= 0) == 0
    X <- matrix(0, nrow(comps), length(q$b), dimnames = list(NULL, names(q$b)))
    X[, "(Intercept)"] <- 1
    X[ok, paste0("ilr", seq_len(length(s$p) - 1))] <- ilr_of(comps[ok, ], s$p)
    for (cv in names(s$cov_x)) X[, cv] <- s$cov_x[[cv]]
    eta <- as.numeric(X %*% q$b)
    se <- sqrt(rowSums((X %*% q$v) * X))
    eta[!ok] <- se[!ok] <- NA
    fam <- if (s$kind == "glm") stats::family(s$fit$model) else stats::gaussian()
    crit <- stats::qt(0.975, q$ddf(NULL))
    hand <- cbind(rbind(data.frame(from = NA, to = NA, delta = 0), h$keys),
                  estimate = fam$linkinv(eta), std.error = fam$mu.eta(eta) * se,
                  df = ifelse(ok, q$ddf(NULL), NA),
                  conf.low = fam$linkinv(eta - crit * se),
                  conf.high = fam$linkinv(eta + crit * se))
    pk <- pkg_1a(nm, s$ref_u, type = "absolute")
    pk <- pk[-which(is.na(pk$from))[-1], ]
    pk <- pk[match(key_of(hand), key_of(pk)), ]
    worst <- max(worst, max_diff(pk[abs_cols], hand[abs_cols]))
  }
  expect_lte(worst, 1e-6)
})

test_that("custom reallocations equal the hand contrast of reference + vector", {
  skip_on_cran()
  num <- reclose <- delta <- 0
  for (nm in names(specs)) {
    s <- specs[[nm]]
    for (reference in list(NULL, s$ref_u)) {
      ref_h <- if (is.null(reference)) hand_ref(s$p, s$total) else reference
      ref_g <- if (is.null(reference)) unlist(s$fit$ref_comp[s$p]) else reference
      full <- lapply(s$scen, full_vec, p = s$p)
      pk <- as.data.frame(predict_realloc(s$fit, realloc = s$scen,
                                          reference = reference,
                                          on_invalid = "na"))
      hand <- hand_contrast(nm, do.call(rbind, lapply(full, `+`, ref_h)), ref_h)
      num <- max(num, max_diff(pk[rel_cols], hand[rel_cols]))
      # Re-closing a balanced reference + vector is the identity.
      gc <- tc_realloc_grid_custom(ref_g, tc_realloc_normalise(s$scen, s$p),
                                   on_invalid = "na")$comps
      reclose <- max(reclose,
                     max_diff(gc, do.call(rbind, lapply(full, `+`, ref_g))))
      delta <- max(delta, max_diff(pk$delta, vapply(full, function(x) {
        sum(x[x > 0])
      }, numeric(1))))
    }
  }
  expect_lte(num, 1e-6)
  expect_lte(reclose, 1e-12)
  expect_identical(delta, 0)
})

test_that("custom vectors reproduce the one-to-one and one-to-all schemes", {
  skip_on_cran()
  worst <- 0
  for (nm in names(specs)) {
    s <- specs[[nm]]
    for (reference in list(NULL, s$ref_u)) {
      ref_g <- if (is.null(reference)) unlist(s$fit$ref_comp[s$p]) else reference
      one <- predict_realloc(s$fit, delta = s$deltas, from = s$pair[1],
                             to = s$pair[2], reference = reference,
                             on_invalid = "na")
      cus <- predict_realloc(s$fit, realloc = lapply(s$deltas, function(d) {
        stats::setNames(c(d, -d), rev(s$pair))
      }), reference = reference, on_invalid = "na")
      worst <- max(worst, max_diff(cus[rel_cols], one[rel_cols]))
      for (f in s$p) {
        o2a <- predict_realloc(s$fit, delta = s$deltas, scheme = "one-to-all",
                               to = f, reference = reference, on_invalid = "na")
        cu <- predict_realloc(s$fit, realloc = lapply(s$deltas, function(d) {
          d * unit_1a(ref_g, f)
        }), reference = reference, on_invalid = "na")
        worst <- max(worst, max_diff(cu[rel_cols], o2a[rel_cols]))
      }
    }
    # A data frame of scenarios gives the list's numbers.
    dfm <- as.data.frame(do.call(rbind, lapply(s$scen, full_vec, p = s$p)))
    worst <- max(worst, max_diff(
      predict_realloc(s$fit, realloc = dfm, reference = s$ref_u,
                      on_invalid = "na")[rel_cols],
      predict_realloc(s$fit, realloc = s$scen, reference = s$ref_u,
                      on_invalid = "na")[rel_cols]))
  }
  # The two routes agree to rounding error, whose size differs by platform.
  expect_lte(worst, 1e-10)
})

test_that("one-to-all lm reallocations match codaredistlm prop-realloc, both directions", {
  skip_on_cran()
  skip_if_not_installed("codaredistlm")
  worst <- 0
  for (nm in c("lm", "lm_cov", "lm_d3", "lm_prop")) {
    s <- specs[[nm]]
    utils::capture.output(suppressMessages(suppressWarnings(
      cr <- codaredistlm::predict_delta_comps(
        dataf = tu, y = "zbmi", comps = s$p,
        covars = if (nm == "lm_cov") "sex" else NULL,
        deltas = s$deltas / s$total, comparisons = "prop-realloc")
    )))
    cr <- as.data.frame(cr)
    d <- round(cr$delta * s$total, 10)
    crit <- stats::qt(0.975, stats::df.residual(s$fit$model))
    want <- cbind(cr$delta_pred, cr$ci_lo, cr$ci_up,
                  (cr$ci_up - cr$delta_pred) / crit)
    # The focal part gains d (a negative d: it gives |d|) ...
    pk <- as.data.frame(predict_realloc(s$fit, delta = s$deltas,
                                        scheme = "one-to-all", to = s$p,
                                        on_invalid = "na"))
    pk <- pk[match(paste(cr$`comp+`, d), paste(pk$to, round(pk$delta, 10))), ]
    worst <- max(worst, max_diff(pk[c("estimate", "conf.low", "conf.high",
                                      "std.error")], want))
    # ... and a negative d is the focal part giving |d| to the rest.
    neg <- cr$delta < 0
    pg <- as.data.frame(predict_realloc(s$fit, delta = -unique(d[neg]),
                                        scheme = "one-to-all", from = s$p,
                                        on_invalid = "na"))
    pg <- pg[match(paste(cr$`comp+`[neg], -d[neg]),
                   paste(pg$from, round(pg$delta, 10))), ]
    worst <- max(worst, max_diff(pg[c("estimate", "conf.low", "conf.high",
                                      "std.error")], want[neg, ]))
  }
  expect_lte(worst, 1e-6)
})

test_that("plot_realloc() at a user reference draws exactly predict_realloc()", {
  skip_on_cran()
  delta_max <- 60
  n_points <- 9
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
  values <- shape <- filtered <- 0
  check <- function(fit, scheme, exclude = NULL, log_y = FALSE) {
    pd <- panels_of(plot_realloc(fit, delta_max = delta_max,
                                 n_points = n_points, scheme = scheme,
                                 exclude = exclude, reference = ref_u))
    kept <- setdiff(parts, exclude)
    shape <<- shape + abs(length(pd) - if (scheme == "one-to-one")
      choose(length(kept), 2) else length(kept))
    inv <- if (log_y) function(y) 10^y else identity
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
                                              scheme = scheme, reference = ref_u,
                                              on_invalid = "na"))
        want <- want[!is.na(want$estimate), ]
        ln <- panel$line[sign * panel$line$x > 0, ]
        rb <- panel$ribbon[sign * panel$ribbon$x > 0, ]
        ml <- match(sign * want$delta, ln$x)
        mr <- match(sign * want$delta, rb$x)
        shape <<- shape + sum(is.na(c(ml, mr))) +
          abs(nrow(ln) - nrow(want)) + abs(nrow(rb) - nrow(want))
        values <<- max(values, max_diff(inv(ln$y[ml]), want$estimate),
                       max_diff(inv(rb$ymin[mr]), want$conf.low),
                       max_diff(inv(rb$ymax[mr]), want$conf.high))
        # The largest magnitude drawn is the largest one that keeps every
        # part positive at the reference.
        feas <- half[vapply(half, function(d) {
          x <- ref_u
          if (scheme == "one-to-one") {
            x[mv$from] <- x[mv$from] - d
            x[mv$to] <- x[mv$to] + d
          } else if (is.null(mv$to)) {
            x <- x - d * unit_1a(ref_u, mv$from)
          } else {
            x <- x + d * unit_1a(ref_u, mv$to)
          }
          all(x > 0)
        }, logical(1))]
        shape <<- shape +
          abs(max(c(0, abs(ln$x))) - max(c(0, feas)))
      }
      zero <- c(panel$line$y[panel$line$x == 0],
                panel$ribbon$ymin[panel$ribbon$x == 0],
                panel$ribbon$ymax[panel$ribbon$x == 0])
      values <<- max(values, max(abs(inv(zero) - if (log_y) 1 else 0)))
    }
    pd
  }
  same_panels <- function(a, b) {
    key <- function(pd) purrr::map_chr(pd, function(x) paste(x$row, x$col))
    m <- match(key(b), key(a))
    if (anyNA(m)) return(Inf)
    max(0, purrr::map2_dbl(a[m], b, function(x, y) {
      max(max_diff(x$line[c("x", "y")], y$line[c("x", "y")]),
          max_diff(x$ribbon[c("x", "ymin", "ymax")],
                   y$ribbon[c("x", "ymin", "ymax")]))
    }))
  }
  for (nm in c("lm", "glm_b")) {
    for (scheme in c("one-to-one", "one-to-all")) {
      full <- check(specs[[nm]]$fit, scheme, log_y = nm == "glm_b")
      excl <- check(specs[[nm]]$fit, scheme, exclude = "MVPA",
                    log_y = nm == "glm_b")
      filtered <- max(filtered, same_panels(full, excl))
    }
  }
  expect_identical(shape, 0)
  expect_lte(values, 1e-12)
  expect_identical(filtered, 0)
})

test_that("summary() and pub_ready() of one-to-all and custom results carry the numbers", {
  skip_on_cran()
  skip_if_not_installed("gt")
  s4 <- specs$lm
  cases <- list(
    predict_realloc(s4$fit, delta = c(5, 15, 30), scheme = "one-to-all",
                    reference = ref_u, on_invalid = "na"),
    predict_realloc(s4$fit, delta = c(5, 15, 30), scheme = "one-to-all",
                    from = parts, reference = ref_u, on_invalid = "na"),
    predict_realloc(s4$fit, delta = c(15, 30), scheme = "one-to-all",
                    type = "absolute", reference = ref_u, on_invalid = "na"),
    predict_realloc(specs$glm_b$fit, realloc = scen4, on_invalid = "na"),
    predict_realloc(specs$cox$fit, realloc = scen4[1:3], reference = ref_u,
                    on_invalid = "na"),
    predict_realloc(specs$lm_prop$fit, realloc = specs$lm_prop$scen,
                    on_invalid = "na")
  )
  # ?summary.realloc: the pooled side of a one-to-all reallocation is shown
  # as "All other parts"; part labels replace part names, also inside a
  # custom reallocation's compact strings.
  rest_label <- "All other parts"
  labs <- c(MVPA = "Exercise", SB = "Sitting")
  want_label <- function(x, labels) {
    out <- ifelse(x == "(rest)", rest_label, x)
    for (pt in names(labels)) {
      out[x == pt] <- labels[[pt]]
      out <- gsub(paste0("(^|, )", pt, " ([+-])"),
                  paste0("\\1", labels[[pt]], " \\2"), out)
    }
    out
  }
  values <- lab_bad <- cells <- 0
  for (r in cases) {
    s <- summary(r)
    rr <- as.data.frame(r)
    rr <- rr[!is.na(rr$from) & !is.na(rr$to) & rr$delta != 0, ]
    st <- as.data.frame(s$table)
    m <- match(paste(st$to, st$from, st$delta), paste(rr$to, rr$from, rr$delta))
    cols <- intersect(c("estimate", "conf.low", "conf.high", "p.value"),
                      names(st))
    values <- max(values, abs(nrow(st) - nrow(rr)) + anyDuplicated(m),
                  max_diff(st[cols], rr[m, cols]))
    sl <- as.data.frame(summary(r, labels = labs)$table)
    lab_bad <- lab_bad +
      sum(st$to_label != want_label(st$to, NULL)) +
      sum(st$from_label != want_label(st$from, NULL)) +
      sum(sl$to_label != want_label(sl$to, labs)) +
      sum(sl$from_label != want_label(sl$from, labs))
    for (dg in c(2, 10)) {
      gd <- as.data.frame(pub_ready(s, digits = dg)$effects[["_data"]])
      gi <- match(paste(st$to_label, st$from_label),
                  paste(gd$to_label, gd$from_label))
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
      cells <- cells + sum(is.na(gi)) +
        abs(nrow(gd) - length(unique(paste(st$to_label, st$from_label)))) +
        sum(!mapply(identical, got_est, want_est)) +
        sum(!mapply(identical, got_ci, want_ci))
    }
  }
  expect_identical(values, 0)
  expect_identical(lab_bad, 0)
  expect_identical(cells, 0)
})

test_that("one-to-all and custom estimates do not depend on the part order", {
  skip_on_cran()
  canon <- function(tab) tab[order(key_of(tab)), ]
  orders <- list(rev(parts), c(parts[-1], parts[1]), c("SB", "MVPA", "SPT", "LPA"))
  cols <- c("estimate", "std.error", "conf.low", "conf.high")
  req <- function(fit, reference) {
    list(o2a = canon(as.data.frame(predict_realloc(
      fit, delta = c(5, 15, 30, -15), scheme = "one-to-all", to = "MVPA",
      reference = reference, on_invalid = "na"))),
      custom = as.data.frame(predict_realloc(
        fit, realloc = scen4, reference = reference, on_invalid = "na")))
  }
  worst <- c(lm = 0, glm_b = 0)
  for (nm in names(worst)) {
    args <- if (nm == "lm") list("zbmi", model = "lm") else
      list("obesity", model = "glm", family = "binomial")
    for (o in orders) {
      refit <- do.call(comp_iv, c(list(declare_composition(tu, o,
                                                           closure = closure)),
                                  args))
      for (reference in list(NULL, ref_u)) {
        a <- req(specs[[nm]]$fit, reference)
        b <- req(refit, reference)
        same_keys <- identical(key_of(a$o2a), key_of(b$o2a))
        worst[[nm]] <- max(worst[[nm]], if (same_keys) 0 else Inf,
                           max_diff(b$o2a[cols], a$o2a[cols]),
                           max_diff(b$custom[cols], a$custom[cols]))
      }
    }
  }
  expect_lte(worst[["lm"]], 1e-8)
  # The glm is refitted iteratively in each basis.
  expect_lte(worst[["glm_b"]], 1e-3)
})
