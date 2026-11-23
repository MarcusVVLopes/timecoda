# Numerical benchmark of plot.comp_iv()'s predicted-outcome surface: the grid
# geometry, and the value at every grid composition against
# predict_realloc(type = "absolute", reference =), against predictions and
# contrasts computed by hand from the fitted model on
# robCompositions::pivotCoord() coordinates, and across part orders (ILR
# bases). A squares-only formula is asserted to be order-dependent, so the
# test documents that rather than hides it.

parts <- c("MVPA", "LPA", "SB", "SPT")
tu <- as.data.frame(timeuse)
tu$count <- as.integer(round(exp(pmin(pmax(tu$zbmi, -3), 3))))

fit_set <- function(cd) {
  list(
    g_bin = comp_iv(cd, "obesity", model = "glm", family = "binomial"),
    l_plain = comp_iv(cd, "zbmi", model = "lm"),
    cox = comp_iv(cd, time = "time", event = "event", model = "cox"),
    q_lm = comp_iv(cd, zbmi ~ ilr()^2 + I(ilr()^2), model = "lm"),
    s_lm = comp_iv(cd, zbmi ~ ilr() + I(ilr()^2), model = "lm")
  )
}

# The surface tibble over every part triplet: each panel's data, not drawn.
surface <- function(fit, r, scale) {
  ref <- unlist(fit$ref_comp)[fit$parts]
  purrr::map(utils::combn(fit$parts, 3, simplify = FALSE), function(tri) {
    tc_surface_panel(fit, tri, ref, r, scale, "value")$data
  }) |> purrr::list_rbind()
}
ref_of <- function(fit) unlist(fit$ref_comp)[fit$parts]
# The full composition of each surface row: the triplet shares on their parts,
# every other part at the reference.
comp_raw <- function(fit, s) {
  m <- matrix(ref_of(fit), nrow(s), length(fit$parts), byrow = TRUE,
              dimnames = list(NULL, fit$parts))
  rows <- seq_len(nrow(s))
  m[cbind(rows, match(s$xpart, fit$parts))] <- s$x
  m[cbind(rows, match(s$ypart, fit$parts))] <- s$y
  m[cbind(rows, match(s$zpart, fit$parts))] <- s$z
  m
}
comp_of <- function(fit, s) {
  m <- comp_raw(fit, s)
  m / rowSums(m) * fit$total
}
ilr_of <- function(fit, comps) {
  z <- as.matrix(robCompositions::pivotCoord(
    as.data.frame(comps[, fit$parts, drop = FALSE]), pivotvar = fit$pivot))
  colnames(z) <- paste0("ilr", seq_len(ncol(z)))
  as.data.frame(z)
}
# Hand linear predictor, averaged with equal weights over the cells of the
# model's factor covariates (emmeans' rule); lmer fixed effects only.
eta_of <- function(fit, comps) {
  m <- fit$model
  nd <- ilr_of(fit, comps)
  xl <- if (isS4(m)) NULL else m$xlevels
  cells <- if (length(xl)) expand.grid(xl, stringsAsFactors = FALSE)
  if (!is.null(cells)) {
    nd <- cbind(nd[rep(seq_len(nrow(nd)), each = nrow(cells)), , drop = FALSE],
                cells[rep(seq_len(nrow(cells)), times = nrow(comps)), , drop = FALSE])
  }
  eta <- if (inherits(m, "merMod")) {
    stats::predict(m, newdata = nd, re.form = NA)
  } else if (inherits(m, "coxph")) {
    stats::predict(m, newdata = nd, type = "lp")
  } else if (inherits(m, "glm")) {
    stats::predict(m, newdata = nd, type = "link")
  } else {
    stats::predict(m, newdata = nd)
  }
  if (is.null(cells)) return(as.numeric(eta))
  rowMeans(matrix(as.numeric(eta), nrow = nrow(comps), byrow = TRUE))
}
is_ratio <- function(fit) {
  inherits(fit$model, "coxph") ||
    (inherits(fit$model, "glm") && stats::family(fit$model)$link != "identity")
}
hand_abs <- function(fit, comps) {
  inv <- if (inherits(fit$model, "glm")) stats::family(fit$model)$linkinv else identity
  inv(eta_of(fit, comps))
}
hand_rel <- function(fit, comps) {
  d <- eta_of(fit, comps) - eta_of(fit, rbind(ref_of(fit)))
  if (is_ratio(fit)) exp(d) else d
}
max_diff <- function(a, b) {
  a <- as.numeric(a)
  b <- as.numeric(b)
  if (length(a) != length(b) || anyNA(c(a, b))) return(Inf)
  max(abs(a - b))
}

test_that("the surface grid covers every triplet on the i/r lattice", {
  skip_on_cran()
  cd <- declare_composition(tu, parts, closure = 1440)
  sim <- sim_parts(5)
  pn <- paste0("p", 1:5)
  cases <- list(
    list(comp_iv(cd, "obesity", model = "glm", family = "binomial"), 6),
    list(comp_iv(cd, "obesity", model = "glm", family = "binomial"), 13),
    list(comp_iv(declare_composition(tu, c("MVPA", "LPA", "SB")), "zbmi",
                 model = "lm"), 6),
    list(comp_iv(declare_composition(sim, pn), "y", model = "lm"), 6)
  )
  worst <- max(purrr::map_dbl(cases, function(cs) {
    fit <- cs[[1]]
    r <- cs[[2]]
    s <- surface(fit, r, "relative")
    ref <- ref_of(fit)
    raw <- comp_raw(fit, s)
    held <- purrr::map_dbl(seq_len(nrow(s)), function(k) {
      sum(ref[setdiff(fit$parts, c(s$xpart[k], s$ypart[k], s$zpart[k]))])
    })
    t3 <- s$x + s$y + s$z
    ii <- s$x / t3 * r
    jj <- s$y / t3 * r
    g <- expand.grid(i = seq_len(r - 1), j = seq_len(r - 1))
    g <- g[g$i + g$j < r, ]
    lattice_ok <- purrr::map_lgl(split(seq_len(nrow(s)), s$triplet), function(k) {
      identical(sort(paste(round(ii[k]), round(jj[k]))), sort(paste(g$i, g$j)))
    })
    if (length(unique(s$triplet)) != choose(length(fit$parts), 3) ||
        !all(lattice_ok) || any(raw <= 0)) {
      return(Inf)
    }
    max(abs(ii - round(ii)), abs(jj - round(jj)),
        abs(rowSums(raw) - fit$total), abs(t3 - (fit$total - held)))
  }))
  expect_lte(worst, 1e-10)
})

test_that("the default band draws exactly the grid compositions inside it", {
  skip_on_cran()
  cd <- declare_composition(tu, parts, closure = 1440)
  fit <- comp_iv(cd, "obesity", model = "glm", family = "binomial")
  q <- purrr::map(as.data.frame(tc_comp_df(fit$comp)), stats::quantile,
                  probs = c(0.05, 0.95), names = FALSE)
  within <- function(s) {
    ok <- function(part, v) v >= q[[part]][1] & v <= q[[part]][2]
    purrr::map2_lgl(s$xpart, s$x, ok) & purrr::map2_lgl(s$ypart, s$y, ok) &
      purrr::map2_lgl(s$zpart, s$z, ok)
  }
  drawn <- purrr::map(plot(fit, resolution = 13, smooth = FALSE), "data") |>
    purrr::list_rbind()
  full <- surface(fit, 13, "relative")
  expect_true(all(within(drawn)))
  expect_identical(paste(drawn$triplet, drawn$x, drawn$y),
                   paste(full$triplet, full$x, full$y)[within(full)])
  expect_lte(max_diff(drawn$.pred, full$.pred[within(full)]), 1e-12)
})

test_that("the smoothed panel carries the computed cells unchanged", {
  skip_on_cran()
  cd <- declare_composition(tu, parts, closure = 1440)
  fit <- comp_iv(cd, "obesity", model = "glm", family = "binomial")
  ref <- ref_of(fit)
  tiled <- tc_surface_panel(fit, parts[1:3], ref, 13, "relative", "value")
  smoothed <- tc_surface_panel(fit, parts[1:3], ref, 13, "relative", "value",
                               smooth = TRUE, link = log, inv = exp)
  expect_identical(smoothed$data, tiled$data)
})

test_that("the surface values equal predict_realloc() and hand predictions", {
  skip_on_cran()
  skip_if_not_installed("robCompositions")
  cd <- declare_composition(tu, parts, closure = 1440)
  fits <- list(
    g_bin = comp_iv(cd, "obesity", model = "glm", family = "binomial"),
    g_cov = comp_iv(cd, "obesity", model = "glm", family = "binomial",
                    covariates = c("sex", "region")),
    l_plain = comp_iv(cd, "zbmi", model = "lm"),
    l_cov = comp_iv(cd, "zbmi", model = "lm", covariates = "sex"),
    g_pois = comp_iv(cd, "count", model = "glm", family = "poisson"),
    lmer = comp_iv(cd, "zbmi", model = "lmer", random = "school"),
    q_lm = comp_iv(cd, zbmi ~ ilr()^2 + I(ilr()^2), model = "lm"),
    cox = comp_iv(cd, time = "time", event = "event", model = "cox")
  )
  # Absolute: predict_realloc()'s reference row at each grid composition, and
  # the hand prediction (1e-8 / 1e-6 for lmer, whose fit is iterative).
  abs_err <- purrr::imap_dbl(fits[names(fits) != "cox"], function(fit, nm) {
    s <- surface(fit, 6, "absolute")
    comps <- comp_of(fit, s)
    pr <- purrr::map_dbl(seq_len(nrow(comps)), function(k) {
      predict_realloc(fit, delta = 10, from = "SB", to = "MVPA",
                      type = "absolute", reference = comps[k, ])$estimate[1]
    })
    max(max_diff(s$.pred, pr) / if (nm == "lmer") 1e-8 else 1e-10,
        max_diff(s$.pred, hand_abs(fit, comps)) / if (nm == "lmer") 1e-6 else 1e-10)
  })
  expect_lte(max(abs_err), 1)
  # Relative: the hand contrast exp(eta - eta_ref) or eta - eta_ref (1e-8 on a
  # ratio scale); the Cox lp is centred, which the difference cancels.
  rel_fits <- fits[c("g_bin", "l_plain", "g_pois", "cox", "lmer", "q_lm")]
  rel_err <- purrr::map_dbl(rel_fits, function(fit) {
    s <- surface(fit, 6, "relative")
    max_diff(s$.pred, hand_rel(fit, comp_of(fit, s))) /
      if (is_ratio(fit)) 1e-8 else 1e-10
  })
  expect_lte(max(rel_err), 1)
  # The two scales are one emmeans grid read through the link.
  link_err <- purrr::map_dbl(fits[c("g_bin", "l_plain")], function(fit) {
    sr <- surface(fit, 6, "relative")
    sa <- surface(fit, 6, "absolute")
    a_ref <- predict_realloc(fit, delta = 10, from = "SB", to = "MVPA",
                             type = "absolute")$estimate[1]
    via <- if (inherits(fit$model, "glm")) {
      fam <- stats::family(fit$model)
      exp(fam$linkfun(sa$.pred) - fam$linkfun(a_ref))
    } else {
      sa$.pred - a_ref
    }
    max_diff(sr$.pred, via)
  })
  expect_lte(max(link_err), 1e-10)
})

test_that("the surface on a 5-part composition equals the hand values", {
  skip_on_cran()
  skip_if_not_installed("robCompositions")
  sim <- sim_parts(5)
  pn <- paste0("p", 1:5)
  fit <- comp_iv(declare_composition(sim, pn), "y", model = "lm")
  sa <- surface(fit, 6, "absolute")
  sa <- sa[sa$triplet == sa$triplet[1], ]
  sr <- surface(fit, 6, "relative")
  sr <- sr[sr$triplet == sr$triplet[1], ]
  comps <- comp_of(fit, sa)
  pr <- purrr::map_dbl(seq_len(nrow(comps)), function(k) {
    predict_realloc(fit, delta = fit$total / 144, from = "p3", to = "p1",
                    type = "absolute", reference = comps[k, ])$estimate[1]
  })
  expect_lte(max(max_diff(sa$.pred, pr),
                 max_diff(sr$.pred, hand_rel(fit, comp_of(fit, sr)))), 1e-10)
})

test_that("the surface is invariant to the part order, except squares only", {
  skip_on_cran()
  orders <- list(parts, c("SB", "SPT", "MVPA", "LPA"), c("SPT", "MVPA", "LPA", "SB"))
  by_order <- purrr::map(orders, function(o) {
    fit_set(declare_composition(tu, o, closure = 1440))
  })
  # Spread of .pred across the orders, rows matched by composition; Inf when
  # a row finds no partner.
  spread <- function(nm, scale) {
    tabs <- purrr::map(by_order, function(fs) {
      fit <- fs[[nm]]
      s <- surface(fit, 6, scale)
      key <- apply(round(comp_of(fit, s)[, sort(parts)] / fit$total, 8), 1,
                   paste, collapse = "|")
      list(key = key, pred = s$.pred)
    })
    vals <- purrr::map(tabs, function(t) t$pred[match(tabs[[1]]$key, t$key)])
    vals <- do.call(cbind, vals)
    if (anyNA(vals) || any(lengths(purrr::map(tabs, "key")) != nrow(vals))) {
      return(Inf)
    }
    max(apply(vals, 1, function(v) max(v) - min(v)))
  }
  cases <- list(c("l_plain", "relative", 1e-10), c("l_plain", "absolute", 1e-10),
                c("q_lm", "relative", 1e-10), c("q_lm", "absolute", 1e-10),
                c("g_bin", "relative", 1e-8), c("g_bin", "absolute", 1e-8),
                c("cox", "relative", 1e-8))
  invariant <- purrr::map_dbl(cases, function(cs) {
    spread(cs[1], cs[2]) / as.numeric(cs[3])
  })
  expect_lte(max(invariant), 1)
  s_spread <- c(spread("s_lm", "relative"), spread("s_lm", "absolute"))
  expect_true(all(is.finite(s_spread) & s_spread > 1e-6))
})
