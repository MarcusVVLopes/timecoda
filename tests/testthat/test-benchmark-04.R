# Numerical benchmark of fit_check() for comp_iv() lm and glm fits against
# DHARMa called directly on the same fitted model: every statistic and
# p-value of the uniformity, dispersion and outlier rows, and the scaled
# residuals, for an lm (250 and 50 simulations), a binomial glm and a Poisson
# glm. Also checks the returned convergence state against the fitted model, and
# that the diagnostics are reproducible: DHARMa::simulateResiduals() fixes
# its own seed (123 by default) and restores the caller's random state, so the
# simulation rows do not move with set.seed(), while DHARMa with seed = NULL
# does. Reference calls are namespace-qualified so each number is visibly
# computed outside timecoda.

skip_if_not_installed("DHARMa")
skip_on_cran()

parts <- c("MVPA", "LPA", "SB", "SPT")
seed <- 20260923
df <- timeuse
df$count <- as.integer(round(exp(pmin(pmax(df$zbmi, -3), 3))))
cd <- declare_composition(df, parts, closure = 1440)
fits <- list(
  lm = comp_iv(cd, outcome = "zbmi", model = "lm"),
  glm_binomial = comp_iv(cd, outcome = "obesity", model = "glm",
                         family = "binomial"),
  glm_poisson = comp_iv(cd, outcome = "count", model = "glm",
                        family = "poisson")
)
dharma_labels <- c("KS uniformity", "Dispersion", "Outliers")

# DHARMa run by hand in the order fit_check() runs it: testOutliers() draws
# bootstrap samples for an integer response, so the order fixes the stream.
direct_dharma <- function(model, n) {
  set.seed(seed)
  sr <- DHARMa::simulateResiduals(model, n = n, plot = FALSE)
  u <- DHARMa::testUniformity(sr, plot = FALSE)
  d <- DHARMa::testDispersion(sr, plot = FALSE)
  o <- DHARMa::testOutliers(sr, plot = FALSE)
  list(sr = sr,
       statistic = unname(c(u$statistic, d$statistic, o$statistic)),
       p.value = unname(c(u$p.value, d$p.value, o$p.value)))
}

test_that("the DHARMa rows match DHARMa called directly", {
  cases <- list(list(fits$lm, 250), list(fits$lm, 50),
                list(fits$glm_binomial, 250), list(fits$glm_poisson, 250))
  for (cs in cases) {
    fit <- cs[[1]]
    n <- cs[[2]]
    set.seed(seed)
    fc <- fit_check(fit, plot = FALSE, n = n)
    ref <- direct_dharma(fit$model, n)
    expect_identical(fc$tests$check, dharma_labels)
    rows <- fc$tests
    expect_equal(rows$statistic, ref$statistic, tolerance = 1e-6)
    expect_equal(rows$p.value, ref$p.value, tolerance = 1e-6)
    expect_true(all(is.na(rows$df)))
    expect_equal(fc$raw$scaledResiduals, ref$sr$scaledResiduals,
                 tolerance = 1e-6)
    expect_identical(fc$raw$nSim, n)
  }
})

test_that("the convergence state reflects the fitted model", {
  for (nm in names(fits)) {
    fit <- fits[[nm]]
    conv <- fit_check(fit, plot = FALSE)$convergence
    if (nm == "lm") {
      expect_true(is.na(conv$converged))
    } else {
      expect_identical(conv$converged, isTRUE(fit$model$converged))
    }
    expect_equal(conv$n_warnings, nrow(fit$fit_warnings))
  }
})

test_that("the diagnostics are reproducible and DHARMa's own seed fixes them", {
  for (fit in fits) {
    set.seed(seed)
    a <- fit_check(fit, plot = FALSE)
    set.seed(seed)
    b <- fit_check(fit, plot = FALSE)
    expect_identical(a$tests, b$tests)
    expect_identical(a$raw$scaledResiduals, b$raw$scaledResiduals)
    # Another global seed leaves the simulation, and so the uniformity and
    # dispersion rows and the outlier count, unchanged.
    set.seed(seed + 1)
    c2 <- fit_check(fit, plot = FALSE)
    keep <- match(c("KS uniformity", "Dispersion"), a$tests$check)
    expect_identical(a$tests[keep, ], c2$tests[keep, ])
    outliers <- match("Outliers", a$tests$check)
    expect_identical(a$tests$statistic[outliers], c2$tests$statistic[outliers])
    expect_identical(a$raw$scaledResiduals, c2$raw$scaledResiduals)
    sr123 <- DHARMa::simulateResiduals(fit$model, n = 250, plot = FALSE,
                                       seed = 123)
    expect_equal(a$raw$scaledResiduals, sr123$scaledResiduals,
                 tolerance = 1e-6)
    # Without the fixed seed the simulation is Monte Carlo.
    set.seed(seed)
    s1 <- DHARMa::simulateResiduals(fit$model, n = 250, plot = FALSE,
                                    seed = NULL)
    set.seed(seed + 1)
    s2 <- DHARMa::simulateResiduals(fit$model, n = 250, plot = FALSE,
                                    seed = NULL)
    expect_gt(max(abs(s1$scaledResiduals - s2$scaledResiduals)), 0)
  }
})
