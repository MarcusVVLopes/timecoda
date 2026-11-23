# data-raw/timeuse.R
#
# Generate `timeuse`, a simulated 24-hour time-use dataset used as a fixture
# for timecoda examples and tests, so they never depend on real data.
#
# Design:
#   * 300 individuals.
#   * 4 movement-behaviour parts summing to 1440 min/day:
#       MVPA  moderate-to-vigorous physical activity
#       LPA   light physical activity
#       SB    sedentary behaviour
#       SPT   sleep period time
#   * Outcomes:
#       obesity  binary factor ("no"/"yes")
#       zbmi   continuous (BMI z-score)
#       time   follow-up time to onset of obesity (survival)
#       event  obesity-onset indicator (1 = onset, 0 = censored)
#   * sex     grouping factor ("male"/"female") for colouring / as a covariate
#   * school  clustering factor (10 schools) for multilevel (lme4) models
#   * region  clustering factor (3 regions) deterministically nesting the 10
#             schools, so `random = c("region", "school")` demonstrates several
#             random intercepts
#   * grade   ordered academic-achievement band ("below" < "expected" < "above")
#             for the ordinal (and, coerced unordered, multinomial) model paths
#
# RNG discipline: the seed and the order of every existing random draw are held
# fixed so that all previously shipped columns stay byte-identical. The only new
# random draw (the `grade` latent-score noise) is appended AFTER every existing
# draw; `region` is a deterministic mapping and costs no RNG.
#
# Run with:  source("data-raw/timeuse.R")
# Regenerate the .rda with usethis when developing the package.

set.seed(2026)

n <- 300L

# --- Grouping and clustering -------------------------------------------------
sex    <- factor(sample(c("male", "female"), n, replace = TRUE),
                 levels = c("female", "male"))
school <- factor(sample(sprintf("school%02d", 1:10), n, replace = TRUE),
                 levels = sprintf("school%02d", 1:10))

# School-level random intercept on the latent BMI scale.
school_re <- rnorm(nlevels(school), mean = 0, sd = 0.35)
names(school_re) <- levels(school)

# --- Composition (4 parts, sum = 1440 min) -----------------------------------
# Dirichlet draws via independent gammas; mean proportions chosen to mimic a
# plausible 24-h movement composition, then rescaled to a 1440-min total.
total      <- 1440
mean_prop  <- c(MVPA = 0.03, LPA = 0.17, SB = 0.45, SPT = 0.35)
conc       <- 35                       # concentration: higher = less variable
alpha      <- conc * mean_prop

g    <- sapply(alpha, function(a) rgamma(n, shape = a, rate = 1))
prop <- g / rowSums(g)
comp <- prop * total

MVPA <- comp[, "MVPA"]
LPA  <- comp[, "LPA"]
SB   <- comp[, "SB"]
SPT  <- comp[, "SPT"]

# --- Latent BMI score driven by the composition -----------------------------
# Use centred log-ratios so the signal is genuinely compositional: more MVPA
# and more sleep lower BMI; more sedentary time raises it.
clr <- log(prop) - rowMeans(log(prop))

lp <-
  0.00 +
  -0.55 * scale(clr[, "MVPA"])[, 1] +
  -0.25 * scale(clr[, "SPT"])[, 1] +
   0.45 * scale(clr[, "SB"])[, 1]  +
   0.20 * (sex == "male") +
   school_re[as.character(school)]

zbmi <- as.numeric(lp + rnorm(n, sd = 0.8))

# --- Binary outcome: obesity --------------------------------------------------
p_obesity <- plogis(-1.1 + 1.0 * lp)
obesity   <- factor(ifelse(rbinom(n, 1, p_obesity) == 1L, "yes", "no"),
                  levels = c("no", "yes"))

# --- Survival outcome: time to onset of obesity ------------------------------
# Weibull event times; higher BMI risk score -> earlier onset. Administrative
# censoring at the study end plus light random dropout.
shape    <- 1.4
base_rate <- 0.015
event_t  <- rweibull(n, shape = shape,
                     scale = (base_rate * exp(0.8 * lp))^(-1 / shape))
cens_t   <- pmin(60, rexp(n, rate = 1 / 90))     # study ends at 60 months
time     <- pmin(event_t, cens_t)
event    <- as.integer(event_t <= cens_t)

# --- Region: a second, coarser cluster level nesting the 10 schools ----------
# Deterministic (no RNG): schools 01-03 -> A, 04-06 -> B, 07-10 -> C. This makes
# `random = c("region", "school")` an honest two-level (region > school) demo.
school_idx <- as.integer(sub("school", "", as.character(school)))
region     <- factor(
  dplyr::case_when(
    school_idx <= 3L ~ "A",
    school_idx <= 6L ~ "B",
    TRUE             ~ "C"
  ),
  levels = c("A", "B", "C")
)

# --- Ordinal outcome: academic-achievement band -----------------------------
# Reuse the composition-driven linear predictor `lp` (negated, so a higher
# achievement band corresponds to the healthier composition: more MVPA/sleep,
# less sedentary time) plus fresh Gaussian noise, then cut at the tertiles into
# an ordered 3-band factor. The noise SD (1.8) is tuned so the ordinal
# composition omnibus is clearly significant (p ~ 3e-5) without a degenerate
# fit. This is the LAST random draw in the script, appended after every existing
# draw so all previously shipped columns stay byte-identical.
grade_lat <- -lp + rnorm(n, sd = 1.8)
grade_q   <- stats::quantile(grade_lat, probs = c(1 / 3, 2 / 3))
grade     <- cut(
  grade_lat,
  breaks = c(-Inf, grade_q, Inf),
  labels = c("below", "expected", "above"),
  ordered_result = TRUE
)

timeuse <- tibble::tibble(
  id     = sprintf("p%03d", 1:n),
  sex    = sex,
  school = school,
  region = region,
  MVPA   = round(MVPA, 1),
  LPA    = round(LPA, 1),
  SB     = round(SB, 1),
  SPT    = round(SPT, 1),
  zbmi   = round(zbmi, 3),
  obesity = obesity,
  time   = round(time, 2),
  event  = event,
  grade  = grade
)

# Guarantee the parts sum exactly to 1440 after rounding by absorbing the
# rounding residual into the largest part (SB).
parts_sum <- timeuse$MVPA + timeuse$LPA + timeuse$SB + timeuse$SPT
timeuse$SB <- timeuse$SB + (total - parts_sum)

stopifnot(
  nrow(timeuse) == n,
  all(abs(timeuse$MVPA + timeuse$LPA + timeuse$SB + timeuse$SPT - total) < 1e-6),
  all(timeuse$MVPA > 0, timeuse$LPA > 0, timeuse$SB > 0, timeuse$SPT > 0),
  # region deterministically nests school (each school maps to one region)
  is.factor(timeuse$region), nlevels(timeuse$region) == 3L,
  all(rowSums(table(timeuse$school, timeuse$region) > 0) == 1L),
  # grade is an ordered 3-band factor with every band populated
  is.ordered(timeuse$grade),
  identical(levels(timeuse$grade), c("below", "expected", "above")),
  all(table(timeuse$grade) > 0)
)

usethis::use_data(timeuse, overwrite = TRUE)
