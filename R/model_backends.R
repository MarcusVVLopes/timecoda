# Model-backend registry ------------------------------------------------------
#
# comp_iv() supports several model classes. The per-type behaviour its helpers
# need (fit label, effect scale, coefficient extraction, emmeans mode) lives in
# one place, the `tc_model_registry` list below, read through
# `tc_model_meta()`, so no helper carries a `switch(model_type, ...)`. Adding a
# model family is then a new registry row.
#
# `model_type` is the engine, and it is also exactly what the user passes as
# `model`: "lm", "glm", "cox"/"coxme", "lmer", "glmer", "clm"/"clmm" (ordinal),
# "mblogit" (multinomial), "lmrob"/"rlmer" (robust linear). There is one
# vocabulary, so a registry row can be looked up straight from user input,
# which lets comp_iv() validate its arguments against the registry before
# fitting.

# Registry of per-type metadata. Fields:
#   label        human-readable model name for print/summary headers
#   show_family  append "(<family>)" to the label, plus ", <link> link" when
#                the link is not that family's default
#   scale        effect scale, drives exponentiation + the effect label:
#                "difference" (identity), "family" (glm rule), "hazard",
#                "odds" (ordinal cumulative OR), "rrr" (multinomial)
#   emm_mode     `mode=` passed to emmeans::ref_grid (e.g. "latent" for ordinal)
#   link_tran    transformation tag applied to the emmGrid so `type = "response"`
#                back-transforms the contrast (ordinal's latent contrast -> OR)
#   emm_via      how to build the reference grid: "ref_grid" (default) or
#                "qdrg" (lmrob has no native emmeans support)
#   per_category TRUE for multinomial (one effect set per outcome category)
#   fixef        extract fixed effects with lme4::fixef() rather than coef()
#   absolute_ok  FALSE when the class has no absolute-response scale (cox, ...)
#   absolute_note why it has none, named in the "relative only" messages
#   check_backend  which fit_check() diagnostic backend to use
#   limited_note the explanatory note fit_check() prints for a class whose
#                backend is "limited" (fit statistics only); NULL gives a
#                generic note
#   icc          which between-cluster measure summary() reports: "lmm",
#                "glmm" (latent ICC), "clmm" (latent ICC, link read from the
#                fit's own $link), "coxme" (median hazard ratio) or "none"
#   emm_df       degrees of freedom forced on every emmeans grid of this engine
#                (`Inf` = asymptotic), when the backend's own df method is not
#                usable; NULL leaves emmeans' choice
#   convergence  how tc_iv_fit() probes the fit for convergence: "lme4"
#                (@optinfo messages + isSingular(); these engines also take an
#                lme4 optimizer through `control`), "flag" ($converged),
#                "ordinal" ($convergence code), "optinfo" (robustlmm:
#                @optinfo$conv$opt, 0 = converged) or "none" (reported as NA)
#
# The next six say how the engine is called; tc_iv_fit() reads them, so there
# is no per-engine code. `fitter` is stored quoted so the fitted object's own
# `call` still names the backend (ordinal::nominal_test() re-fits through it, so
# a function value there would break fit_check()).
#   fitter        the backend, quoted: quote(stats::lm)
#   fixed_args    arguments the backend always receives, e.g. list(trace = FALSE)
#                 (comp_iv() refuses a `model_args` element of the same name)
#   response      "outcome" (default) or "surv" (Surv(time, event) LHS)
#   random_style  "none", "bars" ((1 | g) terms in the formula), or "arg"
#                 (the backend's own `random = ~ 1 | g` argument)
#   outcome_check "ordered" / "unordered" when the engine constrains the
#                 outcome's type; NULL when it does not
#   quiet         TRUE when the backend prints to the console during a fit and
#                 no argument turns that off (a clustered mblogit fit prints
#                 blank lines even with `trace = FALSE`); tc_iv_fit() then
#                 discards that output
#
# The last four say what the user may pass alongside this key. They are
# meaningful because `model` is the key: naming `clmm` asks for the
# clustered engine, so `random` becomes required rather than merely permitted.
#   allows_family    `family` is accepted
#   allows_random    `random` is accepted
#   requires_random  `random` is mandatory
#   pkg              Suggested package the engine needs (checked at call time)
tc_model_registry <- list(
  lm = list(label = "linear model", scale = "difference",
            check_backend = "dharma",
            fitter = quote(stats::lm)),
  glm = list(label = "generalized linear model", show_family = TRUE,
             scale = "family", allows_family = TRUE,
             check_backend = "dharma", convergence = "flag",
             fitter = quote(stats::glm)),
  # `model = TRUE` keeps the model frame on the fit, which cox.zph() and emmeans
  # then read directly instead of re-evaluating the stored call.
  cox = list(label = "Cox proportional hazards model", scale = "hazard",
             absolute_ok = FALSE,
             absolute_note = "a Cox model leaves the baseline hazard unestimated",
             check_backend = "coxzph",
             fitter = quote(survival::coxph), response = "surv",
             fixed_args = list(model = TRUE)),
  coxme = list(label = "mixed-effects Cox proportional hazards model",
               scale = "hazard", absolute_ok = FALSE, fixef = TRUE,
               absolute_note = "a Cox model leaves the baseline hazard unestimated",
               allows_random = TRUE, requires_random = TRUE, pkg = "coxme",
               icc = "coxme",
               check_backend = "limited",
               limited_note = "The proportional-hazards test and residual diagnostics are not run for mixed-effects Cox (coxme) models; only fit statistics are reported.",
               fitter = quote(coxme::coxme), response = "surv",
               random_style = "bars"),
  lmer = list(label = "linear mixed model", scale = "difference",
              fixef = TRUE, allows_random = TRUE, requires_random = TRUE,
              check_backend = "dharma", convergence = "lme4", icc = "lmm",
              fitter = quote(lme4::lmer), random_style = "bars"),
  glmer = list(label = "generalized linear mixed model", show_family = TRUE,
               scale = "family", fixef = TRUE, allows_family = TRUE,
               allows_random = TRUE, requires_random = TRUE,
               check_backend = "dharma", convergence = "lme4", icc = "glmm",
               fitter = quote(lme4::glmer), random_style = "bars"),
  clm = list(label = "cumulative-link (ordinal) model", scale = "odds",
             emm_mode = "latent", link_tran = "log", absolute_ok = FALSE,
             absolute_note = "ordinal predictions are per outcome category",
             pkg = "ordinal", check_backend = "ordinal",
             convergence = "ordinal",
             fitter = quote(ordinal::clm),
             outcome_check = "ordered"),
  clmm = list(label = "cumulative-link mixed (ordinal) model", scale = "odds",
              emm_mode = "latent", link_tran = "log", absolute_ok = FALSE,
              absolute_note = "ordinal predictions are per outcome category",
              allows_random = TRUE, requires_random = TRUE,
              pkg = "ordinal", check_backend = "limited",
              limited_note = "Residual diagnostics and the proportional-odds test are not available for mixed ordinal (clmm) models; only fit statistics are reported.",
              convergence = "ordinal", icc = "clmm",
              fitter = quote(ordinal::clmm), fixed_args = list(Hess = TRUE),
              outcome_check = "ordered", random_style = "bars"),
  mblogit = list(label = "multinomial baseline-logit model", scale = "rrr",
                 per_category = TRUE, link_tran = "log", absolute_ok = FALSE,
                 absolute_note = "multinomial predictions are per outcome category",
                 allows_random = TRUE, pkg = "mclogit",
                 check_backend = "limited", convergence = "flag",
                 limited_note = "Formal residual diagnostics are not available for multinomial (mblogit) models; only fit statistics are reported.",
                 fitter = quote(mclogit::mblogit),
                 fixed_args = list(trace = FALSE), quiet = TRUE,
                 outcome_check = "unordered", random_style = "arg"),
  lmrob = list(label = "robust linear model", scale = "difference",
               emm_via = "qdrg", pkg = "robustbase",
               check_backend = "robust", convergence = "flag",
               fitter = quote(robustbase::lmrob)),
  rlmer = list(label = "robust linear mixed model", scale = "difference",
               fixef = TRUE, allows_random = TRUE, requires_random = TRUE,
               pkg = "robustlmm", check_backend = "robust_mixed", icc = "lmm",
               emm_df = Inf, convergence = "optinfo",
               fitter = quote(robustlmm::rlmer), random_style = "bars")
)

# Filled in under every registry row, so callers read any field without an
# is.null() guard and a row states only what makes its engine different.
tc_model_defaults <- list(
  show_family = FALSE, emm_mode = NULL, link_tran = NULL, emm_via = "ref_grid",
  per_category = FALSE, fixef = FALSE,
  allows_family = FALSE, allows_random = FALSE, requires_random = FALSE,
  pkg = NULL,
  absolute_ok = TRUE, absolute_note = NULL,
  check_backend = "limited", limited_note = NULL, convergence = "none",
  icc = "none", emm_df = NULL,
  fitter = NULL, fixed_args = list(), response = "outcome",
  random_style = "none", outcome_check = NULL, quiet = FALSE
)

# One registry row, defaults filled in.
tc_model_meta <- function(model_type) {
  if (!model_type %in% names(tc_model_registry)) {
    cli::cli_abort(
      c("Unknown model type {.val {model_type}}.",
        "i" = "Known types: {.val {names(tc_model_registry)}}.")
    )
  }
  utils::modifyList(tc_model_defaults, tc_model_registry[[model_type]])
}

# The GLM links comp_iv() accepts, each with the effect it reports: identity
# gives a mean difference; logit an odds ratio; log a ratio, whose name depends
# on the family. Any other link has no such reading, so comp_iv() refuses it.
tc_supported_links <- c("identity", "log", "logit")

# How this fit is described, as opposed to how it is fitted: the header label,
# whether the per-part effect is exponentiated for reporting, and what that
# effect is then called. All three are derivations of one registry row plus the
# GLM family, and every consumer that needs one needs at least one other, so
# they are answered together rather than by three separate lookups.
#
# Kept apart from tc_model_meta() deliberately: two of the three answers depend
# on `family`, which most tc_model_meta() call sites do not have. As fields on
# that list they would come back silently wrong wherever the family was not
# passed; here the signature asks for what the answer needs. `family` is the
# family object a fit stores (`fit$family`), or a family name (resolved to its
# default link), or NULL. For a GLM engine given no family, only the label is
# answerable, so `exponentiate` and `measure` come back NA.
#
# Read by print.timecoda_fit(), summary.comp_iv(), fit_check(), plot.comp_iv()
# and plot_realloc() (Modules 1, 2 and 4), which is why it lives beside the
# registry rather than in any one of them.
tc_report_meta <- function(model_type, family = NULL) {
  meta <- tc_model_meta(model_type)
  if (is.character(family)) {
    family <- tc_match_family(family, env = topenv())
  }
  fam <- family$family
  link <- family$link
  if (meta$scale == "family" && !is.null(link) &&
      !link %in% tc_supported_links) {
    cli::cli_abort("Unsupported link {.val {link}} for {.val {model_type}}.")
  }
  measure <- switch(meta$scale,
    difference = "Mean difference",
    hazard = "Hazard ratio",
    odds = "Cumulative odds ratio",
    rrr = "Relative risk ratio",
    family = if (is.null(link)) {
      NA_character_
    } else {
      switch(link,
        identity = "Mean difference",
        logit = "Odds ratio",
        log = if (fam %in% c("poisson", "quasipoisson")) {
          "Rate ratio"
        } else if (fam %in% c("binomial", "quasibinomial")) {
          "Risk ratio"
        } else {
          "Ratio of means"
        }
      )
    },
    cli::cli_abort("Unhandled effect scale for {.val {model_type}}.")
  )
  # A non-default link is named in the label, since it changes what the effect
  # means; the default is read from the stats family constructor of that name,
  # and a family with no such constructor gets no suffix.
  default_link <- if (!is.null(fam)) {
    ctor <- get0(fam, envir = asNamespace("stats"), mode = "function",
                 inherits = FALSE)
    if (!is.null(ctor)) tryCatch(ctor()$link, error = function(e) NULL)
  }
  link_suffix <- if (!is.null(default_link) && !is.null(link) &&
                     !identical(link, default_link)) {
    paste0(", ", link, " link")
  } else {
    ""
  }
  list(
    label = if (isTRUE(meta$show_family) && !is.null(fam)) {
      paste0(meta$label, " (", fam, link_suffix, ")")
    } else {
      meta$label
    },
    # Every effect but a difference is a ratio, reported exponentiated.
    exponentiate = if (is.na(measure)) NA else measure != "Mean difference",
    measure = measure
  )
}
