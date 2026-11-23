# Module 1: composition as an independent variable.
#
# `comp_iv()` models a `compD` composition as a predictor of some outcome. It
# resolves the interfaces (discrete arguments, `formula`, and the `moderator`
# and `curvature` shortcuts) into one spec, checks that spec against the
# engine's `tc_model_registry` row, and hands it to the shared `tc_iv_fit()`
# builder in fit_models.R, once for the canonical fit and then once per rotated
# pivot for the per-part refits. It
# adds the omnibus test of the whole composition's effect. There is no
# per-engine wrapper: every model class is a registry row. The object keeps the
# `timecoda_fit` parent class, so `predict_realloc()` /
# `plot_realloc()` consume it unchanged.

#' Model a composition as an independent variable
#'
#' Fits a model of an outcome on the pivot ILR coordinates of a time-use
#' composition (plus optional covariates), with the model class selected by
#' `model`. The composition comes entirely from a `compD` object created by
#' [declare_composition()]. The coordinates enter linearly, unless
#' `curvature = TRUE` adds their full quadratic or `formula` writes other terms.
#'
#' @param comp A `compD` object from [declare_composition()].
#' @param outcome Name of the outcome column in `comp$data_raw` (not used by
#'   `model = "cox"` or `"coxme"`, which use `time`/`event`). A model
#'   `formula` may also be passed here positionally (e.g.,
#'   `comp_iv(cd, y ~ ilr() * sex)`), the way [lm()]/[glm()] take a formula
#'   first; see `formula`.
#' @param model The fitting engine. `comp_iv()` wraps established backends, and
#'   `model` takes the backend's own name. Clustered and unclustered variants
#'   are separate keys, and `random` is required by the clustered ones and
#'   rejected by the rest. Backends marked Suggests are checked at call time.
#'   \describe{
#'     \item{`"lm"`, `"glm"`}{[stats::lm()], [stats::glm()].}
#'     \item{`"cox"`, `"coxme"`}{Cox proportional hazards models via
#'       [survival::coxph()] and, with a random intercept or other `(1 | g)`
#'       terms, [coxme::coxme()] (`random` required). Reports hazard ratios.}
#'     \item{`"lmer"`, `"glmer"`}{[lme4::lmer()], [lme4::glmer()]; `random`
#'       required.}
#'     \item{`"clm"`, `"clmm"`}{ordered-factor outcome via [ordinal::clm()] and
#'       [ordinal::clmm()] (`random` required). Reports cumulative odds ratios.}
#'     \item{`"mblogit"`}{unordered-factor outcome via [mclogit::mblogit()],
#'       optionally clustered through `random`. Reports one relative-risk ratio
#'       per non-reference outcome category, so [summary()][summary.comp_iv] and
#'       [predict_realloc()] gain a `response_level` column. For a two-level
#'       outcome, a logistic regression with `model = "glm"` and
#'       `family = "binomial"` reports the effect as an odds ratio. With
#'       `random`, the `convergence` element of the fit holds \pkg{mclogit}'s
#'       own convergence flag, which can read as converged although the inner
#'       iterations warned. A clustered fit has no log-likelihood.}
#'     \item{`"lmrob"`, `"rlmer"`}{robust linear models via
#'       [robustbase::lmrob()] and [robustlmm::rlmer()] (`random` required).}
#'   }
#'   `"cox"`, `"coxme"`, `"clm"`, `"clmm"` and `"mblogit"` support
#'   `predict_realloc(type = "relative")` only.
#' @param covariates Optional character vector of covariate terms. An
#'   interaction between the composition and a covariate (a term involving
#'   `ilr` coordinates) makes reallocation contrasts covariate-conditional;
#'   [predict_realloc()] then needs `at` or `by`.
#' @param family A GLM family (object or name, e.g., `"binomial"`), required for
#'   `model = "glm"` and `model = "glmer"` and rejected by the other engines. A
#'   family given by name is resolved in the calling environment, so a locally
#'   defined family function works. Its link must be `"identity"` (the effect is
#'   a mean difference), `"log"` (a rate ratio for Poisson, a risk ratio for
#'   binomial, a ratio of means otherwise) or `"logit"` (an odds ratio). Any
#'   other link, such as `binomial(link = "probit")` or the default inverse link
#'   of [stats::Gamma()], is an error.
#' @param random Random-effects specification. A character vector of grouping
#'   column names (each expanded to `(1 | name)`) and/or full `lme4` terms (any
#'   element containing `|` is used verbatim), such as
#'   `random = c("province", "school")`. Required by `"coxme"`, `"lmer"`,
#'   `"glmer"`, `"clmm"` and `"rlmer"`, optional for `"mblogit"`, and rejected
#'   by the unclustered engines. `"mblogit"` takes one grouping term, passed as
#'   its own `random` argument: a grouping name gives a random intercept
#'   (`~ 1 | school`), and a bar term keeps its left-hand side
#'   (`"(ilr1 | school)"` gives the random slope `~ ilr1 | school`).
#' @param time,event Survival outcome columns for `model = "cox"` and
#'   `"coxme"`: `time` positive follow-up times, `event` 1 = event,
#'   0 = censored.
#' @param formula Optional model formula, an alternative to
#'   `outcome`/`covariates`/`random`/`time`/`event`/`moderator` for models with
#'   further interactions and non-linear terms, fitted as written. It is
#'   mutually exclusive with the discrete arguments, `moderator` and
#'   `curvature`. The left-hand side is the outcome, or a
#'   `survival::Surv(time, event)` call for `model = "cox"` and `"coxme"`. On
#'   the right-hand side, the marker [ilr()] stands for the composition
#'   (expanded to `ilr1 + ... + ilr{D-1}`), and everything else passes through
#'   to the model function, so `poly()`, `splines::ns()`, `I()`, `a*b` and
#'   `(1 | g)` random terms all work.
#'
#'   Non-linear terms may enter the composition too (`ilr()^2 + I(ilr()^2)`,
#'   `poly(ilr(), 2)`), with or without its main effect; [ilr()] gives the
#'   expansion rules and the specifications that depend on the order of the
#'   parts. If `ilr()` is absent, the composition is added as a main effect,
#'   even beside coordinate names (`ilr1`, ...) typed by hand. When `ilr()` is
#'   present but a coordinate's main effect is not, `comp_iv()` warns that
#'   estimates may be biased. A formula fit, even `outcome ~ ilr() * x`, carries
#'   no `moderation` test and no default `by` for [predict_realloc()]; use
#'   `moderator` for those.
#'
#'   A variable the formula names must be a column of the `compD` data or a
#'   scalar in the formula's environment, such as a knots value. A longer
#'   object, such as a vector of the caller's, must be a column instead. A
#'   scalar is replaced by its value when the model is fitted, and the stored
#'   formula holds the value. The expressions in `covariates` follow the same
#'   rule, with or without `moderator` and `curvature`; their scalars are
#'   found in the environment `comp_iv()` is called from.
#' @param include_total Logical (default `FALSE`). When `TRUE`, fits the
#'   total-as-covariate model of Coenders et al. (2017), for compositions whose
#'   row totals vary (waking hours, dietary intake). The normalized log
#'   geometric mean of the zero-free parts on their original, pre-closure scale
#'   (the `parts_imputed` element of the `compD`),
#'   \eqn{t = \sqrt{D}\,\log\big[(x_1 \cdots x_D)^{1/D}\big]}, enters the model
#'   as the main effect `t_total`. The ILR coefficients then estimate the
#'   effects of proportional shifts at a constant total, and `t_total` the
#'   effect of scaling every part up proportionally. A (near-)constant total,
#'   such as the fixed 1440-minute days of `timeuse`, is an error. Reallocations
#'   hold `t_total` at `total_ref`, and the Wald test of `t_total` is
#'   `total_test`, shown by [summary()][summary.comp_iv]. The option combines
#'   with every engine, `covariates`, `formula`, `moderator` and `curvature`;
#'   `t_total` itself cannot be named in `covariates`, `formula`, `moderator`
#'   or `model_args`.
#' @param moderator Optional name of a single covariate to test as a moderator
#'   of the composition's effect. It adds the `ilr() * moderator` interaction,
#'   so the fit equals `formula = outcome ~ ilr() * moderator` (plus any
#'   `covariates` and `random` terms). Unlike that formula fit, the moderated
#'   fit stores the joint Wald test of the moderation in `moderation`. A factor
#'   or character moderator also becomes the default `by` of
#'   [predict_realloc()], which then gives one estimate per level. The test
#'   identifies the moderator's coefficients by the factor's own contrasts or
#'   those set with `options(contrasts =)`; a coding set through the backend's
#'   `contrasts` argument in `model_args` is not recognized. `outcome`, `time`
#'   and `event` follow the same rules as without a moderator. The moderator
#'   must not also appear in `covariates`, and `moderator` is mutually exclusive
#'   with `formula`.
#' @param curvature Logical (default `FALSE`). When `TRUE`, the composition
#'   enters as the full quadratic in the ILR coordinates,
#'   `ilr()^2 + I(ilr()^2)`: the D-1 main effects, every pairwise product and
#'   every square, (D-1)(D+2)/2 coefficients. The fit equals the one from
#'   `formula = outcome ~ ilr()^2 + I(ilr()^2)` with the same covariates and
#'   random terms, and every method reads it as that formula fit. The full
#'   quadratic is the same model in every ILR basis; a squares-only term or a
#'   spline, which changes with the order of the parts (see [ilr()]), is written
#'   through `formula` instead. With `moderator`, the interaction stays on the
#'   main effects (`ilr()^2 + I(ilr()^2) + ilr() * moderator`), so `moderation`
#'   tests the same coefficients as without `curvature`. The joint Wald test of
#'   the squares and products is stored as `curvature_test`. The argument is an
#'   error beside `formula`, which writes its own curvature terms, and when
#'   `covariates` names a coordinate (`ilr1`, ...) or the [ilr()] marker, both
#'   already held by the quadratic.
#' @param model_args A `list(...)` of further arguments passed straight to the
#'   fitting function of the chosen backend (e.g.,
#'   `model_args = list(weights = w, na.action = stats::na.omit)` for
#'   [stats::glm()], or `list(nAGQ = 0)` for [lme4::glmer()]), and the route to
#'   any backend argument `comp_iv()` does not expose itself. The arguments are
#'   captured unevaluated, so a bare column name works (`weights = w` refers to
#'   a column of the composition data, exactly as in a direct [stats::glm()]
#'   call) and the fitted object's stored call stays readable. A name that is
#'   not a column is found in the environment `comp_iv()` is called from (with
#'   or without `moderator`), or, with a `formula`, in that formula's
#'   environment.
#'
#'   An engine's own control object goes here too, such as
#'   `model_args = list(control = lme4::lmerControl(optimizer = "bobyqa"))` to
#'   retry a mixed fit with a different optimizer, or
#'   `list(control = stats::glm.control(maxit = 100))`. Each backend documents
#'   its own arguments, control functions included. An argument the engine
#'   already sets (`trace` for `"mblogit"`, `model` for `"cox"`, `Hess` for
#'   `"clmm"`) is an error. A `link` other than `"logit"` is an error for
#'   `"clm"` and `"clmm"`, whose effects are cumulative odds ratios.
#'
#' @section Total as a covariate:
#' The \eqn{\sqrt{D}} factor makes `t_total` the coordinate of \eqn{\log x}
#' along \eqn{\mathbf{1}/\sqrt{D}}, which completes the ILR basis to an
#' orthonormal one (Coenders et al. 2017, eq. 3.2). Any positive rescaling,
#' such as a \eqn{\sqrt{D-1}} convention, rescales its coefficient inversely
#' and leaves the fit, the predictions and `total_test` unchanged.
#'
#' @section Convergence:
#' Every engine is fitted with its backend warnings and messages captured (and
#' still re-emitted). The fit carries two additive elements: `fit_warnings`, a
#' tibble of the captured conditions (`source`, `type`, `message`; zero rows
#' when clean) whose `source` is `"canonical"` or a rotated part name, and
#' `convergence`, a one-row tibble (`converged`, `singular`, `n_warnings`).
#' [summary()][summary.comp_iv] prints the captured messages and, for a
#' clustered or iterative fit, a `Convergence: OK` line when clean;
#' [print()][print.timecoda_fit] flags the warning count. `"cox"` and `"coxme"`
#' fits are not probed (`converged` is `NA`), so they print no such line. An
#' `"lmer"` or `"glmer"` fit that fails to converge or is singular can be
#' refitted with a different optimizer through `model_args` (e.g.,
#' `model_args = list(control = lme4::lmerControl(optimizer = "bobyqa"))`).
#'
#' @section Degrees of freedom for clustered fits:
#' For an `"lmer"` fit, the per-part standard errors, intervals and p-values
#' that [summary()][summary.comp_iv] reports and the reallocations of
#' [predict_realloc()] use the same degrees-of-freedom method, chosen by the
#' rule `emmeans` applies. The default is the Kenward-Roger approximation
#' through \pkg{pbkrtest}. The Satterthwaite approximation through
#' \pkg{lmerTest} is used when `emmeans::emm_options(lmer.df = "satterthwaite")`
#' is set, or when \pkg{pbkrtest} is absent, disabled or the fit has more
#' observations than `emmeans::emm_options(pbkrtest.limit =)` (3000). Both
#' \pkg{pbkrtest} and \pkg{lmerTest} are optional, and when neither can be
#' used, or `lmer.df = "asymptotic"` is set, the degrees of freedom are
#' asymptotic (normal). `"glmer"`, `"cox"`, `"coxme"` and `"rlmer"` fits are
#' asymptotic in both places (\pkg{robustlmm} reports no degrees of freedom),
#' and `"lmrob"` fits use their residual degrees of freedom.
#'
#' @return An object of class `c("comp_iv", "timecoda_fit")`: the
#'   [`timecoda_fit`][comp_iv] metadata (fitted `model`, `model_type`, `parts`,
#'   `total`, `ref_comp`, `family`, ...) augmented with the elements below.
#'   `family` is the resolved GLM family object (class `"family"`, carrying its
#'   `link`) for `model = "glm"` and `"glmer"`, and `NULL` for every other
#'   engine. The backend object `model`, like each of `models`, carries its
#'   formula and fitting data in its call. R's own tools therefore re-fit it:
#'   [stats::update()], [stats::drop1()] and [stats::anova()] between two such
#'   fits work on it directly for every engine. [survival::cox.zph()] works on
#'   a plain Cox fit (`model = "cox"`). Current versions of \pkg{survival} do
#'   not support it on a `"coxme"` fit, and it should not be called on one.
#'   Printing that object itself (`print(fit$model)`, `summary(fit$model)`)
#'   shows the data inside the call; `print()` and `summary()` of the
#'   `comp_iv()` fit do not.
#'   \describe{
#'     \item{`omnibus`}{A one-row tibble (`term`, `statistic`, `df`, `p.value`):
#'       the joint Wald test that every coefficient of a term in the ILR
#'       coordinates alone is zero (i.e., that the composition as a whole has no
#'       effect). The statistic is referred to a chi-squared distribution on as
#'       many degrees of freedom as there are such coefficients. These are D-1
#'       for the main effects, plus one per square, product or basis column of a
#'       non-linear term ((D-1)(D+2)/2 for the full quadratic
#'       `ilr()^2 + I(ilr()^2)`), and one set per non-reference category for a
#'       multinomial model. A formula with no term in the coordinates alone
#'       (`ilr():sex`) tests every coefficient crossing a coordinate with a
#'       covariate instead. Aliased (`NA`) coefficients, such as a linear term
#'       repeated by `ilr() + poly(ilr(), 2)`, are left out of the test and of
#'       its degrees of freedom. The test is `NA` on 0 degrees of freedom when
#'       no fixed-effect coefficient carries the composition, for example when
#'       it enters only through a random slope
#'       (`zbmi ~ age + (ilr() | school)`). The test is also `NA`, with a
#'       warning, when every composition coefficient is aliased (a composition
#'       that is the same on every row, or collinear with another term). For a
#'       linear model the statistic equals its degrees of freedom times the F
#'       statistic of the same hypothesis. The chi-squared reference ignores
#'       the uncertainty in the residual variance, so in small samples the test
#'       is slightly anti-conservative, with a p-value a little smaller than
#'       that of the exact F test.}
#'     \item{`comp`}{The `compD` object.}
#'     \item{`models`}{A named list (one element per part) of the rotated
#'       fits: `models[[part]]` is the model refitted with that part as the
#'       first pivot coordinate, so its `ilr1` term is the "part vs remaining"
#'       effect summarized by [summary()][summary.comp_iv]. The refits leave the
#'       canonical `model` unchanged, and no prediction path reads them. For a
#'       mixed model, one refit per part can be slow.}
#'     \item{`moderator`}{The moderator column name (or `NULL`).}
#'     \item{`moderation`}{When `moderator` is set, a one-row tibble with the
#'       joint Wald test that every composition-by-moderator coefficient is zero
#'       (`H0`: no moderation). Its degrees of freedom are (L-1)(D-1) for an
#'       L-level factor moderator and D-1 for a numeric one, times K-1 for a
#'       multinomial fit with K outcome categories. `NULL` otherwise.}
#'     \item{`curvature_test`}{A one-row tibble when the terms in the ILR
#'       coordinates alone are the full quadratic: a square written `I(ilrK^2)`
#'       for every coordinate, a product `ilrJ:ilrK` for every pair, and no
#'       other term above degree one. It is `NULL` otherwise, including for
#'       `poly()` and spline terms and for a square or product written in
#'       another form (such as `I(ilr1 * ilr2)`). The tibble holds the joint
#'       Wald test that every square and product coefficient is zero (`H0`: the
#'       effect is linear in the coordinates). The degrees of freedom are
#'       (D-1)D/2, times K-1 for a multinomial fit with K outcome categories,
#'       less any aliased coefficient. The test is the same whether the
#'       quadratic comes from `curvature = TRUE` or from a formula, with or
#'       without a moderator, and in every ILR basis.}
#'     \item{`random`}{The resolved `random` specification (character vector of
#'       grouping terms), or `NULL` for a non-clustered fit.}
#'     \item{`fit_warnings`}{A tibble of the backend warnings/messages captured
#'       across the canonical and rotated fits (columns `source`, `type`,
#'       `message`; zero rows when every fit was clean). See the *Convergence*
#'       section.}
#'     \item{`convergence`}{A one-row tibble (`converged`, `singular`,
#'       `n_warnings`) probing the canonical fit's convergence state.}
#'     \item{`total_test`}{When `include_total = TRUE`, a one-row tibble: the
#'       Wald test of the `t_total` coefficient (`H0`: the total carries no
#'       effect beyond the composition), on 1 degree of freedom, or K-1 for a
#'       multinomial fit with K outcome categories. `NULL` otherwise.}
#'     \item{`total_ref`}{When `include_total = TRUE`, the mean of `t_total`,
#'       which is `t_total` at the geometric-mean composition of the
#'       original-scale parts; every reallocation holds `t_total` at this
#'       value. `NULL` otherwise.}
#'   }
#'
#' @references
#' Coenders, G., Martín-Fernández, J. A., & Ferrer-Rosell, B. (2017).
#'   When relative and absolute information matter: Compositional predictor with
#'   a total in generalized linear models. *Statistical Modelling*, 17(6),
#'   494-512. \doi{10.1177/1471082X17710398}
#'
#' Cox, D. R. (1972). Regression models and life-tables. *Journal of the Royal
#'   Statistical Society: Series B (Methodological)*, 34(2), 187-202.
#'   \doi{10.1111/j.2517-6161.1972.tb00899.x}
#'
#' Kenward, M. G., & Roger, J. H. (1997). Small sample inference for fixed
#'   effects from restricted maximum likelihood. *Biometrics*, 53(3), 983-997.
#'   \doi{10.2307/2533558}
#'
#' McCullagh, P. (1980). Regression models for ordinal data. *Journal of the
#'   Royal Statistical Society: Series B (Methodological)*, 42(2), 109-127.
#'   \doi{10.1111/j.2517-6161.1980.tb01109.x}
#'
#' @examples
#' cd <- declare_composition(timeuse, c("MVPA", "LPA", "SB", "SPT"),
#'                           closure = 1440)
#' m <- comp_iv(cd, outcome = "obesity", model = "glm", family = "binomial")
#' m$omnibus
#'
#' # The full quadratic in the coordinates, a predicted outcome that can curve
#' curved <- comp_iv(cd, outcome = "zbmi", model = "lm", curvature = TRUE)
#' curved$omnibus
#'
#' @seealso [declare_composition()] to build the `compD` input;
#'   [summary.comp_iv()] and [plot.comp_iv()] to summarize and visualize the
#'   fit; [predict_realloc()] for reallocation predictions; [fit_check()] for
#'   residual diagnostics; [ilr()] for the formula-interface composition
#'   marker.
#' @family composition as predictor
#' @export
comp_iv <- function(comp, outcome = NULL,
                    model = c("lm", "glm", "cox", "coxme", "lmer", "glmer",
                              "clm", "clmm", "mblogit", "lmrob", "rlmer"),
                    covariates = NULL, family = NULL, random = NULL,
                    time = NULL, event = NULL, formula = NULL, moderator = NULL,
                    curvature = FALSE, include_total = FALSE,
                    model_args = list()) {
  if (!inherits(comp, "compD")) {
    cli::cli_abort("{.arg comp} must be a {.cls compD} object from {.fn declare_composition}.")
  }
  model <- match.arg(model)
  user_env <- rlang::caller_env()
  # The family as typed ("binomial", binomial(link = "log"), a variable name),
  # captured before it is resolved below, for the fitted object's stored call.
  family_expr <- rlang::enexpr(family)
  # Captured unevaluated: `model_args = list(weights = w)` must reach the
  # backend as the expression `w`, so the backend resolves it in the model frame
  # the way it would in a direct call. Evaluating here would both fail on a bare
  # column name and inline the value into the fit's stored call.
  model_args <- rlang::enexpr(model_args)
  model_args <- if (rlang::is_call(model_args, "list")) {
    rlang::call_args(model_args)
  } else {
    rlang::eval_tidy(model_args, env = user_env)
  }
  if (!is.list(model_args)) {
    cli::cli_abort("{.arg model_args} must be a list of backend arguments.")
  }
  if (length(model_args) > 0 && !rlang::is_named(model_args)) {
    cli::cli_abort(
      c("Every element of {.arg model_args} must be named.",
        "i" = "Name them as the backend does (e.g. {.code model_args = list(weights = w)}).")
    )
  }
  # Every argument rule is a registry field, read here before anything is fit.
  meta <- tc_model_meta(model)
  if (!is.null(family) && !isTRUE(meta$allows_family)) {
    cli::cli_abort(
      c("{.arg family} does not apply to {.code model = \"{model}\"}.",
        "i" = "It is used by {.code \"glm\"} and {.code \"glmer\"}.")
    )
  }
  if (!is.null(meta$pkg)) {
    rlang::check_installed(meta$pkg,
                           reason = paste0("for `model = \"", model, "\"`."))
  }
  fixed <- intersect(names(model_args), names(meta$fixed_args))
  if (length(fixed) > 0) {
    cli::cli_abort(c(
      "{.arg model_args} names {.arg {fixed}}, which {.code model = \"{model}\"} sets itself.",
      "i" = "Leave {.arg {fixed}} out of {.arg model_args}."
    ))
  }
  # An ordinal effect is reported as a cumulative odds ratio, which only the
  # logit link gives; `model_args` holds expressions, so the link is evaluated.
  if (identical(meta$scale, "odds") && !is.null(model_args[["link"]])) {
    link <- rlang::eval_tidy(model_args[["link"]], env = user_env)
    if (!identical(link, "logit")) {
      cli::cli_abort(c(
        "{.code model = \"{model}\"} reports cumulative odds ratios, which need the logit link; got {.val {link}}.",
        "i" = "Leave {.arg link} out of {.arg model_args}, or fit the model directly with {.fn {deparse(meta$fitter)}} for another link."
      ))
    }
  }
  parts <- comp$parts
  total <- comp$closure

  # Allow the formula to be passed positionally as the second argument
  # (comp_iv(cd, y ~ x, ...)), the way lm()/glm() take a formula first.
  if (inherits(outcome, "formula")) {
    if (!is.null(formula)) {
      cli::cli_abort("Supply the model formula once (positionally or via {.arg formula}).")
    }
    formula <- outcome
    outcome <- NULL
  }

  if (!rlang::is_bool(include_total)) {
    cli::cli_abort("{.arg include_total} must be {.code TRUE} or {.code FALSE}.")
  }
  # The total covariate takes the column name `t_total`, so a model that names
  # it anywhere would collide with the term appended below. A covariate string
  # that does not parse is kept as typed, so the spec resolution below refuses
  # it with the same message as without `include_total`.
  if (include_total) {
    named <- c(
      if (!is.null(formula)) all.vars(formula),
      if (length(covariates) > 0) {
        tryCatch(all.vars(stats::reformulate(covariates)),
                 error = function(e) covariates)
      },
      moderator,
      unlist(purrr::map(model_args, all.vars))
    )
    if ("t_total" %in% named) {
      cli::cli_abort(c(
        "{.arg include_total} adds the covariate {.var t_total}, so the model cannot also name {.var t_total}.",
        "i" = "Leave it out of {.arg formula}, {.arg covariates}, {.arg moderator} and {.arg model_args}; the total enters as a main effect on its own."
      ))
    }
  }

  spec <- tc_iv_resolve_spec(
    comp, model, outcome, covariates, random, time, event, formula, moderator,
    curvature, env = user_env
  )
  # A part enters the model through the coordinates only. Named again, as a
  # covariate, the moderator, a formula term or a `model_args` value such as
  # `weights`, it would be bound into the model frame a second time.
  part_vars <- intersect(
    c(spec$need, unlist(purrr::map(model_args, all.vars))), parts
  )
  if (length(part_vars) > 0) {
    cli::cli_abort(c(
      "The model must not name a composition part as a variable; {.val {part_vars}} {?is a part/are parts}.",
      "i" = "The parts enter the model through their coordinates; leave them out of {.arg covariates}, {.arg moderator}, {.arg formula} and {.arg model_args}."
    ))
  }

  # The clustering rules are checked on the resolved spec, because the formula
  # interface carries its random terms in the formula rather than in `random`.
  if (length(spec$random) > 0 && !isTRUE(meta$allows_random)) {
    clustered <- purrr::keep(names(tc_model_registry), function(k) {
      isTRUE(tc_model_meta(k)$allows_random)
    })
    cli::cli_abort(
      c("{.code model = \"{model}\"} does not accept random-effect terms.",
        "i" = "Use a clustered engine instead: {.val {clustered}}.")
    )
  }
  if (length(spec$random) == 0 && isTRUE(meta$requires_random)) {
    cli::cli_abort(
      c("{.code model = \"{model}\"} is a clustered engine and needs a grouping term.",
        "i" = "Supply {.arg random} (e.g. {.code random = \"school\"}) or a {.code (1 | group)} term in {.arg formula}.")
    )
  }
  # An engine that takes its grouping through its own `random =` argument
  # (multinomial) fits one grouping factor, so the spec may carry one term;
  # terms are counted, so one string holding two bar terms is caught too. A
  # `||` term is checked as typed, before the formula interface splits it.
  if (identical(meta$random_style, "arg")) {
    typed <- c(random, if (!is.null(formula)) deparse1(formula))
    doubled <- typed[grepl("||", typed, fixed = TRUE)]
    if (length(doubled) > 0) {
      correlated <- gsub("||", "|", doubled, fixed = TRUE)
      cli::cli_abort(c(
        "{.code model = \"{model}\"} fits correlated random effects only; {.code ||} is not supported by {.fn {deparse(meta$fitter)}}.",
        "i" = "Write the term with {.code |}: {.val {correlated}}."
      ))
    }
    n_terms <- sum(purrr::map_int(spec$random, function(r) {
      length(attr(stats::terms(stats::as.formula(paste("~", r))), "term.labels"))
    }))
    if (n_terms > 1) {
      cli::cli_abort(
        c("{.code model = \"{model}\"} accepts one grouping term, not {n_terms}.",
          "x" = "Got {.val {spec$random}}.")
      )
    }
  }

  # Resolve the family once, here, in the caller's environment, so a locally
  # defined family function works; this is also where a missing one is caught,
  # after the clustering rules, whose message names the right engine.
  if (isTRUE(meta$allows_family)) {
    family <- tc_match_family(family, env = user_env)
    if (!family$link %in% tc_supported_links) {
      cli::cli_abort(
        c("The {.val {family$link}} link of {.arg family} {.val {family$family}} is not supported.",
          "i" = "Supported links: {.val {tc_supported_links}}.",
          "i" = "The per-part effect must read as a mean difference (identity), a rate, risk or mean ratio (log), or an odds ratio (logit).")
      )
    }
  }

  # Assemble the non-part columns needed from the compD data.
  # `model_args` may name data columns the way a direct backend call would
  # (`weights = wt`), so its free variables join the columns carried into the
  # model frame, where the backend finds them.
  need <- unique(c(spec$need, unlist(purrr::map(model_args, all.vars))))
  need <- need[need %in% names(comp$data_raw)]
  df <- dplyr::bind_cols(tc_comp_df(comp), comp$data_raw[, need, drop = FALSE])

  # The Coenders total covariate (their eq. 3.2): sqrt(D) times the mean log of
  # the pre-closure parts, which declare_composition() leaves zero-free and
  # complete. It is an ordinary main-effect covariate, so every engine, rotated
  # refit and formula carries it.
  if (include_total) {
    raw <- as.matrix(comp$parts_imputed)
    totals <- rowSums(raw)
    # A constant total is collinear with the intercept; one row has no SD (NA).
    cv <- stats::sd(totals) / mean(totals)
    if (!is.finite(cv) || cv < 1e-8) {
      problem <- if (is.finite(cv)) {
        "the observed totals are constant (coefficient of variation {signif(cv, 3)})"
      } else {
        "the composition has fewer than two observations"
      }
      cli::cli_abort(
        c(paste0("{.code include_total = TRUE} needs a variable row total, but ", problem, "."),
          "i" = "A constant total is collinear with the intercept and carries no information; total-adjustment is for variable-total compositions (e.g. waking-hours or dietary intake).",
          "i" = "For a closed composition (e.g. the fixed 1440-minute {.val timeuse} data) leave {.arg include_total} at its default {.code FALSE}.")
      )
    }
    df$t_total <- sqrt(ncol(raw)) * rowMeans(log(raw))
    spec$covariates <- c(spec$covariates, "t_total")
    if (!is.null(spec$model_formula)) {
      spec$model_formula <- rlang::new_formula(
        rlang::f_lhs(spec$model_formula),
        call("+", rlang::f_rhs(spec$model_formula), quote(t_total)),
        env = rlang::f_env(spec$model_formula)
      )
    }
  }

  # Everything that does not depend on the pivot is done once, here, outside
  # the D fits below.
  df <- tc_iv_check_response(df, spec, model)
  check_composition(df, parts, total = total)
  ref_comp <- geom_mean_comp(df, parts, total = total)
  # All D pivot frames from one pivotCoord() call: pivot k's coordinates are
  # pivot 1's times a fixed rotation (by D and k only), solved from probe rows.
  d <- length(parts)
  z1 <- tc_pivot_coords(df[parts], pivot = 1L)
  probe <- df[seq_len(min(nrow(df), d + 2L)), parts, drop = FALSE]
  probe1 <- tc_pivot_coords(probe, pivot = 1L)
  degenerate <- qr(probe1)$rank < d - 1L
  frames <- purrr::map(seq_len(d), function(k) {
    if (k > 1L && degenerate) {
      return(ilr_transform(df, parts, pivot = k))
    }
    rot <- if (k == 1L) {
      diag(d - 1L)
    } else {
      qr.solve(probe1, tc_pivot_coords(probe, pivot = k))
    }
    df[colnames(z1)] <- as.data.frame(z1 %*% rot)
    df
  })
  fit <- tc_iv_fit(frames[[1L]], model, spec, parts, total, ref_comp, family,
                   family_expr = family_expr, pivot = 1L, env = user_env,
                   model_args = model_args, source = "canonical")
  # A backend that gave up can return placeholder coefficients and no
  # covariance (lmrob's covariance is then NA), which no test can read.
  if (isFALSE(fit$convergence$converged) &&
      all(is.na(tryCatch(as.matrix(stats::vcov(fit$model)),
                         error = function(e) NA)))) {
    cli::cli_abort(c(
      "The {.code model = \"{model}\"} fit did not converge and returned no covariance matrix.",
      "i" = "Adjust the backend's control settings through {.arg model_args}; see the control options of {.fn {deparse(meta$fitter)}}."
    ))
  }

  # Omnibus test of the whole composition effect: a joint Wald over the
  # coefficients of every term in the coordinates alone (main effects, squares,
  # products, basis columns), one set or one per outcome category of a
  # per-category model. A formula with no such term (`ilr():sex`) tests every
  # coefficient crossing a coordinate with a covariate instead.
  coords <- tc_coord_names(length(parts))
  b <- tc_model_coefs(fit$model, model)$b
  b_names <- names(b)
  comp_terms <- tc_composition_terms(c(coords, spec$covariates), coords)
  crossed_names <- tc_ilr_coef_names(b_names, model, comp_terms$crossed,
                                     levels = TRUE)
  omnibus_names <- tc_ilr_coef_names(b_names, model, comp_terms$only)
  if (length(omnibus_names) == 0) {
    omnibus_names <- crossed_names
  }
  fit$omnibus <- tc_joint_ilr_test(fit$model, model, omnibus_names)
  # Every composition coefficient aliased (NA): the composition is constant
  # across rows or collinear with another term. A backend that drops such
  # columns instead (lme4) leaves no name to test and says so itself.
  comp_names <- union(omnibus_names, crossed_names)
  if (length(comp_names) > 0 && all(is.na(b[comp_names]))) {
    cli::cli_warn(c(
      "No coefficient of the composition can be estimated in this fit.",
      "i" = "The composition may be the same on every row, or collinear with another term of the model.",
      "i" = "The fit has no composition effect to summarize or reallocate, and its omnibus test is {.code NA}."
    ))
  }
  # The curvature test, stored when the terms in the coordinates alone are the
  # full quadratic: a square I(ilrK^2) for every coordinate, a product
  # ilrJ:ilrK for every pair (keyed by its sorted pair) and nothing else above
  # degree one. It tests the squares and products jointly.
  higher <- setdiff(comp_terms$only, coords)
  keys <- purrr::map_chr(higher, function(label) {
    e <- str2lang(label)
    cmp <- strsplit(label, "(?<!:):(?!:)", perl = TRUE)[[1]]
    if (rlang::is_call(e, "I") && rlang::is_call(e[[2]], "^") &&
        is.name(e[[2]][[2]]) && is.numeric(e[[2]][[3]]) &&
        isTRUE(e[[2]][[3]] == 2)) {
      paste(rep(as.character(e[[2]][[2]]), 2), collapse = ":")
    } else if (length(cmp) == 2 && all(cmp %in% coords) && cmp[1] != cmp[2]) {
      paste(sort(cmp), collapse = ":")
    } else {
      NA_character_
    }
  })
  quadratic <- c(
    paste0(coords, ":", coords),
    if (length(coords) > 1) {
      utils::combn(coords, 2, function(p) paste(sort(p), collapse = ":"))
    }
  )
  if (!anyNA(keys) && setequal(keys, quadratic)) {
    fit$curvature_test <- tc_joint_ilr_test(
      fit$model, model, tc_ilr_coef_names(b_names, model, higher),
      term = "curvature"
    )
  }
  fit$comp <- comp
  # The formula path resolves to character(0) when there is no random term;
  # the documented contract is NULL.
  fit$random <- if (length(spec$random) > 0) spec$random else NULL
  fit$moderator <- spec$moderator
  fit$moderation <- if (!is.null(spec$moderator)) {
    tc_moderation_test(fit$model, fit$model_type, spec$moderator, parts)
  } else {
    NULL
  }
  # The Wald test of the t_total coefficient (one "<level>~t_total" per category
  # of a per-category model, tested jointly), and the value every emmeans grid
  # pins it at, so a reallocation is a constant-total contrast.
  if (include_total) {
    fit$total_test <- tc_joint_ilr_test(
      fit$model, model, tc_ilr_coef_names(b_names, model, "t_total"),
      term = "total"
    )
    fit$total_ref <- mean(df$t_total)
  }
  # Per-part "vs remaining" effects need one model per pivot. The canonical
  # pivot-1 fit above is reused as models[[parts[1]]]; the rest are refit. The
  # omnibus / prediction paths keep using the canonical `fit$model` only.
  # `spec` and `model_args` are reused at every pivot (the formula's column
  # names are pivot-stable); each refit is labelled by its part name so its
  # captured conditions stay distinguishable from the canonical fit's.
  rotated <- purrr::map(seq_along(parts)[-1L], function(k) {
    tc_iv_fit(frames[[k]], model, spec, parts, total, ref_comp, family,
              family_expr = family_expr, pivot = k, env = user_env,
              model_args = model_args, source = parts[k])
  })
  fit$models <- rlang::set_names(
    c(list(fit$model), purrr::map(rotated, "model")), parts
  )
  # Roll up the canonical fit's captured conditions with those from the rotated
  # refits so summary() sees every warning.
  fit$fit_warnings <- dplyr::bind_rows(
    fit$fit_warnings,
    purrr::map(rotated, "fit_warnings") |> purrr::list_rbind()
  )
  fit$convergence$n_warnings <- nrow(fit$fit_warnings)
  class(fit) <- c("comp_iv", class(fit))
  fit
}

# Check the response columns once, before any fit: the survival pair, or the
# single outcome and, where the engine constrains it, its type. Returns `data`,
# with the multinomial engine's outcome coerced to a factor.
tc_iv_check_response <- function(data, spec, model, call = rlang::caller_env()) {
  meta <- tc_model_meta(model)
  surv <- identical(meta$response, "surv")
  outcome <- spec$outcome

  # The formula interface has already checked its own left-hand side.
  if (is.null(spec$model_formula)) {
    if (surv) {
      # `outcome` has no role beside the survival pair, so it is refused
      # rather than silently dropped.
      if (!is.null(outcome)) {
        cli::cli_abort(
          c("{.arg outcome} does not apply to {.code model = \"{model}\"}.",
            "i" = "A survival engine takes its response from {.arg time} and {.arg event}."),
          call = call
        )
      }
      # Both columns are checked before either is read, so one message names
      # everything that is missing.
      bad <- purrr::keep(c("time", "event"), function(arg) {
        nm <- spec[[arg]]
        !is.character(nm) || length(nm) != 1 || !nm %in% names(data)
      })
      if (length(bad) > 0) {
        verb <- if (length(bad) == 1) "name" else "each name"
        cli::cli_abort(
          "{.arg {bad}} must {verb} a single column of the {.cls compD} data.",
          call = call
        )
      }
      if (!is.numeric(data[[spec$time]]) ||
          any(data[[spec$time]] <= 0, na.rm = TRUE)) {
        cli::cli_abort(
          "{.arg time} column {.val {spec$time}} must be positive numeric follow-up times.",
          call = call
        )
      }
      if (!all(data[[spec$event]] %in% c(0, 1))) {
        cli::cli_abort(
          "{.arg event} column {.val {spec$event}} must contain only 0 (censored) and 1 (event).",
          call = call
        )
      }
    } else if (is.null(outcome)) {
      cli::cli_abort("{.arg outcome} is required for this model.", call = call)
    }
  }
  if (!surv && !is.null(outcome) &&
      (!is.character(outcome) || length(outcome) != 1 ||
       !outcome %in% names(data))) {
    cli::cli_abort(
      "{.arg outcome} must name a single column of the {.cls compD} data.",
      call = call
    )
  }

  if (identical(meta$outcome_check, "ordered") &&
      !is.ordered(data[[outcome]])) {
    cli::cli_abort(
      c("Outcome {.val {outcome}} must be an ordered factor for {.code model = \"{model}\"}.",
        "i" = "Make it an ordered factor in the data given to {.fn declare_composition} (e.g. {.code d${outcome} <- ordered(d${outcome})})."),
      call = call
    )
  }
  if (identical(meta$outcome_check, "unordered")) {
    if (is.ordered(data[[outcome]])) {
      cli::cli_abort(
        c("Outcome {.val {outcome}} is an ordered factor; use {.code model = \"clm\"}.",
          "i" = "{.code model = \"mblogit\"} is for unordered categorical outcomes."),
        call = call
      )
    }
    # A numeric column would become one category per distinct value.
    if (is.numeric(data[[outcome]])) {
      cli::cli_abort(c(
        "Outcome {.val {outcome}} is numeric; {.code model = \"{model}\"} is for unordered categorical outcomes.",
        "i" = "Give it as a factor, character or logical column of the data passed to {.fn declare_composition}."
      ), call = call)
    }
    data[[outcome]] <- as.factor(data[[outcome]])
  }
  data
}

# Resolve the discrete-argument / formula / moderator / curvature interfaces
# into one internal spec: list(outcome, time, event, covariates, random,
# model_formula, moderator, need). On the discrete-argument path
# `model_formula` is NULL, and tc_iv_fit() builds the formula from
# `covariates`; `env` is comp_iv()'s caller.
tc_iv_resolve_spec <- function(comp, model, outcome, covariates, random, time,
                               event, formula, moderator, curvature,
                               env = rlang::caller_env(),
                               call = rlang::caller_env()) {
  has_discrete <- !is.null(outcome) || !is.null(covariates) ||
    !is.null(random) || !is.null(time) || !is.null(event)

  if (!rlang::is_bool(curvature)) {
    cli::cli_abort("{.arg curvature} must be {.code TRUE} or {.code FALSE}.",
                   call = call)
  }
  if (!is.null(formula) && !is.null(moderator)) {
    cli::cli_abort("Supply either {.arg formula} or {.arg moderator}, not both.",
                   call = call)
  }
  if (!is.null(formula) && curvature) {
    cli::cli_abort(
      "{.arg curvature} applies to the discrete arguments only; write the curvature terms in {.arg formula} ({.code ilr()^2 + I(ilr()^2)}).",
      call = call
    )
  }
  # The curvature and moderator shortcuts paste the covariates into a formula,
  # so a covariate that does not parse is refused first, with the message the
  # discrete path gives it.
  if ((curvature || !is.null(moderator)) && length(covariates) > 0) {
    rhs <- tryCatch(stats::reformulate(covariates), error = function(e) NULL)
    if (is.null(rhs)) {
      unparsed <- purrr::discard(covariates, function(t) {
        rlang::is_formula(tryCatch(stats::as.formula(paste("~", t)),
                                   error = function(e) NULL))
      })
      # A non-character value fails reformulate() although each element parses.
      if (length(unparsed) == 0) unparsed <- covariates
      cli::cli_abort("Column{?s} {.val {unparsed}} not found in the {.cls compD} data.",
                     call = call)
    }
  }
  # The quadratic holds every coordinate term, so a covariate naming one, or
  # the `ilr()` marker that expands to them, would repeat it.
  if (curvature && length(covariates) > 0) {
    offending <- c(intersect(all.vars(rhs), tc_coord_names(length(comp$parts))),
                   if (tc_expr_has_ilr(rhs)) "ilr()")
    if (length(offending) > 0) {
      cli::cli_abort(
        "{.arg curvature} writes every coordinate term itself, so {.arg covariates} cannot also name a coordinate or the {.fn ilr} marker ({.val {offending}}).",
        call = call
      )
    }
  }
  if (!is.null(formula)) {
    if (has_discrete) {
      cli::cli_abort(
        c("{.arg formula} cannot be combined with the discrete model arguments.",
          "i" = "Use either {.arg formula}, or {.arg outcome}/{.arg covariates}/{.arg random}/{.arg time}/{.arg event}."),
        call = call
      )
    }
    return(tc_iv_parse_formula(formula, comp, model, call = call))
  }
  if (!is.null(moderator)) {
    if (!is.character(moderator) || length(moderator) != 1) {
      cli::cli_abort("{.arg moderator} must be a single column name.", call = call)
    }
    if (!moderator %in% names(comp$data_raw)) {
      cli::cli_abort("{.arg moderator} {.val {moderator}} is not a column of the {.cls compD} data.",
                     call = call)
    }
    if (moderator %in% covariates) {
      cli::cli_abort("{.arg moderator} {.val {moderator}} must not also appear in {.arg covariates}; it is adjusted for automatically.",
                     call = call)
    }
  }
  # A covariate or random term may be an expression ("ilr1:sex",
  # "(1 | school)"); its variables name the columns to carry from the compD
  # data into the model frame, and a plain column name is its own.
  term_vars <- function(t) {
    v <- tryCatch(all.vars(stats::as.formula(paste("~", t))),
                  error = function(e) t)
    if (length(v) == 0) t else v
  }
  # As on the formula path, a covariate's variable that is not a data column
  # but a scalar in the caller's environment is replaced by its value, before
  # either the discrete path or the shortcuts below use the covariates. The
  # value replaces the name wherever it is used as a value, so with `c <- 3`,
  # "I(time * c + c(0))" becomes "I(time * 3 + c(0))". Anything else is
  # refused (the shortcuts by the formula parser), so every variable stays
  # row-aligned with the compD.
  columns <- c(tc_coord_names(length(comp$parts)), names(comp$data_raw))
  cov_vars <- purrr::map(covariates, term_vars)
  cv <- unlist(cov_vars)
  scalars <- tc_caller_scalars(setdiff(cv, columns), env)
  # Recursive, so a name is replaced at any depth; a call's function slot is
  # never visited. Loops, not purrr::map(), so an empty argument (`x[, 1]`)
  # is never passed on as a value.
  subst_values <- function(e) {
    if (is.call(e)) {
      for (i in seq_along(e)[-1]) {
        if (!identical(e[[i]], quote(expr = ))) e[[i]] <- subst_values(e[[i]])
      }
      return(e)
    }
    if (is.name(e) && as.character(e) %in% scalars$present) {
      return(scalars$values[[as.character(e)]])
    }
    e
  }
  scalar_covariates <- covariates
  if (length(scalars$present) > 0) {
    scalar_covariates <- purrr::map2_chr(covariates, cov_vars, function(t, v) {
      if (!any(v %in% scalars$present)) {
        return(t)
      }
      deparse1(subst_values(str2lang(t)))
    })
  }

  if (curvature || !is.null(moderator)) {
    # Sugar: build `outcome ~ ilr()^2 + I(ilr()^2) + ilr() * moderator
    # (+ covariates + random)`, with only the terms of the options given.
    # It is routed through the formula path, so the fit is identical to the
    # equivalent explicit formula.
    # The response arguments follow the discrete path's rules (a survival
    # engine refuses `outcome`), checked before the formula is built from them.
    tc_iv_check_response(
      comp$data_raw, list(outcome = outcome, time = time, event = event), model,
      call = call
    )
    rhs_txt <- paste(c(
      if (curvature) "ilr()^2 + I(ilr()^2)",
      if (!is.null(moderator)) paste0("ilr() * ", moderator),
      scalar_covariates,
      if (!is.null(random)) tc_random_terms(random)
    ), collapse = " + ")
    lhs_txt <- if (identical(tc_model_meta(model)$response, "surv")) {
      sprintf("survival::Surv(%s, %s)", time, event)
    } else {
      outcome
    }
    # A small environment parented on the package namespace, as the discrete
    # path gives its formula (this frame holds the compD); `.tc_synth` lets
    # tc_iv_fit() copy caller-local `model_args` objects into it.
    synth <- stats::as.formula(
      paste(lhs_txt, "~", rhs_txt),
      env = rlang::new_environment(list(.tc_synth = TRUE),
                                   parent = topenv(environment()))
    )
    spec <- tc_iv_parse_formula(synth, comp, model, call = call)
    if (!is.null(moderator)) {
      spec$moderator <- moderator
    }
    return(spec)
  }

  # Discrete-argument path. A bar term's variables must be columns; a plain
  # grouping name is checked when the model is built.
  rand_vars <- purrr::map(random, term_vars)
  absent <- setdiff(
    c(cv, unlist(rand_vars[grepl("|", random, fixed = TRUE)])),
    c(columns, scalars$present)
  )
  if (length(absent) > 0) {
    cli::cli_abort("Column{?s} {.val {absent}} not found in the {.cls compD} data.",
                   call = call)
  }
  need <- unique(c(outcome, time, event, setdiff(cv, scalars$present),
                   purrr::list_c(rand_vars)))
  list(outcome = outcome, time = time, event = event,
       covariates = scalar_covariates, random = random, model_formula = NULL,
       moderator = NULL, need = need)
}
