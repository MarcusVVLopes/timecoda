# Model a composition as an independent variable

Fits a model of an outcome on the pivot ILR coordinates of a time-use
composition (plus optional covariates), with the model class selected by
`model`. The composition comes entirely from a `compD` object created by
[`declare_composition()`](https://marcusvvlopes.github.io/timecoda/reference/declare_composition.md).
The coordinates enter linearly, unless `curvature = TRUE` adds their
full quadratic or `formula` writes other terms.

## Usage

``` r
comp_iv(
  comp,
  outcome = NULL,
  model = c("lm", "glm", "cox", "coxme", "lmer", "glmer", "clm", "clmm", "mblogit",
    "lmrob", "rlmer"),
  covariates = NULL,
  family = NULL,
  random = NULL,
  time = NULL,
  event = NULL,
  formula = NULL,
  moderator = NULL,
  curvature = FALSE,
  include_total = FALSE,
  model_args = list()
)
```

## Arguments

- comp:

  A `compD` object from
  [`declare_composition()`](https://marcusvvlopes.github.io/timecoda/reference/declare_composition.md).

- outcome:

  Name of the outcome column in `comp$data_raw` (not used by
  `model = "cox"` or `"coxme"`, which use `time`/`event`). A model
  `formula` may also be passed here positionally (e.g.,
  `comp_iv(cd, y ~ ilr() * sex)`), the way
  [`lm()`](https://rdrr.io/r/stats/lm.html)/[`glm()`](https://rdrr.io/r/stats/glm.html)
  take a formula first; see `formula`.

- model:

  The fitting engine. `comp_iv()` wraps established backends, and
  `model` takes the backend's own name. Clustered and unclustered
  variants are separate keys, and `random` is required by the clustered
  ones and rejected by the rest. Backends marked Suggests are checked at
  call time.

  `"lm"`, `"glm"`

  :   [`stats::lm()`](https://rdrr.io/r/stats/lm.html),
      [`stats::glm()`](https://rdrr.io/r/stats/glm.html).

  `"cox"`, `"coxme"`

  :   Cox proportional hazards models via
      [`survival::coxph()`](https://rdrr.io/pkg/survival/man/coxph.html)
      and, with a random intercept or other `(1 | g)` terms,
      [`coxme::coxme()`](https://rdrr.io/pkg/coxme/man/coxme.html)
      (`random` required). Reports hazard ratios.

  `"lmer"`, `"glmer"`

  :   [`lme4::lmer()`](https://rdrr.io/pkg/lme4/man/lmer.html),
      [`lme4::glmer()`](https://rdrr.io/pkg/lme4/man/glmer.html);
      `random` required.

  `"clm"`, `"clmm"`

  :   ordered-factor outcome via
      [`ordinal::clm()`](https://rdrr.io/pkg/ordinal/man/clm.html) and
      [`ordinal::clmm()`](https://rdrr.io/pkg/ordinal/man/clmm.html)
      (`random` required). Reports cumulative odds ratios.

  `"mblogit"`

  :   unordered-factor outcome via
      [`mclogit::mblogit()`](https://melff.github.io/mclogit/reference/mblogit.html),
      optionally clustered through `random`. Reports one relative-risk
      ratio per non-reference outcome category, so
      [summary()](https://marcusvvlopes.github.io/timecoda/reference/summary.comp_iv.md)
      and
      [`predict_realloc()`](https://marcusvvlopes.github.io/timecoda/reference/predict_realloc.md)
      gain a `response_level` column. For a two-level outcome, a
      logistic regression with `model = "glm"` and `family = "binomial"`
      reports the effect as an odds ratio. With `random`, the
      `convergence` element of the fit holds mclogit's own convergence
      flag, which can read as converged although the inner iterations
      warned. A clustered fit has no log-likelihood.

  `"lmrob"`, `"rlmer"`

  :   robust linear models via
      [`robustbase::lmrob()`](https://rdrr.io/pkg/robustbase/man/lmrob.html)
      and
      [`robustlmm::rlmer()`](https://rdrr.io/pkg/robustlmm/man/rlmer.html)
      (`random` required).

  `"cox"`, `"coxme"`, `"clm"`, `"clmm"` and `"mblogit"` support
  `predict_realloc(type = "relative")` only.

- covariates:

  Optional character vector of covariate terms. An interaction between
  the composition and a covariate (a term involving `ilr` coordinates)
  makes reallocation contrasts covariate-conditional;
  [`predict_realloc()`](https://marcusvvlopes.github.io/timecoda/reference/predict_realloc.md)
  then needs `at` or `by`.

- family:

  A GLM family (object or name, e.g., `"binomial"`), required for
  `model = "glm"` and `model = "glmer"` and rejected by the other
  engines. A family given by name is resolved in the calling
  environment, so a locally defined family function works. Its link must
  be `"identity"` (the effect is a mean difference), `"log"` (a rate
  ratio for Poisson, a risk ratio for binomial, a ratio of means
  otherwise) or `"logit"` (an odds ratio). Any other link, such as
  `binomial(link = "probit")` or the default inverse link of
  [`stats::Gamma()`](https://rdrr.io/r/stats/family.html), is an error.

- random:

  Random-effects specification. A character vector of grouping column
  names (each expanded to `(1 | name)`) and/or full `lme4` terms (any
  element containing `|` is used verbatim), such as
  `random = c("province", "school")`. Required by `"coxme"`, `"lmer"`,
  `"glmer"`, `"clmm"` and `"rlmer"`, optional for `"mblogit"`, and
  rejected by the unclustered engines. `"mblogit"` takes one grouping
  term, passed as its own `random` argument: a grouping name gives a
  random intercept (`~ 1 | school`), and a bar term keeps its left-hand
  side (`"(ilr1 | school)"` gives the random slope `~ ilr1 | school`).

- time, event:

  Survival outcome columns for `model = "cox"` and `"coxme"`: `time`
  positive follow-up times, `event` 1 = event, 0 = censored.

- formula:

  Optional model formula, an alternative to
  `outcome`/`covariates`/`random`/`time`/`event`/`moderator` for models
  with further interactions and non-linear terms, fitted as written. It
  is mutually exclusive with the discrete arguments, `moderator` and
  `curvature`. The left-hand side is the outcome, or a
  `survival::Surv(time, event)` call for `model = "cox"` and `"coxme"`.
  On the right-hand side, the marker
  [`ilr()`](https://marcusvvlopes.github.io/timecoda/reference/ilr.md)
  stands for the composition (expanded to `ilr1 + ... + ilr{D-1}`), and
  everything else passes through to the model function, so
  [`poly()`](https://rdrr.io/r/stats/poly.html),
  [`splines::ns()`](https://rdrr.io/r/splines/ns.html),
  [`I()`](https://rdrr.io/r/base/AsIs.html), `a*b` and `(1 | g)` random
  terms all work.

  Non-linear terms may enter the composition too
  (`ilr()^2 + I(ilr()^2)`, `poly(ilr(), 2)`), with or without its main
  effect;
  [`ilr()`](https://marcusvvlopes.github.io/timecoda/reference/ilr.md)
  gives the expansion rules and the specifications that depend on the
  order of the parts. If
  [`ilr()`](https://marcusvvlopes.github.io/timecoda/reference/ilr.md)
  is absent, the composition is added as a main effect, even beside
  coordinate names (`ilr1`, ...) typed by hand. When
  [`ilr()`](https://marcusvvlopes.github.io/timecoda/reference/ilr.md)
  is present but a coordinate's main effect is not, `comp_iv()` warns
  that estimates may be biased. A formula fit, even
  `outcome ~ ilr() * x`, carries no `moderation` test and no default
  `by` for
  [`predict_realloc()`](https://marcusvvlopes.github.io/timecoda/reference/predict_realloc.md);
  use `moderator` for those.

  A variable the formula names must be a column of the `compD` data or a
  scalar in the formula's environment, such as a knots value. A longer
  object, such as a vector of the caller's, must be a column instead. A
  scalar is replaced by its value when the model is fitted, and the
  stored formula holds the value. The expressions in `covariates` follow
  the same rule, with or without `moderator` and `curvature`; their
  scalars are found in the environment `comp_iv()` is called from.

- moderator:

  Optional name of a single covariate to test as a moderator of the
  composition's effect. It adds the `ilr() * moderator` interaction, so
  the fit equals `formula = outcome ~ ilr() * moderator` (plus any
  `covariates` and `random` terms). Unlike that formula fit, the
  moderated fit stores the joint Wald test of the moderation in
  `moderation`. A factor or character moderator also becomes the default
  `by` of
  [`predict_realloc()`](https://marcusvvlopes.github.io/timecoda/reference/predict_realloc.md),
  which then gives one estimate per level. The test identifies the
  moderator's coefficients by the factor's own contrasts or those set
  with `options(contrasts =)`; a coding set through the backend's
  `contrasts` argument in `model_args` is not recognized. `outcome`,
  `time` and `event` follow the same rules as without a moderator. The
  moderator must not also appear in `covariates`, and `moderator` is
  mutually exclusive with `formula`.

- curvature:

  Logical (default `FALSE`). When `TRUE`, the composition enters as the
  full quadratic in the ILR coordinates, `ilr()^2 + I(ilr()^2)`: the D-1
  main effects, every pairwise product and every square, (D-1)(D+2)/2
  coefficients. The fit equals the one from
  `formula = outcome ~ ilr()^2 + I(ilr()^2)` with the same covariates
  and random terms, and every method reads it as that formula fit. The
  full quadratic is the same model in every ILR basis; a squares-only
  term or a spline, which changes with the order of the parts (see
  [`ilr()`](https://marcusvvlopes.github.io/timecoda/reference/ilr.md)),
  is written through `formula` instead. With `moderator`, the
  interaction stays on the main effects
  (`ilr()^2 + I(ilr()^2) + ilr() * moderator`), so `moderation` tests
  the same coefficients as without `curvature`. The joint Wald test of
  the squares and products is stored as `curvature_test`. The argument
  is an error beside `formula`, which writes its own curvature terms,
  and when `covariates` names a coordinate (`ilr1`, ...) or the
  [`ilr()`](https://marcusvvlopes.github.io/timecoda/reference/ilr.md)
  marker, both already held by the quadratic.

- include_total:

  Logical (default `FALSE`). When `TRUE`, fits the total-as-covariate
  model of Coenders et al. (2017), for compositions whose row totals
  vary (waking hours, dietary intake). The normalized log geometric mean
  of the zero-free parts on their original, pre-closure scale (the
  `parts_imputed` element of the `compD`), \\t =
  \sqrt{D}\\\log\big\[(x_1 \cdots x_D)^{1/D}\big\]\\, enters the model
  as the main effect `t_total`. The ILR coefficients then estimate the
  effects of proportional shifts at a constant total, and `t_total` the
  effect of scaling every part up proportionally. A (near-)constant
  total, such as the fixed 1440-minute days of `timeuse`, is an error.
  Reallocations hold `t_total` at `total_ref`, and the Wald test of
  `t_total` is `total_test`, shown by
  [summary()](https://marcusvvlopes.github.io/timecoda/reference/summary.comp_iv.md).
  The option combines with every engine, `covariates`, `formula`,
  `moderator` and `curvature`; `t_total` itself cannot be named in
  `covariates`, `formula`, `moderator` or `model_args`.

- model_args:

  A `list(...)` of further arguments passed straight to the fitting
  function of the chosen backend (e.g.,
  `model_args = list(weights = w, na.action = stats::na.omit)` for
  [`stats::glm()`](https://rdrr.io/r/stats/glm.html), or
  `list(nAGQ = 0)` for
  [`lme4::glmer()`](https://rdrr.io/pkg/lme4/man/glmer.html)), and the
  route to any backend argument `comp_iv()` does not expose itself. The
  arguments are captured unevaluated, so a bare column name works
  (`weights = w` refers to a column of the composition data, exactly as
  in a direct [`stats::glm()`](https://rdrr.io/r/stats/glm.html) call)
  and the fitted object's stored call stays readable. A name that is not
  a column is found in the environment `comp_iv()` is called from (with
  or without `moderator`), or, with a `formula`, in that formula's
  environment.

  An engine's own control object goes here too, such as
  `model_args = list(control = lme4::lmerControl(optimizer = "bobyqa"))`
  to retry a mixed fit with a different optimizer, or
  `list(control = stats::glm.control(maxit = 100))`. Each backend
  documents its own arguments, control functions included. An argument
  the engine already sets (`trace` for `"mblogit"`, `model` for `"cox"`,
  `Hess` for `"clmm"`) is an error. A `link` other than `"logit"` is an
  error for `"clm"` and `"clmm"`, whose effects are cumulative odds
  ratios.

## Value

An object of class `c("comp_iv", "timecoda_fit")`: the `timecoda_fit`
metadata (fitted `model`, `model_type`, `parts`, `total`, `ref_comp`,
`family`, ...) augmented with the elements below. `family` is the
resolved GLM family object (class `"family"`, carrying its `link`) for
`model = "glm"` and `"glmer"`, and `NULL` for every other engine. The
backend object `model`, like each of `models`, carries its formula and
fitting data in its call. R's own tools therefore re-fit it:
[`stats::update()`](https://rdrr.io/r/stats/update.html),
[`stats::drop1()`](https://rdrr.io/r/stats/add1.html) and
[`stats::anova()`](https://rdrr.io/r/stats/anova.html) between two such
fits work on it directly for every engine.
[`survival::cox.zph()`](https://rdrr.io/pkg/survival/man/cox.zph.html)
works on a plain Cox fit (`model = "cox"`). Current versions of survival
do not support it on a `"coxme"` fit, and it should not be called on
one. Printing that object itself (`print(fit$model)`,
`summary(fit$model)`) shows the data inside the call;
[`print()`](https://rdrr.io/r/base/print.html) and
[`summary()`](https://rdrr.io/r/base/summary.html) of the `comp_iv()`
fit do not.

- `omnibus`:

  A one-row tibble (`term`, `statistic`, `df`, `p.value`): the joint
  Wald test that every coefficient of a term in the ILR coordinates
  alone is zero (i.e., that the composition as a whole has no effect).
  The statistic is referred to a chi-squared distribution on as many
  degrees of freedom as there are such coefficients. These are D-1 for
  the main effects, plus one per square, product or basis column of a
  non-linear term ((D-1)(D+2)/2 for the full quadratic
  `ilr()^2 + I(ilr()^2)`), and one set per non-reference category for a
  multinomial model. A formula with no term in the coordinates alone
  (`ilr():sex`) tests every coefficient crossing a coordinate with a
  covariate instead. Aliased (`NA`) coefficients, such as a linear term
  repeated by `ilr() + poly(ilr(), 2)`, are left out of the test and of
  its degrees of freedom. The test is `NA` on 0 degrees of freedom when
  no fixed-effect coefficient carries the composition, for example when
  it enters only through a random slope
  (`zbmi ~ age + (ilr() | school)`). The test is also `NA`, with a
  warning, when every composition coefficient is aliased (a composition
  that is the same on every row, or collinear with another term). For a
  linear model the statistic equals its degrees of freedom times the F
  statistic of the same hypothesis. The chi-squared reference ignores
  the uncertainty in the residual variance, so in small samples the test
  is slightly anti-conservative, with a p-value a little smaller than
  that of the exact F test.

- `comp`:

  The `compD` object.

- `models`:

  A named list (one element per part) of the rotated fits:
  `models[[part]]` is the model refitted with that part as the first
  pivot coordinate, so its `ilr1` term is the "part vs remaining" effect
  summarized by
  [summary()](https://marcusvvlopes.github.io/timecoda/reference/summary.comp_iv.md).
  The refits leave the canonical `model` unchanged, and no prediction
  path reads them. For a mixed model, one refit per part can be slow.

- `moderator`:

  The moderator column name (or `NULL`).

- `moderation`:

  When `moderator` is set, a one-row tibble with the joint Wald test
  that every composition-by-moderator coefficient is zero (`H0`: no
  moderation). Its degrees of freedom are (L-1)(D-1) for an L-level
  factor moderator and D-1 for a numeric one, times K-1 for a
  multinomial fit with K outcome categories. `NULL` otherwise.

- `curvature_test`:

  A one-row tibble when the terms in the ILR coordinates alone are the
  full quadratic: a square written `I(ilrK^2)` for every coordinate, a
  product `ilrJ:ilrK` for every pair, and no other term above degree
  one. It is `NULL` otherwise, including for
  [`poly()`](https://rdrr.io/r/stats/poly.html) and spline terms and for
  a square or product written in another form (such as
  `I(ilr1 * ilr2)`). The tibble holds the joint Wald test that every
  square and product coefficient is zero (`H0`: the effect is linear in
  the coordinates). The degrees of freedom are (D-1)D/2, times K-1 for a
  multinomial fit with K outcome categories, less any aliased
  coefficient. The test is the same whether the quadratic comes from
  `curvature = TRUE` or from a formula, with or without a moderator, and
  in every ILR basis.

- `random`:

  The resolved `random` specification (character vector of grouping
  terms), or `NULL` for a non-clustered fit.

- `fit_warnings`:

  A tibble of the backend warnings/messages captured across the
  canonical and rotated fits (columns `source`, `type`, `message`; zero
  rows when every fit was clean). See the *Convergence* section.

- `convergence`:

  A one-row tibble (`converged`, `singular`, `n_warnings`) probing the
  canonical fit's convergence state.

- `total_test`:

  When `include_total = TRUE`, a one-row tibble: the Wald test of the
  `t_total` coefficient (`H0`: the total carries no effect beyond the
  composition), on 1 degree of freedom, or K-1 for a multinomial fit
  with K outcome categories. `NULL` otherwise.

- `total_ref`:

  When `include_total = TRUE`, the mean of `t_total`, which is `t_total`
  at the geometric-mean composition of the original-scale parts; every
  reallocation holds `t_total` at this value. `NULL` otherwise.

## Total as a covariate

The \\\sqrt{D}\\ factor makes `t_total` the coordinate of \\\log x\\
along \\\mathbf{1}/\sqrt{D}\\, which completes the ILR basis to an
orthonormal one (Coenders et al. 2017, eq. 3.2). Any positive rescaling,
such as a \\\sqrt{D-1}\\ convention, rescales its coefficient inversely
and leaves the fit, the predictions and `total_test` unchanged.

## Convergence

Every engine is fitted with its backend warnings and messages captured
(and still re-emitted). The fit carries two additive elements:
`fit_warnings`, a tibble of the captured conditions (`source`, `type`,
`message`; zero rows when clean) whose `source` is `"canonical"` or a
rotated part name, and `convergence`, a one-row tibble (`converged`,
`singular`, `n_warnings`).
[summary()](https://marcusvvlopes.github.io/timecoda/reference/summary.comp_iv.md)
prints the captured messages and, for a clustered or iterative fit, a
`Convergence: OK` line when clean;
[print()](https://marcusvvlopes.github.io/timecoda/reference/print.timecoda_fit.md)
flags the warning count. `"cox"` and `"coxme"` fits are not probed
(`converged` is `NA`), so they print no such line. An `"lmer"` or
`"glmer"` fit that fails to converge or is singular can be refitted with
a different optimizer through `model_args` (e.g.,
`model_args = list(control = lme4::lmerControl(optimizer = "bobyqa"))`).

## Degrees of freedom for clustered fits

For an `"lmer"` fit, the per-part standard errors, intervals and
p-values that
[summary()](https://marcusvvlopes.github.io/timecoda/reference/summary.comp_iv.md)
reports and the reallocations of
[`predict_realloc()`](https://marcusvvlopes.github.io/timecoda/reference/predict_realloc.md)
use the same degrees-of-freedom method, chosen by the rule `emmeans`
applies. The default is the Kenward-Roger approximation through
pbkrtest. The Satterthwaite approximation through lmerTest is used when
`emmeans::emm_options(lmer.df = "satterthwaite")` is set, or when
pbkrtest is absent, disabled or the fit has more observations than
`emmeans::emm_options(pbkrtest.limit =)` (3000). Both pbkrtest and
lmerTest are optional, and when neither can be used, or
`lmer.df = "asymptotic"` is set, the degrees of freedom are asymptotic
(normal). `"glmer"`, `"cox"`, `"coxme"` and `"rlmer"` fits are
asymptotic in both places (robustlmm reports no degrees of freedom), and
`"lmrob"` fits use their residual degrees of freedom.

## References

Coenders, G., Martín-Fernández, J. A., & Ferrer-Rosell, B. (2017). When
relative and absolute information matter: Compositional predictor with a
total in generalized linear models. *Statistical Modelling*, 17(6),
494-512.
[doi:10.1177/1471082X17710398](https://doi.org/10.1177/1471082X17710398)

Cox, D. R. (1972). Regression models and life-tables. *Journal of the
Royal Statistical Society: Series B (Methodological)*, 34(2), 187-202.
[doi:10.1111/j.2517-6161.1972.tb00899.x](https://doi.org/10.1111/j.2517-6161.1972.tb00899.x)

Kenward, M. G., & Roger, J. H. (1997). Small sample inference for fixed
effects from restricted maximum likelihood. *Biometrics*, 53(3),
983-997. [doi:10.2307/2533558](https://doi.org/10.2307/2533558)

McCullagh, P. (1980). Regression models for ordinal data. *Journal of
the Royal Statistical Society: Series B (Methodological)*, 42(2),
109-127.
[doi:10.1111/j.2517-6161.1980.tb01109.x](https://doi.org/10.1111/j.2517-6161.1980.tb01109.x)

## See also

[`declare_composition()`](https://marcusvvlopes.github.io/timecoda/reference/declare_composition.md)
to build the `compD` input;
[`summary.comp_iv()`](https://marcusvvlopes.github.io/timecoda/reference/summary.comp_iv.md)
and
[`plot.comp_iv()`](https://marcusvvlopes.github.io/timecoda/reference/plot.comp_iv.md)
to summarize and visualize the fit;
[`predict_realloc()`](https://marcusvvlopes.github.io/timecoda/reference/predict_realloc.md)
for reallocation predictions;
[`fit_check()`](https://marcusvvlopes.github.io/timecoda/reference/fit_check.md)
for residual diagnostics;
[`ilr()`](https://marcusvvlopes.github.io/timecoda/reference/ilr.md) for
the formula-interface composition marker.

Other composition as predictor:
[`plot.comp_iv()`](https://marcusvvlopes.github.io/timecoda/reference/plot.comp_iv.md),
[`summary.comp_iv()`](https://marcusvvlopes.github.io/timecoda/reference/summary.comp_iv.md)

## Examples

``` r
cd <- declare_composition(timeuse, c("MVPA", "LPA", "SB", "SPT"),
                          closure = 1440)
m <- comp_iv(cd, outcome = "obesity", model = "glm", family = "binomial")
m$omnibus
#> # A tibble: 1 × 4
#>   term        statistic    df    p.value
#>   <chr>           <dbl> <int>      <dbl>
#> 1 composition      28.9     3 0.00000233

# The full quadratic in the coordinates, a predicted outcome that can curve
curved <- comp_iv(cd, outcome = "zbmi", model = "lm", curvature = TRUE)
curved$omnibus
#> # A tibble: 1 × 4
#>   term        statistic    df  p.value
#>   <chr>           <dbl> <int>    <dbl>
#> 1 composition      310.     9 2.10e-61
```
