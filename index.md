# timecoda

`timecoda` implements compositional data analysis (CoDA) workflows for
time-use epidemiology in R. It provides a consistent interface for
declaring and screening compositions, fitting regression models on
isometric log-ratio (ILR) coordinates, estimating isotemporal time
reallocations, and searching for optimal compositions (“Goldilocks
zones”).

The statistical framework draws on the foundational time-use
epidemiology literature (Pedišić, 2014; Pedišić et al., 2017; Dumuid,
Stanford, Martín-Fernández, et al., 2018; Dumuid, Stanford, Pedišić, et
al., 2018; Dumuid et al., 2019, 2020). For more on the field, visit the
[International Network of Time-Use Epidemiologists
(INTUE)](https://www.intue.org/).

## Installation

`timecoda` requires R 4.4.0 or newer. Install the development version
from GitHub with `pak`:

``` r

install.packages("pak")
pak::pak("MarcusVVLopes/timecoda")
```

To install suggested packages for optional features (such as `gt`
tables, `DHARMa` diagnostics, and multilevel or robust model engines):

``` r

pak::pak("MarcusVVLopes/timecoda", dependencies = TRUE)
```

## Quick example

Using the bundled `timeuse` dataset, a typical analysis declares the
daily composition, fits a compositional regression model, and translates
the ILR estimates into minute-based time reallocations:

``` r

library(timecoda)

cd <- declare_composition(
  data    = timeuse,
  parts   = c("MVPA", "LPA", "SB", "SPT"),
  closure = 1440
)

fit_ob <- comp_iv(
  comp    = cd,
  outcome = "obesity",
  model   = "glm",
  family  = "binomial"
)

summary(object = fit_ob)
predict_realloc(object = fit_ob, from = "SB", to = "MVPA", delta = 15)
```

## Getting started

For an end-to-end walkthrough on the `timeuse` dataset, start with
[Getting started with
timecoda](https://marcusvvlopes.github.io/timecoda/articles/timecoda.html).
Function documentation is available on the [package
website](https://marcusvvlopes.github.io/timecoda/), alongside
topic-specific articles:

- [Working with a compositional
  variable](https://marcusvvlopes.github.io/timecoda/articles/composition-variable.html)
  — validating compositions, imputing zeroes, flagging outliers, and
  descriptive plots.
- [Modelling an outcome from a
  composition](https://marcusvvlopes.github.io/timecoda/articles/composition-predictor.html)
  — fitting ILR regression models with
  [`comp_iv()`](https://marcusvvlopes.github.io/timecoda/reference/comp_iv.md),
  checking diagnostics, and plotting response surfaces.
- [Predicting time
  reallocations](https://marcusvvlopes.github.io/timecoda/articles/predicting-reallocations.html)
  — pairwise, one-to-all, and custom reallocations, subgroup moderation,
  and transfer curves.
- [Finding the Goldilocks
  zone](https://marcusvvlopes.github.io/timecoda/articles/optimal-composition.html)
  — single- and multi-outcome optimization over feasible compositions.
- [Advanced model engines and the formula
  interface](https://marcusvvlopes.github.io/timecoda/articles/advanced-models.html)
  — multilevel, survival, ordinal, multinomial, and robust engines, plus
  custom
  [`ilr()`](https://marcusvvlopes.github.io/timecoda/reference/ilr.md)
  formulas and varying-total compositions.

## Getting help and contributing

To ask a question or report a bug, open an issue on the [issue
tracker](https://github.com/MarcusVVLopes/timecoda/issues). Before
writing code for a new feature, please discuss it in an issue first. See
the [contributing
guide](https://github.com/MarcusVVLopes/timecoda/blob/main/.github/CONTRIBUTING.md)
for workflow instructions, and note that participation is governed by
the [code of
conduct](https://github.com/MarcusVVLopes/timecoda/blob/main/.github/CODE_OF_CONDUCT.md).

## Development

We build `timecoda` around verifiable statistical accuracy and
transparent development practices. We specify all statistical methods
ourselves and benchmark each estimator before release against
established R packages (such as `robCompositions`, `compositions`, and
`codaredistlm`), hand-fitted models, or published results. These
benchmarks live in `tests/testthat/test-benchmark-*.R`, so anyone can
inspect and rerun them.

We also use AI tools to assist with documentation, unit tests, code
review, and refactoring. We review, test, and take full responsibility
for everything shipped in the package.

## Citation

To cite `timecoda`, run `citation("timecoda")`:

> Lopes MVV, Tozetto WR, de Camargo R, Zhang K (2026). timecoda:
> Compositional Data Analysis for Time-Use Behaviour Research. R package
> version 0.6.1, <https://github.com/MarcusVVLopes/timecoda>.

## Licence

`timecoda` is distributed under the GNU General Public License, version
3 or any later version. For the full text, see
[LICENSE.md](https://github.com/MarcusVVLopes/timecoda/blob/main/LICENSE.md).

## References

Dumuid, D., Stanford, T. E., Martín-Fernández, J. A., Pedišić, Ž.,
Maher, C. A., Lewis, L. K., Hron, K., Katzmarzyk, P. T., Chaput, J. P.,
Fogelholm, M., Hu, G., Lambert, E. V., Maia, J., Sarmiento, O. L.,
Standage, M., Barreira, T. V., Broyles, S. T., Tudor-Locke, C.,
Tremblay, M. S., & Olds, T. (2018). Compositional data analysis for
physical activity, sedentary time and sleep research. *Statistical
Methods in Medical Research*, 27, 3726-3738.
<https://doi.org/10.1177/0962280217710835>

Dumuid, D., Stanford, T. E., Pedišić, Ž., Maher, C., Lewis, L. K.,
Martín-Fernández, J. A., Katzmarzyk, P. T., Chaput, J. P., Fogelholm,
M., Standage, M., Tremblay, M. S., & Olds, T. (2018). Adiposity and the
isotemporal substitution of physical activity, sedentary time and sleep
among school-aged children: A compositional data analysis approach. *BMC
Public Health*, 18, 311. <https://doi.org/10.1186/s12889-018-5207-1>

Dumuid, D., Pedišić, Ž., Stanford, T. E., Martín-Fernández, J. A., Hron,
K., Maher, C. A., Lewis, L. K., & Olds, T. (2019). The compositional
isotemporal substitution model: A method for estimating changes in a
health outcome for reallocation of time between sleep, physical activity
and sedentary behaviour. *Statistical Methods in Medical Research*, 28,
846-857. <https://doi.org/10.1177/0962280217737805>

Dumuid, D., Pedišić, Ž., Palarea-Albaladejo, J., Martín-Fernández, J.
A., Hron, K., & Olds, T. (2020). Compositional data analysis in time-use
epidemiology: What, why, how. *International Journal of Environmental
Research and Public Health*, 17, 2220.
<https://doi.org/10.3390/ijerph17072220>

Pedišić, Ž. (2014). Measurement issues and poor adjustments for physical
activity and sleep undermine sedentary behaviour research: The focus
should shift to the balance between sleep, sedentary behaviour, standing
and activity. *Kinesiology*, 46, 135-146.

Pedišić, Ž., Dumuid, D., & Olds, T. S. (2017). Integrating sleep,
sedentary behaviour, and physical activity research in the emerging
field of time-use epidemiology: Definitions, concepts, statistical
methods, theoretical framework, and future directions. *Kinesiology*,
49, 252-269.
