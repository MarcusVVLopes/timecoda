# timecoda

`timecoda` is an R package for implementing current methods of time-use
epidemiology, in particular compositional data analysis applied to
behavioural data.

The framework behind the package comes from a series of publications.
Together they set out the definitions, concepts, statistical methods and
research areas of time-use epidemiology (Pedišić, 2014; Pedišić et al.,
2017; Dumuid, Stanford, Martín-Fernández, et al., 2018; Dumuid,
Stanford, Pedišić, et al., 2018; Dumuid et al., 2019, 2020). For further
information on the topic, see the [International Network of Time-Use
Epidemiologists (INTUE)](https://www.intue.org/).

## Installation

`timecoda` needs R 4.4.0 or newer. To install it from GitHub, use the
`pak` package:

``` r

install.packages("pak")
pak::pak("MarcusVVLopes/timecoda")
```

To install all the suggested packages for optional features, run the
same call with `dependencies = TRUE`.

``` r

pak::pak("MarcusVVLopes/timecoda", dependencies = TRUE)
```

## Getting started

For a brief tutorial using the bundled `timeuse` data, see [Getting
started with
timecoda](https://marcusvvlopes.github.io/timecoda/articles/timecoda.html).

To find the reference pages, check the [package
website](https://marcusvvlopes.github.io/timecoda/). For the other
articles, see the list below:

- [Working with a compositional
  variable](https://marcusvvlopes.github.io/timecoda/articles/composition-variable.html)
- [Modelling an outcome from a
  composition](https://marcusvvlopes.github.io/timecoda/articles/composition-predictor.html)
- [Predicting time
  reallocations](https://marcusvvlopes.github.io/timecoda/articles/predicting-reallocations.html)
- [Advanced model engines and the formula
  interface](https://marcusvvlopes.github.io/timecoda/articles/advanced-models.html)
- [Finding the Goldilocks
  zone](https://marcusvvlopes.github.io/timecoda/articles/optimal-composition.html)

## Getting help and contributing

To ask a question or report a bug, use the [issue
tracker](https://github.com/MarcusVVLopes/timecoda/issues). Before you
write any code for a new feature, discuss it in an issue. For
instructions on proposing a change, see the [contributing
guide](https://github.com/MarcusVVLopes/timecoda/blob/main/.github/CONTRIBUTING.md),
and follow the [code of
conduct](https://github.com/MarcusVVLopes/timecoda/blob/main/.github/CODE_OF_CONDUCT.md)
when you contribute.

## Development

The developers specify the statistical methods and write the first
implementation of each feature. Before release, each statistical method
is benchmarked against established R packages (such as
`robCompositions`, `compositions` and `codaredistlm`), hand-fitted
models or published results. These benchmarks are part of the package’s
tests (`tests/testthat/test-benchmark-*.R`), so anyone can rerun them.

AI tools (Claude, by Anthropic) help with documentation, tests, code
review and refactoring. The developers review everything these tools
produce and remain responsible for the package.

## Citation

To cite `timecoda`, run `citation("timecoda")`:

> Lopes MVV, Tozetto WR, de Camargo R, Zhang K (2026). timecoda:
> Compositional Data Analysis for Time-Use Behaviour Research. R package
> version 0.6.0, <https://github.com/MarcusVVLopes/timecoda>.

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
