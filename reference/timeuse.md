# Simulated 24-hour time-use dataset

A simulated dataset of daily movement-behaviour compositions for 300
individuals, used in place of real data by the `timecoda` examples and
tests. The four behaviour parts sum exactly to 1440 minutes (24 h) for
every individual, and the outcomes are generated with a compositional
effect (more moderate-to-vigorous activity and sleep lower adiposity;
more sedentary time raises it).

## Usage

``` r
timeuse
```

## Format

A data frame with 300 rows and 13 variables:

- id:

  Participant identifier (`"p001"` to `"p300"`).

- sex:

  Grouping factor, `"female"`/`"male"` (for ternary colouring and as a
  covariate).

- school:

  Clustering factor with 10 levels (`"school01"` to `"school10"`), for
  multilevel/mixed (`lme4`) models.

- region:

  Coarser clustering factor with 3 levels (`"A"`/`"B"`/`"C"`) that nests
  `school` deterministically (schools 01 to 03 in `"A"`, 04 to 06 in
  `"B"`, 07 to 10 in `"C"`), for a model with several random intercepts
  (`random = c("region", "school")`).

- MVPA:

  Moderate-to-vigorous physical activity (minutes/day).

- LPA:

  Light physical activity (minutes/day).

- SB:

  Sedentary behaviour (minutes/day).

- SPT:

  Sleep period time (minutes/day).

- zbmi:

  Continuous outcome: body-mass-index z-score.

- obesity:

  Binary outcome factor, `"no"`/`"yes"`.

- time:

  Survival outcome: follow-up time, in months, to onset of obesity or
  censoring.

- event:

  Survival outcome: obesity-onset indicator (`1` = onset, `0` =
  censored).

- grade:

  Ordinal outcome: an ordered academic-achievement band,
  `"below" < "expected" < "above"` (tertiles of a composition-driven
  latent score). Used by the ordinal model path and, coerced to an
  unordered factor, by the multinomial path.

## Source

Simulated; see `data-raw/timeuse.R` for the generating script and random
seed.

## Details

The data carry no substantive interpretation and serve only to exercise
the package's linear, logistic, Cox, ordinal, multinomial and multilevel
paths.
