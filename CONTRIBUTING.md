# Contributing to timecoda

There are three ways to contribute.

- To report a bug, open an issue with the “Bug report” form and include
  code that reproduces the problem on `timeuse` or on simulated data.
- To propose a feature, open an issue with the “Feature proposal” form.
  If the feature implements a new method, provide the references that
  support it.
- To fix the documentation, open a pull request. A typo or an unclear
  sentence needs no issue.

A new method, a new function, a new model engine or a new package
dependency should be agreed in its issue before any code is written.

## Opening a pull request

1.  Fork the repository and create a branch from `main`.
2.  Make one change per pull request.
3.  Add a test and a line in `NEWS.md` if you can. Both are added during
    the review when they are missing.
4.  Open the pull request against `main`.

## What happens next

- GitHub runs `R CMD check` on five platforms for the pull request and
  again for each update. The first run for a first-time contributor
  waits for the maintainer’s approval.
- The maintainer reviews the change against the design rules of the
  package: how a composition enters a model, how predictions and their
  confidence intervals are computed, and how results are returned.
- The maintainer may push edits to the branch to meet those rules, each
  as a separate commit. Pull those commits before working on the branch
  again.
- `main` takes no direct push. A pull request is merged after the
  maintainer approves it.

## Licence and credit

Contributions are licensed under the package’s licence, GPL (\>= 3).

Every merged pull request is credited in `NEWS.md` with the
contributor’s GitHub handle, and a contributor who adds a feature is
listed in `DESCRIPTION`.

## Code of conduct

Everyone taking part in this project follows the [code of
conduct](https://github.com/MarcusVVLopes/timecoda/blob/main/.github/CODE_OF_CONDUCT.md).
