# darthtools 0.4.1

## Bug fixes

* `gen_wcc(method = "Simpson1/3")` assigned the 4/3 and 2/3 weights by
  position in the 1-based R vector instead of by time t = 0, 1, ..., n. It
  returned `1/3, 2/3, 4/3, 2/3, ..., 1/3` instead of `1/3, 4/3, 2/3, 4/3, ...,
  1/3`, and the weights summed to n - 2/3 (even n) or n - 1/3 (odd n) instead
  of n. Weights are now indexed by time. Reported independently by Lauren N.
  Steimle (Georgia Institute of Technology) and @ambtnhlr
  (DARTH-git/cohort-modeling-tutorial-intro#2). Thank you both.

## Changes in behaviour

* For an odd number of cycles, `gen_wcc(method = "Simpson1/3")` now applies
  the composite Simpson's 1/3 rule over the first n - 3 cycles and Simpson's
  3/8 rule over the last three, so polynomials up to degree three are still
  integrated exactly.
* `gen_wcc()` now errors when `n_cycles` is not a single positive integer, and
  when `n_cycles = 1` is used with Simpson's rule.

## Infrastructure

* Added unit tests with testthat (3rd edition), starting with `gen_wcc()`.
* Added a GitHub Actions workflow that runs `R CMD check` on every push and
  pull request to `main`.
* Fixed pre-existing `R CMD check` warnings: `surv_prob()` no longer calls
  `require("MASS")` (MASS is now in Imports), and the `B` argument of
  `surv_prob()` and the `v_n_states` argument of `plot_trace_microsim()` are
  now documented.

# darthtools 0.4.0

### Breaking change (with deprecation)
- `t` is deprecated and replaced by `time`
- Monthly probabilities now use `time = 1/12`
- Legacy implementation archived on GitHub
- Version 0.3.0 saved as release v.0


## Breaking changes
- `rate_to_prob()` now interprets `time` as the length of the time step.
- Monthly probabilities should use `time = 1/12`.

## Deprecations
- The `t` argument is deprecated and will be removed in a future release.

## Notes
- The legacy implementation (pre–Jan 2026) is available at GitHub at release with tag v.0.3.0.
This can be installed using `remotes::install_github("DARTH-git/darthtools@v0.3.0")`



# darthtools 0.3.0

# darthtools 0.2.2

# darthtools 0.2.1

# darthtools 0.2.0

# darthtools 0.1.0

* Initial CRAN submission.
