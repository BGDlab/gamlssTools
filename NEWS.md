# gamlssTools (development version)

## Use Notes:
To test the development version (`dev` branch), you can install the development version.
However, be aware this will overwrite any stable `gamlssTools` version you may have 
downloaded:
```
remotes::install_github("BGDlab/gamlssTools@dev")
```

You can also clone the repo, but if you make edits, *please do so in a new branch*

## Data-free prediction

* Centile and z-score predictions on gamlss models no longer require the original 
  fitting data to be in scope. A model's distribution parameters are now reconstructed 
  directly from the stored fit: parametric terms from `coef()` and `model.matrix()`,
  `pb()` smooths from their stored linear coefficient plus stored interpolation 
  function, and `random()` effects from their stored per-level BLUPs. 
  This mirrors the approach taken by the gamlss2charts package.

* Data-free prediction is now the default for `gamlss` fits in `score_centiles()`,
  `centile_fan_values()`, `sigma_values()`, `remove_effects()`, `make_centile_fan()`,
  `plot_sigma()` and `trajectory_diff()`; supplying `fit_data` (the original fitting 
  data) forces the exact `gamlss::predictAll()` path instead.

* Models containing a smoother that cannot be reconstructed this way (`cs()`, `ps()`,
  `ga()`, `s()`) are detected up front and raise an informative error asking for
  `fit_data`, rather than failing obscurely.

* Models containing a parametric term whose columns are computed from the data as a
  whole (`poly()`, `ns()`, `bs()`, `scale()`, `cut()`) are likewise detected up front.
  Such a term is a different basis every time it is rebuilt, so applying the stored
  coefficients to a basis rebuilt from `newdata` alone gives silently wrong
  predictions. Precompute the basis as plain columns and refit to stay on the
  data-free path, as for `pb(log(Age))`.

* `gamlss2` fits continue to use gamlss2's own `predict()` and are unaffected.

## Bug fixes

* Data-free prediction returned silently wrong values for models with a smooth
  on a *transformed* covariate, e.g. `pb(log(Age))`. The covariate was resolved
  with `all.vars(...)[1]`, yielding `Age`, so the stored interpolation function
  -- a smooth of `log(Age)` -- was evaluated at raw `Age`. On a converged test
  fit this put fitted `mu` in the range 6.07-25.97 where `predictAll()` gives
  5.02-8.43: a maximum error of 17.5 on a response whose SD is 0.93, with no
  error or warning. `.datafree_eligible_gamlss()` now resolves the covariate
  with `match.call()` and requires it to be a bare column name, so such models
  fall back to the exact `predictAll()` path (or raise an informative error)
  instead. Precomputing the transform as a column, `pb(logAge)`, keeps a model
  on the data-free path and is exact to ~7e-15.

## Sharing models

* `sanitize_gamlss()` strips every per-observation component out of a fitted
  `gamlss` model while preserving data-free prediction, so models can be shared
  with collaborators who should not be able to reconstruct the training data.
  It keeps only what the data-free path reads (coefficients, formulas, links,
  factor levels and smoother summaries), rebuilds each `pb()` interpolation
  function on a regular grid (its closure otherwise stores the sorted covariate
  values), and clears the `random()` grouping column, fitted values and standard
  errors. Optional `xranges` lets you declare covariate ranges rather than
  disclose your data's min/max, and `random_level_map` pseudonymises `random()`
  level names. Models containing a smoother that cannot be reconstructed
  data-free (`cs()`, `ps()`, `ga()`, `s()`) are refused.

* `audit_gamlss()` recursively walks an object -- lists, attributes and closure
  environments -- and reports every atomic vector longer than `max_len`,
  flagging those of length exactly `n`. `pb()` interpolation functions are
  checked for an evenly spaced grid rather than by length, so a rebuilt smooth
  is not confused with surviving covariate data (and vice versa) when `grid_n`
  happens to match the sample size.

* `check_equivalent()` confirms a sanitized model reproduces the original's
  standardized scores on a covariate grid, so you can verify `grid_n` and the
  `pb()` ranges before sharing. It compares z-scores -- `qnorm()` of the centile,
  as in `score_centiles(standardize = TRUE)` -- rather than raw distribution
  parameters: a score is what collaborators consume and is dimensionless,
  whereas an absolute tolerance on `mu` is meaningless without knowing the
  response's scale (multiply the outcome by 10,000 and `mu`'s error multiplies
  by 10,000 too), and no single tolerance suits `mu`, `sigma` and `nu` at once.
  z-scores rather than the centiles they derive from because centile differences
  are compressed in the tails, so a centile tolerance under-weights exactly the
  tail errors that matter for screening; in z units the error runs roughly flat
  across the distribution (measured spread 1.8x, against ~21x in centile units),
  so one threshold means the same thing everywhere. `tol` defaults to 1e-6 in z
  units -- a millionth of a standard deviation. Note that `grid_n = 500` does not
  always reach it: a smooth spending ~6.5 edf measured 1.7e-06 at `grid_n = 500`
  and 5.3e-09 at 2000. Per-parameter differences are attached as the
  `"parameters"` attribute.

* `check_equivalent()` warns when `newdata` coincides with the rebuilt spline's
  own grid. A natural spline reproduces its nodes exactly, so such a comparison
  reports 0 regardless of how good the rebuild is. This is easy to hit by
  accident: `sim_grid()` returns 500 points over the covariate range and
  `grid_n` defaults to 500, computed by the same `seq()`, so passing `sim_grid()`
  output straight in aliases onto every node. Use any other number of points.

* `check_equivalent()` gains `fit_data`. By default the original model is
  predicted data-free, isolating what sanitizing changed; supplying `fit_data`
  takes the reference from `predictAll()` with the fitting data in scope
  instead, validating the whole chain a collaborator depends on (data-based ->
  data-free -> sanitized) in one number. Note that on some fits `predictAll()`
  refits to make "safe" predictions and its own answer differs from the
  reconstruction by more than sanitizing does, so the two references are worth
  running separately.

  Pass a *dense* grid either way: the error is pointwise, and on one ~11 edf fit
  50 observed rows reported 4.6e-07 where a 2000-point grid over the same range
  reported 5.8e-06.

* `sanitize_gamlss()` refuses models whose `pb()` smooths were fitted with
  `pb.control(quantiles = TRUE)`. Such a fit places its knots at sample
  quantiles of the covariate rather than on an even grid, and the knot vector is
  stored in the model: on a default `inter = 20` fit, 18 of the 24 stored knots
  are exact quantiles of the fitting data. That describes the covariate's
  distribution rather than just its range, so the model is rejected rather than
  sanitized. Refit with the default evenly spaced knots to share it.

* `sanitize_gamlss()` warns when `grid_n` looks too coarse for a smooth's
  effective degrees of freedom (fewer than `points_per_edf = 200` grid points
  per edf). How faithful a rebuilt smooth is depends on how wiggly it is, and
  `grid_n` now defaults to 2000 so that the result meets `check_equivalent()`'s
  1e-6 tolerance: a smooth spending 6.5 edf measured 1.7e-06 at `grid_n = 500`
  against 5.3e-09 at 2000, and one spending 12.8 edf needed the full 2000 to
  reach 9.1e-07. Reaching 1e-6 took roughly 100-160 points per edf across the
  fits measured, so the 200 threshold errs toward warning.

## Diagnostics

* `centile_coverage()` gains `centiles`, so pre-calculated centiles can be passed in
  instead of re-scoring the model. Accepts a numeric vector, the name of a column of
  `data`, or the dataframe returned by `score_centiles(standardize=TRUE)`. `gamlssModel`
  and `data` are now optional: supply either a model (to score `data`) or `centiles`,
  and `data` is only needed when grouping with `group`/`interval_var`.

* `centile_coverage()` also gains `batch_term`, passed through to `score_centiles()` so coverage can
  be checked on data containing unseen levels of a site/study/batch variable. It only applies when
  scoring from `gamlssModel` - combining it with `centiles` is an error.

* `centile_coverage()` now errors on an unrecognized argument rather than letting `...` swallow it.
  `...` is only ever passed to `cut_interval()`, so without `interval_var` a misspelled (or
  not-yet-installed) argument used to disappear silently.

## Out-of-Sample prediction draft

* `score_centiles()` gains `batch_term`, for scoring data containing levels of a
  site/study/batch variable the model was never fit on. Offsets for unseen levels are
  estimated and removed via `gamlss2charts::predict_score()`; rows with known levels
  go through the standard path. HOWEVER, this is validated on and requires the 
  `dev` branch of suggested package gamlss2charts
  (`remotes::install_github("andy1764/gamlss2charts@dev")`), which is what
  `Remotes:` now points to. Full implementation is pending the approval and
  update of gamlss2charts.
  
## Renamed functions

All of the below keep a deprecated alias under the old name, which warns via
`.Deprecated()` and forwards to the new function. The exceptions are noted under
Breaking changes above.

| Old | New |
|-----|-----|
| `centile_predict()`     | `centile_fan_values()`   |
| `pred_og_centile()`     | `score_centiles()`       |
| `sigma_predict()`       | `sigma_values()`         |
| `resid_data()`          | `remove_effects()`       |
| `sim_data()`            | `sim_grid()`             |
| `cent_cdf()`            | `centile_coverage()`     |
| `centile_fan_lifespan()`| `centile_fan_brainchart()` |
| `pred_centile()`        | `.centile_value()` (internal) |

The `gamlss` and `gamlss2` S3 methods moved with their generics.

## Renamed arguments

Argument names are now consistent about which data is which: **`data`** is always the
data you are asking about, and **`fit_data`** is always the data the model was fit on.
Old argument names still work and warn.

| Old | New |
|-----|-----|
| `df` (data to score/plot) | `data` |
| `df` / `og_data` (fitting data) | `fit_data` |
| `range_var`         | `x_var`         |
| `desiredCentiles`   | `centiles`      |
| `sim_data_list`     | `sim_grid_list` |

For now, `make_centile_fan()` and its wrappers deliberately keep `df` and `desiredCentiles`;
only `sim_data_list` was renamed there.

## Breaking changes
Breaking changes were restricted to minor unpopular functions and features

* `centile_fan_lifespan()` has been renamed to `centile_fan_brainchart()`. This is
  the only rename without a deprecated alias — existing calls will error.

* `pred_centile()` is no longer part of the public API. It has become the internal
  `.centile_value()`. The old name still works but warns, and will be removed.

* `gamlss2` and `gamlss2charts` moved from `Imports` to `Suggests`. gamlssTools now
  installs and works without them, but code that relied on gamlssTools attaching
  `gamlss2` as a side effect must now install it explicitly:
  `remotes::install_github("gamlss-dev/gamlss2")`. See `?gamlssTools-optional`.

## Minor new features

* `drop1_all()` gains `fit_data`. `gamlss::drop1()` refits each reduced model by
  re-evaluating the model call, which resolves the fitting data *by name* in the
  global environment; passing `fit_data` makes this work from inside functions and
  scripted pipelines, where it previously failed.

* `trajectory_diff()` gains `datafree`, and accepts `df = NULL` when supplied with a
  pre-built grid via `sim_grid_list`.

* `remove_effects()` now warns and returns the data unchanged when neither `rm_terms`
  nor `zero_terms` is supplied, instead of doing unnecessary work.

* New `?gamlssTools-optional` help topic documenting which functions require the
  suggested packages, and how to install them.

## Documentation and internals

* New tests in `tests/testthat/test-datafree.R` covering data-free output against the
  `predictAll()` gold path, data-free eligibility detection, and that every deprecated
  alias still matches its replacement.

* New tests in `tests/testthat/test-list-predictors.R` checking `list_predictors()`
  against edge-cases.

* New tests in `tests/testthat/test-batch-scoring.R` covering `score_centiles()`'s
  `batch_term` pathway: routing of known versus unseen levels, row-order preservation
  when the two are interleaved, all-out-of-sample data, multiple unseen levels,
  `standardize = TRUE`, and whether an unseen study's offset is actually removed.
  The tests that estimate an offset skip unless the dev branch of gamlss2charts is
  installed.

* `make_centile_fan()` documentation substantially expanded, particularly the three
  residualization routes (`remove_point_effect`, `zero_effect`, and `special_term`)
  and when each is appropriate.

* File reorganization: `R/plotting_functions.R` renamed to `R/plotting_funs.R`; new
  `R/datafree_predict.R` and `R/diagnostic_funs.R`.
