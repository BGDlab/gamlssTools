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

* `gamlss2` fits continue to use gamlss2's own `predict()` and are unaffected.

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
