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

* Data-free prediction is now the default for `gamlss` fits in `score_centiles()`,
  `centile_fan_values()`, `sigma_values()`, `remove_effects()`, `make_centile_fan()`,
  `plot_sigma()` and `trajectory_diff()`; supplying `fit_data` (the original fitting 
  data) forces the exact `gamlss::predictAll()` path instead.

* Models containing a smoother that cannot be reconstructed this way (`cs()`, `ps()`,
  `ga()`, `s()`) and/or a parametric term whose columns are computed from the data as a
  whole (`poly()`, `ns()`, `bs()`, `scale()`, `cut()`) are detected up front and ask for `fit_data`.

* `gamlss2` fits continue to use gamlss2's own `predict()` and are unaffected.

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
  environments -- and reports every atomic vector of length exactly `n`. 
  `pb()` interpolation functions are checked for an evenly spaced grid rather than 
  by length, so a rebuilt smooth is not confused with surviving covariate data 
  (and vice versa) when `grid_n` happens to match the sample size.

## Diagnostics

* `centile_coverage()` gains `centiles`, so pre-calculated centiles can be passed in
  instead of re-scoring the model. Accepts a numeric vector, the name of a column of
  `data`, or the dataframe returned by `score_centiles(standardize=TRUE)`. `gamlssModel`
  and `data` are now optional: supply either a model (to score `data`) or `centiles`,
  and `data` is only needed when grouping with `group`/`interval_var`.

* `centile_coverage()` also gains `batch_term`, passed through to `score_centiles()` so coverage can
  be checked on data containing unseen levels of a site/study/batch variable. It only applies when
  scoring from `gamlssModel` - combining it with `centiles` is an error.
  
* added `compare_scores()` to compare the predicted z-scores of a dataframe and/or predicted
  centile values across two models (or the same model predicted datafree and/or 'sanitized' for)
  sharing.

## Out-of-Sample Prediction

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

## Plotting

* New `compare_centile_fans()`, a wrapper around `make_centile_fan()` that overlays the
  centile fans of two models for visual comparison. It calls `make_centile_fan()` once per
  model and transplants the second model's centile (and peak) layers onto the first
  model's plot, so points, x-axis formatting and centile labels are drawn once and every
  shared argument behaves as it does in `make_centile_fan()`. Color carries the model:
  both fans are drawn solid in the two `model_colors`, at `alpha` (0.6 by default) so the
  overlap stays visible. That frees the grouping variable to be shown as facets, one panel
  per level -- so it is named `facet_var` here rather than `color_var` -- with the data
  points, if `show_points = TRUE`, drawn once in a neutral `point_color`. Pairs with
  `compare_scores()`, which quantifies the same difference.

* Fixed `make_centile_fan()` labelling the x-axis `point_df[[x_var]]` instead of naming
  the variable. The averaged branch (`average_over = TRUE`, or `color_var = NULL`) builds
  its layers from vectors rather than a data/mapping pair, so ggplot derived the literal
  expression as the label. The `x_axis` presets, which name the axis themselves via
  `format_x_axis()`, are unaffected.

* Fixed `sim_grid()` erroring on a factor column when no `factor_var` is given. A
  leftover reference to the old `df` argument name resolved to `stats::df` instead of
  the data ("object of type 'closure' is not subsettable"), which also broke
  `make_centile_fan()` and its wrappers whenever `color_var = NULL`.

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

* New tests in `tests/testthat/test-axis-labels.R` pinning both halves of the axis-label
  fix: the label is repaired on `"custom"` axes and left alone on the `x_axis` presets.

* New tests in `tests/testthat/test-compare-centile-fans.R` covering the layer structure
  `compare_centile_fans()` produces: one solid fan per model in its own color and alpha,
  with distinct predictions; no leftover mapped color to rival the model scale; a two-key
  model legend; one panel per level of `facet_var`, and none when there is nothing to
  facet by; points and centile labels drawn only once; and a single `sim_grid_list`
  being reused for the second model.

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
