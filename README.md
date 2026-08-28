# gamlssTools
This package is intended to make interacting with and plotting GAMLSS models easier. It contains a number of auxiliary functions 
that are compatible with both [gamlss()](https://cran.r-project.org/web/packages/gamlss/index.html) and [gamlss2()](https://github.com/gamlss-dev/gamlss2)

There are 3 vignettes that go over these functions in greater detail: 
- [Plotting Centile Fans](vignettes/centile-fan-plots.Rmd)
- [Model Diagnostics & Centile Scores](vignettes/diagnostics-and-scoring.Rmd)
- [Bootstrapping & Confidence Intervals](vignettes/bootstrapping-and-cis.Rmd)

> **Note:** several functions were renamed in the current version, though the old names still work.
> See [NEWS.md](NEWS.md) for the full list and other recent updates.

### Plotting

`make_centile_fan()` is the primary plotting function. It's designed to cleanly visualize centile fans using ggplot. 
It should be compatible with all gamlss models, regardless of the number of covariates, smooths, random factors, distribution family, etc.

Other plotting/visualization functions include:
- `plot_sigma()`: visualizing sigma on the response scale
- `wp.taki()`: an alternate version of the `gamlss::wp()` function
- `centile_coverage()`: return the cumulative distribution of centiles for a set of datapoints (i.e. what percent of subjects have centile scores <= 50%?)

### Prediction and Scoring

`score_centiles()` returns the centile and/or z-score of each observation in a dataset, given a fitted model.
It works on the model's own fitting data or on new data (currently restricted to data with only levels of
categorical variables that were seen in the model fitting - see `Scoring out-of-sample data` below).

Related functions include `centile_fan_values()` and `sigma_values()` (the values behind the plots),
`pheno_at_centiles()` (the phenotype value at a given centile), and `sim_grid()` (build a prediction grid
once and reuse it across plots).

#### Scoring out-of-sample data (in development)

> ⚠️ **This is a draft implementation.** It requires the `dev` branch of
> [gamlss2charts](https://github.com/andy1764/gamlss2charts/tree/dev):
> `remotes::install_github("andy1764/gamlss2charts@dev")`.

`score_centiles()` can estimate out-of-sample scores (i.e. new levels of a categorical variable that
were not seen in the training data) by passing `batch_term`:

``` r
score_centiles(model, new_data, batch_term = "Study")
```

Rows whose `Study` level was in the fitting data are scored normally. For rows with unseen levels, the level's offset is 
estimated and removed (via `gamlss2charts::predict_score()`) before the centile is computed. 
See the [gamlss2charts](https://github.com/andy1764/gamlss2charts/tree/dev) documentation for additional details.

### Sharing models

`sanitize_gamlss()` strips every per-observation component out of a fitted gamlss model while keeping
data-free prediction intact, so a model can be shared with collaborators who should not be able to
reconstruct the training data. The result still works with `score_centiles()`, `centile_fan_values()`,
`make_centile_fan()` and the rest of the prediction functions.

``` r
clean <- sanitize_gamlss(model, xranges = list(age_days = c(0, 36500)))

audit_gamlss(clean)                       # anything per-observation left?
compare_zscores(model, clean, grid[[1]])  # do predictions still match?
```

`audit_gamlss()` walks the object -- lists, attributes and closure environments -- and reports any
surviving long vector. `compare_zscores()` reports the largest difference, in z units, between the scores the two models
assign to the same grid. Either model can be predicted with its fitting data in scope
(`reference_data` / `comparison_data`) or without it, so the same function also compares two
different models of one outcome, or a model against its own data-free prediction. A few small things survive by necessity (factor levels, `random()` level names, the
fitting `N`, and each smooth's covariate range); see `?sanitize_gamlss` for the full list and check them
by eye before sharing.

### Bootstrapping and Confidence Intervals

Several functions borrow from/adapt from the [gamlss-dev suite](https://github.com/gamlss-dev) to easily fit models to bootstrapped samples
(`bootstrap_gamlss()`), use them to provide confidence intervals on trajectories (i.e. 50th centile curve, sigma; `gamlss_ci()`), and test differences the 
trajectories of different factor levels (`get_median_diffs()`).

Visualization is provided in the wrapper function `plot_centile_cis()`.

### Misc.

The other functions in this package are mostly intended to interact with gamlss model objects. Some highlights include:

- `list_predictors()`: lists all covariates in any moment of a gamlss model
- `get_coeff()`: returns beta coefficient for a specific covariate in a gamlss model
- `age_at_peak()`: returns the value of x (e.g. age) at which the 50th percentile trajectory peaks
- `gamlss_try()`: a wrapper function for fitting gamlss()
- `cohens_f2_local()`: calculate effect size (cohen's fsq) of a covariate using the difference in Rsq of full and nested models
- `drop1_all()`: term-dropping tables across moments, including when the fitting data isn't in the global environment
- `remove_effects()`: residualize (set to mean/median) or zero out selected terms' contributions before plotting

## Installation

You can install the latest version of gamlssTools from [GitHub](https://github.com/BGDlab/gamlssTools) with:

``` r
# install.packages("devtools")
devtools::install_github("BGDlab/gamlssTools", build_vignettes = TRUE) #set build_vignettes to FALSE to save time
```

You can install the development version using: 

``` r
# install.packages("devtools")
devtools::install_github("BGDlab/gamlssTools@dev", build_vignettes = TRUE)
```

### Optional dependencies

Two `Suggests` packages are needed only by a subset of functions rather than by the package as a whole. If needed,
you can install them via:

``` r
remotes::install_github("gamlss-dev/gamlss2")     #needed for all gamlss2 model objects, and for cohens_f2_local()
remotes::install_github("andy1764/gamlss2charts@dev") #dev branch; needed only for scoring out-of-sample centiles
```

Functions that need one of these check for it first and fail with an install hint. See `?gamlssTools-optional`.

## Usage

### Centile Plotting

To plot a basic centile fan using the `iris` dataset:

``` r
library(gamlssTools)

#fit gamlss model
iris_model <- gamlss(formula = Sepal.Width ~ pb(Sepal.Length) + Species, sigma.formula = ~ Sepal.Length, data=iris, family=BCCG)

#basic plot
iris_fan_plot <- make_centile_fan(iris_model, iris, "Sepal.Length", "Species")
print(iris_fan_plot)
```
<img width="3000" height="2100" alt="iris_plot1" src="https://github.com/user-attachments/assets/b844fed7-e736-4a1a-b6c3-ccfadf96bb69" />

You can use all the standard ggplot layers to make your plot prettier.

``` r
iris_fan_plot +
  labs(title="Normative Sepal Width by Length",
  x ="Sepal Length", y = "Sepal Width",
  color = "Species", fill="Species") +
  theme_bw() +
  paletteer::scale_color_paletteer_d("MoMAColors::Smith")
  
```
<img width="3000" height="2100" alt="iris_plot2" src="https://github.com/user-attachments/assets/29fb2e47-e3bd-4b1c-b2f8-c74bb3f6ac52" />

There are also many built-in configuration options, including averaging over categorical variables (like `Species`) 
or changing how centile lines are labeled: 

``` r
make_centile_fan(iris_model, iris, "Sepal.Length", "Species", 
      average_over=TRUE, 
      show_points=FALSE,
      label_centiles="legend") +
      labs(title="Normative Sepal Width by Length",
      x ="Sepal Length", y = "Sepal Width") 
```
<img width="3000" height="2100" alt="iris_plot3" src="https://github.com/user-attachments/assets/aeae149f-6712-4e24-b6d6-0f9e7d1f00b2" />

There are built-in formatting options for the x-axis

``` r
#simulate data
n <- 1000
df <- data.frame(
Age = sample(-140:36525, n, replace = TRUE),
Sex = sample(c("Male", "Female"), n, replace = TRUE),
Study = factor(sample(c("Study_A", "Study_B", "Study_C"), n, replace = TRUE)))

#small study-specific offsets
study_offset <- c(Study_A = -3, Study_B = 1, Study_C = 3)

#log-age in days post-conception (birth = 280 days)
df$logAge <- log(df$Age + 280, base = 10)
df$Pheno <- 6 + 4 * (df$logAge - 2.4) - 1.6 * (df$logAge - 3.6)^2 +
  ifelse(df$Sex == "Male", 1, 0) + ifelse(df$Sex == "Male", .3, 0) * df$logAge +
  unname(study_offset[as.character(df$Study)]) +
  rnorm(n, mean = 0, sd = 1)

#keep the raw scale so new data can be mapped onto the same units later
pheno_raw_range <- range(df$Pheno)
df$Pheno <- scales::rescale(df$Pheno, to = c(1, 50))

#fit a model on log-age and plot a lifespan-style centile fan
pheno_model <- gamlss(formula = Pheno ~ pb(logAge) + Sex + random(Study),
                      sigma.formula = ~ pb(logAge), data = df,
                      family = BCTo)
                      
make_centile_fan(pheno_model, df, x_var="logAge", color_var="Sex",
    label_centiles="legend",
    remove_point_effect="Study",
    x_axis="log_lifespan_fetal",
    point_color_manual = c('Female' = "#F4A15BFF", 'Male' = "#8CB3D1FF"),
    color_manual = c("Female" = "#c05f0d", "Male" = "#38688D"))
```
<img width="3000" height="2100" alt="pheno_plt1" src="https://github.com/user-attachments/assets/aeaa1cca-8530-4610-bbf6-d596697b2889" />

There's also a wrapper function for plotting centile fans in the style of [Bethlehem, Seidlitz, & White et al.](https://www.nature.com/articles/s41586-022-04554-y)

``` r
centile_fan_brainchart(pheno_model, df, x_var="logAge", color_var="Sex")
```
<img width="3000" height="2100" alt="pheno_plt2" src="https://github.com/user-attachments/assets/b51079ca-ec0e-4a37-b365-d0833d8d86d7" />

### Reference Scoring

Individuals' centile scores and standardized "pseudo z-scores" can be calculated using `score_centiles()`.

``` r
scores <- score_centiles(pheno_model, df, standardize=TRUE)
df$centile <- scores$centile
df$std_score <- scores$std_score

#centiles
ggplot(df) +
  geom_histogram(aes(x=centile), bins=20) +
  theme_minimal()

#z-scores  
ggplot(df) +
  geom_histogram(aes(x=std_score), bins=20) +
  theme_minimal()
```
<img width="400" height="400" alt="hist1" src="https://github.com/user-attachments/assets/389d16a4-c2e6-448d-8b21-8523219746a8" />
<img width="400" height="400" alt="hist2" src="https://github.com/user-attachments/assets/1dd69967-4942-41f2-b43e-09c4af520bf4" />


Out-of-sample scoring estimates (i.e. data from new batches) are handled using [gamlss2charts](https://github.com/andy1764/gamlss2charts).

``` r
#simulate new data
n <- 500
new_df <- data.frame(
Age = sample(10000:36525, n, replace = TRUE),
Sex = sample(c("Male", "Female"), n, replace = TRUE))
new_df$Study <- factor("NEW")

#log-age in days post-conception (birth = 280 days)
new_df$logAge <- log(new_df$Age + 280, base = 10)
new_df$Pheno <- 6 + 4 * (new_df$logAge - 2.4) - 1.6 * (new_df$logAge - 3.6)^2 +
  ifelse(new_df$Sex == "Male", 1, 0) + ifelse(new_df$Sex == "Male", .3, 0) * new_df$logAge -
  2 + rnorm(n, mean = 0, sd = 1)
new_df$Pheno <- scales::rescale(new_df$Pheno, to = c(1, 50), from = pheno_raw_range)

#score
new_df$centile <- score_centiles(pheno_model, new_df, batch_term = "Study")

ggplot(new_df) +
  geom_histogram(aes(x=centile), bins=20) +
  theme_minimal()
```
<img width="400" height="400" alt="hist3" src="https://github.com/user-attachments/assets/b610e8ec-bf81-47e1-ace3-080a842c9716" />

