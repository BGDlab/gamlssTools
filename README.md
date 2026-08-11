# gamlssTools
This package is intended to make interacting with and plotting GAMLSS models easier. It contains a number of auxiliary functions 
that are compatible with both [gamlss()](https://cran.r-project.org/web/packages/gamlss/index.html) and [gamlss2()](https://github.com/gamlss-dev/gamlss2)

There are 3 vignettes that go over these functions in greater detail: 
- *Plotting Centile Fans*
- *Model Diagnostics & Scoring*
- *Bootstrapping & Confidence Intervals*.

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

> ⚠️ **This is a draft implementation.** It requires the development version of
> [gamlss2charts](https://github.com/andy1764/gamlss2charts), pending merge approval.

`score_centiles()` can estimate out-of-sample scores (i.e. new levels of a categorical variable that
were not seen in the training data) by passing `batch_term`:

``` r
score_centiles(model, new_data, batch_term = "Study")
```

Rows whose `Study` level was in the fitting data are scored normally. For rows with unseen levels, the level's offset is 
estimated and removed (via `gamlss2charts::predict_score()`) before the centile is computed. 
See the [gamlss2charts](https://github.com/andy1764/gamlss2charts) documentation for additional details.

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
You can install the development version of gamlssTools from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("BGDlab/gamlssTools", build_vignettes = TRUE) #set build_vignettes to FALSE to save time
```

### Optional dependencies

Two `Suggests` packages are needed only by a subset of functions rather than by the package as a whole. If needed,
you can install them via:

``` r
remotes::install_github("gamlss-dev/gamlss2")     #needed for all gamlss2 model objects, and for cohens_f2_local()
remotes::install_github("andy1764/gamlss2charts") #needed only scoring out-of-sample centiles
```

Functions that need one of these check for it first and fail with an install hint. See `?gamlssTools-optional`.

## Usage

### Centile Plotting
To plot a basic centile fan using the `iris` dataset:
```
library(gamlssTools)

#fit gamlss model
iris_model <- gamlss(formula = Sepal.Width ~ pb(Sepal.Length) + Species, sigma.formula = ~ Sepal.Length, data=iris, family=BCCG)

#basic plot
iris_fan_plot <- make_centile_fan(iris_model, iris, "Sepal.Length", "Species")
print(iris_fan_plot)
```
![basic_iris_plot](https://github.com/user-attachments/assets/9ae4e535-94b9-4c7d-a0d5-c13331808d81)

You can use all the standard ggplot layers to make your plot prettier.
```
iris_fan_plot +
  labs(title="Normative Sepal Width by Length",
  x ="Sepal Length", y = "Sepal Width",
  color = "Species", fill="Species") +
  theme_bw() +
  paletteer::scale_color_paletteer_d("MoMAColors::Smith")
  
```
![707785a4-43fa-4358-8527-480a5a7f2608](https://github.com/user-attachments/assets/704806df-bd43-4fb5-963b-e9338a94c855)

There are also many built-in configuration options, including averaging over categorical variables (like `Species`) 
or changing how centile lines are labeled: 
```
make_centile_fan(iris_model, iris, "Sepal.Length", "Species", 
      average_over=TRUE, 
      show_points=FALSE,
      label_centiles="legend",
      ) +
      labs(title="Normative Sepal Width by Length",
      x ="Sepal Length", y = "Sepal Width") 
```
![39e44a98-e8cb-4990-9632-79134e5c7f0b](https://github.com/user-attachments/assets/6f048350-c750-49f7-b062-002b34578779)

There are built-in formatting options for the x-axis

```
#simulate data
n <- 1000
df <- data.frame(
Age = sample(-280:36525, n, replace = TRUE),
Sex = sample(c("Male", "Female"), n, replace = TRUE),
Study = factor(sample(c("Study_A", "Study_B", "Study_C"), n, replace = TRUE)))

#log-age in days post-conception (birth = 280 days)
df$logAge <- log(df$Age + 280, base = 10)
df$Pheno <- 6 + 4 * (df$logAge - 2.4) - 1.6 * (df$logAge - 3.6)^2 +
  ifelse(df$Sex == "Male", 1, 0) + ifelse(df$Sex == "Male", .3, 0) * df$logAge + 
  rnorm(n, mean = 0, sd = 1)
df$Pheno <- scales::rescale(df$Pheno, to = c(1, 50))

#fit a model on log-age and plot a lifespan-style centile fan
pheno_model <- gamlss(formula = Pheno ~ pb(logAge) + Sex + random(Study), sigma.formula= ~ pb(logAge), data = df, family=BCCG)

make_centile_fan(pheno_model, df, x_var="logAge", color_var="Sex",
    label_centiles="legend",
    x_axis="log_lifespan_fetal",
    point_color_manual = c('Female' = "#F4A15BFF", 'Male' = "#8CB3D1FF"),
    color_manual = c("Female" = "#c05f0d", "Male" = "#38688D"))
```

There's also a wrapper function for plotting centile fans in the style of [Bethlehem, Seidlitz, & White et al.](https://www.nature.com/articles/s41586-022-04554-y)

```
centile_fan_brainchart(pheno_model, df, x_var="logAge", color_var="Sex")
```

### Reference Scoring
Section pending!


