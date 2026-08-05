# gamlssTools
This package is intended to make interacting with and plotting GAMLSS models easier. It contains a number of auxiliary functions 
that will be useful for those using [gamlss()](https://cran.r-project.org/web/packages/gamlss/index.html) or [gamlss2()](https://github.com/gamlss-dev/gamlss2)

There are two vignettes that go over these functions in greater detail.

### Plotting

`make_centile_fan()` is the primary plotting function. It's designed to cleanly visualize centile fans using ggplot. 
It should be compatible with all gamlss models, regardless of the number of covariates, smooths, random factors, distribution family, etc.

Other plotting/visualization functions include:
- `plot_sigma()`: visualizing sigma on the response scale
- `wp.taki()`: an alternate version of the `gamlss::wp()` function
- `cent_cdf()`: return the cumulative distribution of centiles for a set of datapoints (i.e. what percent of subjects have centile scores <= 50%?)


### Bootstrapping and Confidence Intervals

Several functions borrow from/adapt from the [gamlss-dev suite](https://github.com/gamlss-dev) to easily fit models to bootstrapped samples
(`bootstrap_gamlss()`), use them to provide confidence intervals on trajectories (i.e. 50th centile curve, sigma; `gamlss_ci()`), and test differences the 
trajectories of different factor levels (`get_median_diffs()`).

Visualization is provided in the wrapper function `plot_centile_cis()`.

### Misc.

The other functions in this package are mostly intended to interact with gamlss model objects. Some highlights include:

- `get_coeff()`: returns beta coefficient for a specific covariate in a gamlss model
- `list_predictors()`: lists all covariates in any moment of a gamlss model
- `pred_og_centile()`: returns the centile and/or z-score values for the datapoints used to fit a gamlss model
- `age_at_peak()`: returns the value of x (e.g. age) at which the 50th percentile trajectory peaks
- `gamlss_try()`: a wrapper function for fitting gamlss()
- `cohens_f2_local()`: calculate effect size (cohen's fsq) of a covariate using the difference in Rsq of full and nested models

## Installation
You can install the development version of gamlssTools from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("BGDlab/gamlssTools", build_vignettes = TRUE) #set build_vignettes to FALSE to save time
```

## Usage

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

There are also many built-in configuration options,  including averaging over categorical variables (like `Species`): 
```
make_centile_fan(iris_model, iris, "Sepal.Length", "Species", average_over=TRUE, color_manual="#467326FF")
```
![39e44a98-e8cb-4990-9632-79134e5c7f0b](https://github.com/user-attachments/assets/6f048350-c750-49f7-b062-002b34578779)

You can also extract and manipulate the layers created by `make_centile_fan()`. For instance, if you're fitting separate models for each sex but still want to plot them on the same plot:

```
#simulate data to plot
df <- data.frame(
  Age = sample(0:36525, 10000, replace = TRUE),
  Sex = sample(c("Male", "Female"), 10000, replace = TRUE),
  Study = factor(sample(c("Study_A", "Study_B", "Study_C"), 10000, replace = TRUE)))

 df$Pheno <- ((df$Age)/365)^3 + rnorm(10000, mean = 0, sd = 100000)
 df$Pheno <- scales::rescale(df$Pheno, to = c(1, 10))

#fit separate models for each sex
df_m <- df %>% filter(Sex=="Male")
df_f <- df %>% filter(Sex=="Female")

m_model <- gamlss(formula = Pheno ~ pb(Age) + random(Study), sigma.formula= ~ pb(Age), data = df_m, family=BCCG)
f_model <- gamlss(formula = Pheno ~ pb(Age) + random(Study), sigma.formula= ~ pb(Age), data = df_f, family=BCCG)

#plot by sex
plot_m <- make_centile_fan(m_model, df_m, "Age", show_points=FALSE, x_axis="lifespan", color_manual="orange")
plot_f <- make_centile_fan(f_model, df_f, "Age", show_points=FALSE, x_axis="lifespan", color_manual="green")

plot_combined <- plot_m
plot_combined$layers <- c(plot_m$layers, plot_f$layers)

plot_combined +
 theme_bw() +
 ggtitle("Male and Female Growth Charts")
```
![2b953588-7b75-4106-a683-f4946d2d9021](https://github.com/user-attachments/assets/cd38b611-e15f-4b23-b49e-4fa0c85f0abd)

