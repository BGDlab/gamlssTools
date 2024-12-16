# gamlssTools
This package is intended to make interacting with and plotting gamlss models easier. It contains a number of auxiliary functions 
that will be useful for those using [gamlss()](https://cran.r-project.org/web/packages/gamlss/index.html). 

### Plotting

`make_centile_fan()` is the primary plotting function. It's designed to cleanly visualize centile fans using ggplot. It should be compatible
with all gamlss models, regardless of the number of covariates, smooths, random factors, distribution family, etc.

### Extracting Information

The other auxiliary functions in this package are mostly intended to interact with gamlss model objects. Some hilights include:

- `get_coeff()`: returns beta coefficient for a specific covariate in a gamlss model
- `list_predictors()`: lists all covariates in any moment of a gamlss model
- `pred_og_centile()`: Returns the centile and/or z-score values for the datapoints used to fit a gamlss model

## Installation
You can install the development version of gamlssTools from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("BGDlab/gamlssTools")
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
![average_iris_plot](https://github.com/user-attachments/assets/cea86418-3da5-4e63-a64c-f35f2f3e9f3a)

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

![39e44a98-e8cb-4990-9632-79134e5c7f0b](https://github.com/user-attachments/assets/6f048350-c750-49f7-b062-002b34578779)


