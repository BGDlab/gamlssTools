# gamlssTools
This package is intended to make interacting with and plotting gamlss models easier. It contains a number of auxiliary functions 
that will be useful for those using [gamlss()](https://cran.r-project.org/web/packages/gamlss/index.html). 

### Plotting

`make_centile_fan()` is the primary plotting function. It's designed to cleanly visualize centile fans using ggplot. It should be compatible
with all gamlss models, regardless of the number of covariates, smooths, random factors, distribution family, etc.

### Extracting Information

The other auxilliary functions in this package are mostly intended to interact with gamlss model objects. Some hilights include:

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
![clean_iris_plot](https://github.com/user-attachments/assets/f3b5b473-d302-4752-86ab-6f905baae35b)


You can use ggplot to clean up the axes a bit
```
iris_fan_plot +
  labs(title="Normative Sepal Width by Length",
  x ="Sepal Length", y = "Sepal Width",
  color = "Species", fill="Species")
```
![basic_iris_plot](https://github.com/user-attachments/assets/9ae4e535-94b9-4c7d-a0d5-c13331808d81)


There are also many built-in configuration options,  including averaging over categorical variables (like `Species`): 
```
make_centile_fan(iris_model, iris, "Sepal.Length", "Species", average_over=TRUE)
```
![average_iris_plot](https://github.com/user-attachments/assets/cea86418-3da5-4e63-a64c-f35f2f3e9f3a)
