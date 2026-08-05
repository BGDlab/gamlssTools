# Robust GAMLSS Simulation and Analysis Script
# This script simulates multiple datasets with male and female subjects of varying ages,
# fits GAMLSS models, and visualizes median trajectories and differences

# Load required libraries
library(gamlss2)
library(gamlss)
library(gamlssTools)
library(ggplot2)
library(dplyr)
library(scales)
library(foreach)
library(rlang)

#source updated functions
source("~/Desktop/BGD_Repos/gamlssTools/R/bootstrapping.R")
source("~/Desktop/BGD_Repos/gamlssTools/R/centile_fans.R")
source("~/Desktop/BGD_Repos/gamlssTools/R/plotting_functions.R")
source("~/Desktop/BGD_Repos/gamlssTools/R/centile_fan_helper_funs.R")


# Set seed for reproducibility
set.seed(123)

# Function to simulate a single dataframe
simulate_dataframe <- function(n_subjects = 1000, age_range = c(18, 80), 
                               male_effect = 2, female_effect = 0, 
                               age_effect = 0.5, interaction_effect = 0.1,
                               noise_sd = 1, dataset_id = 1) {
  
  # Generate ages
  ages <- sample(age_range[1]:age_range[2], n_subjects, replace = TRUE)
  
  # Generate sex (50-50 split)
  sex <- sample(c("Male", "Female"), n_subjects, replace = TRUE)
  
  # Create a phenotype that varies by sex and age
  phenotype <- case_when(
    sex == "Male" ~ male_effect + age_effect * (ages - 40) + interaction_effect * (ages - 40) + rnorm(n_subjects, 0, noise_sd),
    sex == "Female" ~ female_effect + age_effect * (ages - 40) + rnorm(n_subjects, 0, noise_sd)
  )
  
  # Ensure phenotype is positive
  phenotype <- pmax(phenotype, 0.1)
  
  # Create dataframe
  df <- data.frame(
    Age = ages,
    Sex = factor(sex, levels = c("Male", "Female")),
    Phenotype = phenotype,
    Dataset = paste0("Dataset_", dataset_id)
  )
  
  return(df)
}

# Function to fit GAMLSS model with better error handling
fit_gamlss_model <- function(df) {
  # Try BCCG first
  model <- tryCatch({
    gamlss2(Phenotype ~ pb(Age) + Sex + pb(Age):Sex, 
           sigma.formula = ~ pb(Age) + Sex,
           data = df, 
           family = BCCGo,
           method = "RS")
  }, error = function(e) {
    cat("BCCG failed, trying CG\n")
    # Try NO (normal distribution)
    tryCatch({
      gamlss2(Phenotype ~ pb(Age) + Sex + pb(Age):Sex, 
             sigma.formula = ~ pb(Age) + Sex,
             data = df, 
             family = BCCGo,
             method = "CG")
    }, error = function(e2) {
      cat("NO failed, trying simpler model\n")
      # Try even simpler model
      tryCatch({
        gamlss2(Phenotype ~ pb(Age) + Sex, 
               sigma.formula = ~ pb(Age),
               data = df, 
               family = BCCGo,
               method = "RS")
      }, error = function(e3) {
        cat("All models failed for dataset:", df$Dataset[1], "\n")
        return(NULL)
      })
    })
  })
  
  return(model)
}

# Function to create visualizations with error handling
create_visualizations <- function(model, df, dataset_name) {
  if (is.null(model)) {
    cat("Skipping visualizations for", dataset_name, "- model fitting failed\n")
    return(NULL)
  }
  
  cat("Creating visualizations for", dataset_name, "\n")
  
  # Create plots directory if it doesn't exist
  if (!dir.exists("plots")) {
    dir.create("plots")
  }
  
  # plot datapoints
  plt_data <- ggplot(df) +
    geom_point(aes(x=Age, y=Phenotype, color=Sex)) +
    theme_minimal()
  
  ggsave(paste0("plots/", dataset_name, "_datapoints.png"), 
         plot = plt_data, width = 10, height = 6, dpi = 300)
  
  # bootstrap models once
  boot_mod_list <- bootstrap_gamlss(model, df=df, B=500, type="resample", stratify=TRUE, group_var="Sex")
  
  # 1. Plot pointwise centile CIs for median trajectories
  tryCatch({
    cat("  Creating median trajectories plot...\n")
    p1 <- plot_centile_cis(
      gamlssModel = model,
      df = df,
      x_var = "Age",
      color_var = "Sex",
      desiredCentiles = c(0.5),
      interval = 0.95,
      ci_type = "pointwise",
      boot_list = boot_mod_list
    )
    
    # Add title and save
    p1 <- p1 + 
      labs(title = paste("Median Trajectories by Sex -", dataset_name),
           subtitle = "50th percentile with 95% pointwise confidence intervals") +
      theme_minimal()
    
    ggsave(paste0("plots/", dataset_name, "_ptwise_median_trajectories.png"), 
           plot = p1, width = 10, height = 6, dpi = 300)
    
    cat("  Saved median trajectories plot\n")
    
  }, error = function(e) {
    cat("  Error creating median trajectories plot:", e$message, "\n")
  })
  
  # 1. Plot simultaneous centile CIs for median trajectories
  tryCatch({
    cat("  Creating median trajectories plot...\n")
    p1 <- plot_centile_cis(
      gamlssModel = model,
      df = df,
      x_var = "Age",
      color_var = "Sex",
      desiredCentiles = c(0.5),
      interval = 0.95,
      ci_type = "simultaneous",
      boot_list = boot_mod_list
    )
    
    # Add title and save
    p1 <- p1 + 
      labs(title = paste("Median Trajectories by Sex -", dataset_name),
           subtitle = "50th percentile with 95% simultaneous confidence intervals") +
      theme_minimal()
    
    ggsave(paste0("plots/", dataset_name, "_sim_median_trajectories.png"), 
           plot = p1, width = 10, height = 6, dpi = 300)
    
    cat("  Saved median trajectories plot\n")
    
  }, error = function(e) {
    cat("  Error creating median trajectories plot:", e$message, "\n")
  })
  
  # 2. Calculate and plot median differences
  tryCatch({
    cat("  Creating sex differences plot...\n")
    diffs <- get_median_diffs(
      gamlssModel = model,
      df = df,
      x_var = "Age",
      factor_var = "Sex",
      boot_list = boot_mod_list,
      ci_type = "pointwise",
      moment = "mu"
    )
    
    cat(names(diffs))
    col_name <- names(diffs)[grep("minus", names(diffs))]
    cat(col_name)
    
    # Create difference plot
    p2 <- ggplot(diffs) +
      geom_line(aes(x = Age, y = !!sym(col_name)), size = 1) +
      geom_ribbon(aes(x = Age, ymin = lower, ymax = upper), alpha = 0.3) +
      labs(title = paste("Sex Differences in Median Trajectories -", dataset_name),
           subtitle = "difference with 95% pointwise confidence intervals",
           x = "Age", 
           y = paste("Difference", col_name)) +
      theme_minimal()
    
    ggsave(paste0("plots/", dataset_name, "_ptwise_sex_differences.png"), 
           plot = p2, width = 10, height = 6, dpi = 300)
    
    cat("  Saved sex differences plot\n")
    
  }, error = function(e) {
    cat("  Error creating sex differences plot:", e$message, "\n")
  })
  
  tryCatch({
    cat("  Creating sex differences plot...\n")
    diffs <- get_median_diffs(
      gamlssModel = model,
      df = df,
      x_var = "Age",
      factor_var = "Sex",
      boot_list = boot_mod_list,
      ci_type = "simultaneous",
      moment = "mu"
    )
    
    cat(names(diffs))
    col_name <- names(diffs)[grep("minus", names(diffs))]
    cat(col_name)
    
    # Create difference plot
    p2 <- ggplot(diffs) +
      geom_line(aes(x = Age, y = !!sym(col_name)), size = 1) +
      geom_ribbon(aes(x = Age, ymin = lower, ymax = upper), alpha = 0.3) +
      labs(title = paste("Sex Differences in Median Trajectories -", dataset_name),
           subtitle = "difference with 95% simultaneous confidence intervals",
           x = "Age", 
           y = paste("Difference", col_name)) +
      theme_minimal()
    
    ggsave(paste0("plots/", dataset_name, "_sim_sex_differences.png"), 
           plot = p2, width = 10, height = 6, dpi = 300)
    
    cat("  Saved sex differences plot\n")
    
  }, error = function(e) {
    cat("  Error creating sex differences plot:", e$message, "\n")
  })
}

# Main simulation and analysis
main_analysis <- function() {
  cat("Starting Robust GAMLSS simulation and analysis\n")
  cat("============================================\n\n")
  
  # Define simulation parameters - same sample size but different effects
  n_subjects <- 1000  # Same sample size for all datasets
  
  # Define different effect scenarios
  effect_scenarios <- list(
    list(
      name = "Small_Sex_Effects",
      male_effect = 1.5,
      female_effect = 1.0,
      age_effect = 0.3,
      interaction_effect = 0.05,
      description = "Small differences between sexes"
    ),
    list(
      name = "Large_Sex_Effects", 
      male_effect = 3.0,
      female_effect = 0.5,
      age_effect = 0.4,
      interaction_effect = 0.15,
      description = "Large differences between sexes"
    ),
    list(
      name = "Age_Dominant",
      male_effect = 2.0,
      female_effect = 2.0,
      age_effect = 0.8,
      interaction_effect = 0.02,
      description = "Strong age effects, minimal sex differences"
    ),
    list(
      name = "Strong_Interaction",
      male_effect = 1.0,
      female_effect = 1.0,
      age_effect = 0.2,
      interaction_effect = 0.25,
      description = "Strong sex-by-age interaction"
    ),
    list(
      name = "No_Effects",
      male_effect = 1.0,
      female_effect = 1.0,
      age_effect = 0.1,
      interaction_effect = 0.01,
      description = "Minimal effects (null scenario)"
    )
  )
  
  # Storage for results
  all_models <- list()
  all_data <- list()
  
  # Simulate datasets and fit models
  for (i in seq_along(effect_scenarios)) {
    scenario <- effect_scenarios[[i]]
    dataset_name <- scenario$name
    cat("Processing", dataset_name, ":", scenario$description, "\n")
    cat("  Sample size:", n_subjects, "subjects\n")
    cat("  Male effect:", scenario$male_effect, "\n")
    cat("  Female effect:", scenario$female_effect, "\n")
    cat("  Age effect:", scenario$age_effect, "\n")
    cat("  Interaction effect:", scenario$interaction_effect, "\n")
    
    # Simulate data with specific parameters
    df <- simulate_dataframe(
      n_subjects = n_subjects,
      male_effect = scenario$male_effect,
      female_effect = scenario$female_effect,
      age_effect = scenario$age_effect,
      interaction_effect = scenario$interaction_effect,
      dataset_id = i
    )
    
    # Store data
    all_data[[dataset_name]] <- df
    
    # Fit GAMLSS model
    cat("  Fitting GAMLSS model...\n")
    model <- fit_gamlss_model(df)
    all_models[[dataset_name]] <- model
    
    # Create visualizations
    create_visualizations(model, df, dataset_name)
    
    cat("  Completed", dataset_name, "\n\n")
  }
  
  cat("Analysis complete! Check the 'plots' directory for all visualizations.\n")
  
  return(list(models = all_models, data = all_data))
}

# Run the analysis
cat("Running robust simulation analysis...\n")
results <- main_analysis()

