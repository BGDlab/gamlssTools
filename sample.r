library(gamlss)    # fit the model and predict centiles
library(ggplot2)   # plotting

# Load some sample data
fn <- "/Users/youngjm/Projects/GAMLSS/data.csv"
sampleData <- read.csv(fn)

# There are four columns: age_in_years, phenotype, SurfaceHoles, and sex
print(colnames(sampleData))

# Sex needs to be a factor
sampleData$sex <- as.factor(sampleData$sex)

# For Lifespan and CLIP, we used the log of post conception age in days. 
sampleData$logAgeDays <- log((sampleData$age_in_years*365.25 + 280), base=10)

# Set up data frames for model predictions
minAge <- min(sampleData$logAgeDays)
maxAge <- max(sampleData$logAgeDays)
ageRange <- seq(minAge, maxAge, 0.1)  # generate an age range with increments of 0.1 

# Set up a list of tick marks to use on log(post-conception age) x-axes
tickMarks <- c()
for (year in c(0, 1, 2, 5, 10, 20)){ # years
  tickMarks <- append(tickMarks, log(year*365.25 + 280, base=10))
}
tickLabels <- c("Birth", "1", "2", "5", "10", "20")

# Use the median SurfaceHoles value
# Make one set of data for each sex
dataToPredictM <- data.frame(logAgeDays=ageRange,
                             SurfaceHoles=c(rep(median(sampleData$SurfaceHoles), length(ageRange))),
                             sex=c(rep(as.factor("M"), length(ageRange))))
dataToPredictF <- data.frame(logAgeDays=ageRange,
                             SurfaceHoles=c(rep(median(sampleData$SurfaceHoles), length(ageRange))),
                             sex=c(rep(as.factor("F"), length(ageRange))))


# Part 1: create the models and train them on the desired data
# Generate GAMLSS model
formula <- as.formula("phenotype ~ fp(logAgeDays) + SurfaceHoles + sex")    # smooth the age
gamlssModel <- gamlss(formula = formula,
                      sigma.formula = formula,
                      nu.formula = as.formula("phenotype ~ 1"),
                      family = GG,
                      data = na.omit(sampleData),       # drop any missing data
                      control = gamlss.control(n.cyc = 200), # from Lifespan
                      trace = F)
print("Finished training the model")

# Predict phenotype values in a set age range
predictedModelM <- predictAll(gamlssModel, newdata=dataToPredictM)
predictedModelF <- predictAll(gamlssModel, newdata=dataToPredictF)

# Part 2: Predict the median centile for the M and F models using the parameters from each model
predictedPhenotypeM <- qGG(c(0.5),
                           mu=predictedModelM$mu,
                           sigma=predictedModelM$sigma,
                           nu=predictedModelM$nu)

predictedPhenotypeF <- qGG(c(0.5),
                           mu=predictedModelF$mu,
                           sigma=predictedModelF$sigma,
                           nu=predictedModelF$nu)
# Average the two median centiles together to get the overall median
medianPhenotype <- (predictedPhenotypeM + predictedPhenotypeF)/2

# Plot the original data and the predicted centiles
medianPlot <- ggplot() +
    geom_point(data=sampleData, aes(x=logAgeDays, y=phenotype, color=sex), alpha=0.4) +
    geom_line(aes(ageRange, medianPhenotype, linetype="Predicted Centile")) +
    scale_x_continuous(breaks=tickMarks, labels=tickLabels) +
    xlab("Age at scan (log(years)") +
    ylab("Phenotype Value")
print(medianPlot) 

# Part 3: Create a full growth chart (aka fan plot)
# Start by setting up the list of centiles for the fan plot
fanCentiles <- c()
desiredCentiles <- c(0.004, 0.02, 0.1, 0.25, 0.5, 0.75, 0.9, 0.98, 0.996)

# For each desired centile
for (i in c(1:length(desiredCentiles))){
    predictionsM <- qGG(desiredCentiles[[i]],
                        mu=predictedModelM$mu,
                        sigma=predictedModelM$sigma,
                        nu=predictedModelM$nu)

    predictionsF <- qGG(desiredCentiles[[i]],
                        mu=predictedModelF$mu,
                        sigma=predictedModelF$sigma,
                        nu=predictedModelF$nu)

    fanCentiles[[i]] <- (predictionsM + predictionsF)/2
}

# Make the fan plot
sampleCentileFan <- ggplot() +
  geom_point(aes(x=sampleData$logAgeDays, sampleData$phenotype, color=sampleData$sex), alpha=0.5) +
  geom_line(aes(x=ageRange, y=fanCentiles[[1]]), alpha=0.2) +
  geom_line(aes(x=ageRange, y=fanCentiles[[2]]), alpha=0.4) +
  geom_line(aes(x=ageRange, y=fanCentiles[[3]]), alpha=0.6) +
  geom_line(aes(x=ageRange, y=fanCentiles[[4]]), alpha=0.8) +
  geom_line(aes(x=ageRange, y=fanCentiles[[5]])) +
  geom_line(aes(x=ageRange, y=fanCentiles[[6]]), alpha=0.8) +
  geom_line(aes(x=ageRange, y=fanCentiles[[7]]), alpha=0.6) +
  geom_line(aes(x=ageRange, y=fanCentiles[[8]]), alpha=0.4) +
  geom_line(aes(x=ageRange, y=fanCentiles[[9]]), alpha=0.2) +
  scale_x_continuous(breaks=tickMarks, labels=tickLabels, 
                     limits=c(tickMarks[[1]], max(sampleData$logAgeDays))) +
  labs(title="Sample Centile Growth Chart") + 
  xlab("Age at Scan (log(years))") +
  ylab("Phenotype Value")

print(sampleCentileFan)
