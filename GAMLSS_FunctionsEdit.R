
##feature 
phenotype <- data.frame(c("VentricleVolume","CerebralWhiteMatterVol","CortexVol","SubCortGrayVol" ))


createMainGamm <- function(df, measure) {
  
  temp <- transform(df, age_in_years2 = age_in_years^2, age_in_years3 = age_in_years^3)
  formula <- as.formula(paste(measure, "age_in_years+age_in_years2+age_in_years3", sep="~"))
  
  # Make a basic linear model accounting for age and surface holes for the given data frame
  model <- gamlss(formula, data=na.omit(df), family = GG)
  # Return the model
  return(model)
}

createMainGamm(analysisDf,analysisDf$VentricleVolume)

df <- transform(analysisDf, age_in_years2 = age_in_years^2, age_in_years3 = age_in_years^3)

cubVen <-gamlss(VentricleVolume~age_in_years+age_in_years2+age_in_years3, data=na.omit(df), family = GG)
plot(VentricleVolume~age_in_years, col = "lightblue", data = df)
lines(fitted(cubVen)[order(df$age_in_years)]~df$age_in_years[order(df$age_in_years)])

library(gamlss)

createMainGamm(analysisDf, analysisDf$VentricleVolume) 







gamlssFun <- function(df, measure) {
  formula <- as.formula(paste(measure, "age_in_years+I(age_in_years^2)+I(age_in_years^3)", sep="~"))
  
  # Make a basic linear model accounting for age and surface holes for the given data frame
  model <- gamlss(formula, data=na.omit(df), family = GG)
  
 # aModel <- gamlss(VentricleVolume~age_in_years+I(age_in_years^2)+I(age_in_years^3), data = analysisDf, family = NO)
  
  # Return the model
  return(model)
}


gamlssFun(analysisDf,analysisDf$VentricleVolume)























#####functions
gamlss_linear <- function(df, measure) {
  formula <- as.formula(paste(measure, "age_in_years", sep="~"))
  model <- gamlss(formula, data=na.omit(df), family = GG)
  # Return the model
  return(model)
}



gamlss_quad <- function(df, measure) {
  formula <- as.formula(paste(measure, "age_in_years+I(age_in_years^2)+I(age_in_years^3)", sep="~"))
  model <- gamlss(formula, data=na.omit(df), family = GG)
  # Return the model
  return(model)
}



gamlss.gam <- function(df, measure) {
  formula <- as.formula(paste(measure, "pb(age_in_years)+sex+SurfaceHoles", sep="~"))
  model <- gamlss(formula, family = GA, data=na.omit(analysisDf), trace = FALSE)
  # Return the model
  return(model)
}




#### In practice 


phenotype <- (c("VentricleVolume","CerebralWhiteMatterVol","CortexVol","SubCortGrayVol" ))
phenotype[4]



#linear plot
par(mfrow=c(2,2))

for (i in 1:length(phenotype)){
feature <- phenotype[i]
modeltest <- gamlss_linear(analysisDf,feature)
formula1 <- as.formula(paste(feature, "age_in_years", sep="~"))
plot(formula1,data = analysisDf,col = "lightblue")
lines(fitted(modeltest)~(analysisDf$age_in_years))

}

dev.off()



#quadratic plot
par(mfrow=c(2,2))


set <- transform(analysisDf, age_in_years2 = age_in_years^2, age_in_years3 = age_in_years^3)#the function does not use the set polynomial, the graph does
for (i in 1:length(phenotype)){
  feature <- phenotype[i]
  modeltest <- gamlss_quad(analysisDf,feature)
  formula1 <- as.formula(paste(feature, "age_in_years", sep="~"))
  plot(formula1,data = analysisDf,col = "lightblue")
  lines(fitted(modeltest)[order(set$age_in_years)]~set$age_in_years[order(set$age_in_years)])
}

dev.off()


par(mfrow=c(3,3))
#generative additive modelspar(mfrow=c(2,2))
for (i in 1:length(phenotype)){
  feature <- phenotype[i]
  gamModel
  gamModel_phenotype[i] <- gamlss_quad(analysisDf,feature) #have it return three different models 
  #term.plot(gamModel, pages = 1, ask = F)
  #wp(gamModel, ylim.all = .5) 
  }





### generative additive models *****
gamVen <-gamlss(VentricleVolume~pb(age_in_years)+sex+SurfaceHoles, family = GA, data=na.omit(analysisDf), trace = FALSE)
wp(gamVen, ylim.all = .5) 
drop1(gamVen) # sex failed 
term.plot(gamVen, pages = 1, ask = F)

gamCeb <-gamlss(CerebralWhiteMatterVol~pb(age_in_years)+sex+SurfaceHoles, family = GA, data=na.omit(analysisDf), trace = FALSE)
wp(gamCeb, ylim.all = .5) # show how far ordered residuals are from expected value

gamCor <-gamlss(CortexVol~pb(age_in_years)+sex+SurfaceHoles, family = GA, data=na.omit(analysisDf), trace = FALSE)
wp(gamCor, ylim.all = .5) # model is inadequate ?

gamSub <-gamlss(SubCortGrayVol~pb(age_in_years)+sex+SurfaceHoles, family = GA, data=na.omit(analysisDf), trace = FALSE)
wp(gamSub, ylim.all = .5) # this as well












