#
#This script given a dataframe will 
# 1. Pass through Jenna's filter 
# 2. gamlss models(linear -> cubic -> gams -> centiles) on 4 phenotype -- 
# (VentricleVolume, CerebralWhiteMatterVol, CortexVol, SubCortGrayVol)
#

#----
#the gamlss framework packages 
library(gamlss) #to fit model
library(gamlss.cens) #fit censored response variables
library(gamlss.dist) #additional distributions 
library(gamlss.mx) # for fitting finite mixture distributions - uses nnet
library(gamlss.tr) # for fitting truncated distributions
#other packages
library(ggplot2)
library(tidyverse)
library(ggpubr)

#load data ----
analysisDf <- read.csv(file.choose())
summary(analysisDf)

#Jenna's filtering ----
# Drop any rows where neurofibromatosis is in the scan_reason_categories column
# analysisDf <- masterDf[!grepl("neurofibromatosis", masterDf$scan_reason_primary), ]
analysisDf <- analysisDf[!grepl("neurofibromatosis", analysisDf$scan_reason_categories), ]

# Make an age in years column from the age in days column
analysisDf$age_in_years <- analysisDf$age_at_scan_days/365.25

# Some of the columns need to be factors - the checks is for troubleshotting gamlss is difficult if the wrong thing is not a factor
toFactor <- c('sex', 'fs_version', 'MagneticFieldStrength', 'scan_id', 'scan_reason_primary')
analysisDf[toFactor] <- lapply(analysisDf[toFactor], factor)

is.factor(analysisDf$sex)
is.factor(analysisDf$fs_version)
is.factor(analysisDf$MagneticFieldStrength)
is.factor(analysisDf$scan_id)
is.factor(analysisDf$scan_reason_primary)

# Drop any 1.5T scans
analysisDf <- analysisDf[analysisDf$MagneticFieldStrength != "1.5",]

# Drop any scans with ratings less than 0 (-1 rated scans were post contrast in Jenna's initial manual qc)
analysisDf <- analysisDf[analysisDf$rawdata_image_grade >= 0, ]

# Sort the dataframe by patient_id and scanner_id --- # We only one scan per subject
analysisDf <- analysisDf[ with(analysisDf, order(analysisDf$patient_id, analysisDf$scan_id)), ]
# Drop all but the first occurrence of each patient_id
analysisDf <- analysisDf[!duplicated(analysisDf$patient_id), ]
# Convert patient_id to a factor - idk why I did this?
analysisDf$patient_id <- droplevels(as.factor(analysisDf$patient_id))

# Add a column for TCV (Total Cerebrum Volume)
analysisDf$TCV <- analysisDf$TotalGrayVol + analysisDf$CerebralWhiteMatterVol
# Add a column: TotalBrainVolume = TotalGrayVol + CerebralWhiteMatterVol + VentricleVolume + SubCortGrayVol
analysisDf$TotalBrainVol <- analysisDf$TotalGrayVol + analysisDf$CerebralWhiteMatterVol + analysisDf$VentricleVolume + analysisDf$SubCortGrayVol

#prelimary plots 
hist(analysisDf$age_in_years)

ggplot(data = analysisDf, aes(x = age_in_years,y = SurfaceHoles, col = fs_version))+
  geom_point()+
  labs(title = "Euler's vs Age")+
  theme(legend.position = "none")

gven <- ggplot(data = analysisDf, aes(x = age_in_years,y = VentricleVolume, col = fs_version))+
  geom_point()+
  labs(title = "Ventricle Volume vs Age")+  
  theme(legend.position = "none")

gwmv <- ggplot(data = analysisDf, aes(x = age_in_years,y = CerebralWhiteMatterVol, col = fs_version))+
  geom_point()+
  labs(title = "WMV vs Age")+
  theme(legend.position = "none")

ggmv <- ggplot(data = analysisDf, aes(x = age_in_years,y = CortexVol, col = fs_version))+
  geom_point()+
  labs(title = "GMV vs Age")+
  theme(legend.position = "none")

gsgmv <- ggplot(data = analysisDf, aes(x = age_in_years,y = SubCortGrayVol, col = fs_version))+
  geom_point()+
  labs(title = "sGMV vs Age")+
  theme(legend.position = "none")

ggarrange(gven, gwmv, ggmv, gsgmv + rremove("x.text"), 
          labels = c("A", "B", "C","D"),
          ncol = 2, nrow = 2)


####gamlss ------
#1. NORMAL LINEAR _ No log
par(mfrow=c(2,2))

lrVen <- gamlss(VentricleVolume~age_in_years+sex, data = na.omit(analysisDf), family = "NO")
plot(VentricleVolume~(age_in_years),data = analysisDf,col = "lightblue")
lines(fitted(lrVen)~(analysisDf$age_in_years))

lrWMV <- gamlss(CerebralWhiteMatterVol~age_in_years, data = na.omit(analysisDf), family = "GG")
plot(CerebralWhiteMatterVol~(age_in_years),data = analysisDf,col = "lightblue")
lines(fitted(lrWMV)~(analysisDf$age_in_years))

lrGMV <- gamlss(CortexVol~age_in_years, data = na.omit(analysisDf), family = "GG")
plot(CortexVol~(age_in_years),data = analysisDf,col = "lightblue")
lines(fitted(lrGMV)~(analysisDf$age_in_years))

lrGMV <- gamlss(SubCortGrayVol~age_in_years, data = na.omit(analysisDf), family = "GG")
plot(SubCortGrayVol~(age_in_years),data = analysisDf,col = "lightblue")
lines(fitted(lrGMV)~(analysisDf$age_in_years))

dev.off()

###FITTING CUBIC POLY
#
#Classic way to set model: 
     #aModel <- gamlss(VentricleVolume~age_in_years+I(age_in_years^2)+I(age_in_years^3), data = df, family = NO)
#Easier to set polynomial before: 
df <- transform(analysisDf, age_in_years2 = age_in_years^2, age_in_years3 = age_in_years^3)

par(mfrow=c(2,2))

cubVen <-gamlss(VentricleVolume~age_in_years+age_in_years2+age_in_years3, data=na.omit(df), family = NO)
plot(VentricleVolume~age_in_years, col = "lightblue", data = df)
lines(fitted(cubVen)[order(df$age_in_years)]~df$age_in_years[order(df$age_in_years)])

cubCeb <-gamlss(CerebralWhiteMatterVol~age_in_years+age_in_years2+age_in_years3, data=na.omit(df), family = NO)
plot(CerebralWhiteMatterVol~age_in_years, col = "lightblue", data = df)
lines(fitted(cubCeb)[order(df$age_in_years)]~df$age_in_years[order(df$age_in_years)])

cubCor <-gamlss(CortexVol~age_in_years+age_in_years2+age_in_years3, data=na.omit(df), family = NO)
plot(CortexVol~age_in_years, col = "lightblue", data = df)
lines(fitted(cubCor)[order(df$age_in_years)]~df$age_in_years[order(df$age_in_years)])

cubSub <-gamlss(SubCortGrayVol~age_in_years+age_in_years2+age_in_years3, data=na.omit(df), family = NO)
plot(SubCortGrayVol~age_in_years, col = "lightblue", data = df)
lines(fitted(cubSub)[order(df$age_in_years)]~df$age_in_years[order(df$age_in_years)])

dev.off()

### Generative additive models *****
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

#The Error --> in while (abs(olddv - dv) > cc && itn < cyc) { : 
#missing value where TRUE/FALSE needed 
# is a distribution error, helps to use dev.off() and run with GA if it shows up

par(mfrow=c(2,2))

gamVen <-gamlss(VentricleVolume~pb(age_in_years)+sex+SurfaceHoles, family = GG, data=na.omit(analysisDf), trace = F)

gamCeb <-gamlss(CerebralWhiteMatterVol~pb(age_in_years)+sex+SurfaceHoles, family = GG, data=na.omit(analysisDf), trace = F)

gamCor <-gamlss(CortexVol~pb(age_in_years)+sex+SurfaceHoles, family = GG, data=na.omit(analysisDf), trace = F)

gamSub <-gamlss(SubCortGrayVol~pb(age_in_years)+sex+SurfaceHoles, family = GG, data=na.omit(analysisDf), trace = F)

dev.off()

###
#4 Centile Curves using BCPE
#

par(mfrow=c(2,2))

Ven1 <- gamlss(VentricleVolume~pb(age_in_years), sigma.formula=~pb(age_in_years),nu.formula=~pb(age_in_years), 
                 data = analysisDf, family = GG)
Ven2 <- gamlss(VentricleVolume~pb(age_in_years), sigma.formula=~pb(age_in_years), nu.formula=~pb(age_in_years),
                 tau.formula=~pb(age_in_years), data = analysisDf, start.from=Ven1, family ="GG")

centiles.fan(Ven2, xvar=analysisDf$age_in_years, cent=c(3,10,25,50,75,90,97), colors = "topo",
             ylab="Ventricle Volume", xlab = "age(years)")


Ceb1 <- gamlss(CerebralWhiteMatterVol~pb(age_in_years), sigma.formula=~pb(age_in_years),nu.formula=~pb(age_in_years), 
                 data = analysisDf, family = GG)
Ceb2 <- gamlss(CerebralWhiteMatterVol~pb(age_in_years), sigma.formula=~pb(age_in_years), nu.formula=~pb(age_in_years),
                 tau.formula=~pb(age_in_years), data = analysisDf, start.from=Ceb1, family ="GG")

centiles.fan(Ceb2, xvar=analysisDf$age_in_years, cent=c(3,10,25,50,75,90,97), colors = "topo",
             ylab="CerebralWhiteMatterVol", xlab = "age(years)")



Cor1 <- gamlss(CortexVol~pb(age_in_years), sigma.formula=~pb(age_in_years),nu.formula=~pb(age_in_years), 
                 data = analysisDf, family = GG)
Cor2 <- gamlss(CortexVol~pb(age_in_years), sigma.formula=~pb(age_in_years), nu.formula=~pb(age_in_years),
                 tau.formula=~pb(age_in_years), data = analysisDf, start.from=Cor1, family ="GG")

centiles.fan(Cor2, xvar=analysisDf$age_in_years, cent=c(3,10,25,50,75,90,97), colors = "topo",
             ylab="CortexVol", xlab = "age(years)")



Gray1 <- gamlss(SubCortGrayVol~pb(age_in_years), sigma.formula=~pb(age_in_years),nu.formula=~pb(age_in_years), 
                 data = analysisDf, family = GG)
Gray2 <- gamlss(SubCortGrayVol~pb(age_in_years), sigma.formula=~pb(age_in_years), nu.formula=~pb(age_in_years),
                 tau.formula=~pb(age_in_years), data = analysisDf, start.from=Gray1, family ="GG")

centiles.fan(Gray2, xvar=analysisDf$age_in_years, cent=c(3,10,25,50,75,90,97), colors = "topo",
             ylab="SubCortGrayVol", xlab = "age(years)")

####
