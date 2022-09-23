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

#factor Synthe v 
#toFactor <- c( 'fs_version', 'scan_id')
#analysisDf[toFactor] <- lapply(analysisDf[toFactor], factor)


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





 
#The features of interist are:  TCV, GMV, WMV, sGMV, and Ventricles ---- 

#graph centiles first 
gtcv <- ggplot(data = analysisDf, aes(x = age_in_years,y = TCV, col = fs_version))+
  geom_point()+
  labs(title = "TCV vs Age")+  
  theme(legend.position = "none")

gven <- ggplot(data = analysisDf, aes(x = age_in_years,y = VentricleVolume, col = fs_version))+
  geom_point()+
  labs(title = "Ventricles vs Age")+  
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

ggarrange(gtcv, gven, gwmv, ggmv, gsgmv + rremove("x.text"), 
          labels = c("A", "B", "C","D"),
          ncol = 2, nrow = 2)







###Centile 
par(mfrow=c(2,3))


TCV1 <- gamlss(TCV~pb(age_in_years), sigma.formula=~pb(age_in_years),nu.formula=~pb(age_in_years), 
               data = analysisDf, family = GG)
TCV2 <- gamlss(TCV~pb(age_in_years), sigma.formula=~pb(age_in_years), nu.formula=~pb(age_in_years),
               tau.formula=~pb(age_in_years), data = analysisDf, start.from=TCV1, family ="GG")

centiles(TCV2, xvar=analysisDf$age_in_years, cent=c(3,10,25,50,75,90,97), colors = "topo",
             ylab="TCV", xlab = "age(years)")


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



Cor1 <- gamlss(TotalGrayVol~pb(age_in_years), sigma.formula=~pb(age_in_years),nu.formula=~pb(age_in_years), 
               data = analysisDf, family = GG)
Cor2 <- gamlss(TotalGrayVol~pb(age_in_years), sigma.formula=~pb(age_in_years), nu.formula=~pb(age_in_years),
               tau.formula=~pb(age_in_years), data = analysisDf, start.from=Cor1, family ="GG")

centiles.fan(Cor2, xvar=analysisDf$age_in_years, cent=c(3,10,25,50,75,90,97), colors = "topo",
             ylab="TotalGrayVol", xlab = "age(years)")



Gray1 <- gamlss(SubCortGrayVol~pb(age_in_years), sigma.formula=~pb(age_in_years),nu.formula=~pb(age_in_years), 
                data = analysisDf, family = GG)
Gray2 <- gamlss(SubCortGrayVol~pb(age_in_years), sigma.formula=~pb(age_in_years), nu.formula=~pb(age_in_years),
                tau.formula=~pb(age_in_years), data = analysisDf, start.from=Gray1, family ="GG")

centiles.fan(Gray2, xvar=analysisDf$age_in_years, cent=c(3,10,25,50,75,90,97), colors = "topo",
             ylab="SubCortGrayVol", xlab = "age(years)")













library(gamlss)
data()
data("lip")
plot(lip$X1.d, lip$X2.d)
head(lip)

Gray1 <- gamlss(X1.d~pb(X2.d), sigma.formula=~pb(X2.d),nu.formula=~pb(X2.d), 
                data = lip, family = NO)
Gray2 <- gamlss(SubCortGrayVol~pb(age_in_years), sigma.formula=~pb(age_in_years), nu.formula=~pb(age_in_years),
                tau.formula=~pb(age_in_years), data = analysisDf, start.from=Gray1, family ="GG")

centiles.fan(Gray2, xvar=analysisDf$age_in_years, cent=c(3,10,25,50,75,90,97), colors = "topo",
             ylab="SubCortGrayVol", xlab = "age(years)")
















