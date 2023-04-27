############################################## 
############ Packages Preparation ############
##############################################


#Following are the packages needed for fitting the runEDNA model 
install.packages("devtools")
library(devtools)
install_github("alexdiana1992/eDNAPlus",force = TRUE)
library(eDNAPlus)


############################################## 
############ Data Loading ####################
##############################################


#find your R working directory, this is the file path on your computer  that sets the default location of any files 
#you read into R, or save out of R
getwd()
#if the getwd() is not the same as the location that you saved the files SimData, TrueCor, datapoly and the BioMass file then first find 
#the location of those files and store them in a known to you location and set that location as your new 
#working directory by clicking on Session/Set Working Directory/Choose Directory and then browsing as usual to locate your preferred working directory
#then load the SimData
data = readRDS("SimData.RData")
#load the species correlation parameters that used to simulate the SimData set 
true_cor = readRDS("TrueCor.RData")
#load the X,Y axis values used later for displaying the biommases of the species over the sites
datapoly = readRDS("datapoly.RData")
#load the true spatial Biomass parameters used to simulate the SimData set
BioMass = readRDS("BioMass.RData")
#check that your data are on the correct form
head(data$infos,3)
head(data$PCR_table,3)
head(data$X_z,3)
head(data$PCR_spike,3)


############################################## 
############ Model Fitting ###########
##############################################

#We put a prior on the PCR bias variance 
priors <- list("sigma_u" = 3)
#we define the Markov Chain Monte Carlo parameters
MCMCparams <- list("nchain" = 1,
                   "nburn" = 1000,
                   "niter" = 2500,
                   "nthin" = 1)
#finally we run the model
modelResults <- runEDNA(data,
                        priors,
                        jointSpecies = T,
                        MCMCparams = MCMCparams)



############################################## 
############ Outputs ####################
##############################################

install.packages("Matrix") #installs the package Matrix from CRAN
library(Matrix)            #reads the package Matrix
install.packages("MASS")
library(MASS)
install.packages("ggcorrplot")
library(ggcorrplot)
install.packages("gridExtra")
library(gridExtra)

#correlation matrix between spices
p1 = plotCorrelationMatrix(modelResults, idxSpecies = 1:10)
p2 = ggcorrplot(true_cor)
grid.arrange(p1, p2, ncol=2)
#covariate coefficients for biomass availability
plotBiomassCoefficients(modelResults,idxSpecies = 1:10)
#Distribution of Biomass for each species independently, below plots for 3 species
q1 = plotBiomassesMap(modelResults,datapoly,1)
q1
q2 = plotBiomassesMap(modelResults,datapoly,2)
q2
q3 = plotBiomassesMap(modelResults,datapoly,3)
q3
#the probability of a positive pcr when the sample includes species biomass 
plotTruePositiveProbability(modelResults, idxSpecies = 1:10)
#the probability of a positive pcr when the sample does not include species biomass
plotFalsePositiveProbability(modelResults, idxSpecies = 1:10)
#we compare the true species-specific biomass across surveyed sites with the estimated,
# for the first two species
logz_species_1 = BioMass[,1] #first species
logz_species_2 = BioMass[,2] #second species
datapoly2 <- datapoly
#plot for first species
datapoly2$LogBiomass <- rep(logz_species_1, each = 4)
f1 = ggplot() + geom_polygon(data = datapoly2, aes(x = x, 
                                                   y = y, group = id, alpha = LogBiomass), color = "black") + 
  scale_alpha_continuous(range = c(0, 1)) + xlab("X") + 
  ylab("Y") + scale_x_continuous(breaks = c()) + scale_y_continuous(breaks = c()) + 
  theme(plot.title = element_text(hjust = 0.5, size = 20), 
        axis.title = element_text(size = 20, face = "bold"), 
        axis.text = element_text(size = 13, face = "bold", 
                                 angle = 90), panel.background = element_rect(fill = "white", 
                                                                              color = "black"))
#comparison of biomasses for first species on the same plot, left is the estimated parameters and right the true
grid.arrange(q1, f1, ncol=2)
#plot for second species
datapoly2$LogBiomass <- rep(logz_species_2, each = 4)
f2 = ggplot() + geom_polygon(data = datapoly2, aes(x = x, 
                                                   y = y, group = id, alpha = LogBiomass), color = "black") + 
  scale_alpha_continuous(range = c(0, 1)) + xlab("X") + 
  ylab("Y") + scale_x_continuous(breaks = c()) + scale_y_continuous(breaks = c()) + 
  theme(plot.title = element_text(hjust = 0.5, size = 20), 
        axis.title = element_text(size = 20, face = "bold"), 
        axis.text = element_text(size = 13, face = "bold", 
                                 angle = 90), panel.background = element_rect(fill = "white", 
                                                                              color = "black"))
#comparison of biomasses for second species on the same plot, left is the estimated parameters and right the true
grid.arrange(q2, f2, ncol=2)
