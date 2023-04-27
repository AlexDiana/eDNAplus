# install.packages("devtools")
# library(devtools)
# install_github("alexdiana1992/eDNAPlus")
library(eDNAPlus); library(Matrix)
library(here)

# READ DATA ------

setwd(here("Dataset/Simulated"))

data <- read.csv(file = "PCR.csv")
data_spike <- read.csv(file = "PCR_spike.csv")

infos <- data[,1:3]
PCR_table <- data[,4:33]

X_z <- data[,c(1, 34),drop = F]
# X_w <- data[,c(2, 14, 15)]

PCR_spike <- data_spike[,1:6]

data <- list("infos" = infos,
             "PCR_table" = PCR_table,
             "X_z" = X_z,
             # "X_w" = X_w,
             "PCR_spike" = PCR_spike)

# MALAISE DATA --------

load("~/eDNAPlus/Dataset/data_malaise.rda")
# load("~/eDNAPlus/Dataset/data_lakes.rda")
# data <- data_lakes

# ----

priors <- list("sigma_u" = 1,
               "l_gp" = .05)

# mcmc params
{
  MCMCparams <- list("nchain" = 1,
                     "nburn" = 10000,
                     "niter" = 10000,
                     "nthin" = 1,
                     "iterToAdapt" = 100000)
}

modelResults <- runEDNA(data,
                        priors,
                        jointSpecies = T,
                        spatialCorr = T,
                        MCMCparams = MCMCparams)


plotBiomassCoefficients(modelResults)
plotCorrelationMatrix(modelResults, 
                      numSigCorrs = 15)

View(FPNP_plot(modelResults, idxSpecies = 2))

load("~/eDNAPlus/data/polygons_info.rda")
plotBiomassesMap(modelResults, datapoly, 4)

#


