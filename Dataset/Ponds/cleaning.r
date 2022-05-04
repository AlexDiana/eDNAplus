
library(dplyr)

OTUtable1_0 <- read.table("OTUtable1.txt", header = T)
OTUtable2_0 <- read.table("OTUtable2.txt", header = T)
OTUtable3_0 <- read.table("OTUtable3.txt", header = T)

OTUtable1 <- OTUtable1_0
rownames(OTUtable1) <- OTUtable1[,1]
OTUtable1 <- OTUtable1[,-1]
OTUtable1 <- t(OTUtable1)
OTUtable1 <- as.data.frame(OTUtable1)
OTUtable1$PCR <- 1

OTUtable2 <- OTUtable2_0
rownames(OTUtable2) <- OTUtable2[,1]
OTUtable2 <- OTUtable2[,-1]
OTUtable2 <- t(OTUtable2)
OTUtable2 <- as.data.frame(OTUtable2)
OTUtable2$PCR <- 2

OTUtable3 <- OTUtable3_0
rownames(OTUtable3) <- OTUtable3[,1]
OTUtable3 <- OTUtable3[,-1]
OTUtable3 <- t(OTUtable3)
OTUtable3 <- as.data.frame(OTUtable3)
OTUtable3$PCR <- 3

labids <- rownames(OTUtable1)

OTUtable1$Sample <- labids
OTUtable2$Sample <- labids
OTUtable3$Sample <- labids

rownames(OTUtable1) <- NULL
rownames(OTUtable2) <- NULL
rownames(OTUtable3) <- NULL

OTUtable <- rbind(OTUtable1, OTUtable2, OTUtable3)
OTUtable$Site <- OTUtable$Sample

OTUtable <- OTUtable[,c(ncol(OTUtable) - 0:2, 1:(ncol(OTUtable) - 3))]

sumreadspecies <- apply(OTUtable[,-c(1:3)], 2, sum)
OTUtable <- OTUtable[,c(1:3, 3 + order(-sumreadspecies))]

OTUtable <- OTUtable %>% arrange(Site, Sample, PCR)

# X_site <- X_covar[match(OTUtable$Site, X_covar$polygon_ID),-1]

data <- OTUtable[,1:3]
OTU <- OTUtable[,-c(1:3)]

# save(OTUtable, file = "OTUtable.rda")

# covariate -------

covariate <- read.csv(file = "fhtnm_20200229.csv")

colnames(covariate)[1] <- "Sample"

save(data, OTU, covariate, file = "pondsdata.rda")

# AFTER CLEANING ----------

library(here)
load(here("Data/Ponds","pondsdata.rda"))

# cut some species
OTU_sub <- OTU[,1:200000]
nonNaPCR1 <- apply(OTU_sub, 2, function(x){
  sum(x > 10)
})
sumreadspecies1 <- apply(OTU_sub, 2, sum)

OTU_sub2 <- OTU[,200001:ncol(OTU)]
nonNaPCR2 <- apply(OTU_sub2, 2, function(x){
  sum(x > 10)
})
sumreadspecies2 <- apply(OTU_sub2, 2, sum)

nonNaPCR <- c(nonNaPCR1, nonNaPCR2)
sumreadspecies <- c(sumreadspecies1, sumreadspecies2)

OTU <- OTU[,nonNaPCR > 100]

sites <- unique(data$Site)
covariate_subset <- covariate[match(sites, covariate$Sample),]
X_z <- matrix(as.numeric(as.matrix(covariate_subset[,c(39, 40),drop = F])), ncol = 2)

# X_w <- numLeeches

data$Empty <- 0

# data$Site[is.na(data$Site)] <- "Unknown"

data$SiteSample <- paste0(data$Site, data$Sample)
K <- sapply(unique(data$SiteSample), function(x){
  length(unique(data$PCR[data$SiteSample == x]))
})

sites <- unique(data$Site)
emptySites <- rep(0, length(sites))
# #sapply(sites, function(x){
# # data$Empty[data$Site == x][1]
# # })
#
n <- length(sites)
samples <- unique(data$Sample)
M_site <- sapply(1:n, function(i){
  length(unique(data$Sample[data$Site == sites[i]]))
})
S <- ncol(OTU)


# convert to y format
y <- array(NA, dim = c(sum(M_site), max(K), S),
           dimnames = list(unique(data$SiteSample), 1:max(K), colnames(OTU)))
data_short <- data.frame(Site = rep(0, sum(M_site)),
                         Sample = rep(0, sum(M_site)))
for (i in 1:n) {
  print(i)
  site <- sites[i]
  sampleSite <- unique(data$Sample[data$Site == site])
  for (m in 1:M_site[i]) {
    sample <- sampleSite[m]
    replicates <- unique(data$PCR[data$Site == site &
                                    data$Sample == sample])
    numRep <- K[m + sum(M_site[seq_len(i-1)])]
    for (k in 1:numRep) {
      rep <- replicates[k]
      for (j in 1:S) {
        y[m + sum(M_site[seq_len(i-1)]),k,j] <- OTU[data$Site == site &
                                                      data$Sample == sample &
                                                      data$PCR == rep, j]
      }
    }
    data_short$Site[m + sum(M_site[seq_len(i-1)])] <- site
    data_short$Sample[m + sum(M_site[seq_len(i-1)])] <- sample
  }
}
#
# X_z <- as.matrix(X_site[match(sites, X_site$polygon_ID),-1])

ncov_z <- ncol(X_z)

X_z <- scale(X_z)
# X_z <- cbind(1,scale(X_z))
# # assign NA to mean (remove later)
X_z[is.na(X_z[,2]),2] <- 0
X_z[is.na(X_z[,1]),1] <- 0

X_w <- matrix(NA, nrow = sum(M_site), ncol = 0)
# X_w <- as.matrix(X_w[match(samples, X_w$Lab_ID),])
# r <- as.numeric(X_w)
## assign NA to mean (remove later)
# r <- scale(r)
# r[is.na(r)] <- 0
#

emptySites <- rep(0, n)

save(OTU, y, X_w, X_z, n, M_site, K, X_w, data_short, data, emptySites, file = here("Data","Ponds.Rdata"))
