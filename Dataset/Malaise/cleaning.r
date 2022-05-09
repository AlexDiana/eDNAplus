library(here)
library(tidyverse)
library(magrittr)
library(here)
library(glue)
library(tictoc)
library(inspectdf)
library(dplyr)

setwd("C:/Users/Alex/Desktop/Malaise")
# OTUtable1_orig <- read.table(file = "OTUtable_12S_1.txt", header = T)
# OTUtable2_orig <- read.table(file = "OTUtable_12S_2.txt", header = T)
# OTUtable3_orig <- read.table(file = "OTUtable_12S_3.txt", header = T)
OTUtable1_orig <- read.table(file = "OTUtable_1_vsearch97.tsv", header = T)
OTUtable2_orig <- read.table(file = "OTUtable_2_vsearch97.tsv", header = T)
OTUtable3_orig <- read.table(file = "OTUtable_3_vsearch97.tsv", header = T)
# OTUtable1_orig <- read.table(file = "OTUtable_1_vsearch97_all.tsv", header = T)
# OTUtable2_orig <- read.table(file = "OTUtable_2_vsearch97.tsv", header = T)
# OTUtable3_orig <- read.table(file = "OTUtable_3_vsearch97.tsv", header = T)

# -------

OTUtable1 <- OTUtable1_orig
rownames_OTUtable1 <- OTUtable1$OTU
OTUtable1 <- OTUtable1[,c(2:ncol(OTUtable1))]
OTUtable1 <- t(OTUtable1)
OTUtable1 <- as.data.frame(OTUtable1)
colnames(OTUtable1) <- rownames_OTUtable1
OTUtable1$PCR <- 1

OTUtable2 <- OTUtable2_orig
rownames_OTUtable2 <- OTUtable2$OTU
OTUtable2 <- OTUtable2[,c(2:ncol(OTUtable2))]
OTUtable2 <- t(OTUtable2)
OTUtable2 <- as.data.frame(OTUtable2)
colnames(OTUtable2) <- rownames_OTUtable2
OTUtable2$PCR <- 2

OTUtable3 <- OTUtable3_orig
rownames_OTUtable3 <- OTUtable3$OTU
OTUtable3 <- OTUtable3[,c(2:ncol(OTUtable3))]
OTUtable3 <- t(OTUtable3)
OTUtable3 <- as.data.frame(OTUtable3)
colnames(OTUtable3) <- rownames_OTUtable3
OTUtable3$PCR <- 3

samplesid <- rownames(OTUtable1)

samplesid <- gsub('\\.', '-', samplesid)

OTUtable1$SampleID <- samplesid
OTUtable2$SampleID <- samplesid
OTUtable3$SampleID <- samplesid

rownames(OTUtable1) <- NULL
rownames(OTUtable2) <- NULL
rownames(OTUtable3) <- NULL

OTUtable <- rbind(OTUtable1, OTUtable2, OTUtable3)

OTUspecies <- sapply(rownames_OTUtable1, function(x){
  secondPart <- strsplit(x, split = "-")[[1]][2]
  strsplit(secondPart, split = "_")[[1]][1]
})

# link with sampling id ------------

envcov_orig <- read.table(file = "env_covariates_20220406.tsv", header = T)

envcov_orig$Site <- sapply(envcov_orig$sitename, function(x){
  if(!grepl("-",x)){
    return("-")
  } else {
    totalBits <- strsplit(x, split = "-")[[1]]
    return(paste(totalBits[-c(length(totalBits) - 0:1)], collapse = "-"))  
  }
})

envcov_orig$Sample <- envcov_orig$labname

OTUtable$Site <- envcov_orig$Site[match(OTUtable$SampleID, envcov_orig$labname)]
OTUtable$Sample <- envcov_orig$Sample[match(OTUtable$SampleID, envcov_orig$labname)]

OTUtable <- OTUtable[!is.na(OTUtable$Site) & OTUtable$Site != "#N/A",]

OTUtable <- OTUtable[,c(ncol(OTUtable) - c(1,0,3), 
                        1:(ncol(OTUtable) - 4))]

OTUtable <- OTUtable %>% arrange(desc(Site), Sample, PCR)

# add site and sample covariate ---------------

setwd("C:/Users/Alex/Desktop/Malaise")
siteinfovars1 <- readRDS("site-info-vars11-S1.rds")
siteinfovars2 <- readRDS("site-info-vars11-S2.rds")

colnames(siteinfovars2)[10] <- colnames(siteinfovars1)[10]
siteinfovars <- rbind(siteinfovars1, siteinfovars2)

siteinfovars <- siteinfovars[!duplicated(siteinfovars$SiteName),]

# uniqueSite_Xz <- setdiff(unique(OTUtable$Site), "-")

# X_z_0 <- as.matrix(envcov_orig[match(uniqueSite_Xz, envcov_orig$Site),-c(1,2,3,165,166)])

X_z <- siteinfovars[,c("SiteName", "be30", "lg_DistRoad")]

uniqueSample_Xw <- unique(OTUtable$Sample)[-c(117:124)]

X_w_0 <- envcov_orig[match(uniqueSample_Xw, envcov_orig$Sample),-c(1,2,3,165,166)]

# assign spike-ins --------------

spikein_OTU_orig <- read.table(file = "spikein_OTUs.tsv", header = T)

OTUspikes <- sapply(OTUspecies, function(x){
  if(x %in% c("14","22","6093","256609")){
    T
  } else {
    F
  }
})

# cut data set ---------

OTUtable0 <- OTUtable[4:ncol(OTUtable)]

# order by species presence
sumreadspecies <- apply(OTUtable0, 2, sum)
nonNAreadspecies <- apply(OTUtable0, 2, function(x){
  sum(x > 3)
})
# OTUtable <- OTUtable[,c(1:3, 3 + order(-sumreadspecies))]

data <- OTUtable[,1:3]
OTU <- OTUtable[,-c(1:3)]
OTU_0 <- OTU

OTU <- OTU[,order(-nonNAreadspecies)]
OTUspecies <- OTUspecies[order(-nonNAreadspecies)]
OTUspikes <- OTUspikes[order(-nonNAreadspecies)]

OTU_nonspike <- OTU[,1:100]
OTU_spike <- OTU[,OTUspikes]
OTU_spike <- OTU_spike[,c(1,2)]

# rebalance with lysis ratio ---------

OTUtable0 <- OTUtable[4:ncol(OTUtable)]

all(data$Sample %in% envcov_orig$Sample)

lysis_ratio_factor <- sapply(1:nrow(OTU_spike), function(i){
  current <- envcov_orig$lysis_ratio[which(envcov_orig$Sample == data$Sample[i])]
  if(!is.na(current)){
    return(current)
  } else {
    return(NA)    
  }
})

lysis_ratio_factor <- lysis_ratio_factor / mean(lysis_ratio_factor, na.rm = T)
lysis_ratio_factor[is.na(lysis_ratio_factor)] <- 1

i <- 0
OTU_spike <- t(apply(OTU_spike, 1, function(x){
  i <<- i + 1
  round(x * lysis_ratio_factor[i])
}))

i <- 0
OTU_nonspike <- t(apply(OTU_nonspike, 1, function(x){
  i <<- i + 1
  round(x * lysis_ratio_factor[i])
}))

# save(data, OTU_nonspike, OTU_spike,
#      uniqueSite_Xz, X_w_0, uniqueSample_Xw, X_z_0, file = "malaise.rda")

# REARRANGING FOR MODEL new -----------

K <- sapply(unique(data$Sample), function(x){
  length(unique(data$PCR[data$Sample == x]))
})

sites <- unique(data$Site)
sites <- setdiff(sites,"-")
n <- length(sites)
samples <- unique(data$Sample[data$Site %in% sites])
M_site <- sapply(1:n, function(i){
  length(unique(data$Sample[data$Site == sites[i]]))
})

emptyTubes <- length(unique(data$Sample[data$Site == "-"]))

S <- ncol(OTU_nonspike)
S_star <- ncol(OTU_spike)

# convert to y format
y <- array(NA, dim = c(sum(M_site) + emptyTubes, max(K), S + S_star),
           dimnames = list(unique(data$Sample), 1:max(K), c(colnames(OTU_nonspike),
                                                            colnames(OTU_spike))))
data_infos <- data.frame(Site = rep(0, sum(M_site) + emptyTubes),
                         Sample = rep(0, sum(M_site) + emptyTubes))
for (i in 1:n) {
  print(i)
  site <- sites[i]
  sampleSite <- unique(data$Sample[data$Site == site])
  for (m in 1:M_site[i]) {
    sample <- sampleSite[m]
    replicates <- unique(data$PCR[data$Site == site &
                                    data$Sample == sample])
    numRep <- length(replicates)
    for (k in 1:numRep) {
      rep <- replicates[k]
      for (j in 1:S) {
        y[m + sum(M_site[seq_len(i-1)]),k,j] <- OTU_nonspike[data$Site == site &
                                                      data$Sample == sample &
                                                      data$PCR == rep, j]
      }
      for (j in 1:S_star) {
        y[m + sum(M_site[seq_len(i-1)]),k,S + j] <- OTU_spike[data$Site == site &
                                                      data$Sample == sample &
                                                      data$PCR == rep, j]
      }
    }
    data_infos$Site[m + sum(M_site[seq_len(i-1)])] <- site
    data_infos$Sample[m + sum(M_site[seq_len(i-1)])] <- sample
  }
}

emptyTubes_sample <- unique(data$Sample[data$Site == "-"])
for (m in 1:emptyTubes) {
  print(m)
  site <- "-"
  sample <- emptyTubes_sample[m]
  replicates <- unique(data$PCR[data$Site == site &
                                  data$Sample == sample])
  numRep <- length(replicates)
  for (k in 1:numRep) {
    rep <- replicates[k]
    for (j in 1:S) {
      y[m + sum(M_site),k,j] <- OTU_nonspike[data$Site == site &
                                                             data$Sample == sample &
                                                             data$PCR == rep, j]
    }
    for (j in 1:S_star) {
      y[m + sum(M_site),k,S + j] <- OTU_spike[data$Site == site &
                                                              data$Sample == sample &
                                                              data$PCR == rep, j]
    }
  }
  data_infos$Site[m + sum(M_site)] <- site
  data_infos$Sample[m + sum(M_site)] <- sample
}

PCR_table <- cbind(y[,1,1:S],y[,2,1:S],y[,3,1:S])
# OTU_names <- colnames(PCR_table)[1:S]
# OTU_names <- sapply(1:length(OTU_names), function(i){
#   strsplit(OTU_names[i], "_")[[1]][1]
# })
# colnames(PCR_table)[1:S] <- OTU_names
# 
PCR_spike <- cbind(y[,1,S + seq_len(S_star)],
                   y[,2,S + seq_len(S_star)],
                   y[,3,S + seq_len(S_star)])
# OTU_names <- colnames(PCR_spike)[1:S_star]
# OTU_names <- sapply(1:length(OTU_names), function(i){
#   strsplit(OTU_names[i], "_")[[1]][1]
# })
# colnames(PCR_spike)[1:S_star] <- OTU_names

infos <- data_infos
infos$Replicates <- K
infos$Site[infos$Site == "-"] <- "empty"

# X_z <- X_z_0[,c(28, 33)]
# X_z <- apply(X_z, 2, function(x){
#   as.numeric(x)
# })
X_z[,-1] <- scale(X_z[,-1])
X_z <- X_z[match(setdiff(unique(infos$Site),"empty"), X_z[,1]),]

# X_w <- X_w_0[,c(29),drop = F]
# X_w <- apply(X_w, 2, function(x){
#   as.numeric(x)
# })
# X_w <- scale(X_w)
# X_w <- X_w[match(setdiff(data_infos$Sample,
#                          data$Sample[data$Site == "-"]), uniqueSample_Xw),,drop=F]

v_spikes <- matrix(0, nrow = sum(M_site) + emptyTubes, S_star)

spikedSample <- rbind(matrix(1, nrow = sum(M_site), S_star),
                      matrix(1, nrow = 5, S_star),
                      matrix(0, nrow = 3, S_star))

data <- list("PCR_table" = PCR_table,
             "infos" = infos,
             "X_z" = X_z,
             "PCR_spike" = PCR_spike,
             "spikedSample " = spikedSample)

setwd("~/eDNAPlus/Dataset")
save(data, file = "data_malaise.rda")

# REARRANGING FOR MODEL 1 -----------

# data$SiteSample <- paste0(data$Site, data$Sample)
K <- sapply(unique(data$Sample), function(x){
  length(unique(data$PCR[data$Sample == x]))
})

sites <- unique(data$Site)
sites <- setdiff(sites,"-")
n <- length(sites)
samples <- unique(data$Sample[data$Site %in% sites])
M_site <- sapply(1:n, function(i){
  length(unique(data$Sample[data$Site == sites[i]]))
})

emptyTubes <- length(unique(data$Sample[data$Site == "-"]))

S <- ncol(OTU_nonspike)
S_star <- ncol(OTU_spike)

# convert to y format
y <- array(NA, dim = c(sum(M_site) + emptyTubes, max(K), S + S_star),
           dimnames = list(unique(data$Sample), 1:max(K), c(colnames(OTU_nonspike),
                                                            colnames(OTU_spike))))
data_infos <- data.frame(Site = rep(0, sum(M_site) + emptyTubes),
                         Sample = rep(0, sum(M_site) + emptyTubes))
for (i in 1:n) {
  print(i)
  site <- sites[i]
  sampleSite <- unique(data$Sample[data$Site == site])
  for (m in 1:M_site[i]) {
    sample <- sampleSite[m]
    replicates <- unique(data$PCR[data$Site == site &
                                    data$Sample == sample])
    numRep <- length(replicates)
    for (k in 1:numRep) {
      rep <- replicates[k]
      for (j in 1:S) {
        y[m + sum(M_site[seq_len(i-1)]),k,j] <- OTU_nonspike[data$Site == site &
                                                      data$Sample == sample &
                                                      data$PCR == rep, j]
      }
      for (j in 1:S_star) {
        y[m + sum(M_site[seq_len(i-1)]),k,S + j] <- OTU_spike[data$Site == site &
                                                      data$Sample == sample &
                                                      data$PCR == rep, j]
      }
    }
    data_infos$Site[m + sum(M_site[seq_len(i-1)])] <- site
    data_infos$Sample[m + sum(M_site[seq_len(i-1)])] <- sample
  }
}

emptyTubes_sample <- unique(data$Sample[data$Site == "-"])
for (m in 1:emptyTubes) {
  print(m)
  site <- "-"
  sample <- emptyTubes_sample[m]
  replicates <- unique(data$PCR[data$Site == site &
                                  data$Sample == sample])
  numRep <- length(replicates)
  for (k in 1:numRep) {
    rep <- replicates[k]
    for (j in 1:S) {
      y[m + sum(M_site),k,j] <- OTU_nonspike[data$Site == site &
                                                             data$Sample == sample &
                                                             data$PCR == rep, j]
    }
    for (j in 1:S_star) {
      y[m + sum(M_site),k,S + j] <- OTU_spike[data$Site == site &
                                                              data$Sample == sample &
                                                              data$PCR == rep, j]
    }
  }
  data_infos$Site[m + sum(M_site)] <- site
  data_infos$Sample[m + sum(M_site)] <- sample
}

all.equal(unique(data_infos$Site), sites)
all.equal(unique(data_infos$Sample), sites)

X_z <- X_z_0[,c(28, 33)]
X_z <- apply(X_z, 2, function(x){
  as.numeric(x)
})
X_z <- scale(X_z)
X_z <- X_z[match(setdiff(unique(data_infos$Site),"-"), uniqueSite_Xz),]

X_w <- X_w_0[,c(29),drop = F]
X_w <- apply(X_w, 2, function(x){
  as.numeric(x)
})
X_w <- scale(X_w)
X_w <- X_w[match(setdiff(data_infos$Sample,
                         data$Sample[data$Site == "-"]), uniqueSample_Xw),,drop=F]

v_spikes <- matrix(0, nrow = sum(M_site) + emptyTubes, S_star)

spikedSample <- rbind(matrix(1, nrow = sum(M_site), S_star),
                      matrix(1, nrow = 5, S_star),
                      matrix(0, nrow = 3, S_star))

save(y, M_site, K, emptyTubes, S_star, X_w, X_z, 
     data_infos, v_spikes, spikedSample, 
     file = here("Data","Malaise.Rdata"))

