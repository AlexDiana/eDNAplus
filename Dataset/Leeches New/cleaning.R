library(tidyverse)
library(magrittr)
library(here)
library(glue)
library(tictoc)
library(inspectdf)

setwd("C:/Users/Alex/Dropbox/R Folder/PostDoc/Data/Leeches New")
# OTUtable1_orig <- read.table(file = "OTUtable_12S_1.txt", header = T)
# OTUtable2_orig <- read.table(file = "OTUtable_12S_2.txt", header = T)
# OTUtable3_orig <- read.table(file = "OTUtable_12S_3.txt", header = T)
OTUtable1_orig <- read.table(file = "OTUtable_12S_1_reduced.txt", header = T)
OTUtable2_orig <- read.table(file = "OTUtable_12S_2_reduced.txt", header = T)
OTUtable3_orig <- read.table(file = "OTUtable_12S_3_reduced.txt", header = T)

# doug filtering -----------

filter1 <- function(otutab){
  otutab %>% 
    filter(order != "root,unk") %>% 
    filter(family != "Primates,unk", family != "Primates,Hominidae") %>% 
    select(-starts_with("HUMAN"))
}

OTUtable1 <- filter1(OTUtable1_orig)
OTUtable2 <- filter1(OTUtable2_orig)
OTUtable3 <- filter1(OTUtable3_orig)

filter2 <- function(otutab){
  otutab %>% 
    filter(is.na(species) | str_detect(genus, 'unk') | str_detect(species, 'unk'))
}

OTUtable1_unk <- filter2(OTUtable1)
OTUtable2_unk <- filter2(OTUtable2)
OTUtable3_unk <- filter2(OTUtable3)

filter3 <- function(otutab){
  otutab %>% 
    filter(!str_detect(genus, 'unk') & !str_detect(species, 'unk') & !is.na(species))
}

OTUtable1_no_unk <- filter3(OTUtable1)
OTUtable2_no_unk <- filter3(OTUtable2)
OTUtable3_no_unk <- filter3(OTUtable3)


summarise1 <- function(otutab){
  otutab %>% 
    group_by(species) %>% 
    summarise(
      Cluster_Number = first(Cluster_Number),
      QueryID = first(QueryID),
      Size = sum(Size), 
      class = first(class),
      order = first(order),
      family = first(family),
      genus = first(genus),
      species = first(species),
      # across(CX01:NC2, sum)
      across(CX03:NHNC, sum)
    ) %>% 
    relocate(species, .after = genus)
}

OTUtable1_no_unk_summ <- summarise1(OTUtable1_no_unk)
OTUtable2_no_unk_summ <- summarise1(OTUtable2_no_unk)
OTUtable3_no_unk_summ <- summarise1(OTUtable3_no_unk)

OTUtable1_reduced <- bind_rows(OTUtable1_no_unk_summ, OTUtable1_unk)
OTUtable2_reduced <- bind_rows(OTUtable2_no_unk_summ, OTUtable2_unk)
OTUtable3_reduced <- bind_rows(OTUtable2_no_unk_summ, OTUtable2_unk)

# -------


OTUtable1 <- OTUtable1_reduced
rownames_OTUtable1 <- OTUtable1$species
OTUtable1 <- OTUtable1[,c(9:770)]
OTUtable1 <- t(OTUtable1)
OTUtable1 <- as.data.frame(OTUtable1)
colnames(OTUtable1) <- rownames_OTUtable1
OTUtable1$PCR <- 1

OTUtable2 <- OTUtable2_reduced
rownames_OTUtable2 <- OTUtable2$species
OTUtable2 <- OTUtable2[,c(9:770)]
OTUtable2 <- t(OTUtable2)
OTUtable2 <- as.data.frame(OTUtable2)
colnames(OTUtable2) <- rownames_OTUtable2
OTUtable2$PCR <- 2

OTUtable3 <- OTUtable3_reduced
rownames_OTUtable3 <- OTUtable3$species
OTUtable3 <- OTUtable3[,c(9:770)]
OTUtable3 <- t(OTUtable3)
OTUtable3 <- as.data.frame(OTUtable3)
colnames(OTUtable3) <- rownames_OTUtable3
OTUtable3$PCR <- 3

labids <- rownames(OTUtable1)
labids_nop <- sapply(seq_along(labids), function(i) strsplit(labids[i], "[.]")[[1]][1])

OTUtable1$Sample <- labids
OTUtable2$Sample <- labids
OTUtable3$Sample <- labids

rownames(OTUtable1) <- NULL
rownames(OTUtable2) <- NULL
rownames(OTUtable3) <- NULL

OTUtable <- rbind(OTUtable1, OTUtable2, OTUtable3)

# colnames(OTUtable)[is.na(colnames(OTUtable))] <- "NA"

# link with areas id ------------

labtopoly <- read.csv(file = "file1.csv")

sites_samples <- rep(NA, length(labids_nop))

polygonid <- labtopoly$Polygon_ID[match(labids_nop, labtopoly$Lab_ID)]
polygonid <- sapply(seq_along(polygonid), function(i){
  strsplit(polygonid[i], split = "or")[[1]][1]
})

rangerid <- labtopoly$Ranger_ID[match(labids_nop, labtopoly$Lab_ID)]

sites_samples <- polygonid

idxNAPolygonid <- which(is.na(sites_samples))

for (idx in idxNAPolygonid) {
  if(!is.na(rangerid[idx])){
    rangerID <- rangerid[idx]
    polygonIDofranger <- labtopoly$Polygon_ID[labtopoly$Ranger_ID == rangerID]
    # print(paste0("ID = ",rangerID))
    # print(polygonIDofranger)
    if(all(is.na(polygonIDofranger))){
      sites_samples[idx] <- rangerID
    } else {
      sites_samples[idx] <- polygonIDofranger[!is.na(polygonIDofranger)][1]
    }
  } 
}

# sampleid <- as.numeric(sites_samples)
# 
# sampleid[is.na(sampleid)] <- labtopoly$Ranger_ID[match(labids_nop[is.na(sampleid)], labtopoly$Lab_ID)]
# 
# unique_rg <- unique(labtopoly$Ranger_ID[match(labids_nop[is.na(sampleid)], labtopoly$Lab_ID)])
# unique_rg %in% labtopoly$Ranger_ID[match(labids_nop[!is.na(sampleid)], labtopoly$Lab_ID)]

# OTUtable$Site <- sites_samples
OTUtable$Site <- rep(sites_samples, times = 3)

# exclude NA and #N/A
OTUtable <- OTUtable[!is.na(OTUtable$Site) & OTUtable$Site != "#N/A",]

OTU_0 <- OTUtable[,ncol(OTUtable) - 0:2]
OTU_0$Indexes <- 1:nrow(OTU_0)

OTU_0 <- OTU_0 %>% arrange(Site, Sample, PCR)

# reorder
OTUtable <- OTUtable[OTU_0$Indexes,]

OTUtable <- OTUtable[,c(ncol(OTUtable) - 0:2, 1:(ncol(OTUtable) - 3))]
# OTUtable <- OTUtable[!is.na(OTUtable$Site),]

# add site covariate ---------------

X_site <- read.csv(file = "covariates.csv")

uniqueSite <- unique(OTUtable$Site)

X_z <- as.matrix(X_site[match(uniqueSite, X_site$polygon_ID),-1])

# X_z <- as.matrix(X_site[match(sampleid, X_site$polygon_ID),-1])

# add sample covariate ----------

uniqueSample <- unique(OTUtable$Sample)

leechesQty <- read.csv(file = "leechesNum.csv")

X_w <- as.matrix(leechesQty$leech_qty[match(uniqueSample, leechesQty$Lab_ID)])

# cut data set ---------

# order by species presence
sumreadspecies <- apply(OTUtable[,-c(1:3)], 2, sum)
OTUtable <- OTUtable[,c(1:3, 3 + order(-sumreadspecies))]

data <- OTUtable[,1:3]
OTU <- OTUtable[,-c(1:3)]
OTU_0 <- OTU

# cut some species
nonNaPCR <- apply(OTU, 2, function(x){
  sum(x > 100)
})
sumreadspecies <- apply(OTU, 2, sum)

# condition

sum(nonNaPCR > 4)

OTU <- OTU[,nonNaPCR > 4]

# reorder

nonNaPCR <- apply(OTU, 2, function(x){
  sum(x > 100)
})

OTU <- OTU[,order(-nonNaPCR)]

save(data, OTU, X_w, X_z, file = "leeches_new.rda")

# REARRANGING FOR MODEL -----------

data$SiteSample <- paste0(data$Site, data$Sample)
K <- sapply(unique(data$SiteSample), function(x){
  length(unique(data$PCR[data$SiteSample == x]))
})

sites <- unique(data$Site)
n <- length(sites)
samples <- unique(data$Sample)
M_site <- sapply(1:n, function(i){
  length(unique(data$Sample[data$Site == sites[i]]))
})
S <- ncol(OTU)

# convert to y format
y <- array(NA, dim = c(sum(M_site), max(K), S),
           dimnames = list(unique(data$SiteSample), 1:max(K), colnames(OTU)))
data_infos <- data.frame(Site = rep(0, sum(M_site)),
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
    data_infos$Site[m + sum(M_site[seq_len(i-1)])] <- site
    data_infos$Sample[m + sum(M_site[seq_len(i-1)])] <- sample
  }
}

X_z <- X_z[,c(5,11)]

ncov_z <- ncol(X_z)

X_z <- scale(X_z)
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

save(OTU, y, X_w, X_z, n, M_site, K, X_w, data_infos, data, emptySites, file = here("Data","Ponds.Rdata"))

