library(dplyr); library(here); library(reshape2)

data <- read.csv(file = here("Dataset/Francesco Ficetola","lakedata.csv"))

# plateinfo <- read.csv(file = here("Data/Francesco Ficetola","sampleplate.csv"), header = F)
# 
# names(plateinfo)[1] <- "plate"
# names(plateinfo)[4] <- "sample"

# plateinfo$sample <- substr(plateinfo$sample, 1, 5)

# isolate columns i need
data_pcr <- 
  # data[,grepl("Family_name",colnames(data)) | 
  data[,grepl("Scientific_name",colnames(data)) |
                   grepl("sample.ANT",colnames(data)) | 
                   grepl("sample.RMT",colnames(data)) | 
                   grepl("sample.TEM",colnames(data)) ]

data_pcr <- data_pcr %>% group_by(Scientific_name) %>% summarise_all(sum)
rownames_data_pcr <- data_pcr$Scientific_name
data_pcr <- t(data_pcr)
data_pcr <- data_pcr[-1,]
colnames(data_pcr) <- rownames_data_pcr
mode(data_pcr) <- "numeric"
data_pcr <- as.data.frame(data_pcr)
#extract site, sample and pcr from the name
rownames(data_pcr) <- substr(rownames(data_pcr), 8, nchar(rownames(data_pcr)))
data_pcr$Sample <- substr(rownames(data_pcr), 1, 5)
data_pcr$EnvSample <- substr(rownames(data_pcr), 6, 6)
data_pcr$Replicate <- substr(rownames(data_pcr), 7, 8)
# reorder columns
data_pcr <- data_pcr[,c(ncol(data_pcr) - 2:0, 1:(ncol(data_pcr) - 3) )]

# data_pcr$Plate <- plateinfo$plate[match(rownames(data_pcr), plateinfo$sample)]

# match year to data --------

{
  sampleage <- read.csv(file =  here("Dataset/Francesco Ficetola","sampleage.csv"))
  
  # sampleage <- sampleage[sampleage$Age_BC.AD != "-",]
  
  colnames(sampleage)[2] <- "Age"
  
  # isolate unique elements of the sample (only first 5 letters matter)
  sampleage$Sample <- substr(sampleage$Sample, 1, 5)
  sampleage <- sampleage[!duplicated(sampleage$Sample),]
  sampleage$Age <- as.numeric(sampleage$Age)
  
  # uniqueYears <- unique(sampleage$Age_BC.AD)
  # uniqueYears <- sort(as.numeric(uniqueYears))
  # 
  # library(ggplot2)
  # ggplot(data = NULL, aes(x = uniqueYears, xend = uniqueYears,
  #                         y = 0, yend = 1)) +
  #   geom_segment(size = 1, color = "black") +  xlab("") + scale_y_continuous(breaks = c(), name = "") + theme_bw() +
  #   scale_x_continuous(breaks =
  #                        c(min(uniqueYears),seq(-7000,1500,by = 1000), max(uniqueYears)))
}

data_pcr$Age <- as.numeric(sampleage$Age[match(data_pcr$Sample, sampleage$Sample)])
data_pcr <- data_pcr[order(data_pcr$Age),]

OTU <- data_pcr[,-c(1,2,3,ncol(data_pcr) - 0:1)]
data <- data_pcr[,c(1,2,3)]

idxEmptyTubes <- grep("TEM",data$Sample)

data_TR <- data[-idxEmptyTubes,]
data_ET <- data[idxEmptyTubes,]

OTU_TR <- OTU[-idxEmptyTubes,]
OTU_ET <- OTU[idxEmptyTubes,]

# covariates ----------------

variables <- read.table(here("Dataset/Francesco Ficetola","independent_variables.txt"))

data_age <- merge(data_TR, sampleage, by = "Sample", all.x = T)

X_sites <- 
  data_pcr[-idxEmptyTubes,c(1,2,3,ncol(data_pcr) - 0)] %>%
  # data_TR %>% 
  group_by(Sample) %>% 
  summarise(Sample = Sample[1],
            Age = Age[1])

X_sites <- as.data.frame(X_sites)

X_sites$Age <- sapply(1:nrow(X_sites), function(i){
  if(is.na(X_sites$Age[i])){
    return(NA)
  } else {
    # return(variables$Age[which.min(abs(variables$Age - X_sites$Age[i]))]  )
    return(variables$Age[which.min(abs(variables$Age - X_sites$Age[i]))])
  }
})

X <- merge(X_sites, variables, by = "Age", all.x = T)

X <- X %>% 
  rename(Site = Sample)

# X <- matrix(1, 
#             nrow = length(unique(data_TR$Sample)),
#             ncol = 1)
# X <- X[,-1]

data_TR <- data_TR %>% 
  rename(Site = Sample,
         Sample = EnvSample,
         PCR = Replicate)

data <- data %>% 
  rename(Site = Sample,
         Sample = EnvSample,
         PCR = Replicate)

data_ET <- data_ET %>% 
  rename(Site = Sample,
         Sample = EnvSample,
         PCR = Replicate)

# save(X, file = "covariates.rda")
# save(data, OTU, X, file = "francescodata.rda")

# REARRANGING FOR MODEL -----------

data$SiteSample <- paste0(data$Site, data$Sample)
K <- sapply(unique(data$SiteSample), function(x){
  length(unique(data$PCR[data$SiteSample == x]))
})

sites <- unique(data_TR$Site)

sampleage$Age[match(sites, sampleage$Sample)]

n <- length(sites)
samples <- unique(data_TR$Sample)
M_site <- sapply(1:n, function(i){
  length(unique(data_TR$Sample[data_TR$Site == sites[i]]))
})
emptyTubes <- length(unique(data$Sample[idxEmptyTubes]))
S <- ncol(OTU)

# convert to y format
y <- array(NA, dim = c(sum(M_site) + emptyTubes, max(K), S),
           dimnames = list(unique(data$SiteSample), 1:max(K), colnames(OTU)))
data_infos <- data.frame(Site = rep(0, sum(M_site) + emptyTubes),
                         Sample = rep(0, sum(M_site) + emptyTubes))
for (i in 1:n) {
  print(i)
  site <- sites[i]
  sampleSite <- unique(data_TR$Sample[data_TR$Site == site])
  for (m in 1:M_site[i]) {
    sample <- sampleSite[m]
    replicates <- unique(data_TR$PCR[data_TR$Site == site &
                                    data_TR$Sample == sample])
    numRep <- K[m + sum(M_site[seq_len(i-1)])]
    for (k in 1:numRep) {
      rep <- replicates[k]
      for (j in 1:S) {
        y[m + sum(M_site[seq_len(i-1)]),k,j] <- OTU_TR[data_TR$Site == site &
                                                      data_TR$Sample == sample &
                                                      data_TR$PCR == rep, j]
      }
    }
    data_infos$Site[m + sum(M_site[seq_len(i-1)])] <- site
    data_infos$Sample[m + sum(M_site[seq_len(i-1)])] <- sample
  }
}

sitesEmpty <- unique(data_ET$Site)

for (m in seq_len(emptyTubes)) {
  site <- sitesEmpty[m]
  replicates <- unique(data_ET$PCR[data_ET$Site == site])
  numRep <- K[sum(M_site) + m]
  for (k in 1:numRep) {
    rep <- replicates[k]
    for (j in 1:S) {
      y[sum(M_site) + m,k,j] <- OTU_ET[data_ET$Site == site &
                                         data_ET$PCR == rep, j]
    }
  }
  data_infos$Site[sum(M_site) + m] <- 0
  data_infos$Sample[sum(M_site) + m] <- sample
}

X_z_noname <- X[,c(4,5,6,7,8,9),drop = F]

ncov_z <- ncol(X_z_noname)

X_z <- scale(X_z_noname)

X_z <- cbind(X[,c(2)], X_z_noname)
# assign NA to mean (remove later)
X_z[24,c(4,5,6,7)] <- X_z[23,c(4,5,6,7)]
X_z[25,c(4,5,6,7)] <- X_z[26,c(4,5,6,7)]
# X_z[is.na(X_z[,1]),1] <- 0

X_t <- X[,c(2,1)]
X_t[,2] <- scale(X_t[,2])

X_w <- matrix(NA, nrow = sum(M_site), ncol = 0)
# X_w <- as.matrix(X_w[match(samples, X_w$Lab_ID),])
# r <- as.numeric(X_w)
## assign NA to mean (remove later)
# r <- scale(r)
# r[is.na(r)] <- 0
#

save(OTU, y, X_w, X_z, n, M_site, K, X_w, data_infos, emptyTubes, file = here("Data","Lakes.Rdata"))

# REARRANGE FOR MODEL 2 -------

data_infos$Sample[length(data_infos$Sample) - 0:2] <- c("A","B","C")
data_infos$Sample <- paste(data_infos$Site, data_infos$Sample, sep = "")

colnames(X_z)[1] <- "Site"

data_infos$Site[data_infos$Site == "0"] <- "empty"
data_lakes <- list()
data_lakes$PCR_table <- cbind(y[,1,], y[,2,], y[,3,], y[,4,])
data_lakes$infos <- data.frame(data_infos, Replicates = K)
data_lakes$X_z <- X_z
data_lakes$X_s <- X_t

data <- data_lakes



setwd("~/eDNAPlus/Dataset")
save(data, file = "data_lakes.rda")

# ------------------------------------------------------------------------
 
# variables$Age <- sapply(1:nrow(variables), function(i){
  # sampleage$Age[which.min(abs(sampleage$Age - variables$Age[i]))]
# })
X <- merge(variables, sampleage)

samplesNotNAAges <- X$Sample[complete.cases(X)]


# X <- X[match(samplesNotNAAges, X$Sample),]

save(X, file = "covariates_francesco.rda")

X <- X[,2:8]

Y <- length(samplesNotNAAges)
M_y <- sapply(seq_len(length(samplesNotNAAges)), function(i){
  length(unique(data_pcr$EnvSample[data_pcr$Sample==samplesNotNAAges[i]]))
})

K_ym <- matrix(NA, nrow = Y, ncol = max(M_y))
for (i in 1:Y) {
  for (j in 1:M_y[i]) {
    K_ym[i,j] <- length(unique(data_pcr$Replicate[
      data_pcr$Sample==unique(data_pcr$Sample)[i] & data_pcr$EnvSample==unique(data_pcr$EnvSample)[j]
      ]))
  }
}

S <- ncol(data_pcr) - 3

data_pcr_nonna <- data_pcr[data_pcr$Sample %in% samplesNotNAAges,]

y_jimk <- array(0, dim = c(Y, max(M_y), max(K_ym), S))
for (i in 1:nrow(data_pcr)) {
  idxSite <- match(data_pcr$Sample[i], samplesNotNAAges)
  idxSample <- match(data_pcr$EnvSample[i], unique(data_pcr$EnvSample[data_pcr$Sample == data_pcr$Sample[i]]))
  idxReplicate <- match(data_pcr$Replicate[i], unique(data_pcr$Replicate[data_pcr$Sample == data_pcr$Sample[i] & 
                                                                        data_pcr$EnvSample == data_pcr$EnvSample[i]]))
  if(!is.na(idxSite)){
    y_jimk[idxSite, idxSample, idxReplicate,] <- as.numeric(data_pcr[i,4:ncol(data_pcr)])  
  }
  
}

speciesNames <- colnames(data_pcr)[4:ncol(data_pcr)]

save(y_jimk, Y, M_y, K_ym, X, speciesNames, file = "francescodata.rda")


# plot --------

data_species <- t(data[1:10,])
data_species <- data_species[17:400,]
rownames(data_species) <- substr(rownames(data_species),8,nchar(rownames(data_species)) - 1)

mode(data_species) <- "numeric"

# data_species1 <- t(data[2,])
# data_species1 <- data_species1[17:400,]

# names(data_species1) <- substr(names(data_species1),8,nchar(names(data_species1)) - 1)

yearSample <- sampleage$Age_BC.AD[match(
  rownames(data_species),
  substr(sampleage$Sample, 1, nchar(sampleage$Sample) - 1),
         )] 

data_ts <- data.frame(Sample = rownames(data_species), 
                      Year = as.numeric(yearSample), 
                      Reads = data_species)

# 

library(dplyr)
# data_grouped <- data_ts %>% group_by(Year) %>% summarise(Reads = mean(Reads))
data_grouped <- data_ts %>% group_by(Year) %>% summarise(across(5:6, ~ mean(.x, na.rm = TRUE)))
data_grouped <- data_ts %>% group_by(Year) %>% summarise(across(5:6, ~ mean(.x, na.rm = TRUE)))

library(reshape2)

data_grouped_long <- melt(data_grouped, id.vars = "Year")

ggplot(data_grouped_long, aes(Year,value, col=variable)) + geom_line(size = 1) +
  geom_point(size = 2) + 
  theme(legend.position="none") + theme_bw() + xlab("Year") + ylab("Reads")

qplot(data_grouped$Year, data_grouped$Reads.3)
qplot(data_grouped$Year, data_grouped$Reads.4)
qplot(data_grouped$Year, data_grouped$Reads.5)


# PLATES 




      