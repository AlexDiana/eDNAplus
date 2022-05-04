setwd("C:/Users/Alex/Dropbox/R Folder/PostDoc/Data/Oregon")

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
