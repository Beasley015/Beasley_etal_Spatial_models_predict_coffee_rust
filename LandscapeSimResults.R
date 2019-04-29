#####################################################
# Modeling coffee rust movement through a landscape #
# N. Aristizabal, E. Beasley, & E. Bueno            #
# Spring 2019                                       #
#####################################################

#Libraries
library(tidyverse)
library(ggplot2)

# Read in model outputs --------------------------------------
# Read all csv's into a list
filenames <- list.files("Outputs", pattern = "*.csv", full.names = T)
shortnames <- list.files("Outputs", pattern = "*.csv")
output.list <- lapply(filenames, read.csv, header = F)

# Rename columns of each dataframe
newnames <- c("Time", "PercInf")
output.list <- lapply(output.list, setNames, newnames)
head(output.list[[1]])

#Add column to denote replicates
replicate <- logical()
for(i in 1:50){
  new <- rep(i, 1000)
  replicate <- append(replicate, new)

output.list <- lapply(output.list, cbind, replicate)    
