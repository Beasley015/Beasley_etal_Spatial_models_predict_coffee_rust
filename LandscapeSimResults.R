#####################################################
# Modeling coffee rust movement through a landscape #
# N. Aristizabal, E. Beasley, & E. Bueno            #
# Spring 2019                                       #
#####################################################

# Read in model outputs --------------------------------------
# Read all csv's into a list
filenames <- list.files("Outputs", pattern = "*.csv", full.names = T)
output.list <- lapply(filenames, read.csv, header = F)

# Rename columns of each dataframe
newnames <- c("Time", "PercInf")
output.list <- lapply(output.list, setNames, newnames)
head(output.list)

#Libraries
library(tidyverse)
library(ggplot2)
#adding columns defors and disp
defors <- c("10")
disp<- c("2")
for( i in seq_along(output.list)){
  output.list[[i]]$defors <- rep(defors[i],nrow(output.list[[i]]))
  output.list[[i]]$disp<-rep(disp[i], nrow(output.list[[i]]))
  }
print(output.list[[1]])

library(stringr)

library(tools)
files<-file_path_sans_ext(list.files(pattern = "*.csv"))
as.vector(files)
files

output.list[[1]]



    