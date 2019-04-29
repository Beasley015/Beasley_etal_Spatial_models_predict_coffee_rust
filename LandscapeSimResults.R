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
is.atomic(as.matrix(as.data.frame(output.list)))
str_which(output.list,"def")
def<- c("10")
disp<- c("2", "3", "4", "5")
str_match_all(output.list, c(def, disp))

str_extract_all(is.atomic(as.matrix(as.data.frame(output.list))), "disp")

?str_extract_all


for( i in seq_along(output.list)){
  output.list[[i]]$def <- str_extract_all(output.list[i],c("def10", "disp"))
  
}
head()


    