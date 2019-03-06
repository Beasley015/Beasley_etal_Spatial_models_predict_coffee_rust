#####################################################
# Modeling coffee rust movement through a landscape #
# N. Aristizabal, E. Beasley, & E. Bueno            #
# Spring 2019                                       #
#####################################################

#load packages, set seed
library(raster)
library(dismo)

#Make the landscape -------------------------------------------------------------
#Create blank raster
raster.size <- 1000

land <- raster(nrows = raster.size, ncols = raster.size, xmn = 0, ymn = 0, 
               xmx = raster.size, ymx = raster.size)
plot(land) #Make sure it's empty

#Fill it with random values
land[] <- runif(raster.size*raster.size, 0, 1)
plot(land) #randomly distributed values, as expected

#"Smooth" the raster using focal()
land <- focal(x = land, w = matrix(1, nrow = 25, ncol = 25), fun = mean) #need a different function; I want to retain the full range of cell values
plot(land)


