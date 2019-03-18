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
land <- setValues(land, runif(raster.size*raster.size, 0, 1))
plot(land) #randomly distributed values, as expected

#"Smooth" the raster using focal()
land <- focal(x = land, w = matrix(1, nrow = 19, ncol = 15), fun = mean)
plot(land)

#Create a binary raster to represent coffee ----------------------------------
#Blank raster
coffee <- raster(nrows = raster.size, ncols = raster.size, xmn = 0, ymn = 0, 
                 xmx = raster.size, ymx = raster.size)

#Use lowest friction values from land raster to represent coffee
coffee <- setValues(coffee, ifelse(getValues(land) <= 0.47, 1, 0))

#plot the new binary raster
plot(coffee) 

#Merge coffee and landscape rasters ------------------------------------
#Change lanscape so coffee areas = NA
land[which(getValues(coffee)==1)] <- NA

#Merge rasters
landscape <- merge(land, coffee)

#Plot the new raster
plot(landscape)
