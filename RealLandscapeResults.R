##########################################################
# Modeling coffee rust movement through a real landscape #
# N. Aristizabal, E. Beasley, E. Bueno, & E. White       #
# Summer 2021- ????                                      #
##########################################################

# Load packages and data ------------------------------
# Packages
library(landscapemetrics)
library(rgdal)
library(raster)
library(tidyverse)

# Raw raster data
x <- new("GDALReadOnlyDataset", "./coffeeRust_landuses")
getDriver(x)
getDriverLongName(getDriver(x))
xx<-asSGDF_GROD(x)
bothlands <- raster(xx)

# Format rasters ------------------------------
# Convert to matrix
rawmat <- is.na(as.matrix(bothlands))

# Get cols that are all NA
colNA <- which(colSums(rawmat) == nrow(bothlands))

# Trim both rasters
land1extent <- extent(bothlands, 1, nrow(rawmat), 1, colNA[1]-1)
land1raw <- crop(bothlands, land1extent)

land2extent <- extent(bothlands, 1, nrow(rawmat), colNA[3]+1,ncol(rawmat))
land2raw <- crop(bothlands, land2extent)

# Extract coffee values
land1 <- land1raw == 5
land1[land1 == 0] <- NA

land2 <- land2raw == 5
land2[land2 == 0] <- NA

# Calculate aggregation index -------------------
lsm_l_ai(land1)
lsm_l_ai(land2)
