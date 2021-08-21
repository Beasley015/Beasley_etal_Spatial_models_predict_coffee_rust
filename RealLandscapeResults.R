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

# Load model results -----------------------
# Smol landscape
land1res <- read.table("land1.csv", sep = ",")
colnames(land1res) <- c("Time", "PercInf", "X", "Y", "Rep")
land1res$Rep <- land1res$Rep+1

# Large landscape
land2res <- read.table("land2.csv", sep = ",")
colnames(land2res) <- c("Time", "PercInf", "X", "Y", "Rep")
land2res$Rep <- land2res$Rep+1

# Clean data -------------------------------
# Pull out data from final time step
land1final <- subset(land1res, land1res$Time == 999)
land2final <- subset(land2res, land2res$Time == 999)

ggplot()+
  geom_histogram(aes(x = land1final$PercInf, fill = "Landscape 1"),
                 binwidth = 0.05, color = "black")+
  geom_histogram(aes(x = land2final$PercInf, fill = "Landscape 2"),
                 alpha = 0.5, binwidth = 0.05, color = "black")+
  scale_fill_manual(values = c("darkgray", "limegreen"))+
  labs(x = "Percent Infected", y = "Frequency")+
  theme_bw(base_size = 12)+
  theme(panel.grid = element_blank(), legend.title = element_blank())

# ggsave(filename = "realdatahists.jpeg", dpi = 600)
