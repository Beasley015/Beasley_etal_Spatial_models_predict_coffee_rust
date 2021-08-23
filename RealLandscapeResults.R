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
library(effectsize)

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

# Histogram of raw results
ggplot()+
  geom_histogram(aes(x = land1final$PercInf, fill = "Landscape 1"),
                 binwidth = 0.05, color = "black")+
  geom_histogram(aes(x = land2final$PercInf, fill = "Landscape 2"),
                 alpha = 0.5, binwidth = 0.05, color = "black")+
  geom_vline(xintercept = mean(land1final$PercInf), linetype = "dashed")+
  geom_vline(xintercept = mean(land2final$PercInf), linetype = "dotted")+
  scale_fill_manual(values = c("darkgray", "limegreen"))+
  labs(x = "Percent Infected", y = "Frequency")+
  theme_bw(base_size = 12)+
  theme(panel.grid = element_blank(), legend.title = element_blank())

# ggsave(filename = "realdatahists.jpeg", dpi = 600)

# Rate of spread --------------------
# Calculate number of cells per landscape
cells1 <- cellStats(land1, "sum")
cells2 <- cellStats(land2, "sum")

# Get % of total repped by 1 cell
perc.cell1 <- 1/cells1
perc.cell2 <- 1/cells2

# Create columns for # of infected cells
land1res$inf.cells <- land1res$PercInf/perc.cell1
land2res$inf.cells <- land2res$PercInf/perc.cell2

# Combine data frames
all.res <- rbind(land1res, land2res)
all.res$landscape <- c(rep('Land1', nrow(land1res)),
                       rep('Land2', nrow(land2res)))

rates <- all.res %>%
  group_by(landscape, Rep) %>%
  mutate(NewInf = inf.cells-lag(inf.cells)) %>%
  summarise(mean.rate = mean(NewInf, na.rm = T))

# Histogram of avg. new infections:
ggplot(data = rates, aes(x = mean.rate, fill = landscape))+
  geom_histogram(color = "black", bins = 20)+
  scale_fill_manual(values = c("darkgray", "limegreen"))+
  labs(x = "Avg. New Infections per Time Step", y = "Frequency")+
  theme_bw(base_size = 12)+
  theme(panel.grid = element_blank(), legend.title = element_blank())

# Boxplots of avg. new infections:
ggplot(data = rates, aes(x = landscape, y = mean.rate))+
  geom_boxplot(fill = "lightgray")+
  labs(x = "Landscape", y = "Avg. Infections per Time Step")+
  theme_bw(base_size = 12)+
  theme(panel.grid = element_blank())

ggsave(filename = "RateNewInf.jpeg", dpi = 600)
  
# Look at locations of outbreaks: small landscape ----------------------  
# Get high and low quantiles
land1highs <- rates %>%
  filter(landscape == "Land1") %>%
  filter(mean.rate >= quantile(mean.rate, 0.75)) %>%
  mutate(Quant = "High")

land1lows <- rates %>%
  filter(landscape == "Land1") %>%
  filter(mean.rate <= quantile(mean.rate, 0.25)) %>%
  mutate(Quant = "Low")

# Combine
land1quants <- rbind(land1highs, land1lows)

# Get starting coords
land1smol <- subset(land1final, Rep %in% land1quants$Rep)[,c(3:5)]

# Combine data frames
land1coords <- full_join(land1quants, land1smol, by = "Rep")

# Look at locations of outbreaks: large landscape ----------------------  
# Get high and low quantiles
land2highs <- rates %>%
  filter(landscape == "Land2") %>%
  filter(mean.rate >= quantile(mean.rate, 0.75)) %>%
  mutate(Quant = "High")

land2lows <- rates %>%
  filter(landscape == "Land2") %>%
  filter(mean.rate <= quantile(mean.rate, 0.25)) %>%
  mutate(Quant = "Low")

# Combine
land2quants <- rbind(land2highs, land2lows)

# Get starting coords
land2smol <- subset(land2final, Rep %in% land2quants$Rep)[,c(3:5)]

# Combine data frames
land2coords <- full_join(land2quants, land2smol, by = "Rep")
