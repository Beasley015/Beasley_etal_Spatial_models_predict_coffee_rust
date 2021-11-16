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
library(patchwork)
library(ggnewscale)

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
land1raw <- trim(land1raw, values = NA)

land2extent <- extent(bothlands, 1, nrow(rawmat), colNA[3]+1,ncol(rawmat))
land2raw <- crop(bothlands, land2extent)
land2raw <- trim(land2raw, values = NA)

# Extract coffee values
land1 <- land1raw == 5
land1[land1 == 0] <- NA

land2 <- land2raw == 5
land2[land2 == 0] <- NA

# Get % coffee per landscape -------------------
land1vals <- freq(land1raw)
land1vals[5,2]/sum(land1vals[1:7,2])

land2vals <- freq(land2raw)
land2vals[5,2]/sum(land2vals[1:7,2])

# Get % deforested per landscape -----------------
sum(land1vals[c(2,4,6),2])/sum(land1vals[-c(5,7),2])
sum(land2vals[c(2,4,6),2])/sum(land2vals[-c(5,7),2])

# Calculate aggregation index -------------------
lsm_l_ai(land1)
lsm_l_ai(land2)

# Load model results -----------------------
# Read in raw numpy arrays
land1mat <- read.table("land1raw.csv", sep = ",")
colnames(land1mat) <- as.character(1:ncol(land1mat))

land2mat <- read.table("land2raw.csv", sep = ",")
colnames(land2mat) <- as.character(1:ncol(land2mat))

# Smol landscape
land1res <- read.table("land1.csv", sep = ",")
colnames(land1res) <- c("Time", "PercInf", "Y", "X", "Rep")
land1res <- land1res %>%
  mutate(Rep = Rep+1, X = X+1, Y = Y+1)

# Large landscape
land2res <- read.table("land2.csv", sep = ",")
colnames(land2res) <- c("Time", "PercInf", "X", "Y", "Rep")
land2res <- land2res %>%
  mutate(Rep = Rep+1, X = X+1, Y=Y+1)

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
  scale_x_discrete(labels = c("Landscape1", "Landscape2"))+
  theme_bw(base_size = 10)+
  theme(panel.grid = element_blank())

# ggsave(filename = "RateNewInf.tiff", dpi = 1200, width = 84, height = 70,
#        units = "mm")
  
# Look at locations of outbreaks: small landscape ----------------------  
# Get high and low quantiles
land1coord <- rates %>%
  filter(landscape == "Land1") %>%
  mutate(Quant = case_when(mean.rate >= quantile(mean.rate, 0.75)~"High",
                           mean.rate <= quantile(mean.rate, 0.25)~"Low",
                           TRUE ~ "Mid"))

# Get starting coords
land1smol <- subset(land1final, Rep %in% land1coord$Rep)[,c(3:5)]

# Combine data frames
land1coords <- full_join(land1coord, land1smol, by = "Rep")

# Look at locations of outbreaks: large landscape ----------------------  
# Assign to quantiles
land2coord <- rates %>%
  filter(landscape == "Land2") %>%
  mutate(Quant = case_when(mean.rate >= quantile(mean.rate, 0.75)~"High",
                           mean.rate <= quantile(mean.rate, 0.25)~"Low",
                           TRUE ~ "Mid"))

# Get starting coords
land2smol <- subset(land2final, Rep %in% land2coord$Rep)[,c(3:5)]

# Combine data frames
land2coords <- full_join(land2coord, land2smol, by = "Rep")

# Plot on raster -------------------------
# Quick glance
qplot(data = land1coords, x = X, y = Y, color = Quant)
qplot(data = land2coords, x = X, y = Y, color = Quant)

# Convert to data frames
land1df <- land1mat %>%
  mutate(row = rownames(land1mat)) %>%
  pivot_longer(!row, names_to = "col", values_to = "cover_type") %>%
  filter(cover_type != "NaN") %>%
  mutate(cover_type = as.character(cover_type), row = as.numeric(row),
         col = as.numeric(col)) %>%
  mutate(cover_type = case_when(cover_type == "5" ~ "Coffee",
                                cover_type == "1" ~ "Forest",
                                TRUE ~ "Other"))
  

land2df <- land2mat %>%
  mutate(row = rownames(land2mat)) %>%
  pivot_longer(!row, names_to = "col", values_to = "cover_type") %>%
  filter(cover_type != "NaN") %>%
  mutate(cover_type = factor(cover_type), row = as.numeric(row),
         col = as.numeric(col)) %>%
  mutate(cover_type = case_when(cover_type == "5" ~ "Coffee",
                                cover_type == "1" ~ "Forest",
                                TRUE ~ "Other"))

# Plots
land1start <- ggplot()+
  geom_raster(data = land1df, mapping = aes(x = col, y = row,
                                            fill = cover_type))+
  scale_fill_manual(values = gray.colors(n = 8, start = 0.5),
                    name = "Cover Type")+
  geom_point(data = land1coords, mapping = aes(x = X, y = Y,
                                               color = Quant),
             size = 1.5)+
  scale_color_viridis_d(begin = 0.5, breaks = c("Low", "Mid", "High"),
                        name = "Rate of Spread")+
  scale_y_reverse()+
  theme_bw(base_size = 10)+
  theme(axis.title = element_blank(), axis.text = element_blank(), 
        panel.grid = element_blank())

# ggsave("land1start.jpeg", dpi = 600)

land2start <- ggplot()+
  geom_raster(data = land2df, mapping = aes(x = row, y = col, 
                                            fill = cover_type))+
  scale_fill_manual(values = gray.colors(n = 8, start = 0.5),
                    name = "Cover Type")+
  geom_point(data = land2coords, mapping = aes(x = X, y = Y,
                                               color = Quant),
             size = 1.5)+
  coord_flip()+
  scale_x_reverse()+
  scale_color_viridis_d(begin = 0.5, breaks = c("Low", "Mid", "High"),
                        name = "Rate of Spread")+
  theme_bw(base_size = 10)+
  theme(axis.title = element_blank(), axis.text = element_blank(), 
        panel.grid = element_blank(), legend.position = 'none')

# ggsave("land2start.jpeg", dpi = 600)

(land1start|land2start)+
  guides(color = guide_legend(title = "Rate of Spread"))+
  plot_layout(guides = "collect")+
  plot_annotation(tag_levels = "a")

# ggsave("landstarts.tiff", dpi = 1200, width = 174, height = 70,
#        units = "mm")

# Draw "buffers" around points & get % cover --------------------------

# Function to get % of each land cover around starting locations
get_surroundings <- function(start_loc, lands){
  vals <- list()
  
  for(i in 1:nrow(start_loc)){
    # Filter values from df
    smol <- lands %>%
      filter(col > start_loc$X[i]-50 & col < start_loc$X[i]+50) %>%
      filter(row > start_loc$Y[i]-50 & row < start_loc$Y[i]+50)
    
    # Get counts of values and convert to % of total
    vals[[i]] <- smol %>%
      group_by(cover_type) %>%
      summarise(count = n()) %>%
      mutate(percent_cover = count/sum(count)) %>%
      select(cover_type, percent_cover) %>%
      pivot_wider(names_from = cover_type, values_from = percent_cover)
    
    if(nrow(vals[[i]]) == 0){
      start_loc[i,] <- rep(NA, 6)
    }
  }
  
  # Convert list to data frame
  vals.df <- bind_rows(vals)
  
  start_loc <- na.omit(start_loc)
  
  # Append to coords object
  lands <- cbind(start_loc, vals.df)
  
  return(lands)
}

around.land1 <- get_surroundings(start_loc = land1coords, lands = land1df)
around.land2 <- get_surroundings(start_loc = land2coords, lands = land2df)

# Quick plots to look at how % coffee compares across categories
qplot(data = around.land1, x = Quant, y = Coffee, geom = "boxplot")
qplot(data = around.land2, x = Quant, y = Coffee, geom = "boxplot")
qplot(data = around.land1, x = Coffee, y = mean.rate)
qplot(data = around.land2, x = Coffee, y = mean.rate)

# Make nice plots
local.coffee1 <- ggplot(data = around.land1, aes(x = Coffee, 
                                                 y = mean.rate))+
  geom_point()+
  labs(x = "% Coffee Cover", y = "Mean Rate of Spread")+
  theme_bw(base_size = 12)+
  theme(panel.grid = element_blank())

cor.test(around.land1$Coffee, around.land1$mean.rate)

local.coffee2 <- ggplot(data = around.land2, aes(x = Coffee,
                                                 y = mean.rate))+
  geom_point()+
  labs(x = "% Coffee Cover", y = "Mean Rate of Spread")+
  theme_bw(base_size = 12)+
  theme(panel.grid = element_blank())

cor.test(around.land2$Coffee, around.land2$mean.rate)

local.coffee <- (local.coffee1|local.coffee2) +
  plot_annotation(tag_levels = "a")

# ggsave(file = "local_coffee.tiff", dpi = 1200, width = 174, units = "mm")

local.other1 <- ggplot(data = around.land1, aes(x = Other, 
                                                 y = mean.rate))+
  geom_point()+
  labs(x = "% Other Land Cover", y = "Mean Rate of Spread")+
  theme_bw(base_size = 12)+
  theme(panel.grid = element_blank())

local.other2 <- ggplot(data = around.land2, aes(x = Other,
                                                 y = mean.rate))+
  geom_point()+
  labs(x = "% Other Land Cover", y = "Mean Rate of Spread")+
  theme_bw(base_size = 12)+
  theme(panel.grid = element_blank())

local.other <- (local.other1|local.other2) +
  plot_annotation(tag_levels = "a")

cor.test(around.land1$Other, around.land1$mean.rate)
cor.test(around.land2$Other, around.land2$mean.rate)

# ggsave(local.other, filename = "LocalOther.tiff", dpi = 1200, width = 174,
#        units = "mm")

cor.test(around.land1$Forest, around.land1$mean.rate)
cor.test(around.land2$Forest, around.land2$mean.rate)

# Look at domain effects in landscape 1 -----------------
# Get distance from edge
get.distance <- function(x){
  min.dist <- logical()
  
  for(i in 1:nrow(x)){
    test.loc <- x[i,]

    boundaries <- which(land1mat == "NaN", arr.ind = T)

    dist <- abs((test.loc$X-boundaries[,2]))+
                         abs((test.loc$Y-boundaries[,1]))
    
    min.dist[i] <- min(dist)
  }
  
  return(min.dist)
}

land1dist <- get.distance(land1coords)
land2dist <- get.distance(land2coords)

# Compare to rate of spread
land1coords$dist <- land1dist
land2coords$dist <- land2dist

qplot(data = land1coords, x = dist, y = mean.rate)
summary(lm(data = land1coords, mean.rate~dist))
cor.test(land1coords$dist, land1coords$mean.rate, method = "spearman")

qplot(data = land2coords, x = dist, y = mean.rate)
summary(lm(data = land2coords, mean.rate~dist))
cor.test(land2coords$dist, land2coords$mean.rate, method = "spearman")

# distance doesn't seem to matter...

# try Moran's i -------------------------
# Landscape 1
land1distmat <- as.matrix(dist(cbind(land1coords$X, land1coords$Y)))

land1inv <- 1/land1distmat

diag(land1inv) <- 0

Moran.I(land1coords$mean.rate, land1inv)

# Landscape 2
land2distmat <- as.matrix(dist(cbind(land2coords$X, land2coords$Y)))

land2inv <- 1/land2distmat

diag(land2inv) <- 0

Moran.I(land2coords$mean.rate, land2inv)

#both landscapes autocorrelated