#####################################################
# Modeling coffee rust movement through a landscape #
# N. Aristizabal, E. Beasley, E. Bueno, & E. White  #
# Spring 2019- ????                                 #
#####################################################
library(tidyverse)
library(rcompanion)
library(viridis)
library(patchwork)
library(agricolae)

# Sample beta distributions (for figures) ----------------------------------
# neighbor1 <- rbeta(10000, 2, 8)
# neighbor2 <- rbeta(10000, 5, 5)
# neighbor3 <- rbeta(10000, 8, 2)
# 
# ggplot(mapping = aes(x = neighbor1))+
#   geom_density(size = 1)+
#   labs(x = "Infection probability", y = "Density")+
#   scale_x_continuous(limits = c(0,1), expand = c(0,0))+
#   scale_y_continuous(limits = c(0,4.5), expand = c(0,0))+
#   theme_bw()+
#   theme(panel.grid = element_blank(), axis.text.x = element_text(size = 16),
#         axis.text.y = element_blank(), axis.title = element_blank())
# 
# ggsave(file = "beta1.jpg")
# 
# ggplot(mapping = aes(x = neighbor2))+
#   geom_density(size = 1)+
#   labs(x = "Infection probability", y = "Density")+
#   scale_x_continuous(limits = c(0,1), expand = c(0,0))+
#   scale_y_continuous(limits = c(0,4.5), expand = c(0,0))+
#   theme_bw()+
#   theme(panel.grid = element_blank(), axis.text.x = element_text(size = 16),
#         axis.text.y = element_blank(), axis.title = element_blank())
# 
# ggsave(file = "beta2.jpg")
# 
# ggplot(mapping = aes(x = neighbor3))+
#   geom_density(size = 1)+
#   labs(x = "Infection probability", y = "Density")+
#   scale_x_continuous(limits = c(0,1), expand = c(0,0))+
#   scale_y_continuous(limits = c(0,4.5), expand = c(0,0))+
#   theme_bw()+
#   theme(panel.grid = element_blank(), axis.text.x = element_text(size = 16),
#         axis.text.y = element_blank(), axis.title = element_blank())
# 
# ggsave(file = "beta3.jpg")


# Read in model outputs --------------------------------------
# Read all csv's into a list
filenames <- list.files("Outputs", pattern = "*.csv", full.names = T)
shortnames <- list.files("Outputs", pattern = "*.csv")
output.list <- lapply(filenames, read.csv, header = F)

# Rename columns of each dataframe
newnames <- c("Time", "PercInf","X","Y", "?")
output.list <- lapply(output.list, setNames, newnames)
output.list <- map(output.list, ~ (.x %>% select(-'?')))

# Set up data -----------------------------------------------
# Add column to denote replicates
replicate <- logical()
for(i in 1:50){
  new <- rep(i, 1000)
  replicate <- append(replicate, new)
}
  
output.list <- lapply(output.list, cbind, replicate)  

# Pull deforestation and dispersion from file names
loop.ready <- c(1:length(shortnames))
def <- list()
disp <- list()
cluster <- list()

for(i in loop.ready) {
  def[[i]] <- strsplit(shortnames[[i]], split = "def|disp|cluster|.csv")[[1]][2]
  disp[[i]] <- strsplit(shortnames[[i]], split = "def|disp|cluster|.csv")[[1]][3]
  cluster[[i]] <- strsplit(shortnames[[i]], split = "def|disp|cluster|.csv")[[1]][4]
}

for(i in 1:length(output.list)){
  deforest <- rep(def[[i]], nrow(output.list[[1]]))
  dispersion <- rep(disp[[i]], nrow(output.list[[1]]))
  clusters <- rep(cluster[[i]], nrow(output.list[[1]]))
  output.list[[i]] <- cbind(output.list[[i]], deforest, dispersion, clusters)
}
head(output.list[[1]])

# Turn list into big-ass data frame
output.mat <- do.call(rbind, output.list)

# Histograms -------------------------------------
# Pull out data from final time step
step.final <- subset(output.mat, output.mat$Time == 999)

# histogram time steps
histo3 <- ggplot(data = step.final[which(step.final$clusters==0.3),], 
                 aes(PercInf*100)) +
  geom_histogram(fill = "darkgrey", bins = 15) +
  facet_grid(vars(deforest), vars(dispersion)) +
  theme_classic(base_size = 18) +
  labs(x="% Rust Infection", y="Frequency")+
  theme(axis.text.y = element_blank())

histo2 <- ggplot(data = step.final[which(step.final$clusters==0.2),], 
                 aes(PercInf*100)) +
  geom_histogram(fill = "darkgrey", bins = 15) +
  facet_grid(vars(deforest), vars(dispersion)) +
  theme_classic(base_size = 18) +
  labs(x="% Rust Infection", y="Frequency")+
  theme(axis.text.y = element_blank())

histo1 <- ggplot(data = step.final[which(step.final$clusters==0.1),], 
                 aes(PercInf*100)) +
  geom_histogram(fill = "darkgrey", bins = 15) +
  facet_grid(vars(deforest), vars(dispersion)) +
  theme_classic(base_size = 18) +
  labs(x="% Rust Infection", y="Frequency")+
  theme(axis.text.y = element_blank())

# ggsave(histo1, filename = "hist1.jpeg", width = 8.5, height = 6)
# ggsave(histo2, filename = "hist2.jpeg", width = 8.5, height = 6)
# ggsave(histo3, filename = "hist3.jpeg", width = 8.5, height = 6)

# Heat Maps ---------------------------------------------
# Calculate Pearson Skewness Coefficient
step.final %>%
  group_by(deforest, dispersion, clusters) %>%
  summarise(skew = skewness(PercInf)) %>%
  {. ->> data.skew}

# Shape data into list
skew.list <- list()
for(i in 1:length(unique(data.skew$clusters))){
  skew.list[[i]] <- data.skew[which(data.skew$clusters==unique(data.skew$clusters[i])),]
}

heatplot1 <- ggplot(skew.list[[1]], aes(deforest, dispersion, fill = skew)) + 
  geom_raster(hjust = 0, vjust = 0)+
  scale_fill_viridis(name = "Skew", limits = c(0, 3))+
  labs(x = "Deforestation (%)", y = "Dispersion")+
  scale_x_discrete(expand = c(0,0))+
  scale_y_discrete(expand = c(0,0))

heatplot2 <- ggplot(skew.list[[2]], aes(deforest, dispersion, fill = skew)) + 
  geom_raster(hjust = 0, vjust = 0)+
  scale_fill_viridis(name = "Skew", limits = c(0, 3))+
  labs(x = "Deforestation (%)", y = "Dispersion")+
  scale_x_discrete(expand = c(0,0))+
  scale_y_discrete(expand = c(0,0))

heatplot3 <- ggplot(skew.list[[3]], aes(deforest, dispersion, fill = skew)) + 
  geom_raster(hjust = 0, vjust = 0)+
  scale_fill_viridis(name = "Skew", limits = c(0, 3))+
  labs(x = "Deforestation (%)", y = "Dispersion")+
  scale_x_discrete(expand = c(0,0))+
  scale_y_discrete(expand = c(0,0))

# Compare heatplots: where do they differ most? -----------------------------
# Change each element of skew.list into matrix

# Use raster functions to compare matrices

# Does starting location matter -------------------------------
quant90 %>%
  filter(X, Y) %>%
  group_by(X, Y) %>%
  summarise(Count = n()) %>%
  {. ->> loc.counts}

step500 %>%
  filter(X, Y) %>%
  group_by(X, Y) %>%
  summarise(Count = n()) %>%
  {. ->> all.loc}

spatial <- ggplot(data = loc.counts, aes(x = X, y = Y, fill = Count))+
  geom_raster()+
  xlim(0,99)+
  ylim(0,99)+
  theme_classic(base_size = 18)+
  theme(legend.position = "none")

spatial_all <- ggplot(data = all.loc, aes(x = X, y = Y, fill = Count))+
  geom_raster()+
  xlim(0,99)+
  ylim(0,99)+
  theme_classic(base_size = 18)+
  theme(legend.position = "none")

# Saving Plots ---------------------------------------------------
# Histogram
ggsave("histrust.jpeg", histo)

# Line plot
ggsave("rustinfec.jpeg", avg.lines)

# Boxplots
ggsave("rust_infection_all.png", infection.all)
ggsave("rust_infection.png", infection)
ggsave("rust_dispersion.png", dispersion)
ggsave("rust_deforestation.png", deforestation)

ggsave("lotsaboxes.jpeg", many.boxes, width = 13, height = 7, units = "in")

# Heat plots
ggsave("heatplotmean.jpeg", heatplotmean)
ggsave("heatplotmedian.jpeg", heatplotmedian)

# Raster of infection start
ggsave("infecstart.jpeg", spatial)
  
  
