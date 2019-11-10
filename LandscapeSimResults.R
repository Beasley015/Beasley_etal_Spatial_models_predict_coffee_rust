#####################################################
# Modeling coffee rust movement through a landscape #
# N. Aristizabal, E. Beasley, & E. Bueno            #
# Spring 2019- ????                                 #
#####################################################
library(ggplot2)
library(dplyr)
library(tidyverse)
library(rcompanion)
library(dunn.test)
library(viridis)

# Sample beta distributions (for figures) ----------------------------------
neighbor1 <- rbeta(10000, 2, 8)
neighbor2 <- rbeta(10000, 5, 5)
neighbor3 <- rbeta(10000, 8, 2)

ggplot(mapping = aes(x = neighbor1))+
  geom_density(size = 1)+
  labs(x = "Infection probability", y = "Density")+
  scale_x_continuous(limits = c(0,1), expand = c(0,0))+
  scale_y_continuous(limits = c(0,4.5), expand = c(0,0))+
  theme_bw()+
  theme(panel.grid = element_blank(), axis.text.x = element_text(size = 16),
        axis.text.y = element_blank(), axis.title = element_blank())

ggsave(file = "beta1.jpg")

ggplot(mapping = aes(x = neighbor2))+
  geom_density(size = 1)+
  labs(x = "Infection probability", y = "Density")+
  scale_x_continuous(limits = c(0,1), expand = c(0,0))+
  scale_y_continuous(limits = c(0,4.5), expand = c(0,0))+
  theme_bw()+
  theme(panel.grid = element_blank(), axis.text.x = element_text(size = 16),
        axis.text.y = element_blank(), axis.title = element_blank())

ggsave(file = "beta2.jpg")

ggplot(mapping = aes(x = neighbor3))+
  geom_density(size = 1)+
  labs(x = "Infection probability", y = "Density")+
  scale_x_continuous(limits = c(0,1), expand = c(0,0))+
  scale_y_continuous(limits = c(0,4.5), expand = c(0,0))+
  theme_bw()+
  theme(panel.grid = element_blank(), axis.text.x = element_text(size = 16),
        axis.text.y = element_blank(), axis.title = element_blank())

ggsave(file = "beta3.jpg")


# Read in model outputs --------------------------------------
# Read all csv's into a list
filenames <- list.files("Outputs", pattern = "*.csv", full.names = T)
shortnames <- list.files("Outputs", pattern = "*.csv")
output.list <- lapply(filenames, read.csv, header = F)

# Rename columns of each dataframe
newnames <- c("Time", "PercInf","X","ResistProb")
output.list <- lapply(output.list, setNames, newnames)

# Set up data -----------------------------------------------
# Add column to denote replicates
replicate <- logical()
for(i in 1:50){
  new <- rep(i, 500)
  replicate <- append(replicate, new)
}
  
output.list <- lapply(output.list, cbind, replicate)  

# Pull deforestation and dispersion from file names
loop.ready <- c(1:length(shortnames))

# Need to fix this for loop- not calling the correct characters
for(i in loop.ready) {
  def[[i]] <- strsplit(shortnames[[i]], split = "def|disp|prob|.csv")[[1]][2]
  disp[[i]] <- strsplit(shortnames[[i]], split = "def|disp|prob|.csv")[[1]][3]
}

for(i in 1:length(output.list)){
  deforest <- rep(def[[i]], nrow(output.list[[1]]))
  dispersion <- rep(disp[[i]], nrow(output.list[[1]]))
  output.list[[i]] <- cbind(output.list[[i]], deforest, dispersion)
}
head(output.list[[1]])

# Turn list into big-ass data frame
output.mat <- do.call(rbind, output.list)

# Plots-------------------------------------

full.infec <- subset(output.mat, output.mat$PercInf == 1)
  
full.infec %>%
  group_by(ResistProb, replicate, deforest, dispersion) %>%
  summarise(newTime = min(Time)) %>%
  {. ->> full.infec2}

half.infec <- subset(output.mat, output.mat$PercInf <= 0.5)

half.infec %>%
  group_by(ResistProb, replicate, deforest, dispersion) %>%
  summarise(newTime = max(Time)) %>%
  {. ->> half.infec2}
  
# plot variation in Percentage Infestation 
deforestation <- ggplot(data1000, aes(x = deforest, y= newTime)) +
  geom_boxplot(fill = "forestgreen") +
  labs(x="Deforestation (%)",
       y="Leaf rust infection (%)") + 
  theme_classic()

# individual percent infestation ~ dispersion
dispersion <- ggplot(data1000, aes(x = dispersion, y= PercInf)) +
  geom_boxplot(fill = "deepskyblue3") +
  labs(x="Degree of dispersion",
       y="Leaf rust infection (%)") + 
  theme_classic()


# individual percent infestation ~ deforestation among replicates
# ggplot(output.mat, aes(x = deforest, y= PercInf)) +
  #geom_boxplot(aes(fill=factor(deforest)))

full.infec2 %>%
  group_by(deforest, dispersion) %>%
  summarise(mean = mean(newTime), median = median(newTime)) %>%
  {. ->> datameans}

# heat map
heatplotmean <- ggplot(datameans, aes(deforest, dispersion, fill = mean)) + 
  geom_raster(hjust = 0, vjust = 0)+
  scale_fill_viridis(name = "Leaf Rust Infection (%)")+
  labs(x = "Deforestation (%)", y = "Degree of Dispersion")+
  scale_x_discrete(expand = c(0,0))+
  scale_y_discrete(expand = c(0,0))

heatplotmedian <- ggplot(datameans, aes(deforest, dispersion, fill = median)) + 
  geom_raster(hjust = 0, vjust = 0)+
  scale_fill_viridis(name = "Leaf Rust Infection (%)")+
  labs(x = "Deforestation (%)", y = "Degree of Dispersion")+
  scale_x_discrete(expand = c(0,0))+
  scale_y_discrete(expand = c(0,0))

# Infection through time -------------------------------
# percent infestation through time steps
ggplot(output.mat, aes(x = Time, y = PercInf, group = as.factor(replicate), 
                       alpha = 0.3)) + 
  geom_line() +
  facet_grid(vars(deforest), vars(dispersion)) +
  theme_classic() +
  labs(x="Time",
     y="Leaf rust infection (%)") +
  theme_classic() +
  theme(legend.position = "none")
  
 
# ploting % infestation through time steps averaging replicates
dat <- as_tibble(output.mat) %>%
  group_by(Time, deforest, dispersion) %>%
  filter(Time <= 150) %>%
  summarise(m = mean(PercInf), sd = sd(PercInf)) 

# check that it worked
glimpse(dat)

ggplot(data = dat) +
  geom_ribbon(aes(x = Time, ymin = (m-sd), ymax = (m+sd)), fill = "grey70", 
              alpha = 0.6) +
  geom_line(aes(x = Time, y = m)) +
  facet_grid(vars(deforest), vars(dispersion)) +
  theme_classic() +
  labs(x="Time",
       y="Leaf rust infection (%)") +
  theme_classic() +
  theme(legend.position = "none")
                
# histograms time steps
ggplot(data = full.infec2, aes(newTime)) +
  geom_histogram(fill = "darkgrey") +
  facet_grid(vars(deforest), vars(dispersion)) +
  theme_classic() +
  labs(x="Time", y="")

ggplot(data = half.infec2, aes(newTime)) +
  geom_histogram(fill = "darkgrey") +
  facet_grid(vars(deforest), vars(dispersion)) +
  theme_classic() +
  labs(x="Time", y="")

# saving plots
ggsave("rust_infection_all.png", infection.all)
ggsave("rust_infection.png", infection)
ggsave("rust_dispersion.png", dispersion)
ggsave("rust_deforestation.png", deforestation)

ggsave("heatplotmean.jpeg", heatplotmean)
ggsave("heatplotmedian.jpeg", heatplotmedian)

ggsave("hist.rust.jpeg", hist.rust)
  
  
  
