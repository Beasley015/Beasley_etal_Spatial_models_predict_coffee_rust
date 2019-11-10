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
mid <- function(text, start_num, num_char){
  substr(text, start_num, start_num + num_char - 1)
}

def <- list()
disp <- list()

loop.ready <- c(1:length(shortnames))

# Need to fix this for loop- not calling the correct characters
# for(i in loop.ready) {
#   def[[i]] <- mid(shortnames[[i]],4,2) 
#   disp[[i]] <- mid(shortnames[[i]],10,1)
# }

for(i in 1:length(output.list)){
  deforest <- rep(def[[i]], nrow(output.list[[1]]))
  dispersion <- rep(disp[[i]], nrow(output.list[[1]]))
  output.list[[i]] <- cbind(output.list[[i]], deforest, dispersion)
}
head(output.list[[1]])

# Turn list into big-ass data frame
output.mat <- do.call(rbind, output.list)

# Plots (this is where the big edits are)-------------------------------------

data500 <- subset(output.mat, output.mat$Time == 500)
  
# plot variation in Percentage Infestation 
deforestation <- ggplot(data1000, aes(x = deforest, y= PercInf)) +
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

#Try an anova
anov <- aov(data = data1000, PercInf ~ deforest + dispersion + deforest*dispersion)
summary(anov)
#deforestation is significant, nothing else is.

#Nonparametric test
nonpar <- scheirerRayHare(PercInf ~ deforest + dispersion, data = data1000)
#same as anova

dunn <- dunn.test(x = data1000$PercInf, g = data1000$deforest)

# individual percent infestation ~ deforestation among replicates
# ggplot(output.mat, aes(x = deforest, y= PercInf)) +
  #geom_boxplot(aes(fill=factor(deforest)))

data1000 %>%
  group_by(deforest, dispersion) %>%
  summarise(mean = mean(PercInf), median = median(PercInf)) %>%
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

# percent infestation through time steps
infection.all <- ggplot(output.mat, aes(x = Time, y = PercInf, group = as.factor(replicate), alpha = 0.3)) + 
  geom_line() +
  facet_grid(vars(deforest), vars(dispersion)) +
  theme_classic() +
  labs(x="Time",
     y="Leaf rust infection (%)") +
  theme_classic() +
  theme(legend.position = "none")
  
 
# ploting % infestation through time steps averaging replicates
data1 <- as.tibble(output.mat) %>%
  group_by(Time, deforest, dispersion) %>%
  summarise(m = median(PercInf), sd = sd(PercInf)) 

# check that it worked
glimpse(data)

infection <- ggplot(data = data1) +
  geom_ribbon(aes(x = Time, ymin = (m-sd), ymax = (m+sd)), fill = "grey70", alpha = 0.6) +
  geom_line(aes(x = Time, y = m)) +
  facet_grid(vars(deforest), vars(dispersion)) +
  theme_classic() +
  labs(x="Time",
       y="Leaf rust infection (%)") +
  theme_classic() +
  theme(legend.position = "none")
                
# histograms time steps
hist.rust <- ggplot(data = output.mat, aes(PercInf)) +
  geom_histogram(binwidth = 0.15, fill = "darkgrey") +
  facet_grid(vars(deforest), vars(dispersion)) +
  theme_classic() +
  labs(x="Leaf rust infection (%)",
       y="")

# saving plots
ggsave("rust_infection_all.png", infection.all)
ggsave("rust_infection.png", infection)
ggsave("rust_dispersion.png", dispersion)
ggsave("rust_deforestation.png", deforestation)

ggsave("heatplotmean.jpeg", heatplotmean)
ggsave("heatplotmedian.jpeg", heatplotmedian)

ggsave("hist.rust.jpeg", hist.rust)
  
  
  
