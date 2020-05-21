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
library(fitdistrplus)

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
# # Read all csv's into a list
# filenames <- list.files("Outputs", pattern = "*.csv", full.names = T)
# shortnames <- list.files("Outputs", pattern = "*.csv")
# output.list <- lapply(filenames, read.csv, header = F)
# 
# # Rename columns of each dataframe
# newnames <- c("Time", "PercInf","X","Y","Nope")
# output.list <- lapply(output.list, setNames, newnames)
# output.list <- lapply(output.list, function(x) x[!(names(x)) %in% "Nope"])
# 
# # Set up data -----------------------------------------------
# # Add column to denote replicates
# replicate <- logical()
# for(i in 1:50){
#   new <- rep(i, 1000)
#   replicate <- append(replicate, new)
# }
# 
# output.list <- lapply(output.list, cbind, replicate)
# 
# # Pull deforestation and dispersion from file names
# loop.ready <- c(1:length(shortnames))
# def <- list()
# disp <- list()
# kaffee <- list()
# 
# for(i in loop.ready) {
#   def[[i]] <- strsplit(shortnames[[i]], split = "def|disp|cluster|.csv")[[1]][2]
#   disp[[i]] <- strsplit(shortnames[[i]], split = "def|disp|cluster|.csv")[[1]][3]
#   kaffee[[i]] <- strsplit(shortnames[[i]], split = "def|disp|cluster|.csv")[[1]][4]
# }
# 
# for(i in 1:length(output.list)){
#   deforest <- rep(def[[i]], nrow(output.list[[1]]))
#   dispersion <- rep(disp[[i]], nrow(output.list[[1]]))
#   coff <- rep(kaffee[[i]], nrow(output.list[[1]]))
#   output.list[[i]] <- cbind(output.list[[i]], deforest, dispersion, coff)
# }
# head(output.list[[1]])
# 
# # Turn list into big-ass data frame
# output.mat <- do.call(rbind, output.list)
# 
# # Save it as a compressed file
# saveRDS(output.mat, file = "outputmat.rds")

# ------------------------

output.mat <- readRDS("outputmat.rds")

# Pull out data from final time step
step.final <- subset(output.mat, output.mat$Time == 999)

# Summary stats ---------------------------------------
# Histogram of all outcomes
hist.all <- ggplot(data = step.final)+
  geom_histogram(aes(x = PercInf), fill = 'lightgray', color = 'black', 
                 boundary = 0)+
  scale_x_continuous(expand = c(0,0))+
  scale_y_continuous(expand = c(0,0), limits = c(0, 1100))+
  labs(x = "Rust Prevalence", y = "Count")+
  theme_bw(base_size = 18)+
  theme(panel.grid = element_blank())

# ggsave(hist.all, filename = 'histall.jpeg')
  
# Range of prevalence values
range(step.final$PercInf)

# Get parameters of distribution across all parameter combinations
beta.all <- fitdist(step.final$PercInf, "beta")

# Calculate expected value
betaexpec <- function(x){
  a <- x$estimate[1]
  b <- x$estimate[2]
  
  expec <- (a)/(a+b)
  
  return(expec)
}

betaexpec(beta.all)

# Calculate skew
betaskew <- function(x){
  a <- x$estimate[1]
  b <- x$estimate[2]
  
  skew <- (2*(b-a)*sqrt(a+b+1))/((a+b+2)*sqrt(a*b))
  
  return(skew)
}

betaskew(beta.all)

# Calculate concentration parameter
betakappa <- function(x){
  a <- x$estimate[1]
  b <- x$estimate[2]
  
  kap <- a + b
  
  return(kap)
}

betakappa(beta.all)

# Histograms -------------------------------------
# histogram time steps
histo3 <- ggplot(data = step.final[which(step.final$coff==0.3),], 
                 aes(PercInf*100)) +
  geom_histogram(fill = "darkgrey", bins = 15) +
  facet_grid(vars(deforest), vars(dispersion)) +
  lims(x = c(NA,55))+
  theme_classic(base_size = 18) +
  labs(x="% Rust Infection", y="Frequency")+
  theme(axis.text.y = element_blank())

histo2 <- ggplot(data = step.final[which(step.final$coff==0.2),], 
                 aes(PercInf*100)) +
  geom_histogram(fill = "darkgrey", bins = 15) +
  facet_grid(vars(deforest), vars(dispersion)) +
  lims(x = c(NA,55))+
  theme_classic(base_size = 18) +
  labs(x="% Rust Infection", y="Frequency")+
  theme(axis.text.y = element_blank())

histo1 <- ggplot(data = step.final[which(step.final$coff==0.1),], 
                 aes(PercInf*100)) +
  geom_histogram(fill = "darkgrey", bins = 15) +
  facet_grid(vars(deforest), vars(dispersion)) +
  lims(x = c(NA,55))+
  theme_classic(base_size = 18) +
  labs(x="% Rust Infection", y="Frequency")+
  theme(axis.text.y = element_blank())

# ggsave(histo1, filename = "hist1.jpeg", width = 8.5, height = 6)
# ggsave(histo2, filename = "hist2.jpeg", width = 8.5, height = 6)
# ggsave(histo3, filename = "hist3.jpeg", width = 8.5, height = 6)



# Boxplots for each clustering value --------------
# Max infection
step.final %>%
  group_by(deforest, dispersion, coff) %>%
  summarise(max = max(PercInf)) %>%
  {. ->> max.final}

maxbox <- ggplot(data = max.final, aes(x = factor(coff), y = max))+
  geom_boxplot(fill = 'lightgray')+
  labs(x = "Coffee Clustering", y = "Maximum % Infected")+
  scale_x_discrete(labels = c('0.1' = "Low", '0.2' = "Mid", '0.3' = "High"))+
  theme_bw(base_size = 18)+
  theme(panel.grid = element_blank())

# ggsave(maxbox, file = "MaxInfBoxes.jpeg")

# Expected value of the distribution
beta.final <- step.final %>%
  group_by(deforest, dispersion, coff) %>%
  mutate(alpha = fitdist(PercInf, "beta")$estimate[1]) %>%
  mutate(beta = fitdist(PercInf, "beta")$estimate[2]) %>%
  dplyr::select(deforest:beta) %>%
  distinct()

exp.final <- beta.final %>%
  mutate(Expected.Value = (alpha)/(alpha+beta))
  
expected.boxes <- ggplot(exp.final, aes(x = coff, y = Expected.Value))+
  geom_boxplot(fill = 'lightgray')+
  labs(y = "Expected Value", x = "Coffee Clustering")+
  scale_x_discrete(labels = c('0.1' = "Low", '0.2' = "Mid", '0.3' = "High"))+
  theme_bw(base_size = 18)+
  theme(panel.grid = element_blank())

# ggsave(expected.boxes, filename = 'expecplot.jpeg')

# Calculate skew
skew.final <- beta.final %>%
  mutate(skew=((2*(beta-alpha)*sqrt(alpha+beta+1))/
                 ((alpha+beta+2)*sqrt(alpha*beta))))

skew.plot <- ggplot(data = skew.final, aes(x = coff, y = skew))+
  geom_boxplot(fill = 'lightgray')+
  labs(x = "Coffee Clustering", y = "Skew")+
  scale_x_discrete(labels = c('0.1' = "Low", '0.2' = "Mid", '0.3' = "High"))+
  theme_bw(base_size = 18)+
  theme(panel.grid = element_blank())

# ggsave(skew.plot, filename = 'skewplot.jpeg')

# Concentration parameter
conc.final <- beta.final %>%
  mutate(kappa = alpha + beta)

conc.plot <- ggplot(data = conc.final, aes(x = coff, y = kappa))+
  geom_boxplot(fill = 'lightgray')+
  labs(x = "Coffee Clustering", y = "Kappa")+
  scale_x_discrete(labels = c('0.1' = "Low", '0.2' = "Mid", '0.3' = "High"))+
  theme_bw(base_size = 18)+
  theme(panel.grid = element_blank())

# ggsave(conc.plot, filename = "spreadplot.jpeg")

# Put all boxplots together
megabox <- expected.boxes+ggtitle("A)")+
  maxbox+ggtitle("B)")+
  skew.plot+ggtitle("C)")+
  conc.plot+ggtitle("D)")

# ggsave(megabox, file = 'megabox.jpeg', width = 8, height = 6.5, units = 'in')

# Heat Maps ---------------------------------------------
# Write master function
# Which is currently not working
makin.heatmaps <- function(x, nplots, param){
  plotlist <- list()
  cluster <- c(0.1, 0.2, 0.3)
  for(i in 1:nplots){
    plotlist[[i]] <- ggplot(subset(x, coff %in% cluster[i]), 
           aes(deforest, dispersion, fill = param)) + 
      geom_raster(hjust = 0, vjust = 0)+
      labs(x = "Deforestation (%)", y = "Dispersion")+
      scale_x_discrete(expand = c(0,0))+
      scale_y_discrete(expand = c(0,0))+
      theme_classic(base_size = 18)
  }
  return(plotlist)
}

# Layout for patchwork grouping
layout <- 
  "AABB
   #CC#"

# Expected values
exp.heats <- makin.heatmaps(x = exp.final, nplots = 3, param = Expected.Value)

exp.heats[[1]]+ggtitle('A)')+
  exp.heats[[2]]+ggtitle('B)')+
  exp.heats[[3]]+ggtitle('C)')+
  plot_layout(design = layout, guides = 'collect')&
  scale_fill_viridis_c(limits = range(exp.final$Expected.Value),
                       name = "Expected Value")

# Max values
max.heats <- makin.heatmaps(x = max.final, nplots = 3, param = max.final$max)

max.heats[[1]]+ggtitle('A)')+
  max.heats[[2]]+ggtitle('B)')+
  max.heats[[3]]+ggtitle('C)')+
  plot_layout(design = layout, guides = 'collect')&
  scale_fill_viridis_c(limits = range(max.final$max),
                       name = "Max % Infction")

# Make skew heatplots
heatplot1 <- ggplot(subset(data.skew, coff %in% 0.1), 
                    aes(deforest, dispersion, fill = skew)) + 
  geom_raster(hjust = 0, vjust = 0)+
  labs(x = "Deforestation (%)", y = "Dispersion")+
  scale_x_discrete(expand = c(0,0))+
  scale_y_discrete(expand = c(0,0))+
  theme_classic(base_size = 18)

heatplot2 <- ggplot(subset(data.skew, coff %in% 0.2), 
                    aes(deforest, dispersion, fill = skew)) + 
  geom_raster(hjust = 0, vjust = 0)+
  labs(x = "Deforestation (%)", y = "Dispersion")+
  scale_x_discrete(expand = c(0,0))+
  scale_y_discrete(expand = c(0,0))+
  theme_classic(base_size = 18)

heatplot3 <- ggplot(subset(data.skew, coff %in% 0.3), 
                    aes(deforest, dispersion, fill = skew)) + 
  geom_raster(hjust = 0, vjust = 0)+
  labs(x = "Deforestation (%)", y = "Dispersion")+
  scale_x_discrete(expand = c(0,0))+
  scale_y_discrete(expand = c(0,0))+
  theme_classic(base_size =  18)


layout <- "
#AA#
BBCC
"

heats <- heatplot1 + ggtitle('A)') + 
  heatplot2 + ggtitle('B)') +
  heatplot3 + ggtitle('C)') +
  plot_layout(design = layout, guides = 'collect')&
  scale_fill_viridis_c(limits = range(data.skew$skew), name = "Skew")

# ggsave(heats, filename = 'heatmaps.jpeg', height = 6.5, width = 9.5)

# Heat plots for max values
maxheat1 <- ggplot(subset(max.final, coff %in% 0.1), 
                   aes(deforest, dispersion, fill = max)) + 
  geom_raster(hjust = 0, vjust = 0)+
  labs(x = "Deforestation (%)", y = "Dispersion")+
  scale_x_discrete(expand = c(0,0))+
  scale_y_discrete(expand = c(0,0))+
  theme_classic(base_size = 18)

maxheat2 <- ggplot(subset(max.final, coff %in% 0.2), 
                   aes(deforest, dispersion, fill = max)) + 
  geom_raster(hjust = 0, vjust = 0)+
  labs(x = "Deforestation (%)", y = "Dispersion")+
  scale_x_discrete(expand = c(0,0))+
  scale_y_discrete(expand = c(0,0))+
  theme_classic(base_size = 18)

maxheat3 <- ggplot(subset(max.final, coff %in% 0.3), 
                   aes(deforest, dispersion, fill = max)) + 
  geom_raster(hjust = 0, vjust = 0)+
  labs(x = "Deforestation (%)", y = "Dispersion")+
  scale_x_discrete(expand = c(0,0))+
  scale_y_discrete(expand = c(0,0))+
  theme_classic(base_size = 18)

maxheats <- maxheat1 + ggtitle('A)') + 
  maxheat2 + ggtitle('B)') +
  maxheat3 + ggtitle('C)') +
  plot_layout(design = layout, guides = 'collect') &
  scale_fill_viridis_c(name = "Max Infection", limits = range(max.final$max))

# Closer look at lowest clustering value ----------------------
lowest.skew <- skew.list[[1]]

ggplot(data = lowest.skew, aes(x = deforest, y = skew, group = dispersion))+
  geom_line(aes(color = dispersion), size = 2)+
  theme_bw(base_size = 18)+
  theme(panel.grid = element_blank())

ggplot(data = lowest.skew, aes(x = deforest, y = skew))+
  geom_violin()+
  theme_bw(base_size = 18)+
  theme(panel.grid = element_blank())

ggplot(data = lowest.skew, aes(x = dispersion, y = skew))+
  geom_violin()+
  theme_bw(base_size = 18)+
  theme(panel.grid = element_blank())

lowest.max <- max.final %>%
  filter(coff == 0.1)

ggplot(data = lowest.max, aes(x = deforest, y = max, group = dispersion))+
  geom_line(aes(color = dispersion), size = 2)+
  theme_bw(base_size = 18)+
  theme(panel.grid = element_blank())

ggplot(data = lowest.max, aes(x = deforest, y = max))+
  geom_violin()+
  theme_bw(base_size = 18)+
  theme(panel.grid = element_blank())

ggplot(data = lowest.max, aes(x = dispersion, y = max))+
  geom_violin()+
  theme_bw(base_size = 18)+
  theme(panel.grid = element_blank())

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
  
  
