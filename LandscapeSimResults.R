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
library(raster)

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
# filenames <- list.files("Outputs_buffered", pattern = "*.csv", 
#                         full.names = T)
# shortnames <- list.files("Outputs_buffered", pattern = "*.csv")
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

# Read in output matrix ------------------------

output.mat <- readRDS("outputmat.rds")

# Pull out data from final time step
step.final <- subset(output.mat, output.mat$Time == 999)

# Get average rate of spread (new cells per time step)
avg.coffee <- data.frame(coff = as.character(c(0.1, 0.2, 0.3)),
                         ncells = c(3605, 3610, 3630))

rates <- output.mat %>%
  left_join(y = avg.coffee, by = "coff") %>%
  group_by(replicate, deforest, dispersion, coff) %>%
  mutate(inf.cells = PercInf/(1/ncells)) %>%
  mutate(NewInf = inf.cells-lag(inf.cells)) %>%
  summarise(mean.rate = mean(NewInf, na.rm = T))

# Make 0's a ridiculously small number so functions work
step.final$PercInf[which(step.final$PercInf == 0)] <- 0.00000001
rates$mean.rate[which(rates$mean.rate == 0)] <- 0.000000000000001

# Summary stats ---------------------------------------
# Histogram of all outcomes
hist.all <- ggplot(data = rates)+
  geom_histogram(aes(x = mean.rate), fill = 'lightgray', color = 'black', 
                 boundary = 0)+
  scale_x_continuous(expand = c(0,0))+
  scale_y_continuous(expand = c(0,0), limits = c(0,650))+
  labs(x = "Avg. Infections per Time Step", y = "Count")+
  theme_bw(base_size = 12)+
  theme(panel.grid = element_blank())

# ggsave(hist.all, filename = 'histall.tiff', dpi = 1200, width = 84,
#        units = 'mm', height = 60)

# Range of prevalence values
range(rates$mean.rate)

# Get parameters of distribution across all parameter combinations
gamma.all <- fitdist(rates$mean.rate, "gamma")

# Calculate expected value
gammaexpec <- function(x){
  a <- x$estimate[1]
  b <- x$estimate[2]
  
  expec <- a/b
  
  return(expec)
}

gammaexpec(gamma.all)

# Calculate skew
gammaskew <- function(x){
  a <- x$estimate[1]
  b <- x$estimate[2]
  
  skew <- 2/sqrt(a)
  
  return(skew)
}

gammaskew(gamma.all)

# Calculate kurtosis
gammakurtosis <- function(x){
  a <- x$estimate[1]
  b <- x$estimate[2]
  
  k <- 6/a
  
  return(k)
}

gammakurtosis(gamma.all)

# Histograms -------------------------------------
# histogram time steps
histo3 <- ggplot(data = rates[which(rates$coff==0.3),], 
                 aes(mean.rate)) +
  geom_histogram(fill = "darkgrey", bins = 15) +
  lims(x = c(NA,2))+
  facet_grid(vars(deforest), vars(dispersion)) +
  theme_classic(base_size = 18) +
  labs(x="New Infections per Time Step", y="Frequency")+
  theme(axis.text.y = element_blank())

histo2 <- ggplot(data = rates[which(rates$coff==0.2),], 
                 aes(mean.rate)) +
  geom_histogram(fill = "darkgrey", bins = 15) +
  lims(x = c(NA,2))+
  facet_grid(vars(deforest), vars(dispersion)) +
  theme_classic(base_size = 18) +
  labs(x="New Infections per Time Step", y="Frequency")+
  theme(axis.text.y = element_blank())

histo1 <- ggplot(data = rates[which(rates$coff==0.1),], 
                 aes(mean.rate)) +
  geom_histogram(fill = "darkgrey", bins = 15) +
  lims(x = c(NA,2))+
  facet_grid(vars(deforest), vars(dispersion)) +
  theme_classic(base_size = 18) +
  labs(x="New Infections per Time Step", y="Frequency")+
  theme(axis.text.y = element_blank())

# ggsave(histo1, filename = "hist1.jpeg", width = 8.5, height = 6)
# ggsave(histo2, filename = "hist2.jpeg", width = 8.5, height = 6)
# ggsave(histo3, filename = "hist3.jpeg", width = 8.5, height = 6)

# Boxplots for each clustering value --------------
# Max infection
rates %>%
  group_by(deforest, dispersion, coff) %>%
  summarise(max = max(mean.rate)) %>%
  {. ->> max.final}

maxbox <- ggplot(data = max.final, aes(x = factor(coff), y = max))+
  geom_boxplot(fill = 'lightgray')+
  labs(x = "Coffee Clustering", y = "Maximum Infection Rate")+
  scale_x_discrete(labels = c('0.1' = "Low", '0.2' = "Mid", '0.3' = "High"))+
  theme_bw(base_size = 12)+
  theme(panel.grid = element_blank())

# ggsave(maxbox, file = "MaxInfBoxes.jpeg")

# Expected value of the distribution
gamma.final <- rates %>%
  group_by(deforest, dispersion, coff) %>%
  mutate(alpha = fitdist(mean.rate, "gamma")$estimate[1]) %>%
  mutate(beta = fitdist(mean.rate, "gamma")$estimate[2]) %>%
  dplyr::select(deforest:beta) %>%
  distinct()

exp.final <- gamma.final %>%
  mutate(Expected.Value = alpha/beta)
  
expected.boxes <- ggplot(exp.final, aes(x = coff, y = Expected.Value))+
  geom_boxplot(fill = 'lightgray')+
  labs(y = "Expected Value", x = "Coffee Clustering")+
  scale_x_discrete(labels = c('0.1' = "Low", '0.2' = "Mid", 
                              '0.3' = "High"))+
  theme_bw(base_size = 12)+
  theme(panel.grid = element_blank())

# ggsave(expected.boxes, filename = 'expecplot.jpeg')

# Calculate skew
skew.final <- gamma.final %>%
  mutate(skew= 2/sqrt(alpha))

skew.plot <- ggplot(data = skew.final, aes(x = coff, y = skew))+
  geom_boxplot(fill = 'lightgray')+
  labs(x = "Coffee Clustering", y = "Skew")+
  scale_x_discrete(labels = c('0.1' = "Low", '0.2' = "Mid", '0.3' = "High"))+
  theme_bw(base_size = 12)+
  theme(panel.grid = element_blank())

# ggsave(skew.plot, filename = 'skewplot.jpeg')

# Kurtosis parameter
k.final <- gamma.final %>%
  mutate(k=6/alpha)

k.plot <- ggplot(data = k.final, aes(x = coff, y = k))+
  geom_boxplot(fill = 'lightgray')+
  labs(x = "Coffee Clustering", y = "Kurtosis")+
  scale_x_discrete(labels = c('0.1' = "Low", '0.2' = "Mid", '0.3' = "High"))+
  theme_bw(base_size = 12)+
  theme(panel.grid = element_blank())

# ggsave(k.plot, filename = "kurtosisplot.jpeg")

# Put all boxplots together
megabox <- expected.boxes+ggtitle("a)")+
  maxbox+ggtitle("b)")+
  skew.plot+ggtitle("c)")+
  k.plot+ggtitle("d)")

# ggsave(megabox, file = 'megabox.tiff', width = 174, height = 130,
#        units = 'mm', dpi = 1200)

# Closer look at deforestation/dispersion ----------------------
# Expected values ##
exp.hi <- filter(exp.final, coff == '0.3')
exp.mid <- filter(exp.final, coff == '0.2')
exp.lo <- filter(exp.final, coff == '0.1')

# No clear patterns w/expected value at high clustering
ggplot(data = exp.hi, aes(x = deforest, y = Expected.Value))+
  geom_point(aes(color = dispersion), size = 2)+
  geom_smooth(se = F, method = 'lm')+
  labs(x = '% Deforestation', y = 'Expected Value')+
  scale_color_viridis_d(name = "Dispersion")+
  theme_bw(base_size = 16)+
  theme(panel.grid = element_blank())

cor.test(x = as.numeric(exp.hi$deforest), y = exp.hi$Expected.Value, 
         method = "spearman")
cor.test(x = as.numeric(exp.hi$dispersion), y = exp.hi$Expected.Value, 
         method = "spearman")

# No clear patterns at mid clustering
ggplot(data = exp.mid, aes(x = deforest, y = Expected.Value, 
                           color = dispersion))+
  geom_point(size = 2)+
  labs(x = '% Deforestation', y = 'Expected Value')+
  scale_color_viridis_d(name = "Dispersion")+
  theme_bw(base_size = 18)+
  theme(panel.grid = element_blank())

cor.test(x = as.numeric(exp.mid$deforest), y = exp.mid$Expected.Value, 
         method = "spearman")
cor.test(x = as.numeric(exp.mid$dispersion), y = exp.mid$Expected.Value, 
         method = "spearman")

# Ditto low clustering
ggplot(data = exp.lo, aes(x = deforest, y = Expected.Value, color = dispersion))+
  geom_point(size = 2)+
  labs(x = '% Deforestation', y = 'Expected Value')+
  scale_color_viridis_d(name = "Dispersion")+
  theme_bw(base_size = 18)+
  theme(panel.grid = element_blank())

cor.test(x = as.numeric(exp.lo$deforest), y = exp.lo$Expected.Value, 
         method = "spearman")
cor.test(x = as.numeric(exp.lo$dispersion), y = exp.lo$Expected.Value, 
         method = "spearman")

# Maximum infection ##
max.lo <- filter(max.final, coff == '0.1')
max.mid <- filter(max.final, coff == '0.2')
max.hi <- filter(max.final, coff == '0.3')

# No clear patterns at low clustering
ggplot(data = max.lo, aes(x = deforest, y = max, color = dispersion))+
  geom_point(size = 2)+
  labs(x = "% Deforestation", y = "Maximum Infection")+
  scale_color_viridis_d(name = "Dispersion")+
  theme_bw(base_size = 18)+
  theme(panel.grid = element_blank())

cor.test(x = as.numeric(max.lo$deforest), y = max.lo$max, method = "spearman")
cor.test(x = as.numeric(max.lo$dispersion), y = max.lo$max, method = "spearman")

# No clear patterns at mid clustering
ggplot(data = max.mid, aes(x = deforest, y = max, color = dispersion))+
  geom_point(size = 2)+
  labs(x = "% Deforestation", y = "Maximum Infection")+
  scale_color_viridis_d(name = "Dispersion")+
  theme_bw(base_size = 18)+
  theme(panel.grid = element_blank())

cor.test(x = as.numeric(max.mid$deforest), y = max.mid$max, method = "spearman")
cor.test(x = as.numeric(max.mid$dispersion), y = max.mid$max, method = "spearman")

# High dispersion = hi max rate at high clustering
max.hi.plot <- ggplot(data = max.hi, aes(x = dispersion, y = max, 
                                         color = deforest))+
  geom_point(size = 1.5)+
  labs(x = "Dispersion", y = "Maximum Rate of Spread")+
  scale_color_viridis_d(name = "Deforestation")+
  theme_bw(base_size = 10)+
  theme(panel.grid = element_blank())

cor.test(x = as.numeric(max.hi$deforest), y = max.hi$max, method = "spearman")
cor.test(x = as.numeric(max.hi$dispersion), y = max.hi$max, method = "spearman")

# ggsave(max.hi.plot, filename = "maxhiplot.tiff", width = 84, height = 60,
#        units = 'mm', dpi = 1200)

# Skew ##
skew.hi <- filter(skew.final, coff == '0.3')
skew.mid <- filter(skew.final, coff == '0.2')
skew.lo <- filter(skew.final, coff == '0.1')

# No clear pattern at high clustering
ggplot(data = skew.hi, aes(x = deforest, y = skew, color = dispersion))+
  geom_point(size = 2)+
  labs(x = "% Deforestation", y = "Skew")+
  scale_color_viridis_d(name = "Dispersion")+
  theme_bw(base_size = 16)+
  theme(panel.grid = element_blank())

cor.test(x = as.numeric(skew.hi$deforest), y = skew.hi$skew, method = "spearman")
cor.test(x = as.numeric(skew.hi$dispersion), y = skew.hi$skew, method = "spearman")

# No clear patterns at mid and low clustering
ggplot(data = skew.mid, aes(x = deforest, y = skew, 
                            color = dispersion))+
  geom_point(size = 2)+
  labs(x = "% Deforestation", y = "Skew")+
  scale_color_viridis_d(name = "Dispersion")+
  theme_bw(base_size = 18)+
  theme(panel.grid = element_blank())

cor.test(x = as.numeric(skew.mid$deforest), y = skew.mid$skew, method = "spearman")
cor.test(x = as.numeric(skew.mid$dispersion), y = skew.mid$skew, method = "spearman")

ggplot(data = skew.lo, aes(x = deforest, y = skew, color = dispersion))+
  geom_point(size = 2)+
  labs(x = "% Deforestation", y = "Skew")+
  scale_color_viridis_d(name = "Dispersion")+
  theme_bw(base_size = 18)+
  theme(panel.grid = element_blank())

cor.test(x = as.numeric(skew.lo$deforest), y = skew.lo$skew, method = "spearman")
cor.test(x = as.numeric(skew.lo$dispersion), y = skew.lo$skew, method = "spearman")

# Kurtosis ##
k.hi <- filter(k.final, coff == '0.3')
k.mid <- filter(k.final, coff == '0.2')
k.lo <- filter(k.final, coff == '0.1')

# No clear patterns at high clustering
ggplot(data = k.hi, aes(x = as.numeric(dispersion), y = k,
                            color = deforest))+
  geom_point(size = 2)+
  labs(x = "Dispersion", y = "Precision (Kappa)")+
  scale_color_viridis_d(name = "Deforestation")+
  theme_bw(base_size = 16)+
  theme(panel.grid = element_blank())

cor.test(x = as.numeric(k.hi$deforest), y = k.hi$k, method = "spearman")
cor.test(x = as.numeric(k.hi$dispersion), y = k.hi$k, method = "spearman")

# No at mid clustering
ggplot(data = k.mid, aes(x = as.numeric(dispersion), y = k, 
                            color = deforest))+
  geom_point(size = 2)+
  labs(x = "Dispersion", y = "Kurtosis")+
  scale_color_viridis_d(name = "% Deforestation")+
  theme_bw(base_size = 18)+
  theme(panel.grid = element_blank())

cor.test(x = as.numeric(k.mid$deforest), y = k.mid$k, method = "spearman")
cor.test(x = as.numeric(k.mid$dispersion), y = k.mid$k, 
         method = "spearman")

# No clear patterns at low clustering
ggplot(data = k.lo, aes(x = as.numeric(deforest), y = k, 
                           color = dispersion))+
  geom_point(size = 2)+
  labs(x = "% Deforestation", y = "Kurtosis")+
  scale_color_viridis_d(name = "Dispersion")+
  theme_bw(base_size = 12)+
  theme(panel.grid = element_blank())

cor.test(x = as.numeric(k.lo$deforest), y = k.lo$k, method = "spearman")
cor.test(x = as.numeric(k.lo$dispersion), y = k.lo$k, method = "spearman")

# Does starting location matter -------------------------------
all.loc <- step.final %>%
  dplyr::select(X,Y) %>%
  group_by(X,Y) %>%
  mutate(count = n())

hi.loc <- step.final %>%
  group_by(coff) %>%
  filter(PercInf > quantile(PercInf, 0.75)) %>%
  mutate(quant = "Top25")

lo.loc <- step.final %>%
  group_by(coff) %>%
  filter(PercInf < quantile(PercInf, 0.25)) %>%
  mutate(quant = "Bottom25")

test.loc <- rbind(hi.loc, lo.loc) %>%
  group_by(X,Y,quant) %>%
  summarise(count = n()) %>%
  mutate(count = ifelse(quant == "Bottom25", -count, count)) %>%
  ungroup() %>%
  group_by(X,Y) %>%
  mutate(val = sum(count))

# Highest outbreaks tend to start in the middle of the landscape
start.coords <- ggplot(data = test.loc, aes(x = X, y = Y, fill = factor(val)))+
  geom_raster()+
  scale_fill_brewer(type = 'div', name = "")+
  scale_x_discrete(expand = c(0,0))+
  scale_y_discrete(expand = c(0,0))+
  theme_classic(base_size = 18)

# ggsave(start.coords, filename = "startcoords.jpeg")

# Is there any clustering among all starting locations?
all.start <- ggplot(data = all.loc, aes(x = X, y = Y, fill = count))+
  geom_raster()+
  scale_x_discrete(expand = c(0,0))+
  scale_y_discrete(expand = c(0,0))+
  theme_classic(base_size = 18)
# Looks fairly even

# Convert data frame to raster
loc.raster <- rasterFromXYZ(all.loc)

# Look for autocorrelation
Moran(loc.raster$count)
# Close to 0
