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
  theme_bw(base_size = 16)+
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
  scale_x_discrete(labels = c('0.1' = "Low", '0.2' = "Mid", 
                              '0.3' = "High"))+
  theme_bw(base_size = 16)+
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
  theme_bw(base_size = 16)+
  theme(panel.grid = element_blank())

# ggsave(skew.plot, filename = 'skewplot.jpeg')

# Concentration parameter
conc.final <- beta.final %>%
  mutate(kappa = alpha + beta)

conc.plot <- ggplot(data = conc.final, aes(x = coff, y = kappa))+
  geom_boxplot(fill = 'lightgray')+
  labs(x = "Coffee Clustering", y = "Precision (Kappa)")+
  scale_x_discrete(labels = c('0.1' = "Low", '0.2' = "Mid", '0.3' = "High"))+
  theme_bw(base_size = 16)+
  theme(panel.grid = element_blank())

# ggsave(conc.plot, filename = "spreadplot.jpeg")

# Put all boxplots together
megabox <- expected.boxes+ggtitle("A)")+
  maxbox+ggtitle("B)")+
  skew.plot+ggtitle("C)")+
  conc.plot+ggtitle("D)")

# ggsave(megabox, file = 'megabox.jpeg', width = 8, height = 6.5,
#        units = 'in')

# Heat Maps ---------------------------------------------
# Write master function
# Which is currently not working
makin.heatmaps <- function(x, nplots, param){
  plotlist <- list()
  cluster <- c(0.1, 0.2, 0.3)
  for(i in 1:nplots){
    dat <- filter(x, coff == cluster[i])
    plotlist[[i]] <- ggplot(dat, aes_string('deforest', 'dispersion', 
                                            fill = param)) + 
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
exp.heats <- makin.heatmaps(x = exp.final, nplots = 3, param = 'Expected.Value')

exp.heats[[1]]+ggtitle('A)')+
  exp.heats[[2]]+ggtitle('B)')+
  exp.heats[[3]]+ggtitle('C)')+
  plot_layout(design = layout, guides = 'collect')&
  scale_fill_viridis_c(limits = range(exp.final$Expected.Value),
                       name = "Expected Value")

# Max values
max.heats <- makin.heatmaps(x = max.final, nplots = 3, param = 'max')

max.heats[[1]]+ggtitle('A)')+
  max.heats[[2]]+ggtitle('B)')+
  max.heats[[3]]+ggtitle('C)')+
  plot_layout(design = layout, guides = 'collect')&
  scale_fill_viridis_c(limits = range(max.final$max),
                       name = "Max % Infection")

# Skew
skew.heats <- makin.heatmaps(x = skew.final, nplots = 3, param = 'skew')

skew.heats[[1]]+ggtitle('A)')+
  skew.heats[[2]]+ggtitle('B)')+
  skew.heats[[3]]+ggtitle('C)')+
  plot_layout(design = layout, guides = 'collect')&
  scale_fill_viridis_c(limits = range(skew.final$skew),
                       name = "Skew")

# Concentration
conc.heats <- makin.heatmaps(x = conc.final, nplots = 3, param = 'kappa')

conc.heats[[1]]+ggtitle('A)')+
  conc.heats[[2]]+ggtitle('B)')+
  conc.heats[[3]]+ggtitle('C)')+
  plot_layout(design = layout, guides = 'collect')&
  scale_fill_viridis_c(limits = range(conc.final$kappa),
                       name = "Kappa")

# Closer look at promising patterns ----------------------
# Expected values ##
exp.hi <- filter(exp.final, coff == '0.3')
exp.mid <- filter(exp.final, coff == '0.2')
exp.lo <- filter(exp.final, coff == '0.1')

# Expected value increases with deforestation at hi clustering
exp.hi.plot <- ggplot(data = exp.hi, aes(x = deforest, y = Expected.Value))+
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

# ggsave(exp.hi.plot, file = 'exphiplot.jpeg')

# No clear patterns at mid clustering
ggplot(data = exp.mid, aes(x = deforest, y = Expected.Value, color = dispersion))+
  geom_jitter(size = 2)+
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
  geom_jitter(size = 2)+
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

# Hi deforest = hi max infec at mid clustering
ggplot(data = max.mid, aes(x = deforest, y = max, color = dispersion))+
  geom_point(size = 2)+
  labs(x = "% Deforestation", y = "Maximum Infection")+
  scale_color_viridis_d(name = "Dispersion")+
  theme_bw(base_size = 18)+
  theme(panel.grid = element_blank())

cor.test(x = as.numeric(max.mid$deforest), y = max.mid$max, method = "spearman")
cor.test(x = as.numeric(max.mid$dispersion), y = max.mid$max, method = "spearman")

# No clear patterns at high clustering values
max.hi.plot <- ggplot(data = max.hi, aes(x = dispersion, y = max, 
                                         color = deforest))+
  geom_point(size = 2)+
  labs(x = "Dispersion", y = "Maximum Infection")+
  scale_color_viridis_d(name = "Deforestation")+
  theme_bw(base_size = 16)+
  theme(panel.grid = element_blank())

cor.test(x = as.numeric(max.hi$deforest), y = max.hi$max, method = "spearman")
cor.test(x = as.numeric(max.hi$dispersion), y = max.hi$max, method = "spearman")

# Skew ##
skew.hi <- filter(skew.final, coff == '0.3')
skew.mid <- filter(skew.final, coff == '0.2')
skew.lo <- filter(skew.final, coff == '0.1')

# Weird 'cutoff' point at high clustering
skew.hi.plot <- ggplot(data = skew.hi, aes(x = deforest, y = skew, 
                                           color = dispersion))+
  geom_point(size = 2)+
  labs(x = "% Deforestation", y = "Skew")+
  scale_color_viridis_d(name = "Dispersion")+
  theme_bw(base_size = 16)+
  theme(panel.grid = element_blank())

cor.test(x = as.numeric(skew.hi$deforest), y = skew.hi$skew, method = "spearman")
cor.test(x = as.numeric(skew.hi$dispersion), y = skew.hi$skew, method = "spearman")

# ggsave(skew.hi.plot, filename = 'skewhiplot.jpeg')

# No clear patterns at mid and low clustering
ggplot(data = skew.mid, aes(x = deforest, y = skew, color = dispersion))+
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

# Concentration param ##
conc.hi <- filter(conc.final, coff == '0.3')
conc.mid <- filter(conc.final, coff == '0.2')
conc.lo <- filter(conc.final, coff == '0.1')

# Highest dispersion value seems to have higher kappa
kappa.hi.plot <- ggplot(data = conc.hi,
                        aes(x = as.numeric(dispersion), y = kappa,
                            color = deforest))+
  geom_point(size = 2)+
  labs(x = "Dispersion", y = "Precision (Kappa)")+
  scale_color_viridis_d(name = "Deforestation")+
  theme_bw(base_size = 16)+
  theme(panel.grid = element_blank())

# ggsave(kappa.hi.plot, filename = "kappahiplot.jpeg")

cor.test(x = as.numeric(conc.hi$deforest), y = conc.hi$kappa, method = "spearman")
cor.test(x = as.numeric(conc.hi$dispersion), y = conc.hi$kappa, method = "spearman")

# Slight pattern at mid clustering
kappa.mid.plot <- ggplot(data = conc.mid, aes(x = as.numeric(dispersion), y = kappa, 
                            color = deforest))+
  geom_point(size = 2)+
  labs(x = "Dispersion", y = "Precision (Kappa)")+
  scale_color_viridis_d(name = "% Deforestation")+
  theme_bw(base_size = 18)+
  theme(panel.grid = element_blank())

ggsave(kappa.mid.plot, filename = "kappamidplot.jpeg")

cor.test(x = as.numeric(conc.mid$deforest), y = conc.mid$kappa, method = "spearman")
cor.test(x = as.numeric(conc.mid$dispersion), y = conc.mid$kappa, 
         method = "spearman")

# No clear patterns at low clustering
ggplot(data = conc.lo, aes(x = as.numeric(deforest), y = kappa, color = dispersion))+
  geom_line(size = 1.5)+
  labs(x = "% Deforestation", y = "Kappa")+
  scale_color_viridis_d(name = "Dispersion")+
  theme_bw(base_size = 1)+
  theme(panel.grid = element_blank())

cor.test(x = as.numeric(conc.lo$deforest), y = conc.lo$kappa, method = "spearman")
cor.test(x = as.numeric(conc.lo$dispersion), y = conc.lo$kappa, method = "spearman")

# Put plots at high clustering together
megadot <- exp.hi.plot+ggtitle("A)")+
  max.hi.plot+ggtitle("B)")+
  skew.hi.plot+ggtitle("C)")+
  kappa.hi.plot+ggtitle("D)")

ggsave(megadot, file = 'megadot.jpeg', width = 10, height = 7.5,
       units = 'in')

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

# Well this is weird
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
Moran(loc.raster$val)
# Close to 0
