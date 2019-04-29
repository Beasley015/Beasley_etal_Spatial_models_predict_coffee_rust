#####################################################
# Modeling coffee rust movement through a landscape #
# N. Aristizabal, E. Beasley, & E. Bueno            #
# Spring 2019                                       #
#####################################################

# Read in model outputs --------------------------------------
# Read all csv's into a list
filenames <- list.files("Outputs", pattern = "*.csv", full.names = T)
shortnames <- list.files("Outputs", pattern = "*.csv")
output.list <- lapply(filenames, read.csv, header = F)

# Rename columns of each dataframe
newnames <- c("Time", "PercInf")
output.list <- lapply(output.list, setNames, newnames)

# Set up data -----------------------------------------------
# Add column to denote replicates
replicate <- logical()
for(i in 1:50){
  new <- rep(i, 1000)
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

for(i in loop.ready) {
  def[[i]] <- mid(shortnames[[i]],4,2)
  disp[[i]] <- mid(shortnames[[i]],10,1)
}

for(i in 1:length(output.list)){
  deforest <- rep(def[[i]], nrow(output.list[[1]]))
  dispersion <- rep(disp[[i]], nrow(output.list[[1]]))
  output.list[[i]] <- cbind(output.list[[i]], deforest, dispersion)
}
head(output.list[[1]])

# Turn list into big-ass data frame
output.mat <- do.call(rbind, output.list)

# plot variation in Percentage Infestation among replicates
ggplot(output.mat, aes(x = deforest, y= PercInf)) +
  geom_boxplot(aes(fill=factor(deforest)))

# individual percent infestation ~ deforestation
ggplot(output.mat, aes(x = dispersion, y= PercInf)) +
  geom_boxplot(aes(fill=factor(dispersion)))

# percent infestation through time steps
ggplot(output.mat, aes(x = Time, y = PercInf, group = as.factor(replicate), alpha = 0.3)) + 
  geom_line() +
  facet_grid(vars(deforest), vars(dispersion)) +
  theme_classic()
  legend
  
# heat map
ggplot(output.mat, aes(deforest, dispersion, fill = PercInf)) + geom_raster(hjust = 0, vjust = 0)
  
