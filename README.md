# Spatially explicit models predict coffee rust spread in fragmented landscapes
### Beasley, E.M., Aristiz&aacute;bal, N., Bueno, E., & White, E.R. 2022. 

For questions about data or code, please contact the corresponding author: Emily Beasley (ebeasley{at}uvm{dot}edu).

## Abstract
*Context*. Landscape structure influences the spread of plant pathogens, including coffee leaf rust, a fungal disease affecting the coffee industry. Rust transmission is likely affected by landscape structure through the dispersal of wind-borne spores. Previous studies found positive associations between rust incidence and the proportion of pasture cover, suggesting deforestation may facilitate spore dispersal.

*Objectives*. We explored the links between landscape structure and coffee rust by modeling the spread of rust transmission. We investigated how 1) spatial clustering of coffee farms, 2) proportion of landscape deforestation, and 3) clustering of deforestation affects the speed of rust transmission.

*Methods*. We developed a probabilistic model to simulate within-patch and between-patch transmission in simulated and real landscapes. We modeled within-patch transmission using a probabilistic cellular automata model and between-patch transmission using a random walk with spore movement inhibited by canopy cover. 

*Results*. Clustering of coffee farms is the primary driver of rust transmission. Deforestation is a secondary driver of rust spread: outbreaks spread more rapidly in landscapes where deforested areas are evenly dispersed throughout the landscape. When applied to real landscapes in Costa Rica, the model yields the same trends as simulated landscapes and suggests increased amounts of coffee near the starting location of the outbreak are correlated with more rapid rust spread.

*Conclusions*. It is essential to consider landscape structure when managing the spread of crop diseases. Increasing the spacing between coffee farms and reducing forest fragmentation in coffee-growing regions can benefit biodiversity conservation and reduce the economic impacts of coffee rust.

## Code
LandscapeSim.py - generates simulated landscapes and simulates rust spread
RealLandscape.py - simulates rust spread on rasters representing real landscapes

LandscapeSimResults.R - creates figures from results of simulated landscapes
RealLandscapeResults.R - creates figures from results of real landscapes

## Data and outputs
folder coffeeRust_landuses - contains raster data for real landscapes
folder Outputs_buffered - model outputs for simulated landscapes
land1raw.csv, land2raw.csv, land1.csv, land2.csv - model outputs for real landscapes
