import rasterio
import numpy
import random
from matplotlib import pyplot as plt

# Read in land cover data #############
# Read full raster file
filepath = "coffeeRust_landuses/w001001.adf"
full = rasterio.open(filepath)

# Extract array
land = full.read(1)

# Clean rasters #############
# Add null values
land = land.astype('float64')
land[numpy.where(land == 255)] = numpy.nan

# Find cols with all nan to split raster
cols = []
for i in range(land.shape[1]):
    if numpy.isnan(land[:, i]).all():
        cols.append(i)

# Split the raster sections
land1 = land[:, 1:cols[0] - 1]
land2 = land[:, cols[2] + 1:]

# Trim rows with all nan from first raster section
rows = []
for i in range(land1.shape[0]):
    if numpy.isnan(land1[i, :]).all():
        rows.append(i)

land1 = numpy.delete(land1, rows, axis=0)

###################################################
#    Function for creating coffee arrays          #
###################################################
def MakeCoffee(raster_in):
    # Create array of blanks
    coffee = numpy.empty(shape=raster_in.shape)
    # Fill with values representing uninfected coffee
    coffee[numpy.where(raster_in == 3)] = 0
    coffee[numpy.where(raster_in != 3)] = numpy.nan

    # Initialize infection
    coffee_zeros = numpy.where(coffee == 0)
    randrow = random.randint(0, numpy.size(coffee_zeros, 1) - 1)

    coffee[coffee_zeros[0][randrow], coffee_zeros[1][randrow]] = 1

    start = numpy.where(coffee == 1)

    return(coffee, start)

###################################################
#        Creating landscape arrays                #
###################################################

# Create landscape array
land1[numpy.where(land1 == 3)] = numpy.nan