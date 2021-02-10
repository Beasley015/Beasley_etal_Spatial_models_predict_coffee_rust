import rasterio
import numpy
from matplotlib import pyplot

# Read full raster file
filepath = "coffeeRust_landuses/w001001.adf"
full = rasterio.open(filepath)

# Extract array
land = full.read(1)

# Add null values
land = land.astype('float64')
land[numpy.where(land == 255)] = numpy.nan

# Split the two sections



