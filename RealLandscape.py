import rasterio
import numpy
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
    if numpy.isnan(land1[i,:]).all():
        rows.append(i)

land1 = numpy.delete(land1, rows, axis=0)
