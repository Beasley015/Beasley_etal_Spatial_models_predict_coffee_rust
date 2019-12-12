#------------------------------------------------------------------------------
# ABOUT Fig2b.py
#------------------------------------------------------------------------------

# This script is intended to demonstrate the use of the NLMpy Python package,
# by creating Fig. 2b in the paper.

#------------------------------------------------------------------------------
# LICENSING
#------------------------------------------------------------------------------

# The MIT License (MIT)

# Copyright (c) 2014 Thomas R. Etherington, E. Penelope Holland, and
# David O'Sullivan.

# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:

# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

#------------------------------------------------------------------------------
# CREATE GIS BASED HIERARCHICAL NEUTRAL LANDSCAPE MODEL
#------------------------------------------------------------------------------

import nlmpy
import numpy as np

# Read in metadata
inFile = open('DEM.asc', 'r')
textLine = inFile.readline()
nCol = int(textLine.split()[1])
textLine = inFile.readline()
nRow = int(textLine.split()[1])
textLine = inFile.readline()
xll = int(textLine.split()[1])
textLine = inFile.readline()
yll = int(textLine.split()[1])
textLine = inFile.readline()
cellSize = int(textLine.split()[1])
inFile.close()
# Import the digital elevation model (DEM) ASCII file
demArray = np.loadtxt('DEM.asc', skiprows=6)

# Create the NLMs for the different hierarchical levels
np.random.seed(0) # So that the same NLMs are produced each time
# Represent the agricultural land in the valley bottom using a random element
# nearest-neighbour NLM, which is then classified
nlm1 = nlmpy.randomElementNN(nRow, nCol, 5000)
nlm1 = nlmpy.classifyArray(nlm1, [1,1,1])
# Represent the pastoral and forestry land on the mountain slopes using a
# random cluster nearest-neighbour NLM, which is then classified
nlm2 = nlmpy.randomClusterNN(nRow, nCol, 0.58)
nlm2 = nlmpy.classifyArray(nlm2, [1,2])
# Create a NlM that is a combination of these two NLMs, using an elevation of
# 380 as a threshold between the two NLMs
nlm = np.where(demArray<=380, nlm1, nlm2 + 3)
# Replace all values with an elevation of 348 with a no data value to ignore
# the lake.
np.place(nlm, demArray == 348, np.nan)
# Replace any values above an elevation of 1100 with a sixth class representing
# the rough grassland on the mountain tops
np.place(nlm, demArray >= 1100, 5)
# Export the NLM as an ASCII raster grid
nlmpy.exportASCIIGrid("fig2bNLM.asc", nlm, xll, yll, cellSize)

#------------------------------------------------------------------------------
# PLOT HIERARCHICAL NLM
#------------------------------------------------------------------------------

import matplotlib.pyplot as plt
import matplotlib as mpl

# Create figure
fig = plt.figure(1, figsize=(72.5/25.4, 90/25.4))
mpl.rc('axes', linewidth=0.5) # set all axes line widths

# Make classified colour map
cmapClas = mpl.colors.ListedColormap(['#8dd3c7', '#ffffb3', '#bebada', '#fb8072', '#80b1d3', '#fdb462'])
bounds=[0,1,2,3,4,5,6]
norm = mpl.colors.BoundaryNorm(bounds, cmapClas.N)

# Plot hierarchical NLM
plt.subplot2grid((5,1),(0,0), rowspan=4)
plt.xticks(np.arange(0))
plt.yticks(np.arange(0))
np.place(nlm, nlm == -9999, np.nan)
plt.imshow(nlm, interpolation='none', cmap=cmapClas, norm=norm)
plt.title("(b)", fontsize=8)

# Plot classified legend
plt.subplot2grid((5,1), (4,0), colspan=1)
x = np.array([np.array(np.repeat(range(6), 5))])
plt.imshow(x, interpolation='none', aspect=1, cmap=cmapClas, norm=norm)
plt.yticks(np.arange(0))
plt.xticks([2,7,12,17,22,27], [0,1,2,3,4,5], fontsize=8)
plt.tick_params(direction='out', length=3, width=0.5, top='off')
plt.title("Classed value", fontsize=8)

plt.show()

#------------------------------------------------------------------------------
