#------------------------------------------------------------------------------
# ABOUT Fig1.py
#------------------------------------------------------------------------------

# This script is intended to demonstrate the use of the NLMpy Python package,
# by creating Fig. 1 in the paper.

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
# CREATE NEUTRAL LANDSCAPE MODELS
#------------------------------------------------------------------------------

import nlmpy
import numpy as np

np.random.seed(0) # So that the same NLMs are produced each time
nRow = 50 # Number of rows
nCol = 50 # Number of columns

# Random NLM
Fig1a = nlmpy.random(nRow, nCol)
# Planar gradient NLM
Fig1b = nlmpy.planarGradient(nRow, nCol)
# Edge gradient NLM
Fig1c = nlmpy.edgeGradient(nRow, nCol)
# Mask example
Fig1d = np.zeros((nRow, nCol))
Fig1d[10:25, 10:25] = 1
# Distance gradient NLM
Fig1e = nlmpy.distanceGradient(Fig1d)
# Midpoint displacement NLM
Fig1f = nlmpy.mpd(nRow, nCol, 0.75)
# Random rectangular cluster NLM
Fig1g = nlmpy.randomRectangularCluster(nRow, nCol, 4, 8)
# Random element nearest-neighbour NLM
Fig1h = nlmpy.randomElementNN(nRow, nCol, 200)
# Random cluster nearest-neighbour NLM
Fig1i = nlmpy.randomClusterNN(nRow, nCol, 0.4)
# Blended NLM
Fig1j = nlmpy.blendArray(Fig1f, [Fig1c])
# Patch blended NLM
Fig1k = nlmpy.blendClusterArray(Fig1h, [Fig1e], [1.5])
# Classified random cluster nearest-neighbour NLM
Fig1l = nlmpy.classifyArray(Fig1i, [1,1,1,1])
# Percolation NLM
Fig1m = nlmpy.classifyArray(Fig1a, [1 - 0.5, 0.5])
# Binary random rectangular cluster NLM
Fig1n = nlmpy.classifyArray(Fig1g, [1 - 0.75, 0.75])
# Classified midpoint displacement NLM
Fig1o = nlmpy.classifyArray(Fig1f, [1,1,1])
# Classified midpoint displacement NLM, with limited classification extent
Fig1p = nlmpy.classifyArray(Fig1f, [1,1,1], classifyMask=Fig1d)
# Masked planar gradient NLM
Fig1q = nlmpy.planarGradient(nRow, nCol, 90, mask=Fig1n)
# Hierarchical NLM
Fig1r = np.where(Fig1o==2, Fig1m + 2, Fig1o)
# Rotated NLM
Fig1s = np.rot90(Fig1l)
# Transposed NLM
Fig1t = np.transpose(Fig1o)

#------------------------------------------------------------------------------
# PLOT NEUTRAL LANDSCAPE MODELS
#------------------------------------------------------------------------------

import matplotlib.pyplot as plt
import matplotlib as mpl

# Make continuous colour map
cmapCont = mpl.colors.LinearSegmentedColormap.from_list('heat', ['#000000', '#bd0026', '#bd0026', '#fd8d3c', '#fecc5c'])

# Make classified colour map
cmapClas = mpl.colors.ListedColormap(['#8dd3c7', '#ffffb3', '#bebada', '#fb8072'])
bounds=[0,1,2,3,4]
norm = mpl.colors.BoundaryNorm(bounds, cmapClas.N)

# Create figure
fig = plt.figure(1, figsize=(145/25.4, 140/25.4))
mpl.rc('axes', linewidth=0.5) # set all axes line widths

# Plot continuous NLMs
contNLMs = [Fig1a, Fig1b, Fig1c, Fig1e, Fig1f, Fig1g, Fig1h, Fig1i, Fig1j, Fig1k, Fig1q]
labels = ['(a)', '(b)', '(c)', '(e)', '(f)', '(g)', '(h)', '(i)', '(j)', '(k)', '(q)']
subplot = [(0,0),(0,1),(0,2),(0,4),(2,0),(2,1),(2,2),(2,3),(2,4),(4,0),(6,1)]
for n in range(len(contNLMs)):
    plt.subplot2grid((9,5), subplot[n], rowspan=2)
    plt.xticks(np.arange(0))
    plt.yticks(np.arange(0))
    plt.imshow(contNLMs[n], interpolation='none', cmap=cmapCont)
    plt.title(labels[n], fontsize=8)
# Plot classified NLMs
qualNLMs = [Fig1d, Fig1l, Fig1m, Fig1n, Fig1o, Fig1p, Fig1r, Fig1s, Fig1t]
labels = ['(d)', '(l)', '(m)', '(n)', '(o)', '(p)', '(r)', '(s)', '(t)']
subplot = [(0,3),(4,1),(4,2),(4,3),(4,4),(6,0),(6,2),(6,3),(6,4)]
for n in range(len(qualNLMs)):
    plt.subplot2grid((9,5), subplot[n], rowspan=2)
    plt.xticks(np.arange(0))
    plt.yticks(np.arange(0))
    plt.imshow(qualNLMs[n], interpolation='none', cmap=cmapClas, norm=norm)
    plt.title(labels[n], fontsize=8)

# Plot continuous legend
plt.subplot2grid((9,14), (8,1), colspan=5)
x = np.array([range(21)]).astype('float') / 21
plt.imshow(x, aspect=1, cmap=cmapCont)
plt.yticks(np.arange(0))
plt.xticks(np.arange(0, 24, 4), [0.0, 0.2, 0.4, 0.6, 0.8, 1.0], fontsize=8)
plt.tick_params(direction='out', length=3, width=0.5, top='off')
plt.title("Continuous value", fontsize=8)
# Plot classified legend
plt.subplot2grid((9,14), (8,8), colspan=5)
x = np.array([np.array(np.repeat(range(4), 5))])
plt.imshow(x, interpolation='none', aspect=1, cmap=cmapClas, norm=norm)
plt.yticks(np.arange(0))
plt.xticks([2,7,12,17], [0,1,2,3], fontsize=8)
plt.tick_params(direction='out', length=3, width=0.5, top='off')
plt.title("Classed value", fontsize=8)

plt.show()

#------------------------------------------------------------------------------
