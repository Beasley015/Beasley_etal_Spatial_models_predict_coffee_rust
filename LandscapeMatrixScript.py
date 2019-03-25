import numpy
import random
import matplotlib.pyplot as plt

##############################################################################
################          Construct Landscape          #######################
##############################################################################

#Specify matrix size, patch size, and number of patches
matrix_size = 50
n_patches = 2
patch_size = 10

#Create blank landscape
landscape = numpy.zeros(shape = (matrix_size, matrix_size), dtype = dict)

#Create dictionary for coffee
coffee_clean = {'coffee':0}

#Seed two random locations in the landscape with uninfected coffee
randoms = [None]*n_patches
for cells in range(0, n_patches):
    randoms[cells] = (random.randint(1, matrix_size-1), random.randint(1, matrix_size-1))

for coords in randoms:
    landscape[coords] = coffee_clean

#"Grow" coffee patches to desired size by filling in spiral
N, S, W, E = (0, -1), (0, 1), (-1, 0), (1, 0) # directions
turn_right = {N: E, E: S, S: W, W: N} # old -> new direction

def spiral(npatches, ncells):
    for patch in range(0,npatches):
        for step in range(0,ncells):
            x, y = randoms[patch] # start at the patch seed
            dx, dy = N # initial direction
            count = coffee_clean
            while True:
                landscape[x,y] = coffee_clean #change cell value
                # try to turn right
                new_dx, new_dy = turn_right[dx,dy]
                new_x, new_y = x + new_dx, y + new_dy
                if landscape[new_x, new_y] is 0: #can turn right
                    x, y = new_x, new_y
                    dx, dy = new_dx, new_dy
                else: # try to move straight
                    x, y = x + dx, y + dy
                    if not (0 <= x < matrix_size and 0 <= y < matrix_size):
                        return(landscape)


patchy_landscape = spiral(npatches = n_patches, ncells = patch_size)

plt.imshow(patchy_landscape)

#Change zeroes to dictionaries with key "matrix" and blank values

#Fill matrix values with friction values