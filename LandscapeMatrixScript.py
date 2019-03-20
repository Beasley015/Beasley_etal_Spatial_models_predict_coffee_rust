import numpy

##############################################################################
################          Construct Landscape          #######################
##############################################################################

#Specify matrix size, patch size, and number of patches
matrix_size = 50
n_patches = 2
patch_size = 10

#Create blank landscape
landscape = numpy.zeros(shape = (matrix_size, matrix_size))

#Create dictionary for coffee
coffee_clean = {'coffee':0}

#Seed two random locations in the landscape with uninfected coffee

#"Grow" coffee patches to desired size by changing neighbors of coffee cells

#Change zeroes to dictionaries with key "matrix" and blank values

#Fill matrix values with friction values