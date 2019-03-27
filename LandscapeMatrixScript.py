import numpy
import random
import matplotlib.pyplot as plt

##############################################################################
################          Construct Landscape          #######################
##############################################################################

#Specify matrix size, patch size, and number of patches
matrix_size = 50
n_patches = 2
n_draws = 20

#Create blank landscape
coffee = numpy.empty(shape = (matrix_size, matrix_size))
coffee[:] = numpy.nan

#Seed two random locations in the landscape with coffee
randoms = [None]*n_patches
for cells in range(0, n_patches):
    randoms[cells] = (random.randint(1, matrix_size-1), random.randint(1, matrix_size-1))

for coords in randoms:
    coffee[coords] = 0

#Draw from beta dists centered around coords
betavals = numpy.empty((len(randoms), n_draws, 2))

for patch in range(0, n_patches):
    coords = numpy.array(randoms[patch])
    mu = coords/50
    stdev = 0.003
    alpha = ((1-mu)/stdev)-(1/mu)*numpy.square(mu)
    beta = alpha * (1 / mu - 1)
    x = numpy.random.beta(alpha[0], beta[0], n_draws)
    y = numpy.random.beta(alpha[1], beta[1], n_draws)
    betavals[patch,:,0] = x
    betavals[patch,:,1] = y

#Convert beta values to new coords
betavals = betavals*50
betavals = betavals.round()

#Grow patches from coords
for patch in range(0, n_patches):
    coords = betavals[patch,:,:]
    for cell in range(0, len(coords)):
        i,j = coords[cell,:]
        coffee[int(i), int(j)] = 0

plt.matshow(coffee)

#Use neighbors function to remove solitary points







###################################################################
##############    Cellular Automata    ############################
###################################################################

###################################################################
##############    Propagule Release    ############################
###################################################################

###################################################################
#################    Random Walk    ###############################
###################################################################

###################################################################
##############    Propagule infections    #########################
###################################################################