import numpy
import random
import matplotlib.pyplot as plt

##############################################################################
################        Neighbors Function             #######################
##############################################################################

def neighbors_base(mat, row, col, radius=1):
    # neighbors: finds nearest neighbors in a matrix
    # inputs: mat = matrix to look at, row and col = indexes, radius=number of cells around
    # output: a list of the contents of cells around the index (out of bounds returned as 0)
    rows, cols = len(mat), len(mat[0])
    out = [] #out is a list of the neighbors
    for i in range(row - radius - 1, row + radius):
        row = []
        for j in range(col - radius - 1, col + radius):
            if 0 <= i < rows and 0 <= j < cols:
                row.append(mat[i][j])
            else:
                row.append(None)
        out.append(row)

    # make into a flat list
    flat_list = []
    for sublist in out:
        for item in sublist:
            flat_list.append(item)

    return flat_list

##############################################################################
################          Construct Landscape          #######################
##############################################################################

#Specify matrix size, patch size, and number of patches
matrix_size = 50
n_patches = 2
n_draws = 40

def MakeLandscape(size, patches, draws):
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

    #Use neighbors function to remove solitary points
    neighbor_array = numpy.empty(shape=(1,11), dtype = "float")
    neighbor_array[:] = numpy.nan
    for i in range(0, matrix_size):
        for j in range(0, matrix_size):
            if coffee[i,j] == 0:
                neighbor_out = neighbors_base(mat=coffee, row=i, col=j, radius=1)
                neighbor_out.append(i)
                neighbor_out.append(j)
                neighbor_array = numpy.vstack([neighbor_array, neighbor_out])


    for row in range(1, numpy.size(neighbor_array, 0)):
        if numpy.any(neighbor_array[row,0:8] == 0):
            coffee[int(neighbor_array[row,9]),int(neighbor_array[row,10])] = 0
        else:
            coffee[int(neighbor_array[row,9]),int(neighbor_array[row,10])] = None

    indices = numpy.argwhere(coffee)
    landscape = numpy.empty(shape=(matrix_size,matrix_size))
    landscape[coffee == 0] = numpy.nan

    for i in range(0, len(indices)):
        landscape[indices[i][0], indices[i][1]] = numpy.random.uniform(0.65, 0.95, 1)

    # "Deforestation" bit would go here. Would be similar to how coffee patches were drawn

    return coffee, landscape

(coffee, landscape) = MakeLandscape(size=matrix_size, patches=n_patches, draws=n_draws)

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