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
    for i in range(row - radius, row + radius + 1):
        row = []
        for j in range(col - radius, col + radius + 1):
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

    del flat_list[4]
    return flat_list

##############################################################################
################          Construct Landscape          #######################
##############################################################################

#Specify matrix size, patch size, and number of patches
matrix_size = 50
n_patches = 2
n_draws = 50

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
        mu = coords/49
        stdev = 0.003
        alpha = ((1-mu)/stdev)-(1/mu)*numpy.square(mu)
        beta = alpha * (1 / mu - 1)
        x = numpy.random.beta(alpha[0], beta[0], n_draws)
        y = numpy.random.beta(alpha[1], beta[1], n_draws)
        betavals[patch,:,0] = x
        betavals[patch,:,1] = y

    #Convert beta values to new coords
    betavals = betavals*49
    betavals = betavals.round()

    #Grow patches from coords
    for patch in range(0, n_patches):
        coords = betavals[patch,:,:]
        for cell in range(0, len(coords)):
            i,j = coords[cell,:]
            coffee[int(i), int(j)] = 0

    #Use neighbors function to remove solitary points
    neighbor_array = numpy.empty(shape=(1,10), dtype = "float")
    neighbor_array[:] = numpy.nan
    for i in range(matrix_size):
        for j in range(matrix_size):
            if coffee[i,j] == 0:
                neighbor_out = neighbors_base(mat=coffee, row=i, col=j, radius=1)
                neighbor_out.append(i)
                neighbor_out.append(j)
                neighbor_array = numpy.vstack([neighbor_array, neighbor_out])

    neighbor_array = neighbor_array[1:,:]

    coffeenew = coffee
    for row in range(0, numpy.size(neighbor_array, 0)):
        if numpy.any(neighbor_array[row,0:8] == 0) == False:
            coffeenew[int(neighbor_array[row,8]),int(neighbor_array[row,9])] = numpy.nan

    indices = numpy.argwhere(coffee)
    landscape = numpy.empty(shape=(matrix_size,matrix_size))
    landscape[coffee == 0] = numpy.nan

    for i in range(0, len(indices)):
        landscape[indices[i][0], indices[i][1]] = numpy.random.uniform(0.65, 0.95, 1)

    # "Deforestation" bit would go here. Would be similar to how coffee patches were drawn

    return coffee, landscape

(coffee, landscape) = MakeLandscape(size=matrix_size, patches=n_patches, draws=n_draws)
plt.matshow(coffee)
plt.matshow(landscape, vmin=0, vmax=1)

#Pick one coffee cell to initialize infection
coffee_zeros = numpy.where(coffee == 0)
randrow = random.randint(0, numpy.size(coffee_zeros,1))

coffee[coffee_zeros[0][randrow], coffee_zeros[1][randrow]] = 1

plt.matshow(coffee)

###################################################################
##############    Cellular Automata    ############################
###################################################################
def cellaut():
    # Pull index of all coffee cells with value 0
    coffee_zeros = numpy.where(coffee == 0)

    # Get neighborhood of that cell
    neighbors_clean = numpy.empty(shape=(1,8))
    neighbors_clean[:] = numpy.nan

    for i in range(0, numpy.size(coffee_zeros, 1)):
        clean_row = neighbors_base(mat=coffee, row=coffee_zeros[0][i], col=coffee_zeros[1][i], radius=1)
        neighbors_clean = numpy.vstack([neighbors_clean, clean_row])

    neighbors_clean = neighbors_clean[1:,:]
    neighbors_clean = neighbors_clean.astype(numpy.float64)

    # Sum the values of the neighborhood to get number of infected neighbors
    rowsums = []
    for i in range(0, numpy.size(neighbors_clean, 0)):
        row = neighbors_clean[i,:]
        row = row[~numpy.isnan(row)]
        add_row = row.sum()
        rowsums.append(add_row)

    # Use number of infected neighbors to get success probability in a bernoulli trial
    bern_out = []
    for i in range(0,len(rowsums)):
        out = numpy.random.binomial(n = 1, p = (0.1*rowsums[i]))
        bern_out.append(out)

    # If bern trial is a success, change focal cell to infected
    for i in range(0, len(coffee_zeros[1])):
        coffee[coffee_zeros[0][i], coffee_zeros[1][i]] = bern_out[i]

    return(coffee)

coffee = cellaut()
plt.matshow(coffee)

###################################################################
##############    Propagule Release    ############################
###################################################################

#Create list of tuples for changing coords
coord_change = [(-1,-1), (-1,0), (-1, 1),
                (0,-1), (0, 1),
                (1, -1), (1,0),(1,1)]

def new_spore():
    #Get neighbors for each infected cell
    coffee_inf = numpy.where(coffee == 1)

    #Get neighbors of infected cells
    land_neighbors = numpy.empty(shape = (1, 8))
    land_neighbors[:] = numpy.nan

    for i in range(0, numpy.size(coffee_inf, 1)):
        land_row = neighbors_base(mat=coffee, row=coffee_inf[0][i], col=coffee_inf[1][i], radius=1)
        land_neighbors = numpy.vstack([land_neighbors, land_row])

    land_neighbors = land_neighbors[1:,:]

    #Get col numbers of landscape cells
    land_pos = []
    for i in range(0, numpy.size(land_neighbors,0)):
        row = land_neighbors[i,:]
        if numpy.isnan(row).any: #There's an issue here
            nans = numpy.isnan(row)
            pos = numpy.where(nans == True)
            pos = pos[0].tolist()
            land_pos.append(pos)
        else:
            land_pos.append([None])

    #Create sparse matrix of propagules
    spores = {}

    for i in range(0, len(land_pos)):
        place = random.choice(land_pos[i])
        release = coord_change[place]
        new_coord = (coffee_inf[0][i]+release[0], coffee_inf[1][i]+release[1])
        spores.update({new_coord:1})

    return spores

walkers = new_spore()

###################################################################
#################    Random Walk    ###############################
###################################################################

###################################################################
##############    Propagule infections    #########################
###################################################################