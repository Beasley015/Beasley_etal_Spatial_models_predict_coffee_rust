import numpy
import random
import matplotlib.pyplot as plt
import scipy.stats as stats

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

def MakeLandscape(size, patches, draws, deforest, disp, ddraws):
    #Create blank landscape
    coffee = numpy.empty(shape = (size, size))
    coffee[:] = numpy.nan

    #Seed two random locations in the landscape with coffee
    randoms = [None]*patches
    for cells in range(0, patches):
        randoms[cells] = (random.randint(1, size-1), random.randint(1, size-1))

    for coords in randoms:
        coffee[coords] = 0

    #Draw from beta dists centered around coords
    betavals = numpy.empty((draws, 2, patches))

    for patch in range(0, patches):
        coords = numpy.array(randoms[patch])
        mu = coords
        stdev = 1.5
        a1, b1 = (0-mu[0])/stdev, ((size-1)-mu[0])/stdev
        a2, b2 = (0 - mu[1]) / stdev, ((size-1) - mu[1]) / stdev
        x = stats.truncnorm.rvs(a1, b1, size=draws, loc = mu[0], scale = stdev)
        y = stats.truncnorm.rvs(a2, b2, size=draws, loc = mu[1], scale = stdev)
        betavals[:,0,patch] = x
        betavals[:,1,patch] = y

    betavals = betavals.round()

    #Grow patches from coords
    for patch in range(0, patches):
        coords = betavals[:, :, patch]
        for cell in range(0, len(coords)):
            i,j = coords[cell,:]
            coffee[int(i), int(j)] = 0

    #Use neighbors function to remove solitary points
    neighbor_array = numpy.empty(shape=(1,10), dtype = "float")
    neighbor_array[:] = numpy.nan
    for i in range(size):
        for j in range(size):
            if coffee[i,j] == 0:
                neighbor_out = neighbors_base(mat=coffee, row=i, col=j, radius=1)
                neighbor_out.append(i)
                neighbor_out.append(j)
                neighbor_array = numpy.vstack([neighbor_array, neighbor_out])

    neighbor_array = neighbor_array[1:,:]

    for row in range(0, numpy.size(neighbor_array, 0)):
        if numpy.any(neighbor_array[row,0:8] == 0) == False:
            coffee[int(neighbor_array[row,8]),int(neighbor_array[row,9])] = numpy.nan

    #Create fully forested landscape
    indices = numpy.argwhere(coffee)
    landscape = numpy.empty(shape=(size,size))
    landscape[coffee == 0] = numpy.nan

    for i in range(0, len(indices)):
        landscape[indices[i][0], indices[i][1]] = numpy.random.uniform(0.3, 0.95, 1)

    # Simulate deforestation
    while (len(numpy.where(landscape < 0.3)[0])/numpy.count_nonzero(~numpy.isnan(landscape)) < deforest):
        def_seed = random.choice(indices)
        mu = def_seed
        stdev = disp
        a1, b1 = (0 - mu[0]) / stdev, ((size-1) - mu[0]) / stdev
        a2, b2 = (0 - mu[1]) / stdev, ((size-1) - mu[1]) / stdev
        x = stats.truncnorm.rvs(a1, b1, size=ddraws, loc=mu[0], scale=stdev)
        y = stats.truncnorm.rvs(a2, b2, size=ddraws, loc=mu[1], scale=stdev)

        clear_coords = numpy.column_stack((x,y))
        clear_coords = clear_coords.round().astype("int")

        for i in range(0, len(clear_coords[:,0])):
            if coffee[clear_coords[i,0],clear_coords[i,1]] != 0:
                landscape[clear_coords[i,0],clear_coords[i,1]] = numpy.random.uniform(0.05, 0.3, 1)

    # Pick one coffee cell to initialize infection
    coffee_zeros = numpy.where(coffee == 0)
    randrow = random.randint(0, numpy.size(coffee_zeros, 1))

    coffee[coffee_zeros[0][randrow], coffee_zeros[1][randrow]] = 1

    return (coffee, landscape)

###################################################################
##############    Cellular Automata    ############################
###################################################################

def cellaut(mat, land):
    #Create new array to store next time step
    coffee_new = mat

    # Pull index of all coffee cells with value 0
    coffee_zeros = numpy.where(mat == 0)

    # Get neighborhood of that cell
    neighbors_clean = numpy.empty(shape=(1,8))
    neighbors_clean[:] = numpy.nan

    for i in range(0, numpy.size(coffee_zeros, 1)):
        clean_row = neighbors_base(mat=mat, row=coffee_zeros[0][i], col=coffee_zeros[1][i], radius=1)
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
        coffee_new[coffee_zeros[0][i], coffee_zeros[1][i]] = bern_out[i]

    #Replace old array with new one
    mat = coffee_new

    return(mat)

###################################################################
##############    Propagule Release    ############################
###################################################################

def new_spore(mat, coord):
    #Get coords for each infected cell
    coffee_inf = numpy.where(mat == 1)

    #Get neighbors of infected cells
    land_neighbors = numpy.empty(shape = (1, 8))
    land_neighbors[:] = numpy.nan

    for i in range(0, numpy.size(coffee_inf, 1)):
        land_row = neighbors_base(mat=mat, row=coffee_inf[0][i], col=coffee_inf[1][i], radius=1)
        land_neighbors = numpy.vstack([land_neighbors, land_row])

    land_neighbors = land_neighbors[1:,:]

    #Get col numbers of landscape cells
    land_pos = []
    for i in range(0, numpy.size(land_neighbors,0)):
        row = land_neighbors[i,:]
        row = numpy.array(row, dtype = numpy.float64)
        if numpy.isnan(row).any():
            nans = numpy.isnan(row)
            pos = numpy.where(nans == True)
            pos = pos[0].tolist()
            land_pos.append(pos)
        else:
            land_pos.append(None)

    #Create sparse matrix of propagules
    spores = []

    for i in range(0, len(land_pos)):
        if land_pos[i] != None:
            place = random.choice(land_pos[i])
            release = coord[place]
            new_coord = (coffee_inf[0][i]+release[0], coffee_inf[1][i]+release[1])
            if new_coord[0] < matrix_size & new_coord[1] < matrix_size:
                spores.update(new_coord)

    return spores

###################################################################
#################    Random Walk    ###############################
###################################################################

#Wind effects: new array that affects probability of movement, but not
#cost? Or both prob and cost?

#Humidity effects: Change infection prob based on neighborhood of
#target coffee cell

#Clean up this function
#Change walkers object into list of tuples
#so multiple spores can occupy a cell

def spore_walk(spores, land, mat, coord):
    for i in range(0, len(spores)):
        step_credit = 15
        old_coords = list(spores)[i]

        while step_credit > 0:
            land_neighbors = numpy.array(neighbors_base(mat=land, row=old_coords[0], col=old_coords[1]),
                                         dtype=numpy.float64)

            if numpy.isnan(land_neighbors).any():
                land_neighbors[numpy.isnan(land_neighbors)] = 0
                land_neighbors[land_neighbors == None] = 0

            land_neighbors[land_neighbors > 0] = 1
            movement = random.choices(population = coord, weights=land_neighbors, k = 1)
            new_coords = (old_coords[0] + movement[0][0], old_coords[1] + movement[0][1])

            spores.update({old_coords: None})
            spores.update({new_coords: 1})

            old_coords = new_coords

            if new_coords[0] > 99 or new_coords[1] > 99:
                break

            else:
                step_credit = step_credit-land[new_coords[0], new_coords[1]]

                new_neighbors = neighbors_base(mat=mat, row=new_coords[0], col=new_coords[1])
                new_neighbors = numpy.array(new_neighbors)

                if numpy.any(new_neighbors == 0):
                    infec_newcoord = random.choice([coord[i] for i in numpy.where(new_neighbors == 0)[0]])
                    infec_target = [new_coords[0] + infec_newcoord[0], new_coords[1] + infec_newcoord[1]]
                    if all(v < 50 for v in infec_target):
                        infec_prob = numpy.random.binomial(n=1, p=0.75)
                        if infec_prob == 1 and mat[infec_target[0], infec_target[1]] == 0:
                            mat[infec_target[0], infec_target[1]] = 1
                            spores.update({new_coords: None})
                            break

    new_walkers = {k: v for k, v in spores.items() if v is not None}
    spores.clear()
    spores.update(new_walkers)

    return mat, spores

###################################################################
################ Put it all together ##############################
###################################################################

#Specify landscape parameters
matrix_size = 100
n_patches = 30
n_draws = 50
deforest = [10, 35, 50, 65, 80, 95]
deforest_disp = [2, 2.5, 3, 3.5, 4, 4.5]
deforest_draws = 30

#Specify number of landscapes and time steps
n = 50
t = 1000

#Automate this function so it loops through all scenarios
def THE_FUNCTION(nlandscape = n):
    # Create blank array to store results
    perc_inf = numpy.empty((t, 2, n))
    for i in range(n):
        # Create landscapes
        (coffee, landscape) = MakeLandscape(size=matrix_size, patches=n_patches, draws=n_draws, deforest=deforest[0],
                                            disp=[0], ddraws=deforest_draws)

        # Create list of tuples for changing coords
        coord_change = [(-1, -1), (-1, 0), (-1, 1),
                        (0, -1), (0, 1),
                        (1, -1), (1, 0), (1, 1)]

        for j in range(t):
            coffee = cellaut(mat=coffee, land=landscape)
            walkers = new_spore(mat=coffee, coord=coord_change)
            (coffee, walkers) = spore_walk(mat=coffee, land=landscape, spores=walkers, coord=coord_change)
            perc_inf[j, 0, i] = j
            perc_inf[j, 1, i] = numpy.count_nonzero(coffee == 1) / (numpy.count_nonzero(coffee == 1) + numpy.count_nonzero(coffee == 0))
            print("j = " + str(j))

        print("i = " + str(i))

    return perc_inf

perc_inf = THE_FUNCTION(nlandscape=n)

perc_inf2 = perc_inf.transpose(2,0,1).reshape(-1, perc_inf.shape[1])

numpy.savetxt("def10disp2.csv", perc_inf2, delimiter=",")

