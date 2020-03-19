import numpy
import random
import matplotlib.pyplot as plt
import scipy.stats as stats
from nlmpy import nlmpy as nlm

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

def MakeLandscape(size, deforest, disp, cluster):
    #Create coffee matrix
    coffee = nlm.randomClusterNN(size, size, cluster, n = '8-neighbourhood')
    coffee = nlm.classifyArray(coffee, [0.25, 0.75])

    #Create landscape matrix
    landscape_clustered = nlm.randomElementNN(size, size, n=disp*2000, mask=coffee)
    landscape_clustered = nlm.classifyArray(landscape_clustered, [deforest, 1-deforest])

    ones = numpy.where(landscape_clustered == 1)
    zeroes = numpy.where(landscape_clustered == 0)

    landscape_clustered[ones[0], ones[1]] = numpy.random.uniform(0, 0.3, len(ones[1]))
    landscape_clustered[zeroes[0], zeroes[1]] = numpy.random.uniform(0.6, 1, len(zeroes[1]))

    landscape = landscape_clustered

    # Remove 1's from coffee
    coffee_ones = numpy.where(coffee == 1)
    coffee[coffee_ones[0], coffee_ones[1]] = numpy.nan

    #Initialize infection
    coffee_zeros = numpy.where(coffee == 0)
    randrow = random.randint(0, numpy.size(coffee_zeros, 1) - 1)

    coffee[coffee_zeros[0][randrow], coffee_zeros[1][randrow]] = 1

    start = numpy.where(coffee == 1)

    return (coffee, landscape, start)

(coffee, landscape, size) = MakeLandscape(size=100, deforest=0.45, disp=3, cluster=0.075)
plt.matshow(coffee)




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
    rowprobs = []
    for i in range(0, numpy.size(neighbors_clean, 0)):
        row = neighbors_clean[i,:]
        row = row[~numpy.isnan(row)]
        rowsums = row.sum()
        if rowsums > 0:
            add_row = numpy.random.beta(a=rowsums+1, b=8-rowsums+1, size=1)
            rowprobs.append(add_row)
        else:
            rowprobs.append(0)

    # Use number of infected neighbors to get success probability in a bernoulli trial
    bern_out = []
    for i in range(0,len(rowprobs)):
        out = numpy.random.binomial(n = 1, p = (0.1*rowprobs[i]))
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
                spores.append(new_coord)

    return spores

###################################################################
#################    Random Walk    ###############################
###################################################################

def spore_walk(spores, land, mat, coord):
    for i in range(0, len(spores)):
        step_credit = 10
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

            spores[i] = new_coords

            old_coords = new_coords

            if new_coords[0] > 99 or new_coords[1] > 99:
                break

            else:
                step_credit = step_credit-land[new_coords[0], new_coords[1]]

    for i in range(0, len(spores)):
        coords = list(spores[i])
        new_neighbors = neighbors_base(mat=mat, row=coords[0], col=coords[1])
        new_neighbors = numpy.array(new_neighbors)

        if numpy.any(new_neighbors == 0):
            infec_newcoord = random.choice([coord[i] for i in numpy.where(new_neighbors == 0)[0]])
            infec_target = [new_coords[0] + infec_newcoord[0], new_coords[1] + infec_newcoord[1]]
            if all(v < 50 for v in infec_target):
                infec_prob = numpy.random.binomial(n=1, p=0.5)
                if infec_prob == 1 and mat[infec_target[0], infec_target[1]] == 0:
                    mat[infec_target[0], infec_target[1]] = 1
                    spores[i] = None

    spores = list(filter(None, spores))

    return mat, spores

###################################################################
################ Put it all together ##############################
###################################################################

#Specify landscape parameters
matrix_size = 100
deforest = [0.15, 0.3, 0.45, 0.6, 0.75]
deforest_disp = [1, 2, 3, 4, 5]
disag = [0.1, 0.2, 0.3]

#Specify number of landscapes and time steps
n = 50
t = 1000

#Write the master function
def base_function(nlandscape = n):
    # Create blank array to store results
    perc_inf = numpy.empty((t, 5, n))

    for i in range(n):
        # Create landscapes
        (coffee, landscape, start) = MakeLandscape(size=matrix_size, deforest=defor, disp=disp, cluster=disagg)

        # Create list of tuples for changing coords
        coord_change = [(-1, -1), (-1, 0), (-1, 1),
                        (0, -1), (0, 1),
                        (1, -1), (1, 0), (1, 1)]

        for j in range(t):
            coffee = cellaut(mat=coffee, land=landscape)
            walkers = new_spore(mat=coffee, coord=coord_change)
            (coffee, walkers) = spore_walk(mat=coffee, land=landscape, spores=walkers, coord=coord_change)
            perc_inf[j, 0, i] = j
            perc_inf[j, 1, i] = numpy.count_nonzero(coffee == 1) / \
                                (numpy.count_nonzero(coffee == 1) + numpy.count_nonzero(coffee == 0))
            perc_inf[j, 2, i] = start[0]; perc_inf[j, 3, i] = start[1]
            print("j = " + str(j))

            #Add stopping point if landscape is fully infected
            if len(numpy.where(coffee == 0)[0]) == 0:
                break

        print("i = " + str(i))

    return perc_inf

for i in range(len(deforest)):
    for j in range(len(deforest_disp)):
        for k in range(len(disag)):
            defor = deforest[i]
            disp = deforest_disp[j]
            disagg = disag[k]

            perc_inf = base_function(nlandscape=n)

            perc_inf2 = perc_inf.transpose(2,0,1).reshape(-1, perc_inf.shape[1])

            filename = "def"+str(deforest[i])+"disp"+str(deforest_disp[j])+"disagg"+str(disag[k])+".csv"

            numpy.savetxt(filename, perc_inf2, delimiter=",")

