#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
with this script it is possible to generate a (disk) model for Mol3D

You only have to chose a suitable grid and define your model

author: Florian Ober
email: fober@astrophysik.uni-kiel.de

"""
#-------------------------------------------------------------------------------
# import some packages
#
import os
import sys
sys.path.append('./visualize/')
import time
import math
from scipy.interpolate import griddata

import numpy as np
import helper as hlp


#-------------------------------------------------------------------------------
# global paths
PATH_INPUT_GRID = 'input/grid/'

#-------------------------------------------------------------------------------
# global model definitions
MODEL_NAME = 'example' # must not be the Mol3D project name
R_IN = 1.0
R_OU = 100.0

GRID_NAME = 'spherical'                             # a = r,   b = th, c = ph
#~ GRID_NAME = 'cylindrical'                           # a = rho, b = ph, c = z
#~ GRID_NAME = 'cartesian'                             # a = x,   b = y,  c = z

N_A = 100              # no of coordinate a
N_B = 101              # no of coordinate b
N_C = 1                # no of coordinate c

DENSITY_DISTRIBUTION = 'Disk'     # density for a disk
#~ DENSITY_DISTRIBUTION = 'Sphere'   # density for a sphere
#~ DENSITY_DISTRIBUTION = 'User'     # density for user input
#~ DENSITY_DISTRIBUTION = 'Fosite'     # density for user input

FOSITE_FILE = 'bindisk_alphavis_cooling_0010.bin'
if DENSITY_DISTRIBUTION == 'Fosite':
    f = hlp.read_fosite(FOSITE_FILE)
    # make sure, if the units are correct
    # assumption here: we need to convert from 'm' to 'AU'
    #
    X1 = (f['/mesh/bary_centers'][:,:,0]/hlp.AU).flatten()
    Y1 = (f['/mesh/bary_centers'][:,:,1]/hlp.AU).flatten()

    density1 = f['/timedisc/density'].flatten()

    height1 = f['/sources/grav/height'].flatten()/hlp.AU


LINK = True # link the model and boundaries automatically

#-------------------------------------------------------------------------------
# density definitions (feel free to include more models)

def get_density_disk(pos_xyz):
    """
    a shakura & sunyaev disk
    """
    z = pos_xyz[2,:]
    r = np.sqrt(pos_xyz[0,:]**2 + pos_xyz[1,:]**2)
    scale_h = 10.0
    alpha = -2.625
    beta = 1.125

    density = np.zeros_like(r)
    ind = (r >= R_IN) & (r < R_OU)

    h = scale_h * (r/100.0)**beta
    
    density[ind] = (r[ind]/100.0)**(alpha) * np.exp(-0.5 * (z[ind]/h[ind])**2, dtype=np.float64)
    

    return density

def get_density_sphere(pos_xyz):
    """
    a simple sphere with inner an outer radius and constant
    density
    """

    r = np.sqrt(pos_xyz[0, :]**2 +
                pos_xyz[1, :]**2 +
                pos_xyz[2, :]**2)
    density = np.zeros_like(r)
    ind = (r >= R_IN) & (r < R_OU)
    density[ind] = 1.0
    return density

def get_density_fosite(pos_xyz):
    """
    a very(!) simple but working Fosite interface 
    """
    r = np.sqrt(pos_xyz[0, :]**2 + pos_xyz[1, :]**2)
    z = pos_xyz[2,:]
    density = np.zeros_like(r)
    ind = (r >= R_IN) & (r < R_OU)
    xi = pos_xyz[0, :]
    yi = pos_xyz[1, :]
    h = griddata((X1, Y1), height1, (xi[ind], yi[ind]), method='linear')
    
    density[ind] = griddata((X1, Y1), density1, (xi[ind], yi[ind]), method='linear') * \
                   np.exp(-0.5 * (z[ind]/h)**2, dtype=np.float64)


    return density

def get_density_user(pos_xyz):
    """
    a user interface 
    """

    density = 0.0

    return density

FIRST_CALL = True
def get_density(pos_xyz):
    """
    stearing routine to get the density at the given cartesian coordinate 
    """
    global FIRST_CALL

    if DENSITY_DISTRIBUTION == 'Disk':
        
        density = get_density_disk(pos_xyz)

    elif DENSITY_DISTRIBUTION == 'Sphere':
        if FIRST_CALL and GRID_NAME != 'spherical':
            print("Warning, a %s grid might not be the best solution for a sphere" %(GRID_NAME))
        density = get_density_sphere(pos_xyz)

    elif DENSITY_DISTRIBUTION == 'Fosite':
        #TbD
        if FIRST_CALL:
            print('Loading the density from Fosite')
        density = get_density_fosite(pos_xyz)

    elif DENSITY_DISTRIBUTION == 'User':
        #TbD
        #density = get_density_user(pos_xyz) at position 'mid_point_cart'
        print("ERROR: user density definition is not implemeted")
        sys.exit()
    else:
        print("ERROR: density model unknown")
        sys.exit()
    FIRST_CALL = False

    return density

#-------------------------------------------------------------------------------
# main program

def main():
    """
    set up the grid (boundaries) and density at the cell center
    """
    global N_A, N_B, N_C
    print("Hello")
    print("")
    
    # setting boundaries (including some error checks)
    # R = (a, b, c)
    if GRID_NAME == 'spherical':
        # position vector R:
        # R = (a=r, b=theta, c=phi)
        
        # log r coordinate (R_IN,...,R_OU)
        if N_A < 20:
            print("ERROR, not enough cells in r-direction: %d" %(N_A))
            sys.exit()
        a_bounds = np.logspace(np.log10(R_IN), np.log10(R_OU), N_A + 1, dtype=np.float64)

        # linear theta coordinate (-PI/2,...,+PI/2)
        if not(N_B%2):
            N_B += 1
            print("Warning, adjusting theta cells to %d" %(N_B))
        if N_B < 21:
            print("ERROR, not enough cells in theta-direction: %d < 21" %(N_B))
            sys.exit()
        b_bounds = np.linspace(-np.pi/2.0, np.pi/2.0, N_B + 1, dtype=np.float64)

        # linear phi coordinate (0,...,+2*PI)
        if N_C%2 and not(N_C == 1):
            N_C += 1
            print("Warning, adjusting phi cells to %d" %(N_C))
        c_bounds = np.linspace(0.0, 2.0*np.pi, N_C + 1, dtype=np.float64)
        
    elif GRID_NAME == 'cylindrical':
        # position vector R:
        # R = (a=rho, b=phi, c=z)
        
        # log rho coordinate (R_IN,...,R_OU)
        if N_A < 20:
            print("ERROR, not enough cells in rho-direction: %d < 20" %(N_A))
            sys.exit()
        a_bounds = np.logspace(np.log10(R_IN), np.log10(R_OU), N_A + 1, dtype=np.float64)

        # linear phi coordinate (0,...,+2*PI)
        if N_B%2 and not(N_B == 1):
            N_B += 1
            print("Warning, adjusting phi cells to %d" %(N_B))
        b_bounds = np.linspace(0.0, 2.0*np.pi, N_B + 1, dtype=np.float64)

        # sinh z corrdinate (-R_OU,...,R_OU)
        if not(N_C%2):
            N_C += 1
            print("Warning, adjusting number of z cells to %d" %(N_C))
        if N_C < 21:
            print("ERROR, not enough cells in z-direction: %d < 21" %(N_C))
            sys.exit()
        c_bounds = np.zeros(N_C + 1, dtype=np.float64)
        k = 5 # this is a scaling parameter, adjust if necessary
        for i in range(-int((N_C-1)/2), int((N_C-1)/2)+2):
            c_bounds[i+int((N_C-1)/2)] = R_OU*np.sinh(k*(i-0.5)/((N_C-1)/2+0.5))/np.sinh(k)
    elif GRID_NAME == 'cartesian':
        # position vector R:
        # R = (a=x, b=y, c=z)
        print ("Warning: the cartesian grid is not yet tested in Mol3D")

        # sinh x corrdinate (-R_OU,...,R_OU)
        if not(N_A%2):
            N_A += 1
            print("Warning, adjusting number of x cells to %d" %(N_A))
        if N_A < 21:
            print("ERROR, not enough cells in x-direction: %d < 21" %(N_A))
            sys.exit()
        a_bounds = np.zeros(N_A + 1, dtype=np.float64)
        k = 5 # this is a scaling parameter, adjust if necessary
        for i in range(-int((N_A-1)/2), int((N_A-1)/2)+2):
            a_bounds[i+int((N_A-1)/2)] = R_OU*np.sinh(k*(i-0.5)/((N_A-1)/2+0.5))/np.sinh(k)
            
        # sinh y corrdinate (-R_OU,...,R_OU)
        if not(N_B%2):
            N_B += 1
            print("Warning, adjusting number of y cells to %d" %(N_B))
        if N_B < 21:
            print("ERROR, not enough cells in y-direction: %d < 21" %(N_B))
            sys.exit()
        b_bounds = np.zeros(N_B + 1, dtype=np.float64)
        k = 5 # this is a scaling parameter, adjust if necessary
        for i in range(-int((N_B-1)/2), int((N_B-1)/2)+2):
            b_bounds[i+int((N_B-1)/2)] = R_OU*np.sinh(k*(i-0.5)/((N_B-1)/2+0.5))/np.sinh(k)

        # sinh z corrdinate (-R_OU,...,R_OU)
        if not(N_C%2):
            N_C += 1
            print("Warning, adjusting number of z cells to %d" %(N_C))
        if N_C < 21:
            print("ERROR, not enough cells in z-direction: %d < 21" %(N_C))
            sys.exit()   
        c_bounds = np.zeros(N_C + 1, dtype=np.float64)
        k = 5 # this is a scaling parameter, adjust if necessary
        for i in range(-int((N_C-1)/2), int((N_C-1)/2)+2):
            c_bounds[i+int((N_C-1)/2)] = R_OU*np.sinh(k*(i-0.5)/((N_C-1)/2+0.5))/np.sinh(k)

    else:
        print("ERROR: Coordinate system not defined")
        sys.exit()

    # save boundaries
    a_bounds_file = open(PATH_INPUT_GRID + MODEL_NAME + '_' +
                                           DENSITY_DISTRIBUTION + '_' +
                                           GRID_NAME + '_a_boundaries.dat', 'w')
    b_bounds_file = open(PATH_INPUT_GRID + MODEL_NAME + '_' +
                                           DENSITY_DISTRIBUTION + '_' +
                                           GRID_NAME + '_b_boundaries.dat', 'w')
    c_bounds_file = open(PATH_INPUT_GRID + MODEL_NAME + '_' +
                                           DENSITY_DISTRIBUTION + '_' +
                                           GRID_NAME + '_c_boundaries.dat', 'w')

    # write headers
    a_bounds_file.write('%d \n' %(N_A))
    a_bounds_file.write('noscale \n')
    for i_a in range(N_A+1):
        a_bounds_file.write('%.12e \n' %a_bounds[i_a])

    b_bounds_file.write('%d \n' %(N_B))
    b_bounds_file.write('noscale \n')
    for i_b in range(N_B+1):
        b_bounds_file.write('%.12e \n' %b_bounds[i_b])
        
    c_bounds_file.write('%d \n' %(N_C))
    c_bounds_file.write('noscale \n')
    for i_c in range(N_C+1):
        c_bounds_file.write('%.12e \n' %c_bounds[i_c])

    a_bounds_file.close()
    b_bounds_file.close()
    c_bounds_file.close()
    
    N_Cells = N_A * N_B * N_C
    data = np.zeros((N_Cells, 2), dtype=np.float64)

    # create density distribution:
    print("Calculating the mitpoint coordinates for %d cells" %N_Cells)
    i_cell = 0
    mid_point_coord = np.zeros((3, N_Cells), dtype=np.float64)
    mid_point_cart = np.zeros((3, N_Cells), dtype=np.float64)
    for i_a in range(N_A):
        for i_b in range(N_B):
            for i_c in range(N_C):
                data[i_cell, 0] = i_cell+1 # Mol3D/Fortran counts from 1

                # first get the cell midpoint coordinate

                mid_point_coord[0, i_cell] = (a_bounds[i_a] + a_bounds[i_a+1])/2.0
                mid_point_coord[1, i_cell] = (b_bounds[i_b] + b_bounds[i_b+1])/2.0
                mid_point_coord[2, i_cell] = (c_bounds[i_c] + c_bounds[i_c+1])/2.0

                # convert to cartesian
                if GRID_NAME == 'spherical':
                    mid_point_cart[:, i_cell] = hlp.sp2ca(mid_point_coord[:, i_cell])

                elif GRID_NAME == 'cylindrical':
                    mid_point_cart[:, i_cell] = hlp.cy2ca(mid_point_coord[:, i_cell])

                elif GRID_NAME == 'cartesian':
                    mid_point_cart[:, i_cell] = 1.0* mid_point_coord[:, i_cell]

                # increase the cell counter
                i_cell += 1
    print("Calculating the density distribution")
    data[:, 1] = get_density(mid_point_cart)
    
    # save results
    print("Saving the model for the use with Mol3D")
    model_file = open(PATH_INPUT_GRID + MODEL_NAME +  '_' +
                                        DENSITY_DISTRIBUTION + '_' +
                                        GRID_NAME + '_model.dat', 'w')

    # write header
    model_file.write('%d \n' %(N_Cells))

    i_cell = 0
    for i_a in range(N_A):
        for i_b in range(N_B):
            for i_c in range(N_C):
                model_file.write('%d\t%.12e \n' %(data[i_cell, 0], data[i_cell, 1]))
                i_cell += 1
    model_file.close()

    print("")
    print("%s model '%s' succesfully generated on a %s grid"
                                %(DENSITY_DISTRIBUTION, MODEL_NAME, GRID_NAME))
    print("")
    if not(LINK):
        print("")
        print("Don't forget to symlink to the gernerated files: i.e., ")
        print("      ln -sf %s_%s_%s_a_boundaries.dat a_boundaries.dat"
                                %(MODEL_NAME, DENSITY_DISTRIBUTION, GRID_NAME))
        print("      ln -sf %s_%s_%s_b_boundaries.dat b_boundaries.dat"
                                %(MODEL_NAME, DENSITY_DISTRIBUTION, GRID_NAME))
        print("      ln -sf %s_%s_%s_c_boundaries.dat c_boundaries.dat"
                                %(MODEL_NAME, DENSITY_DISTRIBUTION, GRID_NAME))
        print("      ln -sf %s_%s_%s_model.dat model.dat"
                                %(MODEL_NAME, DENSITY_DISTRIBUTION, GRID_NAME))
        print("      (in the 'input/grid/' directory)")
    else:
        os.system("ln -sf %s_%s_%s_a_boundaries.dat %sa_boundaries.dat"
                %(MODEL_NAME, DENSITY_DISTRIBUTION, GRID_NAME, PATH_INPUT_GRID))
        os.system("ln -sf %s_%s_%s_b_boundaries.dat %sb_boundaries.dat"
                %(MODEL_NAME, DENSITY_DISTRIBUTION, GRID_NAME, PATH_INPUT_GRID))
        os.system("ln -sf %s_%s_%s_c_boundaries.dat %sc_boundaries.dat"
                %(MODEL_NAME, DENSITY_DISTRIBUTION, GRID_NAME, PATH_INPUT_GRID))
        os.system("ln -sf %s_%s_%s_model.dat %smodel.dat"
                %(MODEL_NAME, DENSITY_DISTRIBUTION, GRID_NAME, PATH_INPUT_GRID))
        print("Symlinks created")

    print("Use GRID_NAME: '%s' and grid_type: '9' in the Mol3D input file"
                                %(GRID_NAME))
    print("")
    print("bye bye")

main()
