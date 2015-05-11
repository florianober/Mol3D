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
from astropy.io import fits as pf

import numpy as np
import helper as hlp
import matplotlib.pyplot as plt


#-------------------------------------------------------------------------------
# global paths
PATH_INPUT_GRID = 'input/grid/'

#-------------------------------------------------------------------------------
# global model definitions
MODEL_NAME = 'example' # must not be the Mol3D project name
R_IN = 1.0
R_OU = 120.0

GRID_NAME = 'spherical'                             # a = r,   b = th, c = ph
#~ GRID_NAME = 'cylindrical'                           # a = rho, b = ph, c = z
#~ GRID_NAME = 'cartesian'                             # a = x,   b = y,  c = z

N_A = 100              # no of coordinate a
N_B = 101              # no of coordinate b
N_C = 100                # no of coordinate c

#~ DENSITY_DISTRIBUTION = 'Disk'     # density for a disk
#~ DENSITY_DISTRIBUTION = 'Sphere'   # density for a sphere
#~ DENSITY_DISTRIBUTION = 'User'     # density for user input
DENSITY_DISTRIBUTION = 'fosite'     # density for user input

#~ VELOCITY_FIELD = 'Keplerian'
VELOCITY_FIELD = 'fosite'

FOSITE_FILE = '../Mol3d/output/planet-kley_0030.bin'

if DENSITY_DISTRIBUTION == 'fosite':
    f = hlp.read_fosite(FOSITE_FILE)
    # make sure, that the units are correct
    #
    print ('Fosite data:')

    MODEL_NAME += '_' + FOSITE_FILE[-8:-4]
    if f['/config/physics/units'] == 1:
        print('Found SI units')
        # SI units -> convert m to AU
        convert_unit = hlp.AU
        
        R_IN = round(f['/config/mesh/xmin']/convert_unit, 2)
        R_OU = round(f['/config/mesh/xmax']/convert_unit, 2)
        
        print('Rescaling the inner edge: %2.2f AU' %R_IN)
        print('Rescaling the outer edge: %2.2f AU' %R_OU)
        
    elif f['/config/physics/units'] == 3:
        print('Found geometric units')
        # scale free -> convert maximum to R_ou in AU
        convert_unit = f['/config/mesh/xmax']/R_OU
        R_IN = f['/config/mesh/xmin']/ convert_unit
        print('Rescaling the inner edge:', R_IN)
        #~ print(R_OU, f['/config/mesh/xmax']/convert_unit)
    else:
        print('Units are unknown')
        convert_unit = 1.
    
    X1 = (f['/mesh/bary_centers'][:,:,0]/convert_unit).flatten()
    Y1 = (f['/mesh/bary_centers'][:,:,1]/convert_unit).flatten()
    surface_density = f['/timedisc/density']
    
    if (f['/config/mesh/output/rotation']): 
        vxi = f['/timedisc/xvelocity']
        #~ omega = f['/config/timedisc/omega']
        omega = f['/config/sources/grav/pbinary/omega_rot']
        R = np.sqrt(f['/mesh/bary_centers'][:,:,0]**2 + f['/mesh/bary_centers'][:,:,1]**2)
        veta = f['/timedisc/yvelocity'] + omega * R
        rot = f['/mesh/rotation']
        VELOCITY_X = (np.cos(rot)*vxi + np.sin(rot)*veta).flatten()
        VELOCITY_Y = (-np.sin(rot)*vxi + np.cos(rot)*veta).flatten()
    else:
        VELOCITY_X = f['/timedisc/xvelocity']
        VELOCITY_Y = f['/timedisc/yvelocity']
        
    # make a plot of the fosite data
    #x_grid = f['/mesh/grid_x']/convert_unit
    #y_grid = f['/mesh/grid_y']/convert_unit
    #plt.figure('v_x')
    #plt.pcolormesh(x_grid, y_grid, vx)
    #plt.colorbar()
    #plt.figure('v_y')
    #plt.pcolormesh(x_grid, y_grid, vy)
    #plt.colorbar()
    #plt.figure('abs(v)')
    #plt.pcolormesh(x_grid, y_grid, np.sqrt(vy**2+vx**2))
    #plt.colorbar()
    #plt.figure('surface-density')
    #plt.pcolormesh(x_grid, y_grid, surface_density)
    #plt.colorbar()

    if f['/config/mesh/output/volume']:
        M_disk = (surface_density * f['/mesh/volume']).sum()/hlp.M_sun
        print('M_disk: %2.3g M_sun' %(M_disk))
    
    FOSITE_DENSITY = surface_density.flatten()

    if f['/config/sources/grav/output/height']:
        FOSITE_HEIGHT = f['/sources/grav/height'].flatten()/convert_unit
    else:
        FOSITE_HEIGHT = []
    print( '')

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

    xi = pos_xyz[0, :]
    yi = pos_xyz[1, :]
    if FOSITE_HEIGHT != []:
        h = griddata((X1, Y1), FOSITE_HEIGHT, (xi, yi), method='linear')
    else:
        scale_h = 10.0
        beta = 1.125
        h = scale_h * (r/100.0)**beta
    
    density = griddata((X1, Y1), FOSITE_DENSITY, (xi, yi),
                            method='linear') / h * \
                            np.exp(-0.5 * (z/h)**2, dtype=np.float64)


    return density

def get_density_user(pos_xyz):
    """
    a user interface 
    """

    density = 0.0

    return density

def get_density(pos_xyz):
    """
    stearing routine to get the density at the given cartesian coordinate 
    """

    if DENSITY_DISTRIBUTION == 'Disk':
        
        density = get_density_disk(pos_xyz)

    elif DENSITY_DISTRIBUTION == 'Sphere':
        if GRID_NAME != 'spherical':
            print("Warning, a %s grid might not be the best solution for a sphere" %(GRID_NAME))
        density = get_density_sphere(pos_xyz)

    elif DENSITY_DISTRIBUTION == 'fosite':

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

    return density
#-------------------------------------------------------------------------------
# velocity field (feel free to include more models)

def get_velocity_Keplerian(pos_xyz):
    """
    TBD!
    Keplerian rotation for one central star
    """
    z = pos_xyz[2,:]
    r = np.sqrt(pos_xyz[0,:]**2 + pos_xyz[1,:]**2)
    velocity_xyz = np.zeros((len(r),3))
    velocity_xyz[:, 0] = pos_xyz[0,:]
    velocity_xyz[:, 1] = pos_xyz[1,:]
    velocity_xyz[:, 2] = 0.0 

    return velocity_xyz

def get_velocity_fosite(pos_xyz):
    """
    a very(!) simple but working Fosite interface 
    """
    r = np.sqrt(pos_xyz[0, :]**2 + pos_xyz[1, :]**2)
    z = pos_xyz[2,:]

    xi = pos_xyz[0, :]
    yi = pos_xyz[1, :]

    velocity_xyz = np.zeros((len(r),3))
    velocity_xyz[:, 0] = griddata((X1, Y1), VELOCITY_X, (xi, yi), method='linear')
    velocity_xyz[:, 1] = griddata((X1, Y1), VELOCITY_Y, (xi, yi), method='linear')
    velocity_xyz[:, 2] = 0.0
    
    return velocity_xyz

def get_velocity_user(pos_xyz):
    """
    a user interface 
    """
    velocity_xyz = np.zeros_like(pos_xyz)


    return velocity_xyz

def get_velocity(pos_xyz):
    """
    stearing routine to get the density at the given cartesian coordinate 
    """

    if VELOCITY_FIELD == 'Keplerian':
        
        velocity = get_velocity_Keplerian(pos_xyz)

    elif VELOCITY_FIELD == 'fosite':

        print('Loading the velocity field from Fosite')
        velocity = get_velocity_fosite(pos_xyz)

    elif VELOCITY_FIELD == 'User':
        #TbD
        #velocity = get_velocity_user(pos_xyz) at position 'mid_point_cart'
        print("ERROR: user velocity definition is not implemeted")
        sys.exit()
    else:
        print("ERROR: velocity model unknown")
        sys.exit()

    return velocity

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
    data = np.zeros((N_Cells, 13), dtype=np.float64)

    # create density distribution:
    print("Calculating the midpoint coordinates for %d cells" %N_Cells)
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
                    #~ mid_point_cart[:, i_cell] = hlp.sp2ca(mid_point_coord[:, i_cell])
                    data[i_cell, 0:3] = hlp.sp2ca(mid_point_coord[:, i_cell])

                elif GRID_NAME == 'cylindrical':
                    #~ mid_point_cart[:, i_cell] = hlp.cy2ca(mid_point_coord[:, i_cell])
                    data[i_cell, 0:3] = hlp.cy2ca(mid_point_coord[:, i_cell])
                elif GRID_NAME == 'cartesian':
                    #~ mid_point_cart[:, i_cell] = 1.0* mid_point_coord[:, i_cell]
                    data[i_cell, 0:3] =mid_point_coord[:, i_cell]
                # increase the cell counter
                i_cell += 1
    print("Calculating the density distribution")
    data[:, 3] = get_density(data[:, 0:3])
    
    print("Calculating the velocity distribution")
    data[:, 4:7] = get_velocity(data[:, 0:3])
    # save results
    print("Saving the model for the use with Mol3D")

    file_name = PATH_INPUT_GRID + MODEL_NAME +  '_' + \
                DENSITY_DISTRIBUTION + '_' + \
                GRID_NAME
    model_file = open(file_name + '_model.dat', 'w')

    # write header
    model_file.write('%d \n' %(N_Cells))

    i_cell = 0
    for i_a in range(N_A):
        for i_b in range(N_B):
            for i_c in range(N_C):
                model_file.write('%d\t%.12e\t%.12e\t%.12e\t%.12e \n' %(data[i_cell, 0], data[i_cell, 1], data[i_cell, 2], data[i_cell, 3], data[i_cell, 4]))
                i_cell += 1
    model_file.close()

    # save results
    print("Saving the model for the use with Mol3D -> Fits")
    header = pf.Header()
    # write data array in header
    header['HIERARCH /config/dust_density'] = (1, 1)
    header['HIERARCH /config/velocity_x'] = (1, 2)
    header['HIERARCH /config/velocity_y'] = (1, 3)
    header['HIERARCH /config/velocity_z'] = (1, 4)
    header['N_DUST'] = 1
    hdu = pf.PrimaryHDU(data, header=header)
    hdulist = pf.HDUList(hdu)
    hdulist.writeto(file_name + '_model.fits', clobber=True)

    
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
plt.show()
