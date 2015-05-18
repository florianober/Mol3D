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
from matplotlib.colors import LogNorm


#-------------------------------------------------------------------------------
# global paths
PATH_INPUT_GRID = 'input/grid/'

#-------------------------------------------------------------------------------
# global model definitions
MODEL_NAME = 'example' # must not be the Mol3D project name
R_IN = 10
R_OU = 200.0
M_STAR = 0.7*hlp.M_sun
M_DISK = 1e-2*hlp.M_sun

COORDINATE_SYSTEM = 'spherical'                             # a = r,   b = th, c = ph
#~ COORDINATE_SYSTEM = 'cylindrical'                           # a = rho, b = ph, c = z
#~ COORDINATE_SYSTEM = 'cartesian'                             # a = x,   b = y,  c = z

N_A = 100              # no of coordinate a
N_B = 101              # no of coordinate b
N_C = 100                # no of coordinate c

DENSITY_DISTRIBUTION = 'Disk'     # density for a disk
#~ DENSITY_DISTRIBUTION = 'Sphere'   # density for a sphere
#~ DENSITY_DISTRIBUTION = 'User'     # density for user input
#~ DENSITY_DISTRIBUTION = 'fosite'     # density for user input

VELOCITY_FIELD = 'Keplerian'
#~ VELOCITY_FIELD = 'fosite'

TEMPERATURE_DISTRIBUTION = 'no'
#~ TEMPERATURE_DISTRIBUTION = 'power'
#~ TEMPERATURE_DISTRIBUTION = 'fosite'

FOSITE_FILE = ''

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

        # rescale if necessary
        #~ R_IN = round(f['/config/mesh/xmin']/convert_unit, 2)
        #~ R_OU = round(f['/config/mesh/xmax']/convert_unit, 2)
        #~ 
        #~ print('Rescaling the inner edge: %2.2f AU' %R_IN)
        #~ print('Rescaling the outer edge: %2.2f AU' %R_OU)
        
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
        # get velocity
        curv_vx = f['/timedisc/xvelocity']
        curv_vy = f['/timedisc/yvelocity']
        # now convert to cartesian coordinates
        rot = f['/mesh/rotation']
        # routine from Fosite
        cart_vx = np.cos(rot) * curv_vx + np.sin(rot) * curv_vy
        cart_vy = -np.sin(rot) * curv_vx + np.cos(rot) * curv_vy

    else:
        cart_vx = f['/timedisc/xvelocity']
        cart_vy = f['/timedisc/yvelocity']

    # now correct for the rotation frame if necessary and available
    # first we calculate the velocity vector of the rotating frame
    # Just check if an omega value is given somewhere, it is not
    # generalized in fosite at the moment
    # 
    omega_array =['/config/sources/grav/pbinary/omega_rot',
                   '/config/sources/grav/binary/omega_rot',
                   '/config/sources/rotframe/omega']

    omega = 0.0
    for item in omega_array:
        if item in f.keys():
            omega = max(omega, f[item])

    if omega > 0.0:
        print('Rotating Frame found, omega: %2.2g 1/s' %(omega) )
        
        rot_vx = - f['/mesh/bary_centers'][:,:,1] * omega
        rot_vy =   f['/mesh/bary_centers'][:,:,0] * omega
    else:
        rot_vx = 0.0
        rot_vy = 0.0
    # now we only need a vector addition
    VELOCITY_X = (cart_vx + rot_vx).flatten()
    VELOCITY_Y = (cart_vy + rot_vy).flatten()
        
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
        print('M_disk (fosite): %2.3g M_sun' %(M_disk))

    FOSITE_DENSITY = surface_density.flatten()

    if f['/config/sources/grav/output/height']:
        FOSITE_HEIGHT = f['/sources/grav/height'].flatten()/convert_unit
    else:
        FOSITE_HEIGHT = []
        print('No height found fosite data set')
        
    print( '')
SHOW_MODEL = True # show model distribution
LINK = True # link the model and boundaries automatically

#-------------------------------------------------------------------------------
def in_model_space(pos_xyz):
    """
    Test if the given position is inside model space
    returns a mask 'ind'
    """
    if COORDINATE_SYSTEM == 'spherical':
        r = np.sqrt(pos_xyz[:, 0]**2 + pos_xyz[:, 1]**2 +pos_xyz[:, 2]**2)
        ind = (r >= R_IN) & (r < R_OU)
    elif COORDINATE_SYSTEM == 'cylindrical':
        r = np.sqrt(pos_xyz[:, 0]**2 + pos_xyz[:, 1]**2)
        z = np.abs(pos_xyz[:, 2])
        ind = (r >= R_IN) & (r < R_OU) & (z < R_OU)
    elif COORDINATE_SYSTEM == 'cartesian':
        x = np.abs(pos_xyz[:, 0])
        y = np.abs(pos_xyz[:, 1])
        z = np.abs(pos_xyz[:, 2])
        ind =  (x < R_OU) & (y < R_OU) & (z < R_OU)
    else:
        print("ERROR: could not find coordinate system (in_model_space)")
    return ind

#-------------------------------------------------------------------------------
# density definitions (feel free to include more models)

def get_density_disk(pos_xyz):
    """
    a shakura & sunyaev disk
    """
    z = pos_xyz[:, 2]
    r = np.sqrt(pos_xyz[:, 0]**2 + pos_xyz[:, 1]**2)
    scale_h = 10.0
    alpha = -2.625
    beta = 1.125

    density = np.zeros_like(r)

    h = scale_h * (r/100.0)**beta
    ind = (h > 0) & (r < R_OU) & (r >= R_IN)
    
    density[ind] = (r[ind]/100.0)**(alpha) * np.exp(-0.5 * (z[ind]/h[ind])**2, dtype=np.float64)

    return density

def get_density_sphere(pos_xyz):
    """
    a simple sphere with inner an outer radius and constant
    density
    """
    r = np.sqrt(pos_xyz[:, 0]**2 + pos_xyz[:, 1]**2 + pos_xyz[:, 2]**2)
    ind = (r < R_OU) & (r >= R_IN)
    density = np.zeros_like(r)
    density[ind] = 1.0
    return density

def get_density_fosite(pos_xyz):
    """
    a very(!) simple but working Fosite interface 
    """
    z = pos_xyz[:, 2]
    r = np.sqrt(pos_xyz[:, 0]**2 + pos_xyz[:, 1]**2)
    density = np.zeros_like(r)

    xi = pos_xyz[:, 0]
    yi = pos_xyz[:, 1]
    if FOSITE_HEIGHT != []:
        h = griddata((X1, Y1), FOSITE_HEIGHT, (xi, yi), method='linear', fill_value=0)
    else:
        scale_h = 10.0
        beta = 1.125
        h = scale_h * (r/100.0)**beta
    ind = (h > 0)
    density[ind] = griddata((X1, Y1), FOSITE_DENSITY, (xi[ind], yi[ind]),
                            method='linear', fill_value=0) / h[ind] * \
                            np.exp(-0.5 * (z[ind]/h[ind])**2, dtype=np.float64)

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
    ind = in_model_space(pos_xyz)
    density = np.zeros(len(pos_xyz[:, 0]))
    if DENSITY_DISTRIBUTION == 'Disk':
        density[ind] = get_density_disk(pos_xyz[ind])

    elif DENSITY_DISTRIBUTION == 'Sphere':
        if COORDINATE_SYSTEM != 'spherical':
            print("Warning, a %s grid might not be the best solution for a sphere" %(COORDINATE_SYSTEM))
        density[ind] = get_density_sphere(pos_xyz[ind])

    elif DENSITY_DISTRIBUTION == 'fosite':

        density[ind] = get_density_fosite(pos_xyz[ind])

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
    Keplerian rotation profile for one central star with given Mass
    at given location (TdB, for now we assume the star is in the centrum of the
    model space)
    """

    r = (pos_xyz[:, 0]**2+pos_xyz[:, 1]**2 + 1e-16)**(-1.5*0.5)
    konst = (hlp.gamma*M_STAR/hlp.AU)**0.5
    
    velocity_xyz = np.zeros((len(r),3))
    velocity_xyz[:, 0] = - pos_xyz[:, 1] * r * konst
    velocity_xyz[:, 1] =   pos_xyz[:, 0] * r * konst
    velocity_xyz[:, 2] =   0.0 

    return velocity_xyz

def get_velocity_fosite(pos_xyz):
    """
    a very(!) simple but working Fosite interface 
    """

    xi = pos_xyz[:, 0]
    yi = pos_xyz[:, 1]

    velocity_xyz = np.zeros((len(xi),3))
    velocity_xyz[:, 0] = griddata((X1, Y1), VELOCITY_X, (xi, yi), method='linear', fill_value=0)
    velocity_xyz[:, 1] = griddata((X1, Y1), VELOCITY_Y, (xi, yi), method='linear', fill_value=0)
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

    ind = in_model_space(pos_xyz)
    velocity = np.zeros_like(pos_xyz)

    if VELOCITY_FIELD == 'Keplerian':
        
        velocity[ind] = get_velocity_Keplerian(pos_xyz[ind])

    elif VELOCITY_FIELD == 'fosite':

        velocity[ind] = get_velocity_fosite(pos_xyz[ind])

    elif VELOCITY_FIELD == 'User':
        #TbD
        #velocity[ind] = get_velocity_user(pos_xyz[ind]) at position 'mid_point_cart'
        print("ERROR: user velocity definition is not implemeted")
        sys.exit()
    else:
        print("ERROR: velocity model unknown")
        sys.exit()

    return velocity
    
#-------------------------------------------------------------------------------
# temperature distribution (feel free to include more models)

def get_temperature_power(pos_xyz):
    """
    the temperature is following a powerlaw
    """

    r = np.sqrt(pos_xyz[:, 0]**2+pos_xyz[:, 1]**2)
    exp = (-1.5)
    ind = (r > 0)
    temp = np.zeros((len(r)))
    temp[ind] = r[ind]**(exp)

    return temp

def get_temperature_fosite(pos_xyz):
    """
    a very(!) simple but working Fosite interface 
    """

    xi = pos_xyz[:, 0]
    yi = pos_xyz[:, 1]

    temp = np.zeros_like(xi)
    temp = griddata((X1, Y1), FOSITE_TEMPERATURE, (xi, yi), method='linear', fill_value=0)
    
    return temp

def get_temperature_user(pos_xyz):
    """
    a user interface 
    """
    temp = np.zeros_like(pos_xyz)

    return temp

def get_temperature(pos_xyz):
    """
    stearing routine to get the density at the given cartesian coordinate 
    """

    ind = in_model_space(pos_xyz)
    temperature = np.zeros(len(pos_xyz[:, 0]))

    if TEMPERATURE_DISTRIBUTION == 'no':
        pass
    elif TEMPERATURE_DISTRIBUTION == 'power':

        temperature[ind] = get_temperature_power(pos_xyz[ind])

    elif TEMPERATURE_DISTRIBUTION == 'fosite':

        temperature[ind] = get_temperature_fosite(pos_xyz[ind])

    elif TEMPERATURE_DISTRIBUTION == 'User':
        #TbD
        #temperature[ind] = get_temperature_user(pos_xyz[ind]) at position 'mid_point_cart'
        print("ERROR: user temperature definition is not implemeted")
        sys.exit()
    else:
        print("ERROR: temperature model unknown")
        sys.exit()

    return temperature

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
    if COORDINATE_SYSTEM == 'spherical':
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
        
    elif COORDINATE_SYSTEM == 'cylindrical':
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
    elif COORDINATE_SYSTEM == 'cartesian':
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
                                           COORDINATE_SYSTEM + '_a_boundaries.dat', 'w')
    b_bounds_file = open(PATH_INPUT_GRID + MODEL_NAME + '_' +
                                           DENSITY_DISTRIBUTION + '_' +
                                           COORDINATE_SYSTEM + '_b_boundaries.dat', 'w')
    c_bounds_file = open(PATH_INPUT_GRID + MODEL_NAME + '_' +
                                           DENSITY_DISTRIBUTION + '_' +
                                           COORDINATE_SYSTEM + '_c_boundaries.dat', 'w')

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
                if COORDINATE_SYSTEM == 'spherical':
                    #~ mid_point_cart[:, i_cell] = hlp.sp2ca(mid_point_coord[:, i_cell])
                    data[i_cell, 0:3] = hlp.sp2ca(mid_point_coord[:, i_cell])

                elif COORDINATE_SYSTEM == 'cylindrical':
                    #~ mid_point_cart[:, i_cell] = hlp.cy2ca(mid_point_coord[:, i_cell])
                    data[i_cell, 0:3] = hlp.cy2ca(mid_point_coord[:, i_cell])
                elif COORDINATE_SYSTEM == 'cartesian':
                    #~ mid_point_cart[:, i_cell] = 1.0* mid_point_coord[:, i_cell]
                    data[i_cell, 0:3] = mid_point_coord[:, i_cell]
                # increase the cell counter
                i_cell += 1
    header = pf.Header()

    print("Calculating the (dust) density distribution")
    data[:, 3] = get_density(data[:, 0:3])
    # Example for the data arangement:
    #      |        KEYWORD             |    | Position in data cube
    header['HIERARCH Mol3D_dust_density1'] = 4

    # H2 density distribution (total)
    data[:, 5] = data[:, 3] * 100.0
    header['HIERARCH Mol3D_col_density1'] = 6
    # H2 density distribution (para)
    data[:, 6] = 0.25 * data[:, 4]
    header['HIERARCH Mol3D_col_density2'] = 7
    # H2 density distribution (ortho)
    data[:, 7] = 0.75 * data[:, 4]
    header['HIERARCH Mol3D_col_density3'] = 8

    # molecule density distribution
    abundance = 1e-5
    data[:, 4] = data[:, 5] * abundance
    header['HIERARCH Mol3D_mol_density'] = 5

    # dust temperature (if given, no temperature calculation will be done TbD)
    data[:, 8] = get_temperature(data[:, 0:3])
    header['HIERARCH Mol3D_dust_temp1'] = 9
    
    # gas temperature
    data[:, 9] = get_temperature(data[:, 0:3])
    header['HIERARCH Mol3D_gas_temp'] = 10

    print("Calculating the velocity distribution")
    
    data[:, 10:13] = get_velocity(data[:, 0:3])
    header['HIERARCH Mol3D_velocity_x'] = 11
    header['HIERARCH Mol3D_velocity_y'] = 12
    header['HIERARCH Mol3D_velocity_z'] = 13

    # save results
    print("Saving the model for the use with Mol3D -> Fits")
    file_name = PATH_INPUT_GRID + MODEL_NAME +  '_' + \
                DENSITY_DISTRIBUTION + '_' + \
                COORDINATE_SYSTEM

    # write data array in header and fits file
    header['N_DUST'] = 1
    hdu = pf.PrimaryHDU(data, header=header)
    hdulist = pf.HDUList(hdu)
    hdulist.writeto(file_name + '_model.fits', clobber=True)

    print("")
    print("%s model '%s' succesfully generated on a %s grid"
                                %(DENSITY_DISTRIBUTION, MODEL_NAME, COORDINATE_SYSTEM))
    print("")

    # eye candy: show distributions
    if SHOW_MODEL:
        print("Visualisation of the model")
        # xy-plane
        N = 401
        xx = np.linspace(-R_OU, R_OU, N)
        yy = np.linspace(-R_OU, R_OU, N)
        data = np.zeros((N**2))
        extent = [-R_OU, R_OU, -R_OU, R_OU]
        counter = 0
        pos_xyz = np.zeros((N*N, 3))
        for i in range(N):
            for j in range(N):
                pos_xyz[counter, 0] = xx[j]
                pos_xyz[counter, 1] = yy[i]
                pos_xyz[counter, 2] = 0.0
                counter += 1
        
        plt.figure('Normalize Dust density xy-plane')
        data = get_density(pos_xyz)
        data /= data.max()
        vmax = data.max()
        vmin = vmax*1e-4
        plt.imshow(data.reshape(N,N),
                   origin='lower',
                   norm=LogNorm(vmin=vmin,vmax=vmax),
                   extent=extent, interpolation='None')
        plt.colorbar().set_label('Normalized dust-density')
        
        plt.figure('abs(velocity) xy-plane')
        velo = get_velocity(pos_xyz)
        data = np.sqrt(velo[:, 0]**2 + velo[:, 1]**2 + velo[:, 2]**2)
        plt.imshow(data.reshape(N,N)/1000, origin='lower',
                   extent=extent, interpolation='None')
        plt.colorbar().set_label('velocity [km/s]')


        # xz-plane
        data = np.zeros((N**2))
        counter = 0
        pos_xyz = np.zeros((N*N, 3))
        for i in range(N):
            for j in range(N):
                pos_xyz[counter, 0] = xx[j]
                pos_xyz[counter, 1] = 0.0
                pos_xyz[counter, 2] = yy[i]
                counter += 1
        
        plt.figure('Normalize Dust density xz-plane')
        data = get_density(pos_xyz)
        data /= data.max()
        vmax = data.max()
        vmin = vmax*1e-4
        plt.imshow(data.reshape(N,N),
                   origin='lower',
                   norm=LogNorm(vmin=vmin,vmax=vmax),
                   extent=extent, interpolation='None')
        plt.colorbar().set_label('Normalized dust-density')
        
        plt.figure('abs(velocity) xz-plane')
        velo = get_velocity(pos_xyz)
        data = np.sqrt(velo[:, 0]**2 + velo[:, 1]**2 + velo[:, 2]**2)
        plt.imshow(data.reshape(N,N)/1000, origin='lower',
                   extent=extent, interpolation='None')
        plt.colorbar().set_label('velocity [km/s]')

    if not(LINK):
        print("")
        print("Don't forget to symlink to the gernerated files: i.e., ")
        print("      ln -sf %s_%s_%s_a_boundaries.dat a_boundaries.dat"
                                %(MODEL_NAME, DENSITY_DISTRIBUTION, COORDINATE_SYSTEM))
        print("      ln -sf %s_%s_%s_b_boundaries.dat b_boundaries.dat"
                                %(MODEL_NAME, DENSITY_DISTRIBUTION, COORDINATE_SYSTEM))
        print("      ln -sf %s_%s_%s_c_boundaries.dat c_boundaries.dat"
                                %(MODEL_NAME, DENSITY_DISTRIBUTION, COORDINATE_SYSTEM))
        print("      ln -sf %s_%s_%s_model.fits model.fits"
                                %(MODEL_NAME, DENSITY_DISTRIBUTION, COORDINATE_SYSTEM))
        print("      (in the 'input/grid/' directory)")
    else:
        os.system("ln -sf %s_%s_%s_a_boundaries.dat %sa_boundaries.dat"
                %(MODEL_NAME, DENSITY_DISTRIBUTION, COORDINATE_SYSTEM, PATH_INPUT_GRID))
        os.system("ln -sf %s_%s_%s_b_boundaries.dat %sb_boundaries.dat"
                %(MODEL_NAME, DENSITY_DISTRIBUTION, COORDINATE_SYSTEM, PATH_INPUT_GRID))
        os.system("ln -sf %s_%s_%s_c_boundaries.dat %sc_boundaries.dat"
                %(MODEL_NAME, DENSITY_DISTRIBUTION, COORDINATE_SYSTEM, PATH_INPUT_GRID))
        os.system("ln -sf %s_%s_%s_model.fits %smodel.fits"
                %(MODEL_NAME, DENSITY_DISTRIBUTION, COORDINATE_SYSTEM, PATH_INPUT_GRID))
        print("Symlinks created")

    print("Use COORDINATE_SYSTEM: '%s' and grid_type: '9' in the Mol3D input file"
                                %(COORDINATE_SYSTEM))
    print("")
    print("bye bye")

main()
plt.show()
