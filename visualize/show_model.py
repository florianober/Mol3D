#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
show the model set up
ver.: 0.2
routines marked with 'XX' belong to older versions and are not
used at the moment.

date   : 12/09/2014
author : F. Ober
email  : fober@astrophysik.uni-kiel.de
"""

import numpy as np
import matplotlib.pyplot as plt
import os
import sys
from scipy.interpolate import griddata
import matplotlib.patches as mpatches
from astropy.io import fits as pf
import mol3d_routines as m

if len(sys.argv) > 1:
    P_NAME = sys.argv[1]
else:
    P_NAME = 'example'

if len(sys.argv) > 2:
    PATH_RESULTS = sys.argv[2]
else:
    FILE_IN = open('path_result.dat')
    PATH_RESULTS = FILE_IN.readline().split()[0]
    FILE_IN.close()

def make_model(path_results, p_name):
    """ make plot from model file XX  """

    #open model file for the given p_name

    model = pf.open(path_results + p_name+'_model.fits')[0].data
    # tbd: assuming rotational symmetry
    x = (model[:, 0]**2+model[:, 1]**2)**0.5 * \
        np.sign(model[:, 0]) * np.sign(model[:, 1])
    y = model[:, 2]

    # gas Temperature

    z = model[:, 9]
    create_plot(x, y, z, 'Gas Temperature')

    # dust distribution

    z = model[:, 3]
    create_plot(x, y, np.log10(z*1e-6), 'Dust number density distribution')

    # H2 distribution
    z = model[:, 5]
    create_plot(x, y, np.log10(z*1e-6), 'H2 number density distribution')

    # molecule distribution
    z = model[:, 4]
    create_plot(x, y, np.log10(z*1e-6), 'Molecule number density distribution')

    # molecule distribution in relation to H2
    z = model[:, 4]/model[:, 5]
    create_plot(x, y, z, 'Molecule/H2 density distribution')

    #~ # velocity distribution
    z = model[:, 10] / 1000
    create_plot(x, y, z, 'Velocity distribution')

def create_plot(x, y, z, name):
    """ make a plot from arbitary spaced values using griddata XX """
    N = 501
    xi = np.linspace(np.max(x), np.min(x), N)
    yi = np.linspace(np.max(y), np.min(y), N)

    zi = griddata((x, y), z, (xi[None, :], yi[:, None]), method='linear')

    plt.figure(name)
    #cont = [30, 35, 40, 45, 50]
    bar = True
    colors = 151
    if 'Temperature' in name:
        ext = ' [K]'
        cont = [10, 20, 25, 30, 35, 40]
        #~ cont = 1
        #~ colors = 151
        colors = np.round(np.linspace(0, 100, 151))
    elif 'density' in name:
        ext = ' lg [cm^-3]'
        cont = np.round(np.linspace(-5, np.nanmax(zi), 10))
        colors = np.round(np.linspace(-12, np.nanmax(zi), 151))
        if 'H2' in name:
            ext = ' lg [cm^-3]'
            cont = np.round(np.linspace(-5, np.nanmax(zi), 10))
            cont = [5, 5.5, 6, 6.5, 7, 7.5, 8]
            colors = np.round(np.linspace(-12, np.nanmax(zi), 151))
    elif 'Velocity' in name:
        ext = ' [km/s]'
        #~ cont = [zi[,1.5, 2, 3, 5]
        cont = [1, 1.5, 2, 3, 5]
        #~ colors = np.round(np.linspace(0, 16, 151), 2)
    if 'Molecule/H2' in name:
        colors = 1
        bar = False

    if not 'Molecule/H2' in name:

        CS = plt.contour(xi, yi, zi, cont, linewidths=1, colors='k')
        plt.clabel(CS, inline=1, fmt='%2.1f', fontsize=10)

    CS = plt.contourf(xi, yi, zi, colors, cmap=plt.cm.jet)

    plt.xlim(np.min(x), np.max(x))
    plt.ylim(np.min(y), np.max(y))

    plt.title(name)
    plt.xlabel('d [AU]')
    plt.ylabel('d [AU]')
    if bar:
        plt.colorbar().set_label(name + ext)


def show_maps(path_results, p_name):
    """ search visualisation file for different planes """

    project = m.mol3d(p_name, path_results)
    # show xz-plane
    present_plane(path_results+p_name, 'xz')

    # show xy-plane
    #~ present_plane(path_results+p_name, 'xy')

    # show yz-plane
    #~ present_plane(path_results+p_name, 'yz')


def present_plane(file_in, plane):
    """ open visualisation file and make desired plots of that plane """

    if os.path.isfile(file_in + '_visual_' + plane + '.fits.gz'):
        file_path = file_in + '_visual_' + plane + '.fits.gz' 
        fits = pf.open(file_path)
        data_in = fits[0].data
        n_dust = fits[0].header['N_DUST']
        pix = fits[0].header['NAXIS1']
        r_ou = (pix-1)/2 * fits[0].header['CDELT1']
        xlab = 'distance [%s]' %(fits[0].header['CUNIT1'])
        ylab = xlab

        
        fits.close()

        pic = np.zeros((pix, pix, 8+2*n_dust))
        for i in range(8+2*n_dust):
            pic[:, :, i] = data_in[i, :, :]
        
        ext = file_path[-7: -5]+'-plane'
        m_range = [-r_ou,r_ou,-r_ou,r_ou]
       
    else:
        file_path = file_in + '_visual_' + plane + '.dat' 
        map_in = open(file_path)
        row = map_in.readline().split()
        map_size = int(row[0])
        row = map_in.readline()     # empty row
        n_dust = 1 # fixed value here

        for j in range(map_size):
            for k in range(map_size):
                row = map_in.readline().split()
                if j == 0 and j == 0:
                    len_pic = len(row)
                    pic = np.zeros((map_size, map_size, len_pic-2))

                if j == 0 and k == 0:
                    j_min = float(row[0])
                    k_min = float(row[1])
                elif j == map_size-1 and k == map_size-1:
                    j_max = float(row[0])
                    k_max = float(row[1])

                for l in range(2, len_pic):
                    if row[l] == 'NaN':
                        pic[k, j, l-2] = 0
                    else:
                        pic[k, j, l-2] = row[l]

        m_range = [k_min, k_max, j_min, j_max]
        ext = file_path[-6: -4]+'-plane'
        
        
    pic += 1e-250
    #-------------------------------------------------
    #  Dust Temp
    for i in range(n_dust):
        plt.figure('Dust %2.0d Temperature, %s' %(i+1, ext))
        plt.title('Dust %2.0d Temperature, %s' %(i+1, ext))
        plt.xlabel(xlab)
        plt.ylabel(ylab)
        cont = [10, 20, 25, 30, 35, 40] * 2
        data = pic[:, :, 4+n_dust+i]
        CS = plt.contour(data, cont, linewidths=1,
                         colors='k', extent=m_range)
        plt.clabel(CS, inline=1, fmt='%2.1f', fontsize=10)
        plt.imshow(data, extent=m_range, origin='lower',
                   interpolation='None', cmap=plt.cm.jet)
        plt.clim(0, 350)
        plt.colorbar().set_label('Temperature [K]')
    #-------------------------------------------------
    #  Gas Temp

    plt.figure('Gas Temperature, ' + ext)
    plt.title('Gas Temperature, ' + ext)
    plt.xlabel(xlab)
    plt.ylabel(ylab)
    cont = [10, 20, 25, 30, 35, 40] * 2
    data = pic[:, :, 4+2*n_dust]
    CS = plt.contour(data, cont, linewidths=1,
                     colors='k', extent=m_range)
    plt.clabel(CS, inline=1, fmt='%2.1f', fontsize=10)
    plt.imshow(data, extent=m_range, origin='lower',
               interpolation='None', cmap=plt.cm.jet)
    plt.clim(0, 350)
    plt.colorbar().set_label('Temperature [K]')

    #-------------------------------------------------
    #  Velocity
    plt.figure('abs(Velocity), ' + ext)
    plt.title('abs(Velocity), ' + ext)
    plt.xlabel(xlab)
    plt.ylabel(ylab)
    if pic.shape[2] > 9:
        data = np.sqrt(pic[:, :, 5+2*n_dust]**2 +
                       pic[:, :, 6+2*n_dust]**2 +
                       pic[:, :, 7+2*n_dust]**2)
    else:
        data = pic[:, :, 5+2*n_dust]
    CS = plt.contour(data, linewidths=1, colors='k', extent=m_range)
    plt.clabel(CS, inline=1, fmt='%2.1f', fontsize=10)
    plt.imshow(data, extent=m_range, origin='lower',
               interpolation='None', cmap=plt.cm.jet)
    plt.colorbar().set_label('Velocity [m/s]')

    #-------------------------------------------------
    #  molecule density

    data = np.log10(pic[:, :, 0+n_dust] * 1e-6)
    plt.figure('Molecule number density distribution, ' + ext)
    plt.title('Molecule number density distribution, ' + ext)
    plt.xlabel(xlab)
    plt.ylabel(ylab)
    cont = np.round(np.linspace(0, 7, 10))
    CS = plt.contour(data, cont, linewidths=1, colors='k', extent=m_range)
    plt.clabel(CS, inline=1, fmt='%2.1f', fontsize=10)
    plt.imshow(data, extent=m_range, origin='lower',
               interpolation='None', cmap=plt.cm.jet)
    plt.clim(2, 8)
    plt.colorbar().set_label('molecule density lg [cm^-3]')

    #-------------------------------------------------
    #  dust density 
    for i in range(n_dust):
        data = np.log10(pic[:, :, i] * 1e-6)
        plt.figure('Dust %2.0d number density distribution, %s' %(i+1,ext))
        plt.title('Dust %2.0d number density distribution, %s' %(i+1,ext))
        plt.xlabel(xlab)
        plt.ylabel(ylab)
        cont = np.round(np.linspace(-5, 0, 10))

        CS = plt.contour(data, cont, linewidths=1, colors='k', extent=m_range)
        plt.clabel(CS, inline=1, fmt='%2.1f', fontsize=10)
        plt.imshow(data, extent=m_range, origin='lower',
                   interpolation='None', cmap=plt.cm.jet)
        plt.clim(-5, 0)
        plt.colorbar().set_label('dust density lg [cm^-3]')

    #-------------------------------------------------
    #  H2 density

    data = np.log10(pic[:, :, 1+n_dust] * 1e-6)
    plt.figure('H2 number density distribution, ' + ext)
    #~ plt.title('H2 number density distribution, ' + ext)
    plt.xlabel(xlab)
    plt.ylabel(ylab)
    #~ cont = np.round(np.linspace(5.0, 8.0 ,5.), 1)
    cont = [8.0, 7.5, 7.0, 6.0, 5.5]

    CS = plt.contour(data, cont, linewidths=1, colors='k', extent=m_range)
    plt.clabel(CS, inline=1, fmt='%2.1f', fontsize=10)
    plt.imshow(data, extent=m_range, origin='lower',
               interpolation='None', cmap=plt.cm.jet)
    #~ plt.annotate('gap', xy=(80, 10), xycoords='data', color='white',
                 #~ xytext=(0, 80), textcoords='data',
                 #~ arrowprops=dict(width=3, facecolor='white', shrink=0.05),
                 #~ fontsize=24,
                 #~ horizontalalignment='center', verticalalignment='top',)
    #~ plt.annotate('gap', xy=(-80, 10), xycoords='data',
                 #~ xytext=(0, 80), textcoords='data', color='white',
                 #~ arrowprops=dict(width=3, facecolor='white', shrink=0.05),
                 #~ fontsize=23,
                 #~ horizontalalignment='center', verticalalignment='top',)
    plt.clim(5, 9)
    plt.colorbar().set_label('H2 number density lg [cm^-3]')

if __name__ == "__main__":
    """ main routine """
    
    show_maps(PATH_RESULTS, P_NAME)

    plt.show()
