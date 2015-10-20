#!/usr/bin/python
# -*- coding: utf-8 -*-
#------------------------------------------------------------------------------#
# This file is part of Mol3D.
#
#    Mol3D is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    Mol3D is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with Mol3D.  If not, see <http://www.gnu.org/licenses/>.
#
#    Diese Datei ist Teil von Mol3D.
#
#    Mol3D ist Freie Software: Sie können es unter den Bedingungen
#    der GNU General Public License, wie von der Free Software Foundation,
#    Version 3 der Lizenz oder (nach Ihrer Wahl) jeder späteren
#    veröffentlichten Version, weiterverbreiten und/oder modifizieren.
#
#    Mol3D wird in der Hoffnung, dass es nützlich sein wird, aber
#    OHNE JEDE GEWÄHRLEISTUNG, bereitgestellt; sogar ohne die implizite
#    Gewährleistung der MARKTFÄHIGKEIT oder EIGNUNG FÜR EINEN BESTIMMTEN ZWECK.
#    Siehe die GNU General Public License für weitere Details.
#
#    Sie sollten eine Kopie der GNU General Public License zusammen mit diesem
#    Programm erhalten haben. Wenn nicht, siehe <http://www.gnu.org/licenses/>.
#------------------------------------------------------------------------------#
"""
show the model set up
ver.: 0.5

date   : 20/10/2015
author : F. Ober
email  : fober@astrophysik.uni-kiel.de
"""

import numpy as np
import matplotlib.pyplot as plt
import os
import sys

from astropy.io import fits as pf
from matplotlib.colors import LogNorm

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

def make_cuts(path_results, p_name, axis='x'):
    # only x-axis is working/implemented at the moment
    #
    data_in = np.genfromtxt(path_results+ p_name + '_cut_'+axis +'.dat')

    n_dust = int((data_in.shape[1]-9)/2)

    plt.figure('Number density midplane')
    for i_dust in range(n_dust):
        plt.plot(data_in[:, 0], data_in[:, i_dust+1]*1e-6, '+-', label='Dust %d' %(i_dust+1))
    plt.xlabel('distance [AU]')
    plt.ylabel(r'number density cm$^{-3}$')
    plt.plot(data_in[:, 0], data_in[:, n_dust+2]*1e-6, '+-', label=r'H$_2$ %d')
    vmax = data_in[:, 1:n_dust+3].max()*1e-6 *2
    ind = (data_in[:, 1] > 1e-100)
    vmin = data_in[ind, 1:n_dust+3].min()*1e-6 * 0.5
    plt.ylim(vmin, vmax)
    plt.xlim(data_in[0, 0], data_in[-1, 0])
    plt.yscale('log')

    plt.legend()

    plt.figure('Temperatures midplane')
    for i_dust in range(n_dust):
        plt.plot(data_in[:, 0], data_in[:, n_dust + 5 + i_dust], '+-', label='Dust %d' %(i_dust+1))
    plt.plot(data_in[:, 0], data_in[:, 2*n_dust + 5], '+-', label='Gas')
    vmax = data_in[:, 6].max()*1.05
    vmin = 0
    plt.ylim(vmin, vmax)
    plt.xlim(data_in[0, 0], data_in[-1, 0])
    plt.xlabel('distance [AU]')
    plt.ylabel('temperature [K]')

    plt.legend()

def show_maps(path_results, p_name):
    """ search visualisation file for different planes """

    # show xz-plane
    present_plane(path_results+p_name, 'xz')

    # show xy-plane
    present_plane(path_results+p_name, 'xy')

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
        
        ext = file_path[-10: -8]+'-plane'
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
        xlab = 'distance [AU]'
        ylab = xlab
        
    #~ pic += 1e-250
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
                       pic[:, :, 7+2*n_dust]**2)/1000
    else:
        data = pic[:, :, 5+2*n_dust]
    CS = plt.contour(data, linewidths=1, colors='k', extent=m_range)
    plt.clabel(CS, inline=1, fmt='%2.1f', fontsize=10)
    plt.imshow(data, extent=m_range, origin='lower',
               interpolation='None', cmap=plt.cm.jet)
    plt.colorbar().set_label('Velocity [km/s]')

    #-------------------------------------------------
    #  molecule density

    data = pic[:, :, 0+n_dust] * 1e-6
    plt.figure('Molecule distribution, ' + ext)
    plt.title('Molecule distribution, ' + ext)
    plt.xlabel(xlab)
    plt.ylabel(ylab)
    vmax = np.max(data)
    vmin = vmax*1e-6
    cont = 10.0**np.linspace(np.log10(vmax), np.log10(vmin), 7)
    CS = plt.contour(data, levels=cont, linewidths=1, colors='k', extent=m_range)
    plt.clabel(CS, inline=1, fmt='%2.1e', fontsize=10)
    plt.imshow(data, extent=m_range, origin='lower',
               interpolation='None', cmap=plt.cm.jet,
               norm=LogNorm(vmin=vmin, vmax=vmax))

    plt.colorbar().set_label('molecule density [cm$^{-3}$]')

    #-------------------------------------------------
    #  dust density 
    for i in range(n_dust):
        data = pic[:, :, i] * 1e-6
        plt.figure('Dust %2.0d distribution, %s' %(i+1,ext))
        plt.title('Dust %2.0d distribution, %s' %(i+1,ext))
        plt.xlabel(xlab)
        plt.ylabel(ylab)
        vmax = np.max(data)
        vmin = vmax*1e-6
        cont = 10.0**np.linspace(np.log10(vmax), np.log10(vmin), 7)
    
        CS = plt.contour(data, levels=cont, linewidths=1, colors='k', extent=m_range)
        plt.clabel(CS, inline=1, fmt='%2.1e', fontsize=10)
        plt.imshow(data, extent=m_range, origin='lower',
                   interpolation='None', cmap=plt.cm.jet,
                   norm=LogNorm(vmin=vmin, vmax=vmax))
        plt.colorbar().set_label('dust density [cm$^{-3}$]')

    #-------------------------------------------------
    #  H2 density

    data = pic[:, :, 1+n_dust] * 1e-6
    plt.figure('H2 distribution, ' + ext)
    plt.title('H2 distribution, ' + ext)
    plt.xlabel(xlab)
    plt.ylabel(ylab)
    
    vmax = np.max(data)
    vmin = vmax*1e-6
    cont = 10.0**np.linspace(np.log10(vmax), np.log10(vmin), 7)
    
    CS = plt.contour(data, levels=cont, linewidths=1, colors='k', extent=m_range)
    plt.clabel(CS, inline=1, fmt='%2.1e', fontsize=10)
    plt.imshow(data, extent=m_range, origin='lower',
               interpolation='None', cmap=plt.cm.jet,
               norm=LogNorm(vmin=vmin, vmax=vmax))

    plt.colorbar().set_label('H2 number density [cm$^{-3}$]')
    #-------------------------------------------------
    #  N(Mol)/N(H)
    fig = plt.figure('N(mol)/N(H), ' + ext)
    
    H2  = pic[:, :, 1+n_dust]
    Mol = pic[:, :, 0+n_dust]
    data = np.zeros_like(H2)
    idn = np.where(H2 > 0)
    
    data[idn] = Mol[idn]/H2[idn]*0.5
    plt.xlabel(xlab)
    plt.ylabel(ylab)
    vmax = np.max(data)
    vmin = vmax*1e-6
    #~ m_range = [-100,100,-40,40]
    cont = 10.0**np.linspace(np.log10(vmax), np.log10(vmin), 7)
    CS = plt.contour(data[:,:], levels=cont, linewidths=1, colors='k', extent=m_range)
    plt.clabel(CS, inline=1, fmt='%2.1e', fontsize=10)


    plt.imshow(data[:,:], extent=m_range, origin='lower',
               interpolation='None', cmap=plt.cm.jet, norm=LogNorm(vmin=vmin, vmax=vmax))
    plt.colorbar().set_label(r'N(HCO$^{+}$)/N(H)')

if __name__ == "__main__":
    """ main routine """

    show_maps(PATH_RESULTS, P_NAME)
    make_cuts(PATH_RESULTS, P_NAME, axis='x')

    plt.show()
