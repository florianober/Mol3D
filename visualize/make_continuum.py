#!/usr/bin/python
# -*- coding: utf-8 -*
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
""" script to visualize continuum maps and sed's

date: 26/02/2015
author: Florian Ober
email: fober@astrophysik.uni-kiel.de


"""
########################################################
### import some packages
#

import sys
import matplotlib.pyplot as plt
import numpy as np
from mol3d_routines import mol3d
from astropy.io import fits as pf
import helper as hlp
import os
from matplotlib.colors import LogNorm

if len(sys.argv) > 1:
    P_NAME = sys.argv[1]
else:
    P_NAME = 'example8'

if len(sys.argv) > 2:
    PATH_RESULTS = sys.argv[2]
else:
    try:
        FILE_IN = open('path_result.dat')
        PATH_RESULTS = FILE_IN.readline().split()[0]
        FILE_IN.close()
    except:
        PATH_RESULTS = '../results'
########################################################
### main


def make_sed(wavelength, data, 
             xlabel='wavelength [micron]', fig_title='Continuum SED',
             ylabel=r'flux density [Jy]', fig_label=''):
    """
    plot sed
    """
    if fig_label == '':
        plt.figure()
    else:
        plt.figure(fig_label)
    plt.title(fig_title)
    plt.plot(wavelength*1e6, data, 'd--')
    plt.xscale('log')
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)

def make_continuum_maps(wlen, map_in, fig_label='', cmap='Set1',
                       interpol='None', extent=[-1, 1, -1, 1]):
    """
    plot maps
    """
    stokes   = [0, 1, 2, 3]
    stokes_l = ['I', 'Q', 'U', 'V']
    if len(map_in.shape) == 4 and len(wlen) == map_in.shape[1]:
        if map_in.shape[1] > 1:
            w = np.argmax(np.nansum(np.nansum(map_in[0,:, :, :], axis=1), axis=1))
        else:
            w = 0

        if fig_label == '':
            plt.figure()
        else:
            plt.figure(fig_label+ ' 1')
        
        plt.title('wavelength %2.2f micron' %(wlen[w]*1e6))
        xlab = 'distance [as]'
        ylab = xlab
        vmin = 1e-12
        vmax= 1e-2
        #~ vmin = 0.01
        #~ vmax= 0.01
        
        I = map_in[0, w, :, :]
        Q = map_in[1, w, :, :]
        U = map_in[2, w, :, :]

        lin_pol = np.sqrt(Q**2 + U**2)
        plt.imshow(lin_pol, origin='lower',
                   interpolation=interpol, extent=extent,
                   cmap=cmap,
                   norm=LogNorm(vmin=vmin, vmax=vmax))
        plt.colorbar().set_label(r'linear polarisation [$\frac{Jy}{px}$]')
        plt.xlabel(xlab)
        plt.ylabel(ylab)

        # now we need to calculate the polarisation degree and angle
        # we merge neighbouring cells to provide better (visualisation) results
        I_interpolate = hlp.interpol_2d(map_in[0, w, :, :], map_in.shape[2]-1, map_in.shape[2]-1)

        shrink_down = int(0.05*I_interpolate.shape[0])

        I = hlp.shrink(I_interpolate, shrink_down, shrink_down)

        Q_interpolate = hlp.interpol_2d(map_in[1, w, :, :], map_in.shape[2]-1, map_in.shape[2]-1)
        Q = hlp.shrink(Q_interpolate, shrink_down, shrink_down)
        
        U_interpolate = hlp.interpol_2d(map_in[2, w, :, :], map_in.shape[2]-1, map_in.shape[2]-1)
        U = hlp.shrink(U_interpolate, shrink_down, shrink_down)

        pol_deg, pol_angle = hlp.poldegang(I, Q, U)
        pol_angle -= np.pi/2 # Quiver and astro definition of the angle

        x = np.linspace(extent[0], extent[1], I.shape[0])
        y = np.linspace(extent[0], extent[1], I.shape[0])
        X, Y = np.meshgrid(x, y)
        
        props = {'units'  : "dots",
                 #~ 'angles' : pol_angle*180/np.pi,
                 'scale'  : .01,
                 'width'  : 2,
                 'headwidth': 1e-24,
                 'headlength': 1e-24,
                 'headaxislength': 0,
                 'pivot': 'mid',
                 'color': 'black',
                }
        ind = (pol_deg > 0.001)
        plt.quiver(X[ind], Y[ind], pol_deg[ind]*np.cos(pol_angle[ind]),
                         pol_deg[ind]*np.sin(pol_angle[ind]), **props)
        #~ plt.quiver(X, Y, 0, pol_deg, **props)

        qk_posX = -0.1 # quiver key position 
        qk_posY = -0.1
        #~ normlen = np.amax(pol_deg)  # length of key vector
        normlen = 0.4  # length of key vector

        props = {'units' : "dots",
                  'angles': (0.),
                  'scale' : .01,
                  'width' : 2,
                  'headwidth': 1e-24,
                  'headlength': 1e-24,
                  'headaxislength': 0,
                  'pivot': 'mid',   # mid of vector on given (x,y)
                  }

        # Set quiver key
        polvec = plt.quiver(0., 0., 0., 0., color='black', **props) 
        
        qk = plt.quiverkey(polvec, qk_posX, qk_posY, normlen,
                          r'40%', labelpos='N',
                          fontproperties={'weight': 'bold', 'size':12},
                          color='black')
        plt.tight_layout()

        if fig_label == '':
            fig = plt.figure(num=None,figsize=(8,5.5))
        else:
            fig = plt.figure(fig_label+ ' 2',(1,17))

        plt.suptitle('wavelength %2.2f micron' %(wlen[w]*1e6), fontsize=14)

        # I
        st = 0
        vmin = 1e-12
        vmax = 1e-2
        plt.subplot(221+st)
        plt.imshow(map_in[st, w, :, :], origin='lower',
                  interpolation=interpol, extent=extent,
                  cmap=cmap,
                  norm=LogNorm(vmin=vmin, vmax=vmax, clip=True))

        plt.text(0.9*extent[0], 0.8*extent[1], 'Stokes %s' %stokes_l[st], fontsize=12,
                     bbox={'facecolor':'white', 'alpha':0.9, 'pad':5})
        plt.ylabel(ylab)
        plt.colorbar().set_label(r'flux density [$\frac{\rm Jy}{\rm px}$]')

        # Q
        st = 1
        vmin = -1.0e-11
        vmax =  1.0e-11
        plt.subplot(221+st)
        plt.imshow(map_in[st, w, :, :], origin='lower',
                  cmap='bwr',
                  interpolation=interpol, extent=extent, vmin=vmin, vmax=vmax)
        plt.text(0.9*extent[0], 0.8*extent[1], 'Stokes %s' %stokes_l[st], fontsize=12,
                     bbox={'facecolor':'white', 'alpha':0.9, 'pad':5})
        plt.colorbar()

        # U
        st = 2
        vmin = -1.0e-11
        vmax =  1.0e-11
        plt.subplot(221+st)
        plt.imshow(map_in[st, w, :, :], origin='lower',
                  cmap='bwr',
                  interpolation=interpol, extent=extent, vmin=vmin, vmax=vmax)
        plt.text(0.9*extent[0], 0.8*extent[1], 'Stokes %s' %stokes_l[st], fontsize=12,
                     bbox={'facecolor':'white', 'alpha':0.9, 'pad':5})
        plt.colorbar().set_label(r'flux density [$\frac{\rm Jy}{\rm px}$]')
        plt.ylabel(ylab)
        plt.xlabel(ylab)
        # V
        st = 3
        vmin = -1.0e-12
        vmax =  1.0e-12
        plt.subplot(221+st)
        plt.imshow(map_in[st, w, :, :], origin='lower',
                  cmap='bwr',
                  interpolation=interpol, extent=extent, vmin=vmin, vmax=vmax)
        plt.text(0.9*extent[0], 0.8*extent[1], 'Stokes %s' %stokes_l[st], fontsize=12,
                     bbox={'facecolor':'white', 'alpha':0.9, 'pad':5})
        plt.colorbar().set_label(r'flux density [$\frac{\rm Jy}{\rm px}$]')
        plt.xlabel(ylab)

def make_continuum_all(path_results, pname, unit='PX'):
    """
    prepare continuum data and pass them to the corresponding
    routines
    """

    project = mol3d(pname, path_results)
    r_ou = project.attr['image_r_ou']
    arcs = project.attr['image_arcs']
    if unit == 'PX':
        extent = [-r_ou, r_ou, -r_ou, r_ou]
    elif unit == 'AS':
        extent = [-arcs, arcs, -arcs, arcs]
    
    # First temperature calculation sed/maps

    stokes = '_I'
    method_list = ['_bin','_mono', '_raytrace']
    for method in method_list:
        file_in = path_results + pname + '_continuum_sed' + method + stokes+'.dat'
        if os.path.isfile(file_in):
            # SED 
            sed_txt = np.genfromtxt(file_in, skip_header=0,
                                    filling_values="0")
            make_sed(sed_txt[:, 0], sed_txt[:, 1], fig_title=file_in)
            
            # Maps
            map_in = getattr(project, 'return_continuum_maps'+method) (unit)
            wlen = sed_txt[:, 0]
            make_continuum_maps(wlen, map_in, extent=extent)

if __name__ == "__main__":

    make_continuum_all(PATH_RESULTS, P_NAME)
    plt.show()
