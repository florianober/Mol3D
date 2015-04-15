#!/usr/bin/python
# -*- coding: utf-8 -*
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
import mol3d_routines as l
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
            #~ w = np.argmax(map_in[0,:, :, :].sum(axis=1).sum(axis=1))
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
        # we merge neighbouring cells to provide a better result
        
        shrink_down = int(0.1*map_in.shape[2])
        
        # this is excluding the first pixel row/column
        # a bit dirty.....but it works for the moment
        I = hlp.shrink(map_in[0, w, 1:, 1:], shrink_down, shrink_down)
        
        Q = hlp.shrink(map_in[1, w, 1:, 1:], shrink_down, shrink_down)
        
        U = hlp.shrink(map_in[2, w, 1:, 1:], shrink_down, shrink_down)

        pol_deg, pol_angle = hlp.poldegang(I, Q, U)
        pol_angle += np.pi/2 # Quiver and astro definition of the angle

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
        
        plt.quiver(X, Y, pol_deg*np.cos(pol_angle),
                         pol_deg*np.sin(pol_angle), **props)
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
            plt.figure()
        else:
            plt.figure(fig_label+ ' 2')

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
        plt.tight_layout()

def make_continuum_all(path_results, pname):
    """
    prepare continuum data and pass them to the corresponding
    routines
    """

    project = l.mol3d(pname, path_results)
    r_ou = project.attr['r_ou']
    arcs = r_ou/project.attr['distance']
    extent = [-arcs, arcs, -arcs, arcs]
    
    # First temperature calculation sed/maps

    stokes = '_I'
    method_list = ['_bin','_mono', '_raytrace']
    for method in method_list:
        file_in = path_results + pname + '_continuum_sed' + method + stokes+'.dat'
        if os.path.isfile(file_in):
            # SED 
            sed_txt = np.genfromtxt(file_in, skiprows=0,
                                    filling_values="0")
            make_sed(sed_txt[:, 0], sed_txt[:, 1], fig_title=file_in)
            
            # Maps
            map_in = getattr(project, 'continuum_maps'+method)
            wlen = sed_txt[:, 0]
            make_continuum_maps(wlen, map_in, extent=extent)

    #  MC3D test cases
    #~ pname2 = 'th90i00a'
    #~ map_in = np.zeros((4, 100, 401, 401))
    #~ mc3d, vch2 =  hlp.read_mc3dmap2('/home/fober/mc3dtests/mc.v4.0.020/results/'+ pname2 +'.stokes_map.I')
    #~ map_in[0, 21, :, :] = mc3d[0,:, ::-1, ::-1]
    #~ 
    #~ mc3d, vch2 =  hlp.read_mc3dmap2('/home/fober/mc3dtests/mc.v4.0.020/results/'+ pname2 +'.stokes_map.Q')
    #~ map_in[1, 21, :, :] = mc3d[0,:, ::-1, ::-1]
    #~ 
    #~ mc3d, vch2 =  hlp.read_mc3dmap2('/home/fober/mc3dtests/mc.v4.0.020/results/'+ pname2 +'.stokes_map.U')
    #~ map_in[2, 21, :, :] = mc3d[0,:, ::-1, ::-1]
    #~ 
    #~ mc3d, vch2 =  hlp.read_mc3dmap2('/home/fober/mc3dtests/mc.v4.0.020/results/'+ pname2 +'.stokes_map.V')
    #~ map_in[3, 21, :, :] = mc3d[0,:, ::-1, ::-1]

if __name__ == "__main__":

    make_continuum_all(PATH_RESULTS, P_NAME)
    plt.show()
