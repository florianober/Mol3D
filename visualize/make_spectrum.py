#!/usr/bin/python
# -*- coding: utf-8 -*
""" script to create local spectra from mol3d velocity channel maps

date: 01/30/2014
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
from mpl_toolkits.axes_grid1 import AxesGrid, host_axes
import matplotlib as mpl
from astropy.io import fits as pf
import helper as hlp

if len(sys.argv) > 1:
    P_NAME = sys.argv[1]
else:
    P_NAME = 'example'

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

mpl.rcParams['font.size'] += 1

TEXT_SIZE = mpl.rcParams['font.size'] -6

    # make velocity channel overview map

def make_velo_int_plot(data, dvelo=1, cmap=plt.cm.nipy_spectral,
                       interpol='None', extent=[-1, 1, -1, 1],
                       label=r'flux [Jy/as$^2$ * km/s]',
                       conv=1, snr=-1):
    """ make velocity integrated intensity map """
    # show velocity integrated map
    fig = plt.figure('velocity integrated map')
    ax1 = host_axes([0.09, 0.11, 0.8, 0.8])
    extent_conv = [extent[0], 0.5* extent[0],0, 0.5*extent[1], extent[1]]
    dvelo /= 1000
    if extent != [-1, 1, -1, 1]:
        ax1.set_ylabel('["]')
        ax1.set_xlabel('["]')
    else:
        ax1.set_ylabel('distance')
        ax1.set_xlabel('distance')
    if snr == -1:
        im = ax1.imshow(data[:, :, :].sum(axis=0)*dvelo,
                   origin='lower', interpolation=interpol, cmap=cmap,
                   extent=extent)
        plt.colorbar(im).set_label(label)
    else:
        im = ax1.imshow(data[:, :, :].sum(axis=0)*dvelo/(snr),
                   origin='lower', interpolation=interpol, cmap=cmap,
                   extent=extent)
        plt.colorbar(im).set_label('SNR')
    
    if conv != 1:
        ax2 = ax1.twin()
        ax2.axis["right"].toggle(ticklabels=False)
        #~ ax2.set_xticklabels(ax1.get_xticks()*conv)
        ax2.set_xticks(extent_conv)
        ax2.set_xticklabels(np.round(ax2.get_xticks()*conv, 2))
        ax2.set_xlabel('distance [AU]')
    return fig
    
def make_velo_ch_plot(data, vch, N1=3, N2=5, snr='', cmap=plt.cm.nipy_spectral,
                      interpol='None', extent=[-1, 1, -1, 1], unit='["]',
                      label=r'flux density [Jy/as$^2$]', fig_label=''):
    """
    create velocity channel overview map
    """

    if len(vch) != data.shape[0] or data.shape[0] != N1*N2:
        print(len(vch), data.shape[0], N1*N2)
        print("ERROR, data and velocity array do not have the same shape")
        return
    if fig_label != '':
        fig = plt.figure(fig_label)
    else:
        fig = plt.figure()
    grid = AxesGrid(fig, 111, # similar to subplot(132)
                    nrows_ncols=(N1, N2),
                    axes_pad=0.0,
                    share_all=True,
                    label_mode="1",
                    cbar_location="right",
                    cbar_mode="single",
                    cbar_pad=0.02,
                    aspect="auto"
                   )
    vmin = 0
    if snr == '':
        sensitivity = 1
    else:
        sensitivity = snr
        #~ vmin = snr
    vmax = np.max(data[int(len(vch)/2.), :, :]/sensitivity)

    text_pos = [-0.9*extent[1], -0.8*extent[1]]
    for t in range(N1*N2):
        im = grid[t].imshow(data[t, :, :]/sensitivity, vmin=vmin, vmax=vmax,
                            cmap=cmap,
                            origin='lower', interpolation=interpol,
                            extent=extent, aspect="auto")

        # change the axes color from black to red
        grid[t].spines['right'].set_color('red')
        grid[t].spines['left'].set_color('red')
        grid[t].spines['top'].set_color('red')
        grid[t].spines['bottom'].set_color('red')

        # add velocity channel label to map
        grid[t].text(text_pos[0], text_pos[1], '%2.2f km/s' %(vch[t]*1e-3),
                     fontsize=TEXT_SIZE,
                     bbox={'facecolor':'white', 'alpha':0.7, 'pad':5})

    if snr == '':
        grid[0].cax.colorbar(im).set_label_text(label)
    else:
        grid[0].cax.colorbar(im).set_label_text(r'Flux density [$\sigma$]')

    grid[0].cax.toggle_label(True)
    if extent != [-1, 1, -1, 1]:
        grid.axes_llc.set_ylabel(unit)
        grid.axes_llc.set_xlabel(unit)
    else:
        grid.axes_llc.set_ylabel('distance')
        grid.axes_llc.set_xlabel('distance')

    #~ grid.cbar_axes[0].toggle_label(True)
    plt.setp(grid.axes_llc.xaxis.get_majorticklabels(), rotation=70)
    return fig

def make_spectra(path_results, pname, unit='PX'):
    """  present the spectral results of the project """

    project = mol3d(pname, path_results)

    map_in = project.return_velo_ch_map('AS')
    if map_in != []:
        
        vch = project.vch
        r_ou = project.attr['image_r_ou']
        arcs = project.attr['image_arcs']
        
        if unit == 'PX':
            extent = [-r_ou, r_ou, -r_ou, r_ou]
        else:
            extent = [-arcs, arcs, -arcs, arcs]

        # make line spectrum
        
        plt.figure('line spectrum')
        y_arr = map_in.sum(axis=1).sum(axis=1)
        plt.plot(vch, y_arr)
        plt.ylim(0, 1.1*np.max(y_arr))
        plt.xlim(np.min(vch), np.max(vch))
        plt.xlabel('velocity [m/s]')
        plt.ylabel('intensity [Jy]')

        # make intensity map
        map_in = project.return_velo_ch_map(unit)
        fig = make_velo_int_plot(map_in, project.attr['dvelo'],
                                 extent=[-arcs, arcs, -arcs, arcs],
                                 conv=project.attr['distance'])
        #~ fig.savefig(pname +'_intensity_map.pdf',bbox_inches='tight')

        # make velocity channel overview map
        mid = project.attr['i_vel_chan']
        incr = 2
        N1 = 3
        N2 = 5
        offset = int((N1*N2-1)/2.)
        fig = make_velo_ch_plot(map_in[mid-(incr*offset): mid + (incr*offset+1): incr],
                            vch[mid-(incr*offset): mid + (incr*offset+1): incr],
                            N1, N2, extent=extent, interpol='spline36')
        #~ fig.savefig(pname + '_velo_ch_map.pdf',bbox_inches='tight')

def make_spectra_fits(fits_file):
    fits = pf.open(fits_file)
    size = fits[0].data.shape
    if len(size) == 4: 
        map_in = fits[0].data[0, :, :, :]
        
    header =  fits[0].header
    fits.close()
    mid = (map_in.shape[0]-1)/2
    d_v = header['CDELT3']*hlp.c/header['CRVAL3']
    if 'ALTRVAL' in header:
        v_offset = header['ALTRVAL']
    else:
        v_offset = 0.0
    
    vch = np.linspace(v_offset, v_offset-d_v*map_in.shape[0], map_in.shape[0])
    arcs_x = hlp.deg2as(header['CDELT1']*((map_in.shape[1]-1)/2))
    arcs_y = header['CDELT2']*((map_in.shape[1]-1)/2)*3600
    extent = [-arcs_y, arcs_y, -arcs_x, arcs_x]
    incr = 2
    N1 = 5
    N2 = 5
    label = 'flux density ['+header['BUNIT']+']'
    offset = int((N1*N2-1)/2.)
    fig = make_velo_ch_plot(map_in[mid-(incr*offset): mid + (incr*offset+1): incr],
                            vch[mid-(incr*offset): mid + (incr*offset+1): incr],
                            N1, N2, extent=extent, interpol='spline36', label=label)
    # make line spectrum
    if 'BMAJ' in header:
        beam2px = hlp.beam2pixel_header(header)
    else:
        beam2px = 1

    #~ print (beam2px)
    plt.figure('line spectrum')
    y_arr = map_in.sum(axis=1).sum(axis=1)/beam2px
    plt.plot(vch/1000, y_arr)
    plt.ylim(0, 1.1*np.max(y_arr))
    plt.xlim(np.min(vch/1000), np.max(vch/1000))
    plt.xlabel('velocity [km/s]')
    plt.ylabel('intensity [Jy]')
    label = 'flux density ['+header['BUNIT']+' m/s]'
    # make intensity map
    fig = make_velo_int_plot(map_in, d_v, extent=extent, label=label)
    #~ fig.savefig(pname +'_intensity_map.pdf',bbox_inches='tight')
    print (y_arr.sum(axis=0) * d_v/1000 - np.median(y_arr) * d_v/1000*len(vch))
    #~ print (y_arr[0] * d_v/1000)
    
if __name__ == "__main__":
    #~ pass
    if '.fits' in P_NAME:
        make_spectra_fits('/home/fober/test.fits')
    else:
        make_spectra(PATH_RESULTS, P_NAME)
    plt.show()

