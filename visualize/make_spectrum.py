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
import mol3d_routines as l
from mpl_toolkits.axes_grid1 import AxesGrid
import matplotlib as mpl

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

########################################################
### main

mpl.rcParams['font.size'] += 1

TEXT_SIZE = mpl.rcParams['font.size'] -2

    # make velocity channel overview map

def make_velo_int_plot(data, dvelo=1, cmap=plt.cm.nipy_spectral,
                       interpol='None', extent=[-1, 1, -1, 1],
                       label='flux [Jy/px * m/s]'):
    """ make velocity integrated intensity map """

    # show velocity integrated map
    fig = plt.figure('velocity integrated map')
    plt.imshow(data[:, :, :].sum(axis=0)*dvelo,
               origin='lower', interpolation=interpol, cmap=cmap,
               extent=extent)

    plt.colorbar().set_label(label)

    if extent != [-1, 1, -1, 1]:
        plt.ylabel('angular distance ["]')
        plt.xlabel('angular distance ["]')
    else:
        plt.ylabel('distance')
        plt.xlabel('distance')

def make_velo_ch_plot(data, vch, N1=3, N2=5, snr=1, cmap=plt.cm.nipy_spectral,
                      interpol='None', extent=[-1, 1, -1, 1],
                      label='flux [mJy/px]'):
    """
    create velocity channel overview map
    """

    if len(vch) != data.shape[0] or data.shape[0] != N1*N2:
        print(len(vch), data.shape[0], N1*N2)
        print("ERROR, data and velocity array do not have the same shape")
        return

    fig = plt.figure('velocity channel overview map')
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
    vmax = np.max(data[int(len(vch)/2.), :, :]/snr*1000)

    text_pos = [-0.9*extent[1], -0.8*extent[1]]
    for t in range(N1*N2):
        im = grid[t].imshow(data[t, :, :]/snr*1000, vmin=vmin, vmax=vmax,
                            cmap=cmap,
                            origin='lower', interpolation=interpol,
                            extent=extent, aspect="auto")
        grid[t].spines['right'].set_color('red')
        grid[t].spines['left'].set_color('red')
        grid[t].spines['top'].set_color('red')
        grid[t].spines['bottom'].set_color('red')
        grid[t].text(text_pos[0], text_pos[1], '%2.2f km/s' %(vch[t]*1e-3),
                     fontsize=TEXT_SIZE,
                     bbox={'facecolor':'white', 'alpha':0.7, 'pad':5})

    if snr == 1:
        grid[0].cax.colorbar(im).set_label_text(label)
    else:
        grid[0].cax.colorbar(im).set_label_text('SNR')

    grid[0].cax.toggle_label(True)
    if extent != [-1, 1, -1, 1]:
        grid.axes_llc.set_ylabel('angular distance ["]')
        grid.axes_llc.set_xlabel('angular distance ["]')
    else:
        grid.axes_llc.set_ylabel('distance')
        grid.axes_llc.set_xlabel('distance')

    #~ grid.cbar_axes[0].toggle_label(True)
    plt.setp(grid.axes_llc.xaxis.get_majorticklabels(), rotation=70)

    #~ fig.savefig(pname+'_syn_velo_ch_map' +str(i)+'.pdf',bbox_inches='tight')

def make_spectra(path_results, pname):
    """  present the spectral results of the project """

    project = l.mol3d(pname, path_results)

    map_in = project.velo_ch_map
    vch = project.vch
    r_ou = project.attr['r_ou']# * project.attr['sf']
    arcs = r_ou/project.attr['distance']
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
    make_velo_int_plot(map_in, project.attr['dvelo'], extent=extent)

    # make velocity channel overview map
    mid = project.attr['i_vel_chan']
    incr = 2
    N1 = 3
    N2 = 5
    offset = int((N1*N2-1)/2.)
    make_velo_ch_plot(map_in[mid-(incr*offset): mid + (incr*offset+1): incr],
                      vch[mid-(incr*offset): mid + (incr*offset+1): incr],
                      N1, N2, extent=extent, interpol='spline36')

if __name__ == "__main__":
    make_spectra(PATH_RESULTS, P_NAME)
    plt.show()
