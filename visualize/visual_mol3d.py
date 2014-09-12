#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
easy visualisation tool for mol3d.
ver.: 0.1

How to get started: ./visual_mol3d.py example

date   : 2014-02-24
author : F. Ober
email  : fober@astrophysik.uni-kiel.de
"""


import numpy as np
import matplotlib.pyplot as plt
import os
import glob
import sys
import show_model as sm

import mol3d_routines as l

try:
    P_NAME = sys.argv[1]
except:
    P_NAME = 'temp2'

if len(sys.argv) > 2:
    PATH_RESULTS = sys.argv[2]

else:
    FILE_IN = open('path_result.dat')
    PATH_RESULTS = FILE_IN.readline().split()[0]
    FILE_IN.close()

SHOW_ALL = True


if len(glob.glob(os.path.join(PATH_RESULTS, P_NAME + '*'))) < 1:
    print('results not found,' +
          'maybe the path in "path_result.dat" is not correct')
    sys.exit()

def main():
    """ Main visualisation routine  """

    project = l.mol3d(P_NAME, PATH_RESULTS)
    # get attributes
    attr = project.attr

    if not attr:
        sys.exit()

    if SHOW_ALL:
        # present the model
        sm.show_maps(PATH_RESULTS, P_NAME)

        # present temperature in x midplane

        one_dim(PATH_RESULTS+P_NAME+'_temp_x.dat', i=0)
    # present line spectrum

    one_dim(PATH_RESULTS+P_NAME+'_spectrum.dat', i=1)

    # present velocity channel integrated map
    # tbd

    # calculate arcseconds for a disk with given extension and distance
    #~ arcs = attr['sf']*attr['r_ou'] / attr['distance'] / attr['zoom_map']
    arcs = attr['r_ou'] / attr['distance'] / attr['zoom_map']
    extent = [-arcs, arcs, -arcs, arcs]
    plt.figure(P_NAME)
    map_in = np.sum(project.velo_ch_map, axis = 0) * attr['dvelo']
    plt.imshow(map_in*1000, origin='lower',
               interpolation='None', extent=extent)
    #~ plt.clim(vmax=np.nanmax(map_in*1000)/2,vmin=None)
    #~ plt.clim(vmax=0.3,vmin=None)
    plt.colorbar().set_label('Flux [mJy/px]')

def one_dim(file_path, i=0):
    """ 1D ploting routine """
    pic = np.loadtxt(file_path)
    #print pic
    plt.figure(file_path)

    if i == 0:
        plt.xlabel('distance r [AU]')
        plt.ylabel('temperature [K]')
        title = 'midplane temperature distribution'
        maxi = 1
    elif i == 1:
        plt.ylim(0, 1.1)
        plt.xlim(np.min(pic[:, 0]), np.max(pic[:, 0]))
        plt.xlabel('velocity [m/s]')
        plt.ylabel('normalized intensity')
        title = 'velocity spectrum'
        maxi = np.max(pic[:, 1])
    else:
        title = ''
        plt.title(title, fontsize=14)

    plt.plot(pic[:, 0], pic[:, 1] / maxi)

main()
plt.show()

print('bye bye')
