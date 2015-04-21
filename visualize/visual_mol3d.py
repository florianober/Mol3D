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
import make_spectrum as mkspec
import make_continuum  as cont
import glob

from mol3d_routines import mol3d

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

SHOW_ALL = True

if len(glob.glob(os.path.join(PATH_RESULTS, P_NAME + '*'))) < 1:
    print('results not found,' +
          'maybe the path in "path_result.dat" is not correct')
    sys.exit()

def main():
    """ Main visualisation routine  """

    project = mol3d(P_NAME, PATH_RESULTS)
    # get attributes
    attr = project.attr

    if not attr:
        print("ERROR, could not find results fot the requested project")
        sys.exit()

    if SHOW_ALL:
        # present the model
        sm.show_maps(PATH_RESULTS, P_NAME)

        # present temperature in x midplane for all dust species
        temp_paths = glob.glob(PATH_RESULTS+P_NAME+'_temp_x_??.dat')
        n_dust = len(temp_paths)
        for i in range(n_dust):
            one_dim(temp_paths[i], i=0)

    # present line spectrum, intensity map and velocity channel maps
    
    mkspec.make_spectra(PATH_RESULTS, P_NAME)

    # present continuum maps, sed's
    # tbd. Roman?
    cont.make_continuum_all(PATH_RESULTS, P_NAME)
    

def one_dim(file_path, i=0):
    """ 1D ploting routine """
    pic = np.loadtxt(file_path)
    plt.figure(file_path)

    if i == 0:
        plt.xlabel('Distance [AU]')
        plt.ylabel('Temperature [K]')
        plt.xlim(0, np.max(pic[:, 0]))
        title = 'midplane temperature distribution'
        maxi = 1
    elif i == 1:
        plt.ylim(0, 1.1*np.max(pic[:, 1]))
        plt.xlim(np.min(pic[:, 0]), np.max(pic[:, 0]))
        plt.xlabel('Velocity [m/s]')
        plt.ylabel('Normalized intensity')
        title = 'velocity spectrum'
        #~ maxi = np.max(pic[:, 1])
        maxi = 1
    else:
        title = ''
        plt.title(title, fontsize=14)

    plt.plot(pic[:, 0], pic[:, 1] / maxi)

main()
plt.show()

print('bye bye')
