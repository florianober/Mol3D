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
easy visualisation tool for mol3d.
ver.: 1.0

How to get started: ./visual_mol3d.py example

date   : Jul 23 2015
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
