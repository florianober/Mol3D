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
from mpl_toolkits.axes_grid1 import AxesGrid
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

def make_continuum_maps():
    """
    plot maps
    """
    pass

def make_continuum_all(path_results, pname):
    """
    prepare continuum data and pass them to the corresponding
    routines
    
    """
    project = l.mol3d(pname, path_results)
    
    # SED 
    # load from file
    stokes = '_I'
    method = '_bin'
    sed_txt = np.genfromtxt(path_results+pname+'_continuum_sed' +
                         method + stokes+'.dat', skiprows=0,
                         filling_values="0")

    make_sed(sed_txt[:, 0], sed_txt[:, 1])

    # Maps

    # TbD!
    make_continuum_maps()

if __name__ == "__main__":

    make_continuum_all(PATH_RESULTS, P_NAME)
    plt.show()

