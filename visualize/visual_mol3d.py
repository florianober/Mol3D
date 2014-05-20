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
    p_name = sys.argv[1]
except:
    p_name = 'temp2'
    
try:
    inp2 = sys.argv[2]
    
except:
    inp2 = ''

if inp2 == '':

    f = open('path_result.dat')
    inp2 = f.readline().split()[0]
    f.close()

path_results = inp2
    
home         = '../'
show_all = False
#~ show_all = True

if len(glob.glob(os.path.join(path_results,p_name+'*'))) < 1:
    print('results not found, maybe the path in "path_result.dat" is not correct')
    sys.exit()
    
def main():
    
    if show_all:
        # present the model

        sm.make_model(path_results,p_name)
        
        # present temperature in x midplane
        
        oneD(path_results+p_name+'_temp_x.dat',i=0)
        
    # present line spectrum
        
    oneD(path_results+p_name+'_velo_ch_mapsum.dat',i=1)
    
    # present velocity channel integrated map
    
    map_in = l.load_mol3d_map(path_results+p_name+'_velo_ch_mapint.dat')
    
    input_file = open(path_results+p_name+'_input_file.dat',"r")
    full = input_file.readlines()
    input_file.close()
    
    r_ou = 200  # standard value
    dist = 140  # standard value
    
    # search for key in input_file  

    attr = l.get_attr(p_name)
    
    
    
    # calculate arcseconds for a disk with given extension and distance
    arcs = attr['sf']*attr['r_ou'] / attr['distance']
    extent = [-arcs,arcs,-arcs,arcs]
    plt.figure(p_name)
    plt.imshow(map_in*1000,origin='lower',interpolation='None',extent=extent)
    #~ plt.clim(vmax=np.nanmax(map_in*1000)/2,vmin=None)
    #~ plt.clim(vmax=0.3,vmin=None)
    plt.colorbar().set_label('Flux [mJy/px]')
    
def oneD(file_path,i=0):

    pic = np.loadtxt(file_path)
    #print pic
    plt.figure(file_path)

    if i == 0:
            plt.xlabel('distance r [AU]')
            plt.ylabel('temperature [K]')
            title='midplane temperature distribution'
            maxi = 1
    elif i == 1:
            plt.ylim(0,1.1)
            plt.xlim(np.min(pic[:,0]),np.max(pic[:,0]))
            plt.xlabel('velocity [m/s]')
            plt.ylabel('normalized intensity')
            title='velocity spectrum'
            maxi = np.max(pic[:,1])
    else:
        title=''
        plt.title(title,fontsize=14)
        
    plt.plot(pic[:,0],pic[:,1]/maxi)

main()
plt.show()

print('bye bye')
