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
import sys
import show_model as sm

import little as l

try:
    p_name = sys.argv[1]
except:
    p_name = 'temp2'
    
try:
    inp2 = sys.argv[2]
    
except:
    inp2 = '/data/fober/mol3dresults/'

path_results = inp2
    
home         = '../'
show_all = False

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
    l.plot_image(map_in*1000,num=path_results+p_name+'_velo_ch_mapint.dat')
    plt.colorbar().set_label('Flux [mJy]')
    
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
