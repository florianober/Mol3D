#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
show the model set up
ver.: 0.1

How to get started: 

date   : 22/02/2014
author : F. Ober 
email  : fober@astrophysik.uni-kiel.de
"""


import numpy as np
import matplotlib.pyplot as plt
import os
import sys
from scipy.interpolate import griddata    


home         = '../'

try:
    p_name = sys.argv[1]
except:
    p_name = 'example'
    
try:
    path_results = sys.argv[2]
    
except:
    #~ path_results = '../results/'
    path_results = '/data/fober/mol3dresults/'
    

def make_model(path_results,p_name):
    
    #open model file for the given p_name
    
    model = np.loadtxt(path_results+p_name+'_model.dat',skiprows=1)
    
    # tbd: assuming at the moment only one cell in phi direction
    x = -model[:,1]
    y = model[:,3]

    # gas Temperature
    
    z = model[:,10]
    create_plot(x,y,z,'Gas Temperature')
    
    # dust distribution
    
    z = model[:,4]
    create_plot(x,y,np.log10(z),'Dust number density distribution')
    
    # H2 distribution
    
    z = model[:,6]
    create_plot(x,y,np.log10(z),'H2 number density distribution')
    
    # molecule distribution
    
    z = model[:,5]
    create_plot(x,y,np.log10(z),'Molecule number density distribution')
    
    # molecule distribution in relation to H2
    
    z = model[:,5]/model[:,6]
    create_plot(x,y,z,'Molecule/H2 density distribution')

    # velocity distribution
    if model.shape[1] > 11:  # just to work with older versions
        z = model[:,11]/1000
        create_plot(x,y,z,'Velocity distribution')


def create_plot(x,y,z,name):
    
    xi = np.linspace(np.max(x),np.min(x),400)
    yi = np.linspace(np.max(y),np.min(y),400)
    
    zi = griddata((x, y), z, (xi[None,:], yi[:,None]), method='linear')
    
    plt.figure(name)
    #cont = [30,35,40,45,50]
    #~ print name.split()
    colors = 150
    bar = True
    if 'Temperature' in name:
        ext = ' [K]'
        cont = [30,35,40,45]
        
    elif 'density' in name:
        ext = ' lg [m^-3]'
        cont = 10
    elif 'Velocity' in name:
        ext = ' [km/s]'
        cont = [40,10,5,3]
        
    if 'Molecule/H2' in name:
        colors = 1
        bar = False

    if not 'Molecule/H2' in name:
        CS = plt.contour(xi,yi,zi,cont,linewidths=0.5,colors='k')
        plt.clabel(CS,inline=1,fmt='%2.1f', fontsize=10)
    CS = plt.contourf(xi,yi,zi,colors,cmap=plt.cm.jet)
    
    plt.xlim(np.min(x),np.max(x))
    plt.ylim(np.min(y),np.max(y))
    
    plt.title(name)
    plt.xlabel('d [AU]')
    plt.ylabel('d [AU]')
    if bar:
        plt.colorbar().set_label(name+ext)

if __name__ == "__main__":
    make_model(path_results,p_name)
    plt.show()
