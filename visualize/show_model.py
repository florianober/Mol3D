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
    
    # tbd: assuming rotational symmetry
    x = (model[:,1]**2+model[:,2]**2)**0.5*np.sign(model[:,1])*np.sign(model[:,2])
    y = model[:,3]

    # gas Temperature
    
    #~ z = model[:,10]
    #~ create_plot(x,y,z,'Gas Temperature')
    
    # dust distribution
    
    #~ z = model[:,4]
    #~ create_plot(x,y,np.log10(z*1e-6),'Dust number density distribution')
    
    # H2 distribution
    #~ 
    z = model[:,6]
    create_plot(x,y,np.log10(z*1e-6),'H2 number density distribution')
    
    # molecule distribution
    
    z = model[:,5]
    create_plot(x,y,np.log10(z*1e-6),'Molecule number density distribution')
    
    # molecule distribution in relation to H2
    #~ 
    #~ z = model[:,5]/model[:,6]
    #~ create_plot(x,y,z,'Molecule/H2 density distribution')

    # velocity distribution
    #~ if model.shape[1] > 11:  # just to work with older versions
        #~ z = model[:,11]/1000
        #~ create_plot(x,y,z,'Velocity distribution')


def create_plot(x,y,z,name):
    
    N = 501
    xi = np.linspace(np.max(x),np.min(x),N)
    yi = np.linspace(np.max(y),np.min(y),N)
    
    zi = griddata((x, y), z, (xi[None,:], yi[:,None]), method='linear')
    
    plt.figure(name)
    #cont = [30,35,40,45,50]
    bar = True
    colors = 151
    if 'Temperature' in name:
        ext = ' [K]'
        cont = [10,20,25,30,35,40]
        #~ colors = 151
        colors = np.round(np.linspace(0,100,151))
    elif 'density' in name:
        ext = ' lg [cm^-3]'
        cont = np.round(np.linspace(-5,np.nanmax(zi),10))
        colors = np.round(np.linspace(-12,np.nanmax(zi),151))
        if 'H2' in name:
            ext = ' lg [cm^-3]'
            cont = np.round(np.linspace(-5,np.nanmax(zi),10))
            cont = [5,5.5,6,6.5,7,7.5,8]
            colors = np.round(np.linspace(-12,np.nanmax(zi),151))
    elif 'Velocity' in name:
        ext = ' [km/s]'
        #~ cont = [zi[,1.5,2,3,5]
        cont = [1,1.5,2,3,5]
        colors = np.round(np.linspace(0,16,151),2)
    if 'Molecule/H2' in name:
        colors = 1
        bar = False

    if not 'Molecule/H2' in name:
        #~ pass
        CS = plt.contour(xi,yi,zi,cont,linewidths=1,colors='k')
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
