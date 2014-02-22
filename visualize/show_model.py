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

#~ pyfits_log = False
#~ import little as l

#~ try:
    #~ import pyfits as pf
    #~ pyfits_log = True
#~ except:
    #~ print ''
    #~ print 'could not find pyfits. Thus it is not possible to create/manipulate fits files.'
    


path_results = '../results/'
home         = '../'


try:
    p_name = sys.argv[1]
    inp = False 
except:
    p_name = 'example'
    inp = True

def main():
    
    #open model file for the given p_name
    
    model = np.loadtxt(path_results+p_name+'_model.dat',skiprows=1)
    
    # tbd: assuming at the moment only one cell in phi direction
    x = -model[:,1]
    y = model[:,3]

    # gas Temperature
    
    z = model[:,10]
    create_plot(x,y,z,'Gas Temperature [K]')
    
    # dust distribution
    
    z = model[:,4]
    create_plot(x,y,np.log10(z),'Dust number density distribution   ')
    
    # H2 distribution
    
    z = model[:,6]
    create_plot(x,y,np.log10(z),'H2 number density distribution   ')
    
    # molecule distribution
    
    z = model[:,5]
    create_plot(x,y,np.log10(z),'molecule number density distribution   ')
    
    # molecule distribution in relation to H2
    
    z = model[:,5]/model[:,6]
    create_plot(x,y,z,'mol/H2 density distribution   ')

def create_plot(x,y,z,name):
    
    xi = np.linspace(np.max(x),np.min(x),400)
    yi = np.linspace(np.max(y),np.min(y),400)
    
    zi = griddata((x, y), z, (xi[None,:], yi[:,None]), method='linear')
    
    plt.figure(name[:-3])
    #cont = [30,35,40,45,50]
    CS = plt.contour(xi,yi,zi,15,linewidths=0.5,colors='k')
    plt.clabel(CS,inline=1,fmt='%2.1f', fontsize=10)
    CS = plt.contourf(xi,yi,zi,100,cmap=plt.cm.jet)
    
    plt.xlim(np.min(x),np.max(x))
    plt.ylim(np.min(y),np.max(y))
    
    plt.title(name[:-3])
    plt.xlabel('d [AU]')
    plt.ylabel('d [AU]')
    #~ plt.clabel('Temperature [K]')
    plt.colorbar().set_label(name)


main()
plt.show()
