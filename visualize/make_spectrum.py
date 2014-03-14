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
import os
import time
import math
import sys
import matplotlib.pyplot as plt
import numpy as np
import little as l
from mpl_toolkits.axes_grid1 import AxesGrid

try:
    pname = sys.argv[1]

except:
    pname = 'CO3to2nogapana'
    

path_mol3d  = '/data/fober/mol3dresults/'
ch_map      = '_velo_ch_map.dat'
int_map      = '_velo_ch_mapint.dat'

########################################################
### main


    
def main():
    print pname
    
    map_in,vch = l.load_mol3d_fullvchmap(path_mol3d+pname+ch_map)
    print  'channel map loaded'
    
    # plot spectrum of the hole map
    no_pix = map_in.shape[2]*map_in.shape[1]
    plt.figure(0)
    plt.title(pname+ ' - full spectrum')
    plt.xlabel('velocity [m/s]')
    plt.ylabel('brightness temperature [K]')
    plt.plot(vch, map_in.sum(axis=2).sum(axis=1)/no_pix)
    
    # plot velocity map

    no_pix = map_in.shape[2]*map_in.shape[1]
    velmap  = l.mol3d_intmap(map_in,vch)
    
    #~ velmap[200,195] = 0
    l.plot_image(velmap,num=2,fig_title=pname,xlab="d [AU]",ylab="d [AU]",map_range=[-120,120,-120,120]
                    ,color_loc=True)
    #~ plt.figure(8)
    
    #~ plt.plot(vch,map_in[:,200,206])
    # plot local spectrum
    N = 3  # Number of subpixel in one direction (full subpixel = N**3)
    
    no_subpix = float(map_in.shape[2])/N
    # tbd solve if no_subpix is not an integer
    
    i = 1
    #~ plt.figure(1)
    #~ plt.xlabel('velocity [m/s]')
    #~ plt.ylabel('brightness temperature [K]')
    #~ plt.title(pname + ' - local spectra')
    #~ map_in[:,0:267,0:267] = 0.0
    vmax =     np.max(map_in.sum(axis=0)*(vch[1]-vch[0])*1e-3)
    vmin =     0
    
    # Make spectra plots
    
    fig = plt.figure(20)
    grid = AxesGrid(fig, 111, # similar to subplot(132)
                    nrows_ncols = (N, N),
                    axes_pad = 0.0,
                    share_all=True,
                    label_mode = "1",
                    cbar_location = "top",
                    cbar_mode="single",
                    aspect="auto"
                    )
                    
    grd_pos = [6,7,8,3,4,5,0,1,2] 
    t = 0
    for k in range(N):
        for j in range(N):


            spec = map_in[:,k*no_subpix:(k+1)*no_subpix,j*no_subpix:(j+1)*no_subpix].sum(axis=2).sum(axis=1)/no_subpix**2
            #~ grid[t].plot(vch,spec,'k',label='Pix no: '+str(i))
            #~ plt.figure(t+101)
            #fig, ax = plt.subplots()
            #~ im = grid[grd_pos[t]].plot(vch,spec,'k')
            ext = [vch[0],-vch[0],0,np.max(spec)]
            im = grid[grd_pos[t]].imshow(map_in.sum(axis=0)[k*no_subpix:(k+1)*no_subpix,j*no_subpix:(j+1)*no_subpix]*(vch[1]-vch[0])*1e-3,
                origin='lower',interpolation='None',vmin=vmin,vmax=vmax,aspect="auto",extent=ext)
            
            #~ im = grid[t].imshow(map_in.sum(axis=0))
            #~ grid.cbar_axes[t].colorbar(im)

            #~ l.plot_image(map_in.sum(axis=0)[k*no_subpix:(k+1)*no_subpix,j*no_subpix:(j+1)*no_subpix],num=i)
            #print np.shape(map_in.sum(axis=0)[k*no_subpix:(k+1)*no_subpix,j*no_subpix:(j+1)*no_subpix])
            t += 1
            i += 1
    #~ grid.cbar_axes[0].colorbar(im)       
    #~ for cax in grid.cbar_axes:
        #~ cax.toggle_label(False)
            
    #grid.axes_llc.set_xticks([-2, 0, 2])
    #grid.axes_llc.set_yticks([-2, 0, 2])   
    
    N = 5
    fig = plt.figure(22)
    grid = AxesGrid(fig, 111, # similar to subplot(132)
                    nrows_ncols = (N, N),
                    axes_pad = 0.0,
                    share_all=True,
                    label_mode = "1",
                    cbar_location = "top",
                    cbar_mode="single",
                    aspect="auto"
                    )
    #~ print dir(grid[0])
    t = 0.0
    dt = float(len(vch)-1)/(N**2)
    # make velocity channel overview plot
    i = 0
    vmax = np.max(map_in)
    for t in range(N**2):

            #~ spec = map_in[t,k*no_subpix:(k+1)*no_subpix,j*no_subpix:(j+1)*no_subpix].sum(axis=2).sum(axis=1)/no_subpix**2
            #~ grid[t].plot(vch,spec,'k',label='Pix no: '+str(i))
            #~ im = grid[t].plot(vch,spec,'k')
            k = int(round(t*dt))
            #~ print k,t
            im = grid[t].imshow(map_in[k,:,:],
                origin='lower',interpolation='None',vmin=vmin,vmax=vmax,extent=[-120,120,-120,120])
            #~ im = grid[t].imshow(map_in.sum(axis=0))
            #~ grid.cbar_axes[t].colorbar(im)

            grid[t].text(-100,90, '%2.2f Km/s' %(vch[k]*1e-3), fontsize=12,
                      bbox={'facecolor':'white', 'alpha':0.7, 'pad':5}) 
            #~ print int(t*(dt+1))
            #~ l.plot_image(map_in.sum(axis=0)[k*no_subpix:(k+1)*no_subpix,j*no_subpix:(j+1)*no_subpix],num=i)
            #print np.shape(map_in.sum(axis=0)[k*no_subpix:(k+1)*no_subpix,j*no_subpix:(j+1)*no_subpix])
            #~ i += 1
    grid.cbar_axes[0].colorbar(im)       
    for cax in grid.cbar_axes:
        cax.toggle_label(False)
            
    #~ plt.legend()
    #~ print np.shape(map_in.sum(axis=0)[:100,:100])
    #~ l.plot_image(map_in.sum(axis=0)[:100,:100],num=2)
    #~ l.plot_image(map_in.sum(axis=0)[100:,:100],num=3)
    #~ l.plot_image(map_in.sum(axis=0)[:100,101:],num=4)
    #~ l.plot_image(map_in.sum(axis=0)[101:,101:],num=5)
    #~ 
    


    #~ l.plot_image(map_in.sum(axis=0)[:,:]*(vch[1]-vch[0])*1e-3,num=15,fig_title=pname,xlab="d [AU]",ylab="d [AU]",map_range=[-120,120,-120,120] )
    #~ plt.colorbar()

    plt.show()
    


main()