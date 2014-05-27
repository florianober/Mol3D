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
import little as li
import mol3d_routines as l
from mpl_toolkits.axes_grid1 import AxesGrid
from astropy.io import fits as pf

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
    print(pname)
    attr = l.get_attr(pname)
    map_in,vch = l.load_mol3d_fullvchmap(path_mol3d+pname+ch_map)
    print('channel map loaded')


    # plot spectrum of the hole map
    #~ no_pix = map_in.shape[2]*map_in.shape[1]
    #~ plt.figure(0)
    #~ plt.title(pname+ ' - full spectrum')
    #~ plt.xlabel('velocity [m/s]')
    #~ plt.ylabel('brightness temperature [K]')
    #~ plt.plot(vch, map_in.sum(axis=2).sum(axis=1)/no_pix)
    
    # plot velocity map

    # no_pix = map_in.shape[2]*map_in.shape[1]
    # velmap  = l.mol3d_intmap(map_in,vch)
    
    #~ velmap[200,195] = 0
    #~ l.plot_image(velmap,num=2,fig_title=pname,xlab="d [AU]",ylab="d [AU]",map_range=[-120,120,-120,120]
                    #~ ,color_loc=True)
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
    
    # Make spectra134, plots
    
    #~ fig = plt.figure(20)
    #~ grid = AxesGrid(fig, 111, # similar to subplot(132)
                    #~ nrows_ncols = (N, N),
                    #~ axes_pad = 0.0,
                    #~ share_all=True,
                    #~ label_mode = "1",
                    #~ cbar_location = "top",
                    #~ cbar_mode="single",
                    #~ aspect="auto"
                    #~ )
                    #~ 
    #~ grd_pos = [6,7,8,3,4,5,0,1,2] 
    #~ t = 0
    #~ for k in range(N):
        #~ for j in range(N):


            #~ spec = map_in[:,k*no_subpix:(k+1)*no_subpix,j*no_subpix:(j+1)*no_subpix].sum(axis=2).sum(axis=1)/no_subpix**2
            #~ grid[t].plot(vch,spec,'k',label='Pix no: '+str(i))
            #~ plt.figure(t+101)
            #fig, ax = plt.subplots()
            #~ im = grid[grd_pos[t]].plot(vch,spec,'k')
            #~ ext = [vch[0],-vch[0],0,np.max(spec)]
            #~ im = grid[grd_pos[t]].imshow(map_in.sum(axis=0)[k*no_subpix:(k+1)*no_subpix,j*no_subpix:(j+1)*no_subpix]*(vch[1]-vch[0])*1e-3,
                #~ origin='lower',interpolation='None',vmin=vmin,vmax=vmax,aspect="auto",extent=ext)
            
            #~ im = grid[t].imshow(map_in.sum(axis=0))
            #~ grid.cbar_axes[t].colorbar(im)

            #~ l.plot_image(map_in.sum(axis=0)[k*no_subpix:(k+1)*no_subpix,j*no_subpix:(j+1)*no_subpix],num=i)
            #print np.shape(map_in.sum(axis=0)[k*no_subpix:(k+1)*no_subpix,j*no_subpix:(j+1)*no_subpix])
            #~ t += 1
            #~ i += 1
    #~ grid.cbar_axes[0].colorbar(im)       
    #~ for cax in grid.cbar_axes:
        #~ cax.toggle_label(False)
            
    #grid.axes_llc.set_xticks([-2, 0, 2])
    #grid.axes_llc.set_yticks([-2, 0, 2])   
    

    
    N = 3
    #~ FWHM = li.ALMA_FWHM(attr['tr_lam'],16)
    FWHM = 0.09
    print(FWHM)
    #~ print('transition wavelength %2.2f mu'%(attr['tr_lam']*1e6))
    beam = np.array([4*FWHM,2*FWHM,FWHM,FWHM*0.5,FWHM*0.05])
    
    r_ou_new = attr['r_ou'] * attr['sf']
    arcs = r_ou_new/attr['distance']
    extent = [-arcs,arcs,-arcs,arcs]
    text_pos = [-0.9*arcs,-0.8*arcs]
    xx = np.linspace(-arcs,arcs,attr['n_bin_map']*2+1)
    dt = float(len(vch)-1)/(N*(N+2))
    
    step = 1 
    pic_vch = np.arange(step*15,step=step)+(len(vch)-1)/2-step*7
    
    try:
        hdulist = pf.open(path_mol3d+pname+'_final.fits')
        map_casa = hdulist[0].data
    except:
        hdulist = ''
        map_casa = ''

    
    #~ for i in range(len(beam)):
    for i in [4]:
        if i == 4:
            fig = plt.figure(pname+' overview map no real ALMA beam ')
        else:
            fig = plt.figure(pname+' overview map beam %1.1f x FWHM (%3.3f")' %(float(beam[i]/FWHM),FWHM))
        
        grid = AxesGrid(fig, 111, # similar to subplot(132)
                        nrows_ncols = (N, N+2),
                        axes_pad = 0.0,
                        share_all=True,
                        label_mode = "1",
                        cbar_location = "top",
                        cbar_mode="single",
                        aspect="auto"
                        )

        
        # make velocity channel overview map
        
        vmin = 0
        vmax =  np.max(l.conv(map_in[(len(vch)-1)/2,:,:],beam=beam[i],r_ou=r_ou_new,dist=attr['distance']))

        #~ for t in range(N*(N+2)):
        for t in range(15):
                #k = int(round(t*dt+dt/2))
                k = pic_vch[t]

                conv_map = l.conv(map_in[k,:,:],beam=beam[i],r_ou=r_ou_new,dist=attr['distance'])
                # plot x cut
                grid[t].plot(xx,conv_map[attr['n_bin_map'],:]/vmax*arcs*0.8,'k')
                # plot y cut
                grid[t].plot(xx,conv_map[:,attr['n_bin_map']]/vmax*arcs*0.8,'w')
                
                im = grid[t].imshow(conv_map,vmin=vmin,vmax=vmax,
                    origin='lower',interpolation='None',extent=extent,aspect="auto")
                grid[t].text(text_pos[0],text_pos[1], '%2.2f Km/s' %(vch[k]*1e-3), fontsize=10,
                          bbox={'facecolor':'white', 'alpha':0.4, 'pad':5})
        
        cbar = grid.cbar_axes[0].colorbar(im)
        #~ cbar.ax.set_label('Flux [mJy/beam]') 
        cbar.set_label_text('Flux [mJy/beam]') 

        grid.cbar_axes[0].colorbar(im)
        for cax in grid.cbar_axes:
            cax.toggle_label(False)
        
        grid.cbar_axes[0].toggle_label(True)
    
    
    #-----------------------------------compare with CASA ALMA simulation
    
    if hdulist != '':
        N_casa = map_casa.shape[2]
        xx = np.linspace(-arcs,arcs,N_casa)
        fig = plt.figure(pname+' overview map simulated with CASA ')

        grid = AxesGrid(fig, 111, # similar to subplot(132)
                        nrows_ncols = (N, N+2),
                        axes_pad = 0.0,
                        share_all=True,
                        label_mode = "1",
                        cbar_location = "top",
                        cbar_mode="single",
                        aspect="auto"
                        )

        
        # make velocity channel overview map
        
        vmin = 0
        vmax =  np.max(l.conv(map_casa[0,(len(vch)-1)/2,:,:],beam=beam[i],r_ou=r_ou_new,dist=attr['distance']))

        #~ for t in range(N*(N+2)):
        for t in range(15):
                #k = int(round(t*dt+dt/2))
                k = pic_vch[t]

                # plot x cut
                grid[t].plot(xx,map_casa[0,k,N_casa/2,:]/vmax*arcs*0.8,'k')
                # plot y cut
                grid[t].plot(xx,map_casa[0,k,:,N_casa/2]/vmax*arcs*0.8,'w')
                
                im = grid[t].imshow(map_casa[0,k,:,:],vmin=vmin,vmax=vmax,
                    origin='lower',interpolation='None',extent=extent,aspect="auto")
                grid[t].text(text_pos[0],text_pos[1], '%2.2f Km/s' %(vch[k]*1e-3), fontsize=10,
                          bbox={'facecolor':'white', 'alpha':0.4, 'pad':5})
        
        cbar = grid.cbar_axes[0].colorbar(im)
        #~ cbar.ax.set_label('Flux [mJy/beam]') 
        cbar.set_label_text('Flux [mJy/beam]') 

        grid.cbar_axes[0].colorbar(im)
        for cax in grid.cbar_axes:
            cax.toggle_label(False)
        
        grid.cbar_axes[0].toggle_label(True)

    #~ if hdulist == '':
        #~ save fits file for ALMA CASA
        #~ res = 2.0*li.as2deg(arcs)/map_in.shape[2]
        #~ print(res,r_ou_new)
        #~ li.create_CASA_fits(map_in,out_name=pname+'.fits',resolution=res)
    
    #~ try:
        #~ hdulist.close()
    #~ except:
        #~ pass
main()
plt.show()
