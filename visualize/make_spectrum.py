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
import helper as hlp
import mol3d_routines as l
from mpl_toolkits.axes_grid1 import AxesGrid
from astropy.io import fits as pf
from mpl_toolkits.axes_grid1.anchored_artists import AnchoredEllipse
import matplotlib as mpl

try:
    pname = sys.argv[1]

except:
    pname = 'CO3to2nogapana'
    

#~ pname = 'k05m05st01i10CO3'
pname = 'k05m04st01i10CO3'
conf  = 20


path_mol3d  = '/data/fober/mol3dresults/'
path_fits_results  = '/mnt/kronos/data3/fober/casa_out/casa_therm/'
#~ path_mol3d  = '../results/'
ch_map      = '_velo_ch_map.dat'
int_map      = '_velo_ch_mapint.dat'

########################################################
### main

mpl.rcParams['font.size']       += 4

def main():

    project  = l.mol3d('s1_'+pname)
    #~ map_in,vch = l.load_mol3d_fullvchmap(path_mol3d+pname+ch_map)
    
    map_s1  = pf.open(path_fits_results + 's1_' + pname +'_velo_ch_map_alma%2.2d.noisy.fits' %(conf))[0].data[0,:,:,:]
    map_sp1 = pf.open(path_fits_results + 'sp1_'+ pname +'_velo_ch_map_alma%2.2d.noisy.fits' %(conf))[0].data[0,:,:,:]
    
    #~ map_s1  = pf.open('/mnt/kronos/data2/fober/fits/s1_k05m05st01i10CO3_velo_ch_map.fits')[0].data[0,:,:,:]
    #~ map_sp1 = pf.open('/mnt/kronos/data2/fober/fits/sp1_k05m05st01i10CO3_velo_ch_map.fits')[0].data[0,:,:,:]

    #~ map_s1  = pf.open(path_fits_results + 's1_' + pname +'_velo_ch_map_alma%2.2d.fits' %(conf))[0].data[0,:,:,:]
    #~ map_sp1 = pf.open(path_fits_results + 'sp1_'+ pname +'_velo_ch_map_alma%2.2d.fits' %(conf))[0].data[0,:,:,:]
    
    
    
    
    map_in  = np.zeros((3,map_s1.shape[0],map_s1.shape[1],map_s1.shape[2]))
    
    
    
    header = pf.open(path_fits_results + 's1_' + pname +'_velo_ch_map_alma%2.2d.fits' %(conf))[0].header


    beam2px = (abs(header['BMAJ'])*abs(header['BMIN'])*np.pi)/(abs(header['CDELT1'])*abs(header['CDELT2']))
    #~ print(header['BMAJ'], header['BMIN'])
    
    #~ sys.exit()
    
    map_in[0,:,:,:] = map_s1
    map_in[1,:,:,:] = map_sp1
    map_in[2,:,:,:] = np.abs(map_s1-map_sp1)
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
    #~ vmax =     np.max(map_in.sum(axis=0)*(vch[1]-vch[0])*1e-3)
    #~ vmin =     0
    
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

    #~ print('transition wavelength %2.2f mu'%(attr['tr_lam']*1e6))
    
    r_ou_new = project.attr['r_ou'] * project.attr['sf']
    arcs = r_ou_new/project.attr['distance']
    
    extent = [-arcs,arcs,-arcs,arcs]
    text_pos = [-0.9*arcs,-0.8*arcs]
    
    #~ xx = np.linspace(-arcs,arcs,project.attr['n_bin_map']*2+1)
    ysize = int((map_in.shape[2]-1)/2)
    xsize = int((map_in.shape[3]-2)/2)
    xx = np.linspace(-arcs,arcs,map_in.shape[2])
    vch = np.linspace(-project.attr['vel_max'],project.attr['vel_max'],(project.attr['i_vel_chan']*2)+1)
    #~ print( vch)
    dt = float(len(vch)-1)/(N*(N+2))
    
    step = 1 
    pic_vch = np.arange(step*15,step=step)+(len(vch)-1)/2-step*7
    #~ for i in range(len(beam)):
    interpol = 'nearest'
    for i in range(3):
        # plot full spectrum
        if i == 0:
            fig = plt.figure(pname+'full spectrum no planet')
            
        elif i == 1:
            fig =  plt.figure(pname+'full spectrum with jupiter mass planet ')
        elif i == 2:
            fig =  plt.figure(pname+'full spectrum differential map')
        
        
        plt.xlabel('velocity [km/s]')
        plt.ylabel('flux [Jy]')
        plt.plot(vch,map_in[i,:,:,:].sum(axis=1).sum(axis=1)/ beam2px )
        
        
        # plot velocity integrated map (full)
        if i == 0:
            fig = plt.figure(pname+'velocity integrated no planet map')
        elif i == 1:
            fig =  plt.figure(pname+'velocity integrated with jupiter mass planet map')
        elif i == 2:
            fig =  plt.figure(pname+'velocity integrated differential map')
        ax = fig.add_subplot(111, aspect='auto')
        plt.imshow(map_in[i,:,:,:].sum(axis=0)*abs(vch[1]-vch[0])*1e-3,origin='lower',
            interpolation=interpol,extent=extent,aspect="auto")
        ae = AnchoredEllipse(ax.transData, width=hlp.deg2as(header['BMIN']), height=hlp.deg2as(header['BMAJ'])
                             , angle=5., loc=3, pad=0.4, frameon=False)
        ae.ellipse.set_facecolor([0.1,0.1,0.1])
        ax.add_artist(ae)
        plt.xlabel('X ["]')
        plt.ylabel('Y ["]')
        plt.colorbar().set_label('Flux [Jy/beam *km/s]')
        
        # plot velocity integrated map (cut)
        if i == 2:
            
            fig =  plt.figure(pname+'velocity integrated differential map 2')
            ax = fig.add_subplot(111, aspect='auto')
            map1 = map_in[0,:,:,:].sum(axis=0)*abs(vch[1]-vch[0])*1e-3
            map2 = map_in[1,:,:,:].sum(axis=0)*abs(vch[1]-vch[0])*1e-3
            map3 = np.abs(map1-map2)
            plt.imshow(map3,origin='lower',
                interpolation=interpol,extent=extent,aspect="auto")
                
            ae = AnchoredEllipse(ax.transData, width=hlp.deg2as(header['BMIN']), height=hlp.deg2as(header['BMAJ'])
                             , angle=5., loc=3, pad=0.4, frameon=False)
            ae.ellipse.set_facecolor([0.1,0.1,0.1])
            ax.add_artist(ae)
                
            plt.colorbar().set_label('Flux [Jy/beam *km/s]')
            plt.xlabel('X ["]')
            plt.ylabel('Y ["]')
        #~ 
        #~ # make velocity channel overview map
        if i == 0:

            fig = plt.figure(pname+' overview with no planet ')
        elif i == 1:
            fig =  plt.figure(pname+' overview with jupiter mass planet ')
        elif i == 2:
            fig =  plt.figure(pname+' differential map')
            
        grid = AxesGrid(fig, 111, # similar to subplot(132)
                        nrows_ncols = (N, N+2),
                        axes_pad = 0.0,
                        share_all=True,
                        label_mode = "1",
                        cbar_location = "right",
                        cbar_mode="single",
                        aspect="auto"
                        )
        
        vmin = 0
        #~ vmax =  np.max(hlp.conv(map_in[(len(vch)-1)/2,:,:],beam=beam[i],r_ou=r_ou_new,dist=attr['distance']))
        #~ vmax =  np.max(map_in[(len(vch)-1)/2,:,:])
        vmax =  np.max(map_in[i,:,:,:])*1.1

        #~ for t in range(N*(N+2)):
        for t in range(15):
                #k = int(round(t*dt+dt/2))
                k = pic_vch[t]

                #~ conv_map = hlp.conv(map_in[k,:,:],beam=beam[i],r_ou=r_ou_new,dist=attr['distance'])
                conv_map = map_in[i,k,:,:]

                # plot x cut
                #~ grid[t].plot(xx,conv_map[ysize,:]/vmax*arcs*0.8,'k')
                # plot y cut
                #~ grid[t].plot(xx,conv_map[:,xsize]/vmax*arcs*0.8,'w')
                
                im = grid[t].imshow(conv_map,vmin=vmin,vmax=vmax,
                    origin='lower',interpolation=interpol,extent=extent,aspect="auto")
                grid[t].text(text_pos[0],text_pos[1], '%2.2f Km/s' %(vch[k]*1e-3), fontsize=13,
                          bbox={'facecolor':'white', 'alpha':0.7, 'pad':5})
        
        cbar = grid.cbar_axes[0].colorbar(im)
        grid[0].cax.colorbar(im).set_label_text('Flux [Jy/beam]')
        grid[0].cax.toggle_label(True)
        #~ grid.axes_llc.set_xticks(rotation=70) 
        #~ print(dir(grid.axes_llc))
        
        #~ grid.cbar_axes[0].colorbar(im)
        #~ for cax in grid.cbar_axes:
            #~ cax.toggle_label(False)
        grid.axes_llc.set_ylabel('Y ["]')
        grid.axes_llc.set_xlabel('X ["]')
        #~ grid.cbar_axes[0].toggle_label(True)
        plt.setp( grid.axes_llc.xaxis.get_majorticklabels(), rotation=70 )
        fig.savefig('overview_map' +str(i)+'.pdf',bbox_inches='tight')
main()
plt.show()
