#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
make fits file from velocity channel map(s)

date: 26/05/2014
author: Florian Ober
email: fober@astrophysik.uni-kiel.de

"""
######################################################## 
### import some packages
#
import os
import sys
import time
import math
import glob
import tarfile
import numpy as np
import helper as hlp
import mol3d_routines as l
import matplotlib.pyplot as plt
from astropy.io import fits as pf

path_fits = '/mnt/kronos/data2/fober/fits/'
path_save = '/mnt/kronos/data3/fober/diff/'
########################################################
### main

def main():
    
    models = glob.glob(path_fits+'sp1*_velo_ch_map_alma??.fits')
    pname_list = []
    for entry in models:
        pname = entry[len(path_fits)+4:-24]
        alma_conf = entry[-7:-5]
        #~ print(entry[len(path_fits):-24])
        if os.path.isfile(path_fits+'s1_'+pname+'_velo_ch_map_alma'+alma_conf+'.fits'):
            pname_list.append(pname)
            print(pname)
            header = pf.open(path_fits+'s1_'+pname+'_velo_ch_map_alma'+alma_conf+'.fits')[0].header
            s1  = pf.open(path_fits+'s1_'+pname+'_velo_ch_map_alma'+alma_conf+'.fits')[0].data
            sp1 = pf.open(path_fits+'sp1_'+pname+'_velo_ch_map_alma'+alma_conf+'.fits')[0].data
            
            s1_noisy  = pf.open(path_fits+'s1_'+pname+'_velo_ch_map_alma'+alma_conf+'.noisy.fits')[0].data
            sp1_noisy = pf.open(path_fits+'sp1_'+pname+'_velo_ch_map_alma'+alma_conf+'.noisy.fits')[0].data
            
            # write differential image into a new fits file
            # first clean image
            diff_map = np.zeros((s1.shape),dtype=np.float32)
            for i in range(s1.shape[1]):
                diff_map[0,i,:,:] = np.abs(s1[0,i,:,:]-sp1[0,i,:,:])
            hlp.write_image2fits(diff_map,file_name=path_save+pname+'.clean.diff_image_alma%2.2d.fits' %(int(alma_conf)),header=header)
            spec = np.zeros(s1.shape[1])
            #for w in range(len(s1.shape[1])):
            #    spec[w] = np.sum(s1[0,w,:,:])
            #np.savetxt(path_save+pname+'.clean.s1_sed_alma%2.2d.dat' %(int(alma_conf)),spec)


            # write differential image into a new fits file
            # now noisy image
            diff_map = np.zeros((s1_noisy.shape),dtype=np.float32)
            for i in range(s1_noisy.shape[1]):
                diff_map[0,i,:,:] = np.abs(s1_noisy[0,i,:,:]-sp1_noisy[0,i,:,:])
            hlp.write_image2fits(diff_map,file_name=path_save+pname+'.noisy.diff_image_alma%2.2d.fits' %(int(alma_conf)),header=header)
            #~ for w in range(50,91,2):
                #~ plt.figure('clean %2.2g' %w)
                #~ plt.imshow(np.abs(s1[0,w,:,:]-sp1[0,w,:,:]),origin='lower')
                #~ plt.colorbar()
                #~ 
                #~ plt.figure('noisy%2.2g' %w)
                #~ plt.imshow(np.abs(s1_noisy[0,w,:,:]-sp1_noisy[0,w,:,:]),origin='lower')
                #~ plt.colorbar()
    
def tocasa():
    #~ models = glob.glob('/data/fober/mol3dresults/hd_sim/*_model.dat')
    models = glob.glob('/data/fober/mol3dresults/hd_sim/*_velo_ch_map.dat')
    pname_list = []
    #~ print(models[0][-15:])
    for entry in models:
        pname = entry[32:-16]
        #~ pname = entry[32:-10]
        #~ print(pname)
        if not(os.path.isfile(path_fits+pname+'_velo_ch_map.fits')):
        #~ if not(os.path.isfile(l.path_results+pname+'_velo_ch_map.dat')):
            if not('m00' in pname):
                #~ print(pname)
                mk_fits(pname)

    
def mk_fits(pname):
    project  = l.mol3d(pname)
    print(project.pname)
    if len(project.attr)>0:
        r_ou_new = project.attr['r_ou'] * project.attr['sf']
        arcs = r_ou_new/project.attr['distance']
        
        res = 2.0*hlp.as2deg(arcs)/project.velo_ch_map.shape[2]
        hlp.create_CASA_fits(project.velo_ch_map,out_name=path_fits+pname+'_velo_ch_map.fits',
                             object_name='IRAS 04302+2247',resolution=res,
                             freq=project.attr['tr_freq'],deltafreq=project.attr['dtr_freq'])
    else:
        print('skipping missing project')
    del project
    
tocasa()
#~ main()
