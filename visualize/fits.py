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


########################################################
### main

    
def main():
    pname = 'sp1_k08m01st00i0CO3'
    models = glob.glob(l.path_results+'*_velo_ch_map.dat')
    #~ models = glob.glob('/data/fober/mol3dresults/hd_sim/*_input_file.dat')
    pname_list = []
    #~ print(models[0][-15:])
    for entry in models:
        pname = entry[32:-16]
        mk_fits(pname)

    
def mk_fits(pname):
    
    project  = l.mol3d(pname)
    print(project.pname)
    if len(project.attr)>0:
        r_ou_new = project.attr['r_ou'] * project.attr['sf']
        arcs = r_ou_new/project.attr['distance']
        
        res = 2.0*hlp.as2deg(arcs)/project.velo_ch_map.shape[2]
        hlp.create_CASA_fits(project.velo_ch_map,out_name=l.path_results+pname+'_velo_ch_map.fits',resolution=res,
                            freq=project.attr['tr_freq'],deltafreq=project.attr['dtr_freq'])
    else:
        print('skipping missing project')
    del project
main()
