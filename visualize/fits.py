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
    mk_fits(pname)
    
def mk_fits(pname):
    
    project = l.mol3d(pname)
    #~ files = glob.glob(project.attr['r_path']+'*'+'model.dat')
    #~ for f in files:
        #~ print(f)
    print(project.attr)
    #~ sys.exit()
    #~ save fits file for ALMA CASA
    r_ou_new = project.attr['r_ou'] * project.attr['sf']
    arcs = r_ou_new/project.attr['distance']
    
    res = 2.0*hlp.as2deg(arcs)/project.velo_ch_map.shape[2]
    hlp.create_CASA_fits(project.velo_ch_map,out_name=project.attr['r_path']+pname+'_velo_ch_map.fits',resolution=res,
                        freq=project.attr['tr_freq'],deltafreq=project.attr['dtr_freq'])
    
    
main()
