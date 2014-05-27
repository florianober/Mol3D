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
import numpy as np
import little as li
import mol3d_routines as l
import matplotlib.pyplot as plt


########################################################
### main

    
def main():
    pname = 'sp1_k07m03st00i0HCO4'
    mk_fits(pname)
    
def mk_fits(pname):

    project = l.mol3d(pname)

    #~ save fits file for ALMA CASA
    r_ou_new = project.attr['r_ou'] * project.attr['sf']
    arcs = r_ou_new/project.attr['distance']
    
    res = 2.0*li.as2deg(arcs)/project.velo_ch_map.shape[2]
    li.create_CASA_fits(project.velo_ch_map,out_name=project.attr['r_path']+pname+'_velo_ch_map.fits',resolution=res,
                        freq=project.attr['tr_freq'],deltafreq=project.attr['dtr_freq'])
    
    
main()
