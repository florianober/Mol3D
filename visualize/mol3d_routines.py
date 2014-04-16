#!/usr/bin/python
# -*- coding: utf-8 -*-
"""

some routines to access Mol3d data output


"""
from numpy import zeros, argmax
from astropy.convolution import convolve_fft
from astropy.convolution import Gaussian2DKernel

def mol3d_intmap(map_in,vch):
    chmap = zeros((map_in.shape[1],map_in.shape[2]))
    
    for j in range(map_in.shape[1]):
        for k in range(map_in.shape[2]):
            chmap[j,k] = vch[argmax(map_in[:,j,k])]
    return chmap
    
f = open('path_result.dat')
path_results = f.readline().split()[0]
f.close()

def load_mol3d_zerovchmap(file_path,ch=-1):
        map_in = open(file_path)
        row = map_in.readline().split()
        vel_ch = int(row[0])+1

        row = map_in.readline().split()
        map_size = int(row[0])
        pic = zeros((map_size,map_size))
        row = map_in.readline()	#empty row
        if ch != -1:
            map_ch = ch
        else:
            map_ch = (vel_ch-2)/2
        for k in range(vel_ch-1):
            row = map_in.readline()	#empty row
            row = map_in.readline().split()
            ch_val = float(row[0])
            row = map_in.readline()	#empty row
            for j in range(map_size**2):
                row = map_in.readline().split()
                pic[int(row[1]),int(row[0])] = row[2]
            if k == map_ch:
                print(('loaded %6.3f [km/s] channel map' %(ch_val*0.001)))
                break
        return pic
        
def load_mol3d_fullvchmap(file_path):
        map_in = open(file_path)
        row = map_in.readline().split()
        vel_ch = int(row[0])
        row = map_in.readline().split()
        map_size = int(row[0])
        
        pic = zeros((vel_ch,map_size,map_size))
        vch = zeros(vel_ch)
        
        row = map_in.readline()	#empty row
        for k in range(vel_ch):
            row = map_in.readline()	#empty row
            row = map_in.readline().split()
            vch[k] = float(row[0])
            row = map_in.readline()	#empty row
            for j in range(map_size**2):
                row = map_in.readline().split()
                pic[k,int(row[1]),int(row[0])] = row[2]
        return pic,vch


def load_mol3d_map(file_path):

    map_in = open(file_path)
    row = map_in.readline().split()
    map_size = int(row[0])
    row = map_in.readline()	#empty row
    pic = zeros((map_size,map_size))
    #~ print map_size
    for j in range(map_size):
        for k in range(map_size):
            row = map_in.readline().split()
            if j == 0 and k == 0:
                j_min = float(row[0])
                k_min = float(row[1])
            elif j == map_size-1 and k == map_size-1:
                j_max = float(row[0])
                k_max = float(row[1])
            #~ if row[2] == 'NaN':
                #~ pic[k,j] = 0
            #~ else:
            pic[k,j] = row[2]
    return pic

def conv(image, beam=0.5,r_ou=200,dist=140,gg=''):
    from numpy import sqrt 
    from numpy import log
    import sys
    
    N = image.shape[0]
    #convert FWHM -> stddev
    beam_x = beam
    beam_y = beam
    beam_stddev_x = beam_x/(2.*sqrt(2*log(2)))
    beam_stddev_y = beam_y/(2.*sqrt(2*log(2)))
    
    
    arcs = r_ou / dist
    pxwitdth = 2.0*arcs/N 
    g_width = beam_stddev_x/pxwitdth

    gauss = Gaussian2DKernel(width=g_width)  #stddev width in px
    
    # if gg is a given Kernel
    if gg != '':
        z1 = convolve_fft(image, gg,normalize_kernel=True)
    else:
        #~ print('here')
        z1 = convolve_fft(image, gauss,normalize_kernel=True)
    return z1

def get_attr(pname):
    
    input_file = open(path_results+pname+'_input_file.dat',"r")
    full = input_file.readlines()
    input_file.close()
    
    # search for key in input_file and add to dictionary
    attr = {}
    for line in full:
        if 'r_ou' in line:
            r_ou = float(line.partition('{')[-1].rpartition('}')[0])
            attr['r_ou'] = r_ou
        if 'distance' in line:
            dist = float(line.partition('{')[-1].rpartition('}')[0])
            attr['distance'] = dist
    return attr
