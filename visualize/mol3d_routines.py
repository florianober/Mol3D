#!/usr/bin/python
# -*- coding: utf-8 -*-
"""

some routines to access Mol3d data output


"""
from numpy import zeros, argmax

def mol3d_intmap(map_in,vch):
    chmap = zeros((map_in.shape[1],map_in.shape[2]))
    
    for j in range(map_in.shape[1]):
        for k in range(map_in.shape[2]):
            chmap[j,k] = vch[argmax(map_in[:,j,k])]
    return chmap


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
        for k in xrange(vel_ch-1):
            row = map_in.readline()	#empty row
            row = map_in.readline().split()
            ch_val = float(row[0])
            row = map_in.readline()	#empty row
            for j in xrange(map_size**2):
                row = map_in.readline().split()
                pic[int(row[1]),int(row[0])] = row[2]
            if k == map_ch:
                print 'loaded %6.3f [km/s] channel map' %(ch_val*0.001)
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
        for k in xrange(vel_ch):
            row = map_in.readline()	#empty row
            row = map_in.readline().split()
            vch[k] = float(row[0])
            row = map_in.readline()	#empty row
            for j in xrange(map_size**2):
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
    for j in xrange(map_size):
        for k in xrange(map_size):
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
