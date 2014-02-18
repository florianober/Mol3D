#!/usr/bin/python
# -*- coding: utf-8 -*
""" script compare URANIA integral intensity maps with mol3d

date: 07/16/2012
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


path_mol3d  = '/data/fober/mcmolresults/'
#~ mol3d_file_name2   = 'firstlight_velo_ch_mapint2.dat'
mol3d_file_name   = 'Yaroslav_compare_velo_ch_mapint.dat'
mol3d_ch_map      = 'Yaroslav_compare_velo_ch_map.dat'

path_URANIA = '/home/fober/Doktorarbeit/Yaroslav/Result/'
URANIA_file_name  = 'gas-intmap.dat'
URANIA_file_name  = 'gas-channel.dat'

#~ URANIA_file_nameLVG  = 'gas-intmapLVG.dat'
#~ URANIA_file_nameLTE  = 'gas-intmapLTE.dat'
#~ URANIA_file_nameFEP  = 'gas-intmapFEP.dat'

########################################################
### main


def read_URANIA(file_name):
    f = open(path_URANIA+file_name,'r')
    
    # read header 
    line = np.str.split(f.readline())
    mask = line[0:2]
    vel_ch = 1
    y = 1
    x = 1
    line = np.str.split(f.readline())
    while len(line) != 0:
        if line[0] == 'mask:':
            mask = line[1]
        elif line[0] == 'line:':
            J_h  = line[1]
        elif line[0] == 'incl:':
            incl = line[1]
        elif line[0] == 'intx:':
            x = int(line[1])
        elif line[0] == 'inty:':
            y = int(line[1])   
        elif line[0] == 'vchan:':
            vel_ch = int(line[1])      
        line = np.str.split(f.readline())
    pic = np.zeros((x,y,vel_ch))
    
    line = f.readline()
    line = f.readline()
    line = f.readline()
    if vel_ch == 1:
        no = 4
    else:
        no = 6
    
    for i in range(x):
        for j in range(y):
            for ch in range(vel_ch):
                line = np.str.split(f.readline())
                if line[no] == 'NaN':
                    line[no] = 0.0
                
                pic[i,j,ch] = float(line[no])
    f.close()
    return pic
    
def read_mol3d(file_name):
    
    map_in = open(path_mol3d+file_name,'r')
    row = np.str.split(map_in.readline())
    vel_ch = int(row[0])

    row = np.str.split(map_in.readline())
    map_size = int(row[0])
    
    pic = np.zeros((map_size,map_size,vel_ch))
    row = map_in.readline()	#empty row

    for k in xrange(vel_ch):
        row = map_in.readline()	#empty row
        row = np.str.split(map_in.readline())
        ch_val = float(row[0])
        row = map_in.readline()	#empty row
        
        for i in xrange(map_size):
            for j in xrange(map_size):
                row = np.str.split(map_in.readline())
                pic[i,j,k] = float(row[2])
    
    return pic
def main():
    
    
    map_range = [-200,200,-200,200]
    
    
    #~ datei_mol3d  = open(path_mol3d+mol3d_file_name,'r')
    #~ counter = 0
    #~ pic_m = np.zeros((401,401))
    #~ datei_mol3d.readline()
    #~ datei_mol3d.readline()
    #~ for j in range(401):
        #~ for i in range(401):
            #~ counter +=1
            #~ line = np.str.split(datei_mol3d.readline())
            #~ if i < 401 and j < 401:
                #~ pic_m[i,j] = line[2]    
    #~ datei_mol3d.close()
    #~ 
    
#~ 
    #~ plt.imshow((pic_m[:,:]),interpolation='none',extent=map_range)
    #~ plt.title('calculated with mol3d (LVG)')
    #~ plt.xlabel('distance [AU]')
    #~ plt.ylabel('distance [AU]')
    #~ plt.colorbar()
    #~ plt.clim(-0.5, 450)
    
    plt.figure('Velocity Spectrum')
    
    pic = read_URANIA('gas-channelART.dat')
    spec =  np.zeros((len(pic[1,1,:]),2))
    for i in range(len(spec[:,1])):
        spec[i,1] = np.sum(pic[:,:,i])/(len(pic[1,:,1])**2)
        spec[i,0] = -1000+i*1000/40
        
    plt.plot(spec[:,0],spec[:,1],'k-',label='URANIA ART')
    
    
    
    pic = read_URANIA('gas-channelFEP.dat')
    spec =  np.zeros((len(pic[1,1,:]),2))
    for i in range(len(spec[:,1])):
        spec[i,1] = np.sum(pic[:,:,i])/(len(pic[1,:,1])**2)
        spec[i,0] = -1000+i*1000/40
        
    plt.plot(spec[:,0],spec[:,1],'b-',label='URANIA FEP')
    
    pic  =  read_mol3d('Yaroslav_compareFEP_velo_ch_map.dat')
    spec =  np.zeros((len(pic[1,1,:]),2))
    for i in range(len(spec[:,1])):
        spec[i,1] = np.sum(pic[:,:,i])/(len(pic[1,:,1])**2)
        spec[i,0] = -1000+i*1000/40
    plt.plot(spec[:,0],spec[:,1],'b+',label='mol3d FEP')
    
    pic = read_URANIA('gas-channelLTE.dat')
    spec =  np.zeros((len(pic[1,1,:]),2))
    for i in range(len(spec[:,1])):
        spec[i,1] = np.sum(pic[:,:,i])/(len(pic[1,:,1])**2)
        spec[i,0] = -1000+i*1000/40
        
    plt.plot(spec[:,0],spec[:,1],'r-',label='URANIA LTE')
    
    pic  =  read_mol3d('Yaroslav_compareLTE_velo_ch_map.dat')
    spec =  np.zeros((len(pic[1,1,:]),2))
    for i in range(len(spec[:,1])):
        spec[i,1] = np.sum(pic[:,:,i])/(len(pic[1,:,1])**2)
        spec[i,0] = -1000+i*1000/40
    plt.plot(spec[:,0],spec[:,1],'r+',label='mol3d LTE')
    
    
    pic = read_URANIA('gas-channelLVG.dat')
    spec =  np.zeros((len(pic[1,1,:]),2))
    for i in range(len(spec[:,1])):
        spec[i,1] = np.sum(pic[:,:,i])/(len(pic[1,:,1])**2)
        spec[i,0] = -1000+i*1000/40
        
    plt.plot(spec[:,0],spec[:,1],'g-',label='URANIA LVG')
    
    #~ pic = read_URANIA('gas-channelLVGold.dat')
    #~ spec =  np.zeros((len(pic[1,1,:]),2))
    #~ for i in range(len(spec[:,1])):
        #~ spec[i,1] = np.sum(pic[:,:,i])/(len(pic[1,:,1])**2)
        #~ spec[i,0] = -1000+i*1000/40
#~ 
    #~ plt.plot(spec[:,0],spec[:,1],'c-',label='URANIA LVG old')

    pic  =  read_mol3d('Yaroslav_compareLVG_velo_ch_map.dat')
    spec =  np.zeros((len(pic[1,1,:]),2))
    for i in range(len(spec[:,1])):
        spec[i,1] = np.sum(pic[:,:,i])/(len(pic[1,:,1])**2)
        spec[i,0] = -1000+i*1000/40
    plt.plot(spec[:,0],spec[:,1],'g+',label='mol3d LVG')
    
    plt.xlabel('V [m/s]')
    plt.ylabel('T_mb [K]')
    #~ plt.colorbar()
    #~ 
    #~ plt.clim(-0.5, 450)
    
    #~ plt.figure('cuts')
    #~ plt.plot(pic[200,:,0],'r+-',label='URANIA')
    #~ plt.plot(pic_m[200,:],'b+-',label='mol3d')
    
    
    
    #~ plt.figure('Mol3d lvl population')
    #~ pic  =  read_mol3d('Yaroslav_compare_velo_ch_map.dat')
    #~ spec =  np.zeros((len(pic[1,1,:]),2))
    #~ for i in range(len(spec[:,1])):
        #~ spec[i,1] = np.sum(pic[:,:,i])/(len(pic[1,:,1])**2)
        #~ spec[i,0] = -1000+i*1000/40
    #~ plt.plot(spec[:,0],spec[:,1],'r+-',label='LVG')
    #~ 
    #~ 
    #~ 

    #~ 
    #~ 
    #~ 
    #~ pic  =  read_mol3d('Yaroslav_compareLTE_velo_ch_map.dat')
    #~ spec =  np.zeros((len(pic[1,1,:]),2))
    #~ for i in range(len(spec[:,1])):
        #~ spec[i,1] = np.sum(pic[:,:,i])/(len(pic[1,:,1])**2)
        #~ spec[i,0] = -1000+i*1000/40
    #~ plt.plot(spec[:,0],spec[:,1],'g+-',label='LTE')
    #~ plt.xlabel('Velocity [m/s]')
    #~ plt.ylabel('T_mb [K]')
    
    plt.legend()
    

    
    plt.show()


main()
