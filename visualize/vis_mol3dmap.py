#!/usr/bin/python
# -*- coding: utf-8 -*
""" script to visalize mol3d output

date: 08/15/2012
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


########################################################
### set some files paths

path_mol3d  = '/data/fober/mcmolresults/'
mol3d_file_name   = 'firstlight'

########################################################
### functions

print ''
print '################################################'
print '####                mol3d                   ####'
print '####          visualisation tool            ####'
print '################################################'



def vis_ch_map(ch_map_file):
    
    
    ch_map  = open(path_mol3d + ch_map_file, 'r')
    no_vch  = int(np.str.split(ch_map.readline())[0])
    map_len = int(np.str.split(ch_map.readline())[0])
    pic_vch = np.zeros((no_vch,map_len,map_len))
    waste   = ch_map.readline()
    vchlist = np.zeros(no_vch)
    
    for vch in range(no_vch):
        waste = ch_map.readline()
        vchlist[vch] = np.str.split(ch_map.readline())[0]
        waste   = ch_map.readline()
        for j in range(map_len+1):
            for i in range(map_len+1):
                line = np.str.split(ch_map.readline())
                if i < map_len and j < map_len:
                    Noload = np.float(line[2])
                    if  Noload == 0.0:
                        pic_vch[vch,i,j] = 0.01e-270
                    else:
                        pic_vch[vch,i,j] = line[2] 
    ch_map.close()
    
    while True:
        print 'Please choose the velocity channel:    '
        print '--------------    '
        for vch in range(no_vch):
            print '[%2.2i] : %5.0f km/s' %(vch,vchlist[vch])

        choose = raw_input('  :')
        goodinput = False
        if choose == 'x':
            break
        try:
            xx = int(choose)
            goodinput = (xx <= no_vch-1) and (xx >= 0)
        except:
            pass
        
        if goodinput:
            plt.imshow(pic_vch[xx,:,:]**0.1,interpolation='none', cmap=plt.cm.hot)
            plt.colorbar()
            plt.clim(0, 5)
            plt.show()

def vis_intmap(intmap_file):
    
    intmap = open(path_mol3d + intmap_file,'r')
    
    map_len = int(np.str.split(intmap.readline())[0])
    pic_m = np.zeros((map_len,map_len))
    waste =  intmap.readline()
    for j in range(map_len+1):
        for i in range(map_len+1):
            line = np.str.split(intmap.readline())
            if i < map_len and j < map_len:
                Noload = np.float(line[2])
                if  Noload == 0.0:
                    pic_m[i,j] = 0.01e-99
                else:
                    pic_m[i,j] = line[2]    
    intmap.close()
    a = 0.1
    while True:
        plt.imshow((pic_m[:,:])**a,interpolation='none', cmap=plt.cm.hot)
        #plt.title(title,fontsize=14)
        plt.show()
        colorad = raw_input('please enter exponent (default = 0.1) (x quits):   ')
        if colorad == 'x':
            break
        else:
            try:
                a = float(colorad)
            except:
                print 'wrong input. please try again'
            print 'scaled to array^%2.2f' %a

    #plt.imshow(pic_m**0.1,interpolation='none', cmap=plt.cm.hot)
    #plt.imshow(pic_m,interpolation='none', cmap=plt.cm.hot)
    #plt.colorbar()
    #plt.clim(-0.5, 0.5)
    #plt.show()
    
    
########################################################
### main program

def main():
    
    print ''
    try:
        p_name = sys.argv[1]
        print 'project given : '+ sys.argv[1]
        got_input = True

    except: 
        print 'no project given, using : '+ mol3d_file_name
        p_name = mol3d_file_name
        got_input = False
    
    print ''

    #search for results
    
    fileexist = []
    fileexist.append([False ,'_velo_ch_map.dat'])
    fileexist.append([False ,'_velo_ch_mapint.dat'])
    fileexist.append([False ,'_velo_ch_mapsum.dat'])
        
    fileexist.append([False ,'densH2_xz.dat'])

    for j in xrange(len(fileexist)):
        fileexist[j][0] = os.path.exists(path_mol3d+p_name+fileexist[j][1])
    
    #file_open = open(path_mol3d+sys.argv[1],'r')
    
    found = False
    for j in xrange(len(fileexist)):            
        if fileexist[j][0] == True:
            found = True
            
        
    if not found:    
        print 'no files for given project found'
        sys.exit()
        
    #show data available for given project
    while True:
        print ''
        print 'Please choose:    '
        print '--------------    '
        if fileexist[0][0]:
            print '[01] show velocity channel maps'
        if fileexist[1][0]:
            print '[02] show integrated intensity map'
        if fileexist[2][0]:
            print '[03] show velocity spectrum'
            
        print '--------------    '
        choose_plot = raw_input(' :')
        if choose_plot == 'x':  # if input x than exit
            print 'bye bye'
            break
        
        try:
            int(choose_plot)
            is_integer = True
        except:
            is_integer = False
            
        if is_integer and int(choose_plot) < 4:
            if int(choose_plot) == 1:
                vis_ch_map(p_name+fileexist[0][1])
                
            elif int(choose_plot) == 2:
                vis_intmap(p_name+fileexist[1][1])
                
            elif int(choose_plot) == 3:
                print 'TbD: not yet implemented'
        else:
            print 'wrong input, try again'
            print ''
main()
