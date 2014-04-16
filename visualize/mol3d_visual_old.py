#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
easy visualisation tool for mol3d.
ver.: 0.1

NOTE: This file can be used and has some advantages, but it is obsolent, please use
      the new file: visual_mol3d.py

How to get started: 

date   : 30/05/2013
author : F. Ober 
email  : fober@astrophysik.uni-kiel.de
"""


import numpy as np
import matplotlib.pyplot as plt
import os
import sys

pyfits_log = False
import little as l

try:
    import pyfits as pf
    pyfits_log = True
except:
    print('')
    print('could not find pyfits. Thus it is not possible to create/manipulate fits files.')
    


path_results = '/data/fober/mol3dresults/'
home         = '../'

    

try:
    p_name = sys.argv[1]
    inp = False 
except:
    p_name = 'example'
    inp = True

##################### edit this area to create your desired plots ##################################

create_plots = dict(Temp_midplane=False,                      #plot temperature in the midplane
                    Spectrum=False,                           #plot the velocity spectrum
                    Velo_ch_map=False,                         #Intensity map for every velocity channel (selectable)
                    Intensity_map=False,                      #Intensity map integrated over all velocity channels
                    Velo_map_xz=False,                        #Plot velocity in the xz (y = 0) plane
                    Temp_map_xz=False,                        #Plot temperature in the xz (y = 0) plane
                    Velo_map_xy=False,                        #Plot velocity in the xy (z = 0) plane
                    Temp_map_xy=False,                        #Plot temperature in the xy (z = 0) plane
                    DensH2_map_xz=True,                      #Plot H2 distributione in the xz (y = 0) plane
                    DensH2_map_xy=True                       #Plot H2 distributione in the xy (z = 0) plane

                    )
####################              end of editable area            ##################################				

file_names = dict(Temp_midplane='_temp_x.dat',Velo_ch_map='_velo_ch_map.dat',
                  Intensity_map='_velo_ch_mapint.dat',Spectrum='_velo_ch_mapsum.dat',
                  Velo_map_xz='_velo_xz.dat',Temp_map_xz='_temp_xz.dat',
                  Velo_map_xy='_velo_xy.dat',Temp_map_xy='_temp_xy.dat',
                  DensH2_map_xy='_densH2_xy.dat',DensH2_map_xz='_densH2_xz.dat',
                  )

def main():
    
    #search for results with given p_name
    results_types = []
    for n,l in create_plots.items():
        if l:
            t_file = path_results+p_name+file_names[n]
            if os.path.exists(t_file):
                results_types.append(n)
            else:
                print('result not found: ' +n)

    
    
    if len(results_types) == 0:
        sys.exit('no results found')

    for i in range(len(results_types)):
        t_file = path_results+p_name+file_names[results_types[i]]

        if results_types[i] == 'Temp_midplane':
            oneD(t_file,i=0)
        elif results_types[i] == 'Spectrum':
            oneD(t_file,i=1)
        elif results_types[i] == 'Intensity_map':
            present_map(t_file,i=0)
        elif results_types[i] == 'Velo_ch_map':
            present_map(t_file,i=1)
        elif results_types[i] == 'Velo_map_xz':
            present_map(t_file,i=2)
        elif results_types[i] == 'Temp_map_xz':
            present_map(t_file,i=3)
        elif results_types[i] == 'Velo_map_xy':
            present_map(t_file,i=2)
        elif results_types[i] == 'Temp_map_xy':
            present_map(t_file,i=3)
        elif results_types[i] == 'DensH2_map_xy':
            present_map(t_file,i=4)
        elif results_types[i] == 'DensH2_map_xz':
            present_map(t_file,i=4)

def present_map(file_path,i=0):

    if i in [0,2,3,4]:
        map_in = open(file_path)
        row = np.str.split(map_in.readline())
        map_size = int(row[0])
        row = map_in.readline()	#empty row
        pic = np.zeros((map_size,map_size))
        for j in range(map_size):
            for k in range(map_size):
                row = np.str.split(map_in.readline())
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
                    
        m_range = [k_min,k_max,j_min,j_max]
        #m_range = [0,0,0,0]
        xlab = 'distance [AU]'
        ylab = xlab
        expo = 1.0
        maxi = np.max(pic**expo)
        cont = ''

        if i == 0:
            fig_name = 'Intensity map'
            m_range = [-120,120,-120,120]
            print((np.sum(pic)))
        elif i == 2:
            fig_name = 'Velocity map ('+file_path[-6:-4]+'-plane)'
        elif i == 3:
            fig_name = 'Temperature map ('+file_path[-6:-4]+'-plane)'
            maxi = 1
            #~ cont = [12,20,30,40]
            cont = [30,35,40,45]
            
        elif i == 4:
            fig_name = 'H2 density map ('+file_path[-6:-4]+'-plane)'
            expo = 0.2
            maxi = np.max(pic**expo)
        l.plot_image((pic**expo)/maxi,num=file_path,ylab=ylab,xlab=xlab,map_range=m_range,color_loc=True,fig_title=fig_name,cont=cont)
        #~ plt.clim(0, 70)
    elif i == 1:
        map_in = open(file_path)
        row = np.str.split(map_in.readline())
        vel_ch = int(row[0])+1

        row = np.str.split(map_in.readline())
        map_size = int(row[0])
        pic = np.zeros((map_size,map_size))
        row = map_in.readline()	#empty row
        map_ch = 1
        for k in range(vel_ch-1):
            row = map_in.readline()	#empty row
            row = np.str.split(map_in.readline())
            ch_val = float(row[0])
            row = map_in.readline()	#empty row
            for j in range(map_size**2):
                row = np.str.split(map_in.readline())
                pic[int(row[1]),int(row[0])] = row[2]
            if k == map_ch:
                m_range =[-120,120,-120,120]
                print(( np.sum(pic)))
                l.plot_image(pic**0.1,num=file_path,color_loc=True,map_range=m_range,fig_title='%3.0f m/s' %ch_val)
                #pic2 = rebin(pic,(100,100))
                #l.plot_image(pic2**0.4,num='map2 '+str(i),color_loc=True,map_range=m_range,fig_title='%3.0f m/s' %ch_val)
                #~ print np.sum(pic)/len(pic)**2
                
                #~ plt.figure('pic cut ' +str(len(pic)))
                #~ x = np.linspace(-200,200,len(pic[:,0]))
                #~ plt.plot(x,pic[(len(pic[:,0])-1)/2,:])
                #plt.figure('pic cut2')
                #x = np.linspace(-200,200,len(pic2[:,0]))
                #plt.plot(x,pic2[(len(pic2[:,0])/2-1),:])
                break

def rebin(a,shape):
    sh = shape[0],a.shape[0]//shape[0],shape[1],a.shape[1]//shape[1]
    return a.reshape(sh).mean(-1).mean(1)

           
def oneD(file_path,i=0):

    pic = np.loadtxt(file_path)
    #print pic
    plt.figure(file_path)

    if i == 0:
            plt.xlabel('distance r [AU]')
            plt.ylabel('temperature [K]')
            title='midplane temperature distribution'
            maxi = 1
    elif i == 1:
            plt.ylim(0,1.1)
            plt.xlabel('velocity [m/s]')
            plt.ylabel('normalized intensity')
            title='velocity spectrum'
            maxi = np.max(pic[:,1])

    else:
        title=''
        plt.title(title,fontsize=14)
        
    plt.plot(pic[:,0],pic[:,1]/maxi)

main()
plt.show()
