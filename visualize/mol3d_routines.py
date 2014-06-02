#!/usr/bin/python
# -*- coding: utf-8 -*-
"""

objects and routines to access Mol3d data (output)


"""
from numpy import zeros, argmax, array
from astropy.convolution import convolve_fft
from astropy.convolution import Gaussian2DKernel
import os.path
import helper as hlp


# some general definitions
f = open('path_result.dat')
path_results = f.readline().split()[0]
f.close()

resulting_files = array(['_cell_boundaries.dat','_input_file.dat','_model.dat','_temp_x.dat',
                        '_velo_ch_map.dat','_velo_ch_mapint.dat','_velo_ch_mapsum.dat',
                        '_visual_xy.dat','_visual_xz.dat','_visual_yz.dat'])

# some routines to read mol3d results from file


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
        
def load_mol3d_fullvchmap(file_path,return_all=False):
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
        if return_all:       
            return pic,vch
        else:
            return pic
            
def load_mol3d_continuum_map(file_path,return_all=False):
        map_in = open(file_path)
        #~ row = map_in.readline().split()
        #~ n_lam = int(row[0])
        n_lam = 100
        row = map_in.readline().split()
        map_size = int(row[0])
        
        pic  = zeros((n_lam,map_size,map_size))
        wlen = zeros(n_lam)
        
        row = map_in.readline()	#empty row
        for k in range(n_lam):
            row = map_in.readline()	#empty row
            row = map_in.readline().split()
            wlen[k] = float(row[0])
            row = map_in.readline()	#empty row
            for j in range(map_size**2):
                row = map_in.readline().split()
                pic[k,int(row[1]),int(row[0])] = row[2]
        if return_all:
            return pic,wlen
        else:
            return pic

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
    


def get_attr(pname):
    print('This function is old and replaced by the mol3d class')
    
    input_file = open(path_results+pname+'_input_file.dat',"r")
    full = input_file.readlines()
    input_file.close()
    
    # search for key in input_file and add to dictionary
    # this has to be extended
    attr = {}
    for line in full:
        if 'r_ou' in line:
            r_ou = float(line.partition('{')[-1].rpartition('}')[0])
            attr['r_ou'] = r_ou
        if 'r_in' in line:
            r_in = float(line.partition('{')[-1].rpartition('}')[0])
            attr['r_in'] = r_in
        if 'sf =' in line:
            sf = float(line.partition('{')[-1].rpartition('}')[0])
            attr['sf'] = sf
        if 'zoom_map =' in line:
            zoom_map = float(line.partition('{')[-1].rpartition('}')[0])
            attr['zoom_map'] = zoom_map
        if 'distance' in line:
            dist = float(line.partition('{')[-1].rpartition('}')[0])
            attr['distance'] = dist
        if 'n_bin_map' in line:
            attr['n_bin_map']= int(line.partition('{')[-1].rpartition('}')[0])
        if line[0:5] == 'line ':
            attr['line'] = int(line.partition('{')[-1].rpartition('}')[0])
            #~ attr['tr_freq'] = CO_lines_freq[attr['line']-1]
            #~ attr['tr_lam'] = lines_lam[attr['line']-1]
    return attr    


# mol3d class definition

class mol3d:
    """A class to handle mol3d objects"""
    
    def __init__(self,pname):
        
        self.__pname = pname
        self.__attr = {}
        self.__velochmap = []
        self.__files = []
        for f in resulting_files:
            # search for files for this pname
            self.__files.append(os.path.isfile(path_results+pname+f))
        # search for key in input_file and add to dictionary
        # this has to be extended
        
        if self.__files[1]:
            self.__attr['r_path']       = self.get_attr_from_file('r_path')
            self.__attr['r_ou']         = float(self.get_attr_from_file('r_ou'))
            self.__attr['r_in']         = float(self.get_attr_from_file('r_in'))
            self.__attr['sf']           = float(self.get_attr_from_file('sf'))
            self.__attr['distance']     = float(self.get_attr_from_file('distance'))
            self.__attr['n_bin_map']    = int(self.get_attr_from_file('n_bin_map'))
            self.__attr['line']         = int(self.get_attr_from_file('line'))
            self.__attr['gas_cat_name'] = self.get_attr_from_file('gas_cat_name')
            self.__attr['i_vel_chan']   = int(self.get_attr_from_file('i_vel_chan'))
            self.__attr['vel_max']      = float(self.get_attr_from_file('vel_max'))
            
            # now calculate some more properties and make some consitent tests
            if self.__attr['r_path'] != path_results:
                #~ print('project has been copied')
                self.__attr['r_path'] = path_results

            self.__attr['dvelo'] = self.__attr['vel_max']/(self.__attr['i_vel_chan'])
            
            if self.__attr['gas_cat_name'] == 'co.dat':
                self.__attr['tr_freq'] = hlp.CO_lines_freq[self.__attr['line']-1]
                
            elif self.__attr['gas_cat_name'] == 'hco+@xpol.dat':
                self.__attr['tr_freq'] = hlp.HCO_lines_freq[self.__attr['line']-1]
            else:
                print("ERROR, gas type unknown")
                self.__attr['tr_freq'] = 1.0
            self.__attr['tr_lam'] = hlp.c/self.__attr['tr_freq']
            self.__attr['dtr_freq'] = self.__attr['tr_freq']/hlp.c*self.__attr['dvelo']
        else:
            print('ERROR: Could not find results')
        
        
    def load_velo_ch_map(self):
        self.__velochmap = load_mol3d_fullvchmap(path_results+self.__pname+'_velo_ch_map.dat')
        
    def return_velo_ch_map(self):
        if self.__velochmap == []:
            
            if self.__files[4]:
                self.load_velo_ch_map()
            else:
                pass
        else:
            pass
        return self.__velochmap
        
    velo_ch_map = property(return_velo_ch_map)
        
    def return_attr(self):
        return self.__attr
    
    attr = property(return_attr)
    

    def return_name(self):
        return self.__pname
    
    pname = property(return_name)
    
    def get_attr_from_file(self,key):
        input_file = open(path_results+self.__pname+'_input_file.dat',"r")
        full = input_file.readlines()
        input_file.close()
        found = False
        # search for key in input_file and add to dictionary
        # this has to be extended
        for line in full:
            val = line.split()
            if len(val) >=1:
                if val[0] == key:
                    found = True
                    key_value = line.partition('{')[-1].rpartition('}')[0]
                    
        if not(found):
            print('ERROR, could not find key')
            key_value = ''
        return key_value

