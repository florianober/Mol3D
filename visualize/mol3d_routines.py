#!/usr/bin/python
# -*- coding: utf-8 -*-
"""

objects and routines to access Mol3d data (output)


"""
from numpy import zeros, array
from numpy import linspace
#~ from astropy.convolution import convolve_fft
#~ from astropy.convolution import Gaussian2DKernel
from astropy.io import fits as pf
import os.path
import helper as hlp


# some general definitions
FILE_IN = open('path_result.dat')
PATH_RESULTS = FILE_IN.readline().split()[0]
FILE_IN.close()

MOL3D_FILES = array(['_cell_boundaries.dat', '_input_file.dat',
                     '_model.dat', '_temp_x.dat',
                     '_velo_ch_map.dat', '_velo_ch_mapint.dat',
                     '_velo_ch_mapsum.dat',
                     '_visual_xy.dat', '_visual_xz.dat',
                     '_visual_yz.dat', '_velo_ch_map.fits.gz'])

# some routines to read mol3d results from file


def load_mol3d_zerovchmap(file_path, ch=-1):
    map_in = open(file_path)
    row = map_in.readline().split()
    vel_ch = int(row[0])+1

    row = map_in.readline().split()
    map_size = int(row[0])
    pic = zeros((map_size, map_size))
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
            pic[int(row[1]), int(row[0])] = row[2]
        if k == map_ch:
            print(('loaded %6.3f [km/s] channel map' %(ch_val*0.001)))
            break
    return pic

def load_mol3d_fullvchmap(file_path, return_all=False):
    map_in = open(file_path)
    row = map_in.readline().split()
    vel_ch = int(row[0])
    row = map_in.readline().split()
    map_size = int(row[0])

    pic = zeros((vel_ch, map_size, map_size))
    vch = zeros(vel_ch)

    row = map_in.readline()	#empty row
    for k in range(vel_ch):
        row = map_in.readline()	#empty row
        row = map_in.readline().split()
        vch[k] = float(row[0])
        row = map_in.readline()	#empty row
        for j in range(map_size**2):
            row = map_in.readline().split()
            pic[k, int(row[1]), int(row[0])] = row[2]

    if return_all:
        return pic, vch
    else:
        return pic

def load_mol3d_vchmap(file_path, return_all=False):
    fits = pf.open(file_path)

    data = fits[0].data
    header = fits[0].header

    fits.close()

    if return_all:
        return data, header
    else:
        return data

def load_mol3d_continuum_map(file_path, return_all=False):
    map_in = open(file_path)
    row = map_in.readline().split()
    n_lam = int(row[0])
    n_lam = 100
    row = map_in.readline().split()
    map_size = int(row[0])

    pic = zeros((n_lam, map_size, map_size))
    wlen = zeros(n_lam)

    row = map_in.readline()	#empty row
    for k in range(n_lam):
        row = map_in.readline()	#empty row
        row = map_in.readline().split()
        wlen[k] = float(row[0])
        row = map_in.readline()	#empty row
        for j in range(map_size**2):
            row = map_in.readline().split()
            pic[k, int(row[1]), int(row[0])] = row[2]
    if return_all:
        return pic, wlen
    else:
        return pic

def load_mol3d_map(file_path):

    map_in = open(file_path)
    row = map_in.readline().split()
    map_size = int(row[0])
    row = map_in.readline()     #empty row
    pic = zeros((map_size, map_size))
    for j in range(map_size):
        for k in range(map_size):
            row = map_in.readline().split()
            #~ if j == 0 and k == 0:
                #~ j_min = float(row[0])
                #~ k_min = float(row[1])
            #~ elif j == map_size-1 and k == map_size-1:
                #~ j_max = float(row[0])
                #~ k_max = float(row[1])
            pic[k, j] = row[2]
    return pic

# mol3d class definition

class mol3d:
    """A class to handle mol3d objects"""

    def __init__(self, pname, path_results=PATH_RESULTS):
        self.__pname = pname
        self.__attr = {}
        self.__velochmap = []
        self.__vch = []
        self.__files = []
        self.__path_results = path_results
        for item in MOL3D_FILES:
            # search for files for this pname
            self.__files.append(os.path.isfile(self.__path_results +
                                               pname + item))
        # search for key in input_file and add to dictionary
        # this has to be extended

        if self.__files[1]:
            self.__attr['r_path'] = self.get_attr_from_file('r_path')
            self.__attr['r_ou'] = float(self.get_attr_from_file('r_ou'))
            self.__attr['r_in'] = float(self.get_attr_from_file('r_in'))
            self.__attr['sf'] = float(self.get_attr_from_file('sf'))
            self.__attr['distance'] = float(self.get_attr_from_file('distance'))
            self.__attr['n_bin_map'] = int(self.get_attr_from_file('n_bin_map'))
            self.__attr['line'] = int(self.get_attr_from_file('line'))
            self.__attr['gas_cat_name'] = self.get_attr_from_file('gas_cat_name')
            self.__attr['i_vel_chan'] = int(self.get_attr_from_file('i_vel_chan'))
            self.__attr['vel_max'] = float(self.get_attr_from_file('vel_max'))
            self.__attr['zoom_map'] = float(self.get_attr_from_file('zoom_map'))

            # now calculate some more properties and make some consistency tests
            #~ if self.__attr['r_path'] != path_results:
                #~ self.__attr['r_path'] = path_results

            self.__attr['dvelo'] = (self.__attr['vel_max'] /
                                    (self.__attr['i_vel_chan']))

            if self.__attr['gas_cat_name'] == 'co.dat':
                self.__attr['tr_freq'] = (hlp.CO_lines_freq[self.__attr['line']-1])

            elif self.__attr['gas_cat_name'] == 'hco+@xpol.dat':
                self.__attr['tr_freq'] = (hlp.HCO_lines_freq[self.__attr['line']-1])
            else:
                print("ERROR, gas type unknown")
                self.__attr['tr_freq'] = 1.0
            self.__attr['tr_lam'] = hlp.c/self.__attr['tr_freq']
            self.__attr['dtr_freq'] = (self.__attr['tr_freq'] /
                                       hlp.c*self.__attr['dvelo'])
        else:
            print('ERROR: Could not find results')

    def return_velo_ch_map(self):
        if self.__velochmap == []:
            # load map into memory, first try open the fits file (much faster)
            if self.__files[10]:
                self.__velochmap = load_mol3d_vchmap(self.__path_results +
                                                     self.__pname +
                                                     MOL3D_FILES[10])
            # if the fits file is not available try open the ascii file
            # (older versions/results)
            elif self.__files[4]:
                self.__velochmap = load_mol3d_fullvchmap(self.__path_results +
                                                         self.__pname +
                                                         MOL3D_FILES[4])
            else:
                print('ERROR: could not find velocity channel maps')
        else:
            pass
        return self.__velochmap

    velo_ch_map = property(return_velo_ch_map)

    def return_vch_array(self):
        if self.__vch == []:
            self.__vch = linspace(-self.attr['vel_max'],
                                  self.attr['vel_max'],
                                  (self.attr['i_vel_chan']*2)+1)
        else:
            pass
        return self.__vch

    vch = property(return_vch_array)

    def return_attr(self):
        return self.__attr

    attr = property(return_attr)


    def return_name(self):
        """ return the project name """
        return self.__pname

    pname = property(return_name)

    def return_path_results(self):
        """ return the project path """
        return self.__path_results

    path_results = property(return_path_results)

    def get_attr_from_file(self, key):
        input_file = open(self.__path_results+self.__pname +
                          '_input_file.dat', "r")
        full = input_file.readlines()
        input_file.close()
        found = False
        # search for key in input_file and add to dictionary
        # this has to be extended
        for line in full:
            val = line.split()
            if len(val) >= 1:
                if val[0] == key:
                    found = True
                    key_value = line.partition('{')[-1].rpartition('}')[0]

        if not found:
            print('ERROR, could not find key')
            key_value = ''
        return key_value
