#!/usr/bin/python
# -*- coding: utf-8 -*-
#------------------------------------------------------------------------------#
# This file is part of Mol3D.
#
#    Mol3D is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    Mol3D is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with Mol3D.  If not, see <http://www.gnu.org/licenses/>.
#
#    Diese Datei ist Teil von Mol3D.
#
#    Mol3D ist Freie Software: Sie können es unter den Bedingungen
#    der GNU General Public License, wie von der Free Software Foundation,
#    Version 3 der Lizenz oder (nach Ihrer Wahl) jeder späteren
#    veröffentlichten Version, weiterverbreiten und/oder modifizieren.
#
#    Mol3D wird in der Hoffnung, dass es nützlich sein wird, aber
#    OHNE JEDE GEWÄHRLEISTUNG, bereitgestellt; sogar ohne die implizite
#    Gewährleistung der MARKTFÄHIGKEIT oder EIGNUNG FÜR EINEN BESTIMMTEN ZWECK.
#    Siehe die GNU General Public License für weitere Details.
#
#    Sie sollten eine Kopie der GNU General Public License zusammen mit diesem
#    Programm erhalten haben. Wenn nicht, siehe <http://www.gnu.org/licenses/>.
#------------------------------------------------------------------------------#
"""

objects and routines to access Mol3d data (output)


"""
from numpy import zeros, array, linspace, pi, genfromtxt
from astropy.io import fits as pf
import os.path
import helper as hlp

# some general definitions
try:
    FILE_IN = open('path_result.dat')
    PATH_RESULTS = FILE_IN.readline().split()[0]
    FILE_IN.close()
except:
    PATH_RESULTS = './'

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

def load_from_fits(file_path, return_all=False):
    """" load a map from a fits file"""
    fits = pf.open(file_path)

    data = fits[0].data
    header = fits[0].header

    fits.close()
    if return_all:
        return data, header
    else:
        return data

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
        self.__continuum_maps_mono = []
        self.__continuum_maps_bin = []
        self.__continuum_maps_raytrace = []
        self.__vch = []
        self.__continuum_wlen= []
        self.__path_results = path_results

        # search for key in input_file and add to dictionary
        # this has to be extended
        
        if os.path.isfile(self.__path_results + pname + '_input_file.dat'):
            self.__attr['r_path'] = self.get_attr_from_file('r_path')
            self.__attr['r_ou'] = float(self.get_attr_from_file('r_ou'))
            self.__attr['r_in'] = float(self.get_attr_from_file('r_in'))
            self.__attr['distance'] = float(self.get_attr_from_file('distance'))
            
            self.__attr['th_map'] = float(self.get_attr_from_file('th_map'))
            self.__attr['ph_map'] = float(self.get_attr_from_file('ph_map'))
            self.__attr['n_bin_map'] = int(self.get_attr_from_file('n_bin_map'))

            self.__attr['gas_cat_name'] = self.get_attr_from_file('gas_cat_name')
            self.__attr['line'] = int(self.get_attr_from_file('line'))
            
            self.__attr['i_vel_chan'] = int(self.get_attr_from_file('i_vel_chan'))
            self.__attr['vel_max'] = float(self.get_attr_from_file('vel_max'))
            self.__attr['zoom_map'] = float(self.get_attr_from_file('zoom_map'))

            # now calculate some more properties and make some consistency tests
            #~ if self.__attr['r_path'] != path_results:
                #~ self.__attr['r_path'] = path_results
            self.__attr['arcs'] = self.__attr['r_ou']/self.__attr['distance']
            self.__attr['image_arcs'] = self.__attr['arcs']/self.__attr['zoom_map']
            self.__attr['image_r_ou'] = self.__attr['r_ou']/self.__attr['zoom_map']
            
            self.__attr['dvelo'] = (self.__attr['vel_max'] /
                                    (self.__attr['i_vel_chan']))

            if self.__attr['gas_cat_name'] == 'co.dat':
                self.__attr['tr_freq'] = (hlp.CO_lines_freq[self.__attr['line']-1])

            elif self.__attr['gas_cat_name'] == 'hco+@xpol.dat':
                self.__attr['tr_freq'] = (hlp.HCO_lines_freq[self.__attr['line']-1])

            elif self.__attr['gas_cat_name'] == 'hcn.dat':
                self.__attr['tr_freq'] = (hlp.HCN_lines_freq[self.__attr['line']-1])
                
            elif self.__attr['gas_cat_name'] == 'hnc.dat':
                self.__attr['tr_freq'] = (hlp.HNC_lines_freq[self.__attr['line']-1])

            elif self.__attr['gas_cat_name'] == 'c18o.dat':
                self.__attr['tr_freq'] = (hlp.C18O_lines_freq[self.__attr['line']-1])

            elif self.__attr['gas_cat_name'] == 'cs@lique.dat':
                self.__attr['tr_freq'] = (hlp.CS_lines_freq[self.__attr['line']-1])

            elif self.__attr['gas_cat_name'] == '13co.dat':
                self.__attr['tr_freq'] = (hlp._13CO_lines_freq[self.__attr['line']-1])
            else:
                print("ERROR, gas type unknown")
                self.__attr['tr_freq'] = 1.0
            self.__attr['tr_lam'] = hlp.c/self.__attr['tr_freq']
            self.__attr['dtr_freq'] = (self.__attr['tr_freq'] /
                                       hlp.c*self.__attr['dvelo'])

            self.set_units()
        else:
            print('ERROR: Could not find results (input file)')

    def __getitem__(self, name):
        try:
            return self.__attr[name]
        except:
            print("ERROR, '%s' not found, available keys are:" %name)
            for key in self.__attr.keys():
                print(key)

    def set_units(self):
        self.__unit = {}
        self.__unit['PX'] = 1.0
        self.__unit['AS'] = ((self.__attr['n_bin_map'] + 0.5)/self.__attr['image_arcs'])**2
        self.__unit['TB'] = self.__unit['AS'] * (3600*180/pi)**2 * 1e-26 * \
                            (hlp.c/self.__attr['tr_freq'])**2/(2.0*hlp.k)

    def __get_unit_value(self, unit):

        if unit in self.__unit:
            unit_value = self.__unit[unit]
        else:
            print('Desired unit "%s" not found' %(unit))
            print('Fallback to Jy/pixel')
            unit_value = self.__unit['PX']
        return unit_value

    def return_velo_ch_map(self, unit='PX'):
        if self.__velochmap == []:
            # load map into memory, try open the fits file
            file_name = self.__path_results + self.__pname + '_velo_ch_map.fits.gz'
            if os.path.isfile(file_name):
                self.__velochmap = load_from_fits(file_name)
            else:
                return []
        else:
            pass
        return self.__velochmap * self.__get_unit_value(unit)

    velo_ch_map = property(return_velo_ch_map)

    def return_continuum_maps_mono(self, unit='PX'):
        if self.__continuum_maps_mono == []:
            # load maps into memory
            file_name = self.__path_results + self.__pname +                   \
                        '_continuum_map_mono.fits.gz'
            if os.path.isfile(file_name):
                self.__continuum_maps_mono = load_from_fits(file_name)
            else:
                return []
        else:
            pass
        return self.__continuum_maps_mono * self.__get_unit_value(unit)

    continuum_maps_mono = property(return_continuum_maps_mono)

    def return_continuum_maps_bin(self, unit='PX'):
        if self.__continuum_maps_bin == []:
            # load maps into memory
            file_name = self.__path_results + self.__pname +                   \
                        '_continuum_map_bin.fits.gz'
            if os.path.isfile(file_name):
                self.__continuum_maps_bin = load_from_fits(file_name)
            else:
                return []
        else:
            pass
        return self.__continuum_maps_bin * self.__get_unit_value(unit)

    continuum_maps_bin = property(return_continuum_maps_bin)

    def return_continuum_maps_raytrace(self, unit='PX'):
        if self.__continuum_maps_raytrace == []:
            # load maps into memory
            file_name = self.__path_results + self.__pname +                   \
                        '_continuum_map_raytrace.fits.gz'
            if os.path.isfile(file_name):
                self.__continuum_maps_raytrace = load_from_fits(file_name)
            else:
                return [] 
        else:
            pass
        return self.__continuum_maps_raytrace * self.__get_unit_value(unit)

    continuum_maps_raytrace = property(return_continuum_maps_raytrace)

    def return_continuum_wlen(self):
        if self.__continuum_wlen == []:

            if os.path.isfile(self.__path_results + self.__pname +             \
                        '_continuum_sed_bin_I.dat'):
                sed_txt = genfromtxt(self.__path_results + self.__pname +      \
                        '_continuum_sed_bin_I.dat', skip_header=0, filling_values="0")
                self.__continuum_wlen = sed_txt[:, 0]
            elif os.path.isfile(self.__path_results + self.__pname +           \
                        '_continuum_sed_mono_I.dat'):
                sed_txt = genfromtxt(self.__path_results + self.__pname +      \
                        '_continuum_sed_mono_I.dat', skip_header=0, filling_values="0")
                self.__continuum_wlen = sed_txt[:, 0]
            elif os.path.isfile(self.__path_results + self.__pname +           \
                        '_continuum_sed_raytrace_I.dat'):
                sed_txt = genfromtxt(self.__path_results + self.__pname +      \
                        '_continuum_sed_raytrace_I.dat', skip_header=0, filling_values="0")
                self.__continuum_wlen = sed_txt[:, 0]
            else:
                print('Warning: Wavelength table not found')
        else:
            pass
        return self.__continuum_wlen
        
    continuum_wlen = property(return_continuum_wlen)
    
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
