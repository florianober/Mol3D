#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
module which contains serveral little helper routines


date : 05/2014
author : F. Ober
email: fober@astrophysik.uni-kiel.de
"""


import math
import numpy as np
from astropy.io import fits as pf
from astropy.coordinates import SkyCoord

from datetime import datetime, timedelta

### some global constants are defined here:

h = 6.62606957e-34
c = 299792458.
k = 1.3806488e-23
sigma = 5.6704e-8

AU    = 1.496e+11
R_sun = 0.6960e+9        # Radius of the sun
L_sun = 3.85e+26
M_sun = 1.9891e+30
parsec = 3.08567758e16          # meters of 1 parsec

CO_lines_freq =  np.array([115.2712018,230.5380000,345.7959899,461.0407682,576.2679305,691.4730763,806.6518060,
                 921.7997000,1036.9123930,1151.9854520,1267.0144860,1381.9951050,1496.9229090])*1e9
HCO_lines_freq = np.array([89.1885230,178.3750650,267.5576190,356.7342880,445.9029960,535.0617755,
                        624.2086733,713.3420900,802.4583290,891.5579242,980.6374000,1069.6938000])*1e9

def conv(image, beam=0.1,r_ou=200,dist=140,gg=''):
    from numpy import sqrt 
    from numpy import log
    from astropy.convolution import convolve_fft, Gaussian2DKernel
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
        
    # rescale flux
    # correct from Jy/px to mJy/beam
    z1 *= (beam_x*beam_y/4.)/pxwitdth**2*1000
    return z1


def padwithzeros(vector, pad_width, iaxis, kwargs):
    vector[:pad_width[0]] = 0
    vector[-pad_width[1]:] = 0
    return vector

def zeropadding(vec,add):
    return np.lib.pad(vec,add,padwithzeros)

def get_next_pow_2(x):
    return 1<<(x-1).bit_length()

def GetTime(sec_in):
    sec = timedelta(seconds=sec_in)
    d = datetime(1,1,1) + sec

    print("DAYS:HOURS:MIN:SEC")
    print("%d:%d:%d:%d" % (d.day-1, d.hour, d.minute, d.second))

def write_image2fits(image,header='',file_name='example.fits',quiet=True):
    hdu = pf.PrimaryHDU(image,header=header)
    hdulist = pf.HDUList(hdu)
    if quiet:
        hdulist.writeto(file_name,clobber=True)
    else:
        hdulist.writeto(file_name,clobber=True)
        print(('fits file succesfully written: %s' %file_name)) 

def create_CASA_fits(map_in,out_name='standard.fits',position='4h33m16.50s 22d53m20.40s',
                     object_name='unknown',resolution=1.0e-7,freq=345.8e9,deltafreq=1.0e6):
    from numpy import zeros
    header = pf.Header()
    map_size = map_in.shape
    if len(map_size) == 4:
        half00   = (map_size[0]-1)/2
        half0    = (map_size[1]-1)/2   # frequency cube
        half1    = (map_size[2]-1)/2
        half2    = (map_size[3]-1)/2
        map2     = map_in
    
    if len(map_size) == 3:
        half0    = (map_size[0]-1)/2   # frequency cube
        half1    = (map_size[1]-1)/2
        half2    = (map_size[2]-1)/2
        map2 = zeros((1,map_size[0],map_size[1],map_size[2]))
        map2[0,:,:,:] = map_in[:,:,:]
        
    elif len(map_size) == 2:
        half1    = (map_size[0]-1)/2
        half2    = (map_size[1]-1)/2
        map2 = zeros((1,1,map_size[0],map_size[1]))
        map2[0,0,:,:] = map_in[:,:]

    c = SkyCoord(position)
    
    header['BUNIT']    = 'JY/PIXEL'
    header['BTYPE']    = 'Intensity'  
    header['OBJECT']   = object_name
    header['TELESCOP'] = 'ALMA    '
    header['OBSERVER'] = 'Florian Ober'
    header['IMGTYPE']  = 'ALMA Line radiative transfer simulations'          
    header['RADESYS'] = 'ICRS    '
    header['EQUINOX'] =  2000.
    header['EPOCH']   =  2000.0 
    
    header['CTYPE1']  = 'RA---SIN'
    header['CRVAL1']  = c.ra.deg
    header['CRPIX1']  = float(half1+1)
    header['CDELT1']  = resolution*(-1.0)
    header['CUNIT1']  = 'deg'   

    header['CTYPE2']  = 'DEC--SIN'
    header['CRVAL2']  = c.dec.deg
    header['CRPIX2']  = float(half2+1)
    header['CDELT2']  = resolution
    header['CUNIT2']  = 'deg' 
    
    header['CTYPE3']  = 'FREQ-LSR'
    header['CRVAL3']  = freq
    if len(map_size) == 3:
        # maps for a set of wavelength is provided
        header['CDELT3']  = deltafreq
        header['CRPIX3']  = float(half0+1)
        
    elif len(map_size) == 2:
        # only one map at one frequency is provided
        header['CDELT3']  = 7.5e9
        header['CRPIX3']  = 1
    header['CUNIT3']  = 'HZ'
        
    header['CTYPE4']  = 'STOKES'
    header['CRVAL4']  = 1.0
    header['CDELT4']  = 1.0
    header['CRPIX4']  = 1.0
    header['CUNIT4']  = 'Jy'
    
    header['TIMESYS'] = 'UTC     '                                                            
    header['CELLSCAL']  = 'CONSTANT'
    
    write_image2fits(map2,header=header,file_name=out_name)
    
def mol3d_intmap(map_in,vch):
    # this routine needs some improvment.
    # it should produce velocity weighted intensity maps
    chmap = np.zeros((map_in.shape[1],map_in.shape[2]))
    
    for j in range(map_in.shape[1]):
        for k in range(map_in.shape[2]):
            chmap[j,k] = vch[np.argmax(map_in[:,j,k])]
    return chmap

def calc_L(T=4000,R=2):
    return (4.0 * math.pi * (R*R_sun)**2 * sigma * T**4)/L_sun
    
def calc_R(T=4000,L=1):
    
    return ((L*L_sun /(4.0 * math.pi * sigma * T**4))**0.5)/R_sun
    
def calc_T(R=2,L=1):
    
    return (L*L_sun /(4.0 * math.pi * (R*R_sun)**2 * sigma ))**0.25
    
def load_mol3d_zerovchmap(file_path,ch=-1):
        map_in = open(file_path)
        row = np.str.split(map_in.readline())
        vel_ch = int(row[0])+1

        row = np.str.split(map_in.readline())
        map_size = int(row[0])
        pic = np.zeros((map_size,map_size))
        row = map_in.readline()	#empty row
        if ch != -1:
            map_ch = ch
        else:
            map_ch = (vel_ch-2)/2
        for k in range(vel_ch-1):
            row = map_in.readline()	#empty row
            row = np.str.split(map_in.readline())
            ch_val = float(row[0])
            row = map_in.readline()	#empty row
            for j in range(map_size**2):
                row = np.str.split(map_in.readline())
                pic[int(row[1]),int(row[0])] = row[2]
            if k == map_ch:
                print('loaded %6.3f [km/s] channel map' %(ch_val*0.001))
                break
        return pic
        
def load_mol3d_fullvchmap(file_path):
        map_in = open(file_path)
        row = np.str.split(map_in.readline())
        vel_ch = int(row[0])
        row = np.str.split(map_in.readline())
        map_size = int(row[0])
        
        pic = np.zeros((vel_ch,map_size,map_size))
        vch = np.zeros(vel_ch)
        
        row = map_in.readline()	#empty row
        for k in range(vel_ch):
            row = map_in.readline()	#empty row
            row = np.str.split(map_in.readline())
            vch[k] = float(row[0])
            row = map_in.readline()	#empty row
            for j in range(map_size**2):
                row = np.str.split(map_in.readline())
                pic[k,int(row[1]),int(row[0])] = row[2]
        return pic,vch


def load_mol3d_map(file_path):

    map_in = open(file_path)
    row = np.str.split(map_in.readline())
    map_size = int(row[0])
    row = map_in.readline()	#empty row
    pic = np.zeros((map_size,map_size))
    #~ print map_size
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
    return pic

def AU2as(AU,distance=140.0):
    # calculate phi (arcsec) of a x AU in 'distance' pc
    # 1pc = 1AU/1''
    arcs = float(AU) / float(distance)
    return arcs
    
def as2AU(arcs,distance=140.0):   
    AU = float(arcs) * float(distance)
    return AU
    
def maxres(wave,D):
    '''Raylight approximation of the angular resolution
       wave := wavelength in m
       D    := diameter of a telescope or baseline of an interferometer in m
       [alpha] = rad
    '''
    alpha = 1.22*wave/D
    return alpha

def res2AU(wave,D,distance):
    
    AU = as2AU(rad2as(maxres(wave,D)),distance)
    return AU
    
def deg2rad(deg):
    
    rad = math.pi/180.0 *deg
    return rad
    

def rad2deg(rad):
    
    deg = 180.0/math.pi * rad
    return deg

def rad2as(rad):
    as_out = deg2as(rad2deg(rad))
    return as_out
    
def deg2as(deg):
    
    mas = deg*(3600.0)
    return mas

def as2deg(mas):
    
    deg = mas/3600.0
    return deg
    
def as2rad(as_in):
    return deg2rad(as2deg(as_in))
    
    
def calc_a(b):
    
    a = 3.0*(b-0.5)
    return a    

def calc_b(a):
    
    b = a/3.0+0.5
    return b

def read_mc3dmap2(path_file,dtype_in='single'):
    from numpy import empty, zeros, loadtxt, float32, float64
    from numpy import str as nstr
    import time
    
    #~ t1 = time.time()
    
    if dtype_in == 'single':
        dtypo = float32
    elif dtype_in == 'double':
        dtypo = float64
    else:
        print('Warning: dtype not supported')
        dtypo = float32
    pic = open(path_file, 'r')
    
    pic.readline()
    pic.readline()
    
    n_map = int(nstr.split(pic.readline())[0])
    n_lam_map = int(nstr.split(pic.readline())[0])
    n_bin_map = int(nstr.split(pic.readline())[0])
    n_max_map = int(n_bin_map*2.0+1.0)
    
    pic.readline()   
    pic.readline()   
    pic.readline()   
    pic.readline()   
    pic.readline()
    
    wavelength = zeros(int(n_lam_map))

    pic_map = empty((n_map,n_lam_map,n_max_map,n_max_map),dtype=dtypo)
    
    for maps in range(int(n_map)):
        for lam in range(int(n_lam_map)):
            pic.readline()
            wavelength[lam] = float32(nstr.split(pic.readline())[0])
            pic.readline()
            for i in range(n_max_map):
                for j in range(n_max_map):
                    pic_map[maps,lam,i,j] = float32(nstr.split(pic.readline())[0])


    pic.close()
    #~ print 'calculations took %4.2f seconds' %(time.time()-t1)
    return pic_map,wavelength


def read_mc3dmap(path_file,dtype_in='single'):
    # read mc3d map and return the map as np array
    # return pic and wavelength
    from numpy import empty, zeros, loadtxt, float32, float64
    import time
    
    t1 = time.time()
    if dtype_in == 'single':
        dtypo = float32
    elif dtype_in == 'double':
        dtypo = float64
    else:
        print('Warning: dtype not supported')
        dtypo = float32
    
    pic = loadtxt(path_file, dtype=dtypo)
    print('loaded file')

    n_map     = pic[0]
    n_lam_map = pic[1]
    n_bin_map = pic[2]
    n_max_map = int(n_bin_map*2.0+1.0)

    start     = 5
    counter = start
    wavelength = zeros(int(n_lam_map))

    pic_map = empty((n_map,n_lam_map,n_max_map,n_max_map),dtype=dtypo)
    print(pic_map.nbytes)

    for maps in range(int(n_map)):
        for lam in range(int(n_lam_map)):
            wavelength[lam] = pic[counter]
            for i in range(n_max_map):
                for j in range(n_max_map):
                    counter = counter +1
                    pic_map[maps,lam,i,j] = pic[counter]

            counter = counter +1
    del pic
    print('calculations took %4.2f seconds' %(time.time()-t1))
    return pic_map,wavelength
    
    
def mc3d_EandS(emi,sca):
    
    #calculate and 
    
    mapE,wE = read_mc3dmap(emi)
    mapS,wS = read_mc3dmap(sca)
    
    #add fluxes (emission + scattering)
    map_ges = mapE + mapS
    
    for w in range(len(wE)):
        if wE[w] != wS[w]:
            print('Error: the Wavelength are not the same')
    return map_ges,wE
    
    
def plot_image(file_in,num=0,xlab='',ylab='',map_range=[0,0,0,0],color_loc=False,cm='',fig_title='',cont=''):
    from matplotlib.pyplot import imshow, colorbar, figure, xlabel, ylabel, xlim, ylim, title, contour, clabel
    import scipy.ndimage
    figure(num)
    
    if fig_title != '':
        title(fig_title)
    
    if xlab != '':
        xlabel(xlab)
    
    if ylab != '':
        ylabel(ylab)
        
    if all([ v == 0 for v in map_range ]):
        map_range = [0,len(file_in),0,len(file_in)]
    
    if cm == '':
        imshow(file_in,extent=map_range,origin='lower',interpolation='None')
    else:
        imshow(file_in,extent=map_range,origin='lower',interpolation='None',cmap=cm)
    if color_loc:
        colorbar()
    if cont != '':
        #file_in = scipy.ndimage.zoom(file_in, 5)
        CS = contour(file_in,cont,extent=map_range,origin='lower',hold='on',colors=('k',))
        clabel(CS,inline=1,fmt='%2.1f', fontsize=10)
        
        
        
def star_rad(wavel,temp=4000,R=2,distance=140):
    # compute stars (temperature, radius) blackbody radiation for a given input wavelength
    # at a given distance 
    # result : flux in Jy
    
    from numpy import pi
    star = blackbody(wavel,temp)*(R*R_sun)**2/(distance*parsec)**2*wavel**2/c*pi*1e27
    return star
    
def blackbody(wavel,temp):
    
    from numpy import exp

    #-----------------------------------
    # Compute the blackbody flux for a temperature 'temp' and a (or an array of)
    # wavelength 'lambda'. Flux in W/m2/wavelength (m)
    #-----------------------------------
    
    bb = (2.*h*c**2./wavel**5.)*(1./(exp((h*c)/(wavel*k*temp))-1)) #=> EN W.m-2.m-1
    
    return bb
    
def blackbody_hz(freq,temp):
    
    from numpy import exp

    #-----------------------------------
    # Compute the blackbody flux for a temperature 'temp' and a (or an array of)
    # frequency. Flux in W/m2/Hz
    #-----------------------------------

    bb = (2.*h*freq**3/c**2.)*(1./(exp((h*freq)/(k*temp))-1)) 
    return bb
    
def ALMA_FWHM(wavelength,max_base=16.0):
    # calculate approximate FWHM for a given wavelength and the max baselength
    return 76./(max_base*(c)*1e-9/wavelength)

def mkfft(image):
    from numpy.fft import fftshift, ifftshift, fftn
    return fftshift(fftn(ifftshift(image)))

def mkbackfft(image):
    from numpy.fft import fftshift, ifftshift, fftn
    return fftshift(ifftn(ifftshift(image)))

    
    
def intrp1D(x1,y1,x2,y2,xcoord):
    interpInX = ((float(x2)-float(xcoord))*float(y1)+(float(xcoord)-float(x1))*float(y2))/(float(x2)-float(x1))
    return interpInX

def get_uv_coords(array_in, station1, station2, ra, dec, ha):
    from numpy import array,inner,cross,sqrt,sum
    
    # created by Paul Boley (boley@mpia-hd.mpg.de)
    
    # Convert hour angle, ra and dec to radians
    ha %= 24.0
    ha *= math.pi / 12.0
    center = array_in.arrxyz
    ra *= math.pi / 180.0
    dec *= math.pi / 180.0
    # Same for lat/long of array center
    latitude = array_in.latitude * math.pi / 180.0
    longitude = array_in.longitude * math.pi / 180.0

    # Calculate an "east" unit vector
    east = array([center[1]*center[2],-center[0]*center[2],0])
    if center[2] > 0: east *= -1
    east /= sqrt(sum(east**2))
    # Calculate a "north" unit vector
    north = cross(center, east)
    north /= sqrt(sum(north**2))
    up = center / sqrt(sum(center**2))

    # Eq. 3 of Segransan 2007
    B = array([inner(station2.staxyz - station1.staxyz, north),
                  inner(station2.staxyz - station1.staxyz, east),
                  inner(station2.staxyz - station1.staxyz, up)])

    # From Lawson's book
    u = B[1] * math.cos(ha) - B[0] * math.sin(latitude) * math.sin(ha)
    v = B[1] * math.sin(dec) * math.sin(ha) + B[0] * (math.sin(latitude) * math.sin(dec) * math.cos(ha) + math.cos(latitude) * math.cos(dec))

    return u, v
