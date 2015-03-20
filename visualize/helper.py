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
from numpy.linalg import norm
from os.path import dirname

from datetime import datetime, timedelta

### some global constants are defined here:

h = 6.62606957e-34
c = 299792458.
k = 1.3806488e-23
sigma = 5.6704e-8

const_AU = 1.496e+11
R_sun = 0.6960e+9        # Radius of the sun
L_sun = 3.85e+26
M_sun = 1.9891e+30
parsec = 3.08567758e16          # meters of 1 parsec

CO_lines_freq = np.array([115.2712018, 230.5380000, 345.7959899, 461.0407682,
                          576.2679305, 691.4730763, 806.6518060,
                          921.7997000, 1036.9123930, 1151.9854520,
                          1267.0144860, 1381.9951050, 1496.9229090])*1e9
HCO_lines_freq = np.array([89.1885230, 178.3750650, 267.5576190, 356.7342880,
                           445.9029960, 535.0617755, 624.2086733, 713.3420900,
                           802.4583290, 891.5579242, 980.6374000,
                           1069.6938000])*1e9
C18O_lines_freq = np.array([109.7821734, 219.5603541, 329.3305525, 439.0887658,
                           548.8310055, 658.5532782, 768.2515933, 877.9219553,
                           987.5603822, 1097.1628753, 1206.7254487, 1316.2441143,
                           1425.7148854])*1e9
_13CO_lines_freq = np.array([110.201354, 220.398684, 330.587965, 440.765173,
                             550.92629, 661.06728, 771.184125, 881.272809,
                             991.32931, 1101.349597, 1211.329662, 1321.265481,
                             1431.15304, 1540.98832, 1650.76730])*1e9
CS_lines_freq = np.array([48.9909549, 97.9809533, 146.9690287, 195.9542109,
                          244.9355565, 293.9120865, 342.8828503, 391.8468898,
                          440.8032320, 489.7509210, 538.6889972, 587.6164850,
                          636.5324600, 685.4359238, 734.3259290, 783.2015140,
                          832.0617177])*1e9

HNC_lines_freq = np.array([90.66356800, 181.32475800, 271.98114200, 362.6303030,
                           453.2699220, 543.8975540, 634.5108260, 725.1073410,
                           815.6846760, 906.2404590, 996.7723280, 1087.27785800,
                           1177.75467250, 1268.20038490])*1e9

HCN_lines_freq = np.array([88.6316023, 177.2611115, 265.8864343, 354.5054779,
                           443.1161493, 531.7163479, 620.3040022, 708.8770051,
                           797.4332623, 885.9706949, 974.4871998, 1062.9806890,
                           1151.4490880])*1e9

def get_sensitivity_from_table(mol, instrument='ALMA'):
    """
    get the sensitivity from for a given mol(ecule) and transition
    identifier
    """
    sensitivity = 0
    
    # read all stored sensitivities from file
    if instrument == 'ALMA':
        file_in = open(dirname(__file__)+'/sensitivities_ALMA.dat','r')
    else:
        print('Sensitivities for this (%s) instrument not found' %(instrument))
        return sensitivity

    sensi = file_in.readlines()
    file_in.close

    for entry in sensi:
        if mol == entry.rsplit()[0]:
            sensitivity = float(entry.rsplit()[1])
            break
    return sensitivity

def ca2sp(p_vec):
    """ convert cartesian to spherical coordinates """
    r = np.zeros(3)
    r[0] = norm(p_vec)
    r[1] = np.arctan2(p_vec[2], math.sqrt(p_vec[0]**2 + p_vec[1]**2))
    r[2] = np.arctan2(p_vec[1], p_vec[0])
    return r

def sp2ca(p):
    """ convert spherical to cartesian coordinates """
    r = np.zeros(3)
    r[0] = p[0] * math.cos(p[1]) * math.cos(p[2])
    r[1] = p[0] * math.cos(p[1]) * math.sin(p[2])
    r[2] = p[0] * math.sin(p[1])
    return r

def cy2ca(p):
    r = np.zeros(3)
    r[0] = p[0] * math.cos(p[1]) 
    r[1] = p[0] * math.sin(p[1])
    r[2] = p[2]
    return r

def conv(image, beam=0.1, r_ou=200, dist=140, gg=''):
    """ convolve an image with given beam """
    from numpy import sqrt
    from numpy import log
    from astropy.convolution import convolve_fft, Gaussian2DKernel

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
        z1 = convolve_fft(image, gg, normalize_kernel=True)
    else:
        z1 = convolve_fft(image, gauss, normalize_kernel=True)

    # rescale flux
    # correct from Jy/px to mJy/beam
    z1 *= (beam_x*beam_y/4.)/pxwitdth**2*1000
    return z1


def padwithzeros(vector, pad_width, iaxis, kwargs):
    """ extent a vector """
    vector[:pad_width[0]] = 0
    vector[-pad_width[1]:] = 0
    return vector

def zeropadding(vec, add):
    """ add zeros around an input vector """
    return np.lib.pad(vec, add, padwithzeros)

def get_next_pow_2(x):
    """ calculate the next power of 2 for a given integer """
    return 1<<(x-1).bit_length()

def GetTime(sec_in):
    """
    Converts seconds into a nicer view
    """
    sec = timedelta(seconds=sec_in)
    d = datetime(1, 1, 1) + sec

    print("DAYS:HOURS:MIN:SEC")
    print("%d:%d:%d:%d" % (d.day-1, d.hour, d.minute, d.second))

def write_image2fits(image, header='', file_name='example.fits', quiet=True):
    """
    write image data into a fits file with header
    """
    hdu = pf.PrimaryHDU(image, header=header)
    hdulist = pf.HDUList(hdu)
    if quiet:
        hdulist.writeto(file_name, clobber=True)
    else:
        hdulist.writeto(file_name, clobber=True)
        print(('fits file succesfully written: %s' %file_name))
def beam2pixel_header(header):
    """
    returns the beam2pixel conversion unit using a fits header
    """
    return beam2pixel(header['BMAJ'],header['BMIN'],
                      header['CDELT1'],header['CDELT2'] )

def beam2pixel(bmaj, bmin, delt1, delt2):
    """
    returns the beam2pixel conversion value
    note this is an approximation for an ellipsoide beam
    """

    pixel_size = deg2as(abs(delt1))*deg2as(abs(delt2))
    beam_size = deg2as(abs(bmaj))*deg2as(abs(bmin))*np.pi/4.0/np.log(2)
    beam2px = beam_size/pixel_size
    
    return beam2px
def create_CASA_fits(map_in, out_name='standard.fits',
                     position='4h33m16.50s 22d53m20.40s',
                     object_name='unknown', resolution=1.0e-7,
                     freq=345.8e9, deltafreq=1.0e6):
    """
    generate a fits file, which can be used in casa
    """
    from numpy import zeros
    header = pf.Header()
    map_size = map_in.shape
    if len(map_size) == 4:
        half00 = (map_size[0]-1)/2
        half0 = (map_size[1]-1)/2   # frequency cube
        half1 = (map_size[2]-1)/2
        half2 = (map_size[3]-1)/2
        map2 = map_in

    if len(map_size) == 3:
        half0 = (map_size[0]-1)/2   # frequency cube
        half1 = (map_size[1]-1)/2
        half2 = (map_size[2]-1)/2
        map2 = zeros((1, map_size[0], map_size[1], map_size[2]))
        map2[0, :, :, :] = map_in[:, :, :]

    elif len(map_size) == 2:
        half1 = (map_size[0]-1)/2
        half2 = (map_size[1]-1)/2
        map2 = zeros((1, 1, map_size[0], map_size[1]))
        map2[0, 0, :, :] = map_in[:, :]

    sky = SkyCoord(position)
    header['BUNIT'] = 'JY/PIXEL'
    header['BTYPE'] = 'Intensity'
    header['OBJECT'] = object_name
    header['TELESCOP'] = 'ALMA    '
    header['OBSERVER'] = 'Florian Ober'
    header['IMGTYPE'] = 'ALMA Line radiative transfer simulations'
    header['RADESYS'] = 'ICRS    '
    header['EQUINOX'] = 2000.
    header['EPOCH'] = 2000.0

    header['CTYPE1'] = 'RA---SIN'
    header['CRVAL1'] = sky.ra.deg
    header['CRPIX1'] = float(half1+1)
    header['CDELT1'] = resolution*(-1.0)
    header['CUNIT1'] = 'deg'

    header['CTYPE2'] = 'DEC--SIN'
    header['CRVAL2'] = sky.dec.deg
    header['CRPIX2'] = float(half2+1)
    header['CDELT2'] = resolution
    header['CUNIT2'] = 'deg'

    header['CTYPE3'] = 'FREQ-LSR'
    header['CRVAL3'] = freq
    if len(map_size) == 3:
        # maps for a set of wavelength is provided
        header['CDELT3'] = deltafreq
        header['CRPIX3'] = float(half0+1)

    elif len(map_size) == 2:
        # only one map at one frequency is provided
        header['CDELT3'] = 7.5e9
        header['CRPIX3'] = 1
    header['CUNIT3'] = 'HZ'

    header['CTYPE4'] = 'STOKES'
    header['CRVAL4'] = 1.0
    header['CDELT4'] = 1.0
    header['CRPIX4'] = 1.0
    header['CUNIT4'] = 'Jy'

    header['TIMESYS'] = 'UTC     '
    header['CELLSCAL'] = 'CONSTANT'

    write_image2fits(map2, header=header, file_name=out_name)

def mol3d_intmap(map_in, vch):
    """
    this routine needs some improvement.
    it should produce velocity weighted intensity maps
    """

    chmap = np.zeros((map_in.shape[1], map_in.shape[2]))

    for j in range(map_in.shape[1]):
        for k in range(map_in.shape[2]):
            chmap[j, k] = vch[np.argmax(map_in[:, j, k])]
    return chmap

def calc_L(T=4000, R=2):
    """
    calculate the star Luminosity for given temperature and radius
    """
    return (4.0 * math.pi * (R*R_sun)**2 * sigma * T**4)/L_sun

def calc_R(T=4000, L=1):
    """
    calculate the star radius for given temperature and luminosity
    """
    return ((L*L_sun /(4.0 * math.pi * sigma * T**4))**0.5)/R_sun

def calc_T(R=2, L=1):
    """
    calculate the star temperature for given radius and luminosity
    """
    return (L*L_sun /(4.0 * math.pi * (R*R_sun)**2 * sigma))**0.25

def AU2as(AU, distance=140.0):
    """
    calculate phi (arcsec) of a x AU in 'distance' pc

    1pc = 1AU/1''
    """
    arcs = float(AU) / float(distance)
    return arcs

def as2AU(arcs, distance=140.0):
    """ max resulution to AU """
    AU = float(arcs) * float(distance)
    return AU

def maxres(wave, D):
    """
    Raylight approximation of the angular resolution
    wave := wavelength in m
    D    := diameter of a telescope or baseline of an interferometer in m
    [alpha] = rad
    """

    alpha = 1.22 * wave / D
    return alpha

def res2AU(wave, D, distance):
    """
    get distance in AU for a given wavelength, distance and telescope diameter
    """
    AU = as2AU(rad2as(maxres(wave, D)), distance)
    return AU

def deg2rad(deg):
    """ convert degree to radians """
    rad = math.pi/180.0 * deg
    return rad


def rad2deg(rad):
    """ convert radians to degree """
    deg = 180.0/math.pi * rad
    return deg

def deg2as(deg):
    """ convert degree to arcsec """
    mas = deg*(3600.0)
    return mas

def rad2as(rad):
    """ convert radians to arcsec """
    as_out = deg2as(rad2deg(rad))
    return as_out

def as2deg(mas):
    """ convert arcsec to degree """
    deg = mas/3600.0
    return deg

def as2rad(as_in):
    """ convert arcsec to radians """
    return deg2rad(as2deg(as_in))

def calc_a(b):
    """ calculate coupled alpha """
    a = 3.0*(b-0.5)
    return a

def calc_b(a):
    """ calculate coupled beta """
    b = a/3.0+0.5
    return b

def read_mc3dmap2(path_file, dtype_in='single'):
    """
    read mc3d map and return the map as np array
    version 2, read each line of the ascii file

    return: pic and wavelength
    """
    from numpy import empty, zeros, float32, float64
    from numpy import str as nstr

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

    pic_map = empty((n_map, n_lam_map, n_max_map, n_max_map), dtype=dtypo)

    for maps in range(int(n_map)):
        for lam in range(int(n_lam_map)):
            pic.readline()
            wavelength[lam] = float32(nstr.split(pic.readline())[0])
            pic.readline()
            for i in range(n_max_map):
                for j in range(n_max_map):
                    pic_map[maps, lam, i, j] = float32(nstr.split(pic.readline())[0])

    pic.close()
    return pic_map, wavelength

def read_mc3dmap(path_file, dtype_in='single'):
    """
    read mc3d map and return the map as np array
    version 1, using numpy.loadtxt

    return: pic and wavelength
    """
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

    n_map = pic[0]
    n_lam_map = pic[1]
    n_bin_map = pic[2]
    n_max_map = int(n_bin_map*2.0+1.0)

    start = 5
    counter = start
    wavelength = zeros(int(n_lam_map))

    pic_map = empty((n_map, n_lam_map, n_max_map, n_max_map), dtype=dtypo)
    print(pic_map.nbytes)

    for maps in range(int(n_map)):
        for lam in range(int(n_lam_map)):
            wavelength[lam] = pic[counter]
            for i in range(n_max_map):
                for j in range(n_max_map):
                    counter = counter + 1
                    pic_map[maps, lam, i, j] = pic[counter]

            counter = counter +1
    del pic
    print('calculations took %4.2f seconds' %(time.time() - t1))
    return pic_map, wavelength


def mc3d_EandS(emi, sca):
    """ MC3D: add reemission and scattered light images """

    mapE, wE = read_mc3dmap(emi)
    mapS, wS = read_mc3dmap(sca)

    #add fluxes (emission + scattering)
    map_ges = mapE + mapS

    for w in range(len(wE)):
        if wE[w] != wS[w]:
            print('Error: the wavelength are not the same')
    return map_ges, wE


def plot_image(file_in, num=0, xlab='', ylab='', map_range=[0, 0, 0, 0],
               color_loc=False, cm='', fig_title='', cont=''):
    """
    plot a given 2d image
    """
    from matplotlib.pyplot import imshow, colorbar, figure
    from matplotlib.pyplot import xlabel, ylabel
    from matplotlib.pyplot import title, contour, clabel
    figure(num)

    if fig_title != '':
        title(fig_title)

    if xlab != '':
        xlabel(xlab)

    if ylab != '':
        ylabel(ylab)

    if all([v == 0 for v in map_range]):
        map_range = [0, len(file_in), 0, len(file_in)]

    if cm == '':
        imshow(file_in, extent=map_range, origin='lower', interpolation='None')
    else:
        imshow(file_in, extent=map_range, origin='lower',
               interpolation='None', cmap=cm)
    if color_loc:
        colorbar()
    if cont != '':
        CS = contour(file_in, cont, extent=map_range, origin='lower',
                     hold='on', colors=('k', ))
        clabel(CS, inline=1, fmt='%2.1f', fontsize=10)

def star_rad(wavel, temp=4000, R=2, distance=140):
    """
    compute stars (temperature, radius) blackbody radiation for a
    given input wavelength at a given distance .
    result : flux in Jy
    """
    from numpy import pi
    star = blackbody(wavel, temp)*(R*R_sun)**2 / \
           (distance*parsec)**2*wavel**2/c*pi*1e27
    return star

def blackbody(wavel, temp):
    """
    Compute the blackbody flux for a temperature 'temp' and a (or an array of)
    wavelength. Flux in W/m2/wavelength (m)
    """
    from numpy import exp

    bb = (2.*h*c**2./wavel**5.)*(1./(exp((h*c)/(wavel*k*temp))-1))
    return bb

def blackbody_hz(freq, temp):
    """
    Compute the blackbody flux for a temperature 'temp' and a (or an array of)
    frequency. Flux in W/m2/Hz
    """

    from numpy import exp
    bb = (2.*h*freq**3/c**2.)*(1./(exp((h*freq)/(k*temp))-1))
    return bb

def alma_fwhm(wavelength, max_base=16.0):
    """
    calculate approximate FWHM for a given wavelength and the max baselength
    """
    return 76./(max_base*(c)*1e-9/wavelength)

def mkfft(image):
    """
    make fourier transform of a 2d image
    """
    from numpy.fft import fftshift, ifftshift, fftn
    return fftshift(fftn(ifftshift(image)))

def mkbackfft(image):
    """
    make backwards fourier transform of a 2d image
    """
    from numpy.fft import fftshift, ifftshift, ifftn
    return fftshift(ifftn(ifftshift(image)))


def intrp1D(x1, y1, x2, y2, xcoord):
    interpinx = ((float(x2) - float(xcoord)) * float(y1) + \
                (float(xcoord) - float(x1)) * float(y2)) / \
                (float(x2) - float(x1))
    return interpinx

def poldegang(I, Q, U):
    # for getting polarisation angle and degree at same time
    # operates on entire arrays -> much faster
    # based on P. Scicluna, 2015 routine

    poldeg = np.zeros_like(I)
    gamma = np.zeros_like(I)
    
    ind = np.where(I > 0.)
    poldeg[ind] = (np.sqrt(Q[ind]*Q[ind] + U[ind]*U[ind]))/I[ind]

    ind = np.where(poldeg > 0.)
    gamma[ind] = 0.5*np.arctan2(U[ind], Q[ind])
    
    return poldeg, gamma

def shrink(data, rows, cols):
    return data.reshape(rows, data.shape[0]/rows,
                        cols, data.shape[1]/cols).sum(axis=1).sum(axis=2)
