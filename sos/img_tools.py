# -*- coding: utf-8 -*-
# --------------------------------------------------------------------------------- #
# Software de Observaciones Sintéticas S.O.S.
# Image tools
#
# Marcial Becerril, @ 28 June 2021
# Latest Revision: 28 Jun 2021, 21:55 GMT-6
#
# For all kind of problems, requests of enhancements and bug reports, please
# write to me at:
#
# mbecerrilt92@gmail.com
# mbecerrilt@inaoep.mx
#
# --------------------------------------------------------------------------------- #

import os
import numpy as np

from astropy import wcs
from astropy.io import fits
from astropy.coordinates import SkyCoord

from .misc.print_msg import *


def load_fits(path):
    """
        Loading fits file
        Parameters
        ----------
        path : string
            Path file of fits file
        ----------
    """
    # Read FITS file
    hdulist = fits.open(path)
    hdu = hdulist[0]

    # Get header
    header = hdu.header
    # Get data
    data = hdu.data

    return data, header


def move_operations(data, inst, **kwargs):
    """
        Matrix move operations: Rotation to right/left and
        transpose
        Parameters
        ----------
        data : np.array
            Matriz to rotate
        inst : string
            Instructions:
            'r': rotate to the right
            'l': rotate to the left
            't': transpose
            The commands must be separeted by '-'
        ----------
    """

    # Get instructions
    chain = inst.split('-')

    # Apply operations
    for c in chain:
        if c == 'r':
            data = np.flip(data.T,1)
        elif c == 'l':
            data = np.flip(A.T,0)
        elif c == 't':
            data = data.T
        else:
            msg('Instruction not recognized', 'warn')

    return data


def trim_data_cube(path, center, dims, **kwargs):
    """
        Trim data cube
        Parameters
        ----------
        path : string
            Path FITS file
        center : tuple
            Image center
        dims : tuple
            Cut dimensions
            w : width [degrees]
            h : heigth [degrees]
        ----------
    """

    # Key arguments
    # ----------------------------------------------
    # File name
    vel = kwargs.pop('vel', None)
    # ----------------------------------------------

    # Trim center
    RA = center[0]
    DEC = center[1]

    # Coordinates
    coor = SkyCoord(RA, DEC, frame='fk4')
    gal = coor.galactic

    l = gal.l.value 
    b = gal.b.value

    # Read the file
    data, header = load_fits(path)

    # Cut Dimensions
    w, h = dims

    # Headers parameters
    ctype = header['CTYPE*']

    cdelt = np.zeros(len(ctype))
    crpix = np.zeros(len(ctype))
    crval = np.zeros(len(ctype))
    naxis = np.zeros(len(ctype), dtype=int)

    for i in ctype:
        f = header[i].lower()
        if ('lon' in f) or ('ra' in f) or ('x' in f):
            idx = 2
        elif ('lat' in f) or ('dec' in f) or ('y' in f):
            idx = 1
        elif ('vel' in f):
            idx = 0

        c = i[5:]

        cdelt[idx] = header['CDELT'+c]
        crpix[idx] = header['CRPIX'+c]
        crval[idx] = header['CRVAL'+c]
        naxis[idx] = header['NAXIS'+c]

    # Conversion from pixels to galactic coordinates
    lon = (1 + np.arange(naxis[2]) - crpix[2])*cdelt[2] + crval[2]
    lat = (1 + np.arange(naxis[1]) - crpix[1])*cdelt[1] + crval[1]

    xmin = l - w/2.
    xmax = l + w/2.
    xlim = np.where((lon >= xmin) & (lon <= xmax))[0]

    ymin = b - h/2.
    ymax = b + h/2.
    ylim = np.where((lat >= ymin) & (lat <= ymax))[0]

    new_data = np.zeros((naxis[0], len(ylim), len(xlim)))
    for s in range(naxis[0]):
        for k, j in enumerate(ylim):
            for h, i in enumerate(xlim):
                new_data[s, k, h] = data[s, j, i]

    header['CRVAL1'] = -1*xlim[0]*np.abs(cdelt[1])
    header['CRVAL2'] = -1*ylim[0]*np.abs(cdelt[2])

    # Check for the velocity dimension
    if vel:
        if len(naxis) > 2:
            if isinstance(vel, float):
                v_data = new_data[vel,:,:]
            elif isinstance(vel, tuple):
                v_data = new_data[vel[0]:vel[1],:,:]
            else:
                v_data = new_data

            return v_data, header

    return new_data, header


def save_img(data, header, **kwargs):


    # Key arguments
    # ----------------------------------------------
    # File name
    name = kwargs.pop('name', 'img_trim')
    # ----------------------------------------------

    # Save as fits file
    hdu = fits.PrimaryHDU(data=data, header=header)

    # Replace file if it already exist
    os.system('rm -rf '+name+'.fits')
    hdu.writeto(name+'.fits')



# Example:
#
# >> path = '/home/marcial/Documentos/Complex S235/s235_13co.fits'
#
# Coordinates definition
# >> RA = '05h37m44.3s'
# >> DEC = '+35d48m53s'
#
# Cut dimensions
# >> w = 1.0
# >> h = 1.0
#
# Trimming
# >> data, header = trim_data_cube(path, (RA, DEC), (w, h))
# >> save_img(data, header, name='s235_co_trim')
