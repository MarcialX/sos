# -*- coding: utf-8 -*-
# --------------------------------------------------------------------------------- #
# Software de Observaciones Sintéticas S.O.S.
# Molecular cloud data
#
# Marcial Becerril, @ 24 August 2020
# Latest Revision: 15 Sep 2021, 01:00 GMT-6
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

from scipy.optimize import curve_fit
from scipy.integrate import simps
from scipy.signal import find_peaks
from scipy.ndimage import gaussian_filter

from astropy import wcs
from astropy.io import fits
from astropy.wcs import WCS

import sos
from datetime import datetime

from .physics import *
from .misc.print_msg import *
from .misc.units_tool import *
from .mc_db_tools import *
from .line_fitting import *
from .init_vals import *
from .misc.constants import *
from .visual_tools.param_visual import *


class mc(object):
    """
        Molecular cloud object
        Parameters
        ----------
        mcID : string
            Molecular cloud ID
        ----------
    """
    def __init__(self, mcID, *args, **kwargs):

        # Key arguments
        # ----------------------------------------------
        # Index of the cube in case it has more than four dimensions
        idx_4d = kwargs.pop('idx_4d', 0)
        # Volume model
        vol_shape = kwargs.pop('vol_shape', 'sphere')
        # Move operations for polarization
        move_pol = kwargs.pop('move_pol', None)
        # Move operations for data cubes
        move_mol = kwargs.pop('move_mol', None)
        # ----------------------------------------------

        # Load molecular cloud database
        self.mc_data = mc_db(sos.DB_PATH)

        # If the cloud is an empty file, force to load 
        self._flag_emptyMC = False
        if mcID in ['', 'zeroMC']:
            # Load default params
            mcID = 'zeroMC'
            self._flag_emptyMC = True

            # Dimensions of empty maps in pixels
            dims = kwargs.pop('dims', (100, 100, 10))
            # Velocity dims
            vel_lims = kwargs.pop('vel_lims', (-1e3, 1e3))
            # Longitud dims
            lon_lims = kwargs.pop('lon_lims', (-1, 1))
            # Longitud dims
            lat_lims = kwargs.pop('lat_lims', (-1, 1))

            if len(dims) != 3:
                msg('Map dimensions have to be equals to 3', 'fail')
                return

        # Loading the molecular cloud files
        self.mcID = mcID
        msg('Loading object parameters', 'info')

        try:
            self.mc_params = self.mc_data.get_mc_params(self.mcID)
            msg('Parameters loaded', 'ok')

        except Exception as e:
            msg('Loading object parameters. Check the parameters and try again.'+str(e), 'fail')
            return

        # Radius of the molecular cloud
        kRe = self.mc_params['Re_factor']
        self.R = kRe*self.mc_params['dist']*np.tan(0.5*self.mc_params['ang_dia']*np.pi/180.)

        # Read molecules available
        #self.mols = sos.thin_mols + sos.thick_mols
        self.mols = []
        for m in mol_params.keys():
            if m in self.mc_params.keys():
                if not self.mc_params[m] is None:
                    self.mols.append(m)

        self.thin = ''
        self.thick = ''

        # Get default thin and thick molecules
        p = np.zeros(len(self.mols))
        for i, n in enumerate(self.mols):
            p[i] = mol_params[n]['P']

        p = np.argsort(p)

        if len(p) >= 2:
            thick_idx = np.where(p==0)[0][0]
            self.thick = self.mols[thick_idx]
            thin_idx = np.where(p==1)[0][0]
            self.thin = self.mols[thin_idx]
        elif len(p) == 1:
            self.thick = self.mols[0]

        # Header and Data of spectroscopy data
        self.Dmol, self.Hmol = {}, {}
        self.noise, self.noise_stat = {}, {}
        for mol in self.mols:
            self.Dmol[mol], self.Hmol[mol] = [], {}
            self.noise[mol], self.noise_stat[mol] = {}, {}
        # Header and Data of polarization data
        self.Dpol, self.Hpol = [], {}

        # General propierties
        # -----------------------------------------------------------
        # The whole molecular cloud
        self.full = {}

        # Baseline definition
        self.baseline = {}
        # Molecule lines initialization
        for mol in self.mols:
            self.full[mol] = {}
            self.full[mol]['T'] = 0.0
            self.full[mol]['vel'] = np.array([])
            self.full[mol]['line'] = np.array([])
            self.full[mol]['freq'] = 0.0
            self.baseline[mol] = np.array([])

        # Polarization
        self.full['I'] = []
        self.full['Q'] = []
        self.full['U'] = []
        self.full['V'] = []
        # Angle and degree of polarization map
        self.full['pol_angle'] = []
        self.full['pol_degree'] = []
        self.full['pol_vector'] = []

        # MC center
        self.xc, self.yc = 0., 0.

        # Binning grid
        self.bin_grid = np.array([])
        # Binning frame 
        self.bin_frames = []

        # Binned molecular cloud
        # -----------------------------------------------------------
        # Initiate the binned dictionary
        self.binned = {}
        self.binned['B0'] = {}
        self.binned['B0']['flag'] = False

        # Boundaries for each bin
        self.bin_grid = np.array([])
        for mol in self.mols:
            self._init_binned_data(mol)
            msg(mol+' binned initiated', 'ok')

        self.nbins = 1

        # Load molecule cloud
        # Empty cloud
        if self._flag_emptyMC:
            # Load molecules data
            for mol in self.mols:
                self.load_empty_mol(mol, dims, vel_lims, lon_lims, lat_lims)
        # Defined cloud
        else:
            # Load molecules data
            for mol in self.mols:
                self.load_molecule(mol, idx=idx_4d, move_cmd=move_mol)
            # Load polarization map
            path_pol = self.mc_params['stokes']
            self.load_polarization(path_pol, move_cmd=move_pol)

        # Calculate volume
        try:
            self.get_full_volume(model=vol_shape, update_bins=False)
            msg('Volume calculated', 'ok')
        
        except Exception as e:
            msg('Volume error calculation\n'+str(e), 'warn')


    def get_full_volume(self, model, update_bins=True, **kwargs):
        """
            Get the volume of the full cloud
            Parameters
            ----------
            model : string
                Volume model: 'sphere' or 'mc_disk'
            update_bins : boolean
                Update the bins values, if they exist 
            ----------
        """
        # Get dimensions
        flag_dims = True
        aux_xp = 0
        aux_yp = 0
        for i, mol in enumerate(self.mols):
            xp = self.Hmol[mol]['NAXIS1']
            yp = self.Hmol[mol]['NAXIS2']
            if (xp == aux_xp and yp == aux_yp) or i == 0:
                aux_xp = xp
                aux_yp = yp
            else:
                flag_dims = False
                break

        if flag_dims:
            x = self.Hmol[self.mols[0]]['NAXIS1'] 
            y = self.Hmol[self.mols[0]]['NAXIS2']
        else:
            msg('Images have not the same spatial size', 'fail')
            return

        # Define volume array
        Vs = np.zeros((x, y))

        # Define dimensions of each pixel
        delta_lon = np.abs(np.mean(np.diff(self.lon_m)))
        delta_lat = np.abs(np.mean(np.diff(self.lat_m)))
        # Get the volume per pixel
        for i in range(x):
            for j in range(y):
                # If MC is an empty file volume is zero
                if self._flag_emptyMC:
                    V = 0.0
                else:
                    # X and Y limits
                    xlims = (self.lon_m[i]-delta_lon/2, self.lon_m[i]+delta_lon/2)
                    ylims = (self.lat_m[j]-delta_lat/2, self.lat_m[j]+delta_lat/2)
                    # Volume
                    V, corners = model_volume(xlims, ylims, self.R, shape=model, **kwargs)
                # Asign volume
                Vs[i][j] = V

        self.full['vol'] = Vs

        msg('Volume updated', 'ok')

        if update_bins:
            if len(self.bin_grid) > 0 and len(self.bin_frames) > 0:
                # Asign volume to bins
                cnt = 0
                for b in self.binned.keys():
                    # Get bounds and frames
                    # Get bin limits
                    xmin, xmax, ymin, ymax = self.bin_grid[cnt]
                    # Get frames
                    frames_pts = self.bin_frames[cnt]
                    # Get fractional volume
                    if len(frames_pts) > 0:
                        add_vol_pix = self._int_frac_pixels(frames_pts, Vs)
                    else:
                        add_vol_pix = 0.
                    # Asign volume per bin
                    self.binned[b]['vol'] = np.sum(Vs[xmin:xmax, ymin:ymax])+add_vol_pix

                    cnt += 1

                msg('Volume binning updated', 'ok')

            else:
                msg('Bin bounds and frames not defined yet', 'warn')


    def load_fits(self, path, idx_4d=0, cube=True):
        """
            Loading fits file
            Parameters
            ----------
            path : string
                Path file of fits file
            cube : boolean
                Is it a FITS data cube (line observations)?
            ----------
        """
        # Read FITS file
        hdulist = fits.open(path)
        hdu = hdulist[0]

        # Get header
        header = hdu.header
        n_header = header.copy()
        # Get data
        data = hdu.data

        # Check if the data is cube
        dims = data.shape
        if len(dims) > 3:
            cnt_ones = 0
            for i in dims:
                if i == 1:
                    cnt_ones = cnt_ones + 1

            if (len(dims) - cnt_ones) > 4:
                raise Exception("There are more data cubes than indexes to choose it. Check the fits data") 

        # Get number of axes
        n_ax = header['NAXIS']

        # Header files to update
        h_fields = ['CTYPE', 'CDELT', 'CUNIT', 'CRPIX', 'CRVAL', 'NAXIS']
        # Clear the h_fields
        for h in h_fields:
            del(n_header[h+'*'])

        dtype = header['CTYPE*']

        if cube:
            # Axis names [Spectroscopy]
            ax_names = { 'x': ['ra', 'lon', 'x'],
                         'y': ['dec', 'lat', 'y'],
                         'vel': ['vel', 'v']
            }
        else:
            # Axis names [Polarimetry]
            ax_names = { 'x': ['ra-', 'lon', 'x'],
                         'y': ['dec-', 'lat', 'y'],
                         'pol': ['stoke', 'pol']
            }
            h_fields_pol = ['CTYPE', 'CUNIT', 'NAXIS']

        ndim = 3
        new_shape = [0]*ndim

        for i, k in enumerate(ax_names.keys()):
            # Check each of the axis
            for j, dt in enumerate(dtype):
                ax = header[dt]
                # Check where axis it belongs
                if self._check_axes(ax, ax_names[k]):
                    name_axes = dt[5:]
                    num_axes = ""
                    for na in name_axes:
                        if na.isnumeric():
                            num_axes += na
                    new_shape[i] = header['NAXIS'+num_axes]
                    if k == 'pol':
                        hdrs_data = h_fields_pol
                    else:
                        hdrs_data = h_fields
                    # Update the header
                    for h in hdrs_data:
                        if h == 'NAXIS':
                            hd_nm = h+num_axes
                        else:
                            hd_nm = h+name_axes
                        n_header[h+str(i+1)] = header[hd_nm]
                    break

        new_shape = new_shape[::-1]

        if len(dims) > 3:
            aux_shape = list(dims)
            for i,dim in enumerate(aux_shape):
                if not dim in new_shape: 
                    aux_shape[i] = idx_4d + 1
            # Selecting the cube
            data = data[0:aux_shape[0], 0:aux_shape[1], 0:aux_shape[2], 0:aux_shape[3]]

        data = data.reshape(tuple(new_shape))

        return data, n_header


    def load_molecule(self, mol, idx, **kwargs):
        """
            Loading Molecule data
            Parameters
            ----------
            mol : string
                Molecule
            idx : int
                4 dimensions index
            ----------
        """

        # Key arguments
        # ----------------------------------------------
        # Figure size
        move_cmd = kwargs.pop('move_cmd', None)
        # ----------------------------------------------

        # Get molecule path
        path_mol = self.mc_params[mol]
        # Get effective radius of the object
        kRe = self.mc_params['Re_factor']

        msg('Loading '+mol+' map', 'info')
        try:
            self.Dmol[mol], self.Hmol[mol] = self.load_fits(path_mol, idx_4d=idx)
            # Transformation operations
            if move_cmd:
                for s in range(self.Dmol[mol].shape[0]):
                    self.Dmol[mol][s] = move_operations(self.Dmol[mol][s], move_cmd) 
            # Get units factors
            ks = []
            uts = []
            us = self.Hmol[mol]['CUNIT*']
            for u in us:
                unit = us[u].lower()

                # Fill power symbols in case they are missing
                # -------------------------------------------
                cnt = 0
                unit_txt = ''
                flag_pow = False
                for ux in unit:
                    if ux.isalpha():
                        flag_pow = True
                    elif ux in ['-', '+'] or ux.isalnum():
                        if flag_pow:
                            unit_txt += '^'
                            flag_pow = False
                    else:
                        flag_pow = False

                    unit_txt += ux
                unit = unit_txt
                # -------------------------------------------

                if unit:
                    k, ut = get_factors_units(unit, su) 
                    s_ut = unit2str(ut)
                    ks.append(k)
                    uts.append(s_ut)
                    # Update all CUNITS
                    self.Hmol[mol][u] = s_ut
                    self.Hmol[mol]['CDELT'+u[-1]] = k*self.Hmol[mol]['CDELT'+u[-1]] 
                    #self.Hmol[mol]['CRPIX'+u[-1]] = k*self.Hmol[mol]['CRPIX'+u[-1]]
                    self.Hmol[mol]['CRVAL'+u[-1]] = k*self.Hmol[mol]['CRVAL'+u[-1]]

            # Rest frequency
            for rf in ['RESTFREQ', 'RESTFRQ']:
                if rf in self.Hmol[mol]:
                    self.full[mol]['freq'] = self.Hmol[mol][rf]
                    break

            # Get line transition
            if 'LINE' in self.Hmol[mol]:
                self.full[mol]['t_line'] = self.Hmol[mol]['LINE'].replace(mol, '')
            else:
                # In case the line is not defined, the transition J=1-0 is taken as the default
                self.full[mol]['t_line'] = '1-0'

            # Check if the frequency is available in the file
            if self.full[mol]['freq'] == 0.0:
                self.full[mol]['freq'] = mol_params[mol][self.full[mol]['t_line']]['f']

            # Start units
            self.full[mol]['units'] = {}
            # Get Intensity units
            int_units = self.Hmol[mol]['BUNIT'].lower()
            k, ut = get_factors_units(int_units, su)

            flux_units = unit2str(ut)
            mags = list(ut.keys())
            if 'temp' in mags:
                n_data = self.Dmol[mol]
            else:

                # Flux units: flux/pixel
                fu = {
                    'f1': {'mag':['flux', 'pixel'],
                           'pow':[1, -1]},
                    'f2': {'mag':['power', 'length', 'freq', 'pixel'],
                           'pow':[1, -2, -1, -1]}
                }

                for f in fu.keys():
                    if fu[f]['mag'].sort() == mags.sort():
                        flux_flag = True
                        for u in mags:
                            idx = fu[f]['mag'].index(u)
                            if ut[u]['pow'] != fu[f]['pow'][idx]:
                                flux_flag = False
                        if flux_flag:
                            # Get factor conversion to mks
                            flux_mks_factor = get_factors_units('jy', 'mks')[0]
                            k *= flux_mks_factor/k
                            # Get Temperature
                            x, y = self.Hmol[mol]['CDELT1'], self.Hmol[mol]['CDELT2']
                            # As it is in MKS or CGS, units are in degrees
                            deg_area = np.abs(x*y)
                            # Get frequency
                            v = self.full[mol]['freq']
                            n_data = brightness_temperature(self.Dmol[mol], v, deg_area)
                            break

            self.full[mol]['units']['flux'] = 'K'
            self.Hmol[mol]['BUNIT'] = 'K'

            # Scale Flux
            self.Dmol[mol] = k*n_data
            
            # Sort velocity cube slices from lower to higher
            if self.Hmol[mol]['CDELT3'] < 0:
                self.Hmol[mol]['CDELT3'] *= -1
                self.Hmol[mol]['CRPIX3'] = self.Hmol[mol]['NAXIS3']-self.Hmol[mol]['CRPIX3']+1

                aux = np.zeros_like(self.Dmol[mol])
                for i in range(self.Hmol[mol]['NAXIS3']):
                    aux[self.Hmol[mol]['NAXIS3']-i-1] = self.Dmol[mol][i]

                self.Dmol[mol] = aux

            # Velocity
            self.full[mol]['vel'] = ((1+np.arange(self.Hmol[mol]['NAXIS3'])-self.Hmol[mol]['CRPIX3'])*self.Hmol[mol]['CDELT3'] + self.Hmol[mol]['CRVAL3'])
            self.full[mol]['units']['vel'] = uts[2]
            # Longitude
            self.full[mol]['lon'] = ks[0]*((1+np.arange(self.Hmol[mol]['NAXIS1'])-self.Hmol[mol]['CRPIX1'])*self.Hmol[mol]['CDELT1'] + self.Hmol[mol]['CRVAL1'])
            self.full[mol]['units']['lon'] = uts[0]
            # Latitude
            self.full[mol]['lat'] = ks[1]*((1+np.arange(self.Hmol[mol]['NAXIS2'])-self.Hmol[mol]['CRPIX2'])*self.Hmol[mol]['CDELT2'] + self.Hmol[mol]['CRVAL2'])
            self.full[mol]['units']['lat'] = uts[1]

            # Centroid
            self.xc, self.yc = self.mc_params['x0'], self.mc_params['y0']

            # Get centre and dimensions, based on 13CO molecule
            self.lat_m = kRe*self.mc_params['dist']*np.tan(np.pi*(self.full[mol]['lat']-self.yc)/180.)
            self.lon_m = kRe*self.mc_params['dist']*np.tan(np.pi*(self.full[mol]['lon']-self.xc)/180.)

            msg(mol+' map loaded', 'ok')
            
        except Exception as e:
            msg(mol+' map not loaded\n'+str(e), 'warn')


    def load_empty_mol(self, mol, dims, vel_lims, lon_lims, lat_lims):
        """
            Loading empty map
            Parameters
            ----------
            mol : string
                Molecule
            dims : tuple
                Dimensions of the empty cube
            vel_dims : tuple
                Velocity dimensions
            lon_dims : tuple
                Longitude limits
            lat_dims : tuple
                Latitude limits
            ----------
        """

        # Create header file
        self.Hmol[mol] = hdr_default

        # Start units
        self.full[mol]['units'] = {}

        # Bright temperature
        self.full[mol]['units']['flux'] = 'K'
        self.Hmol[mol]['BUNIT'] = 'K'

        # Scale Flux
        self.Dmol[mol] = np.zeros(dims)
        
        # Velocity
        self.full[mol]['vel'] = np.linspace(vel_lims[0], vel_lims[1], dims[0])
        self.Hmol[mol]['NAXIS3'] = dims[0]
        self.Hmol[mol]['CDELT3'] = self.full[mol]['vel'][1]-self.full[mol]['vel'][0]
        self.Hmol[mol]['CRVAL3'] = vel_lims[0]
        self.full[mol]['units']['vel'] = hdr_default['CUNIT3']
        # Longitude
        self.full[mol]['lon'] = np.linspace(lon_lims[0], lon_lims[1], dims[2])
        self.Hmol[mol]['NAXIS1'] = dims[2]
        self.Hmol[mol]['CDELT1'] = self.full[mol]['lon'][1]-self.full[mol]['lon'][0]
        self.Hmol[mol]['CRVAL1'] = lon_lims[0]
        self.full[mol]['units']['lon'] = hdr_default['CUNIT1']
        # Latitude
        self.full[mol]['lat'] = np.linspace(lat_lims[0], lat_lims[1], dims[1])
        self.Hmol[mol]['NAXIS2'] = dims[1]
        self.Hmol[mol]['CDELT2'] = self.full[mol]['lat'][1]-self.full[mol]['lat'][0]
        self.Hmol[mol]['CRVAL2'] = lat_lims[0]
        self.full[mol]['units']['lat'] = hdr_default['CUNIT2']

        # Centroid
        self.xc, self.yc = 0, 0

        # Get centre and dimensions, based on 13CO molecule
        kRe = self.mc_params['Re_factor']
        self.lat_m = kRe*self.mc_params['dist']*np.tan(np.pi*(self.full[mol]['lat']-self.xc)/180.)
        self.lon_m = kRe*self.mc_params['dist']*np.tan(np.pi*(self.full[mol]['lon']-self.yc)/180.)

        msg('Map generated for molecule '+mol, 'ok')


    def load_polarization(self, path_pol, **kwargs):
        """
            Loading Polarization data
        """      

        # Key arguments
        # ----------------------------------------------
        # Figure size
        move_cmd = kwargs.pop('move_cmd', None)
        # ----------------------------------------------

        # Load polarization data
        msg('Loading Polarization map', 'info')
        try:
            # Get data and header from polarized data
            self.Dpol, self.Hpol = self.load_fits(path_pol, cube=False)
            # Check if there is enough params
            if len(self.Dpol) >= 3:
                # Extract useful data
                ndata, mags = self.load_stokes_params(self.Dpol, self.Hpol, norm=True)
                self.full['pol_units'] = {}
                for i, k in enumerate(ndata.keys()):
                    if move_cmd:
                        ndata[k] = move_operations(ndata[k], move_cmd)
                    self.full[k] = ndata[k]
                    self.full['pol_units'][k] = mags[i]

                msg('Polarization map loaded', 'ok')

            else:
                msg('Polarization map has not enough data.\nAt least I, Q and U area needed', 'warn')

        except Exception as e:
            msg('Polarization map not loaded\n'+str(e), 'warn')


    def load_stokes_params(self, data, header, norm):
        """
            Loading polaritzation file
            Parameters
            ----------
            data : numpy.array-3d
                Data with all the stokes params
            header : dict
                Data Header
            norm : boolean
                Apply normalization?
            ----------
        """
        # Get dimensions of data
        dims = data.shape
        nparams = 4
        S_PARAMS = ['I','Q','U','V','tau','N']
        # Get Stokes parameters
        params = header['CTYPE3']
        # Get Stokes units
        units = header['CUNIT3']

        for i, p in enumerate(params):
            if p == ' ':
                break

        ndata = {}        
        stk_params = params[i:].split(',')
        if len(stk_params) == 1:
            nparams = dims[0]
            stk_params = S_PARAMS[:nparams]

        # Get values
        for i, s in enumerate(stk_params):
            if 'opt' in s:
                k = 'tau'
            elif 'den' in s:
                k = 'N'
            else:
                k = s.split(' ')[-1].upper()

            ndata[k] = data[i]

        # Get units
        stk_units = units.split(',')
        if len(stk_units) == 1:
            if stk_units[0] == '':
                stk_units = [header['BUNIT']]
            else:
                stk_units = ['K']

        ax_mags = []
        for u in stk_units:
            uc = ''
            for c in u:
                if not c in ['[', ']']:
                    uc += c
            # Get Intensity units
            k, us = get_factors_units(uc.lower(), su)
            ax_mags.append([k, us])

        mags = []
        for m in ax_mags:
            mag = m[1].keys()
            # Flux
            if (('power' in mag) and ('length' in mag) and \
                ('freq' in mag) and ('pixel' in mag)) or \
                (('flux' in mag) and ('pixel' in mag)):
                if not norm:
                    # Get coordinates
                    x, y = header['CDELT1'], header['CDELT2']
                    # As it is in MKS or CGS, units are in degrees
                    deg_area = np.abs(x*y)
                    # Get frequency
                    v = header['RESTFREQ']
                    for p in S_PARAMS[:nparams]:
                        uts = unit2str(m[1])
                        mags.append(uts)
                        # Temperature conversion
                        ndata[p] = brightness_temperature(ndata[p], v, deg_area)
                else:
                    for p in S_PARAMS[:nparams]:
                        mags.append('ND')
            # Temperature
            elif 'temp' in mag:
                if not norm:
                    uts = 'K'
                else:
                    uts = 'ND'
                for p in S_PARAMS[:nparams]:
                    mags.append(uts)
            # Column density
            elif 'length' in mag:
                # If the magntiude is column density
                if m[1]['length']['pow'] == -2:
                    if 'N' in ndata.keys():
                        ndata['N'] = m[0]*ndata['N']
                        uts = unit2str(m[1])
                        #print(uts)
                        mags.append(uts)
            # Adimensional
            elif len(mag) == 0:
                mags.append('ND')

        # Normalization
        if norm:
            I = ndata['I']
            norm = 1./2.**np.int(np.log2(np.nanmean(I[I > 0])))
            if nparams > 4:
                nparams = 4
            for p in S_PARAMS[:nparams]:
                ndata[p] *= norm

        # Transpose
        for d in ndata.keys():
            ndata[d] = ndata[d].T

        return ndata, mags


    def get_pol_params(self, save=False):
        """
            Polarization map: Includes angle and degree of polaritzation
            Parameters
            ----------
            save : boolean
                Save map as fits
            ----------
        """

        if len(self.full['I']) > 0:

            stokes = {}
            for s in ['I', 'Q', 'U', 'V', 'N']:
                if s in self.full.keys():
                    stokes[s] = self.full[s]

            flag_N = False
            if 'N' in stokes.keys():
                flag_N = True
                N = stokes['N']

            msg('Getting Stokes parameters', 'info')

            # Get polarization angle and degree
            pol_degree, pol_angle = get_pol_params_from_stokes(stokes['I'], stokes['Q'], stokes['U'])
            # Get polarization vectors
            pol_vector = pol_degree*np.cos(pol_angle) + 1j*pol_degree*np.sin(pol_angle)

            # Threshold value
            if flag_N:
                N = stokes['N']
                thresh = np.mean(N)
                b = pol_angle[N > 3*thresh]
                stdev_n = np.std(b)

                msg('Mean value: '+str(thresh), 'ok')
                msg('Std dev signal: '+str(stdev_n), 'ok')

            stdev_polangle = np.std(pol_angle)
            msg('Std dev: '+str(stdev_polangle), 'ok')

            self.full['pol_angle'] = pol_angle
            self.full['pol_degree'] = pol_degree
            self.full['pol_vector'] = pol_vector

            # Create header for polaritzation map
            header_pol = self._header_pol(self.Hpol)

            if save:
                # Create the fits map
                os.system('rm -rf pol_angle_'+self.mcID+'.fits')
                hdu = fits.PrimaryHDU(pol_angle, header=header_pol)
                hdu.writeto('pol_angle_'+self.mcID+'.fits')

                os.system('rm -rf pol_degree_'+self.mcID+'.fits')
                hdu = fits.PrimaryHDU(pol_degree, header=header_pol)
                hdu.writeto('pol_degree_'+self.mcID+'.fits')

                if flag_N:
                    os.system('rm -rf col_dens_'+self.mcID+'.fits')
                    hdu = fits.PrimaryHDU(np.log10(N), header=header_pol)
                    hdu.writeto('col_dens_'+self.mcID+'.fits')

                msg('Angle and Degree of polaritzation saved', 'ok')

            msg('Done', 'ok')

        else:
            msg('Polarization data not available', 'fail')


    def binning_pol(self, nbins, rebin, *args, **kwargs):
        """
            Binning polaritzation
            Parameters
            ----------
            nbins : int
                Number of bins
            rebin : boolean
                Renew the binnig?
            ----------
        """

        # Key arguments
        # ----------------------------------------------
        # Figure size
        move_cmd = kwargs.pop('move_cmd', None)
        # ----------------------------------------------

        x, y = np.shape(self.full['I'])
        xlim, ylim = x/nbins, y/nbins

        # Check if the rebinning is possible
        if xlim < 1.0 or ylim < 1.0:
            msg('Binning limit exceeded.\nMax binning is at least 1 pixel per bin',\
            'fail')
            return

        if (xlim > 1 and xlim < 2) or (ylim > 1 and ylim < 2):
            msg('Binning number not valid. It will create bins of zero dimension.',\
            'fail')
            max_bin = np.max((x, y))
            msg('Choose a bin number between 2-'+str(round(max_bin/2))+' or the maximum binning: '+str(max_bin), 'info')
            return

        # Rebinning
        if rebin:
            self.binned = {}
            self.binned['B0'] = {}
            # Initiate molecule lines binning
            for mol in self.mols:
                self._init_binned_data(mol)
            self.bin_grid = np.zeros((nbins**2, 4), dtype=int)
            self.bin_frames = []
            self.nbins = nbins
            # Volume definition
            volume = self.full['vol']
            msg('Starting rebinning', 'ok')
        else:
            if len(self.bin_frames) > 0 and len(self.bin_grid) > 0:
                nbins = self.nbins
            else:
                msg('Grid does not exist. Rebinning needed', 'fail')
                return

        cnt = 0

        if not 'I' in self.full.keys():
            msg('Get polarization angle and degree first', 'fail')
            return

        frames_pts = []
        I = self.full['I']
        Q = self.full['Q']
        U = self.full['U']

        add_vol_pix = 0

        for i in range(nbins):
            for j in range(nbins):
                flag_binning = True
                name = 'B'+str(cnt)
                if not name in self.binned:
                    if rebin:
                        self.binned[name] = {}
                        self.binned[name]['flag'] = False
                        msg('Bin '+str(cnt)+' created', 'ok')
                    else:
                        msg('Bin '+name+' does not exist. Rebinning needed', 'warn')
                        flag_binning = False

                if flag_binning:
                    # Define positions
                    self.binned[name]['pos'] = [i, j]

                    # X-Y limits for each data bin
                    x1 = i*xlim
                    x2 = (i+1)*xlim
                    y1 = j*ylim
                    y2 = (j+1)*ylim
                    pts = ((x1,x2),(y1,y2))
                    # Get max and min limits
                    xmin = int(np.ceil(x1))
                    xmax = int(x2)
                    ymin = int(np.ceil(y1))
                    ymax = int(y2)
                    lims = ((xmin,xmax),(ymin,ymax))
                    # Extract fractional data
                    if (x%nbins > 0) or (y%nbins > 0):
                        if rebin:
                            # Get frames
                            frames_pts = self._fractional_binning(pts, lims, (x,y))
                            # Get volume
                            add_vol_pix = self._int_frac_pixels(frames_pts, volume)
                        else:
                            frames_pts = self.bin_frames[cnt]
                        # Get pol data and frame
                        add_I_pix, frame_I_pix = self._int_frac_pixels(frames_pts, I, get_frame=True)
                        add_Q_pix, frame_Q_pix = self._int_frac_pixels(frames_pts, Q, get_frame=True)
                        add_U_pix, frame_U_pix = self._int_frac_pixels(frames_pts, U, get_frame=True)
                    else:
                        # No extra volume is null
                        add_vol_pix = 0
                        # Frames are empty
                        frames_pts = []
                        frame_I_pix = []
                        frame_Q_pix = []
                        frame_U_pix = []

                    # Asign volume
                    if rebin:
                        self.binned[name]['vol'] = np.sum(volume[xmin:xmax, ymin:ymax])+add_vol_pix
                        # Asign Grid
                        self.bin_grid[cnt] = np.array([xmin, xmax, ymin, ymax])
                        # Asign frames
                        self.bin_frames.append(frames_pts)

                    # Get Stokes sections
                    bin_I = I[xmin:xmax, ymin:ymax]
                    bin_Q = Q[xmin:xmax, ymin:ymax]
                    bin_U = U[xmin:xmax, ymin:ymax]

                    pol_degree, pol_angle = get_pol_params_from_stokes(bin_I, bin_Q, bin_U)
                    pol_vector = pol_degree*np.cos(pol_angle) + 1j*pol_degree*np.sin(pol_angle)
                    if len(frame_I_pix) > 0:
                        frame_pol_degree, frame_pol_angle = get_pol_params_from_stokes(frame_I_pix, frame_Q_pix, frame_U_pix)
                        frame_pol_vector = frame_pol_degree*np.cos(frame_pol_angle) + 1j*frame_pol_degree*np.sin(frame_pol_angle)
                    else:
                        frame_pol_degree, frame_pol_angle = [], []
                        frame_pol_vector = []

                    # Polarization asignation
                    self.binned[name]['pol_angle'] = {'Area': pol_angle, 'Frame': frame_pol_angle}
                    self.binned[name]['pol_degree'] = {'Area': pol_degree, 'Frame': frame_pol_degree}
                    self.binned[name]['pol_vector'] = {'Area': pol_vector, 'Frame': frame_pol_vector}                 

                    # Basic stadistic of polarizated data
                    pol_angle_data = np.concatenate((pol_angle.flatten(), frame_pol_angle))
                    pol_degree_data = np.concatenate((pol_degree.flatten(), frame_pol_degree))
                    # Polarization Angle
                    self.binned[name]['pol_angle_mean'] = np.mean(pol_angle_data) 
                    angle_std = np.std(pol_angle_data)
                    if angle_std < 1e-5:
                        angle_std = np.nan
                    self.binned[name]['pol_angle_std'] = angle_std
                    # Polarization Degree
                    self.binned[name]['pol_degree_mean'] = np.mean(pol_degree_data)
                    self.binned[name]['pol_degree_std'] = np.std(pol_degree_data)
                    # Polarization vector averaged per bin
                    self.binned[name]['pol_vector_avg'] = np.mean(pol_degree_data)*np.cos(np.mean(pol_angle_data))+1j*np.mean(pol_degree_data)*np.sin(np.mean(pol_angle_data))

                    cnt += 1

        msg('Binning complete', 'info')


    # def M0(self, mol, save=True):
    #     """
    #         Zero moment map
    #         Parameters
    #         ----------
    #         mol : string
    #             Map Molecule
    #         save : boolean
    #             Save map as fits image
    #         ----------
    #     """
    #     if mol in self.mols:
    #         if self.Dmol[mol] == []:
    #             msg(mol + ' Molecule data not available', 'fail')
    #             return
    #         # Parameters
    #         cdelt = self.Hmol[mol]['CDELT*']
    #         crpix = self.Hmol[mol]['CRPIX*']
    #         crval = self.Hmol[mol]['CRVAL*']
    #         # Dimensions
    #         _, y, x = np.shape(self.Dmol[mol])
    #         data = self.Dmol[mol]
    #     else:
    #         msg('Molecule not valid', 'fail')
    #         return

    #     # Get Zero Momentum
    #     m0 = cdelt[2]*np.nansum(data, axis=0)

    #     # Create fits header
    #     header = self._header_vel(self.Hmol[mol])

    #     # Save the map
    #     if save:
    #         # Save data as fits file
    #         hdu = fits.PrimaryHDU(m0, header=header)
    #         os.system("rm -rf map_m0_" + mol + "_" + self.mcID + ".fits")
    #         hdu.writeto('map_m0_' + mol + "_" + self.mcID + ".fits")

    #     M0 = {
    #         'data':m0.T,
    #         'header':header
    #     }

    #     return M0


    def get_n_moment(self, mol, n=0, save=False):
        """
            One Moment map
            Parameters
            ----------
            mol : string
                Map Molecule
            n : int
                Moment order
            save : boolean
                Save map as fits image
            ----------
        """
        if mol in self.mols:
            if self.Dmol[mol] == []:
                msg(mol + ' Molecule data not available', 'fail')
                return
            # Parameters
            cdelt = self.Hmol[mol]['CDELT*']
            crpix = self.Hmol[mol]['CRPIX*']
            crval = self.Hmol[mol]['CRVAL*']
            # Dimensions
            _, y, x = np.shape(self.Dmol[mol])
            data = self.Dmol[mol]
            # Get velocity
            vel = ((1+np.arange(self.Hmol[mol]['NAXIS3'])-crpix[2])*cdelt[2] + crval[2])
        else:
            msg('Molecule not valid', 'fail')
            return

        # Get 0th-Moment
        m0 = cdelt[2]*np.nansum(data, axis=0)

        if n <= 0:
            mn = m0
        else:
            # Get 1st-Moment
            m1 = np.zeros((y, x))
            for i in range(x):
                for j in range(y):
                    m1[j][i] = cdelt[2]*np.nansum(data[:, j, i]*(vel))
            m1 = m1/m0

            if n > 1:
                # Get N-Momentu
                mn = np.zeros((y, x))
                for i in range(x):
                    for j in range(y):
                        mn[j][i] = cdelt[2]*np.nansum(data[:, j, i]*(vel-m1[j][i])**n)
                mn = mn/m0
            else:
                mn = m1

        # Create fits header
        header = self._header_vel(self.Hmol[mol])

        # Save the map
        if save:
            # Save data as fits file
            hdu = fits.PrimaryHDU(mn, header=header)
            os.system("rm -rf map_m"+str(n)+"_" + mol + "_" + self.mcID + ".fits")
            hdu.writeto('map_m'+str(n)+'_' + mol + "_" + self.mcID + ".fits")

        # Get units
        if n <= 0:
            units = self.full[mol]['units']['flux']+' '+self.full[mol]['units']['vel']
        else:
            units = self.full[mol]['units']['vel']
            if n > 1:
                units = '('+units+')^'+str(n)

        Mn = {
            'units':units,
            'data':mn.T,
            'header':header
        }

        return Mn


    def get_M0_mols(self, save=False):
        """
            Get all the Momentum Zero maps for all the molecules
        """
        M0s = {}
        for mol in self.mols:
            M0 = self.get_n_moment(mol, n=0, save=save)
            M0s[mol] = {}
            M0s[mol]['data'] = M0['data']
            M0s[mol]['header'] = M0['header']

        return M0s


    def extract_param_from_binned(self, param):
        """
            Extract molecule from binned data. 
            Parameters
            ----------
            mol : string
                Molecule to extract from binned data
            ----------
        """
        params = list(self.binned['B0'].keys())
        if (not param in params) or param == 'pos':
            msg('Parameter not valid', 'fail')
            return

        mc_mol = {}
        for b in self.binned.keys():
            if isinstance(self.binned[b][param], dict):
                mc_mol[b] = self.binned[b][param]
            else:
                mc_mol[b] = {} 
                mc_mol[b][param] = self.binned[b][param]

            mc_mol[b]['pos'] = self.binned[b]['pos']
            mc_mol[b]['vol'] = self.binned[b]['vol']

        return mc_mol


    def map_vel_mol(self, mol, method='avg'):
        """
            Integrated velocity map of the whole region by molecule
            Parameters
            ----------
            mol : string
                Map Molecule
            ----------
        """
        if mol in self.mols:
            if self.Dmol[mol] == []:
                msg(mol + ' Molecule data not available', 'fail')
                return
            # Parameters
            cdelt = self.Hmol[mol]['CDELT*']
            # Dimensions
            v = np.shape(self.Dmol[mol])[-3]
            data = self.Dmol[mol]
        else:
            msg('Molecule not valid', 'fail')
            return

        # Get map velocity
        map_vel = self.map_vel(v, data, cdelt, method)

        return map_vel


    def map_vel(self, v, data, cdelt, method='avg', *args, **kwargs):
        """
            Integrated velocity map of the data region
            Parameters
            ----------
            v : array
                Velocity array
            data : array
            cdelt : float
            ----------
        """
        
        # Key arguments
        # ----------------------------------------------
        # Frame definition
        frame = kwargs.pop('frame', np.zeros(v))
        # Number of frame bits
        n_frames = kwargs.pop('n_frames', 0)
        # ----------------------------------------------

        # Initiate velocity array
        map_vel = np.zeros(v)
        # Get the velocity profile
        for i in range(v):
            # Remove any nan in data
            data_f = data[i,:,:].flatten()
            data_f = data_f[~np.isnan(data_f)]
            # Get the sum
            map_sum = np.sum([np.sum(data_f),frame[i]])

            if method == 'int':
                map_vel[i] = np.abs(cdelt[1]*cdelt[0])*map_sum
            elif method == 'avg':
                map_vel[i] = map_sum/(n_frames+len(data_f))

        return map_vel


    def get_map_vels(self, method='avg'):
        """
            Get integrated map velocities for all the molecules
        """
        msg('Creating integrated velocity maps...', 'info')
        for mol in self.mols:
            # Define the molecular lines
            self.full[mol]['line'] = self.map_vel_mol(mol, method)
            # Define the baselines
            self.baseline[mol] = np.zeros_like(self.full[mol]['line'])

        msg('Done', 'ok')


    def get_noise(self, mol, channels=[0,-1], margins=None, *args, **kwargs):
        """
            Get the background noise
        """

        # Key arguments
        # ----------------------------------------------
        # Verbose
        verbose = kwargs.pop('verbose', False)
        # ----------------------------------------------

        # Number of channels
        nchns, xdim, ydim = self.Dmol[mol].shape
        # Noise vector
        noise = np.nan*np.zeros(nchns)
        # Get cdelt
        cdelt = self.Hmol[mol]['CDELT*']
        # Get noise data from vector
        if isinstance(channels, list) or isinstance(channels, np.ndarray):
            if len(channels) > 0:
                for chn in channels:
                    if chn < nchns:
                        noise[chn] = np.nansum(self.Dmol[mol][chn,:,:])
                    else:
                        msg('Channel '+str(chn)+' does not exist', 'warn')

        # Get noise from margins
        if isinstance(margins, tuple):
            if len(margins) == 2:
                # Get limits in pixels
                xmin = int(margins[0][0]*xdim)
                xmax = int(margins[0][1]*xdim)

                ymin = int(margins[1][0]*ydim)
                ymax = int(margins[1][1]*ydim)

                # Get data for each margin
                for chn in range(nchns):
                    if not chn in channels:
                        lmargin = self.Dmol[mol][chn][:ymin].flatten()
                        rmargin = self.Dmol[mol][chn][ymax:].flatten()
                        
                        tmargin = self.Dmol[mol][chn][ymin:ymax,:xmin].flatten()
                        bmargin = self.Dmol[mol][chn][ymin:ymax,xmax:].flatten()

                        noise_frame = np.concatenate((lmargin, rmargin, tmargin, bmargin))
                        k = (xdim*ydim)/len(noise_frame)
                        noise[chn] = np.nansum(noise_frame)*k

        # Clean from nan arrays
        noise = noise/(xdim*ydim)
        self.noise[mol] = noise

        # Noise statistics
        self.noise_stat[mol] = {
                'mean':np.nanmean(noise),
                'std':np.nanstd(noise),
                'median':np.nanmedian(noise),
                'rms': np.sqrt(np.nanmean(noise**2))
        }

        if verbose:
            print('========== MAP STATISTICS ==========')
            print('RMS noise map: ', self.noise_stat[mol]['rms'])
            print('Mean noise map: ', self.noise_stat[mol]['mean'])
            print('Standard deviation map: ', self.noise_stat[mol]['std'])
            print('Median map: ', self.noise_stat[mol]['median'])


    def get_param_stat(self, param, *args, **kwargs):
        """
            Get the statistic of a parameter
            Parameters
            ----------
            param : str
                Parameter to get the statistic
            ----------
        """

        # Key arguments
        # ----------------------------------------------
        # Verbose
        verbose = kwargs.pop('verbose', True)
        # ----------------------------------------------

        params = []
        for nbin in self.binned.keys():
            if not self.binned[nbin]['flag']:
                params.append( self.binned[nbin][param] )

        params = np.array(params)

        # Get mean
        mean_param = np.nanmean(params)
        # Get median
        med_param = np.nanmedian(params)
        # Get std
        std_param = np.nanstd(params)

        if verbose:
            print('========== PARAM STATISTICS ==========')
            print('Mean:               ', mean_param)
            print('Standard deviation: ', med_param)
            print('Median:             ', std_param)

        return mean_param, med_param, std_param, params


    def substract_baseline(self, mol, inter=False, **kwargs):
        """
            Substract the baseline
            Parameters
            ----------
            mol : string
                Molecule
            inter : boolean
                Activate the interactive mode
            ----------
        """
        if mol in self.mols:
            if self.full[mol]['line'] == []:
                msg(mol + ' Molecule data not available', 'fail')
                return

            if inter:
                # Add canvas to the app
                ioff()
                fig, ax = subplots()
                close(fig)

                # Instance app
                app = QApplication.instance()

                self.SubBaseLine = BaselineSpecWindow()
                units = [self.full[mol]['units']['vel'], self.full[mol]['units']['flux']]
                save = kwargs.pop('save', False)
                name_bin = self.mcID + '-' + mol
                self.SubBaseLine.load_init_params(fig, ax, self.full[mol]['vel'], self.full[mol]['line'], mol, units, name_bin, save=save)
                # Signal connection to extract baseline data
                self.SubBaseLine.signal_baseline.connect(self._get_baseline_from_ui)
                self.SubBaseLine.show()

                app.exec_()

            else:
                # Get the method from the kwargs
                method = kwargs.pop('method', 'linear')
                if method == 'linear':
                    degree = 1
                elif method == 'poly':
                    degree = kwargs.pop('ndeg', 2)
                elif method == 'bls':
                    l = kwargs.pop('l', 105)
                    p = kwargs.pop('p', 0.05)
                    n = kwargs.pop('n', 10)
                else:
                    msg('Define a valid method to fit the baseline:\n1. Linear\n2. Polynomial\n3. BLS method', 'fail')
                    return

                # Length size
                Nsize = len(self.full[mol]['line'])
                # Get preliminar baseline
                edges = int(Nsize/20)
                prem_baseline = poly_baseline(np.concatenate((self.full[mol]['vel'][:edges], self.full[mol]['vel'][Nsize-edges:])),
                               np.concatenate((self.full[mol]['line'][:edges], self.full[mol]['line'][Nsize-edges:])), 1, 
                               self.full[mol]['vel'])
                # Get find peaks params from kwargs
                dist = kwargs.pop('dist', 5.)
                height_div = kwargs.pop('height_div', 12.)
                # Get peaks and widths
                peaks = find_profile_peaks(self.full[mol]['line']-prem_baseline, dist=dist, height_div=height_div)
                # Guess line parameters to remove
                guess = guess_line(self.full[mol]['vel'], self.full[mol]['line']-prem_baseline, peaks)
                baseline = self.full[mol]['line'].copy()
                dv = np.mean(np.diff(self.full[mol]['vel']))
                # Define the line width constant
                klw = kwargs.pop('klw', 6)

                rem_items = np.array([])
                for i, peak in enumerate(peaks):
                    line_width = klw*int(guess[3*i+2]/dv)
                    # Inferior limit
                    lim_inf = peak-line_width/2
                    if lim_inf < 0:
                        lim_inf = 0
                    # Superior limit
                    lim_sup = peak+line_width/2
                    if lim_sup > Nsize:
                        lim_sup = Nsize

                    rem_aux_items = np.arange(lim_inf, lim_sup)
                    rem_items = np.concatenate((rem_items, rem_aux_items))

                baseline = np.delete(baseline, rem_items.astype(int))
                freq_baseline = np.delete(self.full[mol]['vel'], rem_items.astype(int))

                try: 
                    if method == 'bls':
                        fit_baseline = baseline_als_optimized(baseline, l, p, niter=n)
                    elif (method == 'linear') or (method == 'poly'):
                        fit_baseline = poly_baseline(freq_baseline, baseline, degree, self.full[mol]['vel'])
                    
                    self.baseline[mol] = fit_baseline
                except Exception as e:
                    msg('Baseline couldnt be adjusted.'+str(e), 'fail')                    

        else:
            msg('Molecule not available', 'fail')
            return


    def _get_baseline_from_ui(self, kind):
        """
            Get data with baseline substracted from the UI
            Parameters
            ----------
            kind : string
                Type of spectra [molecule]
            ----------
        """
        fields = kind.split('/')
        if len(fields) > 1:
            name = fields[0]
            mol = fields[1]
            self.binned[name][mol]['bl'] = self.SubBaseLine.data_corrected
        else:
            mol = fields[0]
            self.baseline[mol] = self.SubBaseLine.data_corrected


    def substract_baseline_binning(self, mol, inter=False, **kwargs):
        """
            Substract the baseline
            Parameters
            ----------
            mol : string
                Molecule
            inter : boolean
                Activate the interactive mode
            ----------
        """
        if mol in self.mols:
            names = self.binned.keys()
            for name in names:
                if len(self.binned[name][mol]['line']) == 0:
                    msg(mol + ' Molecule data not available', 'fail')
                    return

                if inter:
                    # Add canvas to the app
                    ioff()
                    fig, ax = subplots()
                    close(fig)

                    # Instance app
                    app = QApplication.instance()

                    self.SubBaseLine = BaselineSpecWindow()
                    units = [self.full[mol]['units']['vel'], self.full[mol]['units']['flux']]
                    save = kwargs.pop('save', False)
                    name_bin = self.mcID + '-' + mol + '-' + name
                    self.SubBaseLine.load_init_params(fig, ax, self.binned[name][mol]['vel'], self.binned[name][mol]['line'], name+'/'+mol, units, name_bin, save=save)
                    # Signal connection to extract baseline data
                    self.SubBaseLine.signal_baseline.connect(self._get_baseline_from_ui)
                    self.SubBaseLine.show()

                    app.exec_()

                else:
                    # Get the method from the kwargs
                    method = kwargs.pop('method', 'linear')
                    if method == 'linear':
                        degree = 1
                    elif method == 'poly':
                        degree = kwargs.pop('ndeg', 2)
                    elif method == 'bls':
                        l = kwargs.pop('l', 105)
                        p = kwargs.pop('p', 0.05)
                        n = kwargs.pop('n', 10)
                    else:
                        msg('Define a valid method to fit the baseline:\n1. Linear\n2. Polynomial\n3. BLS method', 'fail')
                        return

                    # Length size
                    Nsize = len(self.binned[name][mol]['line'])
                    # Get preliminar baseline
                    edges = int(Nsize/20)
                    prem_baseline = poly_baseline(np.concatenate((self.binned[name][mol]['vel'][:edges], self.binned[name][mol]['vel'][Nsize-edges:])),
                                   np.concatenate((self.binned[name][mol]['line'][:edges], self.binned[name][mol]['line'][Nsize-edges:])), 1, 
                                   self.binned[name][mol]['vel'])
                    # Get find peaks params from kwargs
                    dist = kwargs.pop('dist', 5.)
                    height_div = kwargs.pop('height_div', 12.)
                    # Get peaks and widths
                    peaks = find_profile_peaks(self.binned[name][mol]['line']-prem_baseline, dist=dist, height_div=height_div)
                    # Guess line parameters to remove
                    guess = guess_line(self.binned[name][mol]['vel'], self.binned[name][mol]['line']-prem_baseline, peaks)
                    baseline = self.binned[name][mol]['line'].copy()
                    dv = np.mean(np.diff(self.binned[name][mol]['vel']))
                    # Define the line width constant
                    klw = kwargs.pop('klw', 6)

                    rem_items = np.array([])
                    for i, peak in enumerate(peaks):
                        line_width = klw*int(guess[3*i+2]/dv)
                        # Inferior limit
                        lim_inf = peak-line_width/2
                        if lim_inf < 0:
                            lim_inf = 0
                        # Superior limit
                        lim_sup = peak+line_width/2
                        if lim_sup > Nsize:
                            lim_sup = Nsize

                        rem_aux_items = np.arange(lim_inf, lim_sup)
                        rem_items = np.concatenate((rem_items, rem_aux_items))

                    baseline = np.delete(baseline, rem_items.astype(int))
                    freq_baseline = np.delete(self.binned[name][mol]['vel'], rem_items.astype(int))

                    try: 
                        if method == 'bls':
                            fit_baseline = baseline_als_optimized(baseline, l, p, niter=n)
                        elif (method == 'linear') or (method == 'poly'):
                            fit_baseline = poly_baseline(freq_baseline, baseline, degree, self.binned[name][mol]['vel'])
                        
                        self.binned[name][mol]['bl'] = fit_baseline
                    except Exception as e:
                        msg('Baseline couldnt be adjusted.'+str(e), 'fail')                    

        else:
            msg('Molecule not available', 'fail')
            return


    def line_fit(self, mol, inter=False, **kwargs):
        """
            Fit the line profile of the whole map
            Parameters
            ----------
            mol : string
                Molecule
            inter : boolean
                Activate the interactive mode
            ----------
        """
        # Force n lines according with the number of lines of the thin molecule
        # This forced flag has higher hierachy than the forced lines
        forced_thin = kwargs.pop('forced_thin', False)
        # Force n lines
        forced_lines = kwargs.pop('forced_lines', None)
        if not forced_lines is None:
            forced_lines = int(forced_lines)
            if forced_lines < 0:
                forced_lines = 0
    
        # Flag lines
        flag_maps = kwargs.pop('flag_maps', True)
        if flag_maps:
            if len(self.noise[mol]) < 2:
                msg('Not enough noise data selected', 'warn')
                flag_maps = False

        # Sigma threshold
        sigma_thresh = kwargs.pop('sigma_thresh', 1)

        # Sigma noise
        sigma_type = kwargs.pop('sigma_type', 'rms')
        if not sigma_type in ['rms', 'std']:
            msg('Sigma type is not valid. RMS assigned by default', 'warm')
            sigma_type = 'rms'

        # verbose
        verbose = kwargs.pop('verbose', False)

        if mol in self.mols:
            if len(self.full[mol]['line']) == 0:
                msg(mol + ' Molecule spectra data not available', 'fail')
                return

            # Spectra without baseline
            spectra = self.full[mol]['line'] - self.baseline[mol]

            # Check if the map is not noise
            proceed = False
            if flag_maps:
                if sigma_type == 'std':
                    mean_noise = self.noise_stat[mol]['mean']
                    std_noise = self.noise_stat[mol]['std']
                    line_thresh = mean_noise+sigma_thresh*std_noise
                else:
                    rms_noise = self.noise_stat[mol]['rms']
                    line_thresh = sigma_thresh*rms_noise
                
                above_thresh_points = np.where(np.abs(spectra) > line_thresh)[0]

                if len(above_thresh_points) > 0:
                    # Detect spikes
                    diff_above_thresh = np.diff(above_thresh_points)
                    next_above_points = np.where(diff_above_thresh == 1)[0]
                    if len(next_above_points) > 0:
                        proceed = True
            else:
                proceed = True

            if proceed:
                # Forced the same number of lines than the thin molecule
                if forced_thin:
                    #forced_lines = 0
                    if 'A' in self.full[self.thin].keys():
                        forced_lines = len(self.full[self.thin]['A'])

                # If the number of lines are forced
                if forced_lines is None:
                    # Find peaks and get an guess estimation
                    peaks = find_profile_peaks(spectra, **kwargs)
                else:
                    # If forced lines is actived
                    if forced_lines == 0:
                        peaks = []
                    elif forced_lines == 1:
                        peaks = [np.argmax(spectra)]
                    else:
                        item_spectra = np.arange(len(spectra))
                        ndiv = int(np.ceil(len(spectra)/(2+forced_lines)))
                        peaks_raw = item_spectra[::ndiv]
                        peaks = peaks_raw[1:-1]
            else:
                peaks = []
           
            if len(peaks) > 0:
                # Verbose for peaks
                if verbose:
                    print("========== Peaks founded ==========")
                    print("Number of peaks founded: ", len(peaks))
                    for pk in peaks:
                        print("-----------------------------------")
                        print("Peak position: ", self.full[mol]['vel'][pk])
                        print("Peak amplitude: ", spectra[pk])

                popt, pcov = solve_line(self.full[mol]['vel'], spectra, peaks, ['G']*len(peaks))
                A, mu, sigma, fwhm = get_params_from_popt(popt, len(peaks))
            else:
                pcov = []
                A, mu, sigma, fwhm = [], [], [], []
                msg('No lines found', 'warn')

            # Sort parameters through the mu position
            mu_sort = np.argsort(mu)
            A = [A[i] for i in mu_sort]
            mu = [mu[i] for i in mu_sort]
            sigma = [sigma[i] for i in mu_sort]
            fwhm = [fwhm[i] for i in mu_sort]

            # Save line parameters
            self.full[mol]['A'] = A
            self.full[mol]['mu'] = mu
            self.full[mol]['sigma'] = sigma
            self.full[mol]['fwhm'] = fwhm
            self.full[mol]['err'] = pcov

            # If active mode is active
            if inter:
                # Add canvas to the app
                ioff()
                fig, ax = subplots()
                close(fig)

                # Instance app
                app = QApplication.instance()

                # Create interactive window
                self.FitSpectra = FitLinesWindow()
                # Define initial params
                units = [self.full[mol]['units']['vel'], self.full[mol]['units']['flux']]
                save = kwargs.pop('save', False)
                name_bin = self.mcID + '-' + mol
                # Load initial params
                self.FitSpectra.load_init_params(fig, ax, self.full[mol]['vel'], spectra, mol, units, name_bin, save=save)
                # Show the data pre-adjusted
                if len(A) > 0:
                    self.FitSpectra.load_fit_params(popt, pcov, peaks, ['G']*len(peaks))
                # Signal connection to extract baseline data
                self.FitSpectra.signal_fitting.connect(self._get_fitting_from_ui)
                self.FitSpectra.show()

                app.exec_()

        else:
            msg('Molecule not available', 'fail')
            return


    def _get_fitting_from_ui(self, kind):
        """
            Get lines fitted to the spectra from the UI
            Parameters
            ----------
            kind : string
                Type of spectra [molecule]
            ----------
        """
        # Assign parameters
        A = []
        mu = []
        sigma = []
        fwhm = []
        fit_data = self.FitSpectra.fit_data_dict
        # Extract parameters
        for ln in fit_data.keys():
            A.append(fit_data[ln][0])
            mu.append(fit_data[ln][1])
            sigma.append(fit_data[ln][2])
            fwhm.append(2*np.sqrt(2*np.log(2))*fit_data[ln][2])
        
        # Error covariance matrix
        pcov = self.FitSpectra.pcov

        # Sort parameters through the mu position
        mu_sort = np.argsort(mu)
        A = [A[i] for i in mu_sort]
        mu = [mu[i] for i in mu_sort]
        sigma = [sigma[i] for i in mu_sort]
        fwhm = [fwhm[i] for i in mu_sort]

        fields = kind.split('/')
        if len(fields) > 1:
            # Get name and molecule
            name = fields[0]
            mol = fields[1]
            # Save line params
            self.binned[name][mol]['A'] = A
            self.binned[name][mol]['mu'] = mu
            self.binned[name][mol]['sigma'] = sigma
            self.binned[name][mol]['fwhm'] = fwhm
            self.binned[name][mol]['err'] = pcov

        else:
            mol = fields[0]
            # Save line params
            self.full[mol]['A'] = A
            self.full[mol]['mu'] = mu
            self.full[mol]['sigma'] = sigma
            self.full[mol]['fwhm'] = fwhm
            self.full[mol]['err'] = pcov


    def bin2array(self, mol):
        """
            Convert from bin format to dict array format
            Parameters
            ----------
            mol : string
                molecule or parameter
            ----------
        """
        aux_param = {}

        for name in self.binned.keys():
            if not mol in aux_param.keys():
                if type(self.binned[name][mol]) is dict:
                    aux_param[mol] = {}
                else:
                    aux_param[mol] = []

            if type(self.binned[name][mol]) is dict:
                for param in self.binned[name][mol].keys():
                    if not param in ['vel', 'line', 'freq', 'bl']:
                        if not param in aux_param[mol].keys():
                            aux_param[mol][param] = []

                        if type(self.binned[name][mol][param]) is list:
                            value_param = self.binned[name][mol][param]
                        else:
                            value_param = [self.binned[name][mol][param]]
                        aux_param[mol][param] += value_param
            else:
                aux_param[mol].append(self.binned[name][mol])

        return aux_param


    def fixing_params(self, mol, vals=1, method='def', **kwargs):
        """
            Fix the parameters of the binned lines
            Parameters
            ----------
            mol : string
                Molecule
            vals : float or array
                Values to asign to the defined params
            method : string
                Method to asign the fixed value:
                'def': Values defined by the user
                'mean': Gets the mean of a param and asigned as fixed value
            params [kwarg] : string or array
                Params to fixed their values. By default 'sigma' param is fixed
            ----------
        """
        # Key arguments
        # ----------------------------------------------
        # Parameters
        params = kwargs.pop('params', 'sigma')
        # ----------------------------------------------

        if type(params) is str:
            params = [params]
            vals = [vals]

        # Check if the field is valid
        if not mol in sos.all_params:
            msg('Parameter not valid', 'fail')
            return

        # Check dimensions
        if not len(params) == len(vals):
            msg('Parameters and values must have the same length', 'fail')
            return

        # Check if the method is valid
        if not method in ['def', 'mean']:
            msg('Method not valid. Use "def" or "mean"', 'fail')
            return

        flag_width = True
        # Fix the lines with the mean values
        for i,p in enumerate(params):
            # Check if the parameter is valid
            if p in ['T', 'A', 'mu', 'sigma', 'fwhm']:
                if method == 'mean':
                    data = self.bin2array(mol)
                    if type(data[mol]) is dict:
                        val = np.mean(data[mol][p])
                    else:
                        val = np.mean(data[mol])
                elif method == 'def':
                    val = vals[i]

                # Fix a value for all the bins
                for name in self.binned.keys():

                    if type(self.binned[name][mol]) is dict:
                        # Check if it is sigma or fwhm
                        if p == 'sigma' and flag_width:
                            self.binned[name][mol]['fwhm'] = [np.sqrt(8*np.log(2))*val]
                        elif p == 'fwhm' and flag_width:
                            self.binned[name][mol]['sigma'] = [val/np.sqrt(8*np.log(2))]

                        # Check for A or mu parameter
                        if p in ['A', 'mu', 'sigma', 'fwhm']:
                            self.binned[name][mol][p] = [val]
                        else:
                            self.binned[name][mol][p] = val

                    else:
                        if p in ['A', 'mu', 'sigma', 'fwhm']:
                            self.binned[name][mol] = [val]
                        else:
                            self.binned[name][mol] = val

                # Check in case of have both parameters sigma/fwhm, choose only sigma
                if params in ['sigma', 'fwhm']:
                    flag_width = False
            else:
                msg('Sub-parameter not valid', 'warn')

        msg('Done', 'ok')


    def get_gral_params(self, **kwargs):
        """
            Get parameters of the whole molecular cloud
        """
        # Key arguments
        # ----------------------------------------------
        # Sigma of Gaussian filter to get the lines temperatures from maps
        sigma_filter = kwargs.pop('sigma_filter', 0.5)
        # Define virial molecule
        mol_vir = kwargs.pop('mol_vir', 'thin')
        # Method to get radiation temperature 
        method = kwargs.pop('method', 'max')
        # Get general parameters by full map 
        full_map = kwargs.pop('full_map', False)
        # ----------------------------------------------

        # Check if fwhm is calculated
        if not 'fwhm' in self.full[self.thin].keys():
            msg('FWHM is not defined', 'fail')
            return

        # Get conversion factor to km/s
        f, _ = get_factors_units('km s^-1', su)

        # Main physical parameters
        dist = self.mc_params['dist']
        if full_map:
            lon_len = np.abs(self.lon_m[-1] - self.lon_m[0])
            lat_len = np.abs(self.lat_m[-1] - self.lat_m[0])
            A_area = lon_len*lat_len
            theta = 2*np.arctan(np.sqrt(A_area/(np.pi*dist**2)))
            theta = theta*180/np.pi
        else:
            theta = self.mc_params['ang_dia']
        
        kRe = self.mc_params['Re_factor']
        # Put params together
        params = [dist, theta, kRe]

        msg('Calculating parameters...', 'info')

        if len(self.full[self.thin]['line']) > 0 and len(self.full[self.thick]['line']) > 0:

            # Get some factors
            k = self.mc_params['temp_factor']           # Instrumental factor
            k_sigma = self.mc_params['width_factor']    # Sigma factor

            # Get FWHM 
            if mol_vir == 'thin':
                fwhms_vir = self.full[self.thin]['fwhm']
            elif mol_vir == 'thick':
                fwhms_vir = self.full[self.thick]['fwhm']

            # Get FWHM for thin molecule
            fwhms = self.full[self.thin]['fwhm']
            # Number of lines
            nline = len(fwhms)

            # Check if there are lines
            if nline == 0:
                msg('No lines found', 'fail')
                return 

            # LTE mass and N13CO flags
            flag_match_lines = False
            if (nline == len(self.full[self.thick]['A'])) and (nline == len(self.full[self.thin]['A'])):
                flag_match_lines = True
            else:
                msg('Number of '+self.thin+' and '+self.thick+' lines are not the same.\nMass by LTE method and densities columns wont be calculated', 'warn')
                return

            # Get frequency molecules
            v_thick = self.full[self.thick]['freq']
            v_thin = self.full[self.thin]['freq']

            # Get Einstein A coeficient
            line_thin = self.full[self.thin]['t_line']
            for m in mol_params[self.thin].keys():
                if line_thin in m:
                    A10_thin = mol_params[self.thin][m]['A10']
                    X_thin = mol_params[self.thin]['X']
                    J = int(line_thin[0])

            # Initialise parameters
            MVIR_a, MLTE_a, MXF_a, N_13CO_a, N_H2_a, Re = 0, 0, 0, 0, 0, 0
            # Excitation temperature
            Texs = np.zeros(nline)

            T13s = np.zeros(nline)
            T12s = np.zeros(nline)

            # Get parameters by line found
            for idx, fwhm in enumerate(fwhms):

                # Line parameters thin molecule
                # ------------------------------------------------------
                A = self.full[self.thin]['A'][idx]
                mu = self.full[self.thin]['mu'][idx]
                sigma = self.full[self.thin]['sigma'][idx]

                # Line parameters thick molecule
                # ------------------------------------------------------
                A_thick = self.full[self.thick]['A'][idx]
                mu_thick = self.full[self.thick]['mu'][idx]
                sigma_thick = self.full[self.thick]['sigma'][idx]

                # Get line temperatures
                # Use 12CO(optically thick molecule)
                # ------------------------------------------------------

                # Code used to get T12 and T13, now it is obsolete
                # ------------------------------------------------------
                # Thick molecule
                # line_left = np.where(mu>=self.full[self.thick]['vel'])[0]
                # if len(line_left) > 0:
                #     idx_left = line_left[-1]
                # else:
                #     idx_left = 0

                # idx_right = idx_left + 1
                # if idx_right >= len(self.full[self.thick]['vel']):
                #     idx_right = idx_left

                # if (mu - self.full[self.thick]['vel'][idx_left])<=(self.full[self.thick]['vel'][idx_right] - mu):
                #     idx_slice_thick = idx_left
                # else:
                #     idx_slice_thick = idx_right

                # # Thin molecule
                # line_left = np.where(mu>=self.full[self.thin]['vel'])[0]
                # if len(line_left) > 0:
                #     idx_left = line_left[-1]
                # else:
                #     idx_left = 0

                # idx_right = idx_left + 1
                # if idx_right >= len(self.full[self.thin]['vel']):
                #     idx_right = idx_left

                # if (mu - self.full[self.thin]['vel'][idx_left])<=(self.full[self.thin]['vel'][idx_right] - mu):
                #     idx_slice_thin = idx_left
                # else:
                #     idx_slice_thin = idx_right

                # img_thin = gaussian_filter(self.Dmol[self.thin][idx_slice_thin], sigma=sigma_filter)
                # img_thick = gaussian_filter(self.Dmol[self.thick][idx_slice_thick], sigma=sigma_filter)
                # T13, T12 = self.get_line_temp(img_thin, img_thick, method=method)
                
                # # Get normalization factor               
                # # Get real line value
                # max_thin = self.full[self.thin]['line'][idx_slice_thin]
                # max_thick = self.full[self.thick]['line'][idx_slice_thick]

                # # Normalise thin molecule
                # T13 = T13*A/max_thin
                # # Normalise thick molecule
                # P_thick = gaussian(self.full[self.thick]['vel'][idx_slice_thick], A_thick, mu_thick, sigma_thick)
                # T12 = T12*P_thick/max_thick
                # ------------------------------------------------------

                # Assign the peak main beam temperature
                T13 = A         # Temperature of thin molecule
                T12 = A_thick   # Temperature of thick molecule

                T13s[idx] = T13 
                T12s[idx] = T12

                # FWHM correction
                # ------------------------------------------------------
                # Virial FWHM
                fwhm_vir = k_sigma*fwhms_vir[idx]
                fwhm_vir_km = fwhm_vir/f
                # LTE FWHM
                fwhm = k_sigma*fwhm
                fwhm_km = fwhm/f

                # Virial mass
                # ------------------------------------------------------
                MVIR, Re = mass_virial(params, fwhm_vir_km)
                # Accumulate masses
                MVIR_a += MVIR

                # LTE mass
                # ------------------------------------------------------
                if flag_match_lines:
                    # Line temperatures
                    T12CO = k*T12
                    T13CO = k*T13
                    # Get mass
                    MLTE, N_13CO, Re = mass_lte(params, T12CO, T13CO, v_thick, v_thin, fwhm_km, A10_thin, X_thin, J)
                    # Accumulate masses
                    MLTE_a += MLTE
                    # Column densities data
                    N_13CO, Tex, tau = column_density(T12CO, T13CO, v_thick, v_thin, fwhm_km, A10_thin, J)
                    Texs[idx] = Tex
                    N_H2 = get_nh2(self.thin, N_13CO)
                    # Accumulate column density
                    N_13CO_a += N_13CO
                    N_H2_a += N_H2

                # X-factor
                # ------------------------------------------------------
                #kn = T13/max_thin      # ????? A/max_thin
                kn = 1
                line_temp = k*kn*gaussian(self.full[self.thin]['vel'], A, mu, sigma)
                # Get integration limits
                nsigma = 1
                lim_min = mu - nsigma*sigma
                lim_max = mu + nsigma*sigma

                # Get index limits
                line_min = np.where(lim_min>=self.full[self.thin]['vel'])[0]
                line_max = np.where(lim_max<=self.full[self.thin]['vel'])[0]

                if len(line_min) > 0:
                    idx_lim_min = line_min[-1]
                else:
                    idx_lim_min = 0

                if len(line_max) > 0:
                    idx_lim_max = line_max[0]
                else:
                    idx_lim_max = -1

                vel_span = self.full[self.thin]['vel'][idx_lim_min:idx_lim_max]/f
                line_span = line_temp[idx_lim_min:idx_lim_max]

                if len(line_span) > 0:
                    MXF, Re = mass_xf(params, vel_span, line_span)
                else:
                    MXF, Re = 0., 0.

                MXF_a += MXF

            self.full[self.thin]['T'] = T13s
            self.full[self.thick]['T'] = T12s

            # Get masses
            # ------------------------------------------------------
            if flag_match_lines:
                self.full['mass_lte'] = MLTE_a
            else:
                msg('Mass LTE fail', 'fail')

            self.full['mass_vir'] = MVIR_a
            self.full['mass_xf'] = MXF_a
            msg('Masses done', 'ok')
            
            # Column density
            # ------------------------------------------------------
            if flag_match_lines:
                # Asign column density
                self.full['N13CO'] = N_13CO_a
                self.full['NH2'] = N_H2_a
                self.full['tau'] = tau
                msg('Column densities done', 'ok')
            else:
                msg('Column densities fail', 'fail')

            # Excitation temperatures
            # ------------------------------------------------------
            if flag_match_lines:
                self.full['Tex'] = Texs
                msg('Excitation temperatures done', 'ok')
            else:
                msg('Excitation temperatures fail', 'fail')

            # Effective radius
            # ------------------------------------------------------
            # Convert to the used Systems Units
            self.full['Re'] = Re
            msg('Excitation Temperature and Effective Radius done', 'ok')

        else:
            msg(self.thin+' or '+self.thick+' line data are not available, needed to get the mass', 'fail')


    def binning_mol(self, mol, nbins, rebin=False):
        """
            Square binnig of cube data with n bins
            Parameters
            ----------
            mol : string
                Molecule
            nbins : int
                Number of bins
            rebin : boolean
                Renew the binnig?
            ----------
        """

        if mol in self.mols:
            if self.Dmol[mol] == []:
                msg(mol + ' Molecule data not available', 'fail')
                return
            # Parameters
            cdelt = self.Hmol[mol]['CDELT*']
            # Dimensions
            data = self.Dmol[mol]
            # Get volume
            volume = self.full['vol']
        else:
            msg('Molecule not valid', 'fail')
            return

        # Velocity maps
        v, y, x = np.shape(self.Dmol[mol])
        xlim, ylim = x/nbins, y/nbins

        # Check if the rebinning is possible
        if xlim < 1.0 or ylim < 1.0:
            msg('Binning limit exceeded.\nMax binning is at least 1 pixel per bin',\
            'fail')
            return

        if (xlim > 1 and xlim < 2) or (ylim > 1 and ylim < 2):
            msg('Binning number not valid. It will create bins of zero dimension.',\
            'fail')
            max_bin = np.max((x, y))
            msg('Choose a bin number between 2-'+str(round(max_bin/2))+' or the maximum binning: '+str(max_bin), 'info')
            return

        msg('Starting binning for molecule '+mol, 'info')

        # Rebinning
        if rebin:
            self.binned = {}
            self.binned['B0'] = {}
            for m in self.mols:
                self._init_binned_data(m)
            # Get the bins information
            bin_data, bin_volume, self.bin_grid, self.bin_frames = self.binning(data, volume, nbins, cdelt)
            msg('Starting rebinning', 'ok')
        else:
            self._init_binned_data(mol)
            # Check if bin frames exists
            if len(self.bin_frames) > 0:
                # Get binned data
                bin_data = self.binning(data, volume, nbins, cdelt, frames_pts=self.bin_frames)
                msg('Binning for molecule:'+mol, 'ok')
            else:
                msg('Bins not defined. Rebin first', 'fail')
                return

        self.nbins = nbins

        # Asign data per bin
        cnt = 0
        for i in range(nbins):
            for j in range(nbins):
                flag_binning = True
                name = 'B'+str(cnt)
                if not name in self.binned:
                    if rebin:
                        self.binned[name] = {}
                        msg('Bin '+str(cnt)+' created', 'ok')
                    else:
                        msg('Bin '+name+' does not exist. Rebinning needed', 'warn')
                        flag_binning = False

                if flag_binning:
                    # Start binning flagging
                    self.binned[name]['flag'] = False

                    # Define bin coordinates
                    self.binned[name]['pos'] = [i, j]

                    # Asign volume if rebinning
                    if rebin:
                        self.binned[name]['vol'] = bin_volume[i][j]

                    # Create the bin for the molecule
                    if not mol in self.binned[name]:
                        self.binned[name][mol] = {}

                    # Get velocity and intensity of the line
                    self.binned[name][mol]['vel'] = self.full[mol]['vel']
                    self.binned[name][mol]['line'] = bin_data[i][j]
                    
                    # Get restframe frequency
                    #rest_freqs_keys = self.Hmol[mol]['RESTFR*']
                    #rf_key = list(rest_freqs_keys.keys())[0]
                    #self.binned[name][mol]['freq'] = self.Hmol[mol][rf_key]
                    self.binned[name][mol]['freq'] = self.full[mol]['freq']

                    # Define the baseline bin_grid = np.zeros((nbins**2, 4), dtype=int)
                    self.binned[name][mol]['bl'] = np.array([])

                    cnt += 1

        msg('Binning complete', 'info')


    def get_line_temp(self, img_thin, img_thick, **kwargs):
        """
            Get the slice with the maximum value of integration
            Parameters
            ----------
            img_thin : np.array
                Optically thin image
            img_thick : np.array
                Optically thick image
            ----------
        """
        # Key arguments
        # ----------------------------------------------
        # Defining some params
        centre_square = kwargs.pop('centre_square', 0.25)
        # Default region for binned maps
        bounds = kwargs.pop('bounds', None)
        # Frame fractional pixels
        frames = kwargs.pop('frames', np.array([0]))
        # Method
        method = kwargs.pop('method', 'max')
        # ----------------------------------------------

        # Full map
        if bounds is None:
            # Thick and thin images
            img_opt_thick = img_thick
            img_opt_thin = img_thin
        # Binned map
        else:
            # Thick and thin images
            try:
                if len(frames) > 0:
                    frame_thick = np.zeros(len(frames)) 
                    frame_thin = np.zeros(len(frames))
                    for i, f in enumerate(frames):
                        xp = f[0]
                        yp = f[1]
                        wp = f[2]
                        # Define the frame thickness
                        frame_thick[i] = img_thick[yp, xp]*wp
                        frame_thin[i] = img_thin[yp, xp]*wp
                else:
                    frame_thick = np.array([0])
                    frame_thin = np.array([0])

                img_opt_thick = img_thick[bounds[2]:bounds[3], bounds[0]:bounds[1]]
                img_opt_thin = img_thin[bounds[2]:bounds[3], bounds[0]:bounds[1]]
            except:
                msg('Boundaries not valid', 'fail')
                return

        # Get the temperatures
        # Optically thin
        # Get temperature
        if bounds is None:
            if method == 'max':
                T13CO = np.nanmax(img_opt_thin)
            elif method == 'mean':
                T13CO = np.nanmean(img_opt_thin)
        else:
            if method == 'max':
                T13CO = np.nanmax((np.nanmax(frame_thin), np.nanmax(img_opt_thin)))
            elif method == 'mean':
                T13CO = np.nanmean((np.nanmean(frame_thin), np.nanmean(img_opt_thin)))

        # Optically thin
        # Get position
        nrows, ncols = img_opt_thick.shape
        # Has to be the maximum form the thick molecule

        # Check if there is data no nan available
        img_opt_thin_no_nan = img_opt_thin[~np.isnan(img_opt_thin)]

        if len(img_opt_thin_no_nan) > 0:
            flat_max_idx = np.nanargmax(img_opt_thin)
           
            xpos, ypos = int(flat_max_idx/ncols), flat_max_idx%ncols
            # Criterion of maximum values
            x_sqr = int(centre_square*img_opt_thick.shape[0])
            y_sqr = int(centre_square*img_opt_thick.shape[1])
            # Define limits
            x_min = xpos - x_sqr
            if x_min < 0:
                x_min = 0
            x_max = xpos + x_sqr
            if x_max >= nrows:
                x_max = nrows
            y_min = ypos - y_sqr
            if y_min < 0:
                y_min = 0
            y_max = ypos + y_sqr
            if y_max >= ncols:
                y_max = ncols

            n_img_thick = img_opt_thick[x_min:x_max, y_min:y_max]

            # Get temperature
            if bounds is None:
                if method == 'max':
                    T12CO = np.nanmax(n_img_thick)
                elif method == 'mean':
                    T12CO = np.nanmean(n_img_thick)
            else:
                if method == 'max':
                    T12CO = np.nanmax((np.nanmax(frame_thick), np.nanmax(n_img_thick)))
                elif method == 'mean':
                    T12CO = np.nanmean((np.nanmean(frame_thick), np.nanmean(n_img_thick)))

        else:
            T13CO, T12CO = np.nan, np.nan

        return T13CO, T12CO


    def line_fit_binning(self, mol, inter=False, **kwargs):
        """
            Fit the line profile of the whole map
            Parameters
            ----------
            mol : string
                Molecule
            inter : boolean
                Interactive window?
            ----------
        """
        # Key arguments
        # ----------------------------------------------
        # Choose bin to fit
        nbin = kwargs.pop('nbin', None)
        # Channels
        chns = kwargs.pop('chns', (5,5))
        # ----------------------------------------------

        # Force n lines according with the number of lines of the thin molecule
        # This forced flag has higher hierachy than the forced lines
        forced_thin = kwargs.pop('forced_thin', False)
        # Force n lines
        forced_lines = kwargs.pop('forced_lines', None)
        if not forced_lines is None:
            forced_lines = int(forced_lines)
            if forced_lines < 0:
                forced_lines = 0

        # Flag lines
        flag_maps = kwargs.pop('flag_maps', True)
        if flag_maps:
            if len(self.noise[mol]) < 2:
                msg('Not enough noise data selected', 'warn')
                flag_maps = False

        # Sigma threshold
        sigma_thresh = kwargs.pop('sigma_thresh', 1)

        # Sigma noise
        sigma_type = kwargs.pop('sigma_type', 'rms')
        if not sigma_type in ['rms', 'std']:
            msg('Sigma type is not valid. RMS assigned by default', 'warm')
            sigma_type = 'rms'

        # Verbose
        verbose = kwargs.pop('verbose', False)

        if mol in self.mols:
            # Check if it is a single bin or the full set
            if nbin:
                names = [nbin]
            else:
                names = self.binned.keys()

            # Get noise constrains
            #if flag_maps:
                #if sigma_type == 'std':
                    #mean_noise = self.noise_stat[mol]['mean']
                    #std_noise = self.noise_stat[mol]['std']
                    #line_thresh = (mean_noise+sigma_thresh*std_noise)*(tot_bins/4)
                #else:
                    #rms_noise = self.noise_stat[mol]['rms']
                    #line_thresh = (sigma_thresh*rms_noise)*(tot_bins/4)

            for name in names:
                if len(self.binned[name][mol]['line']) == 0:
                    msg(mol + ' Molecule data not available', 'fail')
                    return

                # Spectra without baseline
                if len(self.binned[name][mol]['bl']) > 0:
                    spectra = self.binned[name][mol]['line'] - self.binned[name][mol]['bl']
                else:
                    spectra = self.binned[name][mol]['line']

                # Check if the map is not noise
                proceed = False
                if flag_maps:
                    # Get noise signal
                    noise_data = np.concatenate((self.binned[name][mol]['line'][:chns[0]], self.binned[name][mol]['line'][-1*chns[1]:]))
                    if sigma_type == 'std':
                        mean_noise = np.nanmean(noise_data)
                        std_noise = np.nanstd(noise_data)
                        line_thresh = mean_noise+sigma_thresh*std_noise
                    else:
                        rms_noise = np.sqrt(np.nanmean(noise_data**2))
                        line_thresh = sigma_thresh*rms_noise

                    above_thresh_points = np.where(np.abs(spectra) > line_thresh)[0]
                    if len(above_thresh_points) > 0:
                        # Detect spikes
                        diff_above_thresh = np.diff(above_thresh_points)
                        next_above_points = np.where(diff_above_thresh == 1)[0]
                        if len(next_above_points) > 0:
                            proceed = True
                else:
                    proceed = True

                if proceed:
                    if flag_maps:
                        self.binned[name]['flag'] = False
                    else:
                        self.binned[name]['flag'] = self.binned[name]['flag'] or False 
                    # Forced the same number of lines than the thin molecule
                    if forced_thin:
                        #forced_lines = 0
                        if 'A' in self.binned[name][self.thin].keys():
                            forced_lines = len(self.binned[name][self.thin]['A'])

                    # If the number of lines are forced
                    if forced_lines is None:
                        # Find peaks and get an guess estimation
                        peaks = find_profile_peaks(spectra, **kwargs)
                    else:
                        if forced_lines == 0:
                            peaks = []
                        elif forced_lines == 1:
                            peaks = [np.argmax(spectra)]
                        else:
                            item_spectra = np.arange(len(spectra))
                            ndiv = int(np.ceil(len(spectra)/(2+forced_lines)))
                            peaks_raw = item_spectra[::ndiv]
                            peaks = peaks_raw[1:-1]
                else:
                    self.binned[name]['flag'] = self.binned[name]['flag'] or True
                    peaks = []
                           
                if len(peaks) > 0:
                    # Verbose for peaks
                    if verbose:
                        print("========== Peaks founded ==========")
                        print("Peaks founded for ", name, " :", len(peaks))

                    popt, pcov = solve_line(self.binned[name][mol]['vel'], spectra, peaks, ['G']*len(peaks))
                    A, mu, sigma, fwhm = get_params_from_popt(popt, len(peaks))
                else:
                    pcov = []
                    A, mu, sigma, fwhm = [], [], [], []
                    msg('No lines found in '+name, 'warn')

                # Sort parameters through the mu position
                mu_sort = np.argsort(mu)
                A = [A[i] for i in mu_sort]
                mu = [mu[i] for i in mu_sort]
                sigma = [sigma[i] for i in mu_sort]
                fwhm = [fwhm[i] for i in mu_sort]

                # Save line parameters
                self.binned[name][mol]['A'] = A
                self.binned[name][mol]['mu'] = mu
                self.binned[name][mol]['sigma'] = sigma
                self.binned[name][mol]['fwhm'] = fwhm
                self.binned[name][mol]['err'] = pcov

                # If the interactive mode is active
                if inter:
                    # Add canvas to the app
                    ioff()
                    fig, ax = subplots()
                    close(fig)

                    # Instance app
                    app = QApplication.instance()

                    self.FitSpectra = FitLinesWindow()
                    units = [self.full[mol]['units']['vel'], self.full[mol]['units']['flux']]
                    save = kwargs.pop('save', False)
                    name_bin = self.mcID + '-' + mol + '-' + name
                    self.FitSpectra.load_init_params(fig, ax, self.binned[name][mol]['vel'], spectra, name+'/'+mol, units, name_bin, save=save)
                    # Show the data pre-adjusted
                    if len(A) > 0:
                        self.FitSpectra.load_fit_params(popt, pcov, peaks, ['G']*len(peaks))
                    # Signal connection to extract baseline data
                    self.FitSpectra.signal_fitting.connect(self._get_fitting_from_ui)
                    self.FitSpectra.show()

                    app.exec_()

        else:
            msg('Molecule not available', 'fail')
            return


    def build_data_header(self, param, region):
        """
            Build data/header dict of a parameter
            Parameters
            ----------
            param : string
                Get parameter 
            region : string
                Full or binned region
            ----------
        """
        data_dict = {}

        # If param is in the full map
        if region == 'full':
            # Polarization param
            if param in ['I', 'Q', 'U', 'V', 'N', 'tau', 'pol_degree', 'pol_angle']:
                hd_pm = self.Hpol.copy()
            # Volume
            elif param in ['vol']:
                hd_pm = self.Hmol[self.thin].copy()
            else:
                msg('Parameter not valid', 'fail')
                return
            # Data
            data_dict['data'] = self.full[param]
            # Header
            items_to_del = hd_pm['*3']
            for i in items_to_del:
                del(hd_pm[i])
            data_dict['header'] = hd_pm

            return data_dict

        elif region == 'binned':
            if param in self.binned['B0'].keys():
                if param in ['pol_angle_mean', 'pol_angle_std', 'pol_degree_mean', 'pol_degree_std']:
                    hd_pm = self.Hpol.copy()  
                elif param in ['vol', 'mass_lte', 'mass_xf', 'mass_vir', 'B', 'den', 'N13CO', 'NH2']:
                    hd_pm = self.Hmol[self.thin].copy()
                else:
                    msg('Parameter not valid', 'fail')
                    return
            # Header
            items_to_del = hd_pm['*3']
            for i in items_to_del:
                del(hd_pm[i])

            # Get header dimensions
            m0_shape = self.Dmol[self.thin].shape[1:]
            dims_m0 = (m0_shape[1], m0_shape[0])
            # Get nbins dimensions
            bins = self.binned.keys()
            last_bit = 'B' + str(np.max([int(i[1:]) for i in bins]))
            sx, sy = self.binned[last_bit]['pos']
            dims_bins = (sx+1, sy+1)
            data_dict['header'] = [hd_pm, dims_m0, dims_bins]

            # Data
            data_dict['data'] = np.zeros((sx+1, sy+1))
            for b in self.binned.keys():
                pos = self.binned[b]['pos']
                data_dict['data'][pos[0]][pos[1]] = self.binned[b][param]

            return data_dict


    def binning(self, data, volume, nbins, cdelt, **kwargs):
        """
            Binning data
            Parameters
            ----------
            data : 3-d array
            volume: 3-d array
            nbins : int
                Number of bins as powers of two
            cdelt : list
                Step for the 3 axis
            ----------
        """
        # Key arguments
        # ----------------------------------------------
        # Bin frames definition
        frames_pts = kwargs.pop('frames_pts', [])
        # ----------------------------------------------

        # If frame points are defined
        rebin = False
        if len(frames_pts) == 0:
            rebin = True

        # Velocity maps
        v, y, x = np.shape(data)
        xlim, ylim = x/nbins, y/nbins

        # Velocity grid
        vel_bin = np.zeros((nbins, nbins, v))
        # Grid boundaries array
        bin_grid = np.zeros((nbins**2, 4), dtype=int)
        bin_frames = []
        # Volume grid
        bin_volume = np.zeros((nbins, nbins))

        add_vol_pix = 0
        add_data_pix = np.zeros(v)

        cell = 0
        for i in range(nbins):
            for j in range(nbins):
                # X-Y limits for each data bin
                x1 = i*xlim
                x2 = (i+1)*xlim
                y1 = j*ylim
                y2 = (j+1)*ylim
                pts = ((x1,x2),(y1,y2))
                # Get max and min limits
                xmin = int(np.ceil(x1))
                xmax = int(x2)
                ymin = int(np.ceil(y1))
                ymax = int(y2)

                lims = ((xmin,xmax),(ymin,ymax))
                # Extract fractional data
                if (x%nbins > 0) or (y%nbins > 0):
                    # If the frame points are not defined
                    if rebin:
                        frames_pts = self._fractional_binning(pts, lims, (x,y))
                        # Get the fractional volume
                        add_vol_pix = self._int_frac_pixels(frames_pts, volume)#[0]
                        add_data_pix = self._int_frac_pixels(frames_pts, data)
                    else:
                        # Get the fractional data
                        frame_data = frames_pts[cell]
                        add_data_pix = self._int_frac_pixels(frame_data, data)

                # If data is being rebinned
                if rebin:
                    # Get Volume
                    bin_volume[i][j] = np.sum(volume[xmin:xmax, ymin:ymax])+add_vol_pix
                    # Concatenate grid
                    bin_grid[cell] = np.array([xmin, xmax, ymin, ymax])
                    bin_frames.append(frames_pts)

                # Defining the sub-region
                bin_data = data[:, ymin:ymax, xmin:xmax]  

                if len(frames_pts)>0:
                    if rebin:
                        nframes = frames_pts
                    else:
                        nframes = frames_pts[cell]
                    nsum = 0
                    for nf in nframes:
                        nsum += nf[2]
                else:
                    nsum = 0
                # Get the integrated velocity
                vel_bin[i][j][:] = self.map_vel(v, bin_data, cdelt, frame=add_data_pix, n_frames=nsum)
                cell += 1

        if rebin:
            return vel_bin, bin_volume, bin_grid, bin_frames
        else:
            return vel_bin


    def _fractional_binning(self, pts, lims, map_size):
        """
            Do fractional binning
            Parameters
            ----------
            pts : tuple
            lims: tuple
            map_size : array
                Map dimensions
            ----------
        """
        # Extract points
        x1 = pts[0][0]
        x2 = pts[0][1]
        y1 = pts[1][0]
        y2 = pts[1][1]
        # Get full dimensions
        x = map_size[0]
        y = map_size[1]
        # Extract limits
        xmin = lims[0][0]
        xmax = lims[0][1]
        ymin = lims[1][0]
        ymax = lims[1][1]
        # Get rectangle dimensions
        w = xmax-xmin
        h = ymax-ymin
        dims = [w, h]*2
        # Get fractional part
        xdec_min = xmin-x1
        xdec_max = x2-xmax
        ydec_min = ymin-y1
        ydec_max = y2-ymax
        # Define borders
        xvs = [xmin-1, xmax, xmin-2, xmin-1]
        yvs = [ymin-1, ymin-1, ymax, ymin-2]

        wds = [ydec_min, xdec_max, ydec_max, xdec_min]
        # Define frame container
        frames_pts = []

        for m in range(4):
            xv = xvs[m]
            yv = yvs[m]
            # Check x limits
            for p in range(dims[m]+1):
                if m%2 == 0:                        
                    xv += 1
                else:
                    yv += 1 
                # Only the bits inside the map
                if (xv >= 0 and xv < x) and (yv >= 0 and yv < y):
                    if xv == xmax and yv == ymax:
                        wp = wds[1]*wds[2]
                    elif xv == xmax and yv == ymin-1:
                        wp = wds[1]*wds[0]
                    elif xv == xmin-1 and yv == ymin-1:
                        wp = wds[3]*wds[0]
                    elif xv == xmin-1 and yv == ymax:
                        wp = wds[3]*wds[2]
                    else:
                        wp = wds[m]                      
                    # Get the fractional volume
                    frames_pts.append((xv, yv, wp))

        return frames_pts


    def _int_frac_pixels(self, frame, data, get_frame=False):
        """
            Integrate fractional pixels weighted
            Parameters
            ----------
            frame : np.array 
            lims: np.array
            get_frame : boolean (optional)
            ----------
        """
        dims = data.shape
        if len(dims) == 2:
            z = 1
        else:
            z = dims[0]

        if get_frame:
            acum_data = []

        add_pix = np.zeros(z)
        for k in range(z):
            acum_pix = 0
            for p in range(len(frame)):
                x = frame[p][0]
                y = frame[p][1]
                w = frame[p][2]
                if z == 1:
                    val_pix = data[x, y]*w 
                else:
                    val_pix = data[k, y, x]*w

                if get_frame:
                    acum_data.append(val_pix)

                acum_pix += val_pix

            add_pix[k] = acum_pix

        if get_frame:
            return add_pix, acum_data
        else:
            return add_pix


    def get_bins_params(self, **kwargs):
        """
            Get parameters of each bin of the molecular cloud
        """

        # Key arguments
        # ----------------------------------------------
        # Sigma of Gaussian filter to get the lines temperatures from maps
        sigma_filter = kwargs.pop('sigma_filter', 1)
        # No gaussian filter applied
        # no_filter = kwargs.pop('no_filter', False)
        # Temperature threshold to be considered as a valid line
        #thresh = kwargs.pop('thresh', 3)
        # Flag lines
        flag_maps = kwargs.pop('flag_maps', True)
        # Define virial molecule
        mol_vir = kwargs.pop('mol_vir', 'thin')
        # Method to get radiation temperature
        method = kwargs.pop('method', 'max')
        # ----------------------------------------------

        # Main physical parameters
        dist = self.mc_params['dist']
        kRe = self.mc_params['Re_factor']
        # Factor to convert from rectangle to square
        dims_size = [self.lon_m[-1]-self.lon_m[0], \
                     self.lat_m[-1]-self.lat_m[0]]
        b_axis = np.max(np.abs(dims_size))
        a_axis = np.min(np.abs(dims_size))
        dims_factor = b_axis/a_axis
        kRe_square = kRe*np.sqrt(4*dims_factor/np.pi)
        # Angular size
        # theta = self.mc_params['ang_dia']/self.nbins
        min_theta = np.min((np.abs(self.full[self.thin]['lat'][-1]-self.full[self.thin]['lat'][0]),\
                    np.abs(self.full[self.thin]['lon'][-1]-self.full[self.thin]['lon'][0])))

        theta = min_theta/self.nbins

        # Put params together
        params = [dist, theta, kRe_square]

        msg('Calculating binned parameters...', 'info')

        # Get some molecular cloud parameters
        k = self.mc_params['temp_factor']
        k_sigma = self.mc_params['width_factor']

        # Get conversion to km/s
        f, _ = get_factors_units('km s^-1', su) 

        # Apply a filter for all the images
        # slices = self.Dmol[self.thin].shape[0]
        # filter_data_thin = np.zeros_like(self.Dmol[self.thin])
        # filter_data_thick = np.zeros_like(self.Dmol[self.thick])

        # for s in range(slices):
        #     if no_filter:
        #         filter_data_thin[s] = self.Dmol[self.thin][s]
        #         filter_data_thick[s] =self.Dmol[self.thick][s]
        #     else:   
        #         filter_data_thin[s] = gaussian_filter(self.Dmol[self.thin][s], sigma=sigma_filter)
        #         filter_data_thick[s] = gaussian_filter(self.Dmol[self.thick][s], sigma=sigma_filter)

        for name in self.binned.keys():
            # Check if bin is not flagged
            if self.binned[name]['flag']:
                msg('Flagged bin '+name, 'info')
            # Flag to calculate params
            flag_get_params = True
            # If there is a line for the 'name' bin
            if len(self.binned[name][self.thin]['line']) > 0 and len(self.binned[name][self.thick]['line']) > 0:
                
                # Get FWHM 
                if mol_vir == 'thin':
                    fwhms_vir = self.binned[name][self.thin]['fwhm']
                elif mol_vir == 'thick':
                    fwhms_vir = self.binned[name][self.thick]['fwhm']

                # Get FWHM for thin molecule
                fwhms = self.binned[name][self.thin]['fwhm']

                # Number of lines
                nline = len(fwhms)

                # Check if there are lines
                if nline == 0:
                    msg('No lines found in bin '+name, 'fail')
                    flag_get_params = False          

                # LTE mass and N13CO flags
                flag_match_lines = False
                if (nline == len(self.binned[name][self.thick]['A'])) and (nline == len(self.binned[name][self.thin]['A'])):
                    flag_match_lines = True
                else:
                    msg('Number of '+self.thin+' and '+self.thick+' lines are not the same at bin '+name+'.\nMass by LTE method and densities columns wont be calculated', 'warn')
                    flag_get_params = False

                # Get frequency molecules
                v_thick = self.full[self.thick]['freq']
                v_thin = self.full[self.thin]['freq']

                # Get Einstein A coeficient
                line_thin = self.full[self.thin]['t_line']
                for m in mol_params[self.thin].keys():
                    if line_thin in m:
                        A10_thin = mol_params[self.thin][m]['A10']
                        X_thin = mol_params[self.thin]['X']
                        J = int(line_thin[0])

                # Initialise physical parameters
                tau = 0
                MVIR_a, MLTE_a, MXF_a, N_13CO_a, N_H2_a, Re = 0, 0, 0, 0, 0, 0
                # Excitation temperature
                Texs = np.zeros(nline)

                T13s = np.zeros(nline)
                T12s = np.zeros(nline)

                if flag_get_params:
                    # Get bounds of the bin
                    num_bin = int(name[1:])
                    bounds = self.bin_grid[num_bin]
                    frames = self.bin_frames[num_bin]

                    # If there are peaks at the line
                    if len(fwhms) > 0:
                        # Get physical parameters of each line detected
                        for idx, fwhm in enumerate(fwhms):

                            # Line parameters thin molecule
                            # ------------------------------------------------------
                            A = self.binned[name][self.thin]['A'][idx]
                            mu = self.binned[name][self.thin]['mu'][idx]
                            sigma = self.binned[name][self.thin]['sigma'][idx]

                            # Line parameters thick molecule
                            # ------------------------------------------------------
                            A_thick = self.binned[name][self.thick]['A'][idx]
                            mu_thick = self.binned[name][self.thick]['mu'][idx]
                            sigma_thick = self.binned[name][self.thick]['sigma'][idx]

                            # Get line temperatures
                            # Use 12CO (optically thick molecule)
                            # ------------------------------------------------------
                            
                            # Code used to get T12 and T13, now it is obsolete
                            # ------------------------------------------------------
                            # Thin molecules
                            # line_left = np.where(mu>=self.binned[name][self.thin]['vel'])[0]
                            # if len(line_left) > 0:
                            #     idx_left = line_left[-1]
                            # else:
                            #     idx_left = 0

                            # idx_right = idx_left + 1
                            # if idx_right >= len(self.binned[name][self.thin]['vel']):
                            #     idx_right = idx_left

                            # if (mu - self.binned[name][self.thin]['vel'][idx_left]) <= (self.binned[name][self.thin]['vel'][idx_right] - mu):
                            #     idx_slice_thin = idx_left
                            # else:
                            #     idx_slice_thin = idx_right

                            # # Thick molecules
                            # line_left = np.where(mu>=self.binned[name][self.thick]['vel'])[0]
                            # if len(line_left) > 0:
                            #     idx_left = line_left[-1]
                            # else:
                            #     idx_left = 0

                            # idx_right = idx_left + 1
                            # if idx_right >= len(self.binned[name][self.thick]['vel']):
                            #     idx_right = idx_left

                            # if (mu - self.binned[name][self.thick]['vel'][idx_left]) <= (self.binned[name][self.thick]['vel'][idx_right] - mu):
                            #     idx_slice_thick = idx_left
                            # else:
                            #     idx_slice_thick = idx_right

                            # # Asign slice according to their maximum line profile
                            # img_thin = filter_data_thin[idx_slice_thin]
                            # img_thick = filter_data_thick[idx_slice_thick]

                            # # Check if there is data no nan available
                            # img_thin_no_nan = img_thin[~np.isnan(img_thin)]
                            # img_thick_no_nan = img_thick[~np.isnan(img_thick)]

                            # if len(img_thin_no_nan) > 0 and len(img_thick_no_nan) > 0:
                            #     T13, T12 = self.get_line_temp(img_thin, img_thick, bounds=bounds, frames=frames, centre_square=1.0, method=method)
                            # else:
                            #     T13, T12 = np.nan, np.nan

                            # # Get normalization factor               
                            # # Get real line value
                            # max_thin = self.binned[name][self.thin]['line'][idx_slice_thin]
                            # max_thick = self.binned[name][self.thick]['line'][idx_slice_thick]

                            # # Normalise thin molecule
                            # T13 = T13*A/max_thin
                            # # Normalise thick molecule
                            # P_thick = gaussian(self.binned[name][self.thick]['vel'][idx_slice_thick], A_thick, mu_thick, sigma_thick)
                            # T12 = T12*P_thick/max_thick
                            # ------------------------------------------------------

                            # Assign the peak main beam temperature
                            #if len(img_thin_no_nan) > 0 and len(img_thick_no_nan) > 0:
                            #    T13, T12 = self.get_line_temp(img_thin, img_thick, bounds=bounds, frames=frames, centre_square=1.0, method=method)
                            #else:
                            #    T13, T12 = np.nan, np.nan

                            T13 = A         # Temperature of thin molecule
                            T12 = A_thick   # Temperature of thick molecule

                            T12s[idx] = T12 
                            T13s[idx] = T13

                            # Get Max and Min temperatures
                            # ------------------------------------------------------
                            #T13COmax, T13COmin, T13COmean, T13COstd = self._get_slice_stats(self.thin, idx_slice_thin)
                            #T12COmax, T12COmin, T12COmean, T12COstd = self._get_slice_stats(self.thick, idx_slice_thick)

                            # If the line is over the threshold
                            #if (T13 > (T13COmean - thresh*T13COstd)) and (T12 > (T12COmean - thresh*T12COstd)):

                            # FWHM correction
                            # ------------------------------------------------------
                            # Virial FWHM
                            fwhm_vir = k_sigma*fwhms_vir[idx]
                            fwhm_vir_km = fwhm_vir/f
                            # LTE FWHM
                            fwhm = k_sigma*fwhm
                            fwhm_km = fwhm/f

                            # Virial mass
                            # ------------------------------------------------------
                            MVIR, Re = mass_virial(params, fwhm_vir_km)

                            # LTE mass
                            # ------------------------------------------------------
                            if flag_match_lines:
                                # Line temperatures
                                T12CO = k*T12
                                T13CO = k*T13
                                # Get mass
                                MLTE, N_13CO, Re = mass_lte(params, T12CO, T13CO, v_thick, v_thin, fwhm_km, A10_thin, X_thin, J)
                                if np.isnan(MLTE):
                                    print(name+" tau is negative. Thin/Thick temperature relation is too high")
                                # Column densities data
                                N_13CO, Tex, tau = column_density(T12CO, T13CO, v_thick, v_thin, fwhm_km, A10_thin, J)
                                Texs[idx] = Tex
                                N_H2 = get_nh2(self.thin, N_13CO)
                                N_13CO_a += N_13CO
                                N_H2_a += N_H2
                       
                            # X-factor
                            # ------------------------------------------------------
                            #kn = T13/max_thin      # ?????
                            kn = 1
                            line_temp = k*kn*gaussian(self.binned[name][self.thin]['vel'], A, mu, sigma)
                            # Get integration limits
                            nsigma = 1
                            lim_min = mu - nsigma*sigma
                            lim_max = mu + nsigma*sigma

                            # Get index limits
                            line_min = np.where(lim_min>=self.binned[name][self.thin]['vel'])[0]
                            line_max = np.where(lim_max<=self.binned[name][self.thin]['vel'])[0]

                            if len(line_min) > 0:
                                idx_lim_min = line_min[-1]
                            else:
                                idx_lim_min = 0

                            if len(line_max) > 0:
                                idx_lim_max = line_max[0]
                            else:
                                idx_lim_max = -1

                            vel_span = self.binned[name][self.thin]['vel'][idx_lim_min:idx_lim_max]/f
                            line_span = line_temp[idx_lim_min:idx_lim_max]

                            if len(line_span) > 0:
                                MXF, Re = mass_xf(params, vel_span, line_span)
                            else:
                                MXF, Re = 0., 0.

                            # else:
                            #     # If the threshold is not reached, all the parameters are null
                            #     # Masses
                            #     MVIR = 0.0
                            #     MLTE = 0.0
                            #     MXF = 0.0
                            #     # Column density
                            #     N_13CO = 0.0
                            #     N_H2 = 0.0

                            #     mols_lows = ''
                            #     if not T13 > (T13COmean - thresh*T13COstd):
                            #         mols_lows += self.thin
                            #     if not T12 > (T12COmean - thresh*T12COstd):
                            #         if mols_lows:
                            #             mols_lows += ','+self.thick
                            #         else:
                            #             mols_lows += self.thick

                            #     msg(name+' Line '+mols_lows+' too low to be considered', 'warn')

                            MVIR_a += MVIR
                            MLTE_a += MLTE
                            MXF_a += MXF
                            N_13CO_a += N_13CO
                            N_H2_a += N_H2

                        self.binned[name][self.thin]['T'] = T13s
                        self.binned[name][self.thick]['T'] = T12s

                    else:
                        msg(name+' no '+self.thin+' line', 'warn')                    

                else:
                    MVIR_a, MLTE_a, MXF_a, N_13CO_a, N_H2_a, Re = np.nan, np.nan, np.nan, np.nan, np.nan, np.nan
                    if flag_maps:
                        self.binned[name]['flag'] = True

                # Get masses
                # ------------------------------------------------------
                self.binned[name]['mass_lte'] = MLTE_a
                self.binned[name]['mass_vir'] = MVIR_a
                self.binned[name]['mass_xf'] = MXF_a
                # Column density
                # ------------------------------------------------------
                self.binned[name]['N13CO'] = N_13CO_a
                self.binned[name]['NH2'] = N_H2_a
                self.binned[name]['tau'] = tau
                # Effective radius
                # ------------------------------------------------------
                self.binned[name]['Re'] = Re
                # Excitation temperature
                # ------------------------------------------------------
                self.binned[name]['Tex'] = Texs

            else:
                msg('13CO and 12CO lines data not available, needed to get the mass', 'fail')

        msg('Done', 'ok')


    def get_mag_field(self, mol, region='binned', **kwargs):
        """
            Get magnetic field through Chandrasekhar-Fermi method
            Parameters
            ----------
            mol : string
                Molecule
            region : string
                Full or binned region 
            ----------
        """

        # Key arguments
        # ----------------------------------------------
        # Choose bin to fit
        nbin = kwargs.pop('nbin', None)
        # ----------------------------------------------

        msg('Calculating magnetic field (Chandrasekhar-Fermi method)', 'info')

        # Get conversion factor to km/s
        f, _ = get_factors_units('cm s^-1', su)

        # Full map
        if region == 'full':
            # Check if data is enough
            for i in ['den', 'pol_angle']:
                if not i in self.full.keys():
                    msg('Parameter '+i+' missing', 'fail')
                    return
            if len(self.full[mol]['sigma']) == 0:
                msg('Velocity line fitting missing', 'fail')
                return

            # Get density
            den = self.full['den']
            # Get polarization deviation
            dTheta = np.std(self.full['pol_angle'])
            # Apply DCF
            B = 0
            for dV in self.full[mol]['sigma']:
                #dV = self._get_fwhm(self.full[mol]['vel'], self.full[mol]['line'])
                #dV = dV/(2*np.sqrt(2*np.log(2)))
                # Get Magnetic Field
                B += magnetic_field(den, dV/f, dTheta)

            self.full['B'] = B

        # Binned map
        elif region == 'binned':

            if nbin:
                if isinstance(nbin, str):
                    bins = [nbin]
                elif isinstance(nbin, list):
                    bins = nbin
            else:
                bins = self.binned.keys()

            # Check bin by bin
            for b in bins:
                flag_bin = True
                # Check if data is enough
                for i in ['den', 'pol_angle_std']:
                    if not i in self.binned[b].keys():
                        msg(b+' parameter '+i+' missing', 'warn')
                        flag_bin = False 
                if 'sigma' in self.binned[b][mol].keys():
                    if len(self.binned[b][mol]['sigma']) == 0:
                        msg(b+' velocity line fitting missing', 'fail')
                        flag_bin = False
                else:
                    msg(b+' velocity line fitting missing', 'fail')
                    flag_bin = False              

                if flag_bin:
                    # Get density
                    den = self.binned[b]['den']
                    # Get polarization deviation
                    dTheta = self.binned[b]['pol_angle_std']
                    # Apply DCF
                    B = 0
                    for dV in self.binned[b][mol]['sigma']:
                        #dV = self._get_fwhm(self.binned[kid][mol]['vel'], self.binned[kid][mol]['line'])
                        #dV = dV/(2*np.sqrt(2*np.log(2)))
                        # Get Magnetic Field
                        B += magnetic_field(den, dV/f, dTheta)
                else:
                    B = np.nan

                self.binned[b]['B'] = B


        msg('Magnetic field (Perpendicular component) calculated', 'ok')


    def _get_fwhm(self, vel, line):
        """
            Get full-width at half maximum
            Parameters
            ----------
            vel : array
            line : array
            ----------
        """
        # Limits
        lim_min = np.min(line)
        lim_max = np.max(line)
        # Half power
        hm = (lim_max-lim_min)/2.
        # Full-width
        fw = vel[line>hm]
        if len(fw) > 0:
            fwhm = fw[-1]-fw[0]
        else:
            fwhm = 0

        return fwhm


    def sum_binning_param(self, param):
        """
            Sum of the values of one parameter per bin
            Parameters
            ----------
            param : string
                Bin parameter
            ----------
        """
        if not param in sos.all_params:
            msg('Parameter not valid', 'fail')
            return

        m = 0.
        for nbin in self.binned.keys():
            if not self.binned[nbin]['flag']:
                m_bin = self.binned[nbin][param]
                if not np.isnan(m_bin):
                    m += m_bin

        msg('Total bin sum of '+param+': '+str(m), 'ok')

        return m


    def get_vol_den(self, mass_method='mass_lte', region='binned', **kwargs):#, units=None):
        """
            Get volumetric density.
            Parameters
            ----------
            mass_method : string
                Mass method used: LTE, Virial or X-Factor
            region : string
                Molecular cloud region: full or binned
            ----------
        """

        # Key arguments
        # ----------------------------------------------
        # Choose bin to fit
        nbin = kwargs.pop('nbin', None)
        # ----------------------------------------------

        # Factor conversion from pc to m
        PC2M, _ = get_factors_units('pc', su)

        msg('Calculating Volumetric densities', 'info')

        if not mass_method in ['mass_lte', 'mass_xf', 'mass_vir']:
            msg('Method for the mass not recognized', 'fail')
            return

        if not region in ['full', 'binned']:
            msg('Region not valid', 'fail')
            return

        if region == 'full':
            self.full['den'] = self.full[mass_method]/self.full['vol']
            # Convert to g/cm^3 units 
            self.full['den'] = self.full['den']*1e3*Ms/((1e2*PC2M)**3)

        elif region == 'binned':

            # Check bins to map
            if nbin:
                if isinstance(nbin, str):
                    bins = [nbin]
                elif isinstance(nbin, list):
                    bins = nbin
            else:
                bins = self.binned.keys()

            for b in bins:
                self.binned[b]['den'] = self.binned[b][mass_method]/self.binned[b]['vol']
                # Convert to g/cm^3 units 
                self.binned[b]['den'] = self.binned[b]['den']*1e3*Ms/((1e2*PC2M)**3)

        msg ('Done', 'ok')


    def summary(self):
        """
            Show the results for the general analysis
        """
        msg ('------------------------------------------------', 'info')
        msg ('            Summary Molecular Cloud             ', 'info')
        msg ('------------------------------------------------', 'info')
        print ('MC ID: ', self.mcID)
        print ('Coupling telescope factor: ', self.mc_params['temp_factor'])
        print ('Distance: ', self.mc_params['dist'], ' pc')
        print ('Angular diameter: ', self.mc_params['ang_dia'])
        print ('Width line factor: ', self.mc_params['width_factor'])
        print ('Effective radius factor: ', self.mc_params['Re_factor'])
        
        if 'mass_lte' in self.full.keys():
            msg ('------------------------------------------------', 'info')
            msg ('              Physical properties               ', 'info')
            msg ('------------------------------------------------', 'info')
            msg ('MASS', 'ok')
            print ('LTE mass: {0:.2f} Msun'.format(self.full['mass_lte']))
            print ('Virial mass: {0:.2f} Msun'.format(self.full['mass_vir'], 'Ms'))
            print ('XF mass: {0:.2f} Msun'.format(self.full['mass_xf'], 'Ms'))
            msg ('COLUMN DENSITY', 'ok')
            print ('N13CO: {0:1.2e} cm^-2'.format(self.full['N13CO']))
            print ('NH2: {0:1.2e} cm^-2'.format(self.full['NH2']))
            msg ('OPACITY OF THIN MOLECULE', 'ok')
            print ('tau('+self.thin+'): {0:1.2e}'.format(self.full['tau']))
            msg ('EXCITATION TEMPERATURE', 'ok')
            for t in self.full['Tex']:
                print ('T1: {0:.2f} K'.format(t))
            msg ('DIMENSION', 'ok')
            print ('Effective radius: {0:.2f} pc'.format(self.full['Re']))


    def param_filter(self, param, max_lim=None, min_lim=None):
        """
            Parameter filter file
            Parameters
            ----------
            param : string
                Parameter to filter
            max, min : floats
                Minimum and maximum value to filter
            ----------
        """
        if not param in sos.all_params:
            msg('Parameter not valid', 'fail')
            return

        for kid in self.binned.keys():
            if max_lim and min_lim:
                if not (self.binned[kid][param] >= min_lim and self.binned[kid][param] <= max_lim):
                    self.binned[kid]['flag'] = True
            elif max_lim:
                if not self.binned[kid][param] <= max_lim:
                    self.binned[kid]['flag'] = True
            elif min_lim:
                if not self.binned[kid][param] >= min_lim:
                    self.binned[kid]['flag'] = True

        msg('Done', 'ok')


    def backup(self, name=None):
        """
            Backup the data
            Parameters
            ----------
            name : string
                Name of the backup. If none, it is named with the current date and time
            ----------
        """
        if not name:
            now = datetime.now()
            current_time = now.strftime("%d%m%y_%H%M%S")
            name = './'+self.mcID+'_'+current_time+'_bkp'

        os.mkdir(name)

        # Save the full and binned data
        np.save(name+'/mc_data.npy', [self.full, self.binned])
        msg('Saving molecular cloud full and binned data', 'ok')
        # Save header+data available
        np.save(name+'/mc_data_header.npy', [self.Dmol, self.Hmol])
        msg('Saving molecular cloud data and header', 'ok')
        # Save the polarization data
        np.save(name+'/mc_pol_data.npy', [self.Dpol, self.Hpol])
        msg('Saving stokes params for molecular cloud', 'ok')
        # Save bin bounds and frames
        np.save(name+'/mc_bounds.npy', self.bin_grid)
        np.save(name+'/mc_frames.npy', self.bin_frames)
        msg('Saving bin bounds and frames', 'ok')
        # Save parameters
        np.save(name+'/mc_params.npy', [self.nbins, self.xc, self.yc, self.R])
        msg('Saving molecular cloud parameters', 'ok')


    def load_bkp(self, dirfile):
        """
            Load backup data
            Parameters
            ----------
            dirfile : string
                Directory file of the backup to load
            ----------
        """
        # Load full + binned data
        full, binned = np.load(dirfile+'mc_data.npy', allow_pickle=True)
        self.full = full
        self.binned = binned
        msg('Molecular cloud full and binned data loaded', 'ok')
        # Load header+data available
        Dmol, Hmol = np.load(dirfile+'mc_data_header.npy', allow_pickle=True)
        self.Dmol = Dmol
        self.Hmol = Hmol
        msg('Molecular cloud data and header loaded', 'ok')
        # Load polarization data
        Dpol, Hpol = np.load(dirfile+'mc_pol_data.npy', allow_pickle=True)
        self.Dpol = Dpol
        self.Hpol = Hpol
        msg('Molecular cloud polarization data loaded', 'ok')
        # Load bin bounds and frames
        grids = np.load(dirfile+'mc_bounds.npy', allow_pickle=True)
        frames = np.load(dirfile+'mc_frames.npy', allow_pickle=True) 
        self.bin_grid = grids
        self.bin_frames = frames
        msg('Molecular bin bounds and frames loaded', 'ok')
        # Load header+data available
        nbins, xc, yc, R = np.load(dirfile+'mc_params.npy', allow_pickle=True)
        self.nbins = nbins
        self.xc = xc
        self.yc = yc
        self.R = R
        msg('Molecular cloud parameters loaded', 'ok')


    def _init_binned_data(self, mol):
        """
            Initiate binning data
            Parameters
            ----------
            mol : string
                Molecule to bin
            ----------
        """
        # Initiate molecular bin
        self.binned['B0'][mol] = {}
        self.binned['B0'][mol]['T'] = 0.0
        self.binned['B0'][mol]['vel'] = np.array([])
        self.binned['B0'][mol]['line'] = np.array([])
        self.binned['B0'][mol]['freq'] = 0.0


    def _check_axes(self, ax, name_list):
        """
            Check where axes belongs
            Parameters
            ----------
            ax : axes object
            name_list : list
                Parameter list
            ----------
        """
        for i in name_list:
            if i in ax.lower():
                return True

        return False


    def _header_pol(self, header):

        cdelt = header['CDELT*']
        crpix = header['CRPIX*']
        crval = header['CRVAL*']
        ctype = header['CTYPE*']

        # Create fits header
        w = wcs.WCS(naxis=2)
        w.wcs.crpix = [header['NAXIS1']/2, header['NAXIS2']/2]
        w.wcs.crval = [crval[0], crval[1]]
        w.wcs.cdelt = np.array([cdelt[0], cdelt[1]])
        w.wcs.ctype = [ctype[0], ctype[1]]
        header_pol = w.to_header()

        return header_pol


    def _header_vel(self, header):
        """
            Create the velocity header
            Parameters
            ----------
            header : fits header
                Header seed to create the new velocity header
            ----------
        """
        cdelt = header['CDELT*']
        crpix = header['CRPIX*']
        crval = header['CRVAL*']
        ctype = header['CTYPE*']

        # Create fits header
        w = wcs.WCS(naxis=2)
        w.wcs.crpix = [crpix[0], crpix[1]]
        w.wcs.crval = [crval[0], crval[1]]
        w.wcs.cdelt = np.array([cdelt[0], cdelt[1]])
        w.wcs.ctype = [ctype[0], ctype[1]]
        header_vel = w.to_header()

        return header_vel


    def _get_slice_stats(self, mol, n):

        # Define T vector
        T = np.zeros(4)

        T[0] = np.nanmax(self.Dmol[mol][n])
        # Get min
        T[1] = np.nanmin(self.Dmol[mol][n])
        # Get mean
        T[2] = np.nanmean(self.Dmol[mol][n])
        # Get standard deviation
        T[3] = np.nanstd(self.Dmol[mol][n])

        return T


    def param_visualiser(self, mol, n=0):
        """
            Interactive window to visualise parameters calculated
            Parameters
            ----------
            mol : string
                Molecule for the moment maps
            n : int
                n-Moment 
            ----------
        """

        print('Initialising molecular cloud parameters visualiser...')

        # Add canvas to the app
        ioff()
        fig, ax = subplots()
        close(fig)

        # Instance app
        app = QApplication.instance()

        pm_vis = ParamVisualiser()
        mn_mol = self.get_n_moment(mol, n=n)
        pm_vis.load_fits_file(mol, mn_mol, self.binned, self.mcID)
        pm_vis.show()

        app.exec_()