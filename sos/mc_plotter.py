# -*- coding: utf-8 -*-
# --------------------------------------------------------------------------------- #
# Software de Observaciones Sintéticas" S.O.S.
# Visualization tool
#
# Marcial Becerril, @ 28 August 2020
# Latest Revision: 28 Aug 2020, 22.34 GMT
#
# For all kind of problems, requests of enhancements and bug reports, please
# write to me at:
#
# mbecerrilt92@gmail.com
# mbecerrilt@inaoep.mx
#
# --------------------------------------------------------------------------------- #

import numpy as np
from matplotlib.pyplot import *
#ion()

import sos

from astropy import wcs
from astropy.io import fits
from astropy.wcs import WCS

from .line_fitting import *
from .misc.print_msg import *
from .misc.units_tool import *
from .init_vals import *
from .misc.stats import *


def plot_moment_spectra(Mn, mc, contour=True, fit=True, label=False, **kwargs):
    """
        Plot momentum zero map with lines profiles over.
        Parameters
        ----------
        Mn : dict
            n-Momentum map as a dictionary with two keys: 
            header and data.
        mc: dict
            Molecular cloud from a extracted molecule.
        contour: boolean (optional)
            Add M0 contour
        fit : boolean (optional)
            Plot the fitted profile?
        label: string (optional)
            Label each bin
        ----------
    """

    # Key arguments
    # ----------------------------------------------
    # Figure size
    figsize = kwargs.pop('figsize', (9,9))
    # Level contours
    level_contours = kwargs.pop('level_contours', 5)
    # Level contours
    color_contours = kwargs.pop('color_contours', 'black')
    # Alpha contours
    alpha_contours = kwargs.pop('alpha_contours', 0.5)
    # Line color
    line_color = kwargs.pop('line_color', 'c')
    # Line color
    fit_color = kwargs.pop('fit_color', 'white')
    # Text color
    text_color = kwargs.pop('text_color', 'white')
    # ----------------------------------------------

    # Get binned or full dimensions 
    dims, nbins, mc_binned = _get_mc_dims(mc)
    xbins, ybins = dims

    # Check if fit data is available
    flag_fit = True
    if fit:
        if mc_binned:
            if not 'A' in mc['B0'].keys():
                msg('Fitting data not available', 'warn')
                flag_fit = False
        else:
            if not 'A' in mc.keys():
                msg('Fitting data not available', 'warn')
                flag_fit = False

    # Get header
    mn_header = Mn['header']
    # Get data
    mn_data = Mn['data']

    # Get units conversion factors
    uts = Mn['units']

    # Create figure
    fig = figure(figsize=figsize)
    wcs = WCS(mn_header)
    ax = fig.add_subplot(111, projection=wcs)

    # Get dimensions of M0 map 
    x_mn, y_mn = Mn['data'].shape
    # Get size of M0 bin
    x_mn_bin, y_mn_bin = x_mn/xbins, y_mn/ybins

    # Plot M0 map
    img = ax.imshow(Mn['data'].T, origin='lower')
    fig.colorbar(img, orientation='vertical', label=r''+uts)
    
    # Add contour
    if contour:
        ax.contour(Mn['data'].T, colors=color_contours, alpha=alpha_contours, \
                   levels=level_contours)
    
    # Add axis labels
    ctypex = mn_header['CTYPE1']
    ctypey = mn_header['CTYPE2']
    ax.set_xlabel(r''+ctypex)
    ax.set_ylabel(r''+ctypey)

    # Maximum data for all the bins
    if mc_binned:
        spectra_factor = np.nanmax(get_max_spectra_stats(mc))
    else:
        spectra_factor = np.nanmax(mc['line'])

    # Plot spectra lines
    for b in range(nbins):
        # Binned map
        if mc_binned:
            name = 'B'+str(b)
            # Get bin coordinates
            i, j = mc[name]['pos']
            # Get line profile data
            vel  = mc[name]['vel']
            line = mc[name]['line']
            if fit and flag_fit:
                # Line Parameters
                A     = mc[name]['A']
                mu    = mc[name]['mu']
                sigma = mc[name]['sigma']
        # Full map
        else:
            name = "Full"
            # Coordinates
            i, j = 0, 0
            # Get line profile data
            vel  = mc['vel']
            line = mc['line']
            if fit and flag_fit:
                # Line Parameters
                A     = mc['A']
                mu    = mc['mu']
                sigma = mc['sigma']

        # Observed profile
        veln = vel - vel[0]                                 # Align to zero
        veln = i*x_mn_bin + veln*int(x_mn_bin)/veln[-1]           # Resize velocity vector
        line = j*y_mn_bin + line*y_mn_bin/(1.08*spectra_factor)   # Resize line

        # Plot data line
        ax.plot(veln, line, color=line_color, **kwargs)

        # Fitted profile
        if fit and flag_fit:
            popt = _get_line_params_as_list(A, mu, sigma)
            line_fit = j*y_mn_bin + multiple_lines(vel, ['G']*len(A), *popt)*y_mn_bin/(1.08*spectra_factor)
            ax.plot(veln, line_fit, color=fit_color, **kwargs)
        # Show bin label
        if label:
            ax.text(veln[0], (j+0.75)*y_mn_bin, name, color=text_color)

    # SOS logo
    _add_logo(fig)


def plot_line(mc, fit=True, nbin=None, **kwargs):
    """
        Plot line velocity profile
        Parameters
        ----------
        mc: dict
            Molecular cloud from a extracted molecule.
        fit : boolean (optional)
            Plot the fitted profile?
        nbin : string (optional)
            Number of bin. If the molecular cloud dictionary 
            does not contain bins, it just ignore it.
        ----------
    """

    # Key arguments
    # ----------------------------------------------
    # Figure size
    figsize = kwargs.pop('figsize', (9,9))
    # ----------------------------------------------

    # Get binned or full dimensions 
    dims, nbins, mc_binned = _get_mc_dims(mc)

    # Check if fitting data is available
    flag_fit = True
    if fit:
        if mc_binned:
            if not 'A' in mc['B0'].keys():
                msg('Fitting data not available', 'warn')
                flag_fit = False
        else:
            if not 'A' in mc.keys():
                msg('Fitting data not available', 'warn')
                flag_fit = False

    # Get conversion factors
    vf = get_factors_units('km s^-1', su)[0]

    # Check nbins
    if nbin:
        if not nbin in mc.keys():
            msg('Bin '+nbin+' data not available', 'fail')
            return
        # Line profile
        vel  = mc[nbin]['vel']
        line = mc[nbin]['line']
        # Get fitting params
        if fit and flag_fit:
            A     = mc[nbin]['A']
            mu    = mc[nbin]['mu']
            sigma = mc[nbin]['sigma']

    else:
        # Line profile
        vel  = mc['vel']
        line = mc['line']
        # Get fitting params
        if fit and flag_fit:
            A     = mc['A']
            mu    = mc['mu']
            sigma = mc['sigma']

    # Scale
    vel = np.array(vel)/vf
    mu = np.array(mu)/vf
    sigma = np.array(sigma)/vf

    # Define figure
    fig = figure(figsize=figsize)
    ax = fig.add_subplot(111)
   
    # Get axis labels
    ax.set_xlabel(r'Velocity [km/s]')
    ax.set_ylabel(r'Temperture [K]')
    
    # Plot data
    ax.plot(vel, line, label=r'Measured', **kwargs)
    # Fit data
    if fit and flag_fit:    
        popt = _get_line_params_as_list(A, mu, sigma)
        ax.plot(vel, multiple_lines(vel, ['G']*len(A), *popt), 'ro:', label=r'Fitted')
    # Show legend
    ax.legend()

    # SOS logo
    _add_logo(fig)


def map_param(mc, param, Mn, log=True, **kwargs):
    """
        Plot line velocity profile
        Parameters
        ----------
        mc: dict
            Molecular cloud from a extracted molecule.
        param : string
            Physical Parameter to plot
        Mn : dict
            n-Momentum map as a dictionary with two keys: 
            header and data. Contours will be plotted
        log : boolean (optional)
            Log scale?
        ----------
    """

    # Key arguments
    # ----------------------------------------------
    # Figure size
    figsize = kwargs.pop('figsize', (9,9))
    # Level contours
    level_contours = kwargs.pop('level_contours', 5)
    # Level contours
    color_contours = kwargs.pop('color_contours', 'black')
    # Alpha contours
    alpha_contours = kwargs.pop('alpha_contours', 0.5)
    # Color scale
    cmap = kwargs.pop('cmap', 'gist_gray_r')
    # Log contour
    log_contour = kwargs.pop('log_contour', False)
    # ----------------------------------------------

    # Get Parameters
    m_param, dims_bins = _binparams2img(mc, param)
    
    # If there are parameters available
    if len(m_param) > 0:
        # If log scale
        if log:
            m_param = np.log10(m_param)

        # Create figure
        fig = figure(figsize=figsize)

        # Get header dimensions
        dims_mn = Mn['data'].shape
        mn_header = _header_binned(Mn['header'], dims_mn, dims_bins)
        
        # Get world coordinates
        wcs = WCS(mn_header)
        ax = fig.add_subplot(111, projection=wcs)
        
        # Show axis labels
        ax.set_title(param)
        ctypex = mn_header['CTYPE1']
        ctypey = mn_header['CTYPE2']
        ax.set_xlabel(r''+ctypex)
        ax.set_ylabel(r''+ctypey)

        # Get label
        # ----------------------------------------------
        # Logaritmic scale
        label = r''
        if log:
            label += r'log$_{10}$'
        # Parameter
        pm = param.split('_')[0]

        if param in sos.label_params.keys():
            label += sos.label_params[param]
        elif pm in sos.label_params.keys():
            label += sos.label_params[pm]
        elif param[0] in sos.label_params.keys():
            label += sos.label_params[param[0]]
        else:
            label += 'Undefined'
        # ----------------------------------------------

        # Plot parameter       
        img = ax.imshow(m_param.T, cmap=cmap, origin='lower')
        fig.colorbar(img, orientation='vertical', label=label)

        # Plot M0 map
        ins = ax.inset_axes([0.0,0.0,1.0,1.0])
        ins.patch.set_alpha(0.0)
        ins.axis('off')
        ins.set_xticks([])
        ins.set_yticks([])
        if log_contour:
            data = np.log10(Mn['data'])
        else:
            data = Mn['data']
        ins.contour(data.T, colors=color_contours, alpha=alpha_contours,\
                    levels=level_contours)

        # SOS logo
        _add_logo(fig)


def plot_moment(Mn, contour=True, **kwargs):
    """
        Plot Momentum zero map
        Parameters
        ----------
        Mn : dict
            Momentum Zero map as a dictionary with two keys: 
            header and data. Contours will be plotted
        contour : boolean (optional)
            Show contour?
        ----------
    """

    # Key arguments
    # ----------------------------------------------
    # Figure size
    figsize = kwargs.pop('figsize', (9,9))
    # Level contours
    level_contours = kwargs.pop('level_contours', 5)
    # Level contours
    color_contours = kwargs.pop('color_contours', 'black')
    # Alpha contours
    alpha_contours = kwargs.pop('alpha_contours', 0.5)
    # Color scale
    cmap = kwargs.pop('cmap', 'gist_gray_r')
    # Log scale
    log = kwargs.pop('log', False)
    # ----------------------------------------------

    # Get order
    uts = Mn['units']

    # Create figure
    fig = figure(figsize=figsize)
    mn_header = Mn['header']
    wcs = WCS(mn_header)
    ax = fig.add_subplot(111, projection=wcs)

    # Show M0 map
    data_to_plot = Mn['data'].T
    # Log scale?
    if log:
        data_to_plot = np.log10(data_to_plot)
    img = ax.imshow(data_to_plot, cmap=cmap, origin='lower')
    fig.colorbar(img, orientation='vertical', label=uts)
    
    # Add contours
    if contour:
        ax.contour(Mn['data'].T, colors=color_contours, alpha=alpha_contours,\
                   levels=level_contours)

    # Add axis labels
    ctypex = mn_header['CTYPE1']
    ctypey = mn_header['CTYPE2']
    ax.set_xlabel(r''+ctypex)
    ax.set_ylabel(r''+ctypey)

    # SOS logo
    _add_logo(fig)


def plot_pol_vectors(pol_vectors, map_param, step=1, save=True, **kwargs):
    """
        Plot Polarization vectors
        Parameters
        ----------
        pol_vectors : np.array complex
            Polarization vectors
        spatial_coords : array
            Longitude and Latitude spatial coordinates
        map_param : dict
            Parameter map (N, M0, n). It should contains
            data and header keywords
        step : int (optional)
            Vector step
        save : boolean (optional)
            Save figure?
        ----------
    """

    # Key arguments
    # ----------------------------------------------
    # Figure size
    figsize = kwargs.pop('figsize', (9,9))
    # Level contours
    level_contours = kwargs.pop('level_contours', 5)
    # Level contours
    color_contours = kwargs.pop('color_contours', 'black')
    # Color vectors
    color_vectors = kwargs.pop('color_vectors', 'red')
    # Log-scale
    log = kwargs.pop('log', False)
    # Alpha contours
    alpha_contours = kwargs.pop('alpha_contours', 0.5)
    # Color scale
    cmap = kwargs.pop('cmap', 'gist_gray_r')
    # Name figure
    name = kwargs.pop('name', 'vectors')
    # Rotate vector
    rot = kwargs.pop('rot', 0.)
    # Label
    label = kwargs.pop('label', '---')
    # ----------------------------------------------

    if isinstance(pol_vectors, dict):
        
        if not 'B0' in pol_vectors.keys():
            msg('Bin B0 does not exist', 'fail')
            return

        pol_key = ''
        for i in pol_vectors['B0'].keys():
            if isinstance(pol_vectors['B0'][i], complex):
                pol_key = i
                break

        if pol_key == '':
            msg('Polarization vectors not found', 'fail')
            return

        bins = pol_vectors.keys()
        last_bit = 'B' + str(np.max([int(i[1:]) for i in bins]))
        sx, sy = pol_vectors[last_bit]['pos']
        pol_data = np.zeros((sx+1, sy+1), dtype=complex)
        for b in bins:
            xp, yp = pol_vectors[b]['pos']
            pol_data[xp, yp] = pol_vectors[b][pol_key]

        # Full dimensions comes from the zero map
        dims = map_param['data'].shape
        lon = np.linspace(0, dims[0], sx+1, endpoint=False)
        lon += (lon[1]-lon[0])/2.
        lat = np.linspace(0, dims[1], sy+1, endpoint=False)
        lat += (lat[1]-lat[0])/2.

        # Lon/lat dimensions
        x_vecs = np.arange(0, sx+1, step)
        y_vecs = np.arange(0, sy+1, step)

    else:
        if len(pol_vectors) == 0:
            raise Exception("Polarization vectors array is not defined yet.")

        pol_data = pol_vectors
        # Get data according to the step
        dims_vec = pol_data.shape
        dims_map = map_param['data'].shape
        off_px = -0.5
        lon = np.arange(off_px, dims_map[0]+off_px, step*dims_map[0]/dims_vec[0])
        lat = np.arange(off_px, dims_map[1]+off_px, step*dims_map[1]/dims_vec[1])

        # Lon/lat dimensions
        x_vecs = np.arange(0, dims_vec[0], step)
        y_vecs = np.arange(0, dims_vec[1], step)

    # Check if map param is a binned data
    if isinstance(map_param['header'], list):
        map_param['header'] = _header_binned(*map_param['header'])

    # Create figure
    fig = figure(figsize=figsize)
    m0_header = map_param['header']
    wcs = WCS(m0_header)
    ax = fig.add_subplot(111, projection=wcs)

    # Display Map param
    m0_data = map_param['data']
    if log:
        m0_data = np.log10(m0_data)

    img = ax.imshow(m0_data.T, origin='lower', cmap=cmap)
    fig.colorbar(img, orientation='vertical', label=label)
    if contour:
        ax.contour(m0_data.T, colors=color_contours, alpha=alpha_contours,\
                   levels=level_contours)

    # Add axis labels
    ctypex = m0_header['CTYPE1']
    ctypey = m0_header['CTYPE2']
    ax.set_xlabel(r''+ctypex)
    ax.set_ylabel(r''+ctypey)

    pol_vec_step = np.zeros((len(lon), len(lat)), dtype=complex)
    # Get stepped vectors
    for m, i in enumerate(x_vecs):
        for n, j in enumerate(y_vecs):
            pol_vec_step[m,n] = pol_data[i,j]

    if rot > 0:
        xv = pol_vec_step.real*np.cos(rot) - pol_vec_step.imag*np.sin(rot)
        yv = pol_vec_step.real*np.sin(rot) + pol_vec_step.imag*np.cos(rot)
    else:
        xv = pol_vec_step.real
        yv = pol_vec_step.imag

    # Plot vectors
    ax.quiver(lon, lat, xv.T, yv.T, \
              color=color_vectors, headaxislength=0, headlength=0, **kwargs)
    ax.quiver(lon, lat, -1*xv.T, -1*yv.T, \
              color=color_vectors, headaxislength=0, headlength=0, **kwargs)

    # SOS logo
    _add_logo(fig)

    # Save figure
    if save:
        fig.savefig(name+'_pol_map.png')
        msg('Figure saved', 'ok')


def _get_line_params_as_list(A, mu, sigma):
    """
        Convert dictionary to list fitting parameters
        Parameters
        ----------
        A : float
            Amplitude
        mu : float
            Mean
        sigma : float
            Dispersion
        ----------
    """
    peaks = len(A)
    popt = np.zeros(3*peaks)
    for i in range(peaks):
        popt[3*i] = A[i]
        popt[3*i+1] = mu[i]
        popt[3*i+2] = sigma[i]

    return popt


def _get_mc_dims(mc):
    """
        Maximum data of all the bins
        Parameters
        ----------
        mc : dict
            Molecular cloud dimensions
        ----------
    """
    # Binned molecular cloud
    if 'B0' in mc.keys():
        mc_binned = True
        nbins = len(mc.keys())
        # Find the binned dimensions
        xmax, ymax = 0, 0
        for b in range(nbins):
            bin_name = 'B'+str(b)
            i, j = mc[bin_name]['pos']
            if i >= xmax:
                xmax = i 
            if j >= ymax:
                ymax = j
        dims = [xmax+1, ymax+1]
    # Full molecular cloud
    else:
        mc_binned = False
        nbins = 1
        dims = [1, 1]

    return dims, nbins, mc_binned


def _binparams2img(mc, param):
    """
        Maximum data of all the bins
        Parameters
        ----------
        mc : dict
            Molecular cloud dimensions
        param : boolean
            Parameter
        ----------
    """
    if not param in sos.all_params:
        raise Exception('Parameter not valid')

    # Get binned or full dimensions 
    dims, nbins, mc_binned = _get_mc_dims(mc)
    sx, sy = dims

    # Define paremeter matrix
    param_matrix = np.zeros((sx, sy))

    # Scan all the bins
    for b in range(nbins):
        if mc_binned:
            # Get bin name
            name = 'B'+str(b)
            # Get coordinates
            i, j = mc[name]['pos']
            if not mc[name]['flag']:
                # Get parameter value
                m = mc[name][param]
            else:
                m = np.nan
        else:
            # Get coordinate
            i, j = 0, 0
            # Get parameter value
            m = mc[param]

        #if param == 'den':
        #    m = m/(1.28*1e3*1.672622e-27)     # Neutral gas * proton mass [g]

        param_matrix[i][j] = m

    return param_matrix, dims


def _header_binned(m0_header, dims_m0, dims_bin):
    """
        Maximum data of all the bins
        Parameters
        ----------
        m0_header : header
            Momentum zero header
        dims_m0 : tuple
            Momentum zero map dimensions
        dims_bin : tuple
            Binned map dimensions
        ----------
    """
    # Bins dimensions
    sx, sy = dims_bin
    # M0 dimensions
    x, y = dims_m0
    # Get bin size
    stepx, stepy = x/sx, y/sy

    # Backup the header
    mv_header = m0_header.copy()
    # Update header steps
    mv_header['CDELT1'] = mv_header['CDELT1']*stepx
    mv_header['CDELT2'] = mv_header['CDELT2']*stepy
    # Update home pixels
    mv_header['CRPIX1'] = 0.5 + (mv_header['CRPIX1']-0.5)/stepx
    mv_header['CRPIX2'] = 0.5 + (mv_header['CRPIX2']-0.5)/stepy

    return mv_header


def _add_logo(fig):
    """
        Add logo to the figures
    """
    newax = fig.add_axes([0.01, 0.85, 0.25, 0.135], anchor='NE', zorder=-1)
    newax.imshow(sos.LOGO_SOS)
    newax.axis('off')