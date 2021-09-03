# -*- coding: utf-8 -*-
# --------------------------------------------------------------------------------- #
# Software de Observaciones Sintéticas S.O.S.
# Methods to estimate physical propierties of the molecular cloud
#
# Marcial Becerril, @ 24 August 2020
# Latest Revision: 24 Aug 2020, 18.17 GMT
#
# For all kind of problems, requests of enhancements and bug reports, please
# write to me at:
#
# mbecerrilt92@gmail.com
# mbecerrilt@inaoep.mx
#
# --------------------------------------------------------------------------------- #

import numpy as np
from scipy.integrate import simps
from scipy import integrate

from .misc.constants import *


def column_density(T_thick, T_thin, v_thick, v_thin, fwhm, A10, X=5e5):
    """
        Column density based on Retes-Romero et al. (2018)

        Parameters
        ----------
        T_thick : float
            Temperature line of thick molecule [K]
        T_thin : float
            Temperature line of thin molecule [K]
        v_thick : float
            Transition frequency of thick molecule [Hz]
        v_thin : float
            Transition frequency of thin molecule [Hz]
        fwhm : float
            Full-Width at Half Maximum [km/s] of the molecular line [km s^-1]
        A10 : float
            Einstein A coeficient     
        X : float
            Abundance relative to the Molecular H2       
        ----------
    """
    # Excitation temperature of 12C0
    m_thick = h*v_thick/K
    Tex = m_thick/(np.log(1 + m_thick/(T_thick + Jv(v_thick, Tcmb))))
    # 13C0 Column density and optical depth
    m_thin = h*v_thin/K
    tau_thin = -np.log( np.abs(1.0 - (T_thin / (m_thin/(np.exp(m_thin/Tex) - 1.0) - Jv(v_thin, Tcmb)))) )
    # Column density
    alpha = 1.6*np.pi*K*(v_thin**2)/(3*h*(c**3)*A10)
    N_thin = alpha * (fwhm) * Tex * ( tau_thin/(1.0 - np.exp(-m_thin/Tex)) )

    return N_thin, Tex, tau_thin


def Jv(v, T):
    """
        Jv parametrization

        Parameters
        ----------
        v : float
            Frequency of transition [Hz]
        T : float
            Temperature line [K]
        ----------
    """
    Jv = h*v/(K*(np.exp(h*v/(K*T))-1))

    return Jv


def mass_lte(params, T_thick, T_thin, v_thick, v_thin, fwhm, A10, X=5e5):
    """
        Mass through Local Termodynamic Equilibrium method based on
        Retes-Romero et al. (2018)

        Parameters
        ----------
        params : list
            list of molecular cloud parameters:
                * Distance to the MC [pc]
                * Angular diameter [deg]
                * Effective radius factor [pc]
                * Width line factor
        T12CO : float
            Temperature line of thick molecule [K]
        T13CO : float
            Temperature line of thin molecule [K]
        v_thick : float
            Transition frequency of thick molecule [Hz]
        v_thin : float
            Transition frequency of thin molecule [Hz]
        fwhm : float
            Full-Width at Half Maximum of the molecular line [km s^-1]
        A10 : float
            Einstein A coeficient 
        X : float
            Abundance relative to the H2
        ----------
    """
    # Parameters of the molecular cloud
    d = params[0]       # Distance to the cloud
    theta = params[1]   # Angular diameter
    kRe = params[2]     # Effective radius factor

    # Column density and Excitation temperature
    N_thin, Tex, tau = column_density(T_thick, T_thin, v_thick, v_thin, fwhm, A10, X)
    # Effective radius
    Re = kRe*d*np.tan(theta*0.5*np.pi/180.)

    if tau >= 0:
        # LTE Mass
        MLTE = 0.065 * X * ((Re/5.0)**2) * ((N_thin/1e17))
    else:
        MLTE = 0
        N_thin = 0
        print("Tau is negative. Thin/Thick temperature relation is too high")

    return MLTE, N_thin, Re


def mass_virial(params, fwhm, X=5e5):
    """
        Mass through virial theorem method based on
        Retes-Romero et al. (2018)

        Parameters
        ----------
        params : list
            list of molecular cloud parameters:
                * Distance to the MC [pc]
                * Angular diameter [deg]
                * Effective radius factor [pc]
                * Width line factor [dimensionless]
        fwhm : float
            Full-Width at Half Maximum of the molecular line [km s^-1]
        X : float
            Abundance relative to the H2
        ----------
    """
    # Parameters
    d = params[0]       # Distance to the cloud
    theta = params[1]   # Angular diameter
    kRe = params[2]     # Effective radius factor

    # Effective radius
    Re = kRe*d*np.tan(theta*0.5*np.pi/180.)
    # Virial Mass
    Mvir = 0.0316 * X * (Re/5.0) * (fwhm/5.0)**2

    return Mvir, Re


def mass_xf(params, vel, map_v, X=5e5):
    """
        Mass through virial theorem method based on
        Retes-Romero et al. (2018)

        Parameters
        ----------
        params : list
            list of molecular cloud parameters:
                * Distance to the MC
                * Angular diameter
                * Effective radius factor
                * Width line factor
        vel : array
            Array with all the velocities available for the molecular cloud
        map_v : array
            Flux magnitudes for each velocity
        X : float
            Abundance relative to the H2
        ----------
    """
    # Parameters
    d = params[0]       # Distance to the cloud
    theta = params[1]   # Angular diameter
    kRe = params[2]     # Effective radius factor

    # Effective radius
    Re = kRe*d*np.tan(theta*0.5*np.pi/180.)       # Re [pc]
    # X Factor Mass
    W_13CO = simps(map_v, vel)                    # Integral over all the Velocities
    MX = 0.032 * X * (Re/5.0)**2 * (W_13CO/5.0)       # Factor used by Retes-Romero, direct communication
    
    return MX, Re


def get_nh2(mol, N):
    """
        Hydrogen molecular column density
        Parameters
        ----------
        N_13CO : float
            13CO column density
        ----------
    """
    # H2 Column density
    #XX = 5.0e5     # Factor 13CO/H2 (Dickman 1978)
    XX = mol_params[mol]['X']
    N_H2 = XX*N    # H2 column density

    return N_H2


def get_pol_params_from_stokes(I, Q, U):
    """
        Get polarization parameters
        Parameters
        ----------
        I : float or np.array
            Stokes I parameter
        Q : float or np.array
            Stokes Q parameter
        U : float or np.array
            Stokes U parameter
        ----------
    """
    # Polarization degree
    pol_degree = 100*np.sqrt(np.power(Q, 2)+np.power(U, 2))/np.array(I)
    # Polarization angle
    pol_angle = 0.5*(np.arctan2(U, Q))

    return pol_degree, pol_angle


def sphere_int(x, y, r):
    """
        Sphere integration
        Parameters
        ----------
        x, y : floats
            Points where the evaluation is evaluated
        r : float
            Sphere radius
        ----------
    """
    V = (1/3)*(2*x*y*np.sqrt(r**2-x**2-y**2)-y*(y**2-3*r**2)*np.arctan(x/np.sqrt(r**2-x**2-y**2))- \
    x*(x**2-3*r**2)*np.arctan(y/np.sqrt(r**2-x**2-y**2))-2*r**3*np.arctan(x*y/(r*np.sqrt(r**2-x**2-y**2))))

    return V


def model_volume(xlims, ylims, r, *args, **kwargs):
    """
        Calculate sphere volume of a square line of sight
        Parameters
        ----------
        x0, y0 : float, float
            Central Position of the square section
        w : Width of the square section
            Frequency [Hz]
        r : float
            Sphere radius [pc]
        ----------
    """
    # Parameters
    npoints = kwargs.pop('npoints', 20)            # N-points in double-integral
    shape = kwargs.pop('shape', 'sphere')          # Type of volume: 'sphere' or 'disk'
    thickness = kwargs.pop('thickness', 0.2)       # Thickness of the disk [pc]

    # Y limits
    xmin, xmax = xlims[0], xlims[1]
    # Z limits
    ymin, ymax = ylims[0], ylims[1]

    # Verify the corners
    pts_corners = [(xmin, ymin), (xmax, ymin), (xmax, ymax), (xmin, ymax)]
    fx_on_corners = _on_corners(pts_corners, r, 0, 0)

    # Function covers the full pixel
    if fx_on_corners == 4:
        if shape == 'sphere':
            V = (sphere_int(xmax, ymax, r)-sphere_int(xmin, ymax, r))-(sphere_int(xmax, ymin, r)-sphere_int(xmin, ymin, r))
        elif shape == 'mc_disk':
            V = (xmax-xmin)*(ymax-ymin)*thickness
        else:
            msg('Shape not defined', 'warn')
            return

    # The function at least covers one corner
    elif fx_on_corners > 0:

        # X and Y axis
        xs = np.linspace(xmin, xmax, npoints)
        ys = np.linspace(ymin, ymax, npoints)
        
        # Vectors definition
        f = np.zeros_like(ys)
        fx = np.zeros_like(xs)

        # Get double integral
        for i, y in enumerate(ys):
            for j, x in enumerate(xs):
                # Surface
                Axy = r**2-y**2-x**2
                # SPHERE
                if shape == 'sphere':
                    if Axy >= 0:
                        fxy = 2*np.sqrt(Axy)
                    else:
                        fxy = 0.0
                # DISK
                elif shape == 'mc_disk':
                    if Axy >= 0:
                        fxy = thickness
                    else:
                        fxy = 0.0
                fx[j] = fxy

            mask_x = fx != 0
            xv = xs[mask_x]
            fv = fx[mask_x]

            if len(xv) > 0:
                fy = simps(fv, xv)
            else:
                fy = 0.0
            f[i] = fy

        mask_y = f != 0
        ys = ys[mask_y]
        f = f[mask_y]

        if len(ys) > 0:
            # Get volume    
            V = simps(f, ys)
        else:
            V = 0.0

    # No function
    else:
        V = 0.0

    return V, fx_on_corners


def spot_area(x, y, r, x0, y0):
    """
        Projected area over the sky plane
        Parameters
        ----------
        x, y : floats
            Points where the surface is defined
        r : float
            Radius. The spot area is defined as a circle
        ----------
    """
    sq_spot = r**2-(x-x0)**2-(y-y0)**2
    if sq_spot >= 0:
        fz = np.sqrt(sq_spot)
    else:
        fz = -1
    
    return fz


def _on_corners(pts, r, x0, y0):
    """
        Check how many corners are on the function
        Parameters
        ----------
        pts : tuple of points
            Pixel corners
        r : float
            Radius of the area
        x0, y0 : floats
            Circle offset
        ----------
    """
    cnt_cnr = 0

    for i in range(4):
        f = spot_area(pts[i][0], pts[i][1], r, x0, y0)
        if f >= 0:
            cnt_cnr += 1
            
    return cnt_cnr


def brightness_temperature(S, v, one_px_area):
    """
        Brightness temperature from flux density per omega
        Parameters
        ----------
        S : float
            Flux density per omega unit [Jy/pix]
        v : float
            Frequency [Hz]
        one_px_area : float
            Area in square degrees of 1 pix
        ----------
    """
    # Wavelength
    l = c/v
    omega = 180**2/(one_px_area*np.pi**2)
    # Brightness temperature
    Tb = S*omega*(l**2)/(2*K)

    return Tb


def magnetic_field(den, dV, dTheta):
    """
        Perpendicular Magnetic Field through Chandrasekhar-Fermi method [uG]
        Parameters
        ----------
        den : float
            Volumetric density [g cm^-3]
        dV : float
            Velocity Dispersion [cm/s]
        dTheta : float
            Polarization angle dispersion [rad]
        ----------
    """
    # Epsilon factor 0.7 a 1 de acuerdo con Cho & Yoo (2016)
    # E = 0.5 Ostriker 2001
    E = 0.5
    # Velocity to cgs unit system
    #dV = 100*dV # [cm/s]
    # Magnetic field
    B = E*((4*np.pi*den)**(0.5))*dV/dTheta

    return B*1e6

