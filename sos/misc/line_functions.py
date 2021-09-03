# -*- coding: utf-8 -*-
# --------------------------------------------------------------------------------- #
# Software de Observaciones Sintéticas S.O.S.
# Line functions
#
# Marcial Becerril, @ 25 August 2020
# Latest Revision: 25 Aug 2020, 21.33 GMT
#
# For all kind of problems, requests of enhancements and bug reports, please
# write to me at:
#
# mbecerrilt92@gmail.com
# mbecerrilt@inaoep.mx
#
# --------------------------------------------------------------------------------- #

import numpy as np

def gaussian(x, A, mu, sigma):
    """
        Gaussian function
        Parameters
        ----------
        x : int/float/array
        A : float
            Amplitude
        mu : float
            Mean
        sigma : float
            Dispersion
        y0 : float
            Offset
        ----------
    """
    return A*np.exp(-((x-mu)**2)/(2*sigma**2))

def lorentzian(x, A, mu, w):
    """
        Gaussian function
        Parameters
        ----------
        x : int/float/array
        A : float
            Amplitude
        mu : float
            Mean
        w : float
            Width
        y0 : float
            Offset
        ----------
    """
    w = np.abs(w)
    return A*(w/(4*(x-mu)**2 + w**2))
