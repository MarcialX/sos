# -*- coding: utf-8 -*-
# --------------------------------------------------------------------------------- #
# Software de Observaciones Sintéticas S.O.S.
# Methods to get statistics
#
# Marcial Becerril, @ 28 June 2020
# Latest Revision: 28 Jun 2021, 1:25 GMT-6
#
# For all kind of problems, requests of enhancements and bug reports, please
# write to me at:
#
# mbecerrilt92@gmail.com
# mbecerrilt@inaoep.mx
#
# --------------------------------------------------------------------------------- #

import numpy as np


def get_max_spectra_stats(mc, **kwargs):
    """
        Maximum and minimum of spectra lines of the 
        binned data
        Parameters
        ----------
        mol : string
            Molecule
        ----------
    """
    if not 'line' in mc['B0'].keys():
        raise Exception("Molecular cloud format not valid")

    max_lines_vector = np.zeros(len(mc.keys()))
    
    for i, b in enumerate(mc.keys()):
        max_lines_vector[i] = np.max(mc[b]['line'])

    return max_lines_vector

