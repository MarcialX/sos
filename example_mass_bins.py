# -*- coding: utf-8 -*-
# --------------------------------------------------------------------------------- #
# Software de Observaciones Sintéticas" S.O.S.
# Usage example of sos. Getting masses different bins
#
# Marcial Becerril, @ 4 September 2020
# Latest Revision: 15 June 2021, 12:47 GMT-6
##
# For all kind of problems, requests of enhancements and bug reports, please
# write to me at:
#
# mbecerrilt92@gmail.com
# mbecerrilt@inaoep.mx
#
# --------------------------------------------------------------------------------- #

# Import sos library
import sos
import numpy as np
from matplotlib.pyplot import *

# Create molecular cloud object, defining the ID of the cloud
mc = sos.mc('TaurusMC')
# Get the integrated line profiles of the whole map, for all the molecules
# which data is available
mc.get_map_vels()

# FULL REGION
# ===========
# Fit 13CO line profile
mc.line_fit('13CO')
# Fit 12CO line profile
mc.line_fit('12CO')
# Calculate the physical parameters of the whole cloud
mc.get_gral_params()

mass_lte_one_bin = mc.full['mass_lte']
mass_vir_one_bin = mc.full['mass_vir']
mass_xf_one_bin = mc.full['mass_xf']

# BINNING THE MAP
# ===============
bins = np.arange(7)
mass_lte = []
mass_vir = []
mass_xf = []
manual_13CO = [[],[],[],[],[],[],[]]
manual_12CO = [[],[],[],['B36'],[],[],['B1457','B1586','B1839','B1903','B1904','B1967', 'B3186']]
for n in bins:
    print('BIN: ', n)
    # Binning the map of molecule '13CO' and '12CO', dividing the map in 16x16 bins
    mc.binning_mol('13CO', 2**n, rebin=True)
    mc.binning_mol('12CO', 2**n)
    # Fit the line of each bin
    mc.line_fit_binning('13CO', forced_lines=1)
    mc.line_fit_binning('12CO', forced_lines=1)
    # Fitting bad analysis
    if len(manual_12CO[n]) > 0:
        for i in range(len(manual_12CO[n])):
            ioff()
            mc.line_fit_binning('12CO', nbin=manual_12CO[n][i], inter=True)
            ion()
    if len(manual_13CO[n]) > 0:
        for i in range(len(manual_13CO[n])):
            ioff()
            mc.line_fit_binning('13CO', nbin=manual_13CO[n][i], inter=True)
            ion()
    # Get physical parameters of each bin
    mc.get_bins_params()

    mass_lte.append(mc.sum_binning_param('mass_lte'))
    mass_vir.append(mc.sum_binning_param('mass_vir'))
    mass_xf.append(mc.sum_binning_param('mass_xf'))

# PLOT RESULTS
# ============
figure()
plot(bins, mass_lte, 'rs-', label=r'LTE')
plot(bins, mass_vir, 'bs-', label=r'Virial')
plot(bins, mass_xf, 'gs-', label=r'XF')

axhline(mass_lte_one_bin, color='r', lw=1)
axhline(mass_vir_one_bin, color='b', lw=1)
axhline(mass_xf_one_bin, color='g', lw=1)

xlabel(r'log$_2$($\#$ Number of bins$^{1/2}$)')
ylabel(r'Mass M$_{\odot}$')
legend()
