# -*- coding: utf-8 -*-
# --------------------------------------------------------------------------------- #
# Software de Observaciones Sintéticas S.O.S.
# Usage example of sos
#
# Marcial Becerril, @ 4 September 2020
# Latest Revision: 29 June 2021, 10:37 GMT-6
##
# For all kind of problems, requests of enhancements and bug reports, please
# write to me at:
#
# mbecerrilt92@gmail.com
# mbecerrilt@inaoep.mx
#
# --------------------------------------------------------------------------------- #

from matplotlib.pyplot import *

# Import sos library
import sos

# Create molecular cloud object, defining the ID of the cloud
mc = sos.mc('TaurusMC', move_pol='r-r-t')

# FULL REGION
# ===========
# Get the integrated line profiles of the whole map, for all the molecules
# which data is available
mc.get_map_vels()
# Fit 13CO line profile
mc.line_fit('13CO', forced_lines=1)
# Fit 12CO line profile
mc.line_fit('12CO', forced_lines=1)
# Calculate the physical parameters of the whole cloud
mc.get_gral_params()

# BINNING THE MAP
# ===============
# Binning the map of molecule '13CO' and '12CO', dividing the map in 32x32 bins
mc.binning_mol('13CO', nbins=41, rebin=True)
mc.binning_mol('12CO', 41)
# Fit the line of each bin
mc.line_fit_binning('13CO', forced_lines=1)
mc.line_fit_binning('12CO', forced_lines=1)
# Get physical parameters of each bin
mc.get_bins_params()
# Filter mass errors calculations
mc.param_filter('mass_lte', max=1e4)

# PLOT MASS RESULTS
# =================
# Get the M0 map of each molecule
m0_13co = mc.get_n_moment('13CO', n=0)
# Plot momentum zero
sos.plot_moment_spectra(m0_13co, mc.extract_param_from_binned('13CO'), label=False)
# Plot mass lte as a function of position
sos.map_param(mc.binned, 'mass_lte', m0_13co, cmap='rainbow', log=False)

# Display the general results
mc.summary()

# Polarization analysis
# =====================
# Segment the polarization maps
mc.binning_pol(nbins=32, rebin=False)

# Get Volumetric density using the mass from LTE method
mc.get_vol_den('mass_lte')

# Get Magnetic Field, using the dispersion velocity of the 
# 13CO molecule
mc.get_mag_field('13CO')

# Filter magnetic field over-estimations on the edges
mc.param_filter('B', max=120)

# PLOT MAG RESULTS
# =================
sos.map_param(mc.binned, 'B', m0_13co, cmap='rainbow', log=False)

# Polarization vector (one vector each 4 pixels) and column density[N]
# Get polarization vectors
mc.get_pol_params()
# Plot magnetic intensity with polarization vectors rotated 90 degrees
n_pol = mc.build_data_header('B', 'binned')
sos.plot_pol_vectors(mc.full['pol_vector'], n_pol, step=4, rot=np.pi/2., label=r'$\mu$G', level_contours=3)

# Show all the maps
show()