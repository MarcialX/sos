# -*- coding: utf-8 -*-
# --------------------------------------------------------------------------------- #
# Software de Observaciones Sintéticas S.O.S.
# Define empty molecular cloud and apply polarization analysis
#
# Marcial Becerril, @ 02 June 2021
# Latest Revision: 29 June 2021, 10:30 GMT-6
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

# Map dimensions
dims = (50, 256, 256)
# Velocity limits
vel_lims = (-5e3, 5e3)
# Longitude dims
lon_lims = (-1, 1)
# Latitude dims
lat_lims = (-1, 1)

# Create empty molecular cloud
mc = sos.mc('', dims=dims, vel_lims=vel_lims, lon_lims=lon_lims, lat_lims=lat_lims)

# BINNING THE MAP
# ===============
# Number of bins
bins = 16
# Binning the map of molecule '13CO' and '12CO', dividing the map in 32x32 bins
mc.binning_mol('13CO', nbins=bins, rebin=True)
mc.binning_mol('12CO', bins)

# Define Temperature and line parameters
# 13CO molecule
defined_13co = {
				'T': 10,
				'sigma': 250,
				'mu': 0,
				'A': 0.75
}
mc.fixing_params('13CO', vals=list(defined_13co.values()), params=list(defined_13co.keys()))

# 12CO molecule parameters
defined_12co = {
				'T': 15,
				'sigma': 350,
				'mu': 0,
				'A': 1.
}
mc.fixing_params('12CO', vals=list(defined_12co.values()), params=list(defined_12co.keys()))

# Get physical parameters of each bin
mc.get_bins_params()

# Polarization analysis
# =====================
# Load polarization map in the following path
#path_pol = '/home/marcial/Documentos/sos/data_cube/new_polaris_detector_nr0003.fits'
path_pol = '/home/marcial/Descargas/polaris_detector_nr0003.fits'

# As the image is rotated, we rotate 180° the image with 2 rotations of 90° (r) and then
# transpose the array (t)
mc.load_polarization(path_pol, rotate=90, move_cmd='r-r-t')
# Get polarization vectors
mc.get_pol_params()
# Apply binning
mc.binning_pol(nbins=bins, rebin=False)

# Get density map
for bin in mc.binned.keys():
	mc.binned[bin]['B'] = 1/mc.binned[bin]['pol_angle_std']

N = mc.build_data_header('N', 'full')

# Polarization vector (one vector each 4 pixels) and column density[N]
sos.plot_pol_vectors(mc.full['pol_vector'], N, step=4, log=True, level_contours=3)
# Polarization vectors averaged
sos.plot_pol_vectors(mc.extract_param_from_binned('pol_vector_avg'), N, step=1, log=True, level_contours=3)

# Show mean polarization map degree [averaging all the pixels along the bin]
sos.map_param(mc.binned, 'pol_degree_mean',  N, cmap='rainbow', log=False, log_contour=True, level_contours=3)
# Show standard deviation of the polarization angle
sos.map_param(mc.binned, 'pol_angle_std',  N, cmap='rainbow', log=False, log_contour=True, level_contours=3)

show()

