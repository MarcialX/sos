# -*- coding: utf-8 -*-
# --------------------------------------------------------------------------------- #
# Software de Observaciones Sintéticas S.O.S.
# Some initial parameters
#
# Marcial Becerril, @ 2 Jun 2021
# Latest Revision: 2 Jun 2021, 14:00 GMT-6
#
# For all kind of problems, requests of enhancements and bug reports, please
# write to me at:
#
# mbecerrilt92@gmail.com
# mbecerrilt@inaoep.mx
#
# --------------------------------------------------------------------------------- #

from astropy.io.fits import Header


# MC Parameters
all_params = ['mass_lte', 'mass_vir', 'mass_xf', 'N13CO', 'NH2', 'Re', 'Tex', 'vol', 'B',
 'den', 'pol_angle_mean', 'pol_angle_std', 'pol_degree_std', 'pol_degree_mean',
 '13CO','12CO', 'C18O', 'tau']

# Labels units per parameter
label_params = {'mass':r'M$_{\odot}$', 'N':r'cm$^{-2}$', 'Re':r'Mpc',
				'T':r'K', 'vol':r'pc$^3$', 'den':r'cm$^{-3}$',
				'pol_angle_std':r'rad', 'pol_angle_mean':r'rad',
				'pol_degree_std':r'%', 'pol_degree_mean':r'%',
				'B':r'$\mu$G', 'tau':r'$\tau$' }

# Default unit system
su = 'mks'

# Default Header
hdr_default = Header.fromstring("""\
SIMPLE  =                    T / conforms to FITS standard
BITPIX  =                    8 / array data type
NAXIS   =                    0 / number of array dimensions
EXTEND  =                    T 
NAXIS1  =                   10
CTYPE1  =            'RA---SIN'
CDELT1  =                    1
CUNIT1  =                 'deg'
CRPIX1  =                  0.5
CRVAL1  =                  0.0
NAXIS2  =                   10
CTYPE2  =            'DEC--SIN'
CDELT2  =                    1
CUNIT2  =                 'deg'
CRPIX2  =                  0.5
CRVAL2  =                  0.0
NAXIS3  =                   10	
CTYPE3  =            'VELO-LSR'
CDELT3  =                    1
CUNIT3  =              'm s^-1'
CRPIX3  =                  0.5
CRVAL3  =                  0.0
BUNIT   =                  'K'
BSCALE  =                  1.0
RESTFREQ=           1000000000

""", sep='\n')