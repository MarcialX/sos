# -*- coding: utf-8 -*-
# --------------------------------------------------------------------------------- #
# Software de Observaciones Sintéticas S.O.S.
# Physical constants
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

import numpy as np

Ms = 1.9891e30      # Solar mass [kg]
h = 6.62607e-34     # Planck constant
c = 299.792458e6    # Light speed [m/s]
K = 1.38064e-23     # Stefan-Boltzmann constant
# CMB temperature
Tcmb = 2.72548

# Molecules parameters
# SOS available molecules
mol_params = {
	'12CO':{
		'1-0': {
			'f': 115271204000.0,
			'A10': 7.43e-8
			},
		'2-1': {
			'f': 230538001000.0,
			'A10': 7.13e-7
			},
		'X':5e5,
		'P':0
	},
	'13CO':{
		'1-0': {
			'f': 110201370000.0,
			'A10': 6.49e-8
			},
		'2-1': {
			'f': 220398714000.0,
			'A10': 6.23e-7
			},
		'X':5e5,
		'P':1
	},
	'C18O':{
		'1-0': {
			'f': 109782182000.0,
			'A10': 6.42e-8
			},
		'2-1': {
			'f': 219560369000.0,
			'A10': 6.16e-7
			},
		'X':5e5,
		'P':2
	},
}