# -*- coding: utf-8 -*-
# --------------------------------------------------------------------------------- #
# Software de Observaciones Sintéticas S.O.S.
# Units handling
#
# Marcial Becerril, @ 5 May 2021
# Latest Revision: 5 May 2021, 14:00 GMT-6
#
# For all kind of problems, requests of enhancements and bug reports, please
# write to me at:
#
# mbecerrilt92@gmail.com
# mbecerrilt@inaoep.mx
#
# --------------------------------------------------------------------------------- #

import numpy as np

# Fundamental magnitudes for sos
# MKS system [REFERENCE SYSTEM]
mks = {
    'mass'    : {'name': 'kg',
                'units': {'g':1000, 'lb':2.20462, 'Ms':5.0274e-31}},
    'length'  : {'name': 'm',
                'units': {'m':1,  'ft':3.28084, 'pc':3.24077929e-17}},
    'time'    : {'name': 's',
                'units': {'s':1,  'min':1/60,  'hr':1/3600}},
    'angular' : {'name': 'deg',
                'units': {'deg':1, 'degree':1,'rad':np.pi/180, 'arcmin':60, 'arcsec':3600}},
    'temp'    : {'name': 'K',
                'units': {'k':1, 'kcmb':1}},
    'pixel'   : {'name': 'pixel',
                'units': {'pixel':1, 'px':1}},
    'flux'    : {'name': 'W m^-2 Hz^-1',
                'units': {'w m^-2 hz^-1':1, 'jy':1e26}},
    'freq'    : {'name': 'Hz',
                'units': {'hz':1}},
    'power'   : {'name': 'W',
                'units': {'w':1, 'erg s^-1':1e7}},
    'den_mag' : {'name': 'T',
                'units': {'t':1, 'g':1e4}},
    'energy'  : {'name': 'J', 
                'units': {'j':1, 'erg':1e7}},
    'flux_mag': {'name': 'Wb',
                'units': {'wb':1, 'mx':1e8}}
}
# CGS system
cgs = {
    'mass'    : {'name': 'g',
                'units': {'g':1, 'lb':0.00220462, 'Ms':5.0274e-33}},
    'length'  : {'name': 'cm',
                'units': {'m':0.01,  'ft':0.0328084, 'pc':3.24077929e-19}},
    'time'    : {'name': 's',
                'units': {'s':1,  'min':1/60,  'hr':1/3600}},
    'angular' : {'name': 'deg',
                'units': {'deg':1, 'degree':1, 'rad':np.pi/180, 'arcmin':60, 'arcsec':3600}},
    'temp'    : {'name': 'K',
                'units': {'k':1, 'kcmb':1}},
    'pixel'   : {'name': 'pixel',
                'units': {'pixel':1, 'px':1}},
    'flux'    : {'name': 'erg cm^-2 s^-1 Hz^-1',
                'units': {'erg cm^-2 s^-1 hz^-1':1, 'jy':1e23}},
    'freq'    : {'name': 'Hz',
                'units': {'hz':1}},
    'power'   : {'name': 'erg s^-1',
                'units': {'erg s^-1':1, 'w':1e-7}},
    'den_mag' : {'name': 'G',
                'units': {'g':1, 't':1e-4}},
    'energy'  : {'name': 'erg', 
                'units': {'erg':1, 'j':1e-7}},
    'flux_mag': {'name': 'Mx',
                'units': {'mx':1, 'wb':1e-8}}
}

# Prefix-Subfix
sub_pre = {'a':-18,'f':-15,'p':-12,'n':-9,'u':-6,'m':-3,'c':-2,'d':-1,
           'D':1,'H':2,'k':3,'M':6,'G':9,'T':12,'P':15,'E':18}


# Extract subfix or prefix
def _get_presub_fix(unit, sys):

    # Get all the units in terms of the International System of Units
    for mag in sys.keys():
        for u in sys[mag]['units'].keys():
            if unit == u:
                return (1/sys[mag]['units'][u], mag)
            elif len(unit) > 0:
                if unit[1:] == u:
                    if unit[0] in sub_pre:
                        factor = 10**sub_pre[unit[0]]/sys[mag]['units'][u]
                        return (factor, mag)
                    else:
                        print("Prefix/Subfix not defined")
                        return (None, None)
    print("Magnituide is not available")
    return (None, None)


def get_factors_units(units, ref):
    """
        Get the prefix and subfix of the units
    """
    # Identify the system of units
    if ref == 'mks':
        su_ref = mks
    elif ref == 'cgs':
        su_ref = cgs

    # Factors
    factors = {}
    # Get unities
    uts = _parse_units(units)
    # Initiate pows
    pows = 0
    # Get factor
    f = 1
    # Get the multiples
    for u in uts.keys():
        # Search for the param
        k, mag = _get_presub_fix(u, su_ref)
        # Get factors
        if mag in factors:
            factors[mag]['factor'] = factors[mag]['factor']*(k)**uts[u]
            factors[mag]['pow'] += uts[u]
        else:
            factors[mag] = {}
            factors[mag]['factor'] = k**uts[u]
            factors[mag]['pow'] = uts[u]

        # Get factor
        f = f*factors[mag]['factor']
        # Get units
        factors[mag]['units'] = su_ref[mag]['name']

    return f, factors


def unit2str(units):
    """
        Transform units dictionary format to dictionary
    """
    str_units = ""
    for u in units.keys():
        p = units[u]['pow']
        ut = units[u]['units']
        if p == 1:
            str_units += ut+' '
        else:
            str_units += ut+'^'+str(p)+' '

    return str_units[:-1]


def _parse_units(units):
    """
        Get the units factors regarding the cgs system
    """
    # Get all the physical magnitudes
    mags = []
    pwrs = []

    # Sign of the power
    sign_pwr = 1
    val_pwr = 1

    # Multiplications symbols
    mul_signs = [' ', '*']
    # Auxiliar units
    aux_unit = ''
    aux_pwr = ''
    # Some flags
    flag_apply = False
    flag_div = False
    flag_mag = False
    flag_pwr = False
    for u in units:
        if u.isalpha():
            aux_unit += u
            flag_mag = True
        elif u == '^':
            flag_pwr = True 
        elif flag_pwr and (u.isnumeric() or u in [',','.','-']):
            if u == '-' and len(aux_pwr) > 0:
                if aux_pwr[0] == '-':
                    aux_pwr = aux_pwr[1:]
            else:
                aux_pwr += u
        elif u == '/':
            flag_div = True
            flag_apply = True
        else:
            flag_apply = True
        
        # If the character is the last
        if u == units[-1]:
            flag_apply = True

        if flag_apply:
            if flag_pwr or flag_mag:
                if aux_pwr == '':
                    aux_pwr = '1'
                elif aux_pwr == '-':
                    aux_pwr = '-1'
                pwrs.append(int(aux_pwr))
                flag_pwr = False
                if flag_div:
                    aux_pwr = '-'
                    flag_div = False
                else:
                    aux_pwr = ''
            if flag_mag:
                mags.append(aux_unit)
                aux_unit = ''
                flag_mag = False
            flag_apply = False

    # Check if the units are repeated
    uts = {}
    while len(mags) > 0:
        # Masking
        mask = [a==mags[0] for a in mags]
        # Get powers
        pwrs_mag = [pwrs[n] for n in range(len(mask)) if mask[n]]
        # Asign not repeated magnitude and powers
        uts[mags[0]] = sum(pwrs_mag)
        # Update magnitudes and powers
        mags = [mags[m] for m in range(len(mask)) if not mask[m]]
        pwrs = [pwrs[m] for m in range(len(mask)) if not mask[m]]

    return uts


# su = 'cgs'
# units = "W m^-2 Hz^-1"
# factors, uts = _get_factors_units(units, su)