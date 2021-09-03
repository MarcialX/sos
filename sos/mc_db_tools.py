# -*- coding: utf-8 -*-
# --------------------------------------------------------------------------------- #
# Software de Observaciones Sintéticas" S.O.S.
# Tools for the handling of data base of molecular clouds and thier parameters
#
# Marcial Becerril, @ 25 August 2020
# Latest Revision: 25 Aug 2020, 22.33 GMT
#
# For all kind of problems, requests of enhancements and bug reports, please
# write to me at:
#
# mbecerrilt92@gmail.com
# mbecerrilt@inaoep.mx
#
# --------------------------------------------------------------------------------- #

import yaml

from .misc.print_msg import *

class mc_db(object):

    def __init__(self, path):
        """
            Load the molecular cloud data base
            Parameters
            ----------
            path : string
                File path of the data base
            ----------
        """
        # Starting variables
        self._dist = 0.0
        self._ang_dia = 0.0
        self._mc_factor = 1.0
        self._Re_factor = 1.0
        self._inst_factor = 1.0
        self._M13CO = ''
        self._M12CO = ''
        self._MC18O = ''
        self._x0, self._y0 = 0.0, 0.0

        self.mc_db = {}
        self.mc_id = ''
        self.path = path

        # Loading data base
        msg('Loading data base', 'info')
        try:
            with open(path) as file:
                mc_db = yaml.safe_load(file)
            msg('Data base loaded', 'ok')
            self.mc_db = mc_db
        except Exception as e:
            msg('Fail loading data base. '+str(e), 'fail')
            return


    def get_mc_params(self, mcID):
        """
            Get molecular cloud parameters
            Parameters
            ----------
            mcID : string
                Molecular cloud ID
            ----------
        """
        self.mc_id = mcID
        mc_params = self.mc_db[mcID]

        self._dist = mc_params['dist']
        self._ang_dia = mc_params['ang_dia']
        self._width_factor = mc_params['width_factor']
        self._Re_factor = mc_params['Re_factor']
        self._temp_factor = mc_params['temp_factor']
        self._M13CO = mc_params['13CO']
        self._M12CO = mc_params['12CO']
        self._MC18O = mc_params['C18O']
        self._x0, self._y0 = mc_params['x0'], mc_params['y0']

        return mc_params


    def add_mc(self, name, params):
        """
            Add molecular cloud to the data base
            Parameters
            ----------
            name : string
                Molecular cloud ID
            params : list
                Array with the parameters of the molecular cloud
                    dist: float
                    ang_dia: float
                    width_factor: float
                    Re_factor: float
                    temp_factor: float
                    mc_center: tuple
                    13CO: string
                    12CO: string
                    C18O: string
            ----------
        """
        if len(params) == 9:
            params_dict = {'dist'   :   params[0],
                           'ang_dia':   params[1],
                           'width_factor':   params[2],
                           'Re_factor':   params[3],
                           'temp_factor':   params[4],
                           'x0': params[5],
                           'y0': params[6], 
                           '13CO':   params[7],
                           '12CO':   params[8],
                           'C18O':   params[9]
            }
            self.mc_db[name] = params_dict
            msg(name+' added!', 'ok')
        else:
            msg('Adding molecular cloud. Parameters size wrong', 'fail')


    def save_mc_db(self):
        """
            Save molecular cloud data base
        """
        try:
            with open(self.path, 'w') as file:

                self.mc_db[self.mc_id]['dist'] = self._dist
                self.mc_db[self.mc_id]['ang_dia'] = self._ang_dia
                self.mc_db[self.mc_id]['width_factor'] = self._width_factor
                self.mc_db[self.mc_id]['Re_factor'] = self._Re_factor
                self.mc_db[self.mc_id]['temp_factor'] = self._temp_factor
                self.mc_db[self.mc_id]['x0'] = self._x0
                self.mc_db[self.mc_id]['y0'] = self._y0
                self.mc_db[self.mc_id]['13CO'] = self._M13CO
                self.mc_db[self.mc_id]['12CO'] = self._M12CO
                self.mc_db[self.mc_id]['C18O'] = self._MC18O

                documents = yaml.dump(self.mc_db, file)

            msg('Molecular Cloud data base saved', 'ok')
        except:
            msg('Saving data base', 'fail')


    # Distance getter function
    @property
    def dist(self):
        return self._dist

    # Distance setter function
    @dist.setter
    def dist(self, dist):
        msg('Distance to '+self.mc_id+': '+str(dist), 'ok')
        self._dist = dist

    # Angular diameter getter function
    @property
    def ang_dia(self):
        return self._ang_dia

    # Angular diameter setter function
    @ang_dia.setter
    def ang_dia(self, ang):
        msg('Angular diameter to '+self.mc_id+': '+str(ang), 'ok')
        self._ang_dia = ang

    # MC factor getter function
    @property
    def mc_factor(self):
        return self._mc_factor

    # MC factor setter function
    @mc_factor.setter
    def mc_factor(self, mc_factor):
        msg(self.mc_id+' factor: '+str(mc_factor), 'ok')
        self._mc_factor = mc_factor

    # Effective radius getter function
    @property
    def Re_factor(self):
        return self._Re_factor

    # Angular diameter setter function
    @Re_factor.setter
    def Re_factor(self, Re_factor):
        msg(self.mc_id+' effective radius: '+str(Re_factor), 'ok')
        self._Re_factor = Re_factor

    # Telescope Coupling factor getter function
    @property
    def inst_factor(self):
        return self._inst_factor

    # Telescope Coupling factor setter function
    @inst_factor.setter
    def inst_factor(self, inst_factor):
        msg(self.mc_id+' copuling factor: '+str(inst_factor), 'ok')
        self._inst_factor = inst_factor

    # Map 13C0 path getter function
    @property
    def path_13CO(self):
        return self._M13CO

    # Map 13C0 factor setter function
    @path_13CO.setter
    def path_13CO(self, path):
        msg(self.mc_id+' 13CO file:'+str(path), 'ok')
        self._M13CO = path

    # Map 12C0 path getter function
    @property
    def path_12CO(self):
        return self._M12CO

    # Map 12C0 factor setter function
    @path_12CO.setter
    def path_12CO(self, path):
        msg(self.mc_id+' 12CO file:'+str(path), 'ok')
        self._M12CO = path

    # Map C18O path getter function
    @property
    def path_C18O(self):
        return self._MC18O

    # Map C18O factor setter function
    @path_C18O.setter
    def path_C18O(self, path):
        msg(self.mc_id+' C18O file:'+str(path), 'ok')
        self._MC18O = path
