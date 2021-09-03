# -*- coding: utf-8 -*-
# --------------------------------------------------------------------------------- #
# Software de Observaciones Sintéticas S.O.S.
# This is the set of tools and functionalities to analyse the propierties of
# molecular clouds.
#
# SPA group at INAOE, @ 24 August 2020
# Latest Revision: 24 Aug 2020, 17:39 GMT
#
# For all kind of problems, requests of enhancements and bug reports, please
# write to me at:
#
# mbecerrilt92@gmail.com
# mbecerrilt@inaoep.mx
#
# --------------------------------------------------------------------------------- #

def help():
    print (
    """
    Welcome to the Software de Observaciones Sintéticas (SOSPY) for molecular clouds analysis.
    For help use sos.help() at anytime to print available top-level options.
    Alternatively, use pcp.<SUBMODULE>.help() to get detailed help for a given <SUBMODULE>.
    """
    )

from .misc import *

from .misc.print_msg import *
#from .misc.units_tool import *

from .mc import *
from .mc_plotter import *
from .mc_db_tools import *
from .misc.constants import *
from .init_vals import *
from .misc.units_tool import *
from .specs import *


DB_PATH = './sos/data/mc_db.yaml'
LOGO_SOS = imread('./sos/res/logo_sos.png')

