# -*- coding: utf-8 -*-
# --------------------------------------------------------------------------------- #
# Software de Observaciones Sintéticas S.O.S.
# Usage example of spectra module of sos
#
# Marcial Becerril, @ 12 March 2021
# Latest Revision: 12 March 2021, 15:37 GMT-6
##
# For all kind of problems, requests of enhancements and bug reports, please
# write to me at:
#
# mbecerrilt92@gmail.com
# mbecerrilt@inaoep.mx
#
# --------------------------------------------------------------------------------- #

# Importamos la librería SOS
import sos

# Definimos la dirección del archivo con el espectro [csv]
path = './data_specs/MYS1192.csv'

# Creamos el objecto espectro
sp = sos.specs(path)

# Sustraemos la línea base [Modo interactivo]
sp.substract_baseline()

# Ajustamos las líneas disponibles en el espectro libre de contaminación [línea base]
sp.fit_spectra()

# Mostramos los resultados
sp.summary()

# Generamos el reporte
sp.generate_report()


