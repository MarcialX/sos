# -*- coding: utf-8 -*-
# --------------------------------------------------------------------------------- #
# Software de Observaciones Sintéticas S.O.S.
# Ejemplo versión [15-Sep-2021]. Implementación de filtro de ruido
#
# Marcial Becerril, @ 15 September 2020
# Latest Revision: 15 September 2021, 01:00 GMT-6
##
# For all kind of problems, requests of enhancements and bug reports, please
# write to me at:
#
# mbecerrilt92@gmail.com
# mbecerrilt@inaoep.mx
#
# --------------------------------------------------------------------------------- #

# Librería para graficar
from matplotlib.pyplot import *

# Librería numérica
import numpy as np

# Import sos library
import sos

# Creamos el objeto nube molecular
mc = sos.mc('3Orion')

# RUIDO
# ===========
# Calculamos el ruido utilizando los primeros 5 canales y los últimos 5
noise_channels = np.arange(-5,5)
# Calulamos el ruido para todas las moleculas.
# NOTA. El verbose (que significa verborrea) solo le indica a la función
# que ademas de hacer el cálculo, también lo imprima en la pantalla
mc.get_noise('13CO', channels=np.arange(-5, 5), verbose=True)
mc.get_noise('12CO', channels=np.arange(-5, 5), verbose=True)

# REGIÓN COMPLETA
# ===========
# Obtenemos los espectros integrados de todas las moléculas
mc.get_map_vels()
# Ajustamos una sola línea al espectro a las moléculas 13CO y 12CO
mc.line_fit('13CO', forced_lines=1)
mc.line_fit('12CO', forced_lines=1)
# Calculate the physical parameters of the whole cloud
mc.get_gral_params()

# REGIÓN SEGMENTADA
# ===============
# Segmentamos el mapa en 62x62 bines, la máxima posible para esta imagen, 
# es decir, igual a la resolución del mapa
nbins = 62
mc.binning_mol('13CO', nbins=nbins, rebin=True)
mc.binning_mol('12CO', nbins)
# Ajustamos una linea a cada espectro. Se define un umbral de detección 'sigma_thresh',
# si esta por debajo se ignora el pixel. Este valor normalmente es más chico a medida
# que incrementamos la resolución. Sugiero probar con este e irlo cambiando para ver
# los cambios
mc.line_fit_binning('13CO', forced_lines=1, sigma_thresh=0.4)
mc.line_fit_binning('12CO', forced_lines=1, sigma_thresh=0.4)
# Calculamos los parámetros de cada pixel.
# El parámetro 'no_filter' sólo le indica al programa que no aplique suavizado a los
# mapas. Esto se recomienda para cubos de datos que contengan valores vaciós (nan)
# como es el caso de este mapa que es circular y las esquinas son nan.
mc.get_bins_params(no_filter=True)
# Con el criterio de sólo usar espectros con cierto umbral los resultados parecen
# no necesitar de este filtro. Sin embargo, dependerá de la situación.
#mc.param_filter('mass_lte', max=1e4)

# GRÁFICAR LOS RESULTADOS
# =================
# Obtenemos el mapa de momento cero de la molécula 13CO
# 'n' es para indicar el grado del momento
m0_13co = mc.get_n_moment('13CO', n=0)
# Obtenemos el mapa de momento con los espectros de cada bin sobrepuestos
sos.plot_moment_spectra(m0_13co, mc.extract_param_from_binned('13CO'), label=False)
# Graficamos la distribución espacial de la masa calculada
sos.map_param(mc.binned, 'mass_lte', m0_13co, cmap='RdBu_r', log=False)

# Presentamos los resultados finales generales
mc.summary()

# Sumamos la masa de todos los bines activos
mc.sum_binning_param('mass_lte')

# Desplegamos los mapas
show()