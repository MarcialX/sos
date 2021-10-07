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
mc.get_noise('13CO', channels=noise_channels, verbose=True)
mc.get_noise('12CO', channels=noise_channels, verbose=True)

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
# Segmentamos el mapa en 63x63 bines, la máxima posible para esta imagen, 
# es decir, igual a la resolución del mapa
nbins = 8
mc.binning_mol('13CO', nbins=nbins, rebin=True)
mc.binning_mol('12CO', nbins)
# Ajustamos una linea a cada espectro. Se define un umbral de detección 'sigma_thresh',
# si esta por debajo se ignora el pixel. Este valor normalmente es más chico a medida
# que incrementamos la resolución. Sugiero probar con este e irlo cambiando para ver
# los cambios
mc.line_fit_binning('13CO', forced_lines=1, sigma_thresh=1)
mc.line_fit_binning('12CO', forced_lines=1, sigma_thresh=1)
# Calculamos los parámetros de cada pixel.
# El parámetro 'no_filter' sólo le indica al programa que no aplique suavizado a los
# mapas. Esto se recomienda para cubos de datos que contengan valores vacíos (nan)
# como es el caso de este mapa que es circular y las esquinas son nan.
mc.get_bins_params()
# Con el criterio de sólo usar espectros con cierto umbral los resultados parecen
# no necesitar de este filtro. Sin embargo, dependerá de la situación.
#mc.param_filter('mass_lte', max_lim=10)

# GRÁFICAR LOS RESULTADOS
# =================
# Obtenemos el mapa de momento cero de la molécula 13CO
# 'n' es para indicar el grado del momento
m0_13co = mc.get_n_moment('13CO', n=0)
# Obtenemos el mapa de momento con los espectros de cada bin sobrepuestos
# [Esta linea se comenta para que sea más rápido el programa]
#sos.plot_moment_spectra(m0_13co, mc.extract_param_from_binned('13CO'), label=False)

# [FILTRO ADICIONAL]
# Podemos agregar filtros externos, activando y desactivando, según algún criterio, los bines.
# Por ejemplo, en este caso, filtraremos según la intensidad del mapa de momento zero
# Quitaremos del mapa aquellos píxeles con intensidad < 2.5 K km/s
m0_thresh = 2500
data = m0_13co['data']
for name in mc.binned.keys():
	# Full pixels
	nbin = int(name[1:])
	xmin, xmax, ymin, ymax = mc.bin_grid[nbin]
	full_samples = len(data[xmin:xmax,ymin:ymax].flatten())
	full_pixel = np.nansum(data[xmin:xmax,ymin:ymax])
	# Frame pixels
	val = 0
	fm_samples = 0
	for fm in mc.bin_frames[nbin]:
		val += fm[2]*data[fm[0], fm[1]]
		fm_samples += fm[2]

	m0_avg = (full_pixel+val)/(fm_samples+full_samples)

	if m0_avg < m0_thresh:
		mc.binned[name]['flag'] = True

# Graficamos la distribución espacial de la masa calculada
sos.map_param(mc.binned, 'mass_lte', m0_13co, cmap='magma', log=False)

# Presentamos los resultados finales generales
mc.summary()

# Sumamos la masa de todos los bines activos
mc.sum_binning_param('mass_lte')

# Obtenemos la estadística de la masa
mean_mass, med_mass, std_mass, data_mass = mc.get_param_stat('mass_lte')
# Y densidad columnar
mean_nh2, med_nh2, std_nh2, data_nh2 = mc.get_param_stat('NH2')

# Plot the histograms
fig, ax = subplots(1, 2)

rcParams['axes.linewidth'] = 1.25
rc('font', family='serif', size='18')

# LTE Mass
ax[0].hist(data_mass, 50, color='r', linewidth=1.5)

legend_elements = [Line2D([0], [0], marker='o', color='k', label='LTE mass',
              markerfacecolor='r', markersize=18)]

ax[0].legend(handles=legend_elements, loc='best', prop={'size': 18})   
ax[0].set_xlabel(r"Masa LTE [M$\odot$]")
ax[0].set_ylabel(r"Número de muestras")
ax[0].grid()

# NH2
ax[1].hist(data_nh2, 50, color='b', linewidth=1.5)

legend_elements = [Line2D([0], [0], marker='o', color='k', label='NH2',
              markerfacecolor='b', markersize=18)]

ax[1].legend(handles=legend_elements, loc='best', prop={'size': 18})   
ax[1].set_xlabel(r"Densidad columnar NH2 [cm$^{-2}$]")
ax[1].set_ylabel(r"Número de muestras")
ax[1].grid()

# Desplegamos los mapas
show()