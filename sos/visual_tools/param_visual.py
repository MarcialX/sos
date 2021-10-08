# -*- coding: utf-8 -*-
# --------------------------------------------------------------------------------- #
# Software de Observaciones Sintéticas S.O.S.
# Parameter visualiser
#
# Marcial Becerril, @ 28 September 2021
# Latest Revision: 07 Oct 2021, 19:10 GMT-6
#
# For all kind of problems, requests of enhancements and bug reports, please
# write to me at:
#
# mbecerrilt92@gmail.com
# mbecerrilt@inaoep.mx
#
# --------------------------------------------------------------------------------- #

import os

from astropy.io import fits

from PyQt5 import QtCore, QtWidgets, uic, QtGui
from PyQt5.QtCore import Qt, QObject, QThread
from PyQt5.QtWidgets import QApplication, QWidget, QMessageBox
from PyQt5.QtGui import QPixmap, QIcon

from ..fits_funcs.table_header import TableModelHeader
from ..mc_plotter import *


class ParamVisualiser(QWidget):
    """
        Interactive window to fit spectral lines
        Parameters
        ----------
        Not initial params
        ----------
    """
    # Signal to update data
    save = QtCore.pyqtSignal(bool)

    def __init__(self):

        super(ParamVisualiser, self).__init__()

        uic.loadUi("./sos/res/gui/param_visualiser.ui", self)

        # Bin parameters
        self.bin_params = [['Bin name', ''], ['LTE mass [Ms]', ''], ['Virial mass [Ms]', ''],
                           ['XF mass [Ms]', ''], ['H2 col dens [cm^-2]', ''], ['Volume [pc^3]', ''],
                           ['Opacity', ''], ['Tex [K]', ''], ['Re [pc]', ''], ['Position', '']]

        # Get params
        self.GRAL_PARAMS = {'LTE Mass':'mass_lte', 'Virial Mass':'mass_vir', 'XF Mass':'mass_xf',
                            'H2 Column density':'NH2', 'Magnetic field':'B', 'Excitation Temperature': 'Tex',
                            'Volume':'vol', 'Volume Density':'den'}

        # Get initial param to plot
        self.param = 'mass_lte'
        self.paramsComboBox.currentTextChanged.connect(self.on_param_changed)

        self.color = 'magma'
        self.colorBox.currentTextChanged.connect(self.on_color_changed)

        self.scale = 'Lineal'
        self.scaleBox.currentTextChanged.connect(self.on_scale_changed)

        self.contour = True
        self.contoursBox.stateChanged.connect(self.change_map_contours)

        self.edit_flag = False
        self.flag_editButton.mousePressEvent = self.editing_map

        self.levels = self.levelsBox.value()
        self.levelsBox.valueChanged.connect(self.change_contour_levels)

        # Initialise some variables
        self.currentPos = ''
        self.fix = False

        self.flag_editButton.setStyleSheet('QPushButton {\
                        color: white;\
                        background-color:rgb(87, 195, 237);\
                        border-style: inset;\
                    }')


    def editing_map(self, event):
        # Change edit flag
        self.edit_flag = not self.edit_flag
        if self.edit_flag:
            self.flag_editButton.setText('Editing...')
            self.flag_editButton.setStyleSheet('QPushButton {\
                color: white;\
                background-color:red;\
                border-style: inset;\
            }')
        else:
            self.flag_editButton.setText('Flag Bin Edit')
            self.flag_editButton.setStyleSheet('QPushButton {\
                color: white;\
                background-color:rgb(87, 195, 237);\
                border-style: inset;\
            }')          


    def change_map_contours(self, event):
        # Get button status
        self.contour = self.contoursBox.isChecked()

        # Update main graph
        self.update_main_plot()


    def change_contour_levels(self, event):
        # Get number of contour levels
        self.levels = self.levelsBox.value()

        # Update main graph
        self.update_main_plot()


    def on_color_changed(self, event):
        # Get current color selected
        self.color = self.colorBox.currentText()

        # Update main graph
        self.update_main_plot()


    def on_scale_changed(self, event):
        # Get current color selected
        self.scale = self.scaleBox.currentText()

        # Update main graph
        self.update_main_plot()


    def on_param_changed(self, event):
        # Get current param selected
        param = self.paramsComboBox.currentText()
        self.param = self.GRAL_PARAMS[param]

        # Update main graph
        self.update_main_plot()

        # Update histogram
        data_hist = self.get_hist(self.param)
        self.ax_hist.cla()
        self.ax_hist.hist(data_hist, 30, color='b')
        self.fig_hist.canvas.draw_idle()

        self.update_stats(data_hist)


    def update_stats(self, data_hist):
        # Get statistics
        max_data = np.nanmax(data_hist)
        min_data = np.nanmin(data_hist)
        mean_data = np.nanmean(data_hist)
        med_data = np.nanmedian(data_hist)
        sum_data = np.nansum(data_hist)

        if mean_data > 1e4:
            max_data = '{:.2e}'.format(max_data)
            min_data = '{:.2e}'.format(min_data)
            mean_data = '{:.2e}'.format(mean_data)
            med_data = '{:.2e}'.format(med_data)
            sum_data = '{:.2e}'.format(sum_data)
        else:
            max_data = '{:.2f}'.format(max_data)
            min_data = '{:.2f}'.format(min_data)
            mean_data = '{:.2f}'.format(mean_data)
            med_data = '{:.2f}'.format(med_data)
            sum_data = '{:.2f}'.format(sum_data)

        stats_text = 'Sum='+sum_data+'     Max='+max_data+'     Min='+min_data+'     Mean='+mean_data+'     Median='+med_data
        self.statsLabel.setText(stats_text)


    def update_main_plot(self):
        # Show initial graphs
        log = False
        if self.scale == 'Logarithmic':
            log = True

        self.fig, self.ax = map_param(self.binned, self.param, self.moment, cmap=self.color, log=log, \
                            left=0.135, right=0.99, top=0.95, bottom=0.05, return_figure=True, \
                            show_contours=self.contour, level_contours=self.levels)
        self._rmmpl_main()
        self._addmpl_main(self.fig)

        # Events
        self._hover_event = self.fig.canvas.mpl_connect('motion_notify_event', self._onhover)
        self._onclick_xy = self.fig.canvas.mpl_connect('button_press_event', self._onclick)


    def load_fits_file(self, mol, moment, mc_binned, name):
        # Load initial params
        self.mol = mol
        self.moment = moment
        self.binned = mc_binned
        self.name = name

        # Update map name
        self.nameLabel.setText(self.name)

        self.model = TableModelHeader(self.bin_params)
        self.tableHeader.setModel(self.model)

        # Set font size
        font = QtGui.QFont("Courier New", 11)
        self.tableHeader.setFont(font)

        # Hide table headers
        self.tableHeader.horizontalHeader().hide()
        self.tableHeader.verticalHeader().hide()
        
        # Set table dimensions
        self.tableHeader.setColumnWidth(0, 150)   
        self.tableHeader.setColumnWidth(1, 100)

        # Enable some items
        all_params = [self.paramsComboBox.itemText(i) for i in range(self.paramsComboBox.count())]
        for i, a in enumerate(all_params):
            if not self.GRAL_PARAMS[a] in self.binned['B0'].keys():
                self.paramsComboBox.model().item(i).setEnabled(False)

        # Show initial graphs
        self.fig, self.ax = map_param(self.binned, self.param, self.moment, cmap='magma', log=False, \
                            left=0.135, right=0.99, top=0.95, bottom=0.05, return_figure=True)
        self._addmpl_main(self.fig)

        self._hover_event = self.fig.canvas.mpl_connect('motion_notify_event', self._onhover)
        self._onclick_xy = self.fig.canvas.mpl_connect('button_press_event', self._onclick)

        # Initial spectral graph
        self.fig_spectra, self.ax_spectra = subplots()
        self.fig_spectra.subplots_adjust(left=0.125, bottom=0.135, right=0.96,
                                top=0.98, wspace=None, hspace=None)
        self._addmpl_spectra(self.fig_spectra)

        self.ax_spectra.tick_params(axis='both', labelsize=7)

        # Initial historgram graph
        self.fig_hist, self.ax_hist = subplots()
        self.fig_hist.subplots_adjust(left=0.125, bottom=0.195, right=0.96,
                                top=0.98, wspace=None, hspace=None)
        self._addmpl_hist(self.fig_hist)
        self.ax_hist.tick_params(axis='both', labelsize=7)

        # Get histogram
        data_hist = self.get_hist(self.param)
        self.ax_hist.hist(data_hist, 30, color='b')
        self.fig_hist.canvas.draw_idle()

        self.update_stats(data_hist)


    def get_hist(self, param):
        # Get data
        hist_data = []
        for b in self.binned.keys():
            if not self.binned[b]['flag']:
                vals = self.binned[b][param]
                if isinstance(vals, np.ndarray):
                    vals = np.mean(vals)
                hist_data.append(vals)

        return np.array(hist_data)


    def _onhover(self, event):
        """
            On hover event to select bins
        """
        if event.inaxes == self.ax and not self.fix:
            ix, iy = event.xdata, event.ydata
            if not((ix is None) or (iy is None)):

                # Get bin
                bin_name = self.get_bin_name(ix, iy)

                if bin_name != '' and self.currentPos != bin_name:
                    # Plot spectra
                    velx = self.binned[bin_name][self.mol]['vel']
                    line = self.binned[bin_name][self.mol]['line']

                    self.ax_spectra.cla()
                    self.ax_spectra.plot(velx, line, 'r', lw=1.5)

                    self.fig_spectra.canvas.draw_idle()

                    # Update table data
                    self.bin_params[0][1] = bin_name
                    params = ['mass_lte', 'mass_vir', 'mass_xf', 'NH2',
                              'vol', 'tau', 'Tex', 'Re', 'pos']
                    for i in range(len(params)):
                        vals = self.binned[bin_name][params[i]]
                        if isinstance(vals, np.ndarray):
                            if len(vals) > 0:
                                vals = np.mean(vals)
                                if vals > 1e4:
                                    vals = '{:.2e}'.format(vals)
                                else:
                                    vals = '{0:.2f}'.format(vals)
                            else:
                                vals = np.nan
                        elif isinstance(vals, list):
                            vals = str(vals[0])+','+str(vals[1])
                        else:
                            if vals > 1e4:
                                vals = '{:.2e}'.format(vals)
                            else:
                                vals = '{0:.2f}'.format(vals)

                        self.bin_params[i+1][1] = vals

                    self.model = TableModelHeader(self.bin_params)
                    self.tableHeader.setModel(self.model)

                    self.currentPos = bin_name


    def get_bin_name(self, ix, iy):
        # Get position
        xpos, ypos = int(np.round(ix)), int(np.round(iy))

        # Get bin
        bin_name = ''
        for b in self.binned.keys():
            xbin, ybin = self.binned[b]['pos']
            if xbin == xpos and ybin == ypos:
                bin_name = b
                break

        return bin_name


    def _onclick(self, event):
        """
            On click event to select lines
        """
        if event.inaxes == self.ax:
            ix, iy = event.xdata, event.ydata
            # Left-click
            if event.button == 1:
                if self.edit_flag:
                    # Get bin name
                    bin_name = self.get_bin_name(ix, iy)
                    # Change flag
                    self.binned[bin_name]['flag'] = not self.binned[bin_name]['flag']
                    # Update main graph
                    self.update_main_plot()

                    # Update histogram
                    data_hist = self.get_hist(self.param)
                    self.ax_hist.cla()
                    self.ax_hist.hist(data_hist, 30, color='b')
                    self.fig_hist.canvas.draw_idle()

                    self.update_stats(data_hist)

                else:
                    if not self.fix:
                        # Get position
                        xpos, ypos = int(np.round(ix)), int(np.round(iy))
                        self.ax.plot(xpos, ypos, 'k*')
                        self.fig.canvas.draw_idle()
                    # Add peaks
                    self.fix = True

            # Right-click
            elif event.button == 3:
                #ix, iy = event.xdata, event.ydata
                self.fix = False
                self.update_main_plot()


    def _addmpl_hist(self, fig):
        
        self.canvas_hist = FigureCanvas(fig)
        self.histPlotLayout.addWidget(self.canvas_hist)
        self.canvas_hist.draw()
        #self.toolbar_hist = NavigationToolbar(self.canvas_hist,
        #   self, coordinates=True)
        #self.histPlotLayout.addWidget(self.toolbar_hist)


    def _rmmpl_hist(self):

        self.histPlotLayout.removeWidget(self.canvas_hist)
        self.canvas_hist.close()
        #self.histPlotLayout.removeWidget(self.toolbar_hist)
        #self.toolbar_hist.close()


    def _addmpl_spectra(self, fig):
        
        self.canvas_spec = FigureCanvas(fig)
        self.spectraPlotLayout.addWidget(self.canvas_spec)
        self.canvas_spec.draw()
        #self.toolbar_spec = NavigationToolbar(self.canvas_spec,
        #   self, coordinates=True)
        #self.spectraPlotLayout.addWidget(self.toolbar)


    def _rmmpl_spectra(self):

        self.spectraPlotLayout.removeWidget(self.canvas_spec)
        self.canvas_spec.close()
        self.spectraPlotLayout.removeWidget(self.toolbar_spec)
        self.toolbar_spec.close()


    def _addmpl_main(self, fig):
        
        self.canvas_main = FigureCanvas(fig)
        self.mainPlotLayout.addWidget(self.canvas_main)
        self.canvas_main.draw()
        self.toolbar_main = NavigationToolbar(self.canvas_main,
           self, coordinates=True)
        self.mainPlotLayout.addWidget(self.toolbar_main)


    def _rmmpl_main(self):

        self.mainPlotLayout.removeWidget(self.canvas_main)
        self.canvas_main.close()
        self.mainPlotLayout.removeWidget(self.toolbar_main)
        self.toolbar_main.close()
