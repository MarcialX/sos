# -*- coding: utf-8 -*-
# --------------------------------------------------------------------------------- #
# Software de Observaciones Sintéticas S.O.S.
# specs. Module to read, analyse (substract baseline, fit gauss/lorentzian, etc.), save parameters
#
# Marcial Becerril, @ 25 January 2021
# Latest Revision: 25 Jan 2021, 13:38 GMT
#
# For all kind of problems, requests of enhancements and bug reports, please
# write to me at:
#
# mbecerrilt92@gmail.com
# mbecerrilt@inaoep.mx
#
# --------------------------------------------------------------------------------- #

import csv

from matplotlib.pyplot import *
matplotlib.use('Qt5agg')
#ion()

from matplotlib.backends.backend_qt5agg import(
    FigureCanvasQTAgg as FigureCanvas,
    NavigationToolbar2QT as NavigationToolbar)


import numpy as np
from scipy.integrate import simps

from PyQt5 import QtWidgets

from .misc.print_msg import *
from .line_fitting import *

import time
import sys


class specs(object):
    """
        Spectra object
        Parameters
        ----------
        filename : string
            File name path of the spec file (As .csv file, so far...)
        show :  boolean
            Show spectra during the analysis?
        ----------
    """
    def __init__(self, filename, delimiter=' ', show=False):

        # Load molecular cloud database
        self.pathSpec = filename

        # Initialise fitting container
        self.fit_data = {}

        self.x, self.y = 0, 0
        self.units = ['um', '']
        self.x, self.y = self.load_data(filename, delimiter)
        
        self.baseline_substracted = None

        self.show = show
        if show:
            self.plot_spectra()


    def load_data(self, path, delimiter=' '):
        """
            Load spectra data
            Parameters
            ----------
            filename : string
                Read spectra data (as .csv)
            ----------
        """
        try:
            with open(path, newline='') as csvfile:
                reader = csv.reader(csvfile, delimiter=delimiter)
                x = []
                y = []
                for i, row in enumerate(reader):
                    nrow = []
                    for r in row:
                        if r != '': 
                            nrow.append(r)
                    x.append(float(nrow[0]))
                    y.append(float(nrow[1]))

            msg("Spectra file read!", 'ok')
            return x, y
        except:
            msg("File couldn't be read. Check the format", 'fail')
            return 

    
    def substract_baseline(self, inter=True, *args, **kwargs):
        """
            Substract baseline
            Parameters
            ----------
            inter: boolean
                Activate interactive mode
            *args
            *kwargs:
                ----- Baseline fitting -----
                method [string]: Method to fit baseline:
                    linear, poly, bls
                ndeg [int]: Degree of polynom. Only for linear and poly methods
                If method is 'bls' the following params have to be defined
                    l [int]: Smoothing
                    p [float]: Asymmetry
                    n [int]: Number of iterations
                ----- Find peaks -----
                dist [float]: distance between peaks
                height_div [float]: height factor to set minimum limit
                ----- Line width -----
                klw [float]: Linewidth constant factor
            ----------
        """
        if inter:
            # Add canvas to the app
            ioff()
            fig, ax = subplots()
            close(fig)

            # Instance app
            app = QApplication.instance()

            self.SubBaseLine = BaselineSpecWindow()
            save = kwargs.pop('save', False)
            name = self.pathSpec.split('/')[-1]
            self.SubBaseLine.load_init_params(fig, ax, self.x, self.y, 'spectra', self.units, name, save=save)
            # Signal connection to extract baseline data
            self.SubBaseLine.signal_baseline.connect(self._get_baseline_from_ui)
            self.SubBaseLine.show()

            app.exec_()

        else:
            # Get the method from the kwargs
            method = kwargs.pop('method', 'linear')
            if method == 'linear':
                degree = 1
            elif method == 'poly':
                degree = kwargs.pop('ndeg', 2)
            elif method == 'bls':
                l = kwargs.pop('l', 105)
                p = kwargs.pop('p', 0.05)
                n = kwargs.pop('n', 10)
            else:
                msg('Define a valid method to fit the baseline:\n1. Linear\n2. Polynomial\n3. BLS method', 'fail')
                return

            # Length size
            Nsize = len(self.y)
            # Get preliminar baseline
            edges = int(Nsize/20)
            prem_baseline = poly_baseline(np.concatenate((self.x[:edges], self.x[Nsize-edges:])),
                           np.concatenate((self.y[:edges], self.y[Nsize-edges:])), 1, self.x)
            # Get find peaks params from kwargs
            dist = kwargs.pop('dist', 5.)
            height_div = kwargs.pop('height_div', 12.)
            # Get peaks and widths
            peaks = find_profile_peaks(self.y-prem_baseline, dist=dist, height_div=height_div)
            # Guess line parameters to remove
            guess = guess_line(np.array(self.x), np.array(self.y-prem_baseline), peaks)
            baseline = self.y.copy()
            dv = np.mean(np.diff(self.x))
            # Define the line width constant
            klw = kwargs.pop('klw', 6)

            rem_items = np.array([])
            for i, peak in enumerate(peaks):
                line_width = klw*int(guess[3*i+2]/dv)
                # Inferior limit
                lim_inf = peak-line_width/2
                if lim_inf < 0:
                    lim_inf = 0
                # Superior limit
                lim_sup = peak+line_width/2
                if lim_sup > Nsize:
                    lim_sup = Nsize

                rem_aux_items = np.arange(lim_inf, lim_sup)
                rem_items = np.concatenate((rem_items, rem_aux_items))

            baseline = np.delete(baseline, rem_items.astype(int))
            freq_baseline = np.delete(self.x, rem_items.astype(int))

            try: 
                if method == 'bls':
                    fit_baseline = baseline_als_optimized(baseline, l, p, niter=n)
                elif (method == 'linear') or (method == 'poly'):
                    fit_baseline = poly_baseline(freq_baseline, baseline, degree, self.x)
                
                self.baseline_substracted = fit_baseline
            except Exception as e:
                msg('Baseline couldnt be adjusted.'+str(e), 'fail')    

    
    def _get_baseline_from_ui(self, kind):
        """
            Get data with baseline substracted from the UI
        """
        self.baseline_substracted = self.SubBaseLine.data_corrected


    def fit_spectra(self, data=None, inter=True, *args, **kwargs):
        """
            Fit lines into a spectrum
            Parameters
            ----------
            data : array
                Y component of the spectra.
                If it is None, it uses the data from the
                initially loaded file.
            inter : boolean
                Activate interactive mode?
            ----------
        """

        if not data:
            if self.baseline_substracted is None:
                spectra = self.y
            else:
                spectra = self.baseline_substracted
        else:
            spectra = data

        if inter:
            # Add canvas to the app
            ioff()
            fig, ax = subplots()
            close(fig)

            # Instance app
            app = QApplication.instance()

            self.FitSpectra = FitLinesWindow()
            save = kwargs.pop('save', False)
            name = self.pathSpec.split('/')[-1]
            self.FitSpectra.load_init_params(fig, ax, self.x, spectra, 'spectra', self.units, name, save=save)
            # Enable selection button
            self.FitSpectra.inter = True
            self.FitSpectra._change_inter_icon()
            # Signal connection to extract baseline data
            self.FitSpectra.signal_fitting.connect(self._get_fitting_from_ui)
            self.FitSpectra.show()

            app.exec_()

        else:

            # Get find peaks params from kwargs
            dist = kwargs.pop('dist', 100.)
            height_div = kwargs.pop('height_div', 5.)
            # Find peaks
            peaks = find_profile_peaks(spectra, dist=dist, height_div=height_div)
                       
            # The automatic mode only works with purely Gaussian funtions
            lines_method = ['G']*len(peaks)

            if len(peaks) > 0:
                popt, pcov = solve_line(np.array(self.x), spectra, peaks, lines_method)
                A, mu, sigma, fwhm = get_params_from_popt(popt, len(peaks))
            else:
                msg('No lines found', 'warn')
                return
            
            # Assign to the fitting results dictionary
            self.fit_data = {}
            for i in range(len(A)):
                param = [A[i], mu[i], sigma[i]]
                curve, area = get_fit_curve(self.x, lines_method[i], param)
                if lines_method[i] == 'G':
                    line_func = 'Gaussian'
                elif lines_method[i] == 'L':
                    line_func = 'Lorentzian'
                self.fit_data['L'+str(i+1)] = [A[i], mu[i], sigma[i], area, line_func]


    def _get_fitting_from_ui(self):
        """
            Get lines fitted to the spectra from the UI
        """
        self.fit_data = self.FitSpectra.fit_data_dict


    def summary(self):
        """
            Show the results of the fitting
        """
        msg ('------------------------------------------------', 'info')
        msg ('            Summary Spectra Fitting             ', 'info')
        msg ('------------------------------------------------', 'info')
        for n in self.fit_data.keys():
            print ('Line: ', n)
            print ('Method: ',self.fit_data[n][4] )
            print ('Amplitude: ', np.round(self.fit_data[n][0],2), self.units[1])
            print ('Position: ', np.round(self.fit_data[n][1],2), self.units[0])
            print ('Width: ', np.round(self.fit_data[n][2],2), self.units[0])
            print ('Area: ', np.round(self.fit_data[n][3],2))
            print('------------------------------------------------')


    def generate_report(self, path=None):
        """
            Generate report of results
            Parameters
            ----------
            path: string
                Path of the report
            ----------
        """
        if path is None:
            name = self.pathSpec.split('/')[-1]
            name = name.split('.')[0]
            path = './'+name+'_fit_report.txt'

        try:
            with open(path, 'w') as file:
                header = "Name\tMethod\tAmplitude "+self.units[1]+"\tPosition "+self.units[0]+"\tWidth "+self.units[0]+"\tArea\n"
                file.write(header)
                for n in self.fit_data.keys():         
                    row = n+"\t"+self.fit_data[n][4]+"\t"+str(np.round(self.fit_data[n][0],4))+ \
                          "\t"+str(np.round(self.fit_data[n][1],4))+"\t"+str(np.round(self.fit_data[n][2],4))+ \
                          "\t"+str(np.round(self.fit_data[n][3],4))+"\n"
                    file.write(row)
            msg('Done :)', 'ok')
        except Exception as e:
            msg("File couldn't be write."+str(e), 'fail')
            return 


    def plot_spectra(self):
        """
            Plot the spectra loaded
        """
        figure()
        plot(self.x, self.y)
        show()