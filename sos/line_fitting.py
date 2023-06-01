# -*- coding: utf-8 -*-
# --------------------------------------------------------------------------------- #
# Software de Observaciones Sintéticas S.O.S.
# Line fitting functions
#
# Marcial Becerril, @ 26 August 2020
# Latest Revision: 26 Aug 2020, 22.58 GMT
#
# For all kind of problems, requests of enhancements and bug reports, please
# write to me at:
#
# mbecerrilt92@gmail.com
# mbecerrilt@inaoep.mx
#
# --------------------------------------------------------------------------------- #


import numpy as np

from matplotlib import colors
from matplotlib.pyplot import *
from scipy.optimize import curve_fit

from PyQt5 import QtCore, QtWidgets, uic, QtGui
from PyQt5.QtCore import Qt, QObject, QThread
from PyQt5.QtWidgets import QApplication, QWidget, QMessageBox
from PyQt5.QtGui import QPixmap, QIcon

import sos
from .misc.line_functions import *
from .misc.print_msg import *
from .misc.table_model import *

from scipy.signal import find_peaks
from scipy.signal import savgol_filter
from scipy import sparse
from scipy.sparse.linalg import spsolve
from scipy.integrate import simps

from datetime import datetime

from matplotlib.backends.backend_qt5agg import(
    FigureCanvasQTAgg as FigureCanvas,
    NavigationToolbar2QT as NavigationToolbar)


def find_profile_peaks(data_array, **kwargs):
    """
        Find profile peaks
        Parameters
        ----------
            data_array : array
            dist : float
                Minimum distance between peaks
            height_div : float
                Factor which divide the maximum value of the array, to define the
                minimum peak detectable
        ----------
    """

    # Find peaks keyword parameters
    # Distance between lines
    dist = kwargs.pop('dist', 10.0)
    # Minimum height for the lines
    height_div = kwargs.pop('height_div', 5.0)

    # Height division
    height = np.max(data_array)/height_div
    peaks, _ = find_peaks(data_array, distance=dist, height=height)

    return peaks


def get_fit_curve(x, func, params):
    """
        Calculate the area under the curve of the line function
        Parameters
        ----------
        x : x-axis data
        params : array
            Line function parameters
        func : string
            Function: 'G' Gaussian, or 'L' Lorentzian
        ----------
    """
    if func == 'G':
        curve = gaussian(x, *params)
    elif func == 'L':
        curve = lorentzian(x, *params)
    
    area = simps(curve, x)

    return curve, area


def multiple_lines(x, func, *params):
    """
        Multiple gaussian profiles
        Parameters
        ----------
            x : array
            params : array
                Parameters of the line
            func : function
                Line function model
        ----------
    """
    k = 0
    nparams = 3
    y = np.zeros_like(x)
    for i in range(0, len(params), nparams):
        a = params[i]
        b = params[i+1]
        c = params[i+2]

        if func[k] == 'G':
            f = gaussian
        elif func[k] == 'L':
            f = lorentzian

        y = y + f(x, a, b, c)
        k = k + 1

    #offset = params[-1]
    #y = y + offset

    return y


def guess_line(x, data, peaks, func=None):
    """
        Guess line parameters
        Parameters
        ----------
            x : array
            data : array
                Line profile
            peaks : array
                Peaks positions
            func : array of functions
                Functions of line profiles
        ----------
    """
    guess = []
    step = np.mean(np.diff(x))
    std_step = np.std(np.diff(x))

    # Smoothing the signal
    data = savgol_filter(data, 5, 2)

    for i in range(len(peaks)):
        # Amplitude
        A = data[peaks[i]]
        # Mean
        mu = x[peaks[i]]
        # Dispersion
        hm = np.abs(A/2.)
        hwm = np.where(data>hm)[0]

        if peaks[i] > 0:
            left_reg = np.where(np.abs(data[:peaks[i]])>hm)[0]
            if len(left_reg) == 0:
                left_edge = 0
            else:
                left_edge = left_reg[np.append(False, np.abs(np.diff(x[left_reg]))>(2*std_step+step))]
                if len(left_edge) == 0:
                    left_edge = np.where(np.abs(data[:peaks[i]])>hm)[0][0]
                else:
                    left_edge = left_edge[0]
        else:
            left_edge = 0

        if peaks[i] < len(x)-1:
            right_reg = np.where(np.abs(data[peaks[i]:])>hm)[0]
            if len(right_reg) == 0:
                right_edge = 0
            else:
                right_edge = right_reg[np.append(False, np.abs(np.diff(x[right_reg]))>(2*std_step+step))]
                if len(right_edge) == 0:
                    right_edge = np.where(np.abs(data[peaks[i]:])>hm)[0][-1]+peaks[i]
                else:
                    right_edge = right_edge[0]
        else:
            right_edge = 0

        fwhm = x[right_edge] - x[left_edge]
        if func is None:
            s = fwhm/(2*np.sqrt(2*np.log(2)))
        else:
            if func[i] == 'G':
                s = fwhm/(2*np.sqrt(2*np.log(2)))   
            elif func[i] == 'L':
                s = fwhm

        guess += [A, mu, np.abs(s)]

    # Normalise the guess [according the number of lines]
    # xmax_idx = np.argmax(data[peaks])
    # max_comb = multiple_lines(x[peaks[xmax_idx]], peaks, *guess)

    # Amax = np.max(guess[::3])
    # k = Amax/max_comb
    
    # guess[::3] = list(k*np.array(guess[::3]))

    # Get the bounds
    bounds_min = -1*np.ones_like(guess)*np.inf
    bounds_min[::3] = np.zeros_like(bounds_min[::3]) 
    bounds_max = np.ones_like(guess)*np.inf

    # Medians in the range
    ndata = np.abs(np.mean(np.diff(x)))

    vel_bounds = [x[0] - ndata/2., x[-1] - ndata/2.]
    min_bound = np.min(vel_bounds)
    max_bound = np.max(vel_bounds)

    bounds_min[1::3] = min_bound*np.ones_like(bounds_min[1::3])

    bounds_max = np.ones_like(guess)*np.inf
    bounds_max[1::3] = max_bound*np.ones_like(bounds_min[1::3])

    return guess, bounds_min, bounds_max


def solve_line(x, spectra, peaks, lines_methods):
    """
        Line fitting solver
        Parameters
        ----------
            x : array
            spectra : array
                Line profile
            peaks : array
                Peaks positions
            lines_methods : array of functions
                Functions of line profiles
        ----------
    """
    guess, bounds_min, bounds_max = guess_line(x, spectra, peaks, func=lines_methods)

    try:
        popt, pcov = curve_fit(lambda x, *params: multiple_lines(x, ['G']*len(peaks), *params), 
                     x, spectra, p0=guess, bounds=(bounds_min, bounds_max))

        return popt, pcov

    except Exception as e:
        msg('Fitting error. Check the data or the initial guess. '+str(e), 'warn')
        return None, None


def get_params_from_popt(popt, n_peaks):
    """
        Get parameters from the potp array
        Parameters
        ----------
            popt : array
                Array solutions
            n_peaks : array
                Number of peaks
        ----------
    """
    A = []
    mu = []
    sigma = []
    fwhm = []

    if not popt is None:
        for peak in range(n_peaks):
            popt[2+3*peak] = np.abs(popt[2+3*peak])

            A.append(popt[3*peak])
            mu.append(popt[3*peak+1])
            sigma.append(popt[3*peak+2])
            fwhm.append(2*np.sqrt(2*np.log(2))*popt[3*peak+2])

    return A, mu, sigma, fwhm


def poly_baseline(x, y, degree, base):
    """
        Fit baseline to a polynomial function
        Parameters
        ----------
            x : x-data to fit
            y : y-data to fit
            degree : polynomial degree
            base : x-axis base array
        ----------
    """
    z = np.polyfit(x, y, degree)
    p = np.poly1d(z)

    return p(base)


def baseline_als_optimized(y, l, p, niter=10):
    """
        Extract baseline using Asymmetric Least Squares Smoothing
        Eliers and Boelens 2005. Taken from:
        https://stackoverflow.com/questions/29156532/python-baseline-correction-library
        Parameters
        ----------
            y : y-data to fit
            l : smoothness. 102 <= l <= 109.
            p : asymmetry. 0.001 <= p <= 0.1 good choice for positive peaks.
            niter : number of iterations. 10 by default
        ----------
    """
    L = len(y)
    D = sparse.diags([1,-2,1],[0,-1,-2], shape=(L,L-2))
    D = l * D.dot(D.transpose())

    w = np.ones(L)
    W = sparse.spdiags(w, 0, L, L)
    for i in range(niter):
        W.setdiag(w)
        Z = W + D
        z = spsolve(Z, w*y)
        w = p*(y > z) + (1 - p)*(y < z)

    return z


class BaselineSpecWindow(QWidget):
    """
        Substract baselines
        Parameters
        ----------
            x : array
            y : array
            method : array
                Method of baseline substraction:
                1. line. Get the best fit line that fit several points
                2. polynomial N-degree. Best fit with a polynomial function of N-degree
                3. bls. BLS Algorithm
        ----------
    """
    # Signal to update data
    signal_baseline = QtCore.pyqtSignal(str)

    def __init__(self):

        super(BaselineSpecWindow, self).__init__()

        uic.loadUi("./sos/res/gui/baseline.ui", self)

        # Initialisation of variables
        self.inter = False
        self.selPointsFlag = False
        self.selRegionFlag = False
        self.method = None
        self.flag_apply = False
        self.iter_points = []

        # Assign buttons
        self.cancelButton.mousePressEvent = self.close_widget
        # Activate interaction
        self.interactButton.mousePressEvent = self.activate_interactive
        # Type of selection
        self.choosePointsButton.mousePressEvent = self.act_points_selection
        self.removePointsButton.mousePressEvent = self.act_region_selection
        # Selection method
        self.linealButton.mousePressEvent = self.linear_selection
        self.polyButton.mousePressEvent = self.poly_selection
        self.blsButton.mousePressEvent = self.bls_selection
        # Clear button
        self.clearButton.mousePressEvent = self.reset_canvas
        # Apply button
        self.applyButton.mousePressEvent = self.apply_baseline_correction
        # Accept button
        self.acceptButton.mousePressEvent = self.accept_baseline_correction

        #cid = self.f1.figure.canvas.mpl_connect('button_press_event', self.resDraw)


    def load_init_params(self, fig, ax, x, y, kind, units, name, save):
        # Load initial params
        self.x = x
        self.y = y

        # File name
        self.nameLabel.setText(name)

        # Initialise corrected data array
        self.data_corrected = self.y.copy()

        # Load units
        self.ux = units[0]
        self.uy = units[1]

        # Data type
        self.kind = kind

        # Get figure
        self.fig = fig
        self.fig.subplots_adjust(left=0.12, bottom=0.12, right=0.98,
                                top=0.98, wspace=None, hspace=None)
        self.ax = ax 
        
        # Save figure?
        self.save = save

        # Update plot
        self._addmpl(self.fig)

        # Initial plot
        self.initial_plot()


    def close_widget(self, event):
        # Disable graphs
        self.close()


    def activate_interactive(self, event):
        # Interactive activation
        self.inter = not self.inter
        if self.inter:
            if self.selPointsFlag or self.selRegionFlag:
                self._onclick_xy = self.fig.canvas.mpl_connect('button_press_event', self._onclick)
            else:
                msg('Choose one selection mode', 'warn')
                self.inter = not self.inter
                return

            icon_img = './sos/res/icons/int_sel.png'
        else:
            if self.selPointsFlag or self.selRegionFlag:
                self.fig.canvas.mpl_disconnect(self._onclick_xy)

            icon_img = './sos/res/icons/int.png'

        self.interactButton.setIcon(QIcon(icon_img))


    def act_points_selection(self, event):
        # Points selection activated
        self.selection_settings(True)
        self._update_selected_plot(self.iter_points)


    def act_region_selection(self, event):
        # Region selection activated
        self.selection_settings(False)
        self._update_selected_plot(self.iter_points)


    def linear_selection(self, event):
        # Linear baseline method
        self.baseline_method('linear')


    def poly_selection(self, event):
        # Linear baseline method
        self.baseline_method('poly')


    def bls_selection(self, event):
        # Linear baseline method
        self.baseline_method('bls')


    def baseline_method(self, method):

        self.method = method

        if self.method == 'linear':
            linear = './sos/res/icons/lineal_icon_sel.png'
            poly = './sos/res/icons/poly_curve.png' 
            bls = './sos/res/icons/bls_icon.png'    
            # Disable the other functions
            self.nDegreeBox.setEnabled(False)
            self.lamdbaBLSEdit.setEnabled(False)
            self.pBLSEdit.setEnabled(False)
            self.iterBLSEdit.setEnabled(False)

        elif self.method == 'poly':
            linear = './sos/res/icons/lineal_icon.png'
            poly = './sos/res/icons/poly_curve_sel.png' 
            bls = './sos/res/icons/bls_icon.png' 
            # Disable the other functions
            self.nDegreeBox.setEnabled(True)
            self.lamdbaBLSEdit.setEnabled(False)
            self.pBLSEdit.setEnabled(False)
            self.iterBLSEdit.setEnabled(False)

        elif self.method == 'bls':
            linear = './sos/res/icons/lineal_icon.png'
            poly = './sos/res/icons/poly_curve.png' 
            bls = './sos/res/icons/bls_icon_sel.png' 
            # Disable the other functions
            self.nDegreeBox.setEnabled(False)
            self.lamdbaBLSEdit.setEnabled(True)
            self.pBLSEdit.setEnabled(True)
            self.iterBLSEdit.setEnabled(True)

        self.linealButton.setIcon(QIcon(linear))
        self.polyButton.setIcon(QIcon(poly))
        self.blsButton.setIcon(QIcon(bls))


    def selection_settings(self, ptsFlag):
        # Grpah Selection Configuration
        self.selPointsFlag = ptsFlag
        self.selRegionFlag = not self.selPointsFlag

        if self.selPointsFlag:
            points = './sos/res/icons/choosePoints_sel.png'
            region = './sos/res/icons/removePoints.png'
        else:
            points = './sos/res/icons/choosePoints.png'
            region = './sos/res/icons/removePoints_sel.png'

        self.choosePointsButton.setIcon(QIcon(points))
        self.removePointsButton.setIcon(QIcon(region))


    def apply_baseline_correction(self, event):
        # Baseline correction
        if self.method == 'bls':
            l = float(self.lamdbaBLSEdit.text())
            p = float(self.pBLSEdit.text())
            n = int(float((self.iterBLSEdit.text())))
            self.iterBLSEdit.setText(str(n))
            y_baseline = baseline_als_optimized(self.y, l, p, niter=n)
        
        else:
            points = self.iter_points
            points = np.sort(points)

            if self.selPointsFlag:
                x_filtered = []
                y_filtered = []
                # Extracting data [points]
                for i in range(len(points)):
                    x_filtered.append(self.x[points[i]])
                    y_filtered.append(self.y[points[i]])

            elif self.selRegionFlag:
                # Extracting data [region]
                adquire_data = False
                x_mask = [True]*len(self.x)
                y_mask = [True]*len(self.y)
                for i in range(len(points)):
                    if adquire_data:
                        x_mask[points[i-1]:points[i]] = [False]*(points[i]-points[i-1])
                        y_mask[points[i-1]:points[i]] = [False]*(points[i]-points[i-1]) 
                    adquire_data = not adquire_data

                x_filtered = np.array(self.x)[x_mask]
                y_filtered = np.array(self.y)[y_mask]

            else:
                msg('Choose one selection mode', 'warn')
                return

            x_baseline = self.x

            if self.method == 'linear':
                y_baseline = poly_baseline(x_filtered, y_filtered, 1, x_baseline)

            elif self.method == 'poly':
                degree = self.nDegreeBox.value()
                y_baseline = poly_baseline(x_filtered, y_filtered, degree, x_baseline)
            
        # Set flag
        self.flag_apply = True

        data_corrected = self.y - y_baseline

        # Update data with baseline substrated
        self.data_corrected = data_corrected

        self._update_plot(y_baseline, data_corrected)


    def accept_baseline_correction(self, event):
        # Applying baseline correction
        if not self.flag_apply:
            self.apply_baseline_correction(event)
        
        self.signal_baseline.emit(self.kind)

        if self.save:
            now = datetime.now()
            name = now.strftime("%d-%m-%Y_%H-%M-%S")
            self.fig.savefig('fig_'+name+'_bl.png')

        self.close()


    def _onclick(self, event):
        """
            On click event to select lines
        """
        if event.inaxes == self.ax:
            # Left-click
            if event.button == 1:
                ix, iy = event.xdata, event.ydata
                # Add peaks
                xarray = np.where(self.x>ix)[0]
                if len(xarray) > 0:
                    xpos = xarray[0]
                else:
                    xpos = len(self.x)-1
                self.iter_points.append(xpos)

                self.flag_apply = False

            # Right-click
            elif event.button == 3:
                ix, iy = event.xdata, event.ydata
                popt = []
                # Remove points
                # Define threshold
                thresh = 5*np.mean(np.diff(self.x))
                xlines = np.where((ix >= (np.array(self.x)[self.iter_points] - thresh)) & 
                               (ix < (np.array(self.x)[self.iter_points] + thresh)))[0]

                try:
                    if len(xlines) > 0:
                        x_min = np.argmin(np.abs((np.array(self.x)[np.array(self.iter_points)[xlines]] - ix)))
                        if self.selRegionFlag:
                            self.iter_points.remove(self.iter_points[xlines[x_min]])
                        else:
                            ylines = np.where((iy >= (np.array(self.y)[self.iter_points] - thresh)) & 
                                              (iy < (np.array(self.y)[self.iter_points] + thresh)))

                            if len(ylines) > 0:
                                y_min = np.argmin(np.abs((np.array(self.y)[np.array(self.iter_points)[ylines]] - iy)))
                                self.iter_points.remove(self.iter_points[xlines[y_min]])
                        
                        self.flag_apply = False
                except:
                    pass

            # Update plot
            self._update_selected_plot(self.iter_points)


    def _update_selected_plot(self, points):
        """
            Update selected Points/Region in the canvas
        """
        self.ax.clear()

        # Label axes
        ux_label = self.ux
        if self.ux:
            ux_label = '['+ux_label+']'
        uy_label = self.uy
        if self.uy:
            uy_label = '['+uy_label+']'

        self.ax.set_xlabel(r''+ux_label)
        self.ax.set_ylabel(r'Temperature '+uy_label)

        self.ax.plot(self.x, self.y, 'k')
        for i in range(len(points)):
            if self.selPointsFlag:
                self.ax.plot(self.x[points[i]], self.y[points[i]], 'r+')
            elif self.selRegionFlag:
                self.ax.axvline(self.x[points[i]], color='r', linewidth=1)

        self.ax.grid()

        self.fig.canvas.draw_idle()


    def _update_plot(self, baseline, data_corrected):
        """
            Update baseline in the canvas
        """
        self.ax.clear()

        # Label axes
        ux_label = self.ux
        if self.ux:
            ux_label = '['+ux_label+']'
        uy_label = self.uy
        if self.uy:
            uy_label = '['+uy_label+']'

        self.ax.set_xlabel(r''+ux_label)
        self.ax.set_ylabel(r'Temperature '+uy_label)

        self.ax.plot(self.x, self.y, 'k', linewidth=0.75)
        self.ax.plot(self.x, baseline, 'c-.', linewidth=0.75)
        self.ax.plot(self.x, data_corrected, 'r')
        self.ax.grid()

        self.fig.canvas.draw_idle()


    def initial_plot(self):
        """
            Initial plot
        """
        self.ax.clear()
        self.ax.plot(self.x, self.y, 'k')

        # Label axes
        ux_label = self.ux
        if self.ux:
            ux_label = '['+ux_label+']'
        uy_label = self.uy
        if self.uy:
            uy_label = '['+uy_label+']'

        self.ax.set_xlabel(r''+ux_label)
        self.ax.set_ylabel(r'Temperature '+uy_label)

        self.ax.grid()

        self.fig.canvas.draw_idle()


    def reset_canvas(self, event):
        # Restart to initial plot
        self.initial_plot()

        self.iter_points = []


    def _addmpl(self, fig):
        
        self.canvas = FigureCanvas(fig)
        self.plotLayout.addWidget(self.canvas)
        self.canvas.draw()
        self.toolbar = NavigationToolbar(self.canvas,
           self, coordinates=True)
        self.plotLayout.addWidget(self.toolbar)


    def _rmmpl(self):
        self.plotLayout.removeWidget(self.canvas)
        self.canvas.close()
        self.plotLayout.removeWidget(self.toolbar)
        self.toolbar.close()



class FitLinesWindow(QWidget):
    """
        Interactive window to fit spectral lines
        Parameters
        ----------
            x : array
            y : array
            method : array
                Method of baseline substraction:
                1. line. Get the best fit line that fit several points
                2. polynomial N-degree. Best fit with a polynomial function of N-degree
                3. bls. BLS Algorithm
        ----------
    """
    # Signal to update data
    signal_fitting = QtCore.pyqtSignal(str)

    def __init__(self):

        super(FitLinesWindow, self).__init__()

        uic.loadUi("./sos/res/gui/fit_line.ui", self)

        # Initialisation of variables
        self.inter = False

        self.fit_method('gauss')
        
        self.flag_apply = False
        
        self.iter_points = []
        self.lines_method = []
        self.areas = []

        self.fit_data_dict = {}

        # Assign buttons
        self.cancelButton.mousePressEvent = self.close_widget
        # Activate interaction
        self.interactButton.mousePressEvent = self.activate_interactive
        # Automatic button
        self.autoButton.mousePressEvent = self.auto_fitting
        # Function to fit
        self.gaussButton.mousePressEvent = self.gauss_line_active
        self.lorentzButton.mousePressEvent = self.lorentz_line_active
        # Clear button
        self.clearButton.mousePressEvent = self.reset_canvas
        # Apply button
        self.applyButton.mousePressEvent = self.apply_line_fitting
        # Accept button
        self.acceptButton.mousePressEvent = self.accept_line_fitting

        #cid = self.f1.figure.canvas.mpl_connect('button_press_event', self.resDraw)


    def load_init_params(self, fig, ax, x, y, kind, units, name, save):
        # Load initial params
        self.x = np.array(x)
        self.y = np.array(y)

        # File name
        self.nameLabel.setText(name)

        # Initialise corrected data array
        self.data_corrected = self.y.copy()

        # Load units
        self.ux = units[0]
        self.uy = units[1]

        # Data type
        self.kind = kind

        # Get figure
        self.fig = fig
        self.fig.subplots_adjust(left=0.1, bottom=0.125, right=0.98,
                                top=0.98, wspace=None, hspace=None)
        self.ax = ax 
        
        # Save figure?
        self.save = save

        # Update plot
        self._addmpl(self.fig)

        # Initial plot
        self.initial_plot()


    def load_fit_params(self, popt, pcov, peaks, lines_method):
        # Load predefined data

        self.lines_method = lines_method

        # Update params
        self._update_fit_params(popt, pcov)
        # Update plot
        self._update_plot(popt, lines_method)
        # Fill the table
        self._fill_table()

        self.iter_points = list(peaks)

        # Enable selection button
        self.inter = True
        self._change_inter_icon()

        self.flag_apply = True


    def _update_fit_params(self, popt, pcov):
        # New parameters
        self.A = popt[::3]
        self.mu = popt[1::3]
        self.sigma = popt[2::3]
        self.pcov = pcov


    def _update_peaks_points(self):
        # Update peaks points 

        new_iter_points = []
        for i in self.mu:
            idx_array = np.where(i<=self.x)
            if len(idx_array) > 0:
                idx = idx_array[0][0]
            else:
                idx = len(self.x)-1
            new_iter_points.append(idx)

        self.iter_points = new_iter_points


    def _fill_table(self):

        # Create table
        data = np.zeros((4, len(self.A)))
        data[0] = self.A
        data[1] = self.mu
        data[2] = self.sigma
        data[3] = self.areas

        self._createLineTable(data.T, self.lines_method, self.colors)


    def close_widget(self, event):
        # Disbale graphs
        self.close()


    def auto_fitting(self, event):
        # Get distance and heigth limits for peaks
        dist = float(self.d_peaks_Edit.text())
        height_div = float(self.p_peaks_Edit.text())
        # Get peaks and widths
        peaks = find_profile_peaks(self.y, dist=dist, height_div=height_div)
        # Apply the adjust [Only Gaussian functions]
        if len(peaks) > 0:
            popt, pcov = solve_line(self.x, self.y, peaks, ['G']*len(peaks))
            if len(popt) > 0:
                self.load_fit_params(popt, pcov, peaks, ['G']*len(peaks))
        else:
            msg('No lines found', 'warn')


    def gauss_line_active(self, event):
        # Gaussian method fit
        self.fit_method('gauss')


    def lorentz_line_active(self, event):
        # Lorentzian method fit
        self.fit_method('loren')


    def activate_interactive(self, event):
        # Interactive activation
        self.inter = not self.inter
        self._change_inter_icon()


    def _change_inter_icon(self):
        # Change iteraction icon
        if self.inter:
            self._onclick_xy = self.fig.canvas.mpl_connect('button_press_event', self._onclick)

            icon_img = './sos/res/icons/int_sel.png'
        else:
            self.fig.canvas.mpl_disconnect(self._onclick_xy)

            icon_img = './sos/res/icons/int.png'

        self.interactButton.setIcon(QIcon(icon_img))


    def fit_method(self, method):

        self.method = method

        if self.method == 'gauss':
            gauss = './sos/res/icons/gaussian_sel.png'
            loren = './sos/res/icons/lorentzian.png'   

        elif self.method == 'loren':
            gauss = './sos/res/icons/gaussian.png'
            loren = './sos/res/icons/lorentzian_sel.png'  

        self.gaussButton.setIcon(QIcon(gauss))
        self.lorentzButton.setIcon(QIcon(loren))


    def apply_line_fitting(self, event):

        if not self.flag_apply:
            # Guess parameters
            guess, bounds_min, bounds_max = guess_line(self.x, self.y, self.iter_points, func=self.lines_method)

            #self._update_plot(guess, self.iter_points, self.lines_method)

            # Variable backup, in case that the fitting wouldn't be optimal
            iter_bkp = np.copy(self.iter_points)
            line_bkp = np.copy(self.lines_method)
            
            try:
                # Fit gaussian
                popt, pcov = curve_fit(lambda x, *params: multiple_lines(x, self.lines_method, *params), 
                             self.x, self.y, p0=guess, bounds=(bounds_min, bounds_max))

                # Update params
                self._update_fit_params(popt, pcov)
                # Update plot
                self._update_plot(popt, self.lines_method)
                # Fill the table
                self._fill_table()
                # Update peaks positions
                self._update_peaks_points()

                # Set flag
                self.flag_apply = True

            except Exception as e:
                self.iter_points = list(iter_bkp)
                self.lines_method = list(line_bkp)
                msg('Optimal parameters not found. Select other point\n.'+str(e), 'fail')


    def accept_line_fitting(self, event):
        # Applying baseline correction
        if not self.flag_apply:
            self.apply_line_fitting(event)
        
        self.signal_fitting.emit(self.kind)
        
        if self.save:
            now = datetime.now()
            name = now.strftime("%d-%m-%Y_%H-%M-%S")
            self.fig.savefig('fig_'+name+'_fit_lines.png')

        self.close()


    def _onclick(self, event):
        """
            On click event to select lines
        """
        if event.inaxes == self.ax:
            # Left-click
            if event.button == 1:
                ix, iy = event.xdata, event.ydata
                # Add peaks
                xarray = np.where(self.x>ix)[0]
                if len(xarray) > 0:
                    xpos = xarray[0]
                else:
                    xpos = len(self.x)-1
                self.iter_points.append(xpos)

                if self.method == 'gauss':
                    self.lines_method.append('G')
                elif self.method == 'loren':
                    self.lines_method.append('L')

                self.flag_apply = False

            # Right-click
            elif event.button == 3:
                ix, iy = event.xdata, event.ydata
                popt = []
                # Remove points
                # Define threshold
                thresh = 5*np.mean(np.diff(self.x))
                xlines = np.where((ix >= (np.array(self.x)[self.iter_points] - thresh)) & 
                               (ix < (np.array(self.x)[self.iter_points] + thresh)))[0]

                try:
                    if len(xlines) > 0:
                        x_min = np.argmin(np.abs((np.array(self.x)[np.array(self.iter_points)[xlines]] - ix)))
                        self.iter_points.remove(self.iter_points[xlines[x_min]])

                        if self.method == 'gauss':
                            self.lines_method.remove(self.lines_method[xlines[x_min]])
                        elif self.method == 'loren':
                            self.lines_method.remove(self.lines_method[xlines[x_min]])

                        self.flag_apply = False

                except Exception as e:
                    pass

            # Update plot
            self._update_selected_lines(self.iter_points)


    def _update_selected_lines(self, points):
        """
            Update selected Points/Region in the canvas
        """
        self.ax.clear()

        # Label axes
        ux_label = self.ux
        if self.ux:
            ux_label = '['+ux_label+']'
        uy_label = self.uy
        if self.uy:
            uy_label = '['+uy_label+']'

        self.ax.set_xlabel(r''+ux_label)
        self.ax.set_ylabel(r'Temperature '+uy_label)

        maxSpec = np.max(self.y)/20.
        self.ax.plot(self.x, self.y, 'k')

        for i in range(len(points)):
            if self.lines_method[i] == 'G':
                color = 'b'
            elif self.lines_method[i] == 'L':
                color = 'r'

            self.ax.plot([self.x[points[i]], self.x[points[i]]], [self.y[points[i]]+maxSpec, self.y[points[i]]+2*maxSpec], color, linewidth=1)
            #self.ax.plot(self.x[points[i]], self.y[points[i]], 'r+')

        self.ax.grid()

        self.fig.canvas.draw_idle()


    def _update_plot(self, params, funcs):
        """
            Update canvas
        """
        self.ax.clear()

        # Label axes
        ux_label = self.ux
        if self.ux:
            ux_label = '['+ux_label+']'
        uy_label = self.uy
        if self.uy:
            uy_label = '['+uy_label+']'

        self.ax.set_xlabel(r''+ux_label)
        self.ax.set_ylabel(r'Temperature '+uy_label)

        # Get colors
        n = len(funcs)
        clrs = cm.get_cmap('viridis', n)
        maxSpec = np.max(self.y)/20.

        ck = []
        areas = []
        for i in range(n):
            param = params[3*i:3*(i+1)]
            curve, area = get_fit_curve(self.x, funcs[i], param)
            areas.append(area)

            ck.append(colors.to_hex(clrs(i)))

            self.ax.plot(self.x, curve, '-.', linewidth=0.75, color=clrs(i), label=r'Line '+str(i+1))
            self.ax.plot([param[1], param[1]], [maxSpec+param[0], 3*maxSpec+param[0]], color=clrs(i), linewidth=1.5)
            self.ax.text(param[1], 3.2*maxSpec+param[0], 'L'+str(i+1))

        self.colors = ck
        self.areas = areas

        self.ax.plot(self.x, self.y, 'k')
        self.ax.plot(self.x, multiple_lines(self.x, funcs, *params), 'r')

        self.ax.grid()
        #self.ax.legend()

        self.fig.canvas.draw_idle()


    def _createLineTable(self, data, line_funcs, color):
        
        # Create header for the data
        # Label axes
        ux_label = self.ux
        if self.ux:
            ux_label = '['+ux_label+']'
        uy_label = self.uy
        if self.uy:
            uy_label = '['+uy_label+']'

        header = ["Line", "Color", "Amplitude"+uy_label, 
                  "Position"+ux_label, "Width"+ux_label, 
                  "Area", "Method"]
        
        data_table = []
        data_table.append(header)

        self.fit_data_dict = {}
        for i in range(len(data)):
            # Kind of line
            if line_funcs[i] == 'G':
                line_label = "Gaussian"
            elif line_funcs[i] == 'L':
                line_label = "Lorentzian"
            # Extra parameters
            row = ['L'+str(i+1), color[i], data[i][0], data[i][1], 
                    np.abs(data[i][2]), data[i][3], line_label]
            row_data_to_write = [data[i][0], data[i][1], 
                    np.abs(data[i][2]), data[i][3], line_label]

            self.fit_data_dict['L'+str(i+1)] = row_data_to_write

            data_table.append(row)
        
        # Build the table for Historial Panel
        self.model = TableLinesModel(data_table)
        self.tableFitLines.setModel(self.model)

        self.tableFitLines.horizontalHeader().hide()
        self.tableFitLines.verticalHeader().hide()
        

    def initial_plot(self):
        """
            Initial plot
        """
        self.ax.clear()
        self.ax.plot(self.x, self.y, 'k')

        # Label axes
        ux_label = self.ux
        if self.ux:
            ux_label = '['+ux_label+']'
        uy_label = self.uy
        if self.uy:
            uy_label = '['+uy_label+']'

        self.ax.set_xlabel(r''+ux_label)
        self.ax.set_ylabel(r'Temperature '+uy_label)

        self.ax.grid()

        self.fig.canvas.draw_idle()


    def reset_canvas(self, event):
        # Restart to initial plot
        self.initial_plot()

        self.iter_points = []
        self.lines_method = []
        self.areas = []
        self.colors = []


    def _addmpl(self, fig):
        
        self.canvas = FigureCanvas(fig)
        self.plotLayout.addWidget(self.canvas)
        self.canvas.draw()
        self.toolbar = NavigationToolbar(self.canvas,
           self, coordinates=True)
        self.plotLayout.addWidget(self.toolbar)


    def _rmmpl(self):
        self.plotLayout.removeWidget(self.canvas)
        self.canvas.close()
        self.plotLayout.removeWidget(self.toolbar)
        self.toolbar.close()


class HeaderEditWindow(QWidget):
    """
        Interactive window to fit spectral lines
        Parameters
        ----------

        ----------
    """
    # Signal to update data
    signal_fitting = QtCore.pyqtSignal(str)

    def __init__(self):

        super(HeaderEditWindow, self).__init__()

        uic.loadUi("./sos/res/gui/header.ui", self)

        # Add button
        self.addButton.mousePressEvent = self.addKey
        # Edit button
        self.editButton.mousePressEvent = self.editKey
        # Delete button
        self.deleteButton.mousePressEvent = self.deleteKey
        # Apply button
        self.applyButton.mousePressEvent = self.apply
        # Cancel button
        self.cancelButton.mousePressEvent = self.cancel
        # Accept button
        self.acceptButton.mousePressEvent = self.accept

        #cid = self.f1.figure.canvas.mpl_connect('button_press_event', self.resDraw)


    def load_init_params(self, path):
        # Load initial params
        self.x = np.array(x)
        self.y = np.array(y)

        # File name
        self.nameLabel.setText(name)

        # Initialise corrected data array
        self.data_corrected = self.y.copy()

        # Load units
        self.ux = units[0]
        self.uy = units[1]

        # Data type
        self.kind = kind

        # Get figure
        self.fig = fig
        self.fig.subplots_adjust(left=0.1, bottom=0.125, right=0.98,
                                top=0.98, wspace=None, hspace=None)
        self.ax = ax 
        
        # Save figure?
        self.save = save

        # Update plot
        self._addmpl(self.fig)

        # Initial plot
        self.initial_plot()


    def load_fit_params(self, popt, pcov, peaks, lines_method):
        # Load predefined data

        self.lines_method = lines_method

        # Update params
        self._update_fit_params(popt, pcov)
        # Update plot
        self._update_plot(popt, lines_method)
        # Fill the table
        self._fill_table()

        self.iter_points = list(peaks)

        # Enable selection button
        self.inter = True
        self._change_inter_icon()

        self.flag_apply = True


    def _update_fit_params(self, popt, pcov):
        # New parameters
        self.A = popt[::3]
        self.mu = popt[1::3]
        self.sigma = popt[2::3]
        self.pcov = pcov


    def _update_peaks_points(self):
        # Update peaks points 

        new_iter_points = []
        for i in self.mu:
            idx_array = np.where(i<=self.x)
            if len(idx_array) > 0:
                idx = idx_array[0][0]
            else:
                idx = len(self.x)-1
            new_iter_points.append(idx)

        self.iter_points = new_iter_points


    def _fill_table(self):

        # Create table
        data = np.zeros((4, len(self.A)))
        data[0] = self.A
        data[1] = self.mu
        data[2] = self.sigma
        data[3] = self.areas

        self._createLineTable(data.T, self.lines_method, self.colors)


    def close_widget(self, event):
        # Disbale graphs
        self.close()


    def auto_fitting(self, event):
        # Get distance and heigth limits for peaks
        dist = float(self.d_peaks_Edit.text())
        height_div = float(self.p_peaks_Edit.text())
        # Get peaks and widths
        peaks = find_profile_peaks(self.y, dist=dist, height_div=height_div)
        # Apply the adjust [Only Gaussian functions]
        if len(peaks) > 0:
            popt, pcov = solve_line(self.x, self.y, peaks, ['G']*len(peaks))
            if len(popt) > 0:
                self.load_fit_params(popt, pcov, peaks, ['G']*len(peaks))
        else:
            msg('No lines found', 'warn')


    def gauss_line_active(self, event):
        # Gaussian method fit
        self.fit_method('gauss')


    def lorentz_line_active(self, event):
        # Lorentzian method fit
        self.fit_method('loren')


    def activate_interactive(self, event):
        # Interactive activation
        self.inter = not self.inter
        self._change_inter_icon()


    def _change_inter_icon(self):
        # Change iteraction icon
        if self.inter:
            self._onclick_xy = self.fig.canvas.mpl_connect('button_press_event', self._onclick)

            icon_img = './sos/res/icons/int_sel.png'
        else:
            self.fig.canvas.mpl_disconnect(self._onclick_xy)

            icon_img = './sos/res/icons/int.png'

        self.interactButton.setIcon(QIcon(icon_img))


    def fit_method(self, method):

        self.method = method

        if self.method == 'gauss':
            gauss = './sos/res/icons/gaussian_sel.png'
            loren = './sos/res/icons/lorentzian.png'   

        elif self.method == 'loren':
            gauss = './sos/res/icons/gaussian.png'
            loren = './sos/res/icons/lorentzian_sel.png'  

        self.gaussButton.setIcon(QIcon(gauss))
        self.lorentzButton.setIcon(QIcon(loren))


    def apply_line_fitting(self, event):

        if not self.flag_apply:
            # Guess parameters
            guess, bounds_min, bounds_max = guess_line(self.x, self.y, self.iter_points, func=self.lines_method)

            #self._update_plot(guess, self.iter_points, self.lines_method)

            # Variable backup, in case that the fitting wouldn't be optimal
            iter_bkp = np.copy(self.iter_points)
            line_bkp = np.copy(self.lines_method)
            
            try:
                # Fit gaussian
                popt, pcov = curve_fit(lambda x, *params: multiple_lines(x, self.lines_method, *params), 
                             self.x, self.y, p0=guess, bounds=(bounds_min, bounds_max))

                # Update params
                self._update_fit_params(popt, pcov)
                # Update plot
                self._update_plot(popt, self.lines_method)
                # Fill the table
                self._fill_table()
                # Update peaks positions
                self._update_peaks_points()

                # Set flag
                self.flag_apply = True

            except Exception as e:
                self.iter_points = list(iter_bkp)
                self.lines_method = list(line_bkp)
                msg('Optimal parameters not found. Select other point\n.'+str(e), 'fail')


    def accept_line_fitting(self, event):
        # Applying baseline correction
        if not self.flag_apply:
            self.apply_line_fitting(event)
        
        self.signal_fitting.emit(self.kind)
        
        if self.save:
            now = datetime.now()
            name = now.strftime("%d-%m-%Y_%H-%M-%S")
            self.fig.savefig('fig_'+name+'_fit_lines.png')

        self.close()


    def _onclick(self, event):
        """
            On click event to select lines
        """
        if event.inaxes == self.ax:
            # Left-click
            if event.button == 1:
                ix, iy = event.xdata, event.ydata
                # Add peaks
                xarray = np.where(self.x>ix)[0]
                if len(xarray) > 0:
                    xpos = xarray[0]
                else:
                    xpos = len(self.x)-1
                self.iter_points.append(xpos)

                if self.method == 'gauss':
                    self.lines_method.append('G')
                elif self.method == 'loren':
                    self.lines_method.append('L')

                self.flag_apply = False

            # Right-click
            elif event.button == 3:
                ix, iy = event.xdata, event.ydata
                popt = []
                # Remove points
                # Define threshold
                thresh = 5*np.mean(np.diff(self.x))
                xlines = np.where((ix >= (np.array(self.x)[self.iter_points] - thresh)) & 
                               (ix < (np.array(self.x)[self.iter_points] + thresh)))[0]

                try:
                    if len(xlines) > 0:
                        x_min = np.argmin(np.abs((np.array(self.x)[np.array(self.iter_points)[xlines]] - ix)))
                        self.iter_points.remove(self.iter_points[xlines[x_min]])

                        if self.method == 'gauss':
                            self.lines_method.remove(self.lines_method[xlines[x_min]])
                        elif self.method == 'loren':
                            self.lines_method.remove(self.lines_method[xlines[x_min]])

                        self.flag_apply = False

                except Exception as e:
                    pass

            # Update plot
            self._update_selected_lines(self.iter_points)


    def _update_selected_lines(self, points):
        """
            Update selected Points/Region in the canvas
        """
        self.ax.clear()

        # Label axes
        ux_label = self.ux
        if self.ux:
            ux_label = '['+ux_label+']'
        uy_label = self.uy
        if self.uy:
            uy_label = '['+uy_label+']'

        self.ax.set_xlabel(r''+ux_label)
        self.ax.set_ylabel(r'Temperature '+uy_label)

        maxSpec = np.max(self.y)/20.
        self.ax.plot(self.x, self.y, 'k')

        for i in range(len(points)):
            if self.lines_method[i] == 'G':
                color = 'b'
            elif self.lines_method[i] == 'L':
                color = 'r'

            self.ax.plot([self.x[points[i]], self.x[points[i]]], [self.y[points[i]]+maxSpec, self.y[points[i]]+2*maxSpec], color, linewidth=1)
            #self.ax.plot(self.x[points[i]], self.y[points[i]], 'r+')

        self.ax.grid()

        self.fig.canvas.draw_idle()


    def _update_plot(self, params, funcs):
        """
            Update canvas
        """
        self.ax.clear()

        # Label axes
        ux_label = self.ux
        if self.ux:
            ux_label = '['+ux_label+']'
        uy_label = self.uy
        if self.uy:
            uy_label = '['+uy_label+']'

        self.ax.set_xlabel(r''+ux_label)
        self.ax.set_ylabel(r'Temperature '+uy_label)

        # Get colors
        n = len(funcs)
        clrs = cm.get_cmap('viridis', n)
        maxSpec = np.max(self.y)/20.

        ck = []
        areas = []
        for i in range(n):
            param = params[3*i:3*(i+1)]
            curve, area = get_fit_curve(self.x, funcs[i], param)
            areas.append(area)

            ck.append(colors.to_hex(clrs(i)))

            self.ax.plot(self.x, curve, '-.', linewidth=0.75, color=clrs(i), label=r'Line '+str(i+1))
            self.ax.plot([param[1], param[1]], [maxSpec+param[0], 3*maxSpec+param[0]], color=clrs(i), linewidth=1.5)
            self.ax.text(param[1], 3.2*maxSpec+param[0], 'L'+str(i+1))

        self.colors = ck
        self.areas = areas

        self.ax.plot(self.x, self.y, 'k')
        self.ax.plot(self.x, multiple_lines(self.x, funcs, *params), 'r')

        self.ax.grid()
        #self.ax.legend()

        self.fig.canvas.draw_idle()


    def _createLineTable(self, data, line_funcs, color):
        
        # Create header for the data
        # Label axes
        ux_label = self.ux
        if self.ux:
            ux_label = '['+ux_label+']'
        uy_label = self.uy
        if self.uy:
            uy_label = '['+uy_label+']'

        header = ["Line", "Color", "Amplitude"+uy_label, 
                  "Position"+ux_label, "Width"+ux_label, 
                  "Area", "Method"]
        
        data_table = []
        data_table.append(header)

        self.fit_data_dict = {}
        for i in range(len(data)):
            # Kind of line
            if line_funcs[i] == 'G':
                line_label = "Gaussian"
            elif line_funcs[i] == 'L':
                line_label = "Lorentzian"
            # Extra parameters
            row = ['L'+str(i+1), color[i], data[i][0], data[i][1], 
                    np.abs(data[i][2]), data[i][3], line_label]
            row_data_to_write = [data[i][0], data[i][1], 
                    np.abs(data[i][2]), data[i][3], line_label]

            self.fit_data_dict['L'+str(i+1)] = row_data_to_write

            data_table.append(row)
        
        # Build the table for Historial Panel
        self.model = TableLinesModel(data_table)
        self.tableFitLines.setModel(self.model)

        self.tableFitLines.horizontalHeader().hide()
        self.tableFitLines.verticalHeader().hide()
        

    def initial_plot(self):
        """
            Initial plot
        """
        self.ax.clear()
        self.ax.plot(self.x, self.y, 'k')

        # Label axes
        ux_label = self.ux
        if self.ux:
            ux_label = '['+ux_label+']'
        uy_label = self.uy
        if self.uy:
            uy_label = '['+uy_label+']'

        self.ax.set_xlabel(r''+ux_label)
        self.ax.set_ylabel(r'Temperature '+uy_label)

        self.ax.grid()

        self.fig.canvas.draw_idle()


    def reset_canvas(self, event):
        # Restart to initial plot
        self.initial_plot()

        self.iter_points = []
        self.lines_method = []
        self.areas = []
        self.colors = []


    def _addmpl(self, fig):
        
        self.canvas = FigureCanvas(fig)
        self.plotLayout.addWidget(self.canvas)
        self.canvas.draw()
        self.toolbar = NavigationToolbar(self.canvas,
           self, coordinates=True)
        self.plotLayout.addWidget(self.toolbar)


    def _rmmpl(self):
        self.plotLayout.removeWidget(self.canvas)
        self.canvas.close()
        self.plotLayout.removeWidget(self.toolbar)
        self.toolbar.close()
