# -*- coding: utf-8 -*-
# --------------------------------------------------------------------------------- #
# Software de Observaciones Sintéticas S.O.S.
# Header edit
#
# Marcial Becerril, @ 24 August 2020
# Latest Revision: 6 Sep 2021, 19:00 GMT-6
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

from .table_header import TableModelHeader, HeaderDelegate


class HeaderEditWindow(QWidget):
    """
        Interactive window to fit spectral lines
        Parameters
        ----------

        ----------
    """
    # Signal to update data
    save = QtCore.pyqtSignal(bool)

    def __init__(self):

        super(HeaderEditWindow, self).__init__()

        uic.loadUi("./sos/res/gui/header.ui", self)

        # Add button
        self.addButton.mousePressEvent = self.addKey
        # Delete button
        self.deleteButton.mousePressEvent = self.deleteKey
        # Cancel button
        self.cancelButton.mousePressEvent = self.cancel
        # Accept button
        self.acceptButton.mousePressEvent = self.accept


    def _header2list(self, header, maxlen=True):

        len_list = []
        header_list = []
        for key in header:
            value = header[key]
            header_list.append([key, str(value)])
            len_list.append(len(str(value)))

        if maxlen:
            return header_list, max(len_list)

        return header_list


    def load_fits_file(self, data, header, name):
        # Load initial params
        self.data = data
        self.header = header

        # Update FITS name
        self.nameEdit.setText(name)
        if len(name) > 33:
            name = name[:30]+'...'
        self.nameLabel.setText(name)

        # Create and customise Table
        header_list, maxlen = self._header2list(header)

        header_delegate = HeaderDelegate()
        self.tableHeader.setItemDelegate(header_delegate)

        self.model = TableModelHeader(header_list)
        self.tableHeader.setModel(self.model)

        # Set font size
        font = QtGui.QFont("Courier New", 11)
        self.tableHeader.setFont(font)

        # Hide table headers
        self.tableHeader.horizontalHeader().hide()
        self.tableHeader.verticalHeader().hide()
        
        # Set table dimensions
        self.tableHeader.setColumnWidth(0, 90)   
        self.tableHeader.setColumnWidth(1, 180)
        #self.tableHeader.setColumnWidth(1, 3.75*maxlen/2)


    def addKey(self, event):
        # Get position of selected row
        x, y = self._get_pos_selected()

        # Remove the row
        if len(x)>0:
            min_row = min(x)    
        else:
            min_row = self.tableHeader.model().rowCount(0)
        self.tableHeader.model().insertRow(min_row)


    def deleteKey(self, event):
        # Get position of selected row
        x, y = self._get_pos_selected()

        # Remove the row
        if len(x)>0:
            max_row = max(x)
            min_row = min(x)       
            self.tableHeader.model().removeRows(min_row, nrows=max_row-min_row+1)


    def _get_pos_selected(self):
        #print(self.tableHeader.selectionModel().selectedRows())
        indexes = self.tableHeader.selectionModel().selectedIndexes()

        x, y = [], []
        for index in indexes:
            x.append(index.row())
            y.append(index.column())

        return x, y


    def accept(self, event):
        # If changes has been accepted
        # Read the table
        data = {}
        for row in range(self.tableHeader.model().rowCount(0)):
            cell = str(self.tableHeader.model().index(row, 0).data())
            data[cell] = {}
            data[cell] = str(self.tableHeader.model().index(row, 1).data())

        # Get name FITS file
        name = self.nameEdit.text()

        # If some keys were modified or deleted
        hd_list = []
        for key in self.header:
            hd_list.append(key)
            if key in data.keys():
                if data[key] != str(self.header[key]):
                    try:
                        if isinstance(self.header[key], int):
                            data[key] = int(data[key])
                        elif isinstance(self.header[key], float):
                            data[key] = float(data[key])
                        elif isinstance(self.header[key], bool):
                            if data[key] == 'True':
                                data[key] = True
                            else:
                                data[key] = False
                            data[key] = bool(data[key])
                        else:
                            data[key] = str(data[key])
                    except:
                        data[key] = str(data[key])

                    self.header[key] = data[key]
            else:
                del self.header[key]
        
        # If some keys were added
        for key in data.keys():
            if not key in hd_list:
                if data[key].isnumeric():
                    data[key] = float(data[key])
                elif data[key] == 'True':
                    data[key] = True
                elif data[key] == 'False':
                    data[key] = False
                self.header[key] = (data[key], 'New parameter')

        # Save FITS file
        # Overwrite the FITS file
        name = name.split('.fits')[0]
        os.system('rm -rf '+name)
        hdu = fits.PrimaryHDU(self.data, header=self.header)
        hdu.writeto(name+'.fits')

        self.close()


    def cancel(self, event):
        # Disable graphs
        self.close()
