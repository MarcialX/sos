# -*- coding: utf-8 -*-
# --------------------------------------------------------------------------------- #
# Software de Observaciones Sintéticas S.O.S.
# Table model
#
# Marcial Becerril, @ 03 March 2021
# Latest Revision: 03 Mar 2021, 01:58 GMT-6
#
# For all kind of problems, requests of enhancements and bug reports, please
# write to me at:
#
# mbecerrilt92@gmail.com
# mbecerrilt@inaoep.mx
#
# --------------------------------------------------------------------------------- #

import numpy as np

from PyQt5 import QtCore, QtWidgets, uic, QtGui
from PyQt5.QtCore import Qt, QObject
from PyQt5.QtWidgets import QApplication, QWidget


# Table model for Line fitting tool
# =========================================================
class TableLinesModel(QtCore.QAbstractTableModel):

    def __init__(self, data):
        super(TableLinesModel, self).__init__()

        self._data = data

    def data(self, index, role):
        if role == Qt.DisplayRole:
            if index.row() > 0 and index.column() == 1:
                pass
            else:
                value = self._data[index.row()][index.column()]
                if index.row() > 0 and index.column() > 1 and index.column() != 6:
                    value = str(np.round(value,2)) 
                return value

        if role == Qt.TextAlignmentRole:
            value = self._data[index.row()][index.column()]
            if isinstance(value, int) or isinstance(value, float):
                return Qt.AlignVCenter + Qt.AlignRight
            else:
                return Qt.AlignVCenter + Qt.AlignHCenter

        if role == Qt.FontRole:
            if index.row() < 1:
                font = QtGui.QFont()
                font.setBold(True)
                return font

        if role == Qt.BackgroundRole:
            if index.row() < 1:
                return QtGui.QColor("#4CE3E0")
            elif index.row() > 0 and index.column() == 1:
                value = self._data[index.row()][index.column()]
                return QtGui.QColor(value)


    def rowCount(self, index):
        # The length of the outer list.
        return len(self._data)


    def columnCount(self, index):
        # The following takes the first sub-list, and returns
        # the length (only works if all rows are an equal length)
        return len(self._data[0])