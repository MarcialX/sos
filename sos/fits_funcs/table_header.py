# -*- coding: utf-8 -*-
# --------------------------------------------------------------------------------- #
# Software de Observaciones Sintéticas S.O.S.
# Header Table Model
#
# Marcial Becerril, @ 24 August 2020
# Latest Revision: 8 Sep 2021, 01:08 GMT-6
#
# For all kind of problems, requests of enhancements and bug reports, please
# write to me at:
#
# mbecerrilt92@gmail.com
# mbecerrilt@inaoep.mx
#
# --------------------------------------------------------------------------------- #


from PyQt5 import QtCore, QtWidgets, QtGui
from PyQt5.QtCore import Qt


class TableModelHeader(QtCore.QAbstractTableModel):

    def __init__(self, data):
        super(TableModelHeader, self).__init__()
        self._data = data

    def data(self, index, role):
        if role == Qt.DisplayRole:
            return self._data[index.row()][index.column()]

        if role == Qt.FontRole:
            if index.column() == 0:
                font = QtGui.QFont()
                font.setBold(True)
                return font

        if role == Qt.BackgroundRole:
            if index.column() == 0:
                return QtGui.QColor("#BBBBBB")


    def flags(self, index): 
        if not index.isValid():
            return Qt.ItemIsEnabled
        return Qt.ItemFlags(QtCore.QAbstractTableModel.flags(self, index) | Qt.ItemIsEditable)


    def rowCount(self, index):
        # The length of the outer list.
        return len(self._data)


    def columnCount(self, index):
        # The following takes the first sub-list, and returns
        # the length (only works if all rows are an equal length)
        return len(self._data[0])


    def removeRows(self, init_row, nrows=1):
        self.beginRemoveRows(QtCore.QModelIndex(), init_row, init_row + nrows - 1)       
        self._data = self._data[:init_row] + self._data[init_row + nrows:]
        self.endRemoveRows()


    def insertRow(self, row):
        self.beginInsertRows(QtCore.QModelIndex(), row, row)
        self._data.insert(row+1,  ["[NEW KEY]", None])
        self.endInsertRows()
        return True


    def setData(self, index, value, role):
        if not index.isValid():
            return False
        if role != QtCore.Qt.EditRole:
            return False
        row = index.row()
        if row < 0 or row >= len(self._data):
            return False
        column = index.column()
        if column < 0 or column >= len(self._data[0]):
            return False
        self._data[row][column] = value
        self.dataChanged.emit(index, index)
        return True


class HeaderDelegate(QtWidgets.QItemDelegate):

    def createEditor(self, parent, option, index):
        if index.column() <= 1:
            return super(HeaderDelegate, self).createEditor(parent, option, index)
        return None

    def setEditorData(self, editor, index):
        if index.column() <= 1:
            # Gets display text if edit data hasn't been set.
            text = index.data(Qt.EditRole) or index.data(Qt.DisplayRole)
            editor.setText(text)

