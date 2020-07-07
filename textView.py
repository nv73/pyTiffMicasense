# -*- coding: utf-8 -*-
"""
Created on Fri May 15 14:50:30 2020

@author: nick.viner
"""
from PyQt5 import QtWidgets

import textView_ui


class textView_Form(QtWidgets.QWidget, textView_ui.Ui_viewText):

    def __init__(self, parent=None):
        super().__init__()

        self.setupUi(self)

        self.actionClose.clicked.connect(self.close)

    def close(self):
        self.close()
