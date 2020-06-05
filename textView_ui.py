# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'textView.ui'
#
# Created by: PyQt5 UI code generator 5.9.2
#
# WARNING! All changes made in this file will be lost!

from PyQt5 import QtCore, QtGui, QtWidgets

class Ui_viewText(object):
    def setupUi(self, viewText):
        viewText.setObjectName("viewText")
        viewText.resize(492, 302)
        self.horizontalLayout = QtWidgets.QHBoxLayout(viewText)
        self.horizontalLayout.setObjectName("horizontalLayout")
        self.actionClose = QtWidgets.QPushButton(viewText)
        self.actionClose.setObjectName("actionClose")
        self.horizontalLayout.addWidget(self.actionClose)
        self.textEdit = QtWidgets.QTextEdit(viewText)
        self.textEdit.setObjectName("textEdit")
        self.horizontalLayout.addWidget(self.textEdit)

        self.retranslateUi(viewText)
        QtCore.QMetaObject.connectSlotsByName(viewText)

    def retranslateUi(self, viewText):
        _translate = QtCore.QCoreApplication.translate
        viewText.setWindowTitle(_translate("viewText", "View Text"))
        self.actionClose.setText(_translate("viewText", "Close"))


if __name__ == "__main__":
    import sys
    app = QtWidgets.QApplication(sys.argv)
    viewText = QtWidgets.QWidget()
    ui = Ui_viewText()
    ui.setupUi(viewText)
    viewText.show()
    sys.exit(app.exec_())

