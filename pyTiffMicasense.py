# -*- coding: utf-8 -*-
"""
Created on Fri May 15 13:55:38 2020

@author: nick.viner
"""
import pyTiff
from PyQt5 import QtCore, QtWidgets, QtGui
import pyTiffMicasense_ui
import textView
import matplotlib.pyplot as plt
from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvasQt
from matplotlib.figure import Figure
import sys
from PIL import Image, ImageQt
from geoCanvas import geoCanvas
import gdal
import numpy as np
import osr
from datetime import datetime
import os
import tkinter
import FileDialog

class pyTiffMicasense_Form(QtWidgets.QMainWindow, pyTiffMicasense_ui.Ui_MainWindow):
    
    def __init__(self, parent=None):
        
        super().__init__()
        
        self.setupUi(self)
        
        self.logFileName = "%s_log.txt" % str(datetime.now())[:10]
        
        self.log("--------------------------------------------------------------")
        
        self.log("Running pyTiff: Micasense on %s" % os.environ['COMPUTERNAME'])
        
        self.log("Written by Nick Viner (Hakai Institute, 2020)")
        
        self.log("Contact nick.viner@hakai.org for questions or assistance")
        
        self.log("Current Date / Time: %s" % str(datetime.now()))
        
        #Connecting UI objects to their corresponding functions
        self.actionImport.triggered.connect(self.importImage)
        
        self.actionAdd_Band.triggered.connect(self.addBand)
        
        self.actionMerge_Bands_to_Raster.triggered.connect(self.mergeBands)
        
        self.actionRadiance_to_Reflectance.triggered.connect(self.radianceToReflectance)
        
        self.actionViewGeodetics.clicked.connect(self.showGeodetics)
        
        self.actionViewCameraInfo.clicked.connect(self.showCameraInfo)
        
        self.bandSelectBox.activated.connect(self.updateBandView)
        
        self.actionImport_Calibration_Panel.triggered.connect(self.importCalibrationPanel)
    
        #Container for image. Using dictionary allows for multiple bands to be added
        self.image = {}
                
        ### Matplotlib images of the input data
        self.inputImage = None
        self.inputCalibration = None
        
        #Add plot of raw image to UI
        self.imageFigure = Figure()
        self.imageFigureCanvas = FigureCanvasQt(self.imageFigure)
        self.graphicsView.addWidget(self.imageFigureCanvas)
        self.imageAx = None
        
        #Add plot of panel image to UI
        self.panelFigure = Figure()
        self.panelFigureCanvas = FigureCanvasQt(self.panelFigure)
        self.calibrationGraphicsView.addWidget(self.panelFigureCanvas)
        self.panelAx = None
        
        ####
        
        self.reflectance = {}
        
        self.calibrationPanelPath = None
        
        self.resolution = 0.0839
        
        self.bandCount = 0
        
        self.conversionValue = 0
        
        self.checkForReflectance = False
        
    ### FILE DROP DOWN ###
    
    def importImage(self):
        
        if self.bandCount != 0:
            
            self.log("Image already loaded.")
            
            return
        
        #Get filepath for image file
        f = QtWidgets.QFileDialog.getOpenFileName(
            self, 'Select image to open...', '', ".tif(*.tif)")
        
        #Get filepath for calibration panel file
        c = QtWidgets.QFileDialog.getOpenFileName(
            self, 'Select calibration panel image...', '', ".tif(*.tif)")
        
        self.calibrationPanelPath = c[0]
        
        #Load files to memory
        self.inputImage = plt.imread(f[0])
        self.inputCalibration = plt.imread(c[0])
        
        # Get user input resolution value
        value, ok = QtWidgets.QInputDialog().getText(self, "Select Resolution",
                                     "Resolution (can be changed later):", QtWidgets.QLineEdit.Normal)
        
        #Check if value was entered, if not, set to 1
        if ok and value:
            self.resolution = float(value)
        else:
            self.resolution = 1
        
        #Generate pyTiff object and assign to dictionary
        self.image[1] = pyTiff.pyTiff(filename=f[0], calibPanel=c[0], micasense=True, overrideResolution=self.resolution)
        
        self.bandCount += 1
                        
        self.log("Imported image %s" % f[0])
        
        self.log("Imported panel %s" % c[0])
        
        ### Update the data displays
        
        #Calibration Panel
        self.panelAx = self.panelFigure.add_subplot(111)
        self.panelAx.imshow(self.inputCalibration)
        self.panelFigureCanvas.draw()
        
        #Input Image
        self.imageAx = self.imageFigure.add_subplot(111)
        self.imageAx.imshow(self.inputImage)
        self.imageFigureCanvas.draw()
        
        #Run the panel calibration
        self.conversionValue = self.image[1].get_calibration()
        
        self.fileGeoLocationX.setText(str(self.image[1].x))
        self.fileGeoLocationY.setText(str(self.image[1].y))
        self.fileGeoLocationLat.setText(str(self.image[1].lat))
        self.fileGeoLocationLon.setText(str(self.image[1].lon))
        self.filePathText.setText(str(self.image[1].filename))
        self.calPathText.setText(str(self.image[1].calibrationImagePath))
        self.redValText.setText(str(self.image[1].calibrationValues["Red"]))
        self.blueValText.setText(str(self.image[1].calibrationValues["Blue"]))
        self.greenValText.setText(str(self.image[1].calibrationValues["Green"]))
        self.redEdgeValText.setText(str(self.image[1].calibrationValues["Red edge"]))
        self.nearInfraredValText.setText(str(self.image[1].calibrationValues["NIR"]))
        self.panelRadianceText.setText(str(self.image[1].calibrationValues["meanPanelRadiance"]))
        self.panelCornerText.setText(str(self.image[1].panel))
        
    def importCalibrationPanel(self):
        
        c = QtWidgets.QFileDialog.getOpenFileName(
            self, 'Select panel image to open...', '', ".tif(*.tif)")
        
        self.calibrationPanelPath = c[0]
        
    ### PROCESS DROP DOWN ###
    
    def addBand(self):
        
        f = QtWidgets.QFileDialog.getOpenFileName(
            self, 'Select image to open...', '', ".tif(*.tif)")
        
        self.bandCount += 1
        
        self.image[self.bandCount] = pyTiff.pyTiff(filename=f[0], calibPanel=self.calibrationPanelPath,
                                                   micasense=True, overrideResolution=self.resolution)
        
        self.image[self.bandCount].get_calibration()
        
        self.bandSelectBox.addItem("Band %i" % self.bandCount)
        
        self.log("Added image %s to band %i" % (f[0], self.bandCount))
        
    def mergeBands(self):
                                 
        fileList = []
        
        for x in range(1, self.bandCount + 1):
            
            fileList.append(self.image[x].tiff.ReadAsArray())
        
        outRaster = gdal.GetDriverByName("GTiff").Create("Test.tif", self.image[1].reflectanceImage.shape[1], 
                                                         self.image[1].reflectanceImage.shape[0], self.bandCount, gdal.GDT_Float32)
        
        for x in range(0, self.bandCount):
            
            outRaster.GetRasterBand(x + 1).WriteArray(fileList[x])
        
        outRaster.SetGeoTransform(self.image[1].geoTransform)
        
        srs = osr.SpatialReference()
        
        srs.ImportFromEPSG(self.image[1].epsg)
        
        outRaster.SetProjection(srs.ExportToWkt())
        
    def radianceToReflectance(self):
        
        print('radiance to reflectance')
        
        saveImage = QtWidgets.QMessageBox.question(self, 'Save images to file?', 'Save image to file?', QtWidgets.QMessageBox.Yes, QtWidgets.QMessageBox.No)
        
        if self.bandCount >= 1:
            
            if saveImage == QtWidgets.QMessageBox.No:
                
                for x in range(1, self.bandCount + 1):
                
                    self.reflectance[x] = self.image[x].get_reflectance(write_to_file=False)
                
            else:
                
                for x in range(1, self.bandCount + 1):
                    
                    self.reflectance[x] = self.image[x].get_reflectance(write_to_file=True)
            
            self.imageAx.clear()
            
            self.imageAx.imshow(self.reflectance[1])
            
            self.imageFigureCanvas.draw()
            
            self.checkForReflectance = True
            
    def updateBandView(self):
        
        self.imageAx.clear()
        
        band = self.bandSelectBox.currentText()[-1]
        
        if self.checkForReflectance == True:
        
            self.imageAx.imshow(self.reflectance[int(band)])
            
        else:
            
            self.imageAx.imshow(plt.imread(self.image[int(band)].filename))
        
        self.imageFigureCanvas.draw()
                   
    ### MAIN TAB ###
    
    def showGeodetics(self):
        
        if len(self.image) == 0:
            
            self.log("No file loaded")
            
        else:
            
            self.log("\n---------GEODETICS---------\n")
            
            ypr = self.image[1].meta.dls_pose()
            pos = self.image[1].geoTransform
            
            self.log("yaw, pitch, roll: %f, %f, %f" % (ypr[0], ypr[1], ypr[2]))
            self.log("x, y: %f, %f" % (pos[0], pos[3]))
            self.log("lat, lon: %f, %f" % (self.image[1].lat, self.image[1].lon))
            self.log("EPSG Code: %s" % self.image[1].epsg)
        
    def showCameraInfo(self):
        
        if len(self.image) == 0:
            
            self.log("No file loaded")
            
        else:
        
            activeBand = int(self.bandSelectBox.currentText()[-1])
        
            xmpTags = [x for x in self.image[activeBand].meta.get_all() if "XMP" in x]
        
            exifTags = [x for x in self.image[activeBand].meta.get_all() if "EXIF" in x]
        
            self.log("\n---------CAMERA INFO---------\n")
        
            self.log("\n----XMP TAGS----\n")
        
            for x in xmpTags:
            
                self.log("%s: %s" % (x[4:], self.image[activeBand].meta.get_all()[x]))
            
            self.log("\n----EXIF TAGS----\n")
            
            for x in exifTags:
                
                if x == "EXIF:GDALMetadata":
                    
                    pass
                
                else:
            
                    self.log("%s: %s" % (x[5:], self.image[activeBand].meta.get_all()[x]))
        
    def log(self, inputText, saveLog = True):
        """Input text to be displayed in the GUI logging window"""
        
        self.logView.insertPlainText(">> %s \n" % inputText)
        
        self.logView.moveCursor(QtGui.QTextCursor.End)
        
        if saveLog == True:
            
            f = open(self.logFileName, "a")
            
            f.write(inputText + '\n')
            
            f.close()
        

app = QtWidgets.QApplication(sys.argv)
form = pyTiffMicasense_Form()
form.show()
app.exec_()