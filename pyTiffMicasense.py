# -*- coding: utf-8 -*-
"""
Created on Fri May 15 13:55:38 2020

@author: nick.viner
"""
import os
import sys
from datetime import datetime

import gdal
import matplotlib.pyplot as plt
import osr
from PyQt5 import QtWidgets, QtGui
from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvasQt
from matplotlib.figure import Figure

import pyTiff
import pyTiffMicasense_ui


class pyTiffMicasense_Form(QtWidgets.QMainWindow, pyTiffMicasense_ui.Ui_MainWindow):

    def __init__(self, parent=None):

        super().__init__()

        self.setupUi(self)

        self.logFileName = "%s_log.txt" % str(datetime.now())[:10]

        self.log("--------------------------------------------------------------")

        self.log("Running pyTiff: Micasense on %s" % os.environ.get('COMPUTERNAME', 'unknown'))

        self.log("Written by Nick Viner (Hakai Institute, 2020)")

        self.log("Contact nick.viner@hakai.org for questions or assistance")

        self.log("Current Date / Time: %s" % str(datetime.now()))
        
        self.log("--------------------------------------------------------------")

        # Connecting UI objects to their corresponding functions
        self.actionImport.triggered.connect(self.importImage)
        self.actionAdd_Band.triggered.connect(self.addBand)
        self.actionMerge_Bands_to_Raster.triggered.connect(self.mergeBands)
        self.actionRadiance_to_Reflectance.triggered.connect(self.radianceToReflectance)
        self.actionViewGeodetics.clicked.connect(self.showGeodetics)
        self.actionViewCameraInfo.clicked.connect(self.showCameraInfo)
        self.bandSelectBox.activated.connect(self.updateBandView)
        self.actionImport_Calibration_Panel.triggered.connect(self.importCalibrationPanel)
        self.menuProcess.setEnabled(False)
        self.actionImport_Calibration_Panel.setEnabled(False)      
        self.actionMerge_Bands_to_Raster.setEnabled(False)        
        self.blueValText.editingFinished.connect(self.updateCalibrationValues)       
        self.greenValText.editingFinished.connect(self.updateCalibrationValues)        
        self.redValText.editingFinished.connect(self.updateCalibrationValues)        
        self.redEdgeValText.editingFinished.connect(self.updateCalibrationValues)        
        self.nearInfraredValText.editingFinished.connect(self.updateCalibrationValues)

        # Container for image. Using dictionary allows for multiple bands to be added
        self.image = {}

        ### Matplotlib images of the input data
        self.inputImage = None
        self.inputCalibration = None

        # Add plot of raw image to UI
        self.imageFigure = Figure()
        self.imageFigureCanvas = FigureCanvasQt(self.imageFigure)
        self.graphicsView.addWidget(self.imageFigureCanvas)
        self.imageAx = None

        # Add plot of panel image to UI
        self.panelFigure = Figure()
        self.panelFigureCanvas = FigureCanvasQt(self.panelFigure)
        self.calibrationGraphicsView.addWidget(self.panelFigureCanvas)
        self.panelAx = None
        
        # Zero the initial calibration panel values
        self.redValText.setText("0.0")
        self.greenValText.setText("0.0")
        self.blueValText.setText("0.0")
        self.redEdgeValText.setText("0.0")
        self.nearInfraredValText.setText("0.0")

        ####

        self.reflectance = {}

        self.calibrationPanelPath = None
        
        self.calibrationPanels = {}
        
        self.calibrationValues = {"Blue": 0.54,"Green": 0.55,"Red": 0.54,"Red edge": 0.52,"NIR": 0.48}

        self.resolution = 0.0839

        self.bandCount = 0

        self.conversionValue = 0

        self.checkForReflectance = False
        
    def updateCalibrationValues(self):

        self.calibrationValues["Red"] = float(self.redValText.text())
        self.calibrationValues["Blue"] = float(self.blueValText.text())
        self.calibrationValues["Green"] = float(self.greenValText.text())
        self.calibrationValues["Red edge"] = float(self.redEdgeValText.text())
        self.calibrationValues["NIR"] = float(self.nearInfraredValText.text())
        
        print(self.calibrationValues)

    ### FILE DROP DOWN ###

    def importImage(self):

        if self.bandCount != 0:
            self.log("Image already loaded.")

            return

        # Get filepath for image file
        f = QtWidgets.QFileDialog.getOpenFileName(
            self, 'Select image to open...', '', ".tif(*.tif)")

        # Get filepath for calibration panel file
        c = QtWidgets.QFileDialog.getOpenFileName(
            self, 'Select calibration panel image...', '', ".tif(*.tif)")

        self.calibrationPanelPath = c[0]
        
        self.calibrationPanels[1] = c[0]
        
        # Load files to memory
        self.inputImage = plt.imread(f[0])
        self.inputCalibration = plt.imread(c[0])

        # Get user input resolution value
        #value, ok = QtWidgets.QInputDialog().getText(self, "Select Resolution",
         #                                            "Resolution (can be changed later):", QtWidgets.QLineEdit.Normal)

        # Check if value was entered, if not, set to 1
        #if ok and value:
            #self.resolution = float(value)
        #else:
            #self.resolution = 1

        # Generate pyTiff object and assign to dictionary
        self.image[1] = pyTiff.pyTiff(filename=f[0], calibPanel=c[0], micasense=True)
        
        imageMetadata = self.image[1].get_micaMetaData()
        self.cameraNameText.setText(imageMetadata.get_item("EXIF:Model"))
        
        if imageMetadata.get_item("EXIF:Model") != "RedEdge-M":
            
            self.log("WARNING: CAMERA MODEL NOT REDEDGE-M, USING DEFAULT IMAGER SIZE DIMENSIONS.")
            self.log("CHECK CAMERA TAB. CAN UPDATE VALUES IF NEED BE.")
            self.log("CHECK YOUR MICASENSE CAMERA USER MANUAL FOR MORE INFO.")
        
        self.serialNumberText.setText(imageMetadata.get_item("XMP:Serial"))
        w = imageMetadata.get_item("EXIF:ImageWidth")
        h = imageMetadata.get_item("EXIF:ImageHeight")
        wh = str(w) + " , " + str(h)
        self.pixelResText.setText(wh)
        self.focalLengthText.setText(str(imageMetadata.get_item("EXIF:FocalLength")))
        self.fovText.setText(str(imageMetadata.get_item("Composite:FOV")))
        self.imageResText.setText(str(self.image[1].cellSize[0]))
       
        self.bandCount += 1

        self.log("Imported image %s" % f[0])

        self.log("Imported panel %s" % c[0])

        ### Update the data displays

        # Calibration Panel
        self.panelAx = self.panelFigure.add_subplot(111)
        self.panelAx.imshow(self.inputCalibration)
        self.panelFigureCanvas.draw()

        # Input Image
        self.imageAx = self.imageFigure.add_subplot(111)
        self.imageAx.imshow(self.inputImage)
        self.imageFigureCanvas.draw()

        # Run the panel calibration
        #self.conversionValue = self.image[1].get_calibration()

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
        self.menuProcess.setEnabled(True)
        self.actionImport_Calibration_Panel.setEnabled(False)
        
        self.log("CALIBRATION PANEL VALUES ARE UNIQUE TO YOUR PANEL")
        self.log("PLEASE ENSURE THE CALIBRATION VALUES MATCH THOSE PROVIDED BY THE MANUFACTURER")
        
    def importCalibrationPanel(self):

        c = QtWidgets.QFileDialog.getOpenFileName(
            self, 'Select panel image to open...', '', ".tif(*.tif)")

        self.calibrationPanelPath = c[0]
        
    ### PROCESS DROP DOWN ###

    def addBand(self):
        
        f = QtWidgets.QFileDialog.getOpenFileName(
            self, 'Select image to open...', '', ".tif(*.tif)")
        
        c = QtWidgets.QFileDialog.getOpenFileName(
            self, 'Select calibration panel...', '', ".tif(*.tif)")
        
        self.bandCount += 1
        self.calibrationPanels[self.bandCount] = c[0]
        self.image[self.bandCount] = pyTiff.pyTiff(filename=f[0], calibPanel=self.calibrationPanels[self.bandCount], micasense=True)
        self.bandSelectBox.addItem("Band %i" % self.bandCount)
        self.log("Added image %s to band %i" % (f[0], self.bandCount))

    def mergeBands(self):

        fileList = []

        fileName = "merged_%s.tif" % str(datetime.now())[:10]

        for x in range(1, self.bandCount + 1):
            
            fileList.append(self.image[x].tiff.ReadAsArray())

        outRaster = gdal.GetDriverByName("GTiff").Create(fileName, self.image[1].reflectanceImage.shape[1],
                                                         self.image[1].reflectanceImage.shape[0], self.bandCount,
                                                         gdal.GDT_Float32)

        for x in range(0, self.bandCount):
            
            outRaster.GetRasterBand(x + 1).WriteArray(fileList[x])

        outRaster.SetGeoTransform(self.image[1].geoTransform)

        srs = osr.SpatialReference()

        srs.ImportFromEPSG(self.image[1].epsg)

        outRaster.SetProjection(srs.ExportToWkt())

        self.log("File %s written to %s" % (fileName, str(os.getcwd())))

    def radianceToReflectance(self):
        
        checkCal = QtWidgets.QMessageBox()
        
        checkCalAnswer = checkCal.question(self,'',
                                         "Quick Reminder: Have you double-checked your Calibration Panel values? Do they match the manufacture reported values?",
                                         checkCal.Yes | checkCal.No)
        
        if checkCalAnswer == checkCal.No:
            
            return
        
        for i in self.image:
                
            self.image[i].get_calibration()
            
        self.panelRadianceText.setText(str(self.image[1].calibrationValues["meanPanelRadiance"]))
        self.panelCornerText.setText(str(self.image[1].panel)) 
        
        print('radiance to reflectance')
        saveImage = QtWidgets.QMessageBox.question(self, 'Save images to file?', 'Save image to file?',
                                                   QtWidgets.QMessageBox.Yes, QtWidgets.QMessageBox.No)

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
            
        self.actionMerge_Bands_to_Raster.setEnabled(True)

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
                    
    #WORK ON THIS
    def estimateCellSize(self, metadata):
        
        meta = metadata.get_all()
        
        width = meta.get_item("EXIF:ImageWidth")
        height = meta.get_item("EXIF:ImageHeight")
        focalLength = meta.get_item("EXIF:FocalLength")
        fov = meta.get_item("Composite:FOV")
        imgerSize = self.imagerWidthText.Text()
        
        cs = (imgerSize * height * 100) / (focalLength * width)
        
        print(cs)
        
        return(cs)
        
    def log(self, inputText, saveLog=True):
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
