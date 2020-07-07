# -*- coding: utf-8 -*-
"""
Created on Fri Feb 14 22:51:30 2020

@author: nickv
"""
import gdal
# from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as figureCanvas
# from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as navigationToolbar
from PIL import Image, ImageQt
from PyQt5 import QtCore, QtWidgets, QtGui


# Extension of the Qt QGraphicsView
# Enables streamlined handling and viewing of georeferenced, multiband imagery
class geoCanvas(QtWidgets.QGraphicsView):
    photoClicked = QtCore.pyqtSignal(QtCore.QPoint)

    def __init__(self):

        super(geoCanvas, self).__init__()

        # Initialize the scene for the graphics view
        self.scene = QtWidgets.QGraphicsScene(self)
        self.setScene(self.scene)

        self._QtImage = QtWidgets.QGraphicsPixmapItem()

        # QGraphicsView properties
        self.setVerticalScrollBarPolicy(QtCore.Qt.ScrollBarAlwaysOff)
        self.setHorizontalScrollBarPolicy(QtCore.Qt.ScrollBarAlwaysOff)
        self.setTransformationAnchor(QtWidgets.QGraphicsView.AnchorUnderMouse)
        self.setResizeAnchor(QtWidgets.QGraphicsView.AnchorUnderMouse)
        self.setMouseTracking(True)

        # some helper variables
        self._empty = True
        self._zoom = 0

        # Holder for geoImages
        self.geoImageIndex = {}

        # Coordinate values for the mouse cursor
        self.mouse_coordinates = None
        self.selected_coordinates = None

        # Graphical coordinate indicators
        self.displayed_coordinates = QtWidgets.QGraphicsTextItem()
        self.displayed_coordinates.setTransformOriginPoint(self.displayed_coordinates.boundingRect().topLeft())
        self.displayed_coordinates_font = self.displayed_coordinates.font()
        self.displayed_coordinates_font.setPointSize(2)
        self.displayed_coordinates.setFont(self.displayed_coordinates_font)
        self.displayed_coordinates_scale = 1

        # Add the initialized image and text to the graphics view
        self.scene.addItem(self._QtImage)
        self.scene.addItem(self.displayed_coordinates)

    # Taken from https://stackoverflow.com/questions/35508711/how-to-enable-pan-and-zoom-in-a-qgraphicsvie
    def setQtImage(self, pixmap=None):

        self._zoom = 0

        if pixmap and not pixmap.isNull():

            self._empty = False

            self.setDragMode(QtWidgets.QGraphicsView.ScrollHandDrag)

            self.viewport().setCursor(QtCore.Qt.CrossCursor)

            self._QtImage.setPixmap(pixmap)

        else:

            self._empty = True

            self.setDragMode(QtWidgets.QGraphicsView.NoDrag)

            self._QtImage.setPixmap(QtGui.QPixmap())

        self.fitInView()

    def hasQtImage(self):

        return not self._empty

    # Taken from https://stackoverflow.com/questions/35508711/how-to-enable-pan-and-zoom-in-a-qgraphicsview
    def fitInView(self, scale=True):

        rect = QtCore.QRectF(self._QtImage.pixmap().rect())

        if not rect.isNull():

            self.setSceneRect(rect)

            if self.hasQtImage():
                unity = self.transform().mapRect(QtCore.QRectF(0, 0, 1, 1))

                self.scale(1 / unity.width(), 1 / unity.height())

                viewrect = self.viewport().rect()

                scenerect = self.transform().mapRect(rect)

                factor = min(viewrect.width() / scenerect.width(),
                             viewrect.height() / scenerect.height())

                self.scale(factor, factor)

            self._zoom = 0

    def wheelEvent(self, event):

        if not self._empty:

            if event.angleDelta().y() > 0:

                factor = 1.25

                self._zoom += 1

                self.displayed_coordinates.setScale(self.displayed_coordinates_scale)

                self.displayed_coordinates_scale = self.displayed_coordinates_scale * 0.8

            else:

                if self.displayed_coordinates_scale < 1.0:
                    self.displayed_coordinates_scale = self.displayed_coordinates_scale * 1.25

                factor = 0.8

                self._zoom -= 1

                self.displayed_coordinates.setScale(self.displayed_coordinates_scale)

            if self._zoom > 0:

                self.scale(factor, factor)

            elif self._zoom == 0:

                self.fitInView()

            else:

                self._zoom = 0

    def enterEvent(self, event):

        self.viewport().setCursor(QtCore.Qt.CrossCursor)

        super(geoCanvas, self).enterEvent(event)

    def mousePressEvent(self, event):

        if self._QtImage.isUnderMouse():
            self.photoClicked.emit(QtCore.QPoint(event.pos()))

        selected_coordinates = self.mapToScene(event.x(), event.y())

        self.selected_coordinates = (selected_coordinates.x(), selected_coordinates.y())

        self.viewport().setCursor(QtCore.Qt.CrossCursor)

        super(geoCanvas, self).mousePressEvent(event)

    def mouseReleaseEvent(self, event):

        super(geoCanvas, self).mouseReleaseEvent(event)

        self.viewport().setCursor(QtCore.Qt.CrossCursor)

    def mouseMoveEvent(self, event):

        mouse_coords = self.mapToScene(event.x(), event.y())

        self.mouse_coordinates = (mouse_coords.x(), mouse_coords.y())

        self.displayed_coordinates.setPlainText("X: %i, Y: %i" % (mouse_coords.x(), mouse_coords.y()))

        self.displayed_coordinates.setPos(mouse_coords.x(), mouse_coords.y())

        super(geoCanvas, self).mouseMoveEvent(event)

    def toggleDragMode(self):

        if self.dragMode() == QtWidgets.QGraphicsView.ScrollHandDrag:

            self.setDragMode(QtWidgets.QGraphicsView.NoDrag)

        elif not self.QtImage.pixmap().isNull():

            self.setDragMode(QtWidgets.QGraphicsView.ScrollHandDrag)

    # Convert a PIL image to a QPixmap
    def imageToQPixmap(self, image):

        imgQt = QtGui.QImage(ImageQt.ImageQt(image))

        qPixMap = QtGui.QPixmap.fromImage(imgQt)

        return (qPixMap)

    # Display a 2d array within the viewport (no georeferencing)
    def displayArray(self, arr):

        if len(arr.shape) != 2:

            print("Input must be 2-dimensional array.")

        else:

            img = Image.fromarray((arr * 255).astype('uint8'), mode='L')

            imgPixMap = self.imageToQPixmap(img)

            self.setQtImage(imgPixMap)

    # Load a geoimage in gdal format into a dictionary
    # Using the dictionary allows for multiple rasters to be loaded at once
    # Currently loading them directly into memory which will be a limiting factor.
    def importGeoImage(self, geoImagePath=None):

        if geoImagePath == None:

            geoImagePath = QtWidgets.QFileDialog.getOpenFileName(self, "Import GeoImage", "", ".tif(*.tif)")

            self.geoImageIndex[str(geoImagePath)] = gdal.Open(geoImagePath)

        else:

            self.geoImageIndex[str(geoImagePath)] = gdal.Open(geoImagePath)

    def displayGeoImage(self, geoImage, band=0):
        pass

# GetGeoTransform
# 0 == x coordinate
# 1 == grid size
# 3 == y coordinate
