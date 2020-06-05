from osgeo import gdal
import numpy as np
from PIL import Image, ImageQt
import os
import osr
import exifread
import pyproj
import spectral.io.envi as envi
import matplotlib.pyplot as plt
import micasense.metadata as micaMetadata
import micasense.utils as utils
from micasense.image import Image
import pyzbar.pyzbar as pyzbar
import mapboxgl
from micasense.panel import Panel
import micasense.plotutils as plotutils
import exiftool
import cv2
import math

######Scroll to the bottom to see how to use this! :)

class pyTiff(object):

    def __init__(self, filename, QtTools = False, hyperspectral=False, overrideResolution=0, 
                 calibPanel=None, micasense=False, exiftoolpath=r"C:\exiftool\exiftool.exe", epsg=26909):

        self.filename = filename

        self.geodetics = None

        self.x = 0

        self.y = 0

        self.bands = []

        self.numberOfBands = 0

        self.hasWorldFile = True

        self.image = None

        self.cellSize = None

        self.geoTransform = None

        self.rotations = None

        self.worldfile = None
        
        self.meta = None

        self.epsg = epsg

        self.imageSize = None
                
        self.calibrationImagePath = calibPanel
        
        self.calibrationValues = None
        
        self.panel = None
        
        self.isMicasense = micasense
        
        self.exiftoolpath = exiftoolpath
        
        self.reflectanceImage = None
        
        self.orientation = None
        
        self.lat = None
        
        self.lon = None
        
        self.altitude = None
        
        self.scaleFactor = 0
        
        self.procFile = None
        
        if hyperspectral == False:

            self.tiff = gdal.Open(self.filename)

            self._load_tif()

            if overrideResolution > 0:

                self.cellSize = (overrideResolution, overrideResolution)

            if self.geoTransform[0] == 0.0:

                self.meta = self.get_metadata(self.filename)

                self.georeference_from_metadata()

        else:

            headerfilename = filename + '.hdr'

            self._load_hyperspectral(headerfilename, filename)

            self._decode_envi_header(headerfilename)
    
###The following 4 functions are specific to images taken by micasense cameras.###

    #Get calibration values from an image of the micasense calibration panel
    def get_calibration(self, plot=False):

        imgPath = self.calibrationImagePath
        img = plt.imread(self.calibrationImagePath)
        meta = micaMetadata.Metadata(imgPath, self.exiftoolpath)
        
        #Generate radiance image
        calib_rads = self.raw_to_rad(meta, img)
        
        panel = Panel(Image(imgPath))
        
        self.panel = panel.get_corners()
        corners = self.panel
        panel.plot(figsize=(8,8))
        
        x1 = corners[3][0] #upper left column coord (x) ulx
        y1 = corners[3][1] #upper left row coord (y) uly
        x2 = corners[1][0] #lower right column coord (x) lrx
        y2 = corners[1][1] #lower right row coord (y) lry

        #Just the values from within the calibration panel
        panelRegion = calib_rads[y1:y2, x1:x2]
        print("Panel Corners: %f, %f, %f, %f" % (x1, y1, x2, y2))
        
        #Mean radiance value within the panel
        meanRadiance = panelRegion.mean()
        print("Mean Radiance in panel region: %f" % meanRadiance)
        
        #Calibration values
        bandName = meta.get_item('XMP:BandName')
        
        panelCalibration = {
                "Blue": 0.54,
                "Green": 0.55,
                "Red": 0.54,
                "Red edge": 0.52,
                "NIR": 0.48
                }
        
        
        panelReflectance = panelCalibration[bandName]

        radianceToReflectance = panelReflectance / meanRadiance
        print("Radiance to Reflection conversion factor: {:1.3f}".format(radianceToReflectance))

        reflectanceImage = calib_rads * radianceToReflectance
    
        #Check quality of calibration image

        gaussianPlot = cv2.GaussianBlur(reflectanceImage[y1:y2, x1:x2],(55,55),5)
    
        print('Min Reflectance in panel region: {:1.2f}'.format(gaussianPlot.min()))
        print('Max Reflectance in panel region: {:1.2f}'.format(gaussianPlot.max()))
        print('Mean Reflectance in panel region: {:1.2f}'.format(gaussianPlot.mean()))
        print('Standard Deviation in region: {:1.4f}'.format(gaussianPlot.std()))

        undistortedReflectance = utils.correct_lens_distortion(meta, reflectanceImage)
        
        if plot == True:
        
            plotutils.plotwithcolorbar(panel, 'Raw image')
            plotutils.plotwithcolorbar(reflectanceImage, 'Reflectance Image')
            plotutils.plotwithcolorbar(gaussianPlot, 'Panel Region with Gaussian Blur')
            plotutils.plotwithcolorbar(undistortedReflectance, 'Undistorted Reflectance')
            
        self.calibrationValues = panelCalibration
        self.calibrationValues["meanPanelRadiance"] = meanRadiance
        
        self.scaleFactor = radianceToReflectance
    
    #Convert a raw micasense image to radiance
    def raw_to_rad(self, metadata, image):
        
        radiance, L, V, R = utils.raw_image_to_radiance(metadata, image)
        
        return(radiance)
    
    #Convert a radiance image to reflectance. 
    def get_reflectance(self, write_to_file = True):

        #I plan on using bands = 0 as a keyword later. For now, I am just coding it here.
        bands = 0

        scaleFactor = self.scaleFactor
            
        img_raw = plt.imread(self.filename)
            
        meta = micaMetadata.Metadata(self.filename, exiftoolPath = self.exiftoolpath)
    
        img_rad,_,_,_ = utils.raw_image_to_radiance(meta, img_raw)
        
        img_ref = img_rad * scaleFactor

        img_ref_undist = utils.correct_lens_distortion(meta, img_ref)
        
        self.reflectanceImage = img_ref_undist
        
        if write_to_file == True:
            
            filename = str(os.path.splitext(self.filename)[0]) + "proc.ref.tif"
            
            self.procFile = filename
            
            outRaster = gdal.GetDriverByName("GTiff").Create(filename, self.reflectanceImage.shape[1],
                                                             self.reflectanceImage.shape[0], 1, gdal.GDT_Float32)
            
            outRaster.SetGeoTransform(self.geoTransform)
            
            srs = osr.SpatialReference()
            
            srs.ImportFromEPSG(self.epsg)
            
            outRaster.SetProjection(srs.ExportToWkt())
            
            print("\nNumber of bands in source image: %i" % self.numberOfBands)
            
            print("Raster export dimensions: %f, %f" % (self.reflectanceImage.shape[1], self.reflectanceImage.shape[0]))
            
            print("Calibration image path: %s" % self.calibrationImagePath)
            
            print("Calibration values: %s" % self.calibrationValues)
            
            outRaster.GetRasterBand(1).WriteArray(self.reflectanceImage)
                    
            outRaster.FlushCache()
            
        return(self.reflectanceImage)
             
    #Extract the metadata from a file.
    #Used in the case that no tfw file exists with the data,
    #but EXIF info is provided.
    def get_metadata(self, filename):
        
        if self.isMicasense == False:

            image = open(filename, 'rb')

            tags = exifread.process_file(image)

            return(tags)
            
        else:
            
            meta = micaMetadata.Metadata(filename, exiftoolPath=self.exiftoolpath)
            
            return(meta)

###End micasense functions#######################################################

    #Loads the tif data into class variables.
    #Note that the file itself is loaded once the object is instantiated,
    #this function merely extracts the required data and saves it as
    #class variables.
    def _load_tif(self):

        tiff = self.tiff

        #Extract all the information required to position the data
        geodetics = tiff.GetProjection()[6:]
        
        for p in ['[', ']', '"']:

            geodetics = geodetics.replace(p, '')

        geodetics = geodetics.split(',')

        self.geodetics = geodetics

        self.geoTransform = tiff.GetGeoTransform()

        self.x = self.geoTransform[0]

        self.y = self.geoTransform[3]

        self.cellSize = (self.geoTransform[1], self.geoTransform[5])

        self.rotations = (self.geoTransform[2], self.geoTransform[4])

        self.numberOfBands = tiff.RasterCount

        #Iterate through each band and save them to a class variable.
        #This results in a list of lists.
        for b in range(self.numberOfBands):

            b += 1

            activeBand = tiff.GetRasterBand(b)

            activeBand = np.array(activeBand.ReadAsArray())

            self.bands.append(activeBand)

    #Rescale the data to values between 0 and 1. This allows
    #them to be easily converted later to uint8 (0-255) for the generation
    #of images.
    def _rescale_image(self, arr):

        arr_range = (arr.min(), arr.max())

        return((arr - arr_range[0])/(arr_range[1] - arr_range[0]))

    #Takes the exifread GPS position and converts it to decimal degrees.
    def _to_degrees(self, coord):

        #The GPS position data has weird ratio values, so we have to convert
        #them to float values (.num == numerator, .den == denominator
        degrees = float(coord.values[0].num) / float(coord.values[0].den)
        
        minutes = float(coord.values[1].num) / float(coord.values[1].den)
        
        seconds = float(coord.values[2].num) / float(coord.values[2].den)
        
        return(degrees + (minutes / 60.0) + (seconds / 3600))

    #Takes in a set of geographic coordinates (lat, lon) and converts them
    #to x, y coordinates (meters).
    #Defaults to WGS84, UTM 9N, Western hemisphere (Vancouver Island, BC)
    def geo_to_utm(self, lat, lon, utmzone=9, ellipse='WGS84', hemi='W'):

        if hemi == 'W':
            
            longitude = lon * -1

        else:

            longitude = lon
            
        latitude = lat

        proj = pyproj.Proj(proj='utm',zone=str(utmzone),ellps=ellipse, preserve_units=False)

        x2,y2 = proj(longitude, latitude)

        return((x2, y2))

    #Generates a Pillow Image from the input data. This is primarily for preparing
    #the dataset to be displayed in a QT application.
    #Supports 1 band (greyscale) and 3 band (RGB) images.
    def image_from_bands(self, band1, band2 = None, band3 = None, save=True):

        #Get the dimensions of the array to be written to the image
        band_w = self.bands[band1].shape[0]
        
        band_h = self.bands[band1].shape[1]      

        #Set up the filename for the output file.
        #Filetype is locked into jpeg for now.
        filename = str(os.path.splitext(self.filename)[0]) + ".jpeg"

        #If only one band is entered, create a greyscale image
        if band2 == None or band3 == None:

            image = self._rescale_image(self.bands[band1]) * 255
            
            image = image.astype(np.uint8)

            image = Image.fromarray(image, mode='L')

            self.image = image

        #If all three bands are entered, create an RGB image
        else:

            rgb_image = np.zeros((band_w, band_h, 3))
            rgb_image[:,:,0] = self.bands[band1]
            rgb_image[:,:,1] = self.bands[band2]
            rgb_image[:,:,2] = self.bands[band3]

            image = Image.fromarray((self._rescale_image(rgb_image) * 255).astype(np.uint8), 'RGB')

            self.image = image

        #The data has been stored to a class variable, but can
        #be output to a file if so desired.
        if save == True:

            self.image.save(filename)

    def show_data(self):

        self.image.show()
       
    #Takes a list containing the desired bands and outputs a tiff containing those bands
        #ie. bands = [1,3,4,7] would output a 4 band tiff with bands 1,3,4, and 7.
        #if bands == 0, all bands will be output 
    def write_bands_to_tiff(self, bands, epsg=4326):

        if bands == 0:

            bandLength = 1

        else:

            bandLength = len(bands)

        #The filename which will be assigned to the new output file
        filename = str(os.path.splitext(self.filename)[0]) + ".proc.tif"

        #Initialize the GeoTiff file to which the data will be written to
        outRaster = gdal.GetDriverByName("GTiff").Create(filename,
                                                         self.bands[0].shape[1],
                                                         self.bands[0].shape[0],
                                                         bandLength, gdal.GDT_UInt16)

        #Apply offsets, rotations, and resolution to the dataset
        outRaster.SetGeoTransform(self.geoTransform)
        
        #Initialize a spatial reference system
        srs = osr.SpatialReference()

        #Using the EPSG code, add a reference system to the object
        #This will default to WGS84
        srs.ImportFromEPSG(epsg)

        #Apply the reference system to the tiff file
        outRaster.SetProjection(srs.ExportToWkt())
        
        #Add the array dataset to the tiff file
        #If bands == 0, add all bands to the output tiff
        #Otherwise, add the input bands only to the output tiff
        if bands == 0:

            for x in range(self.numberOfBands):

                x+=1

                outRaster.GetRasterBand(x).WriteArray(self.bands[x-1])

        else:

            for x in range(len(bands)):

                x += 1

                outRaster.GetRasterBand(x).WriteArray(self.bands[bands[x-1]])

        #Dump the tiff data to a file.
        outRaster.FlushCache()

    #Extract geo data from image EXIF and use it for Tiff georeferencing
    def georeference_from_metadata(self):

        meta = self.meta
        
        lat = meta.get_item("EXIF:GPSLatitude")
        lon = meta.get_item("EXIF:GPSLongitude")
        alt = meta.get_item("EXIF:GPSAltitude")
        
        self.lat = lat
        self.lon = lon
        self.altitude = alt
        
        #Convert lat lon into x, y coordinates
        x,y = self.geo_to_utm(lat, lon)
        
        #We need to adjust image position for camera orientation        
        orientation = meta.dls_pose()
        yaw = orientation[0]
        pitch = orientation[1]
        roll = orientation[2]
        self.orientation = (roll, pitch, yaw)
        print("Altitude, Pitch, Roll (Radians) = %f, %f, %f, %f" % (alt, pitch, roll, yaw))
        print("Raw Position: X, Y = %f, %f" % (x, y))
        x = x - (alt * math.tan(roll))
        y = y + (alt * math.tan(pitch))
        
        print("Orientation: X, Y = %f, %f" % (x, y))
        
        #Adjust coordinates to represent the top left corner of image, as the currently
        #represent the centre of the image
        width = meta.get_item("EXIF:ImageWidth")
        height = meta.get_item("EXIF:ImageHeight")
        print("Width, Height: %f, %f" % (width, height))
        
        x = x - ((width * (self.cellSize[0]) / 2))
        y = y - ((height * (self.cellSize[1]) / 2))
        
        self.x = round(x, 3)
        self.y = round(y, 3)
        
        print("Offset to Corner: X, Y = %f, %f" % (self.x, self.y))

        #Build the geotransform based off the image EXIF data.
        #self.geoTransform = (self.x, self.cellSize[0], 0.0, self.y, 0.0, self.cellSize[1])  
        self.geoTransform = (self.x, math.cos(yaw) * self.cellSize[0], -math.sin(yaw) * self.cellSize[0], 
                             self.y, math.sin(yaw) * self.cellSize[1], math.cos(yaw) * self.cellSize[1])

    #Import hyperspectral data (ENVI format)
    def _load_hyperspectral(self, headerfilepath, imagefilepath):
        
        '''Load Hyperspectral Data'''
        #Open the header and image files
        hsData = envi.open(headerfilepath, imagefilepath)

        #Load the data into an array
        hsArray = hsData.load()

        #Get image dimensions
        x,y,bands = hsArray.shape

        self.numberOfBands = bands
        self.imageSize = (x,y)

        #Store each band to a class variable
        for x in range(bands):

            #Reshape the hyperspectral data slice into a 2d array
            band2D = hsArray[:,:,x].reshape(self.imageSize[0], self.imageSize[1])

            #add the data slice to the band list
            self.bands.append(band2D)

    #Load and decode the hyperspectral header file
    def _decode_envi_header(self, headerfilepath):

        #open the header file and load into a dictionary
        header = envi.read_envi_header(headerfilepath)

        #Extract the map info and coordinate system tags.
        #The coordinate system tag is not currently used, but I will
        #likely remove it in the future.
        geodetics = header['map info']
        coord_system = header['coordinate system string']

        #values required for geoTransform
        self.x = float(geodetics[3])
        self.y = float(geodetics[4])
        self.cellSize = (float(geodetics[5]),float(geodetics[6]))

        #Build the geoTransform required for exporting georeferenced tiffs
        self.geoTransform = (self.x, self.cellSize[0], 0.0, self.y, 0.0, self.cellSize[1])

#Right now the code relies on an accurate epsg code being used. Need to have auto projection available

if __name__ == "__main__":
    
    cal = r"C:\Users\nick.viner\Documents\GitHub\pyTiff\data\IMG_0001_5.tif"
    img = r"C:\Users\nick.viner\Documents\GitHub\pyTiff\data\IMG_0006_1.tif"
    b = pyTiff(img, calibPanel=cal, micasense=True, overrideResolution=0.0839)
    
    for item in b.meta.get_all():
        print(item)
     #   if "XMP" in item:
      #      print("{}: {}".format(item, b.meta.get_item(item)))

    
    #plt.imshow(reflectance)
    #b = pyTiff("OWK-BASIN2-2M-AVG-BATHY1.tif")
    #b = pyTiff(".\\data\\KoeyeImagerySubset.tif")
    #b = pyTiff(".\data\samson_1.img", hyperspectral=True)
    #b = pyTiff(".\data\IMG_0012_1.tif",overrideResolution=0.25)
    #b = pyTiff('test.tif')
    #b.image_from_bands(0,1,2)
    #b.write_bands_to_tiff(0, epsg=26909)
    #For testing purposes:
        #32610 = wgs84 utm10N
        #26909 = nad83 utm9N
    #b.write_bands_to_tiff(bands=[15,55,22,43,2,78],epsg=32610)

