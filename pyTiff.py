from osgeo import gdal
import numpy as np
from PIL import Image, ImageQt
import os
import osr
import exifread
import pyproj
import spectral.io.envi as envi

class pyTiff(object):

    def __init__(self, filename, QtTools = False, hyperspectral=False):

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

        self.epsg = 4326

        self.imageSize = None

        #Qt tools take up a lot of space and can use
        #a lot of runtime importing. 
        if QtTools == True:

            try:

                from PyQt4 import QtCore, QtGui

            except e as Exception:

                print("PyQt4 Packages not installed.")

        if hyperspectral == False:

            self.tiff = gdal.Open(self.filename)

            self.load_tif()

            self.get_world_file(self.filename)

        else:

            headerfilename = filename + '.hdr'

            self.load_hyperspectral(headerfilename, filename)

            self.decode_envi_header(headerfilename)

    #checks for a world file (.tfw). Creates one if none exist.
    def get_world_file(self, filename):

        fname = str(filename)
        
        self.worldfile = os.path.splitext(filename[0] + ".tfw")

        self.meta = self.get_metadata(filename)

        filepath = '\\'.join(filename.split('\\')[0:-1])

        if os.path.basename(str(self.worldfile)) not in os.listdir(filepath):

            try:

                self.make_world_file()

            except:

                pass
            
    #Extract the metadata from a file.
    #Used in the case that no tfw file exists with the data,
    #but EXIF info is provided.
    def get_metadata(self, filename):

        image = open(filename, 'rb')

        tags = exifread.process_file(image)

        return(tags)

    #Loads the tif data into class variables.
    #Note that the file itself is loaded once the object is instantiated,
    #this function merely extracts the required data and saves it as
    #class variables.
    def load_tif(self):

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

        longitude = lon * -1
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

        #The filename which will be assigned to the new output file
        filename = str(os.path.splitext(self.filename)[0]) + ".proc.tif"

        #Initialize the GeoTiff file to which the data will be written to
        outRaster = gdal.GetDriverByName("GTiff").Create(filename,
                                                         self.bands[0].shape[1],
                                                         self.bands[0].shape[0],
                                                         len(bands), gdal.GDT_UInt16)

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

    #Outputs a .tfw file to accompany a tif and provide georeferencing
    def make_world_file(self):

        #Load the exif metadata
        metadata = self.meta

        #Extract lat and lon from the metadata
        lat = self._to_degrees(metadata["GPS GPSLatitude"])
        
        lon = self._to_degrees(metadata["GPS GPSLongitude"])

        #Convert the lat, lon to x, y coordinates
        x,y = self.geo_to_utm(lat, lon)

        #A .tfw requires the following ordered info: x resolution, rotation1, rotation2, y resolution, x, y
        alt = metadata["GPS GPSAltitude"]
        alt = float(alt.values[0].num) / float(alt.values[0].den)
        resx = self.cellSize[0]
        resy = self.cellSize[1]
        rot1 = self.rotations[0]
        rot2 = self.rotations[1]

        #create the file and write the requried values to it
        with open(self.worldfile, 'a') as file:

            file.write(str(resx) + '\n' + str(rot1) + '\n' + str(rot2) + '\n' +
                       str(resy) + '\n' + str(x) + '\n' + str(y))           
        
    #Import hyperspectral data (ENVI format)
    def load_hyperspectral(self, headerfilepath, imagefilepath):

        #Open the header and image files
        hsData = envi.open(headerfilepath, imagefilepath)

        #Load the data into an array
        hsArray = hsData.load()

        #Get image dimensions
        x,y,bands = hsArray.shape

        self.numberOfBands = bands
        self.imageSize = (x,y)

        print(self.imageSize)
        #Store each band to a class variable
        for x in range(bands):

            band2D = hsArray[:,:,x].reshape(self.imageSize[0], self.imageSize[1])

            self.bands.append(band2D)

    def decode_envi_header(self, headerfilepath):

        hdr = envi.open(headerfilepath)
        
        header = envi.read_envi_header(headerfilepath)

        geodetics = header['map info']
        coord_system = header['coordinate system string']

        #values required for geoTransform
        self.x = float(geodetics[3])
        self.y = float(geodetics[4])
        self.cellSize = (float(geodetics[5]),float(geodetics[6]))       
        
        zone = geodetics[7]
        datum = geodetics[9]
        units = geodetics[10].split('=')[1]

        self.geoTransform = (self.x, self.cellSize[0], 0.0, self.y, 0.0, self.cellSize[1])
        
        print(self.geoTransform)


#b = pyTiff("OWK-BASIN2-2M-AVG-BATHY1.tif")
#b = pyTiff(".\\data\\KoeyeImagerySubset.tif")
b = pyTiff(".\data\samson_1.img", hyperspectral=True)
#b = pyTiff("IMG_0012_1.tif")
#b = pyTiff('test.tif')
#b.image_from_bands(0)
#b.raster_from_bands()
b.write_bands_to_tiff(bands=[15,55,154,43,2,78],epsg=26910)
