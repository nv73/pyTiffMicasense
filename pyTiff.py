from osgeo import gdal
import numpy as np
from PIL import Image, ImageQt
import os
import osr
import exifread
import pyproj

class pyTiff(object):

    def __init__(self, filename, QtTools = False):

        self.filename = filename

        self.tiff = gdal.Open(filename)

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

        #leaving a default of wgs84
        self.epsg = 4326

        if QtTools == True:

            try:

                from PyQt4 import QtCore, QtGui

            except e as Exception:

                print("PyQt4 Packages not installed.")

        self.load_tif()

        #Generate a world file if none is found

        self.worldfile = os.path.splitext(filename)[0] + ".tfw"

        self.meta = self.get_metadata(self.filename)

        if self.worldfile not in os.listdir('.'):

            print("World file not found. Attempting to create...")

            try:
                
                self.make_world_file()

            except:

                print("This image has no metadata.")

        #print(self.geodetics)
        
    def get_metadata(self, filename):

        image = open(filename, 'rb')

        tags = exifread.process_file(image)

        return(tags)

    def load_tif(self):

        tiff = self.tiff
        
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

        for b in range(self.numberOfBands):

            b += 1

            activeBand = tiff.GetRasterBand(b)

            activeBand = np.array(activeBand.ReadAsArray())

            self.bands.append(activeBand)

        print(self.bands[0])

    def _rescale_image(self, arr):

        arr_range = (arr.min(), arr.max())

        return((arr - arr_range[0])/(arr_range[1] - arr_range[0]))

    def _to_degrees(self, coord):

        degrees = float(coord.values[0].num) / float(coord.values[0].den)
        
        minutes = float(coord.values[1].num) / float(coord.values[1].den)
        
        seconds = float(coord.values[2].num) / float(coord.values[2].den)
        
        return(degrees + (minutes / 60.0) + (seconds / 3600))

    def geo_to_utm(self, lat, lon, utmzone=9, ellipse='WGS84', hemi='W'):

        longitude = lon * -1
        latitude = lat

        proj = pyproj.Proj(proj='utm',zone=str(utmzone),ellps=ellipse, preserve_units=False)

        x2,y2 = proj(longitude, latitude)

        return((x2, y2))

    def image_from_bands(self, band1, band2 = None, band3 = None, save=True):
            
        band_w = self.bands[band1].shape[0]
        
        band_h = self.bands[band1].shape[1]      

        filename = str(os.path.splitext(self.filename)[0]) + ".jpeg"

        print("Height: %i, Width: %i" % (band_w, band_h))
        print("Number of bands: %i" % self.numberOfBands)
        print("X Corner: %i Y Corner: %i" % (self.x, self.y))
        
        if band2 == None or band3 == None:

            image = self._rescale_image(self.bands[band1]) * 255
            
            image = image.astype(np.uint8)

            image = Image.fromarray(image, mode='L')

            self.image = image

        else:

            rgb_image = np.zeros((band_w, band_h, 3))
            rgb_image[:,:,0] = self.bands[band1]
            rgb_image[:,:,1] = self.bands[band2]
            rgb_image[:,:,2] = self.bands[band3]

            image = Image.fromarray((self._rescale_image(rgb_image) * 255).astype(np.uint8), 'RGB')

            self.image = image

        if save == True:

            self.image.save(filename)

    def show_data(self):

        self.image.show()

    def write_array_to_tiff(self, band1, band2 = None, band3 = None, epsg=4326):

        b = self.numberOfBands

        filename = str(os.path.splitext(self.filename)[0]) + ".proc.tif"
        
        outRaster = gdal.GetDriverByName("GTiff").Create(filename,
                                                         self.bands[0].shape[1],
                                                         self.bands[0].shape[0],
                                                         b, gdal.GDT_Float32)

        outRaster.SetGeoTransform(self.geoTransform)

        srs = osr.SpatialReference()

        srs.ImportFromEPSG(epsg)

        outRaster.SetProjection(srs.ExportToWkt())

        #the if statement is unneccesary 
        if band2 == None or band3 == None:

            outRaster.GetRasterBand(1).WriteArray(band1)

        else:

            for x in range(b):

                x += 1

                outRaster.GetRasterBand(x).WriteArray(self.bands[x-1])
                
        outRaster.FlushCache()

    def write_bands_to_tiff(self, bands, epsg=4326):

        b = len(bands)

        filename = str(os.path.splitext(self.filename)[0]) + ".proc.tif"

        outRaster = gdal.GetDriverByName("GTiff").Create(filename,
                                                         self.bands[0].shape[1],
                                                         self.bands[0].shape[1],
                                                         b, gdal.GDT_Float32)

        outRaster.SetGeoTransform(self.geoTransform)

        srs = osr.SpatialReference()

        srs.ImportFromEPSG(epsg)

        outRaster.SetProjection(srs.ExportToWkt())

        for x in range(b):

            x += 1

            outRaster.GetRasterBand(x).WriteArray(self.bands[bands[x-1]])

        outRaster.FlushCache()

    def make_world_file(self):

        metadata = self.meta

        lat = self._to_degrees(metadata["GPS GPSLatitude"])
        
        lon = self._to_degrees(metadata["GPS GPSLongitude"])

        x,y = self.geo_to_utm(lat, lon)
        
        alt = metadata["GPS GPSAltitude"]
        alt = float(alt.values[0].num) / float(alt.values[0].den)
        resx = self.cellSize[0]
        resy = self.cellSize[1]
        rot1 = self.rotations[0]
        rot2 = self.rotations[1]

        with open(self.worldfile, 'a') as file:

            file.write(str(resx) + '\n' + str(rot1) + '\n' + str(rot2) + '\n' +
                       str(resy) + '\n' + str(x) + '\n' + str(y))           
        
    def load_hyperspectral(self):

        pass

#b = pyTiff("OWK-BASIN2-2M-AVG-BATHY1.tif")
b = pyTiff("KoeyeImagerySubset.tif")
#b = pyTiff("IMG_0012_1.tif")
#b = pyTiff('test.tif')
#b.image_from_bands(0)
#b.raster_from_bands()
b.write_array_to_tiff(0,1,2,epsg=3005)
