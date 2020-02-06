# pyTiff
Streamlined handling of geotiff images in Python. Mostly keeping in mind compatibility with PyQT and multi-band rasters.

Requires the following packages:

	GDAL
	numpy
	pillow
	osr
	exifread
	pyproj
	spectral
	
All packages can be installed in Anaconda (open as administrator) via:

	conda install -c conda-forge package_name
	

pyTiff Class

	An object which loads and holds geoTiff or ENVI Hyyperspectral data.
	
	pyTiff(filename, QtTools=False, hyperspectral=False)
	
	filename = The filepath to the desired tiff / hyperspectral data file
	
	QtTools = Only for dev purposes. Ignore this
	
	hyperspectral = Boolean value. False if input data is not hyperspectral. True if it is.
	
	Example:
	
		tiff = pyTiff(".\tiffImage.tiff")
		
		hyperspec = pyTiff(".\hsData.img", hyperspectral=True)
		
	#### Class Functions ####
			
		image_from_bands(band1, band2=None, band3=None, save=True)
		
			Generates a pillow image object from the loaded tiff / hyperspectral data
			
			Will save the image to a .jpeg file if save == True (default)
			
		write_bands_to_tiff(bands, epsg=4326)
		
			arg bands can be either a list of integer values (starting at 0) or 0.
			
				if 0, then all bands will be written to a new tiff file
				
				otherwise, only the bands specified in the list will be written
				
					eg. pyTiffObject.write_bands_to_tiff([0,4,8], epsg=26909) will write bands 1, 5, and 9 to the tiff in NAD83 UTM 9N
					
			for proper georeferencing, the epsg code for the data projection must be specified.
			
			** Note that with the tiff files generated via this method, the georeferencing is encoded within the file. No .tfw is needed.
			
		show_data()
		
			image_from_bands has been run, this function can display the data temporarily.
	
	
