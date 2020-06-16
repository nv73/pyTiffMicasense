-------Install Instructions-------

1. Open “Anaconda Prompt” (already installed on Mustang), you can download it from here: https://www.anaconda.com/products/individual

2. Run the command “activate micasense”

    •	There is a high likelihood this won’t immediately work (especially on any computer other than mustang)
    •	If so, use the prompt to navigate to the pyTiff folder in which the yaml file exists.
    •	Run the command “conda env create –f micasense.yaml”
             o	This could take a while
    •	Run the command “activate micasense”

3. Try the command “spyder”, if this doesn’t work, use “idle” instead

4. When the gui pops up, click anywhere within, then press F5 to run the application

5. If any errors arise, take a screenshot and send it to Nick 


------Processing Steps-------

File > Import Single File

	3 Windows will appear in the following order: Select Image, Select Calibration Panel Image, Enter Image Resolution 

		I would select the red band of the images you are working with to begin with

Process > Add Band

	Run this to add a new band to the image. Run this multiple times to select multiple bands
	
		Note that this will order the bands in whichever order you have selected them.

Process > Radiance To Reflectance

	This will convert all added bands to a reflectance image

	You will be asked if you wish to save to a file, click "Yes"

		This is needed to merge the bands to a single raster image, for simple data checks, the files needn't be created.

Process > Merge Bands to Raster

	Will merge all bands to a single georeferenced image and save it to the same directory where this script was ran.