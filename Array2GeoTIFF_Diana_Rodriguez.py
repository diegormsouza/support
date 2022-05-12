#---------------------------------------------------------------------------------------------------------------------------
# Required modules
from osgeo import gdal, osr, ogr # Python bindings for GDAL
import numpy as np # Import the Numpy package
#---------------------------------------------------------------------------------------------------------------------------

def getGeoTransform(extent, nlines, ncols):
    resx = (extent[2] - extent[0]) / ncols
    resy = (extent[3] - extent[1]) / nlines
    return [extent[0], resx, 0, extent[3] , 0, -resy]

# Create a test 2D array (randon numbers between 0 and 9)
nlines = 30
ncolumns = 30
data = np.random.randint(0, 10, size=(nlines,ncolumns))

# Define the data extent (min. lon, min. lat, max. lon, max. lat)
extent = [-93.0, -60.00, -25.00, 18.00] # South America

# Export the test array to GeoTIFF ================================================

# Get GDAL driver GeoTiff
driver = gdal.GetDriverByName('GTiff')

# Get dimensions
nlines = data.shape[0]
ncols = data.shape[1]
nbands = len(data.shape)
data_type = gdal.GDT_Int16 # gdal.GDT_Float32

# Create a temp grid
#options = ['COMPRESS=JPEG', 'JPEG_QUALITY=80', 'TILED=YES']
grid_data = driver.Create('grid_data', ncols, nlines, 1, data_type)#, options)

# Write data for each bands
grid_data.GetRasterBand(1).WriteArray(data)

# Lat/Lon WSG84 Spatial Reference System
srs = osr.SpatialReference()
srs.ImportFromProj4('+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs')

# Setup projection and geo-transform
grid_data.SetProjection(srs.ExportToWkt())
grid_data.SetGeoTransform(getGeoTransform(extent, nlines, ncols))

# Save the file
file_name = 'my_test_data.tif'
print(f'Generated GeoTIFF: {file_name}')
driver.CreateCopy(file_name, grid_data, 0)	

# Close the file
driver = None
grid_data = None

# Delete the temp grid
import os                
os.remove('grid_data')

#==============================================================================