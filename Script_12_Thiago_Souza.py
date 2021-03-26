# Training: Python and GOES-R Imagery: Script 12 - Cropping the Full Disk
#-----------------------------------------------------------------------------------------------------------
# Required modules
from netCDF4 import Dataset              # Read / Write NetCDF4 files
import matplotlib.pyplot as plt          # Plotting library
from datetime import datetime            # Basic Dates and time types
import cartopy, cartopy.crs as ccrs      # Plot maps
import os                                # Miscellaneous operating system interfaces
from utilities import download_CMI       # Our own utilities
from utilities import geo2grid, convertExtent2GOESProjection      # Our own utilities
#-----------------------------------------------------------------------------------------------------------
# Input and output directories
input = "Samples"; os.makedirs(input, exist_ok=True)
output = "Output"; os.makedirs(output, exist_ok=True)

# Desired extent
extent = [-64.0, -36.0, -40.0, -15.0]  # Min lon, Max lon, Min lat, Max lat
'''
# AMAZON repository information 
# https://noaa-goes16.s3.amazonaws.com/index.html
bucket_name = 'noaa-goes16'
product_name = 'ABI-L2-CMIPF'
yyyymmddhhmn = '202102181800'
band = '13'

# Download the file
file_name = download_CMI(yyyymmddhhmn, band, input)
'''
#-----------------------------------------------------------------------------------------------------------
# Open the GOES-R image
file = Dataset(f'OR_ABI-L2-DSIF-M6_G16_s20210632010164_e20210632019472_c20210632020547.nc')
                   
# Convert lat/lon to grid-coordinates
lly, llx = geo2grid(extent[1], extent[0], file)
ury, urx = geo2grid(extent[3], extent[2], file)
        
# Get the pixel values
data = file.variables['CAPE'][ury:lly, llx:urx]          
#-----------------------------------------------------------------------------------------------------------
# Compute data-extent in GOES projection-coordinates
img_extent = convertExtent2GOESProjection(extent)
#-----------------------------------------------------------------------------------------------------------
# Choose the plot size (width x height, in inches)
plt.figure(figsize=(10,10))

# Choose the plot size (data size, in pixels)
# Print the dimension of our data
#print("Array dimension: ", data.shape)
#dpi = 150
#fig = plt.figure(figsize=(data.shape[1]/dpi, data.shape[0]/dpi), dpi=dpi)

# Use the Geostationary projection in cartopy
#ax = plt.axes([0, 0, 1, 1], projection=ccrs.Geostationary(central_longitude=-75.0, satellite_height=35786023.0))
ax = plt.axes(projection=ccrs.Geostationary(central_longitude=-75.0, satellite_height=35786023.0))

# Define the color scale based on the channel
colormap = "jet" # White to black for IR channels
        
# Plot the image
img = ax.imshow(data, origin='upper', extent=img_extent, cmap=colormap)

# Add coastlines, borders and gridlines
ax.coastlines(resolution='10m', color='black', linewidth=0.8)
ax.add_feature(cartopy.feature.BORDERS, edgecolor='black', linewidth=0.5)
ax.gridlines(color='black', alpha=0.5, linestyle='--', linewidth=0.5)

# Add a colorbar
plt.colorbar(img, label='CAPE', extend='both', orientation='horizontal', pad=0.05, fraction=0.05)

# Extract the date
date = (datetime.strptime(file.time_coverage_start, '%Y-%m-%dT%H:%M:%S.%fZ'))

# Add a title
plt.title('GOES-16 CAPE ' + date.strftime('%Y-%m-%d %H:%M') + ' UTC', fontweight='bold', fontsize=10, loc='left')
plt.title('Reg.: ' + str(extent) , fontsize=10, loc='right')
#-----------------------------------------------------------------------------------------------------------
# Save the image
plt.savefig(f'OR_ABI-L2-DSIF-M6_G16_s20210632010164_e20210632019472_c20210632020547.png', bbox_inches='tight', pad_inches=0, dpi=300)

# Show the image
plt.show()
