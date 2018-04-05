"""Script to create a set of test profiles given an input raster and a line shapefile. 
Created by Taylor Smith on 5.4.2018."""

import numpy as np
import os, json
from osgeo import gdalnumeric

import matplotlib.pyplot as plt

from create_swath import *

#Create boxes from shapefile
line = #Input line feature
boxes = #Output boxes along that line
if os.path.exists(boxes):
    for fid in glob.glob(line.split('.')[0] + '*'):
        os.remove(fid)
slice_width =  #In map units
overlap = #In map units
edge_length = #In map units
write_shapefile(boxes, create_boxes(line, slice_width, overlap, edge_length, 1000), cs)
print 'Clipper created...'

#Set the input raster
input_raster = #For example, a DEM
output_image_directory = #Where to save the profile images
tmpdir = #Set a directory to store temporary files

#Load the boxes and use them to create clips of the underlying raster
driver = ogr.GetDriverByName("ESRI Shapefile")
dataSource = driver.Open(boxes, 0)
layer = dataSource.GetLayer()
cs = layer.GetSpatialRef()

for feature in layer:
    #Get the id of the slice
    slice_id = feature.GetField('id')
    
    #Get corner coords and centroid
    ul, ur, lr, ll = json.loads(feature.ExportToJson())['geometry']['coordinates'][0][0:4]
    xcent = np.nanmean([ul[0],ur[0]])
    
    #Use bottom of box to compute rotation angle
    alpha = get_angle(ll, lr, [lr[0],ll[1]])
    
    #Use the box to clip the underlying raster
    subarr_fid = Subset_Array(feature, cs, input_raster, tmpdir)
    subarr = gdalnumeric.LoadFile(subarr_fid).astype(float)
    subarr[subarr < 0] = np.nan #Useful to flag NaNs before rotation of the image

    #Create 'flat' profiles from the underlying data
    profile, profile_std = create_profile(subarr, alpha)
    
    #Create a swath profile image
    fig, ax = plt.subplots(1)
    ax.plot(range(len(profile)), profile, 'k-', linewidth=2)
    ax.fill_between(range(len(profile)), profile - profile_std, profile + profile_std, facecolor='grey', alpha=0.5)
    
    fig.savefig(output_image_directory + str(slice_id) + '_profile.jpg', dpi=300)