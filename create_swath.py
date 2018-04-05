"""Set of helper functions to create profiles along a given line shapefile. 
Requires gdal, PIL and shapely.
Created by Taylor Smith on 5.4.2018"""

import numpy as np
import shapely.geometry as shpgeo
import json, math, os, time, subprocess

from PIL import Image
from osgeo import ogr, osr, gdal

def cut(line, distance):
    from shapely.geometry import LineString, Point
    # Cuts a line in two at a distance from its starting point
    if distance <= 0.0 or distance >= line.length:
        return [LineString(line)]
    coords = list(line.coords)
    for i, p in enumerate(coords):
        pdist = line.project(Point(p))
        if pdist == distance:
            return [
                LineString(coords[:i+1]),
                LineString(coords[i:])]
        if pdist > distance:
            cp = line.interpolate(distance)
            return [
                LineString(coords[:i] + [(cp.x, cp.y)]),
                LineString([(cp.x, cp.y)] + coords[i:])]

def line_slice(line, slice_width, max_slice):
    fid = ogr.Open(line)
    lyr = fid.GetLayer(0)
    shape = lyr.GetFeature(0)
    lstring = shpgeo.shape(json.loads(shape.ExportToJson())['geometry'])
    fts = {}
    i = 0
    while i < max_slice:
        if i == 0:
            ft, newbase = cut(lstring, slice_width) #This gives ft of size 0-slice_width, and newbase of size slice_width-end
            fts[i] = ft.to_wkt()
        else:
            try:
                ft, newbase = cut(newbase, slice_width)
                fts[i] = ft.to_wkt()
            except:
                i = max_slice #Break the loop if max_slice is bigger than the maximum number of slices from the line
        i += 1
    return fts
    
def create_polygon(key, cd):          
    ring = ogr.Geometry(ogr.wkbLinearRing)
    coords = cd[key]
    for c in coords:
        ring.AddPoint(float(c[0]), float(c[1]))

    # Create polygon
    poly = ogr.Geometry(ogr.wkbPolygon)
    poly.AddGeometry(ring)
    return poly.ExportToWkt()

def write_shapefile(out_shp, cd, cs):
    """
    https://gis.stackexchange.com/a/52708/8104
    """
    # Now convert it to a shapefile with OGR    
    driver = ogr.GetDriverByName('Esri Shapefile')
    ds = driver.CreateDataSource(out_shp)
    srs = osr.SpatialReference()
    srs.ImportFromWkt(cs)
    #layer = ds.CreateLayer('', srs, ogr.wkbMultiLineString)
    layer = ds.CreateLayer('', srs, ogr.wkbPolygon)
    # Add one attribute
    layer.CreateField(ogr.FieldDefn('id', ogr.OFTInteger))
    defn = layer.GetLayerDefn()

    ## If there are multiple geometries, put the "for" loop here

    # Create a new feature (attribute and geometry)
    for key in cd.keys():
        feat = ogr.Feature(defn)
        feat.SetField('id', key)
    
        # Make a geometry, from Shapely object
        geom = ogr.CreateGeometryFromWkt(create_polygon(key, cd))
        feat.SetGeometry(geom)
    
        layer.CreateFeature(feat)
        feat = geom = None  # destroy these

    # Save and close everything
    ds = layer = feat = geom = None   
 
def perp_pts(x, y, m, edge_length, edges):
    x0, y0, x1, y1 = edges
    if m == 0:
        #Flat edge cases
        if y0 == y1:
            y1 = y + edge_length
            y2 = y - edge_length
            x1 = x
            x2 = x
        if x0 == x1:
            y1 = y
            y2 = y
            x1 = x + edge_length
            x2 = x - edge_length
    else:
        #A line perpendicular to x-y will have slope (-1/m)
        m_perp = (-1/m)
        
        if m > 0:
            x1 = x + (edge_length / math.sqrt(1 + m_perp**2))
            y1 = y + ((edge_length * m_perp) / math.sqrt(1+m_perp**2))
            x2 = x - (edge_length / math.sqrt(1 + m_perp**2))
            y2 = y - ((edge_length * m_perp) / math.sqrt(1+m_perp**2))
        
        if m < 0:
            x1 = x - (edge_length / math.sqrt(1 + m_perp**2))
            y1 = y - ((edge_length * m_perp) / math.sqrt(1+m_perp**2))
            x2 = x + (edge_length / math.sqrt(1 + m_perp**2))
            y2 = y + ((edge_length * m_perp) / math.sqrt(1+m_perp**2))
            
    return x1, y1, x2, y2
 
def create_boxes(line, slice_width, overlap, edge_length, max_slice):
    linedict = line_slice(line, slice_width, max_slice)
    
    #Get first and last point on the line, use them to build the box
    corner_dict = {}
    for key in linedict.keys():
        l = linedict[key]
        x0, y0 = l.replace('LINESTRING (', '').replace(')', '').split(',')[0].split(' ')
        x1, y1 = l.replace('LINESTRING (', '').replace(')', '').split(',')[-1].split(' ')[1:]
        x0, y0, x1, y1 = float(x0), float(y0), float(x1), float(y1)
        #x0, y0, x1, y1  are the leftmost and rightmost vertices of the line segment you want to make a box around
        
        #Get the slope of the line
        m = (y1 - y0) / (x1 - x0)
        #xdif, ydif = np.abs(x0-x1),np.abs(y0-y1)
        
        #Add the buffer onto the ends...
        if not overlap == 0:                
            s = (overlap / slice_width) + 1 #Get scaling factor
            x0_n = x0 * (1+s)/2 + x1 * (1-s)/2
            y0_n = y0 * (1+s)/2 + y1 * (1-s)/2

            x1_n = x1 * (1+s)/2 + x0 * (1-s)/2
            y1_n = y1 * (1+s)/2 + y0 * (1-s)/2

            x0, y0, x1, y1 = x0_n, y0_n, x1_n, y1_n

        #Get the bounding box using the slope of the line
        topx_l, topy_l, botx_l, boty_l = perp_pts(x0, y0, m, edge_length, [x0, y0, x1, y1])
        topx_r, topy_r, botx_r, boty_r = perp_pts(x1, y1, m, edge_length, [x0, y0, x1, y1])
        
        corner_dict[key] = [[topx_l, topy_l], [topx_r, topy_r], [botx_r, boty_r], [botx_l, boty_l], [topx_l, topy_l]]

    return corner_dict
    
def Rotate2D(pts,cnt,ang):
    '''pts = {} Rotates points(nx2) about center cnt(2) by angle ang(1) in radian'''
    return np.dot(pts-cnt,np.array([[math.cos(ang),math.sin(ang)],[-math.sin(ang),math.cos(ang)]]))+cnt

def Load_And_Rotate(arr, alpha):
    arr[arr == 0] = np.nan
    arr[arr == -9999] = np.nan
    pix = Image.fromarray(arr)
    del arr
    arr = np.array(pix.rotate(alpha, resample=Image.BICUBIC, expand=True))
    arr[arr == 0] = np.nan
    return arr

def create_profile(arr, alpha):
    smallindex, firstval, lastval, firstval_s, lastval_s = whichaxis(arr)
    xvar = Load_And_Rotate(arr, alpha)
    if smallindex == 0:
        profile = np.nanmean(xvar[firstval_s:lastval_s,firstval:lastval], 0)[::-1]
    elif smallindex == 1:
        profile = np.nanmean(xvar[firstval:lastval,firstval_s:lastval_s], 1)[::-1]
    if smallindex == 0:
        pstd = np.nanstd(xvar[firstval_s:lastval_s,firstval:lastval], 0)[::-1]
    elif smallindex == 1:
        pstd = np.nanstd(xvar[firstval:lastval,firstval_s:lastval_s], 1)[::-1]

    return profile, pstd   

def whichaxis(rot):
    p1 = np.nanmean(rot, 0)
    p2 = np.nanmean(rot, 1)

    firstval = np.where(~np.isnan(p1))[0][0] + 2
    lastval = np.where(~np.isnan(p1))[0][-1] - 2
    
    firstval2 = np.where(~np.isnan(p2))[0][0] + 2
    lastval2 = np.where(~np.isnan(p2))[0][-1] - 2
    
    if lastval - firstval > lastval2 - firstval2:
        return 0, firstval, lastval, firstval2, lastval2
    else:
        return 1, firstval2, lastval2, firstval, lastval
    
def get_angle(p0, p1, p2):
    ''' compute angle (in degrees) for p0p1p2 corner
    Inputs:
        p0,p1,p2 - points in the form of [x,y]
    '''
    if p2 is None:
        p2 = p1 + np.array([1, 0])
    v0 = np.array(p0) - np.array(p1)
    v1 = np.array(p2) - np.array(p1)

    angle = np.math.atan2(np.linalg.det([v0,v1]),np.dot(v0,v1))
    return np.degrees(angle)
    
def Subset_Array(feature_object, cs_in, tifgrid, tmpdir):
    """ 
    Subset_Array(feature_object, cs_in, tifgrid, tmpdir)
    Takes a feature object and coordinate system
    from GDAL and returns a clip of the input tif 
    """
    randint = str(np.random.randint(0,1000000))
    geom = feature_object.geometry()
    #geom2 = geom.Buffer(0) #Repair faulty geometry
    tmp = tmpdir + randint + 'tmp.shp'
    drv = ogr.GetDriverByName('ESRI Shapefile')

    dst_ds = drv.CreateDataSource(tmp)
    dst_layer = dst_ds.CreateLayer('Clipper', srs=cs_in)
    f = ogr.Feature(dst_layer.GetLayerDefn())
    f.SetGeometry(geom)
    f.SetFID(0)
    dst_layer.CreateFeature(f)
    dst_ds.Destroy()
    
    c = gdal.Open(tifgrid)
    b = c.GetRasterBand(1)
    if not b.GetNoDataValue() == None:
        ndataval = str(b.GetNoDataValue())
    else:
        if b.GetMinimum() == -9999:
            ndataval = '-9999'
        else:
            ndataval = 'nan'

    out = tmpdir + randint + 'tmp.tif'
    FNULL = open(os.devnull, 'w')
    gdal_comm = ['gdalwarp', '-cutline', tmp, '-crop_to_cutline', '-srcnodata', ndataval, tifgrid, out]
    retcode = subprocess.call(gdal_comm, stdout=FNULL, stderr=subprocess.STDOUT)
    while os.path.exists(out) == False:
        time.sleep(1)
        
    os.remove(tmp)
    os.remove(tmpdir + randint + 'tmp.shx')
    os.remove(tmpdir + randint + 'tmp.prj')
    os.remove(tmpdir + randint + 'tmp.dbf')

    return out     