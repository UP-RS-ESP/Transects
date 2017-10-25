'''This program takes an input line shapefile and creates a set of boxes of arbitary width and height at a given spacing along the line. This is 
useful for many applications in earth science, when transects need to be made. For example, profiles along a mountain range.

Created by Taylor Smith on 25.10.17'''

import math, sys, os, glob
from osgeo import ogr, osr
import shapely.geometry as shpgeo
import json

def cut(line, distance):
    # Cuts a line in two at a distance from its starting point
    # See here: https://toblerity.org/shapely/manual.html
    if distance <= 0.0 or distance >= line.length:
        return [shpgeo.LineString(line)]
    coords = list(line.coords)
    for i, p in enumerate(coords):
        pdist = line.project(shpgeo.Point(p))
        if pdist == distance:
            return [
                shpgeo.LineString(coords[:i+1]),
                shpgeo.LineString(coords[i:])]
        if pdist > distance:
            cp = line.interpolate(distance)
            return [
                shpgeo.LineString(coords[:i] + [(cp.x, cp.y)]),
                shpgeo.LineString([(cp.x, cp.y)] + coords[i:])]

def line_slice(line, slice_width, max_slice):
    ''' Take an input shapefile and slice it into pieces based on a given width and maximum number of slices'''
    fid = ogr.Open(line)
    lyr = fid.GetLayer(0)
    shape = lyr.GetFeature(0)
    lstring = shpgeo.shape(json.loads(shape.ExportToJson())['geometry']) #Get the line shapefile as a WKT LineString
    fts = {}
    i = 0
    print 'Input loaded...'
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
    '''Create a rectangle polygon out of a given set of points'''
    ring = ogr.Geometry(ogr.wkbLinearRing)
    coords = cd[key]
    for c in coords:
        ring.AddPoint(float(c[0]), float(c[1]))

    poly = ogr.Geometry(ogr.wkbPolygon)
    poly.AddGeometry(ring)
    return poly.ExportToWkt()

def write_shapefile(out_shp, cd, cs):
    '''Write the results to a polygon shapefile, based on the 'cd' coordinate dictionary'''
    driver = ogr.GetDriverByName('Esri Shapefile')
    ds = driver.CreateDataSource(out_shp)
    srs = osr.SpatialReference()
    srs.ImportFromWkt(cs) #Get the spatial reference from the input line
    layer = ds.CreateLayer('', srs, ogr.wkbMultiLineString)
    layer.CreateField(ogr.FieldDefn('ID', ogr.OFTInteger))
    defn = layer.GetLayerDefn()

    for key in cd.keys():
        feat = ogr.Feature(defn)
        feat.SetField('ID', key)
        geom = ogr.CreateGeometryFromWkt(create_polygon(key, cd))
        feat.SetGeometry(geom)
    
        layer.CreateFeature(feat)
        del feat, geom
    del ds, layer
 
def perp_pts(x, y, m, edge_length, edges):
    '''Create points perpendicular to the line segment at a set distance'''
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
        
        #Use vector math to get points along perpendicular line
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
    '''Main function, takes input line and saves out output file'''
    #Slice the line into chunks...
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
        xdif, ydif = abs(x0-x1),abs(y0-y1)
        #Add the buffer onto the ends...
        if not overlap == 0:
            if y0 == y1:
                x0, x1 = x0 - overlap, x1 + overlap
            if x0 == x1:
                y0, y1 = y0 - overlap, y1 + overlap
            if m < 0:
                alpha = math.tan(ydif/xdif)
                ydist = math.sin(alpha) * overlap
                xdist = math.cos(alpha) * overlap
                x0, y0, x1, y1 = x0 - xdist, y0 + ydist, x1 + xdist, y1 - ydist
            if m > 0:
                alpha = math.tan(ydif/xdif)
                ydist = math.sin(alpha) * overlap
                xdist = math.cos(alpha) * overlap
                x0, y0, x1, y1 = x0 - xdist, y0 - ydist, x1 + xdist, y1 + ydist

        #Get the bounding box using the slope of the line
        topx_l, topy_l, botx_l, boty_l = perp_pts(x0, y0, m, edge_length, [x0, y0, x1, y1])
        topx_r, topy_r, botx_r, boty_r = perp_pts(x1, y1, m, edge_length, [x0, y0, x1, y1])
        
        corner_dict[key] = [[topx_l, topy_l], [topx_r, topy_r], [botx_r, boty_r], [botx_l, boty_l], [topx_l, topy_l]]

    return corner_dict
    
line = sys.argv[1]
output = sys.argv[2]        
slice_width = float(sys.argv[3])
overlap = float(sys.argv[4])
edge_length = float(sys.argv[5])
max_poly = int(sys.argv[6])
print 'Input:', line, 'Output:', output, 'Width:', slice_width, 'Overlap:', overlap, 'Edge Length:', edge_length, 'Max Polygons:', max_poly

fid = ogr.Open(line)
cs = fid.GetLayer().GetSpatialRef().ExportToWkt() #Get the coordinate system from the input for the output

#Clean up the shapefile if it already exists
if os.path.exists(output):
    print 'Deleting old shapefile...'
    for fid in glob.glob(output.split('.')[0] + '*'):
        os.remove(fid)

write_shapefile(output, create_boxes(line, slice_width, overlap, edge_length, max_poly), cs)
print 'Data Written.'
