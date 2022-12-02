#-- my_code_hw02.py
#-- hw02 GEO1015.2022
#-- [Malte Schade]
#-- [5850282] 

from datetime import datetime
from pytz import timezone
import suncalc
import pyproj
import math
import numpy as np
import rasterio
from rasterio import features


def lin(a, b, c, d):
    """
    Calculate Line Intersection.
     
    Input:
        a:  Point A (Vector A-B)
        b:  Point B (Vector A-B)
        c:  Point C (Vector C-D)
        d:  Point D (Vector C-D)
    Output:
        [double, double]: Intersection Point of Lines [x,y] OR
        False: No intersection between lines
    """

    # calculate two-coordinate intersection point of two lines 
    t = ((a[0]-c[0])*(c[1]-d[1]) - (a[1]-c[1])*(c[0]-d[0])) / ((a[0]-b[0])*(c[1]-d[1]) - (a[1]-b[1])*(c[0]-d[0]))
    u = ((a[0]-c[0])*(a[1]-b[1]) - (a[1]-c[1])*(a[0]-b[0])) / ((a[0]-b[0])*(c[1]-d[1]) - (a[1]-b[1])*(c[0]-d[0]))
    
    # check if intersection point is between line-defining points
    if (0<=t and t<=1 and 0<=u and u<=1):
        return [a[0] + t*(b[0]-a[0]), a[1] + t*(b[1]-a[1])]
    else: 
        return False


def is_sunny(dataset, px, py, dt):
    """
    Does the sun shine there at that moment?
     
    Input:
        dataset: the rasterio input dataset
        px:  x-coordinate of p
        py:  y-coordinate of p
        dt:  ISO-formatted datetime ('YYYY-MM-DD HH:MM'), eg '2022-08-12 13:32'
    Output:
        True/False: True = ground is illuminated; False = ground is not illuminated
           (raise Exception if p is outside the extent of the dataset)
           (raise Exception if p is no_data value)
    """
    # get index of raster cell that contains view point
    pr = dataset.index(px, py)

    # check if view point is inside raster extend
    if (pr[0] >= dataset.height) | (pr[1] >= dataset.width) | (pr[0] < 0) | (pr[1] < 0):
        raise Exception("Point given is outside the extent of the dataset")
    
    # raster has no data at view point
    if dataset.dataset_mask()[pr] == 0:
        raise Exception("Point given has no_data value")
    
    # transform time to utc
    ams_tz = timezone('Europe/Amsterdam')
    dto = ams_tz.localize(datetime.fromisoformat(dt))
    time_utc = dto.astimezone(timezone('UTC'))

    # transform view point crs
    transfo = pyproj.Transformer.from_crs("EPSG:28992", "EPSG:4326")
    px_t, py_t = transfo.transform(px, py)

    # calculate sun position
    pos_sun = suncalc.get_position(time_utc, px_t, py_t)
    
    # read z-data into numpy array
    data = dataset.read(1)

    # read z value for viewpoint
    z0 = data[pr]

    # altitude <= 0  => night
    if pos_sun['altitude'] <= 0:
        return False
        
    # calculate ray vector end point with dataset bounds (a+b >= d(c))
    b = dataset.bounds
    v_len = abs(b[0]-b[2])+abs(b[1]-b[3])
    px2 = px-np.sin(pos_sun['azimuth'])*v_len
    py2 = py-np.cos(pos_sun['azimuth'])*v_len

    # intersect vector with raster and find all intersection points (GEOJSON format)
    shapes = [({'type':'LineString', 'coordinates':[(px, py), (px2,py2)]}, 1)]
    re = features.rasterize(shapes, out_shape=dataset.shape, all_touched=True, transform=dataset.transform)

    # find all point indices for intersected raster blocks; filter view point block
    pi = np.argwhere(re)
    pi = pi[~np.all(pi==pr, axis=1)]

    # no information outside boundary => no ray obstruction
    if not pi.size:
        return True
        
    # get dtm height for intersection blocks and replace no-data values by minimum of grid
    zi = np.array([data[p[0], p[1]] for p in pi])
    zi[zi > 10000] = data.min() # no object/terrain on earth higher than 10.000m; no-data values higher

    # get intersection block midpoints from intersection indices
    pix, piy = dataset.xy(pi[:,0], pi[:,1])

    # get raster resolution
    res = dataset.res

    # calculate intersection points of in-line-blocks with ray
    pin = []
    for ix, iy in zip(pix, piy):
        # get box corners
        p1 = (ix-res[0]/2, iy-res[1]/2)
        p2 = (ix-res[0]/2, iy+res[1]/2)
        p3 = (ix+res[0]/2, iy-res[1]/2)
        p4 = (ix+res[0]/2, iy+res[1]/2)
        
        # calculate intersection points with user defined function
        li1 = lin((px, py), (px2, py2), p1, p2)
        li2 = lin((px, py), (px2, py2), p2, p4)
        li3 = lin((px, py), (px2, py2), p4, p3)
        li4 = lin((px, py), (px2, py2), p3, p1)
        
        # filter non-intersecting edges
        li = np.array([li1, li2, li3, li4], dtype=object)
        li = li[np.argwhere(li)]
        
        # select minimum distance intersection point
        li_d = [np.linalg.norm(np.subtract([px,py], i[0])) for i in li]
        li = li[np.argwhere(li_d==min(li_d))][0,0,0]
        pin.append(li)
        
    # calculate distance between intersection points and view point
    p_dist = [np.linalg.norm(np.subtract([px,py], p)) for p in pin]

    # calculate ray z-values at intersection point location 
    zp = np.add(z0, np.multiply(np.tan(pos_sun['altitude']), p_dist))

    # compare ray height with dtm height
    comp = np.greater(zp, zi)

    # check if all values True
    result = np.all(comp)

    return result

