#-- my_code_hw02.py
#-- hw02 GEO1015.2022
#-- [YOUR NAME]
#-- [YOUR STUDENT NUMBER] 

from datetime import datetime
from pytz import timezone
import suncalc #-- https://pypi.org/project/suncalc/
import pyproj
import math
import numpy as np
import rasterio
from rasterio import features

#input for function is_sunny:
#dataset = rasterio.open('C:/Users/selin/OneDrive/Dokumente/Daten/Geophysik/Matlab/geo1015.2022/hw/01/dtm/hw2/ahn3_data/ahn3_dsm50cm_bk_small.asc')
#px = 85173.64
#py = 446877.28
#dt = '2022-12-01 14:45'

def is_sunny(dataset, px, py, dt):
    """
    !!! TO BE COMPLETED !!!
     
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
    
    # raise Exception("Point given is outside the extent of the dataset")
    b = dataset.bounds
    print(b)
    if b[0]>px or px >b[2] or b[1]>py or py> b[3]:
        raise Exception("Point given is outside the extent of the dataset")
    # raise Exception("Point given has no_data value")
    pxy = dataset.index(px, py) #row and column index of coordinates(x+y) of p
    if dataset.dataset_mask()[pxy] == 0:
        raise Exception("Point given has no_data value")


    #-- define the timezone for where Delft is located
    ams_tz = timezone('Europe/Amsterdam')
    #-- get the datetime object for a given local time at a given date
    dto = ams_tz.localize(datetime.fromisoformat(dt))
    #print(dto) 
    #-- now let's get the time UTC time, which must be used with suncalc 
    #-- notice that this gives us either +01:00 in winter 
    #-- and +02:00 because of EU summer time the local time
    time_utc = dto.astimezone(timezone('UTC'))
    #print(time_utc)

    #-- convert from EPSG:28992 to EPSG:4326
    transfo = pyproj.Transformer.from_crs("EPSG:28992", "EPSG:4326")
    latlon = transfo.transform(px,py)
    #print(latlon) 

    #calculate position of sun
    possun = suncalc.get_position(time_utc, latlon[1], latlon[0])
    #print(possun)

    ##finding a length(long enough) of the horizontal line on the 2D plane between the point and the sun
    b= dataset.bounds
    l= abs(b[0]-b[2])+abs(b[1]-b[3])
    #print(l)
    #finding x-y-position of sun: 
    sx= px-np.sin(possun['azimuth'])*l
    sy= py-np.cos(possun['azimuth'])*l
    #print(sx,sy)

    #make a line:
    v = {}
    v["type"] = "LineString"
    v["coordinates"] = []
    v["coordinates"].append([px,py]) #hier füge ich Angangswert ein also in Pixel-Indexierung
    v["coordinates"].append([sx,sy]) #hier Endpunkt
    shapes = [(v, 1)] #WOFÜR STEHT 1????
    #print(shapes)

    #rasterize our line
    #wir spezifizieren keinen Wert und nehmen an dass alle Pixel, welche von der Linie geschnitten werden 1 erhalten und die andren 0
    # re is a numpy array with d.shape where the line is rasterised (values != 0)
    re = features.rasterize(shapes, 
                        out_shape=dataset.shape,  #im Raster file habe ich dann eine Linie
                        all_touched=True,
                        transform=dataset.transform) # If the source dataset was georeferenced, you would get similarly georeferenced geometries like this:
    
    #die Pixel von rasterized line bekommen z-value of dataset:
    n1  = dataset.read(1)
    for i in range(re.shape[0]): #shape gives size of array in (n,m), N= nr. rows & m = nr. columns, shape[0] = nr.rows
        for j in range(re.shape[1]):#shape[1]= nr. of columns
            if re[i][j] == 1:
                re[i][j]=n1[i][j]

    #aus numpy array wieder raster machen:
    output_file = 'testing.tiff' 
    with rasterio.open(output_file, 'w', 
                       driver='GTiff', 
                       height=re.shape[0],
                       width=re.shape[1], 
                       count=1, 
                       dtype=rasterio.uint8,
                       crs=dataset.crs, 
                       transform=dataset.transform) as dst: #und hier sage ich dass, es ein Raster-dsm file werden soll??
        dst.write(re.astype(rasterio.uint8), 1)#Values for the height, width, and dtype keyword arguments are taken directly from attributes of the 2-D array (re)
    #1 ist die bandnr./Anzahl Variablen???
    print("File written to '%s'" % output_file)

    #above generated file öffnen:
    test = rasterio.open('testing.tiff')


    #height of sun for every middle point of pixel touching the line
    ##hmmm MAYBE WE NEED TO project the middle points of the touched pixel onto the line in the xy-plane
    sun_nosun_list= []
    for i in range(re.shape[0]): #shape gives size of array in (n,m), N= nr. rows & m = nr. columns, shape[0] = nr.rows
        for j in range(re.shape[1]):#shape[1]= nr. of columns
            if re[i][j] != 0:
                coordinates= test.xy([i], [j]); #we get the coordinates of the central point of the rasterized pixels
                c= np.array(coordinates, dtype=float) # c = x und y-coordinate, need to change to float numbers to do calculation below
                adj_side = np.sqrt(np.power(c[0]-px,2) + np.power(c[1]-py,2))
                z= np.multiply(np.tan(possun['altitude']),adj_side)
                sun = z > re[i][j] #sun position higher than dsm, logical value, true or false
                sun_nosun_list.append(sun) #append results to list
    answer = False in sun_nosun_list #is there one no-sun = false in the list?
    #print(sun_nosun_list)
    return answer #if there is one false element in the list it returns false 
    #but does it also return automatically true if there are only true (sun) elements??


