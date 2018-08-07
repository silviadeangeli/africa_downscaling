__author__ = 'silvia'

from Geotiff import readFile, readFile_withNoData, readFile_withNoData_band2, writeGeotiffSingleBand
import numpy as np
from os.path import join
import os
from rasterRegrid import rasterRegrid
import matplotlib.pylab as plt
from osgeo import gdal, gdalconst, ogr, gdal_translate

#########################################
############## FILE PATH  ###############
#########################################

# original input
#shape =  "E:/africa_downscaling/admin0GAUL.shp"
#raster = "E:/africa_downscaling/ESACCI-LC-L4-LC10-Map-20m-P1Y-2016-v1/ESACCI-LC-L4-LC10-Map-20m-P1Y-2016-v1.0.tif"
#output = "E:/africa_downscaling/LC_cut_Angola.tif"

#$shape =  "E:\africa_downscaling\admin0GAUL.shp"
#$raster = "E:\africa_downscaling\ESACCI-LC-L4-LC10-Map-20m-P1Y-2016-v1\ESACCI-LC-L4-LC10-Map-20m-P1Y-2016-v1.0.tif"
#$output = "E:\africa_downscaling\LC_cut_Angola.tif"


#gdalwarp -cutline shape raster output

gdal_translate -of GTIFF -projwin 9 -2.7 25 -19.2 -projwin_srs EPSG:4326 "E:\africa_downscaling\ESACCI-LC-L4-LC10-Map-20m-P1Y-2016-v1\ESACCI-LC-L4-LC10-Map-20m-P1Y-2016-v1.0.tif" "E:\africa_downscaling\LC_cut_Angola.tif"