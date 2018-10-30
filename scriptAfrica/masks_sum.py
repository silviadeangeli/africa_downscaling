__author__ = 'silvia'

from Geotiff_Silvia import readFile, readFile_withNoData, writeGeotiffSingleBand, readFileBand
import numpy as np
from os.path import join
import os
from rasterRegrid import rasterRegrid
import matplotlib.pylab as plt
from africa_tools import log_print

            
path_1 = '/Users/silvia/Documents/AFRICA_DATA/Guinea_Bissau/integrazione_buiA/gw_buia_01_box.tif'
path_2 = '/Users/silvia/Documents/AFRICA_DATA/Guinea_Bissau/regrid_outputs_20181018/GW_buiA_GHSL_GUF_LC_OSM_wp3_90m.tif'
#path_3 = '/Users/silvia/Documents/AFRICA_DATA/Namibia/OSM.tif'
path_output =  '/Users/silvia/Documents/AFRICA_DATA/Guinea_Bissau/integrazione_buiA/GW_buiA_GHSL_GUF_LC_OSM_wp3_90m_integration.tif'

xsize, ysize, geotransform, geoproj, data_1   = readFile(path_1)
xsize, ysize, geotransform, geoproj, data_2   = readFile(path_2)
#xsize, ysize, geotransform, geoproj, data_3   = readFile(path_3)

#data_tot = data_1 + data_2 + data_3
data_tot = data_1 + data_2
data_tot [data_tot > 0] = 1
data_tot [data_tot <= 0] = 0
writeGeotiffSingleBand(path_output, geotransform, geoproj, data_tot)