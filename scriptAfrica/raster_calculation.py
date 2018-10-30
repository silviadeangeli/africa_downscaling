__author__ = 'silvia'

from Geotiff_Silvia import readFile, readFile_withNoData, writeGeotiffSingleBand, readFileBand
import numpy as np
from os.path import join
import os
from rasterRegrid import rasterRegrid
import matplotlib.pylab as plt
from africa_tools import log_print
from high_res_mask import write_hres_mask
from raster_sum import raster_sum

            
#inputs
path = '/Users/silvia/Documents/AFRICA_DATA/Sz/Crops'
path_file = join(path, 'sz_agrV_sugar_cane_HighRes.tif')
value = 40.0

    
# outputs
path_output = join(path, 'sz_agrP_sugar_cane_other.tif')

xsize, ysize, geotransform, geoproj, data   = readFile(path_file)

data_new = data/value
writeGeotiffSingleBand(path_output, geotransform, geoproj, data_new)



