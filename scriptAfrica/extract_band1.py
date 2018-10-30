__author__ = 'silvia'

from Geotiff import readFile, readFile_withNoData, readFile_withNoData_band2, writeGeotiffSingleBand
import numpy as np
from os.path import join
import os
from rasterRegrid import rasterRegrid
import matplotlib.pylab as plt


#########################################s
############## FILE PATH  ###############
#########################################

# original input
sFileMask =  '/Users/silvia/Documents/AFRICA_DATA/Ghana/hrsl_gha.tif'

# output directory
dir_out= '/Users/silvia/Documents/AFRICA_DATA/Ghana/'

if not os.path.exists (dir_out):
    os.mkdir(dir_out)

# regridded outputs
#sFileWorldPop_regrid = join(dir_out, "rw_mask_facebook.tif")

#########################################
#########################################
#########################################

# START


[xsize, ysize, geotransform, geoproj, data_mask]   = readFile_withNoData(sFileMask)

writeGeotiffSingleBand(sFileMask+'_pop.tif', geotransform, geoproj, data_mask)
