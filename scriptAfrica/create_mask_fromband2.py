__author__ = 'lauro'

from Geotiff import readFile, readFile_withNoData, readFile_withNoData_band2, writeGeotiffSingleBand
import numpy as np
from os.path import join
import os
from rasterRegrid import rasterRegrid
import matplotlib.pylab as plt


#########################################
############## FILE PATH  ###############
#########################################

# original input
sFileMask =  "E:/africa_downscaling/Rwanda/hrsl_rwa.tif"
sFileFinalGrid = "E:/africa_downscaling/Rwanda/rw_box.tif"

# output directory
dir_out="E:/africa_downscaling/rw_2/"

if not os.path.exists (dir_out):
    os.mkdir(dir_out)

# regridded outputs
sFileWorldPop_regrid = join(dir_out, "rw_mask_facebook.tif")

#########################################
#########################################
#########################################

# START


[xsize, ysize, geotransform, geoproj, data_mask]   = readFile_withNoData_band2(sFileMask)

writeGeotiffSingleBand(sFileMask+'_2.tif', geotransform, geoproj, data_mask)



[xsize, ysize, geotransform, geoproj, data_grid]   = readFile_withNoData(sFileFinalGrid)
#regrid world pop to 90m
match_geotrans, match_proj = rasterRegrid(sFileMask+'_2.tif', sFileFinalGrid, sFileWorldPop_regrid ,"nearest")


