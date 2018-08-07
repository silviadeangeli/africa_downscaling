__author__ = 'lauro'

from Geotiff import readFile, readFile_withNoData, writeGeotiffSingleBand
import numpy as np
from os.path import join
import os
from rasterRegrid import rasterRegrid
import matplotlib.pylab as plt


#########################################
############## FILE PATH  ###############
#########################################

# original input
sFileWorldPop =  "E:/africa_downscaling/Rwanda/hrsl_rwa.tif"
sFileFinalGrid = "E:/africa_downscaling/Rwanda/rw_box.tif"

# output directory
dir_out="E:/africa_downscaling/rw_2/"

if not os.path.exists (dir_out):
    os.mkdir(dir_out)

# regridded outputs
sFileWorldPop_regrid = join(dir_out, "rw_facebook_pop_90m.tif")

#########################################
#########################################
#########################################

# START

# regrid and read GRID
[xsize, ysize, geotransform, geoproj, data_grid]   = readFile_withNoData(sFileFinalGrid)

#regrid world pop to 90m
match_geotrans, match_proj = rasterRegrid(sFileWorldPop, sFileFinalGrid, sFileWorldPop_regrid ,"nearest")


