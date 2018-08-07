__author__ = 'silvia'

from Geotiff import readFile, readFile_withNoData, writeGeotiffSingleBand
import numpy as np
from os.path import join
import os
from rasterRegrid import rasterRegrid
import matplotlib.pylab as plt


#########################################
############## FILE PATH  ###############
#########################################

# original inputs
sFileOriginalRaster =  "E:/africa_downscaling/crops/ESA_CCI_LC_cut_crops.tif"
sFileFinalGrid = "E:/africa_downscaling/crops/zm_GAUL.tif"

# output directory
dir_out="E:/africa_downscaling/crops/cuts"

if not os.path.exists (dir_out):
    os.mkdir(dir_out)

# regridded outputs
sFileOriginalRaster_regrid = join(dir_out, "zm_crops.tif")
sFileOriginalRaster_masked = join(dir_out, "zm_crops_mask.tif")


#########################################
#########################################

# START

# regrid and read GRID
[xsize, ysize, geotransform, geoproj, data_grid]   = readFile_withNoData(sFileFinalGrid)

# regrid and read OriginalRaster
match_geotrans, match_proj = rasterRegrid(sFileOriginalRaster, sFileFinalGrid, sFileOriginalRaster_regrid ,"nearest")

[xsize, ysize, geotransform, geoproj, data_OriginalRaster_regrid]   = readFile_withNoData(sFileOriginalRaster_regrid)

data_OriginalRaster_masked = data_OriginalRaster_regrid * data_grid
writeGeotiffSingleBand(sFileOriginalRaster_masked,geotransform,geoproj, data_OriginalRaster_masked)









