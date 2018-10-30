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
sFileOriginalRaster = "/Users/silvia/Documents/AFRICA_DATA/Botswana/BW_pop_HRSL_90m_sum.tif"  
sFileFinalGrid = "/Users/silvia/Documents/AFRICA_DATA/Botswana/original_layers/bw_box.tif"  

# output directory
dir_out= "/Users/silvia/Documents/AFRICA_DATA/Botswana/regridded_hrsl"

if not os.path.exists (dir_out):
    os.mkdir(dir_out)

# regridded outputs
sFileOriginalRaster_regrid = join(dir_out, "hrsl_bw_pop_sum_90m_regridded.tif")
sFileOriginalRaster_regrid_m = join(dir_out, "hrsl_bw_buiA_mask_90m_regridded.tif")


#########################################
#########################################

# START

# regrid and read GRID
[xsize, ysize, geotransform, geoproj, data_grid]   = readFile_withNoData(sFileFinalGrid)

# regrid and read OriginalRaster
match_geotrans, match_proj = rasterRegrid(sFileOriginalRaster, sFileFinalGrid, sFileOriginalRaster_regrid ,"nearest")

[xsize, ysize, geotransform, geoproj, data_OriginalRaster_regrid]   = readFile_withNoData(sFileOriginalRaster_regrid)

writeGeotiffSingleBand(sFileOriginalRaster_regrid,geotransform,geoproj, data_OriginalRaster_regrid)

data_OriginalRaster_regrid [data_OriginalRaster_regrid > 0 ] = 1
writeGeotiffSingleBand(sFileOriginalRaster_regrid_m,geotransform,geoproj, data_OriginalRaster_regrid)







