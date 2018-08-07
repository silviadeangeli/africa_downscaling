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

# original inputs
sFileFP =        "E:/africa_downscaling/zm/zm_FP.tif"
sFileGUF =       "E:/africa_downscaling/zm/zm_GUF.tif"
sFileGHSL =      "E:/africa_downscaling/zm/zm_GHSL.tif"
sFileWorldPop =  "E:/africa_downscaling/zm/zm_pop_WP.tif"
sFileFinalGrid = "E:/africa_downscaling/zm/zm_GAUL.tif"

# output directory
dir_out="E:/africa_downscaling/zm/masks_out"

if not os.path.exists (dir_out):
    os.mkdir(dir_out)

# regridded outputs
sFileFP_regrid = join(dir_out, "zm_buiA_FP.tif")
sFileGUF_regrid = join(dir_out, "zm_buiA_GUF.tif")
sFileGHSL_regrid = join(dir_out, "zm_buiA_GHSL.tif")
sFileWorldPop_regrid = join(dir_out, "zm_WP_90m.tif")

# regridded merged outputs
sFileGHSL_GUF_regrid = join(dir_out, "zm_buiA_GUF_GHSL.tif")
sFileGHSL_GUF_FP_regrid = join(dir_out, "zm_buiA_GUF_GHSL_FP.tif")
sFileMaskWPop2  = join(dir_out, "zm_buiA_GUF_GHSL_FP_WP2.tif")
sFileMaskWPop3  = join(dir_out, "zm_buiA_GUF_GHSL_FP_WP3.tif")
sFileMaskWPop5  = join(dir_out, "zm_buiA_GUF_GHSL_FP_WP5.tif")
sFileMaskWPop01 = join(dir_out, "zm_buiA_GUF_GHSL_FP_WP01.tif")

#########################################
#########################################
#########################################

# START

# regrid and read GRID
[xsize, ysize, geotransform, geoproj, data_grid]   = readFile_withNoData(sFileFinalGrid)

# regrid and read FP
match_geotrans, match_proj = rasterRegrid(sFileFP, sFileFinalGrid, sFileFP_regrid ,"nearest")
[xsize, ysize, geotransform, geoproj, data_FP]   = readFile(sFileFP_regrid)

# regrid and read GUF and GHSL
match_geotrans, match_proj = rasterRegrid(sFileGUF, sFileFinalGrid, sFileGUF_regrid ,"nearest")
[xsize, ysize, geotransform, geoproj, data_GUF]   = readFile(sFileGUF_regrid)


#os.system ("gdal_translate -of GTIFF -projwin 33.8 5.2 42.3 -5 -projwin_srs EPSG:4326 GHS_BUILT_LDSMT_GLOBE_R2015B_3857_38_v1_0.vrt Kenya.tif")
sFileGHSL4326 = sFileGHSL.split(".")[0]+"4326.tif"
os.system ("gdalwarp -t_srs EPSG:4326 "+sFileGHSL+" "+sFileGHSL4326)
match_geotrans, match_proj = rasterRegrid(sFileGHSL4326, sFileFinalGrid, sFileGHSL_regrid ,"nearest")
[xsize, ysize, geotransform, geoproj, data_GHSL]   = readFile(sFileGHSL_regrid)

data_GHSL [data_GHSL <3 ] =0
data_GHSL [data_GHSL >6 ] =0
data_GHSL [data_GHSL ==3] =1
data_GHSL [data_GHSL ==4] =1
data_GHSL [data_GHSL ==5] =1
data_GHSL [data_GHSL ==6] =1

#Save GHSL [0,1]
data_GHSL = data_GHSL * data_grid
writeGeotiffSingleBand(sFileGHSL_regrid, geotransform, geoproj, data_GHSL)

data_GUF_GHSL = (data_GUF + data_GHSL)
data_GUF_GHSL [data_GUF_GHSL >= 1] = 1
data_GUF_GHSL [data_GUF_GHSL < 1] = 0

#Save merge GUF + GHSL
data_GUF_GHSL = data_GUF_GHSL * data_grid
writeGeotiffSingleBand(sFileGHSL_GUF_regrid, geotransform, geoproj, data_GUF_GHSL)
del data_GUF, data_GHSL


data_GUF_GHSL_FP = (data_GUF_GHSL + data_FP)
data_GUF_GHSL_FP [data_GUF_GHSL_FP >= 1] = 1
data_GUF_GHSL_FP [data_GUF_GHSL_FP < 1] = 0

#Save merge GUF + GHSL + FP
data_GUF_GHSL_FP = data_GUF_GHSL_FP * data_grid
writeGeotiffSingleBand(sFileGHSL_GUF_FP_regrid, geotransform, geoproj, data_GUF_GHSL_FP)

#regrid world pop to 90m
match_geotrans, match_proj = rasterRegrid(sFileWorldPop, sFileFinalGrid, sFileWorldPop_regrid ,"nearest")
#read world pop to 90m

[xsize, ysize, geotransform, geoproj, data_WorldPop]   = readFile(sFileWorldPop_regrid)

data_merge_WorldPop2 = np.copy(data_WorldPop)
data_merge_WorldPop2 [data_merge_WorldPop2 < 2]= 0
data_merge_WorldPop2 [data_merge_WorldPop2 >= 2]= 1
data_merge_WorldPop2 = (data_GUF_GHSL_FP + data_merge_WorldPop2)
data_merge_WorldPop2 [data_merge_WorldPop2 >= 1]= 1

#Save merge World Pop 2
data_merge_WorldPop2 = data_merge_WorldPop2 * data_grid
writeGeotiffSingleBand(sFileMaskWPop2, geotransform, geoproj, data_merge_WorldPop2)
del data_merge_WorldPop2

data_merge_WorldPop3 = np.copy(data_WorldPop)
data_merge_WorldPop3 [data_merge_WorldPop3 < 3]= 0
data_merge_WorldPop3 [data_merge_WorldPop3 >= 3]= 1
data_merge_WorldPop3 = (data_GUF_GHSL_FP+ data_merge_WorldPop3)
data_merge_WorldPop3 [data_merge_WorldPop3 >= 1]= 1

#Save merge World Pop 3
data_merge_WorldPop3 = data_merge_WorldPop3 * data_grid
writeGeotiffSingleBand(sFileMaskWPop3, geotransform, geoproj, data_merge_WorldPop3)
del data_merge_WorldPop3

data_merge_WorldPop5 = np.copy(data_WorldPop)
data_merge_WorldPop5 [data_merge_WorldPop5 < 5]= 0
data_merge_WorldPop5 [data_merge_WorldPop5 >= 5]= 1
data_merge_WorldPop5 = (data_GUF_GHSL_FP+ data_merge_WorldPop5)
data_merge_WorldPop5 [data_merge_WorldPop5 >= 1]= 1

#Save merge World Pop 5
data_merge_WorldPop5 = data_merge_WorldPop5 * data_grid
writeGeotiffSingleBand(sFileMaskWPop5, geotransform, geoproj, data_merge_WorldPop5)
del data_merge_WorldPop5

data_merge_WorldPop01 = np.copy(data_WorldPop)
data_merge_WorldPop01 [data_merge_WorldPop01 < 0.1]= 0
data_merge_WorldPop01 [data_merge_WorldPop01 >= 0.1]= 1
data_merge_WorldPop01 = (data_GUF_GHSL_FP+ data_merge_WorldPop01)
data_merge_WorldPop01 [data_merge_WorldPop01 >= 1]= 1

#Save merge World Pop 10
data_merge_WorldPop01 = data_merge_WorldPop01 * data_grid
writeGeotiffSingleBand(sFileMaskWPop01, geotransform, geoproj, data_merge_WorldPop01)
