__author__ = 'silvia'

from Geotiff_Silvia import readFile, readFile_withNoData, writeGeotiffSingleBand, readFileBand
import numpy as np
from os.path import join
import os
from rasterRegrid import rasterRegrid
import matplotlib.pylab as plt
from africa_tools import log_print


#########################################
############## FILE PATH  ###############
#########################################

path_valfis = "E:/africa_downscaling/ke/KE_valfis/"
sFile_C1 = path_valfis + "C1.tif"
sFile_M1 = path_valfis + "M1.tif"
sFile_T1 = path_valfis + "T1.tif"
sFile_M2 = path_valfis + "M2.tif"
sFile_W1 = path_valfis + "W1.tif"

path_p = "E:/africa_downscaling/ke/percent_int_20180805/"
sFile_C1_p = path_p + "C1.tiff"
sFile_M1_p = path_p + "M1.tiff"
sFile_T1_p = path_p + "T1.tiff"
sFile_M2_p = path_p + "M2.tiff"
sFile_W1_p = path_p + "W1.tiff"

# outputs
sFile_hous = "E:/africa_downscaling/ke/new_sectors_0805/ke_housing.tif"
sFile_ind = "E:/africa_downscaling/ke/new_sectors_0805/ke_industrial.tif"
sFile_serv = "E:/africa_downscaling/ke/new_sectors_0805/ke_services.tif"


[xsize, ysize, geotransform, geoproj, data_C1]   = readFile(sFile_C1)
[xsize, ysize, geotransform, geoproj, data_M1]   = readFile(sFile_M1)
[xsize, ysize, geotransform, geoproj, data_M2]   = readFile(sFile_M2)
[xsize, ysize, geotransform, geoproj, data_T1]   = readFile(sFile_T1)
[xsize, ysize, geotransform, geoproj, data_W1]   = readFile(sFile_W1)

[xsize, ysize, geotransform, geoproj, data_C1_p1]   = readFileBand(sFile_C1_p, 1)
[xsize, ysize, geotransform, geoproj, data_M1_p1]   = readFileBand(sFile_M1_p, 1)
[xsize, ysize, geotransform, geoproj, data_M2_p1]   = readFileBand(sFile_M2_p, 1)
[xsize, ysize, geotransform, geoproj, data_T1_p1]   = readFileBand(sFile_T1_p, 1)
[xsize, ysize, geotransform, geoproj, data_W1_p1]   = readFileBand(sFile_W1_p, 1)

[xsize, ysize, geotransform, geoproj, data_C1_p2]   = readFileBand(sFile_C1_p, 2)
[xsize, ysize, geotransform, geoproj, data_M1_p2]   = readFileBand(sFile_M1_p, 2)
[xsize, ysize, geotransform, geoproj, data_M2_p2]   = readFileBand(sFile_M2_p, 2)
[xsize, ysize, geotransform, geoproj, data_T1_p2]   = readFileBand(sFile_T1_p, 2)
[xsize, ysize, geotransform, geoproj, data_W1_p2]   = readFileBand(sFile_W1_p, 2)

[xsize, ysize, geotransform, geoproj, data_C1_p3]   = readFileBand(sFile_C1_p, 3)
[xsize, ysize, geotransform, geoproj, data_M1_p3]   = readFileBand(sFile_M1_p, 3)
[xsize, ysize, geotransform, geoproj, data_M2_p3]   = readFileBand(sFile_M2_p, 3)
[xsize, ysize, geotransform, geoproj, data_T1_p3]   = readFileBand(sFile_T1_p, 3)
[xsize, ysize, geotransform, geoproj, data_W1_p3]   = readFileBand(sFile_W1_p, 3)

[xsize, ysize, geotransform, geoproj, data_C1_p5]   = readFileBand(sFile_C1_p, 5)
[xsize, ysize, geotransform, geoproj, data_M1_p5]   = readFileBand(sFile_M1_p, 5)
[xsize, ysize, geotransform, geoproj, data_M2_p5]   = readFileBand(sFile_M2_p, 5)
[xsize, ysize, geotransform, geoproj, data_T1_p5]   = readFileBand(sFile_T1_p, 5)
[xsize, ysize, geotransform, geoproj, data_W1_p5]   = readFileBand(sFile_W1_p, 5)

[xsize, ysize, geotransform, geoproj, data_C1_p6]   = readFileBand(sFile_C1_p, 6)
[xsize, ysize, geotransform, geoproj, data_M1_p6]   = readFileBand(sFile_M1_p, 6)
[xsize, ysize, geotransform, geoproj, data_M2_p6]   = readFileBand(sFile_M2_p, 6)
[xsize, ysize, geotransform, geoproj, data_T1_p6]   = readFileBand(sFile_T1_p, 6)
[xsize, ysize, geotransform, geoproj, data_W1_p6]   = readFileBand(sFile_W1_p, 6)


data_hous = (data_C1*data_C1_p1+data_C1*data_C1_p2+data_M1*data_M1_p1+data_M1*data_M1_p2+data_M2*data_M2_p1+data_M2*data_M2_p2+data_T1*data_T1_p1+data_T1*data_T1_p2+data_W1*data_W1_p1+data_W1*data_W1_p2)/100

hous_tot = np.sum(data_hous)

log_print(f'total housing: [{hous_tot}]')


data_ind = (data_C1*data_C1_p6+data_M1*data_M1_p6+data_M2*data_M2_p6+data_T1*data_T1_p6+data_W1*data_W1_p6)/100

ind_tot = np.sum(data_ind)

log_print(f'total industrial: [{ind_tot}]')

data_serv = (data_C1*data_C1_p3+data_C1*data_C1_p5+data_M1*data_M1_p3+data_M1*data_M1_p5+data_M2*data_M2_p3+data_M2*data_M2_p5+data_T1*data_T1_p3+data_T1*data_T1_p5+data_W1*data_W1_p3+data_W1*data_W1_p5)/100

serv_tot = np.sum(data_serv)

log_print(f'total service: [{serv_tot}]')

writeGeotiffSingleBand(sFile_hous, geotransform, geoproj, data_hous)
writeGeotiffSingleBand(sFile_ind, geotransform, geoproj, data_ind)
writeGeotiffSingleBand(sFile_serv, geotransform, geoproj, data_serv)




