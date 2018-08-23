__author__ = 'silvia'

from Geotiff_Silvia import readFile, readFile_withNoData, writeGeotiffSingleBand, readFileBand
import numpy as np
from os.path import join
import os
from rasterRegrid import rasterRegrid
import matplotlib.pylab as plt
from africa_tools import log_print


import rasterio as rio

def readMultiBandGeotiff(file):
    """
        read a multiband geotiff and returns an array of 2d images
        if the band count is 1, it returns a single image
        
        example (multiband):       
            data_1, data_2, data_3, data_4, data_5, data_6 = readMultiBandGeotiff(file)

        example (singleband):       
            data = readMultiBandGeotiff(file)
        
    """

    with rio.open(file) as f:
        if f.count == 1:
            data = f.read(1)
            return data
        else:
            ret_vals = []
            for band in range(f.count):
                data = f.read(band+1)
                ret_vals.append(data)
            return ret_vals

#########################################
############## FILE PATH  ###############
#########################################

path_valfis = "E:/africa_downscaling/data/Rwanda/RW_valfis/"
path_p = "E:/africa_downscaling/data/Rwanda/percent_int_20180806/"

# outputs
sFile_hous = "E:/africa_downscaling/data/Rwanda/new_sectors_0806/rw_housing.tif"
sFile_ind = "E:/africa_downscaling/data/Rwanda/new_sectors_0806/rw_industrial.tif"
sFile_serv = "E:/africa_downscaling/data/Rwanda/new_sectors_0806/rw_services.tif"


sFile_C1 = join(path_valfis, "C1.tif")
sFile_M1 = join(path_valfis, "M1.tif")
sFile_M2 = join(path_valfis, "M2.tif")
sFile_W1 = join(path_valfis, "W1.tif")
sFile_T1 = join(path_valfis, "T1.tif")


sFile_C1_p = join(path_p, "C1.tiff")
sFile_M1_p = join(path_p, "M1.tiff")
sFile_M2_p = join(path_p, "M2.tiff")
sFile_W1_p = join(path_p, "W1.tiff")
sFile_T1_p = join(path_p, "T1.tiff")

xsize, ysize, geotransform, geoproj, data_C1   = readFile_withNoData(sFile_C1)
log_print(ysize)
log_print(xsize)
log_print(np.sum(data_C1))
xsize, ysize, geotransform, geoproj, data_M1   = readFile(sFile_M1)
xsize, ysize, geotransform, geoproj, data_M2   = readFile(sFile_M2)
xsize, ysize, geotransform, geoproj, data_W1   = readFile(sFile_W1)

data_hous = np.zeros((ysize, xsize))
data_ind = np.zeros((ysize, xsize))
data_serv = np.zeros((ysize, xsize))

if os.path.isfile(sFile_C1_p):
    data_C1_p1, data_C1_p2, data_C1_p3, _, data_C1_p5, data_C1_p6 \
        = readMultiBandGeotiff(sFile_C1_p)
    data_hous += (data_C1 * data_C1_p1 + data_C1 * data_C1_p2)/100
    data_ind  += (data_C1 * data_C1_p6)/100
    data_serv += (data_C1 * data_C1_p3 + data_C1 * data_C1_p5)/100

if os.path.isfile(sFile_M1_p):
    data_M1_p1, data_M1_p2, data_M1_p3, _, data_M1_p5, data_M1_p6 = readMultiBandGeotiff(sFile_M1_p)
    data_hous += (data_M1 * data_M1_p1 + data_M1 * data_M1_p2)/100
    data_ind += (data_M1*data_M1_p6)/100
    data_serv += (data_M1*data_M1_p3 + data_M1*data_M1_p5)/100
    
if os.path.isfile(sFile_M2_p):
    data_M2_p1, data_M2_p2, data_M2_p3, _, data_M2_p5, data_M2_p6 = readMultiBandGeotiff(sFile_M2_p)
    data_hous += (data_M2 * data_M2_p1 + data_M2 * data_M2_p2)/100
    data_ind += (data_M2*data_M2_p6)/100
    data_serv += (data_M2*data_M2_p3 + data_M2*data_M2_p5)/100

if os.path.isfile(sFile_W1_p):
    data_W1_p1, data_W1_p2, data_W1_p3, _, data_W1_p5, data_W1_p6 = readMultiBandGeotiff(sFile_W1_p)
    data_hous += (data_W1 * data_W1_p1 + data_W1 * data_W1_p2)/100
    data_ind += (data_W1*data_W1_p6)/100
    data_serv += (data_W1*data_W1_p3 + data_W1*data_W1_p5)/100

if os.path.isfile(sFile_T1_p):
    data_T1_p1, data_T1_p2, data_T1_p3, _, data_T1_p5, data_T1_p6 = readMultiBandGeotiff(sFile_T1_p)
    data_hous += (data_T1 * data_T1_p1 + data_T1 * data_T1_p2)/100
    data_ind += (data_T1*data_T1_p6)/100
    data_serv += (data_T1*data_T1_p3 + data_T1*data_T1_p5)/100

hous_tot = np.nansum(data_hous)
log_print(f'total housing: [{hous_tot}]')
ind_tot = np.nansum(data_ind)
log_print(f'total industrial: [{ind_tot}]')
serv_tot = np.nansum(data_serv)
log_print(f'total service: [{serv_tot}]')

writeGeotiffSingleBand(sFile_hous, geotransform, geoproj, data_hous)
writeGeotiffSingleBand(sFile_ind, geotransform, geoproj, data_ind)
writeGeotiffSingleBand(sFile_serv, geotransform, geoproj, data_serv)




