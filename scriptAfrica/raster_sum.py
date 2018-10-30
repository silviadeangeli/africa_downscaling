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


def raster_sum(data):
    
    data [data <= 0] = 0

    tot = np.nansum(data)

    log_print(f'total: [{tot}]')
    return tot

if __name__ == '__main__':
    
    sFile = '/Users/silvia/Desktop/GDP/TZ_gdp_2050_90m_normalizedWP_20181001.tif'
    [xsize, ysize, geotransform, geoproj, data]   = readFile(sFile)
    raster_sum(data)
    
    
