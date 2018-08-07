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

sFile = "E:/africa_downscaling/ao/ao_industrial.tif"

[xsize, ysize, geotransform, geoproj, data]   = readFile(sFile)

tot = np.sum(data)

log_print(f'total: [{tot}]')



