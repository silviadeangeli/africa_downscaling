__author__ = 'silvia'

from Geotiff_Silvia import readFile, readFile_withNoData, writeGeotiffSingleBand, readFileBand
import numpy as np
from os.path import join
import os
from rasterRegrid import rasterRegrid
import matplotlib.pylab as plt
from africa_tools import log_print

            
mask = '/Users/silvia/Documents/AFRICA_DATA/Tanzania/regridded_hrsl/hrsl_tza_buiA_mask_90m_regridded.tif'
path_valfis = '/Users/silvia/Documents/AFRICA_DATA/Tanzania/VALFIS_edu'
path_output =  '/Users/silvia/Documents/AFRICA_DATA/Tanzania/TZ_criEDU_GAR_downscaled_VALFIS_90m_20181016'

xsize, ysize, geotransform, geoproj, data_mask   = readFile(mask)

materials = [
    'C1',
    'M1',
    'M2',
    'W1',
    'T1'
]


for m in materials:
    output_file = join(path_output, f'{m}.tif')
    m_file = join(path_valfis, f"{m}.tif")
    if os.path.isfile(m_file):
        print(f'{m} m_file found')
        xsize, ysize, geotransform, geoproj, data_m   = readFile(m_file)
        sum1 = np.nansum(data_m)
        data_m_new = data_m*data_mask
        sum2 = np.nansum(data_m_new)
        lost_p = (sum1-sum2)*100/sum1
        print(f'{m} lost percentage: {lost_p} %') 
        writeGeotiffSingleBand(output_file, geotransform, geoproj, data_m_new)