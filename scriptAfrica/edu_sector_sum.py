__author__ = 'silvia'

from Geotiff_Silvia import readFile, readFile_withNoData, writeGeotiffSingleBand, readFileBand
import numpy as np
from os.path import join
import os
from rasterRegrid import rasterRegrid
import matplotlib.pylab as plt
from africa_tools import log_print

            
output_file =  '/Users/silvia/Documents/AFRICA_DATA/Botswana/BW_education_sector_20181018.tif'
path_valfis = '/Users/silvia/Documents/AFRICA_DATA/Botswana/BW_criEDU_GAR_downscaled_VALFIS_90m_20181018'

materials = [
    'C1',
    'M1',
    'M2',
    'W1',
    'T1'
]

init = False 
for m in materials:
    m_file = join(path_valfis, f"{m}.tif")
    if os.path.isfile(m_file):
        print(f'{m} m_file found')
        xsize, ysize, geotransform, geoproj, data_m = readFile(m_file)
        if not init:
            data_edu = np.zeros((ysize, xsize)) 
            init = True
    data_edu += data_m
writeGeotiffSingleBand(output_file, geotransform, geoproj, data_edu)