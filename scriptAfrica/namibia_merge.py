__author__ = 'silvia'

from Geotiff_Silvia import readFile, readFile_withNoData, writeGeotiffSingleBand, readFileBand
import numpy as np
from os.path import join
import os
from rasterRegrid import rasterRegrid
import matplotlib.pylab as plt
from africa_tools import log_print
from toolbox_mirko.africa_tools import readMultiBandGeotiff


#########################################
############## FILE PATH  ###############
#########################################

path = '/Users/silvia/Documents/AFRICA_DATA/Guinea_Bissau'
path_valfis1 = join(path, 'valfis_20181025')
path_valfis2 = join(path, 'valfis_FILTERED_20181025')
mask = join(path, 'integrazione_buiA/GW_buiA_GHSL_GUF_LC_OSM_wp3_90m_integration.tif')

# outputs
path_output = join(path, 'VALFIS_merge_outputs')

materials = [
    'C1',
    'M1',
    'M2',
    'W1',
    'T1'
]

xsize, ysize, geotransform, geoproj, mask   = readFile(mask)
mask_false = np.zeros((ysize, xsize))
mask_false [mask == 1 ] =0 
mask_false [mask != 1 ] =1


for m in materials:
    output_file = join(path_output, f'{m}.tif')
    m_file1 = join(path_valfis1, f"{m}.tif")
    if os.path.isfile(m_file1):
        print(f'{m} 1st file found')
        xsize, ysize, geotransform, geoproj, data_m1   = readFile(m_file1)
    
        m_file2 = join(path_valfis2, f"{m}.tif")         
        if os.path.isfile(m_file2):
            print(f'{m} 2nd file found')
            xsize, ysize, geotransform, geoproj, data_m2   = readFile(m_file2)
    
            data_m_tot = data_m1*mask + data_m2*mask_false
            print(f'{m} data_tot_calculated')
            writeGeotiffSingleBand(output_file, geotransform, geoproj, data_m_tot)
            print(f'{m} file saved')
        else:
            print(f'second {m} file not found')
    else:
        print(f'{m} file not found')

        
        
