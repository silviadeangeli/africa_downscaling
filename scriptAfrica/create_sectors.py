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
path_valfis = join(path, 'VALFIS_merge_outputs')
path_p = join(path, 'use_percentage')

# outputs
sFile_hous = join(path, 'GW_buiA_housing_sector_20181025.tif')
sFile_ind = join(path, 'GW_buiA_productive_sector_20181025.tif')
sFile_serv = join(path, 'GW_buiA_service_sector_20181025.tif')

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
        xsize, ysize, geotransform, geoproj, data_m   = readFile(m_file)

        if not init:
            data_hous = np.zeros((ysize, xsize))
            data_ind = np.zeros((ysize, xsize))
            data_serv = np.zeros((ysize, xsize))
            init = True
        
        p_file = join(path_p, f'{m}.tiff')
        if os.path.isfile(p_file):
            print(f'{m} p_file found')            
            data_p1, data_p2, data_p3, data_p4, data_p5, data_p6 \
                = readMultiBandGeotiff(p_file)

            data_hous += (data_m * data_p1 + data_m * data_p2)/100
            data_ind  += (data_m * data_p6)/100
            data_serv += (data_m * data_p3 + data_m * data_p5)/100


hous_tot = np.nansum(data_hous)
log_print(f'total housing: [{hous_tot}]')
ind_tot = np.nansum(data_ind)
log_print(f'total industrial: [{ind_tot}]')
serv_tot = np.nansum(data_serv)
log_print(f'total service: [{serv_tot}]')

writeGeotiffSingleBand(sFile_hous, geotransform, geoproj, data_hous)
writeGeotiffSingleBand(sFile_ind, geotransform, geoproj, data_ind)
writeGeotiffSingleBand(sFile_serv, geotransform, geoproj, data_serv)
