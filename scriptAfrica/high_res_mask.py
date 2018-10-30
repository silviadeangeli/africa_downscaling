__author__ = 'silvia'

from Geotiff_Silvia import readFile, readFile_withNoData, writeGeotiffSingleBand, readFileBand
import numpy as np
from os.path import join
import os
from rasterRegrid import rasterRegrid
import matplotlib.pylab as plt
from africa_tools import log_print

            
def write_hres_mask(path_valfis, path_output, materials):
    init = False
    for m in materials:
        m_file = join(path_valfis, f"{m}.tif")
        if os.path.isfile(m_file):
            print(f'{m} m_file found')
            xsize, ysize, geotransform, geoproj, data_m   = readFile(m_file)

            if not init:
                data_tot = np.zeros((ysize, xsize))
                init = True

            data_tot += data_m
    if not init:
        print("no file was loaded")
        return
    
    data_tot [data_tot > 0] = 1
    data_tot [data_tot <= 0] = 0
    writeGeotiffSingleBand(path_output, geotransform, geoproj, data_tot)

if __name__ == '__main__':
    path = '/Users/silvia/Documents/AFRICA_DATA/Guinea_Bissau'
    path_valfis = join(path, 'valfis')
    # outputs
    path_output = join(path, 'gw_buiA_mask.tif')
    materials = [
        'C1',
        'M1',
        'M2',
        'W1',
        'T1'
    ]


    write_hres_mask(path_valfis, path_output, materials)