__author__ = 'silvia'

from Geotiff_Silvia import readFile, readFile_withNoData, writeGeotiffSingleBand, readFileBand
import numpy as np
from os.path import join
import os
from rasterRegrid import rasterRegrid
import matplotlib.pylab as plt
from africa_tools import log_print
from high_res_mask import write_hres_mask
from raster_sum import raster_sum

            
#inputs
path = '/Users/silvia/Documents/AFRICA_DATA/Guinea_Bissau/'
path_valfis_folder = join(path, 'VALFIS_merge_outputs')
path_pop = join(path, 'regrid_outputs/GW_pop_WP_90m.tif')
pop_2016 = 1815698
pop_2050 = 2488000
    
# outputs
path_output_mask = join(path, 'GW_buiA_mask_20181025.tif')
path_output_pop = join(path, 'GW_pop_in_buiA_20181025.tif')
path_output_pop_2016 = join(path, 'GW_pop_2016_in_buiA_20181025.tif')
path_output_pop_2050 = join(path, 'GW_pop_2050_in_buiA_20181025.tif')
    
materials = [
    'C1',
    'M1',
    'M2',
    'W1',
    'T1'
 ]
    
write_hres_mask(path_valfis_folder, path_output_mask, materials)
xsize, ysize, geotransform, geoproj, data_mask   = readFile(path_output_mask)
xsize, ysize, geotransform, geoproj, data_pop   = readFile(path_pop)
    
data_pop_on_mask = data_mask*data_pop
writeGeotiffSingleBand(path_output_pop, geotransform, geoproj, data_pop_on_mask)
    
tot_pop = raster_sum(data_pop_on_mask)

factor_2016 = pop_2016/tot_pop
data_pop_2016 = data_pop_on_mask*factor_2016
writeGeotiffSingleBand(path_output_pop_2016, geotransform, geoproj, data_pop_2016)
#check
raster_sum(data_pop_2016)

factor_2050 = pop_2050/tot_pop
data_pop_2050 = data_pop_on_mask*factor_2050
writeGeotiffSingleBand(path_output_pop_2050, geotransform, geoproj, data_pop_2050)
#check
raster_sum(data_pop_2050)


