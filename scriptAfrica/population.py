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
path = '/Users/silvia/Documents/AFRICA_DATA/Botswana/'
path_pop = join(path, 'BW_gdp_2016_GAR_normalizedWP_90m_20180730.tif')
pop_2016 = 15649000000
pop_2050 = 88597000000
    
# outputs
path_output_pop_2016 = join(path, 'BW_gdp_2016_GAR_normalizedWP_90m_20181024.tif')
path_output_pop_2050 = join(path, 'BW_gdp_2050_GAR_normalizedWP_90m_20181024.tif')
    

xsize, ysize, geotransform, geoproj, data_pop   = readFile(path_pop)
tot_pop = raster_sum(data_pop)

factor_2016 = pop_2016/tot_pop
data_pop_2016 = data_pop*factor_2016
writeGeotiffSingleBand(path_output_pop_2016, geotransform, geoproj, data_pop_2016)
#check
raster_sum(data_pop_2016)

factor_2050 = pop_2050/tot_pop
data_pop_2050 = data_pop*factor_2050
writeGeotiffSingleBand(path_output_pop_2050, geotransform, geoproj, data_pop_2050)
#check
raster_sum(data_pop_2050)


