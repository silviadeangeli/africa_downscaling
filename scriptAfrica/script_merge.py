__author__ = 'lauro'

import os
from Geotiff import writeGeotiffSingleBand, readFile

dir1 = '/share/africa/Global_datasets/Builtup/AFR_buiA_GUF_100m_cuts_20180312'
dir2 = '/share/africa/Global_datasets/Builtup/WRL_buiA_GHSL_100m_nearest_20180312'
#origineDir = '/Users/lauro/Downloads'

dirOut = '/share/africa/Global_datasets/Builtup/AFR_buiA_merged_100m_20180312'

for f in os.listdir(dir1):

    if f.endswith('.tif'):

        #print (f)
        filename, ext = os.path.splitext(f)
        #print filename

        #if filename == "UR_Tanzania_tanga":
        [xsize, ysize, geotransform, geoproj, data1]   = readFile(os.path.join(dir1, f))
        [xsize, ysize, geotransform, geoproj, data2]   = readFile(os.path.join(dir2, filename+"_100m_nearest.tif"))
        #print(filename+"_100m_nearest.tif")
        sFileOut = os.path.join(dirOut, filename + '_GHSL_GUF_merged.tif')
        print ("writing... "+ sFileOut)
        writeGeotiffSingleBand(sFileOut,geotransform,geoproj,data1+data2)
