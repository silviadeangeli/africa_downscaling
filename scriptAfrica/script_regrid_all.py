__author__ = 'lauro'

import os
from rasterRegrid import rasterRegrid

origineDir = '/share/africa/Global_datasets/Builtup/WRL_buiA_GHSL_38m_cuts_20180305/Ritagli_subnazionali'
#origineDir = '/Users/lauro/Downloads'
matchDir = '/share/africa/Global_datasets/Builtup/AFR_buiA_GUF_100m_cuts_20180312'

method = 'nearest'

for f in os.listdir(origineDir):

    if f.endswith('.tif'):
        sFileMatch = os.path.join(matchDir,(os.path.basename(f)))
        #print (sFileMatch)
        filename, ext = os.path.splitext(os.path.basename(f))
        print (filename)
        sFileOut = os.path.join(origineDir + '_90m_'+ method , filename + '_90m_'+ method +'.tif')
        #print (origineDir + '/' + f,sFileMatch, sFileOut)

        #match_geotrans,match_proj= rasterRegrid (origineDir+'/'+f ,sFileMatch, sFileOut, method)