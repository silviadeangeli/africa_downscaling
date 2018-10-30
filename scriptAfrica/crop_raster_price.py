__author__ = 'lauro'


from Geotiff import readFile, writeFile, writeGeotiffSingleBand, readFile_withNoData
from collections import Counter
from osgeo import gdal, gdalconst, ogr
import os, sys, csv, getopt
import matplotlib.pylab as plt
import numpy as np
from rasterRegrid import rasterRegrid
import datetime
import random

#os.system(gdal_translate -a_srs EPSG:32736 -outsize 25% 25% -of GTiff -co COMPRESS=DEFLATE /Users/lauro/Documents/PROJECTS/FP7_H2020/RASOR_2012/MALAWI/crop/annual_crop_MW.tif /Users/lauro/Documents/PROJECTS/FP7_H2020/RASOR_2012/MALAWI/crop/annual_crop_60m.tif)

def log_print(s):
    print(datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S') + " - " + s)

def non_unique (array):
    log_print ("Making unique values...")
    uniq_val = []
    for elem in array:
        if elem != 0 and not np.isnan (elem):
            #print(elem)
            elem = elem + random.random() * 0.01
        uniq_val.append(elem)
    return uniq_val

def plot_image(data):
    fig = plt.figure()
    cax = plt.imshow(data, cmap='Blues')
    plt.colormaps()
    #plt.clim(0,400)
    cbar = fig.colorbar(cax, orientation='vertical')
    #cbar.ax.set_yticklabels(['< -1', '0', 1, 2,'> 10'])# vertically oriented colorbar
    plt.show()

def get_crop (sFileMapSpamLowRes):
    return (os.path.basename(sFileMapSpamLowRes).split('_total')[0]).split('production_')[1]

def get_country (sFileMask):
    return (os.path.basename(sFileMask).split('_crops')[0])

def readPrice (sFileIn):

    #a1sInputName=[]
    #a1sInputValue=[]
    cropPrice = {}

#np.loadtxt(inputfilename, comments='#', delimiter=' ',dtype='string', ndmin=2)

    with open(sFileIn, 'r') as filePar:

        for line in filePar:
            try:
                if "#" in line:
                    (line1,line2) = line.split('#',1)
                else:
                    line1=line
                #print line1
                (sVarName,sVarValue) = line1.replace(' ','\t').split('\t',1)
                if sVarName.strip() =='' and sVarValue.strip()=='':
                    continue
                else:
                    cropPrice[sVarName.strip()] = sVarValue.strip()
                    #a1sInputName.append(sVarName.strip())
                    #a1sInputValue.append(sVarValue.strip())

            except:
                print ('File: '+sFileIn)
                print ('Skipping line: '+line)

    #return a1sInputName, a1sInputValue
    return cropPrice

def downscale_crop (sFileMapSpamLowResTot,sFileMask, sFileClipper, sFilePrice):

    sCropType = get_crop (sFileMapSpamLowResTot)

    sCountry = get_country (sFileMask)

    print ('\n')
    log_print('Crop type: '+sCropType)

    #making dir out 
    dir_out = os.path.join(os.path.dirname(sFileMapSpamLowResTot),sCountry+"_outputs")
    if not os.path.exists(dir_out):
        os.mkdir(dir_out)

    sFileMapSpamLowRes = os.path.join(dir_out, os.path.basename(sFileMapSpamLowResTot.split('.')[0] + "_"+sCountry+".tif"))

    warp = '''gdalwarp -co "COMPRESS=DEFLATE" -dstnodata {noDataValue} -crop_to_cutline -overwrite -cutline {clipper} "{infile}" "{outfile}"'''.format(
                    noDataValue=-9999, clipper=sFileClipper, infile=sFileMapSpamLowResTot, outfile=sFileMapSpamLowRes)
    #print (warp)
    os.system (warp)

    [xsize_orig, ysize_orig, geotransform, geoproj, data_low] = readFile_withNoData(sFileMapSpamLowRes)
    log_print ("Original MAPSPAM map size: ("+str(ysize_orig)+","+str(xsize_orig)+") (rows,cols)")

    #TODO verify
    total_exposure = np.nansum(data_low);
    
    #making unique values
    data_low=np.reshape(non_unique (data_low.ravel()), data_low.shape)

    writeGeotiffSingleBand(sFileMapSpamLowRes, geotransform, geoproj, data_low)

    if bDisplayImages: plot_image(data_low)

    log_print("Total original production for crop type="+sCropType+": %.2f" % np.nansum(data_low))

    del data_low

    sFileMapSpamHiRes  = os.path.join(dir_out, sFileMapSpamLowRes.split('.')[0] + "_regrid.tif")

    match_geotrans, match_proj = rasterRegrid(sFileMapSpamLowRes, sFileMask, sFileMapSpamHiRes ,"nearest")

    [xsize, ysize, geotransform_high, geoproj_high, data_high]   = readFile_withNoData(sFileMapSpamHiRes)
    [xsize, ysize, geotransform_high, geoproj_high, data_mask]   = readFile_withNoData(sFileMask)

    data_high [np.isnan (data_high)] = 0
    data_mask [np.isnan (data_mask)] = 0

    data_high_masked = data_high * data_mask

    if bDisplayImages:
        plot_image(data_mask)
        plot_image(data_high)
        plot_image(data_high_masked)

    data_high_masked [np.isnan (data_high_masked) ] = 0

    #debug
    #plot_image(data_high_masked)

    log_print ("Starting counting...")
    dictCounter = Counter(data_high_masked.ravel())
    log_print ("Data counted")

    data_mask_spam=np.copy(data_high)
    ratio = abs ((geotransform[1]*geotransform[5]) / (geotransform_high[1]*geotransform_high[5])) #number of high res pixels in low res pixel

    i=0
    data_gar_excluded = []
    #non posso leggere data low perche cambia nel salvarli
    data_high_unique=np.unique(data_high.ravel())
    for key in data_high_unique:
            if key not in dictCounter.keys() and key != 0:
                dictCounter[key] = ratio
                data_mask_spam[data_mask_spam==key]=-9999
                data_gar_excluded.append(key)
    excluded_exposure = sum(data_gar_excluded)
    log_print ("MAPSPAM production not overlapping mask: %.2f (%.1f %%)" % (excluded_exposure , excluded_exposure/total_exposure*100))

    #building mask (union of the original mask + non-zero values from MapSpam)
    data_mask_spam[data_mask_spam!=-9999]=0
    data_mask_spam[data_mask_spam==-9999]=1
    data_mask = data_mask_spam + data_mask
    data_mask[data_mask>0] = 1
    data_high_masked = data_high * data_mask

    del data_mask, data_high, data_mask_spam

    if bVerbose: log_print ("Counter length: "+str(dictCounter.__len__()))
    if bVerbose: log_print ("Unique length: "+str(len(np.unique(data_high_masked))))

    for key,value  in dictCounter.items():
        data_high_masked [data_high_masked==key] = key/value
        if bVerbose: print ('Amount of pixel for crop' + str(key) + ': ' + str(value))

    if bDisplayImages:
        fig = plt.figure()
        cax = plt.imshow(data_high_masked[:, :])
        plt.colormaps()
        cbar = fig.colorbar(cax, orientation='vertical')  # vertically oriented colorbar
        #cbar.ax.set_yticklabels(['< -1', '0', 1, 2,'> 10'])
        plt.show()

    log_print("Total downscaled production value for crop type="+sCropType+": %.2f" % np.nansum(data_high_masked))

    #from production to economic value
    cropPrice = readPrice(sFilePrice)
    data_high_masked = data_high_masked * float(cropPrice[sCropType])

    log_print ('Price of '+ sCropType +' = '+cropPrice[sCropType]+' [USD]')

    sFileDownscaled = os.path.join(dir_out, sCountry+"_"+sCropType+"_HighRes.tif")

    description = 'Production value for ' + sCropType + ' [USD]'
    bandMetadata = {}
    bandMetadata[('unit', '1')] = 'USD'
    iNbands=1;

    writeFile(sFileDownscaled, match_geotrans, match_proj, data_high_masked, bandMetadata, description, {'Crop type': sCropType},iNbands)
    log_print ("Exposure downscaled. Final result save in file: " + sFileDownscaled)
    
################################
############# START ############
################################

bDisplayImages = 0
bVerbose = 0 

sFileMask = '/Users/silvia/Documents/AFRICA_DATA/Guinea_Bissau/GW_agrM_ESA_90m_20180628.tif'
sFileBoundary = '/Users/silvia/Documents/AFRICA_DATA/Guinea_Bissau/GW_bou_GAUL_admin0_20180901/GW_bou_GAUL_admin0_20180901.shp'
sFilePrice = '/Users/silvia/Documents/AFRICA_DATA/Guinea_Bissau/CropsPrice.txt'

cropPrice = readPrice (sFilePrice)

# crop production directory
dirIn = '/Users/silvia/Documents/AFRICA_DATA/Guinea_Bissau/Crops'

for file in os.listdir(dirIn):
    if  file.startswith('spam') and file.endswith('tiff'):
        downscale_crop(os.path.join(dirIn,file),sFileMask, sFileBoundary, sFilePrice)
