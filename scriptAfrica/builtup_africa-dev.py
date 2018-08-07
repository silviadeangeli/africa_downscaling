__author__ = 'lauro'

from Geotiff import readFile, readFile_withNoData, writeFile, writeGeotiffSingleBand
from collections import Counter
from osgeo import gdal, gdalconst, ogr
import os, sys, csv, getopt
import matplotlib.pylab as plt
import numpy as np
from rasterRegrid import rasterRegrid
import datetime
import random
from progressbar import ProgressBar as PB

# default nan value
NAN_VALUE = np.nan

#os.system(gdal_translate -a_srs EPSG:32736 -outsize 25% 25% -of GTiff -co COMPRESS=DEFLATE /Users/lauro/Documents/PROJECTS/FP7_H2020/RASOR_2012/MALAWI/crop/annual_crop_MW.tif /Users/lauro/Documents/PROJECTS/FP7_H2020/RASOR_2012/MALAWI/crop/annual_crop_60m.tif)

def log_print(s):
    print(datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S') + " - " + s)

def displayHelp():

    print ('\nScript to downscale GAR input into a mask map')
    print ('Options:')
    print ('         -v | --verbose                          print debug output on the screen')
    print ('         -h | --help                             display this help')
    print ('         -i | --input                            input file or single population value')
    print ('         -m | --mask                             file name of high res mask')
    print ('         -e | --exposure                         exposure type: VALFIS or VALHUM')
    print ('   -n VALUE | --nan VALUE                        replaces nan values with VALUE')
    print ('         -f | --filter                           attribute filter name (filter), ALL takes all attributes')
    print ('\nUse -f ALL is VALHUM is selected')
    print ('\nExample: python builtup_africa.py -i input.shp -m mask.tif -e [VALFIS|VALHUM] -f BUILDING_TYPE')


def createBuffer(inputfn, outputBufferfn, bufferDist):
    inputds = ogr.Open(inputfn)
    inputlyr = inputds.GetLayer()

    shpdriver = ogr.GetDriverByName('ESRI Shapefile')
    if os.path.exists(outputBufferfn):
        shpdriver.DeleteDataSource(outputBufferfn)
    outputBufferds = shpdriver.CreateDataSource(outputBufferfn)
    bufferlyr = outputBufferds.CreateLayer(outputBufferfn, geom_type=ogr.wkbPolygon)
    featureDefn = bufferlyr.GetLayerDefn()

    for feature in inputlyr:
        ingeom = feature.GetGeometryRef()
        geomBuffer = ingeom.Buffer(bufferDist)

        outFeature = ogr.Feature(featureDefn)
        outFeature.SetGeometry(geomBuffer)
        bufferlyr.CreateFeature(outFeature)
        outFeature = None

def non_unique (list):
    log_print ("Making unique values...")
    uniq_val = []
    for elem in list:
        if elem != 0:
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

def downscale_pop(sFileGARPoint, sFileMask, sExposureType):
    #sBuildingType = "UFB"
    #sBuildingType = "UCB"

    bDisplayImages = False

    #sFileGARPoint = '/Users/lauro/Desktop/EU_ACP_UNISDR_Africa/DATA/UR_Tanzania/AFR_buiA_GAR_5km_20180313/africa_tza.shp'
    #sFileMask = '/Users/lauro/Desktop/EU_ACP_UNISDR_Africa/DATA/UR_Tanzania/AFR_buiA_GUF_100m_nations_20180315/UR_Tanzania.tif'

    #sFileGARPoint = '/Users/lauro/Desktop/EU_ACP_UNISDR_Africa/DATA/Swaziland/AFR_buiA_GAR_5km_20180313/africa_swz.shp'
    #sFileMask = '/Users/lauro/Desktop/EU_ACP_UNISDR_Africa/DATA/Swaziland/AFR_buiA_GUF_100m_nations_20180315/Swaziland.tif'

    #os.system('ogr2ogr -overwrite -f "ESRI Shapefile" -where "curve=\''+sBuildingType+'\'" ' + sFileGARPointCurveType + ' ' + sFileGARPoint)

    inDriver = ogr.GetDriverByName("ESRI Shapefile")
    inDataSource = inDriver.Open(sFileGARPoint, 0)
    inLayer = inDataSource.GetLayer()
    inLayerDefn = inLayer.GetLayerDefn()
    log_print("Reading GAR shape file...")
    log_print("Total features: "+str(inLayer.GetFeatureCount()))

    #inLayer.SetAttributeFilter("curve = \''+sBuildingType+'\'")

    (fXmin, fXmax, fYmin, fYmax) = inLayer.GetExtent()
    iXLowRes = 0.04167
    iYLowRes = 0.04167
    fXmin = fXmin - iXLowRes / 2
    fYmin = fYmin - iYLowRes / 2
    fXmax = fXmax + iXLowRes / 2
    fYmax = fYmax + iYLowRes / 2

    sFileGARPointBase= os.path.basename(sFileGARPoint)
    dir_out = os.path.join(os.path.dirname(sFileGARPoint),sFileGARPointBase.split('.')[0]+"_"+"outputs")
    if not os.path.exists(dir_out):
        os.mkdir(dir_out)

    sFileGARPointPop =  os.path.join(dir_out, sFileGARPointBase.split('.')[0]+"_"+sExposureType+".shp")
    outDriver = ogr.GetDriverByName( "ESRI Shapefile")
    outDataSource = outDriver.CreateDataSource(sFileGARPointPop)
    proj = inLayer.GetSpatialRef()
    ##Creating the layer with its fields
    outLayer = outDataSource.CreateLayer(sFileGARPointPop, proj, ogr.wkbPoint )
    field_exposure = ogr.FieldDefn(sExposureType, ogr.OFTReal)
    outLayer.CreateField(field_exposure)

    pop = []
    geoms = []

    for feature in inLayer:

            if feature.GetGeometryRef().GetPoint() not in geoms:
                geoms.append(feature.GetGeometryRef().GetPoint())
                pop.append(feature.GetField(sExposureType))

            else:
                id = geoms.index(feature.GetGeometryRef().GetPoint())
                pop[id]=pop[id]+feature.GetField(sExposureType)

    for i in range(len(geoms)):
        outFeature= ogr.Feature(outLayer.GetLayerDefn())
        outFeature.SetField(sExposureType, pop[i])
        point = ogr.Geometry(ogr.wkbPoint)
        point.AddPoint(*geoms[i])
        outFeature.SetGeometry(point)
        outLayer.CreateFeature(outFeature)

    total_exposure = sum (pop)
    log_print("Total unique features in outLayer: "+str(outLayer.GetFeatureCount()))
    log_print("Total original exposure value (for curve type="+sCurveType+", exposure type="+sExposureType+"): %.2f" % total_exposure)

    inLayer.ResetReading()
    inDataSource = None
    outDataSource = None

#    for i in range (0, len(geoms)):

    log_print("Converting points to raster...")

    sFileGARRasterEXPLowRes = os.path.join(dir_out, sFileGARPointBase.split('.')[0] + "_"+sExposureType+".tif")
    if bVerbose: print('gdal_rasterize -co compress=DEFLATE -a '+sExposureType+' -l ' + os.path.basename(sFileGARPointPop).split('.')[0] + ' -tr ' + str(iXLowRes) + ' ' + str(iYLowRes) + ' -te ' + str(fXmin) + ' ' + str(fYmin) + ' ' + str(fXmax) + ' ' + str(fYmax) + ' ' + sFileGARPointPop + ' ' +sFileGARRasterEXPLowRes)
    os.system('gdal_rasterize -co compress=DEFLATE -a '+sExposureType+' -l ' + os.path.basename(sFileGARPointPop).split('.')[0] + ' -tr ' + str(iXLowRes) + ' ' + str(iYLowRes) + ' -te ' + str(fXmin) + ' ' + str(fYmin) + ' ' + str(fXmax) + ' ' + str(fYmax) + ' ' + sFileGARPointPop + ' ' +sFileGARRasterEXPLowRes)

    # if bVerbose: print('gdal_rasterize -co compress=DEFLATE -a '+sExposureType+' -where "curve=\''+sBuildingType+'\'" -l ' + os.path.basename(sFileGARPoint).split('.')[0] + ' -tr ' + str(iXLowRes) + ' ' + str(iYLowRes) + ' -te ' + str(fXmin) + ' ' + str(fYmin) + ' ' + str(fXmax) + ' ' + str(fYmax) + ' ' + sFileGARPoint + ' ' +sFileGARRasterEXPLowRes)
    # os.system('gdal_rasterize -co compress=DEFLATE -a '+sExposureType+' -where "curve=\''+sBuildingType+'\'" -l ' + os.path.basename(sFileGARPoint).split('.')[0] + ' -tr ' + str(iXLowRes) + ' ' + str(iYLowRes) + ' -te ' + str(fXmin) + ' ' + str(fYmin) + ' ' + str(fXmax) + ' ' + str(fYmax) + ' ' + sFileGARPoint + ' ' +sFileGARRasterEXPLowRes)

    [xsize_orig, ysize_orig, geotransform, geoproj, data_low]   = readFile(sFileGARRasterEXPLowRes)

    log_print ("Original GAR map size: ("+str(ysize_orig)+","+str(xsize_orig)+") (rows,cols)")

    uniq = non_unique (data_low.ravel().tolist())

    data_low=np.reshape(uniq, data_low.shape)

    writeGeotiffSingleBand(sFileGARRasterEXPLowRes, 
                    geotransform, geoproj, 
                    data_low, nan_value=NAN_VALUE)

    if bDisplayImages: plot_image(data_low)
    # DEBUG
    #b = np.unique(data_low.ravel()).tolist()

    del data_low

    sFileGARRasterEXPHiRes  = os.path.join(dir_out, sFileGARPointBase.split('.')[0] + "_"+sCurveType+"_"+sExposureType+"_regrid.tif")

    match_geotrans, match_proj = rasterRegrid(sFileGARRasterEXPLowRes, sFileMask, sFileGARRasterEXPHiRes ,"nearest")

    [xsize, ysize, geotransform_high, geoproj, data_high]   = readFile(sFileGARRasterEXPHiRes)
    [xsize, ysize, geotransform_high, geoproj, data_mask]   = readFile(sFileMask, fix_nan_value=255)

    ## DEBUG
    #c = np.unique(data_high.ravel()).tolist()
    #print (len(b),len(c))

    data_high_masked = data_high * data_mask

    if bDisplayImages:
        plot_image(data_mask)
        plot_image(data_high)
        plot_image(data_high_masked)

    log_print ("Starting counting...")
    dictCounter = Counter(data_high_masked.ravel().tolist())
    log_print ("Data counted")

    data_mask_gar=np.copy(data_high)
    ratio = abs ((geotransform[1]*geotransform[5]) / (geotransform_high[1]*geotransform_high[5])) #number of high res pixels in low res pixel

    data_gar_excluded = []
    for key in data_high.ravel().tolist():    # for name, age in list.items():  (for Python 3.x)
        if key not in dictCounter.keys() and key != 0 :
            dictCounter[key] = ratio
            data_mask_gar[data_mask_gar==key]=-9999
            data_gar_excluded.append(key)
    excluded_exposure = sum(data_gar_excluded)
    log_print ("GAR exposure value not overlapping mask: %.2f (%.1f %%)" % (excluded_exposure , excluded_exposure/total_exposure*100))
    del data_gar_excluded

    #building mask (union of the original mask + non-zero values form GAR)
    data_mask_gar[data_mask_gar!=-9999]=0
    data_mask_gar[data_mask_gar==-9999]=1
    data_mask = data_mask_gar + data_mask
    data_mask[data_mask>0] = 1
    data_high_masked = data_high * data_mask

    del data_mask, data_high, data_mask_gar

    if bVerbose: log_print ("Counter length: "+str(dictCounter.__len__()))
    if bVerbose: log_print ("Unique length: "+str(len(np.unique(data_high_masked))))

    for key,value  in dictCounter.items():
        data_high_masked [data_high_masked==key] = key/value
        if bVerbose: log_print ('Amount of pixel for population ' + str(key) + ': ' + str(value))

    if bDisplayImages:
        fig = plt.figure()
        cax = plt.imshow(data_high_masked[:, :])
        plt.colormaps()
        cbar = fig.colorbar(cax, orientation='vertical')  # vertically oriented colorbar
        #cbar.ax.set_yticklabels(['< -1', '0', 1, 2,'> 10'])
        plt.show()

    sFileDownscaled = os.path.join(dir_out, sFileGARPointBase.split('.')[0] + "_"+sCurveType+"_"+sExposureType+"_HighRes.tif")
    bandMetadata = {}
    bandMetadata[('unit', '1')] = 'people'
    #bandMetadata[('unit', '2')] = 'USD/tons'

    if sExposureType == "VALHUM":
        description = 'VALHUM in ' + sCurveType + ' (downscaling of GAR2015 data)'
    elif sExposureType == "VALFIS":
        description = 'VALFIS in ' + sCurveType + ' (downscaling of GAR2015 data)'

    iNbands = 1
    
    writeFile(sFileDownscaled, match_geotrans, match_proj, 
                data_high_masked, bandMetadata, description, 
                {'Building_type': sCurveType}, iNbands, 
                nan_value=NAN_VALUE)

    log_print("Total downscaled exposure value (for curve type="+sCurveType+", exposure type="+sExposureType+"): %.2f" % np.sum(data_high_masked))
    log_print ("Exposure downscaled. Final result save in file: " + sFileDownscaled)

    
    
def single_value_to_population(value, file_mask):
    """
    read mask file, calculate population count and rescale population/pixel value according to 'value'
    writes population_{value}_highres.tiff file
    """
    
    dir_out = os.path.join(os.path.dirname(file_mask), "outputs")
    os.makedirs(dir_out, exist_ok=True)
    
    xsize, ysize, geotransform, geoproj, mask_data = readFile(file_mask)
    pop_count = np.sum(mask_data)
    scaling_coeff = value/pop_count
    new_data = mask_data * scaling_coeff
    
    file_downscaled = os.path.join(
                                dir_out, 
                                "population_" + str(value) + "_HighRes.tif"
    )
    
    description = 'Downscaled value'
    bandMetadata = {}
    bandMetadata[('unit', '1')] = 'people'
    writeFile(file_downscaled, geotransform, 
            geoproj, new_data, bandMetadata, 
            description, {}, 1, 
            nan_value=NAN_VALUE)


def downscale(sFileGARPoint, sFileMask, sBuildingType, sExposureType):
    #sBuildingType = "UFB"
    #sBuildingType = "UCB"

    bDisplayImages = False

    print (" ")
    log_print ("Curve: " + sBuildingType)

    #sFileGARPoint = '/Users/lauro/Desktop/EU_ACP_UNISDR_Africa/DATA/UR_Tanzania/AFR_buiA_GAR_5km_20180313/africa_tza.shp'
    #sFileMask = '/Users/lauro/Desktop/EU_ACP_UNISDR_Africa/DATA/UR_Tanzania/AFR_buiA_GUF_100m_nations_20180315/UR_Tanzania.tif'

    #sFileGARPoint = '/Users/lauro/Desktop/EU_ACP_UNISDR_Africa/DATA/Swaziland/AFR_buiA_GAR_5km_20180313/africa_swz.shp'
    #sFileMask = '/Users/lauro/Desktop/EU_ACP_UNISDR_Africa/DATA/Swaziland/AFR_buiA_GUF_100m_nations_20180315/Swaziland.tif'

    # # from points 10km to raster 10km
    sFileGARPointBase= os.path.basename(sFileGARPoint)
    dir_out = os.path.join(os.path.dirname(sFileGARPoint),sFileGARPointBase.split('.')[0]+"_"+"outputs")
    if not os.path.exists(dir_out):
        os.mkdir(dir_out)

    sFileGARPointBuildingType =  os.path.join(dir_out, sFileGARPointBase.split('.')[0]+'_'+sBuildingType+".shp")
    #os.system('ogr2ogr -overwrite -f "ESRI Shapefile" -where "curve=\''+sBuildingType+'\'" ' + sFileGARPointBuildingType + ' ' + sFileGARPoint)

    inDriver = ogr.GetDriverByName("ESRI Shapefile")
    inDataSource = inDriver.Open(os.path.join(sFileGARPoint), 0)
    inLayer = inDataSource.GetLayer()
    inLayerDefn = inLayer.GetLayerDefn()
    log_print("Reading GAR shape file...")
    log_print("Total features: "+str(inLayer.GetFeatureCount()))

    #inLayer.SetAttributeFilter("curve = \''+sBuildingType+'\'")

    (fXmin, fXmax, fYmin, fYmax) = inLayer.GetExtent()
    iXLowRes = 0.04167
    iYLowRes = 0.04167
    fXmin = fXmin - iXLowRes / 2
    fYmin = fYmin - iYLowRes / 2
    fXmax = fXmax + iXLowRes / 2
    fYmax = fYmax + iYLowRes / 2

    outDriver = ogr.GetDriverByName( "ESRI Shapefile")
    outDataSource = outDriver.CreateDataSource(sFileGARPointBuildingType)
    proj = inLayer.GetSpatialRef()
    ##Creating the layer with its fields
    outLayer = outDataSource.CreateLayer(sFileGARPointBuildingType, proj, ogr.wkbPoint )
    field_exposure = ogr.FieldDefn(sExposureType, ogr.OFTReal)
    outLayer.CreateField(field_exposure)

    pop = []
    geoms = []
    ii=0
    for feature in inLayer:
        local_curve = feature.GetField("curve")

        if local_curve is not None and '.' in local_curve:

#            ii=ii+1;print (str(ii) +' - '+ local_curve.split('.')[0])

            if local_curve.split('.')[0]==sBuildingType:

                if feature.GetGeometryRef().GetPoint() not in geoms:
                    geoms.append(feature.GetGeometryRef().GetPoint())
                    pop.append(feature.GetField(sExposureType))

                else:
                    id = geoms.index(feature.GetGeometryRef().GetPoint())
                    pop[id]=pop[id]+feature.GetField(sExposureType)

    for i in range(len(geoms)):
        outFeature= ogr.Feature(outLayer.GetLayerDefn())
        outFeature.SetField(sExposureType, pop[i])
        point = ogr.Geometry(ogr.wkbPoint)
        point.AddPoint(*geoms[i])
        outFeature.SetGeometry(point)
        outLayer.CreateFeature(outFeature)

    total_exposure = sum (pop)
    log_print("Total unique features in outLayer: "+str(outLayer.GetFeatureCount()))
    log_print("Total original exposure value (for curve type="+sCurveType+", exposure type="+sExposureType+"): %.2f" % total_exposure)

    inLayer.ResetReading()
    outDriver = None
    outFeature = None
    outLayer = None
    inDataSource = None
    outDataSource = None

#    for i in range (0, len(geoms)):

    log_print("Converting points to raster...")

    sFileGARRasterEXPLowRes = os.path.join(dir_out, sFileGARPointBase.split('.')[0] + "_"+sBuildingType+"_"+sExposureType+".tif")
    if bVerbose: print('gdal_rasterize -co compress=DEFLATE -a '+sExposureType+' -l "' + os.path.basename(sFileGARPointBuildingType).split('.')[0] + '" -tr ' + str(iXLowRes) + ' ' + str(iYLowRes) + ' -te ' + str(fXmin) + ' ' + str(fYmin) + ' ' + str(fXmax) + ' ' + str(fYmax) + ' "' + sFileGARPointBuildingType + '" "'+sFileGARRasterEXPLowRes+'"')
    os.system('gdal_rasterize -co compress=DEFLATE -a '+sExposureType+' -l "' + os.path.basename(sFileGARPointBuildingType).split('.')[0] + '" -tr ' + str(iXLowRes) + ' ' + str(iYLowRes) + ' -te ' + str(fXmin) + ' ' + str(fYmin) + ' ' + str(fXmax) + ' ' + str(fYmax) + ' "' + sFileGARPointBuildingType + '" "'+sFileGARRasterEXPLowRes+'"')

    # if bVerbose: print('gdal_rasterize -co compress=DEFLATE -a '+sExposureType+' -where "curve=\''+sBuildingType+'\'" -l ' + os.path.basename(sFileGARPoint).split('.')[0] + ' -tr ' + str(iXLowRes) + ' ' + str(iYLowRes) + ' -te ' + str(fXmin) + ' ' + str(fYmin) + ' ' + str(fXmax) + ' ' + str(fYmax) + ' ' + sFileGARPoint + ' ' +sFileGARRasterEXPLowRes)
    # os.system('gdal_rasterize -co compress=DEFLATE -a '+sExposureType+' -where "curve=\''+sBuildingType+'\'" -l ' + os.path.basename(sFileGARPoint).split('.')[0] + ' -tr ' + str(iXLowRes) + ' ' + str(iYLowRes) + ' -te ' + str(fXmin) + ' ' + str(fYmin) + ' ' + str(fXmax) + ' ' + str(fYmax) + ' ' + sFileGARPoint + ' ' +sFileGARRasterEXPLowRes)

    [xsize_orig, ysize_orig, geotransform, geoproj, data_low] = readFile(sFileGARRasterEXPLowRes)
    log_print ("Original GAR map size: ("+str(ysize_orig)+","+str(xsize_orig)+") (rows,cols)")

    uniq = non_unique (data_low.ravel().tolist())
    data_low=np.reshape(uniq, data_low.shape)

    writeGeotiffSingleBand(sFileGARRasterEXPLowRes, 
                    geotransform, geoproj, 
                    data_low, nan_value=NAN_VALUE)

    if bDisplayImages: plot_image(data_low)

    del data_low, uniq

    sFileGARRasterEXPHiRes  = os.path.join(dir_out, sFileGARPointBase.split('.')[0] + "_"+sBuildingType+"_"+sExposureType+"_regrid.tif")

    match_geotrans, match_proj = rasterRegrid(sFileGARRasterEXPLowRes, sFileMask, sFileGARRasterEXPHiRes ,"nearest")

    [xsize, ysize, geotransform_high, geoproj_high, data_high]   = readFile(sFileGARRasterEXPHiRes)
    #[xsize, ysize, geotransform_high, geoproj_high, data_mask]   = readFile(sFileMask, fix_nan_value=255)
    [xsize, ysize, geotransform_high, geoproj_high, data_mask]   = readFile_withNoData(sFileMask)

    data_mask_zeroes = data_mask.copy()
    data_mask[data_mask==0] = np.nan
    data_high_masked = data_high * data_mask
    
    if bDisplayImages:
        plot_image(data_mask)
        plot_image(data_high)
        plot_image(data_high_masked)

    log_print ("Starting counting...")
    # take out the nan values from the array
    values = data_high_masked.ravel()
    values = values[~np.isnan(values)]

    dictCounter = Counter()
    # show a progressbar if the array is "big"
    # count 1000 elements at a time
    if len(values)>1000000:
        stride = 1000
        bar = PB()
        for v in bar(range(0, len(values), stride)):
            dictCounter.update(values[v:v+stride])
    else:
        dictCounter.update(values)

    
    log_print ("Data counted")

    data_mask_gar = np.copy(data_high)
    ratio = abs ((geotransform[1]*geotransform[5]) / (geotransform_high[1]*geotransform_high[5])) #number of high res pixels in low res pixel

    data_gar_excluded = []
    for key in data_high.ravel().tolist():    # for name, age in list.items():  (for Python 3.x)
        if key not in dictCounter.keys() and key != 0 :
            dictCounter[key] = ratio
            data_mask_gar[data_mask_gar==key]=-9999
            data_gar_excluded.append(key)
    excluded_exposure = sum(data_gar_excluded)
    log_print ("GAR exposure value not overlapping mask: %.2f (%.1f %%)" % (excluded_exposure , excluded_exposure/total_exposure*100))

    #building mask (union of the original mask + non-zero values form GAR)
    data_mask_gar[data_mask_gar!=-9999]=0
    data_mask_gar[data_mask_gar==-9999]=1
    data_mask = data_mask_gar + data_mask
    data_mask[data_mask>0] = 1
    data_high_masked = data_high * data_mask

    del data_high, data_mask_gar

    if bVerbose: log_print ("Counter length: "+str(dictCounter.__len__()))
    if bVerbose: log_print ("Unique length: "+str(len(np.unique(data_high_masked))))

    for key,value  in dictCounter.items():
        data_high_masked [data_high_masked==key] = key/value
        if bVerbose: print ('Amount of pixel for population ' + str(key) + ': ' + str(value))

    if bDisplayImages:
        fig = plt.figure()
        cax = plt.imshow(data_high_masked[:, :])
        plt.colormaps()
        cbar = fig.colorbar(cax, orientation='vertical')  # vertically oriented colorbar
        #cbar.ax.set_yticklabels(['< -1', '0', 1, 2,'> 10'])
        plt.show()

    sFileDownscaled = os.path.join(dir_out, sFileGARPointBase.split('.')[0] + "_"+sBuildingType+"_"+sExposureType+"_HighRes.tif")
    bandMetadata = {}
    bandMetadata[("unit", "1")] = "USD"

    if sExposureType == "VALHUM":
        description = 'VALHUM in ' + sBuildingType + ' (downscaling of GAR2015 data)'
    elif sExposureType == "VALFIS":
        description = 'VALFIS in ' + sBuildingType + ' (downscaling of GAR2015 data)'

    iNbands = 1
    data_high_masked_USD = data_high_masked*1000000
    
    writeFile(sFileDownscaled, match_geotrans, 
                match_proj, data_high_masked_USD, 
                bandMetadata, description, 
                {'Building_type': sBuildingType}, iNbands,
                nan_value=NAN_VALUE)

    log_print("Total downscaled exposure value (for curve type="+sCurveType+", exposure type="+sExposureType+"): %.2f" % np.sum(data_high_masked))
    log_print ("Exposure downscaled. Final result save in file: " + sFileDownscaled)

    sMaskFileDownscaled = os.path.join(dir_out, sFileGARPointBase.split('.')[0] + "_mask_HighRes.tif")    
    
    writeFile(sMaskFileDownscaled, match_geotrans, 
                match_proj, data_mask, 
                bandMetadata, description, 
                {}, iNbands,
                nan_value=NAN_VALUE)

    log_print ("Mask saved in file: " + sMaskFileDownscaled)
    
if __name__ == "__main__":

    bVerbose = False
    print(" ")
    log_print ('Process started...\n')
    argv=sys.argv[1:]

    try:
        a1Opts, a1Args = getopt.getopt(
                           argv,
                           "vhi:m:e:f:n:",
                           ["verbose","help","input=","mask=","exposure=","filter=","nan="]
        )
    except getopt.GetoptError:
        displayHelp();
        sys.exit(2)
    if a1Opts==[]:
        print ('\nPlease provide an input file...')
        displayHelp();
        sys.exit(2);
        
    for opt, arg in a1Opts:
        print(opt, arg)
        if opt in ('-h',"--help"):
            displayHelp();
            sys.exit()
        elif opt in ("-i", "--input"):
            SINGLE_POP_VALUE = None
            try:
                SINGLE_POP_VALUE = float(arg)
            except:
                if os.path.isfile(arg): 
                    sFileGARPoint = arg
                else:
                    print(f'error, {arg} is neither a value or a file')
                    sys.exit(2)                    
                    
            #sFileGARPoint = arg
        elif opt in ("-m", "--mask"):
            sFileMask = arg
        elif opt in ("-e", "--exposure"):
            sExposureType = arg
        elif opt in ("-f", "--filter"):
            sCurveType = arg
        elif opt in ("-n", "--nan"):
            NAN_VALUE = float(arg)
        elif opt in ("-v", "--verbose"):
            bVerbose=True
        else:
            print ('\nPlease provide an input files...')
            displayHelp();
            sys.exit(2)


    if sExposureType == "VALFIS":
        if sFileGARPoint is None:
            print(f'error, input file required')
            sys.exit(2)
            
        driver = ogr.GetDriverByName('ESRI Shapefile')
        sDataSource = driver.Open(sFileGARPoint, 0) # 0 means read-only. 1 means writeable.
        oLayer = sDataSource.GetLayer()

        if sCurveType == "ALL":
            aBuildingTypes = []
            aBuildingTypesUnique = []
            for i in range(0, oLayer.GetFeatureCount()):
                oFeature = oLayer.GetFeature(i)
                curve = oFeature.GetField("curve")
                if not curve or '.' not in curve:
                    print('Warning: curve for %s is not defined' % oFeature.GetField("ID_5X5"))
                    continue
                type = curve.split('.')[0]
                if type not in aBuildingTypesUnique:
                    aBuildingTypesUnique.append(type)
            log_print("Exposure type: "+sExposureType)
            log_print("Total curve types:"+str(aBuildingTypesUnique))
            for sCurveType in aBuildingTypesUnique:
                downscale (sFileGARPoint, sFileMask, sCurveType, sExposureType)
        else:
            downscale (sFileGARPoint, sFileMask, sCurveType, sExposureType)
            
    elif sExposureType == "VALHUM":
        if SINGLE_POP_VALUE is not None:
            single_value_to_population(SINGLE_POP_VALUE, sFileMask)
        else:
            downscale_pop(sFileGARPoint, sFileMask, sExposureType)   
    else:
        print ("Exposure type incorrect!")
        print ("Please use one of the following values: VALHUM or VALFIS")

    #sFileGARPoint = '/Users/lauro/Desktop/EU_ACP_UNISDR_Africa/DATA/UR_Tanzania/AFR_buiA_GAR_5km_20180313/africa_tza.shp'
    #sFileMask = '/Users/lauro/Desktop/EU_ACP_UNISDR_Africa/DATA/UR_Tanzania/AFR_buiA_GUF_100m_nations_20180315/UR_Tanzania.tif'
    #downscale (sFileGARPoint,sFileMask,"A")

    #sFileGARPoint = '/Users/lauro/Desktop/EU_ACP_UNISDR_Africa/DATA/Swaziland/AFR_buiA_GAR_5km_20180313/africa_swz.shp'
    #sFileMask = '/Users/lauro/Desktop/EU_ACP_UNISDR_Africa/DATA/Swaziland/AFR_buiA_GUF_100m_nations_20180315/Swaziland.tif'
    #downscale (sFileGARPoint,sFileMask,"UCB")

    #clip (sFileRaster, sFileClipper,sAttributeName,'')
    log_print ("Process ended.\n")