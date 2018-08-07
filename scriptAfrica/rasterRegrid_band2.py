__author__ = 'lauro'

from osgeo import gdal, gdalconst

def rasterRegrid (sFileIn ,sFileMatch, sFileOut, method):
    ## resample sFileIn at the same resolution of sFileMatch

    if method == "nearest":
        interpMethod = gdalconst.GRA_NearestNeighbour
    elif method == "max":
        interpMethod = gdalconst.GRA_Max
    else:
        interpMethod = gdalconst.GRA_NearestNeighbour
    #print (interpMethod)

    src_all = gdal.Open(sFileIn, gdalconst.GA_ReadOnly)
    src_proj = src_all.GetProjection()
    src_geotrans = src_all.GetGeoTransform()

    # The source will be converted in order to match this:
    match_ds = gdal.Open(sFileMatch, gdalconst.GA_ReadOnly)
    match_proj = match_ds.GetProjection()
    match_geotrans = match_ds.GetGeoTransform()
    wide = match_ds.RasterXSize
    high = match_ds.RasterYSize
    
    print ("Regridded map size: ("+str(high)+","+str(wide)+")")

    # Output / destination
    dst = gdal.GetDriverByName('GTiff').Create(sFileOut, wide, high, 1, gdalconst.GDT_Float32, options = ['COMPRESS=DEFLATE', 'PREDICTOR=3'])
    dst.SetGeoTransform( match_geotrans )
    dst.SetProjection( match_proj)

    # gdal.ReprojectImage (src, dst, src_proj, src_proj, gdalconst.GRA_NearestNeighbour)
    gdal.ReprojectImage (src, dst, src_proj, match_proj, interpMethod)

    return match_geotrans,match_proj

"""
sFileIn = '/share/africa/Global_datasets/Builtup/WRL_buiA_GHSL_38m_cuts_20180305/Ritagli_subnazionali/UR_Tanzania_morogoro.tif'
sFileMatch = '/share/africa/Global_datasets/Builtup/AFR_buiA_GUF_100m_cuts_20180312/UR_Tanzania_morogoro.tif'
sFileOut = '/share/africa/Global_datasets/Builtup/WRL_buiA_GHSL_38m_cuts_20180305/Ritagli_subnazionali_100m/UR_Tanzania_morogoro_100m.tif'
match_geotrans,match_proj=rasterRegrid (sFileIn ,sFileMatch, sFileOut)
print ('finish')
"""