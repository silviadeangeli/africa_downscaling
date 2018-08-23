from matplotlib.pyplot import plot
import gdal
import numpy as np
# filename= ('/Users/lauro/Documents/PROJECTS/BOLIVIA/Bolivia2/simulazioni/Prova6/Results6b_verbose2_WATERDEPTH_36000.tif');

# readFile(filename):
def readFile(filename):
    filehandle = gdal.Open(filename)
    band1 = filehandle.GetRasterBand(1)
    geotransform = filehandle.GetGeoTransform()
    geoproj = filehandle.GetProjection()
    band1data = band1.ReadAsArray()
    xsize = filehandle.RasterXSize
    ysize = filehandle.RasterYSize
    return xsize,ysize,geotransform,geoproj,band1data

def readFileBand(filename,band):
    filehandle = gdal.Open(filename)
    band1 = filehandle.GetRasterBand(band)
    geotransform = filehandle.GetGeoTransform()
    geoproj = filehandle.GetProjection()
    band1data = band1.ReadAsArray()
    xsize = filehandle.RasterXSize
    ysize = filehandle.RasterYSize
    return xsize,ysize,geotransform,geoproj,band1data



def readFile_withNoData(filename):
    print(filename)
    filehandle = gdal.Open(filename)
    print(filehandle)
    band1 = filehandle.GetRasterBand(1)
    print(band1)
    geotransform = filehandle.GetGeoTransform()
    geoproj = filehandle.GetProjection()
    band1data = band1.ReadAsArray()
    print(np.max(band1data))
    band1dataNP=np.array(band1data).astype(float)
    nodata = band1.GetNoDataValue()
    band1dataNP[band1dataNP==nodata]=np.nan
    print(band1data)
    print(np.max(band1data))
    xsize = filehandle.RasterXSize
    ysize = filehandle.RasterYSize
    return xsize,ysize,geotransform,geoproj,band1dataNP

def writeGeotiffSingleBand(filename, geotransform, 
                            geoprojection, data, 
                            nan_value=np.nan):
    if ~np.isnan(nan_value):
        data = data.copy()
        data[np.isnan(data)] = nan_value

    (x,y) = data.shape
    format = "GTiff"
    driver = gdal.GetDriverByName(format)
    #dst_datatype = gdal.GDT_Byte #byte
    dst_datatype = gdal.GDT_Float32
    dst_ds = driver.Create(filename,y,x,1,dst_datatype,options = [ 'COMPRESS=DEFLATE'])
    #sDATETIME= "2013:04:30 12:00:00"#The format is: "YYYY:MM:DD HH:MM:SS", with hours like those on a 24-hour clock, and one space character between the date and the time. The length of the string, including the terminating NUL, is 20 bytes.
    dst_ds.SetMetadata( {'TIFFTAG_SOFTWARE': 'Telemac2D', 'TIFFTAG_DATETIME': 'sDATETIME'} )
    dst_ds.GetRasterBand(1).WriteArray(data)
    dst_ds.SetGeoTransform(geotransform)
    dst_ds.SetProjection(geoprojection)
    return 1


def writeFile(filename, geotransform, 
                geoprojection, data, 
                bandMetadata, description, 
                name, iNbands, 
                nan_value=np.nan):
    #metedata is like:
    #metadata[('unit')] = 'tons'
    #metedata[('unit')] = 'USD/tons'
    #where the 2nd tag is the number of the band
    
    if ~np.isnan(nan_value):
        data = data.copy()
        data[np.isnan(data)] = nan_value

    if iNbands>1:
        (x,y,iNbands) = data.shape
    else:
        (x,y) = data.shape
    format = "GTiff"
    driver = gdal.GetDriverByName(format)
    dst_datatype = gdal.GDT_Float32
    dst_ds = driver.Create(filename,y,x,iNbands,dst_datatype,options = [ 'COMPRESS=DEFLATE'])
    #VALUES FOR THE WHOLE FILE
    dst_ds.SetDescription(description)
    dst_ds.SetMetadata( name )
    dst_ds.SetGeoTransform(geotransform)
    dst_ds.SetProjection(geoprojection)

    #write data
    if iNbands==1:
        dst_ds.GetRasterBand(1).WriteArray(data)
    else:
        for i in range (0,iNbands):
            dst_ds.GetRasterBand(i+1).WriteArray(data[:,:,i])
    #write metadata
#    for i in range (0,len(bandMetadata.keys())):
    for key, value in bandMetadata.items():
        #print (attributes, values)
        iBand = key[1]
        attribute = key[0]
        dst_ds.GetRasterBand(int(iBand)).SetMetadata({attribute: value} )
    return 1