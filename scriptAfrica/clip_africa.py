__author__ = 'lauro'
from osgeo import ogr
import sys, getopt, os
import datetime
from Geotiff import readFile, writeFile
import matplotlib.pylab as plt
import numpy as np

def simple_plot (a):
    fig,ax=plt.subplots()
    im = plt.imshow(a,interpolation='none')
    plt.colormaps()
    #cbar = fig.colorbar(cax,ticks=[-1, 0, 1, 2, 10],orientation='vertical')
    cbar = fig.colorbar(im,orientation='vertical')
    #plt.clim(0,1)
    #cbar.ax.set_yticklabels(['< -1', '0', 1, 2,'> 10'])# vertically oriented colorbar
    #ax.format_coord = Formatter(im)
    plt.show()


# gdalwarp -q -cutline "/Users/lauro/Documents/PROJECTS/FP7_H2020/RASOR_2012/MALAWI/Administrative limits/AdminLimitsUTM36S.shp" -crop_to_cutline -of GTiff /Users/lauro/Documents/PROJECTS/FP7_H2020/RASOR_2012/MALAWI/crop/rice/malawi_agri_rice_production30m.tif /Users/lauro/Documents/PROJECTS/FP7_H2020/RASOR_2012/MALAWI/crop/rice/malawi_agri_rice_production30mCropped.tif
def displayHelp():
    print ('\nScript to clip a rater in multiple tiles according to a feature attribute')
    print ('Options:')
    print ('         -v | --verbose                          print the output on the screen')
    print ('         -h | --help                             display this help')
    print ('         -c | --clipper                          file name of the shp file (clipper)')
    print ('         -i | --input                                input GeoTIFF file')
    print ('         -o | --output                           output core name')
    print ('         -f | --filter                           attribute field name (filter)')
    print ('         -g | --group                            limit to group of nations called group')
    print ('         -l | --level                            administrative level [0 or 1]')

    print ('\nExample: python clip_africa.py -i FILE.tif -c FILE_CLIPPER.shp -f ATTRIBUTE_FIELD_NAME')


def clip_africa(sFileIn, sFileClipper, sAttributeName, sFileOutCore,sGroup,sLevel):
    sFileName = os.path.splitext(os.path.basename(sFileIn))[0]
    sDirectory = os.path.dirname(sFileIn)
    driver = ogr.GetDriverByName('ESRI Shapefile')
    sDataSource = driver.Open(sFileClipper, 0)  # 0 means read-only. 1 means writeable.
    oLayer = sDataSource.GetLayer()

    for i in range(0, oLayer.GetFeatureCount()):
        #for i in range(0, 1):
        #for oFeature in oLayer:
        # Get the input Feature
        oFeature = oLayer.GetFeature(i)
        tile = oFeature.GetField(sAttributeName)
        sPrefix=sFileOutCore

        if sLevel == "0":
                if tile == sGroup or sGroup == "ALL":
                    if sFileOutCore == '':
                         sFileOut = sDirectory + '/' + tile + '.tif'
                    else:
                        sFileOut = sDirectory + '/' + sFileOutCore + '_' + tile + '.tif'
                    print (tile)
                    warp = '''gdalwarp -co "COMPRESS=DEFLATE" -crop_to_cutline -overwrite -cutline {clipper} -cwhere "{TILE_NAME} = '{tile}'" "{infile}" "{outfile}"'''.format(
                    clipper=sFileClipper, TILE_NAME=sAttributeName, tile=tile, infile=sFileIn, outfile=sFileOut)
                    print (warp)
                    os.system(warp)
                    listFilesOut.append(sFileOut)


        if sLevel == "1":
                sNation = "ADM0_NAME"
                if oFeature.GetField(sNation) == sGroup or sGroup == "ALL":
                    sPrefix = oFeature.GetField(sNation)
                    sFileOut = sDirectory + '/' + sFileOutCore + '_' + sPrefix + '_' + tile.lower() + '.tif'
                    print (sPrefix + '_' + tile.lower())
                    warp = '''gdalwarp -co "COMPRESS=DEFLATE" -crop_to_cutline -overwrite -cutline {clipper} -cwhere "{TILE_NAME} = '{tile}'" "{infile}" "{outfile}"'''.format(
                    clipper=sFileClipper, TILE_NAME=sAttributeName, tile=tile, infile=sFileIn, outfile=sFileOut)
                    print (warp)
                    os.system(warp)
                    listFilesOut.append(sFileOut)

        #else:
        #   sFileOut = sDirectory + '/' + sFileOutCore + '_' + tile.lower() + '.tif'


    return (listFilesOut)



if __name__ == "__main__":
    print ('\n' + str(datetime.datetime.now()) + ' - Process started...\n')
    argv = sys.argv[1:]

    sFileOutCore = ''

    try:
        a1Opts, a1Args = getopt.getopt(argv, "vhi:c:o:f:g:l:", ["verbose", "help", "ifile=", "ofile=", "group=", "filter=", "level="])
    except getopt.GetoptError:
        displayHelp();
        sys.exit(2)
    if a1Opts == []:
        print ('\nPlease provide an input file...')
        displayHelp();
        sys.exit(2);
    for opt, arg in a1Opts:
        if opt in ('-h', "--help"):
            displayHelp();
            sys.exit()
        elif opt in ("-v", "--verbose"):
            bVerbose = True
        elif opt in ("-i", "--ifile"):
            sFileIn = arg
        elif opt in ("-o", "--ofile"):
            sFileOutCore = arg
        elif opt in ("-c", "--clipper"):
            sFileClipper = arg
        elif opt in ("-f", "--filter"):
            sAttributeName = arg
        elif opt in ("-g", "--group"):
            sGroup = arg
        elif opt in ("-l", "--level"):
            sLevel = arg
        else:
            print ('\nPlease provide an input files...')
            displayHelp();
            sys.exit(2)

    listFilesOut = []
    listFilesOut = clip_africa(sFileIn, sFileClipper, sAttributeName, sFileOutCore, sGroup, sLevel)
    #print (a2sFileOut)

    for f in listFilesOut:
        print("masking: "+f)
        [xsize, ysize, geotransform, geoproj, data] = readFile(f)
        data=np.array(data)
        data[data < 255] = 0        #data[data <  3] = 0 ### in case of GHSL
        data[data == 255] = 1       #data[data >= 3] = 1 ### in case of GHSL
        bandMetadata = {}
        bandMetadata[('built-up', '1')] = '1'
        description = 'Built up layer form GHSL'
        #writeFile(os.path.splitext(f)[0]+"_masked.tif", geotransform, geoproj, data, bandMetadata, description, "built-up", 1)
        writeFile(f, geotransform, geoproj, data, bandMetadata, description, "built-up", 1)

    print (str(datetime.datetime.now()) + ' - Process ended.\n')