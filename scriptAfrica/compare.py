__author__ = 'lauro'

from Geotiff import readFile
import matplotlib.pylab as plt
import numpy as np

def plot_image(data):
    fig = plt.figure()
    cax = plt.imshow(data, cmap='Blues')
    plt.colormaps()
    #plt.clim(0,400)
    cbar = fig.colorbar(cax, orientation='vertical')
    #cbar.ax.set_yticklabels(['< -1', '0', 1, 2,'> 10'])# vertically oriented colorbar
    plt.show()

#sFileLowRes  = '/Users/lauro/Desktop/EU_ACP_UNISDR_Africa/DATA/UR_Tanzania/AFR_buiA_GAR_5km_20180313/africa_tza_A_VALHUM.tif'
#sFileHighRes = '/Users/lauro/Desktop/EU_ACP_UNISDR_Africa/DATA/UR_Tanzania/AFR_buiA_GAR_5km_20180313/africa_tza_A_VALHUM_HighRes.tif'

sFileLowRes  = '/Users/lauro/Desktop/EU_ACP_UNISDR_Africa/DATA/Swaziland/AFR_buiA_GAR_5km_20180313/swz_ucb_valhum.tif'
sFileHighRes = '/Users/lauro/Desktop/EU_ACP_UNISDR_Africa/DATA/Swaziland/AFR_buiA_GAR_5km_20180313/africa_swz_UCB_VALHUM_HighRes.tif'


[xsize_orig, ysize_orig, geotransform, geoproj, data_low]   = readFile(sFileLowRes)
[xsize, ysize, geotransform, geoproj, data_high]   = readFile(sFileHighRes)

data_low [data_low<0]=0
#for el in data_low.ravel().tolist():
#    if el < 0:
#        print (el)

sum_low = np.sum(data_low)
sum_high = np.sum(data_high)

print (sum_low,sum_high)

#plot_image(data_low)
#plot_image(data_high)
