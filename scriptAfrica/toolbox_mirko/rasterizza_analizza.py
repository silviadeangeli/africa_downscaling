# coding: utf-8

from __future__ import division
import geopandas as gpd
import pandas as pd
import rasterio as rio
import numpy as np
from rasterio import crs, transform, warp, enums


def write_geotiff(filename, src, bounds):
    west, north, east, south = bounds
    rows, cols = src.shape
    trans = transform.from_bounds(west, south, east, north, cols, rows)
    out_crs = crs.CRS({'init': 'EPSG:4326', 'no_defs': True})
    
    with rio.Env():
        with rio.open(
                filename,
                'w',
                driver='GTiff',
                width=src.shape[1],
                height=src.shape[0],
                count=1,
                dtype=np.float32,
                nodata=0,
                transform=trans,
                crs=out_crs) as f: f.write(src.astype('float32'), indexes=1)


def get_matrix_shape(df):
    minx, miny, maxx, maxy = df.total_bounds
    
    diffx = np.diff(df.geometry.x)
    diffx[diffx==0] = np.NaN
    step = np.nanmin(abs(diffx))

    nrows = int((maxy-miny)/step) +1
    ncols = int((maxx-minx)/step) +1
    
    bounds = (minx-.5*step, miny-.5*step, maxx+.5*step, maxy+.5*step)
    shape = (nrows, ncols)
    return step, bounds, shape

def rasterize(df, column, step, bounds, shape):
    minx, miny, maxx, maxy = bounds

    M = np.ones(shape) * np.nan
    
    next_perc = 0
    for i, (_, p) in enumerate(df.iterrows()):
        perc = 100*(i/df.size)
        if perc>next_perc:
            print(str(round(perc)))
            next_perc += 5

        x, y = p.geometry.x, p.geometry.y
        c = int((x-minx)/step)
        r = int((y-miny)/step)
        try:
            val = float(p[column])
            if np.isnan(M[r,c]):
                M[r,c] = 0
            M[r,c] += val
        except:
            pass
    
    return M


def stat_string(mean_val_over_mask, std_val_over_mask, 
                perc_mapped_cells, perc_mapped_value, perc_missing, value_missing):
    
    final_str = '\n'.join([
        'media valore economico per cella di costruito: %.2f'  %(mean_val_over_mask*1000000),
        'std valore economico per cella di costruito: %.2f'  %(std_val_over_mask*1000000),
        'percentuale celle del GAR mappate dalla maschera di costruito: %.2f'  %(perc_mapped_cells),
        'percentuale valore del GAR mappato dalla maschera di costruito: %.2f' %(perc_mapped_value),
        'percentuale celle 5x5 con costruito non aventi valore GAR: %0.2f' % (perc_missing),
        'percentuale celle maschera di costruito non aventi valore GAR: %0.2f' % (value_missing)    
    ])
    return final_str

def statistics (data_values, data_mask, data_mask_resolution):
    # get the values where we have at least one point
    vd = (data_mask > 0) & (data_values > 0)

    val_over_mask = (data_values[vd])/(data_mask[vd]*data_mask_resolution**2) 
    
    mean_val_over_mask = np.nanmean(val_over_mask)
    std_val_over_mask = np.nanstd(val_over_mask)
    
    # evaluate statistics over positive VALFIS data
    pd = (data_values > 0)
    n_cells = np.sum(pd)
    perc_mapped_cells = 100 * (data_mask[vd].size/n_cells)

    all_gar_value = np.sum(data_values[pd])
    mapped_value = np.sum(data_values[pd&vd])
    perc_mapped_value = mapped_value*100/all_gar_value
    
    # evaluate missing data points
    value_missing =  100*np.sum(data_mask[(data_mask > 0) & ~(data_values > 0)])/  np.sum(data_mask[(data_mask > 0)])
    perc_missing = 100 * np.sum((data_mask > 0) & ~(data_values > 0)) / np.sum((data_mask > 0))
    

    return mean_val_over_mask, std_val_over_mask, perc_mapped_cells, perc_mapped_value, perc_missing, value_missing


def analyse(shp_5km, filename, out_tiff, valfis_tiff, out_txt):
    #shp_5km_out = shp_5km + '_out.pickle'
    
    #if (os.exist(shp_5km_out)):
    #    pass
        #M_VALFIS, M_VALHUM, step, bounds, shape = load(shp_5...)
        #load shp_5k_out
    #else:
    df_5 = gpd.read_file(shp_5km)
    ID_FIELD = 'ID_5X5'
    df_5.VALFIS = pd.to_numeric(df_5.VALFIS, errors='coerce')
    df_5.VALHUM = pd.to_numeric(df_5.VALHUM, errors='coerce')
    df_sum = gpd.GeoDataFrame()
    groups = df_5.groupby(ID_FIELD)
    df_sum['VALFIS'] = groups.VALFIS.sum()
    df_sum['VALHUM'] = groups.VALHUM.sum()
    df_sum.geometry = groups.geometry.first()
    step, bounds, shape  = get_matrix_shape(df_sum)
    M_VALFIS = rasterize(df_sum, 'VALFIS', step, bounds, shape)
    M_VALHUM = rasterize(df_sum, 'VALHUM', step, bounds, shape)
        
        #pickle.dump((M_VALFIS, M_VALHUM, step, bounds, shape), open(shp_5km_out, 'wb'))
        #save shp_5km_out

    data_raster = rio.open(filename)
    data = data_raster.read(1)
    data = data.astype('float')
    data[data==255] = np.NaN
    data[data>1] = 1

    min_lon, min_lat, max_lon, max_lat = bounds

    M2 = np.zeros(shape)
    for r in range(shape[0]):
        for c in range(shape[1]):
            lat = min_lat + (r*step)        
            lon = min_lon + (c*step)            

            pr, pc = data_raster.index(lon, lat)
            nr, nc = data_raster.index(lon+step, lat+step)
            m = data[nr:pr, pc:nc]


            M2[r,c] = np.sum(m>0)

    write_geotiff(out_tiff, M2, bounds)
    write_geotiff(valfis_tiff, M_VALFIS, bounds)

    result = stat_string(*statistics(M_VALFIS, M2, 90))
    f = open(out_txt,'w')
    f.write(result)
    f.close()
    return result

if __name__ == '__main__':    
    shp_5km = sys.argv[1]
    filename = sys.argv[2]
    out_tiff = sys.argv[3]
    valfis_tiff = sys.argv[4]
    out_txt = sys.argv[5]
    print(analyse(shp_5km, filename, out_tiff, valfis_tiff, out_txt))