# coding: utf-8
__author__ = 'mirko'

import geopandas as gpd
import pandas as pd
import numpy as np
import rasterio as rio

import argparse
import json
import os

from contextlib import contextmanager
from rasterio import crs, transform, warp, enums, features
from functools import partial
from datetime import datetime
from progressbar import ProgressBar as PB

CURVE = 'curve'
ID = 'ID_5X5'
USE_SECTOR = 'USE_SECTOR'
VALUE_COLUMN = 'VALFIS'

use_sector_table = {
        "emp_serv":  "SERV",
        "emp_agr":  "AGR",
        "emp_gov":  "GOV",
        "emp_ind":  "IND",
        "ic_low":    "RES_LI", 
        "ic_high":    "RES_MHI", 
        "ic_mhg":    "RES_MHI", 
        "ic_mlow":    "RES_MHI"
}
band_order = ["RES_LI", "RES_MHI", "SERV", "AGR", "GOV", "IND"]

def log_print(s):
    print(f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} - {s}")

def override_config(config_file):
    config_obj = json.load(config_file)
    use_sector_table = config_obj['use_sector_table']
    band_order = config_obj['band_order']

def init_parser():       
    parser = argparse.ArgumentParser(description='''
        Process a shapefile and extracts the % of use by curve.
        Writes a tiff file for each curve''',
        epilog='''
        Usage example:       
        python extract_use.py shapefile.shp config.json mask_file.tiff output_dir''')
    
    parser.add_argument('shapefile', help='shapefile to process',
        type=argparse.FileType('r'))
    
    parser.add_argument('config',
                        help='json file with the "use" configuration',
                        type=argparse.FileType('r'))

    parser.add_argument('mask',
                        help='geotiff file for regrid',
                        type=argparse.FileType('r'))
    
    parser.add_argument('outdir', 
                        help='output directory',
                        type=str)
    return parser
    
def rasterize(gdf_geom, column, shape, trans, dst_crs, dtype=np.uint8, nodata=0):
    """rasterize the column value of the GeoDataFrame to the 
        raster defined by shape and trans
    """
    raster_data = features.rasterize(
        shapes=[
            (r[1].geometry, r[1][column]) 
            for r in gdf_geom.iterrows()
        ],
        out_shape=shape, 
        fill=nodata, 
        transform=trans,
    )
    return raster_data.astype(dtype)

@contextmanager
def geotiff_writer(filename, trans, dst_crs, shape, n_bands, 
                   dtype=np.uint8, nodata=0):
    """writes the raster as a multiband geotiff
    returns an open geotiff writer
    """
    with rio.Env():
        with rio.open(
                filename,
                'w',
                driver='GTiff',
                width=shape[1],
                height=shape[0],
                count=n_bands,
                dtype=dtype, 
                nodata=nodata,
                transform=trans,
                crs=dst_crs) as f:
            yield f


def get_use_on_df(df, idx):
    r = df.iloc[idx]
    if r[USE_SECTOR] in use_sector_table:
        return use_sector_table[r[USE_SECTOR]] 
    else:
        return r[USE_SECTOR]

def extract_geo_info(geometries):
    """extract geographical informations for the raster
    returns transformation and shape
    """
    lon, lat = zip(*list(map(lambda p: (p.x, p.y), geometries.geometry)))
    west, east, south, north = min(lon), max(lon), min(lat), max(lat)

    diff_lat = np.abs(np.diff(lat))
    diff_lon = np.abs(np.diff(lon))
    step_x = np.min(diff_lon[diff_lon>0])
    step_y = np.min(diff_lat[diff_lat>0])
    cols = int((east - west)/step_x) + 1
    rows = int((north - south)/step_y) + 1

    trans = rio.transform.from_bounds(
                            west  -step_x/2, 
                            south -step_y/2, 
                            east  +step_x/2, 
                            north +step_y/2, 
                            cols, rows)
    
    return trans, (rows, cols)

def reproject(src, src_trans, dst_trans, dst_shape, _crs, 
              src_nodata=0, dst_nodata=0, dtype=None):
    
    if dtype==None:
        dtype = src.dtype
        
    dst = np.empty(dst_shape, dtype=dtype) 
    with rio.Env():
        warp.reproject(
            source=np.ascontiguousarray(src), 
            destination=dst,
            src_crs=_crs, 
            dst_crs=_crs,
            dst_transform=dst_trans, 
            src_transform=src_trans,
            src_nodata=src_nodata,
            dst_nodata=src_nodata,
            resampling=enums.Resampling.nearest,
            num_threads=4
        )

    return dst


def check_dataframe(df):
    all_columns = True
    for C in (VALUE_COLUMN, ID, USE_SECTOR, CURVE):
        if C not in df:
            print(f'ERROR: column {C} not in shapefile')
            all_columns = False
    if not all_columns:
        sys.exit(-1)
        
    unique_use = df[USE_SECTOR].unique()
    missing_use = filter(
        lambda u: u not in use_sector_table, 
        unique_use
    )
    if len(list(missing_use))>0:
        log_print(f'warning missing uses: [{", ".join(missing_use)}]')

def process(shapefile, output_dir, mask):
    """processes the shapefile and writes the geotiff files
    """
    log_print('loading file')
    df = gpd.read_file(shapefile)

    check_dataframe(df)

    log_print('getting geometries')    

    geometries = df                                         \
                    .groupby(ID)                            \
                    .agg({'geometry':'first'})
    
    
    
    out_crs = rio.crs.CRS({'init': 'EPSG:4326', 'no_defs': True})            
    trans, shape = extract_geo_info(geometries)

    with rio.open(mask) as ref_raster:
        out_trans, out_shape = ref_raster.affine, ref_raster.shape
    
    os.makedirs(output_dir, exist_ok=True)
    log_print('writing geotiffs')
    
    filename = f'{output_dir}/{curve}.tiff'
    log_print(f'writing {filename}')
    with geotiff_writer(filename, out_trans, out_crs, out_shape, n_bands) as writer:
        # intersect geometries with data 
        gdf_geom = gpd.GeoDataFrame(geometries)

        # create the low res raster for the band
        raster = rasterize(gdf_geom, VALUE_COLUMN, shape, trans, out_crs) 
        # upscale the raster to the mask resolution
        out_raster = reproject(raster, trans, out_trans, out_shape, out_crs)
        # write using the context-manager writer
        writer.write(out_raster, indexes=band+1)

            

if __name__ == '__main__':   
    parser = init_parser()
    args = parser.parse_args()
    log_print(f'processing {args.shapefile.name}')
    
    if args.config:
        override_config(args.config)
        
    if args.outdir is None:
        args.outdir = './'
    
    args.shapefile.close()
    args.mask.close()
    process(args.shapefile.name, args.outdir, args.mask.name)
    log_print(f'finished processing {args.shapefile.name}')