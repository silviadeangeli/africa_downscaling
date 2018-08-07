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
from functools import partial
from rasterio import crs, transform, warp, enums, features
from africa_tools import log_print, rasterize, geotiff_writer, extract_geo_info, reproject, ID

CURVE = 'curve'
USE_SECTOR = 'USE_SECTOR'
VALUE_COLUMN = 'VALFIS'

use_sector_table = {
        "emp_serv":  "SERV",
        "emp_agr":  "AGR",
        "emp_gov":  "GOV",
        "emp_ind":  "IND",
        "ic_low":    "RES_LI", 
        "ic_high":    "RES_MHI", 
        "ic_mhigh":    "RES_MHI", 
        "ic_mlow":    "RES_MHI"
}
band_order = ["RES_LI", "RES_MHI", "SERV", "AGR", "GOV", "IND"]
curve_name_mapping = {
    "I_T1(m).fvu": "T1",
    "I_M1(m).fvu": "M1",
    "I_M2(m).fvu": "M2",
    "I_W1(m).fvu": "W1",
    "I_C1(m).fvu": "C1"
}

def override_config(config_file):
    config_obj = json.load(config_file)
    use_sector_table = config_obj['use_sector_table']
    band_order = config_obj['band_order']
    curve_name_mapping = config_obj['curve_name_mapping']

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

def get_use_on_df(df, idx):
    r = df.iloc[idx]
    if r[USE_SECTOR] in use_sector_table:
        return use_sector_table[r[USE_SECTOR]] 
    else:
        return r[USE_SECTOR]

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

def fix_percentage(x):
    """
        get the int percentage of the values in the dataframe, 
        fixes the highest value if the sum is not 100%
    """
    x_perc = round(100 * x/x.sum())
    if (x_perc.sum()[0])<100:
        err = 100-x_perc.sum()        
        x_perc.loc[x_perc.idxmax()] += err
    return x_perc
                  

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
    
    
    
    log_print('aggregate on use')            
    get_use = partial(get_use_on_df, df)
    df_abs_values = df                                          \
                    .loc[                                       \
                        df[USE_SECTOR].isin(use_sector_table)   \
                    ]                                           \
                    .groupby([CURVE, get_use, ID])              \
                    .agg({VALUE_COLUMN: 'sum'})

    log_print('calculating percentage')
    df_perc_values = df_abs_values                              \
                     .groupby(level=(CURVE,ID))                 \
                     .apply(fix_percentage)   

    out_crs = rio.crs.CRS({'init': 'EPSG:4326', 'no_defs': True})            
    trans, shape = extract_geo_info(geometries)

    with rio.open(mask) as ref_raster:
        out_trans, out_shape = ref_raster.affine, ref_raster.shape
    
    os.makedirs(output_dir, exist_ok=True)
    
    log_print('writing geotiffs')
    for curve, df_curve in df_perc_values.groupby(CURVE):
        curve_name = curve_name_mapping.get(curve, curve)
        filename = f'{output_dir}/{curve_name}.tiff'
        n_bands = len(band_order)        
        log_print(f'writing {filename}')
        with geotiff_writer(filename, 
                            out_trans, 
                            out_crs, 
                            out_shape, 
                            n_bands,
                            dtype=np.uint8
                           ) as writer:
            for use, df_data in df_curve.groupby(level=1):                
                #select the band as configured
                band = band_order.index(use)
                
                log_print(f'   use "{use}" as band {band+1}')
                
                # intersect geometries with data 
                df_geom = df_data.join(geometries)
                gdf_geom = gpd.GeoDataFrame(df_geom)
                
                # create the low res raster for the band
                raster = rasterize(gdf_geom, VALUE_COLUMN, shape, trans, out_crs, 
                                   dtype=np.float32, nodata=0.0) 
                # upscale the raster to the mask resolution
                out_raster = reproject(raster, trans, out_trans, out_shape, out_crs)
                # write using the context-manager writer
                out_raster_int = (np.round(out_raster)).astype(np.uint8)
                writer.write(out_raster_int, indexes=band+1)

            

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