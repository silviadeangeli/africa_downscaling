# coding: utf-8
__author__ = 'mirko'

import geopandas as gpd
import pandas as pd
import numpy as np
import rasterio as rio
import argparse
import json
import os

from africa_tools import log_print, rasterize, geotiff_writer, extract_geo_info, reproject, ID

def init_parser():       
    parser = argparse.ArgumentParser(description='''
        Extract a new mask by clipping the original mask where the GAR points are defined''',
        epilog='''
        Usage example:       
        python clip_on_gar.py shapefile.shp mask_file.tiff output_file.tiff''')
    
    parser.add_argument('shapefile', help='shapefile to process',
        type=argparse.FileType('r'))
    
    parser.add_argument('mask',
                        help='geotiff file for regrid',
                        type=argparse.FileType('r'))
    
    parser.add_argument('output_file', 
                        help='output file',
                        type=str)
    return parser
    

def process(shapefile, mask, output_file):
    df = gpd.read_file(shapefile)
    log_print('getting geometries')    

    geometries = df.groupby(ID).agg({'geometry':'first'})

    out_crs = rio.crs.CRS({'init': 'EPSG:4326', 'no_defs': True})            
    trans, shape = extract_geo_info(geometries)

    with rio.open(mask) as ref_raster:
        out_trans, out_shape = ref_raster.affine, ref_raster.shape
        mask_data = ref_raster.read(1)


    gdf_geom = gpd.GeoDataFrame(geometries)
    gdf_geom['value'] = 1

    log_print(f'writing {output_file}')
    with geotiff_writer(output_file, out_trans, out_crs, out_shape, 1) as writer:
        gdf_geom = gpd.GeoDataFrame(geometries)
        # create the low res raster for the band
        raster = rasterize(gdf_geom, 'value', shape, trans, out_crs) 

        # upscale the raster to the mask resolution
        hr_raster = reproject(raster, trans, out_trans, out_shape, out_crs)
        new_mask = (mask_data * hr_raster)
        # write using the context-manager writer
        writer.write(new_mask.astype('uint8'), indexes=1)

    perc_missing = 100.0 * np.nansum(mask_data - new_mask) / np.nansum(mask_data)        
    log_print('Percentage missing: %.2f'%perc_missing)

            

if __name__ == '__main__':   
    parser = init_parser()
    args = parser.parse_args()
    log_print(f'processing {args.shapefile.name}')
    
    args.shapefile.close()
    args.mask.close()
    process(args.shapefile.name, args.mask.name, args.output_file)
    log_print(f'finished processing {args.shapefile.name}')