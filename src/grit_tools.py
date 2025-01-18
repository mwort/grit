#!/usr/bin/env python
import warnings

import numpy as np
import pandas as pd
from tqdm import tqdm
import os
def create_rectangle(min_x, max_x, min_y, max_y,crs=4326):
    from shapely.geometry import Polygon
    # Create a rectangular polygon using Shapely
    polygon = Polygon([(min_x, min_y), (max_x, min_y), (max_x, max_y), (min_x, max_y)])
    # Create a GeoDataFrame with the polygon
    gdf = gpd.GeoDataFrame(geometry=[polygon],crs=crs)
    return gdf

def gdf_subset(input_file,boundary,output,epsg=None,spatial='intersects'):
    import geopandas as gpd
    import fiona
    from pyproj import CRS
    # get the crs of the input_file
    # import ipdb;ipdb.set_trace()
    with fiona.open(input_file) as src:
        epsg_input = int(CRS.from_string(src.crs['init']).to_epsg())
    if os.path.isfile(boundary):
        gdf_boun = gpd.read_file(boundary)
        if 'reach_id' in gdf_boun.columns:
            gdf_boun = gdf_boun.drop(columns=['reach_id'])
        if int(gdf_boun.crs.to_epsg())!=epsg_input:
            print(gdf_boun.crs.to_epsg())
            gdf_boun=gdf_boun.to_crs(epsg=epsg_input)
    else:
        try: 
            min_x,max_x,min_y,max_y=boundary.split(' ')
            print(min_x,max_x,min_y,max_y)
            # create rectangle from boun
            gdf_boun = create_rectangle(min_x, max_x, min_y, max_y,crs=epsg_input)
        except:
            raise ValueError('Cannot create polygon!')
    gdf_subset = gpd.read_file(input_file,bbox=gdf_boun)

    gdf_subset = gpd.sjoin(gdf_subset, gdf_boun,predicate=spatial)
    # process column names
    # Drop columns copied from the right dataframe
    columns_to_drop = [col for col in gdf_subset.columns if col.endswith('_right')]
    gdf_subset = gdf_subset.drop(columns=columns_to_drop)
    # Rename columns from the left dataframe (remove the '_left' suffix)
    gdf_subset.columns = [col.replace('_left', '') for col in gdf_subset.columns]
    # columns = gdf.columns
    # polygon check

    is_polygon = gdf_subset.geometry.geom_type.unique() == 'MultiPolygon'
    # is_polygon = gdf_subset.geometry.geom_type.unique() == 'Polygon'
    if is_polygon.all():
        print('Clip the input polygon by the boundary!')
        print('The rows before clipping is ',len(gdf_subset))
        gdf_subset = gpd.clip(gdf_subset,gdf_boun)
        print('The rows after clipping is ',len(gdf_subset))
        gdf_subset = gdf_subset.to_crs(epsg=epsg)
        # import ipdb;ipdb.set_trace()    
        gdf_subset['area'] = gdf_subset['geometry'].apply(lambda geom: geom.area)
        area_threshold = 10 # in square meters
        gdf_subset = gdf_subset[gdf_subset['geometry'].apply(lambda geom: geom.area) > area_threshold]
        print('The rows after area threshold is ',len(gdf_subset))
        gdf_subset = gdf_subset.to_crs(epsg_input)
    gdf_subset = gdf_subset.rename(columns={'global_id':'reach_id'})
    print(gdf_subset.columns,len(gdf_subset))
    # if 'catchment_id' not in gdf_subset.columns:
    #     try:
    #         gdf_subset = gdf_subet.set_index('reach_id')
    #         gdf_boun = gdf_boun.set_index('reach_id')
    #         idx = 
    #         gdf_subset['catchment_id']=gdf_boun['catchment_id']
    #     except:

    if len(gdf_subset)>0:
        gdf_subset.to_file(output,driver='GPKG',overwrite=True)
    else:
        raise ValueError('Subset dataframe is empty!')
def domain_search(domain,boundary,output):
    import geopandas as gpd
    gdf = gpd.read_file(domain)
    gdf_boun = gpd.read_file(boundary)
    if gdf_boun.crs!=gdf.crs:
        gdf_boun = gdf_boun.to_crs(gdf.crs)
    gdf_intersect = gpd.sjoin(gdf,gdf_boun,how='inner',predicate="intersects")
    df = gdf_intersect[['cat','domain']]
    df.to_csv(output,index=False)

def segment_index():
    import networkx as nx
    import river_network_utils as ru

    start_segment_id = 1320
    d=ru.river_networkx("segments")
    source_st, source_en, _ = [(s, e, k) for s, e, k, d in d.edges(keys=True, data=True) if d["id"] == start_segment_id][0]
    downstream_segment_ids = [d.get_edge_data(*n, key=0)["id"] for n in nx.bfs_edges(d, source_en)]
    
if __name__ == "__main__":
    from commandline import interface as _cli
    _cli()
