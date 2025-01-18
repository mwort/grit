#!/usr/bin/env python
import os
import os.path as osp
from collections import defaultdict
import warnings

import geopandas as gpd
from geopandas import read_file, GeoDataFrame
import pandas as pd
import fiona
from tqdm import tqdm

def parse_geometry(geometry_string):
    from shapely import wkt
    from shapely.geometry import MultiPoint, Point
    try:
        geometry = wkt.loads(geometry_string)
        # import ipdb;ipdb.set_trace()
        if isinstance(geometry, MultiPoint):
            # Convert MultiPoint to Point
            geometry = Point(geometry.centroid)
        elif isinstance(geometry, Point):
            # If it's already a Point, leave it as is
            pass
        else:
            raise ValueError("Invalid geometry format")
        return geometry
    except Exception as e:
        raise ValueError(f"Error parsing geometry: {str(e)}")

def points_to_geo(input, x_column="lon", y_column="lat", epsg=4326, output=None):
    """Convert a dataframe with point entries to a geopandas dataframe."""
    df = pd.read_csv(input) if type(input) == str else input
    try:
        geometry=df['geometry'].apply(parse_geometry)
    except:
        geometry = gpd.points_from_xy(x=df[x_column], y=df[y_column])
    gdf = GeoDataFrame(df, geometry=geometry, crs=epsg)
    if output:
        gdf.to_file(output)
        return
    return gdf


def geo_to_csv(input, output=None, keep_coor=False, keep_lonlat=False, **readkw):
    """Convert a geopandas readable dataframe/file to a dataframe."""
    gdf = read_file(input, **readkw) if type(input) == str else input
    if 'MultiPoint' in gdf.geometry.type.unique():
        gdf = gdf.explode(index_parts=False)
    df = pd.DataFrame(gdf.drop('geometry', axis=1))
    if keep_coor:
        df["x"] = gdf.geometry.x
        df["y"] = gdf.geometry.y
    if keep_lonlat:
        lonlat = gdf.to_crs(epsg=4326)
        df["lon"] = lonlat.geometry.x
        df["lat"] = lonlat.geometry.y 
    if output:
        df.to_csv(output)
        return
    return df


def merge_gpkg(output, *input, index=None, preserve_index=False,
               src_layer_field_name=None, crs=None):
    """Merge multiple geopandas readable files into a single gpkg file.

    Arguments:
    ----------
    output : gpkg output path (will be overwritten if it exists)
    *input : input file paths
    index : name of existing unique int column to use for the fid index
    preserve_index : if index is given, keep the index if true
    src_layer_field_name : add column with input file basename
    """
    layers_out = defaultdict(lambda: [])
    for f in tqdm(input, "Reading files"):
        # preserve layers in output
        layers = fiona.listlayers(f)
        src = osp.splitext(osp.basename(f))[0]
        for l in layers:
            inp = read_file(f, layer=l)
            if type(src_layer_field_name) is str:
                inp[src_layer_field_name] = src
            layers_out[l].append(inp)
    # make sure output doesnt exist, otherwise layers get appended
    if osp.exists(output):
        os.remove(output)
    # write out
    for i, l in tqdm(enumerate(layers_out), "Writing layers", total=len(layers_out)):
        df = pd.concat(layers_out[l])
        if type(index) is str and (index not in df.columns or df[index].isna().sum()):
            warnings.warn(f"Index column {index} does not exist or includes NAs. Skipping.")
        elif type(index) is str:
            df["fid"] = df[index] if preserve_index else df.pop(index)
        if crs:
            print(f"Reprojecting to {crs}...")
            df = df.to_crs(crs)
        df.to_file(output, layer=l)
    return


if __name__ == "__main__":
    from commandline import interface as _cli
    _cli()
