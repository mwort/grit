#!/usr/bin/env python
import gc

from tqdm import tqdm
from grass_attribute_table import GrassAttributeTable
import grass.script as grass


# category ints to merge with GRWL, first one found
GRWL_CATEGORIES = [
    ("canal", 0),  # GRWL normally 86, but shouldnt be included in DEM conditioning
    ("pond", 0),
    ("lake", 180),
]
GRWL_DEFAULT = 255  # river

MEMORY = 10000


def categorise_grwl(vector, default=None):

    tbl = GrassAttributeTable(vector)
    tbl["grwl_value"] = default or GRWL_DEFAULT
    stats = {k: 0 for k, v in GRWL_CATEGORIES}

    for i, row in tqdm(tbl.iterrows(), total=len(tbl)):
        st = str(row.values).lower()
        for c, v in GRWL_CATEGORIES:
            found = c in st
            if found:
                tbl.loc[i, "grwl_value"] = v
                stats[c] += 1
                break
    print(f"Found these categories: {stats}")
    tbl.write()
    return


def make_grwl_raster(vector, output, default=None):
    osm_vect = vector.split("@")[0]
    grass.run_command("g.copy", vector=(vector, osm_vect))
    categorise_grwl(osm_vect, default=default)
    gc.collect()
    grass.run_command("v.to.rast", input=osm_vect, output=output, use="attr", memory=MEMORY,
                      attribute_column="grwl_value", where="grwl_value > 0")
    return



def high_res():
    """
    g.region -p res=10

    v.to.rast input=$input_river_mask type="area" output=riverbank__raster use="cat" \
        memory=$MEMORY label_column=grwl_value

    # remove all units that are entirely surrounded by others by constraining the area changes 
    r.grow.distance input=riverbank__raster value=riverbank__value__grown
    r.mapcalc exp="riverbank__value__grown=int(riverbank__value__grown)" --o
    r.mapcalc exp="cell__1=int(1)"

    r.stats.zonal base=riverbank__value__grown cover=cell__1 method=sum output=riverbank__grown__area  # in cells
    r.stats.zonal base=riverbank__raster cover=cell__1 method=sum output=riverbank__area  # in cells

    # now use grwl value from the label
    r.mapcalc exp="riverbank__grwl__highres=if(riverbank__grown__area/riverbank__area < 1.5, null(), @riverbank__raster)"

    # width
    r.grow.distance -n input=riverbank__grwl__highres distance=riverbank__halfwidth__highres

    # resample
    g.region -pa res=30
    r.resamp.stats in=riverbank__halfwidth__highres method=maximum out=riverbank_halfwidth
    r.resamp.stats in=riverbank__grwl__highres method=maximum out=riverbank_grwl
    """

if __name__ == "__main__":
    from commandline import interface as _cli
    _cli()