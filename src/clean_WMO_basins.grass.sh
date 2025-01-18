#!/usr/bin/env bash

export GRASS_OVERWRITE=1

input=$1
output=$2

v.in.ogr $input out=wmobb_basins

# clip dateline overlapping polygons
g.region n=90 s=-90 w=-179.9999 e=179.9999
v.in.region out=world_roi
v.overlay ainput=wmobb_basins binput=world_roi operator=and output=wmobb_basins_clean

# dissolve by basin where there are subbasins, otherwise by name
v.db.addcolumn wmobb_basins_clean column="name varchar(256)"
v.db.update wmobb_basins_clean col=name query_col=a_WMOBB_NAME
large=$(python -c """import grass.script as grass
bas=grass.vector_db_select('wmobb_basins_clean', where=\"a_WMOBB_SUBBASIN != '---'\", column='cat,a_WMOBB_BASIN')['values']
print(str(tuple(set([b for i, (cat, b) in bas.items()]))))""")
v.db.update wmobb_basins_clean col=name query_col=a_WMOBB_BASIN where="a_WMOBB_BASIN in $large"
v.dissolve wmobb_basins_clean col=name out=basins_dissolved

# add region back to table
v.db.addcolumn basins_dissolved column="region varchar(128)"
python -c """
import grass.script as grass
from collections import defaultdict
bas=grass.vector_db_select('wmobb_basins_clean', columns='cat,name,a_REGNAME')['values']
regions = {n: r for c, n, r in bas.values()}
basr=grass.vector_db_select('basins_dissolved', columns='cat,name')['values']
reg = defaultdict(lambda: list())
for c, n in basr.values():
    reg[regions[n]].append(c)
for r, cats in reg.items():
    grass.run_command('v.db.update', map='basins_dissolved', where='cat in (%s)' % (','.join(cats)), column='region', value=r)
    print(r)
"""

v.out.ogr basins_dissolved output_layer=basins out=$output
