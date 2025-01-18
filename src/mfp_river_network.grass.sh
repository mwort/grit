#!/usr/bin/env bash
set -e -u
export GRASS_OVERWRITE=1

# MASK FILLING
MIN_BAR_SIZE=1  # sq.km
MIN_LAKE_ISLAND_SIZE=50
MAX_REACH_LENGTH=1000  # m


# network densification/gap filling
DEM_VERT_ACCURACY=15  # dem units (likely meters)
MIN_STREAM_DAREA=50  # km2
LINE_END_ATTRACTION_BUFFER=90  # m


MAX_CLEANING_ITERATIONS=15
# integer ID maximum - 1 within each domain
LOCAL_ID_SPACE=10000000

# for regional drainage density (not working)
SUBBASIN_TARGET_AREA=150000  #km2, will need to be bounded by domain size
MIN_CENTRELINE_SUBBASIN_CELLS=10000
MAX_WIDTH_THRESHOLD=500  # only fit width - darea relationship to rivers smaller that size
MAX_DAREA_THRESHOLD=5000
MIN_DAREA_THRESHOLD=500


MEMORY=7000  # MB
CLEAN_TMP_MAPS="true"

# minimum elevation of intersections of catchment polygons with grwl mask to be interbasin nodes
INTERBASIN_TRANSFER_MIN_ELEVATION=30  # m
INTERBASIN_COAST_BUFFER=60000  # m

# sink discretisation
SINK_MIN_VOLUME=20  # km3
SINK_DEM_THRESHOLD=10  # m

# station snapping
STATION_DAREA_MAX_DEVIATION=20  # %, max difference in reported and derived drainage area
STATION_MAX_SNAP_DISTANCE=2000  # m

# rivgraph
RIVGRAPH_MAINCHANNEL_DAREA=1000 # km2

# tile size of raster output, must be divisible by raster resolution
RASTER_OUTPUT_GRID_SIZE=300000 # m in domain location, i.e. Equal Earth


clean_tmp(){
    input_do_clean=${1:-$CLEAN_TMP_MAPS}
    if $input_do_clean; then
        echo Cleaning tmp maps...
        g.remove type=all pattern=*__* -f
        # g.remove type=all pattern=*__* -b -f
    fi
    # explicitly vacuum sqlite db to avoid leftover journal files
    db=$(g.gisenv GISDBASE,LOCATION_NAME,MAPSET sep=/)/sqlite/sqlite.db
    if [ -f "$db" ]; then
        sqlite3 $db vacuum
        sleep 20
    fi
}

create_empty_vector(){
    output_vector=$1
    v.edit map=$output_vector tool=create
    v.db.addtable map=$output_vector
}


current_resolution(){
    local res=$(g.region -g | grep nsres)
    echo ${res##*=}
}


river_mask_filling(){
    input_river_mask=$1
    output_mapset=$2  # may exist already

    g.mapset -c $output_mapset
    g.region -p rast=$input_river_mask

    # grown value needed to identify islands in lakes and value of all islands
    r.grow.distance input=$input_river_mask value=rivers__grwl__all

    # identify bar/island sizes
    r.mapcalc exp="rivers__bool=if(isnull($input_river_mask), 0, 1)"
    r.clump in=rivers__bool out=rivers__clumped
    r.mapcalc exp="cellarea__=area()/(1000*1000)"
    r.stats.zonal base=rivers__clumped cover=cellarea__ method=sum output=rivers__clumpsize  # in cells
    # islands in lakes
    r.stats.zonal base=rivers__clumped cover=rivers__grwl__all method=average output=clump__grwl__average  # in cells


    # # fill and grow grwl values into gaps
    r.mapcalc exp="rivers__mask=if(rivers__bool, 1, \
                          if(clump__grwl__average == 180 & rivers__clumpsize<$MIN_LAKE_ISLAND_SIZE, 1, \
                          if(rivers__clumpsize<$MIN_BAR_SIZE, 1, null() ) ) )"
    r.mask --o rast=rivers__mask
    r.mapcalc exp="rivers_filled=int(rivers__grwl__all)"
    r.mask -r

    clean_tmp
}


extent_coastal(){
    input_grwl_filled=$1
    input_land_mask=$2
    input_coast_dist=$3
    output_mapset=$4

    g.mapset -c $output_mapset
    g.region -p rast=$input_grwl_filled

    COASTAL_DISTANCE=60000  # m, roughly the widest expected river
    MAX_WIDTH_EXTENT_FACTOR=3

    # if no coast present
    if [ ! -z $(r.info -s $input_coast_dist | grep sum=NULL) ]; then
        g.copy rast=$input_grwl_filled,rivers_grwl
        exit 0
    fi

    r.mapcalc exp="coastal__water=if(not(isnull($input_land_mask)) & \
        $input_coast_dist < $COASTAL_DISTANCE & $input_grwl_filled != 180, \
        $input_grwl_filled, null())"
    # unique component ids and their maximum width/thickness
    r.clump in=coastal__water out=coastal__water__clumped
    r.grow.distance in=coastal__water distance=coastal__bank__dist -n
    r.stats.zonal -r base=coastal__water__clumped cover=coastal__bank__dist method=max out=coastal__bank__dist__max

    # grow unique values out
    r.mapcalc exp="grwl__land=if(isnull($input_land_mask), null(), if(isnull(coastal__water__clumped), 0, coastal__water__clumped))"
    r.grow.distance input=grwl__land value=grwl__land__grown
    # move max dists to grown ids
    tmpf=$(g.tempfile pid=$$)
    r.stats -ln coastal__bank__dist__max sep="=" out=$tmpf
    r.reclass in=grwl__land__grown rules=$tmpf out=coastal__bank__dist__max__grown
    # GRWL values in grown ids
    r.stats -n coastal__water__clumped,$input_grwl_filled sep="=" out=$tmpf
    r.reclass in=grwl__land__grown rules=$tmpf out=grwl__values__grown
    # extend river by factor of max bank distance (1/2 width) and give it the coastal GRWL ID 
    r.mapcalc exp="rivers__extended=if(isnull($input_land_mask) & $input_coast_dist < \
        (coastal__bank__dist__max__grown * $MAX_WIDTH_EXTENT_FACTOR * 2), grwl__values__grown, null())"
    r.mapcalc exp="rivers_grwl=if(isnull($input_land_mask), int(rivers__extended), int($input_grwl_filled))"
    r.colors map=rivers_grwl color=water

    clean_tmp
}


patch_large_grwl_lakes(){
    input_glwd_lakes=$1
    input_grwl_tiles=$2
    input_grwl=$3
    output_mapset=$4

    LARGE_LAKE_MIN_AREA=230  # km2
    GLWD_SHRINK=2000  # m
    MAX_GAP_FILL=500  # m

    g.mapset -c $output_mapset
    g.region -p rast=$input_grwl

    v.to.lines input=$input_grwl_tiles out=tile__lines
    if [[ ! -z $(v.db.select -c $input_glwd_lakes) ]] ; then
        v.extract in=$input_glwd_lakes out=glwd__selected where="AREA_SKM > 230"
        if [[ ! -z $(v.db.select -c glwd__selected) ]] ; then
            v.buffer in=glwd__selected dist=-$GLWD_SHRINK out=glwd__shrunk
            v.overlay bin=glwd__shrunk ain=tile__lines out=tile_lines_clipped op=and
            echo Found $(v.db.select -c tile_lines_clipped | wc -l) GRWL tile lines in large lakes
        else
            create_empty_vector tile_lines_clipped
        fi
    else
        create_empty_vector tile_lines_clipped
    fi
    if [[ ! -z $(v.db.select -c tile_lines_clipped) ]] ; then
        v.buffer in=tile_lines_clipped dist=$MAX_GAP_FILL out=grwl_gap
        v.to.rast grwl_gap out=grwl_large_lake_gap use=val value=180 memory=$MEMORY
    else
        r.mapcalc exp="grwl_large_lake_gap=null()"
    fi
    clean_tmp
}

patch_coastal_grwl(){
    input_occurrence=$1
    input_grwl=$2
    input_land=$3
    output_mapset=$4

    g.mapset -c $output_mapset
    g.region rast=$input_grwl

    COAST_DIST=100000
    MIN_OCCURRENCE=80

    # only consider 100km strip away from ocean
    r.grow.distance -n in=$input_land dist=distance_to_coast maximum_dist=$COAST_DIST
    
    r.mask distance_to_coast

    r.mapcalc exp="water__=if($input_occurrence>=$MIN_OCCURRENCE, 1, null())"
    r.clump in=water__ out=water_clumped

    # count n coastal cells
    r.mapcalc exp="grwl__coastal=if($input_grwl==126, 1, null())"
    r.stats.zonal base=water_clumped cover=grwl__coastal method=count output=n__coastal

    r.mask -r
    r.mapcalc exp="coastal_water=if(n__coastal>1, 126, null())"

    clean_tmp
}


skeletonise_parallel(){
    input_grwl=$1
    input_cpus=$2  # will be split into twice as many tiles as cpus to utilise idle CPUs due to different speeds of tiles
    output_mapset=$3  # may exist already

    g.mapset -c $output_mapset
    g.region rast=$input_grwl

    # upper threshold, will complete successfully after this even if not properly thinned
    iter=1000

    if [ $input_cpus == 1 ]; then
        r.thin input=$input_grwl output='rivers_thinned' iterations=$iter
    else
        # thin entire mask in parallel
        python -c "
import grass.script as gs
from grass.pygrass.modules.grid import GridModule
x, y = (gs.region()[i] for i in ('cols', 'rows'))
wx, wy = (x//$input_cpus + 1, y//2 + 1) if x > y else (x//2 + 1, y//$input_cpus + 1)
rthin = GridModule('r.thin', input='$input_grwl', output='rivers_thinned', iterations=$iter,
    overlap=$iter, processes=$input_cpus, width=int(wx), height=int(wy), mapset_prefix='rthin$output_mapset')
rthin.run(patch=False, clean=False)
# patch only those with existing output to avoid errors because tiles dont have GRWL cells
maps = [m for m in gs.mapsets() if m.startswith('rthin$output_mapset')]
outp = gs.list_strings('rast', 'rivers_thinned', mapset=maps)
gs.message(f'Found {len(outp)} tiles to patch.')
# copy rasters one by one checking if raster exists
gs.run_command('v.mkgrid', map='parallel__grid', grid=[len(rthin.bboxes), len(rthin.bboxes[0])])
outn = [r.replace('@', '_') for r in outp]
gs.use_temp_region()
for r, n in zip(outp, outn):
    # get row/column from mapset
    row, col = [int(n.split('_')[i])+1 for i in (-2, -1)]
    reg = gs.parse_command('v.db.select', map='parallel__grid', flags='r', where=f'row={row} and col={col}')
    gs.run_command('g.region', flags='a', **reg)
    gs.mapcalc(exp=f'{n}={r}')
gs.del_temp_region()
gs.run_command('r.patch', input=outn, output='rivers_thinned')
gs.run_command('g.remove', type='rast', name=outn, flags='f')
rthin.clean_location()
    "
    fi
    
    clean_tmp
}

build_network(){
    input=$1
    output=$2

    v.build.polylines input=$input output=${input%@*}__poly type=line
    v.category in=${input%@*}__poly op=add type=line layer=1 out=${input%@*}__poly__cat
    v.net -c in=${input%@*}__poly__cat out=$output operation=nodes
    v.db.droptable -f map=$output layer=1 2> /dev/null || echo No table to remove.
    v.db.droptable -f map=$output layer=2 2> /dev/null || echo No table to remove.
    v.db.addtable map=$output layer=1 table=${output}_lines
    v.db.addtable map=$output layer=2 table=${output}_nodes
}

snap_points(){
    input_from_points=$1
    input_to_vect=$2
    output_points=$3

    if [ $(v.info $input_from_points -t | grep map3d | cut -d= -f2) == 1 ]; then
        v.to.3d -r $input_from_points out=points__snapp
    else
        g.copy vect=$input_from_points,points__snapp
    fi

    v.db.addcolumn points__snapp \
        column='segment_id int,snapped_x double,snapped_y double,snapped_distance double'
    v.distance from=points__snapp to=$input_to_vect upload=cat,to_x,to_y,dist \
        column='segment_id,snapped_x,snapped_y,snapped_distance'
    # create new points vector with all attributes transfered
    v.db.select points__snapp column="cat,snapped_x,snapped_y" -c | \
        v.in.ascii input=- out=$output_points cat=1 x=2 y=3

    # move table over
    v.db.addtable $output_points
    v.db.join map=$output_points col=cat other_table=points__snapp other_col=cat
}


vectorize(){
    input_thin_rivers=$1
    input_bank_distance=$2
    input_coast=$3
    output_mapset=$4  # may exist already

    g.mapset -c $output_mapset
    g.region rast=$input_bank_distance

    { # try
        r.to.vect $input_thin_rivers output=rivers__ type=line
    } || { # catch
        skeletonise_parallel $input_thin_rivers 2 $output_mapset
        r.to.vect rivers_thinned output=rivers__ type=line
    }

    build_network rivers__ rivers__network

    # clean using network approach
    v.to.db rivers__network type=line option=length column=length units=meters
    v.what.rast map=rivers__network layer=2 raster=$input_bank_distance col=river_halfwidth
    v.what.rast map=rivers__network layer=2 raster=$input_coast col=is_outlet
    v.db.update map=rivers__network layer=2 col=is_outlet where="is_outlet IS NULL" value=0
    # only do basic cleaning, dangles of narrow rivers and lakes are needed for gap filling
    river_network_utils.py clean_network rivers__network rivers_cleaned --factor=1.4
    
    clean_tmp
}


land_ocean_coast(){
    input_land_vector=$1
    output_mapset=$2

    COAST_WIDTH=10  # cells

    g.mapset -c $output_mapset

    # remove small island and merge smaller adjacent areas with larger ones
    v.dissolve in=$input_land_vector out=land__dissolved
    v.clean in=land__dissolved out=land tool=break,bpol,rmdupl,rmarea thresh=0,0,0,$(($MIN_STREAM_DAREA*1000000)) -c

    # land, ocean, coast rasters
    v.to.rast in=land use=val out=land memory=$MEMORY
    r.mapcalc exp="ocean=if(isnull(land), 1, null())"
    r.grow input=ocean out=ocean__with__coast radius=$COAST_WIDTH new=1
    r.mapcalc exp="coast=if(isnull(ocean), ocean__with__coast, null())"
    # distance to coast into both land and ocean
    r.grow.distance in=coast distance=coast_distance__fcell
    r.mapcalc exp="coast_distance=int(coast_distance__fcell)"

    clean_tmp
}


interbasin_nodes(){
    # Identify intersection points of external catchment polygons with the grwl lines as potential
    # interbasin canals. Where the catchment polygons are wrong (mostly in deltas), the points need to be deleted.
    input_catchments=$1
    input_grwl_lines=$2
    input_elevation=$3
    input_coast_dist=$4
    output_mapset=$5

    g.mapset -c $output_mapset

    v.overlay ain=$input_grwl_lines atype=line bin=$input_catchments out=line_polygon__overlay op=and

    # only continue if any overlaps
    if [ -z "$(v.db.select line_polygon__overlay -c)" ]; then
        create_empty_vector interbasin_nodes
        exit 0
    else
        v.net input=line_polygon__overlay output=line_polygon__overlay__net op=nodes -c

        v.db.addtable line_polygon__overlay__net col="dist double" layer=2
        v.distance from=line_polygon__overlay__net from_layer=2 from_type=point to=$input_catchments \
            to_type=boundary dmax=0 upload=dist column=dist
    fi

    # only continue if there are any intersection points
    if [ -z "$(v.db.select line_polygon__overlay__net where="dist=0" layer=2 -c)"]; then
        create_empty_vector interbasin_nodes
    else
        v.extract line_polygon__overlay__net layer=2 type=point where="dist=0" output=line_polygon__nodes
        v.category line_polygon__nodes op=chlayer layer=2,1 output=intersection_nodes
        v.db.connect intersection_nodes layer=1 table=intersection_nodes
        v.db.connect intersection_nodes layer=2 -d

        # interbasin defined by elevation threshold
        v.what.rast intersection_nodes col=elevation rast=$input_elevation
        v.what.rast intersection_nodes col=coast_distance rast=$input_coast_dist
        v.extract intersection_nodes out=interbasin_nodes \
            where="elevation>$INTERBASIN_TRANSFER_MIN_ELEVATION AND coast_distance>$INTERBASIN_COAST_BUFFER" -t
        v.db.addtable interbasin_nodes
    fi
    
    clean_tmp
}


dem_conditioning(){
    input_halfwidth=$1
    input_center_line=$2
    input_osm_centreline_raster=$3
    input_osm_centerline_raster_small=$4
    input_grwl=$5
    output_mapset=$6  # may exist already

    g.mapset -c $output_mapset
    g.region rast=$input_halfwidth -p

    # conversions (mapset and resolution)
    v.to.rast input=$input_center_line use=cat out=centerlines memory=$MEMORY
    #r.resamp.stats input=$input_halfwidth output=river_width method=maximum
    # lower surrounding of line ends/starts to make flow gravitate to them and help bridge small gaps better
    build_network $input_center_line centreline_segments
    river_network_utils.py upload_line_counts centreline_segments
    v.to.rast in=centreline_segments layer=2 where="n_lines=1" out="centreline_ends" use=cat memory=$MEMORY
    r.grow.distance in=centreline_ends distance=line_end_zone max=$LINE_END_ATTRACTION_BUFFER

    # DEM combined lowering
    r.mapcalc << EOF
    eval( \
    bankdistance_gradient = if($input_grwl == 180, log($input_halfwidth/ewres()), $input_halfwidth/ewres()) , \
    bankdistance_lowering = if(isnull(bankdistance_gradient), 0, bankdistance_gradient) , \
    centerline_lowering = if(isnull(centerlines), 0, $input_halfwidth/(ewres()*2) + 5) , \
    osm_lowering = if(isnull($input_halfwidth) & not(isnull($input_osm_centreline_raster)), 2, 0) , \
    osm_lowering_small = if(isnull($input_halfwidth) & not(isnull($input_osm_centerline_raster_small)), 1, 0) , \
    line_end_lowering = if(isnull(line_end_zone), 0, ($LINE_END_ATTRACTION_BUFFER-line_end_zone)/$LINE_END_ATTRACTION_BUFFER) * 3 \
    )
    dem_lowering = (bankdistance_lowering + centerline_lowering + osm_lowering + osm_lowering_small + line_end_lowering) * $DEM_VERT_ACCURACY
EOF

    clean_tmp
}


sinks_merit(){
    input_merit_flow_direction=$1
    input_rivers_grwl=$2
    input_lake_ids=$3
    input_glwd_lakes=$4
    output_mapset=$5  # may exist already

    MERIT_SINKS_BUFFER=180  #m

    g.mapset -c $output_mapset
    g.region rast=$input_rivers_grwl -p

    # take MERIT sinks and buffer
    r.mapcalc exp="sinks__merit=if(($input_merit_flow_direction==255), 1, null())"
    g.message mes="$(echo Found $(r.stats sinks__merit -cn) MERIT sink cells)"

    r.grow -m in=sinks__merit out=sinks__merit__grown radius=$MERIT_SINKS_BUFFER

    # check GRWL lakes that have MERIT sinks
    r.stats.zonal base=$input_lake_ids cover=sinks__merit__grown out=lakes__nsinks method=count

    # GLWD lakes only those with MGLD_TYPE closed or probably closedx
    if [ -z "$(v.db.select $input_glwd_lakes -c)" ]; then
        r.mapcalc exp="glwd__closed=null()"
    else
        v.to.rast in=$input_glwd_lakes out=glwd__closed where="MGLD_TYPE in ('closed', 'closedx')" use=cat
    fi
    # patch
    r.mapcalc exp="sinks = if( (not(isnull(lakes__nsinks)) &&& lakes__nsinks > 0) ||| not(isnull(glwd__closed)), 1, sinks__merit__grown)"

    # vectorise to later break network
    r.to.vect in=sinks out=sinks type=area

    clean_tmp
}


sinks(){
    input_dem=$1
    input_rivers_grwl=$2
    input_lake_ids=$3
    input_glwd_lakes=$4
    input_memory=$5
    output_mapset=$6  # may exist already

    g.mapset -c $output_mapset
    g.region rast=$input_rivers_grwl -p

    # burn grwl and neighbouring cells into DEM (also copy dem raster to this mapset, otherwise r.hydrodem wont read it)
    r.mapcalc exp="dem__=if(isnull($input_rivers_grwl) & \
        isnull($input_rivers_grwl[-1,0]) & isnull($input_rivers_grwl[1,0]) & \
        isnull($input_rivers_grwl[0,-1]) & isnull($input_rivers_grwl[0,1]) & \
        isnull($input_rivers_grwl[-1,-1]) & isnull($input_rivers_grwl[1,1]) & \
        isnull($input_rivers_grwl[-1,1]) & isnull($input_rivers_grwl[1,-1]), $input_dem, 0)"
    # remove all (-a) sinks by filling (-f)
    r.hydrodem -af in=dem__ out=depressionless__dem mem=$input_memory

    # only consider filling above min threshold, take difference against original DEM (not burned)
    r.mapcalc exp="depression__filling=if(depressionless__dem - $input_dem > $SINK_DEM_THRESHOLD, (depressionless__dem - $input_dem), null())"
    r.mapcalc exp="depression__area=if(depression__filling, 1, null())"

    r.clump -d input=depression__area output=depression__ids
    r.stats.zonal base=depression__ids cover=depression__filling out=depression__volume method=sum
    r.stats.quantile base=depression__ids cover=$input_dem percentile=2, out=depression__min__elevation
    r.mapcalc exp="dem__sinks=if(depression__volume*area()/(10^9) > $SINK_MIN_VOLUME &&& $input_dem <= depression__min__elevation, depression__ids, null())"

    # check GRWL lakes that have sinks and add entire lake
    r.stats.zonal base=$input_lake_ids cover=dem__sinks out=lakes__nsinks method=count

    # GLWD lakes only those with MGLD_TYPE closed or probably closedx
    if [ -z "$(v.db.select $input_glwd_lakes -c)" ]; then
        r.mapcalc exp="glwd__closed=null()"
    else
        v.to.rast in=$input_glwd_lakes out=glwd__closed where="MGLD_TYPE in ('closed', 'closedx')" use=cat mem=$input_memory
    fi

    # patch
    r.mapcalc exp="sinks = if( (not(isnull(lakes__nsinks)) &&& lakes__nsinks > 0) ||| not(isnull(glwd__closed)), 1, dem__sinks)"

    # vectorise to later break network
    r.to.vect in=sinks out=sinks type=area
    # add sink volume/depth
    v.what.rast sinks rast=depression__volume col=volume type=centroid
    v.db.addcolumn sinks col="volume_km3 double"
    v.db.update sinks col=volume_km3 query="volume*30*30/1000000000"
    v.db.dropcolumn sinks col=volume
    v.what.rast sinks rast=depression__filling col=depth_m type=centroid

    g.remove rast name=dem__ -f
    clean_tmp
}


r_watershed_sfd(){

    input_dem=$1
    input_dem_lowering=$2
    input_sinks=$3
    memory=$4
    output_mapset=$5  # may exist already

    g.mapset -c $output_mapset
    g.region rast=$input_dem -p

    r.mapcalc exp="dem__burned=if(isnull($input_dem), 0, $input_dem) - $input_dem_lowering"

    subbasin_ncells=$(python -c "import grass.script as g; print(int(float($SUBBASIN_TARGET_AREA)/(g.region()['nsres']**2 * 1e-6)))")
    # decide whether to use in-memory version
    ncells=$(g.region -g | grep cells)
    #inmem=$(if ((${ncells##*=} > 2**31 | ${ncells##*=}*32/10**6 > $memory)); then echo -m memory=$memory; fi)
    inmem="-m memory=$memory"
    #g.gisenv set="DEBUG=2"
    r.watershed -s $inmem elevation=dem__burned drainage=drainage accumulation=accumulation__cells \
        depression=$input_sinks
        #basin=subbasins
    # parallel version
    #r_watershed_parallel.grass.py -k elevation=dem__burned accumulation=accumulation \
    #    drainage=drainage segments=segments overlap=250000 \
    #    memory=$memory processes=$n_processes r_watershed_processes=$(($n_processes/4)) --o
    # accumulation cells to drainage area
    r.mapcalc exp="accumulation_sqkm=abs(accumulation__cells*area()/1000000)"
    r.colors -n -a map=accumulation_sqkm color=inferno
    # report negative accumulation for quality assurance
    r.mapcalc exp="accumulation_negative_sqkm=if(accumulation__cells < 0, accumulation_sqkm, null())"
    
    clean_tmp
}


gap_filled_network(){
    input_accumulation_sqkm=$1
    input_high_res_river_mask=$2
    input_grwl_lines=$3
    input_bankdistance=$4
    input_coast=$5
    input_cpus=$6
    output_mapset=$7  # may exist already

    g.mapset -c $output_mapset
    g.region -p rast=$input_accumulation_sqkm
    res=$(current_resolution)

    ## crucial to resample to have low res river mask that covers a wider mask
    #r.resamp.stats in=$input_high_res_river_mask out=low_res_river_mask
    #r.to.vect $input_high_res_river_mask type=area out=river_mask
    # remove single cells
    #v.clean river_mask tool=rmarea thresh=$(($res*$res)) output=river_mask__clean
    # create streams rasters
    r.mapcalc exp="streams=if(abs($input_accumulation_sqkm) > $MIN_STREAM_DAREA, int(abs($input_accumulation_sqkm)), null())"
    skeletonise_parallel streams $input_cpus $output_mapset  # whipes all tmp maps in this mapset
    r.to.vect rivers_thinned out=streams type=line

    # remove all loops
    tmpf=$(g.tempfile pid=$$)
    echo start > $tmpf
    iteration=0
    loopy=streams
    while [ ! -z "$(cat $tmpf)" ] && (($iteration < $MAX_CLEANING_ITERATIONS)); do
        build_network $loopy streams__network
        v.what.rast streams__network layer=2 raster=$input_accumulation_sqkm column=darea
        v.to.db streams__network layer=1 option=length col=length
        river_network_utils.py remove_cycles streams__network streams__nocycles --drop_lines_file=$tmpf
        echo "Finished loop purge iteration $iteration"
        loopy=streams__nocycles
        iteration=$(($iteration + 1))
    done

    # v.overlay ain=streams bin=river_mask__clean out=streams_no_grwl operator=not
    # v.to.points input=streams_no_grwl out=nodes use=node
    # # fix cats caused by v.to.points
    # v.category nodes option=del cat=-1 out=nodes_nocats
    # v.category nodes_nocats layer=2,1 option=chlayer out=nodes_cats
    # # select nodes laying on boundary of river mask
    # v.select ain=nodes_cats bin=river_mask__clean out=touch_nodes

    # v.db.addcolumn touch_nodes column="segment_cat int,snap_dist double,segment_x double,segment_y double"
    # v.distance from=touch_nodes to=$input_grwl_lines output=missing_bits \
    #     upload=cat,dist,to_x,to_y column=segment_cat,snap_dist,segment_x,segment_y
    # v.extract input=missing_bits output=missing_bits_linesonly type=line
    # # add nodes at snap intersections (done done by v.clean break)
    # # the start of missing_bits is already a node because it's the end of streams_no_grwl
    # v.to.points -t missing_bits_linesonly out=snap_nodes use=end
    # v.net in=$input_grwl_lines operation=connect points=snap_nodes out=grwl_lines_broken thresh=0
    # v.extract in=grwl_lines_broken out=grwl_lines_only type=line
    # v.patch in=grwl_lines_only,streams_no_grwl,missing_bits_linesonly out=rivers__dupl
    # # clean the overlapping nodes
    # v.clean input=rivers__dupl tool=snap,break,rmdupl,rmdangle \
    #     thresh=1,0,0,${res} output=rivers

    # patch all streams with rivers
    v.patch in=$input_grwl_lines,streams__nocycles out=rivers__dupl
    v.clean input=rivers__dupl tool=snap,break,rmdupl,rmdangle \
        thresh=1,0,0,${res} output=rivers__streams

    # make sure lines are really broken where patched lines exactly end along other lines
    v.to.points in=rivers__streams out=rivers__streams__nodes use=node
    v.net in=rivers__streams points=rivers__streams__nodes op=connect out=rivers__streams__net -c thresh=0
    v.category input=rivers__streams__net option=add cat=1 step=1 layer=3 output=rivers__streams__net__recat type=line
    v.extract rivers__streams__net__recat type=line layer=3 out=rivers__streams__net__lines
    v.category input=rivers__streams__net__lines op=chlayer layer=3,1 out=rivers__streams__net__cat1

    # snap small gaps not filled with streams, in two iterations to snap short isolated lines
    netholes=rivers__streams__net__cat1
    for i in $(seq 2); do
        build_network $netholes river__lines__network__$i
        v.what.rast map=river__lines__network__$i layer=2 column=grwl_value raster=$input_high_res_river_mask
        v.what.rast map=river__lines__network__$i layer=2 raster=$input_bankdistance col=river_halfwidth
        # set nodes in coast to different grwl_value to skip those
        v.what.rast map=river__lines__network__$i layer=2 raster=$input_coast col=is_coast
        v.db.update map=river__lines__network__$i layer=2 where="is_coast is not null" col=grwl_value value=0
        # avoid 3 pixel wide river bias by enforcing a minimum width
        v.db.update map=river__lines__network__$i layer=2 col=river_halfwidth \
            where="river_halfwidth < 90 OR river_halfwidth IS NULL" value=90
        v.db.addcolumn river__lines__network__$i col='length double' layer=1
        v.to.db river__lines__network__$i option=length column="length"
        # snap endpoints
        river_network_utils.py snap_network_ends river__lines__network__$i river__lines__snapped__$i
        # save added lines
        v.extract river__lines__snapped__$i where="is_valid=1" out=snap_lines_$i
        netholes=river__lines__snapped__$i
    done
    g.rename vect=$netholes,river__lines__snapped

    # remove weak links created from the merge of vectorised centrelines and D8 streams
    # first filter by length
    v.to.db river__lines__snapped layer=1 option=length col=length
    # min legth is one straight and one diagonal resolution
    minl=$(python -c "import math; print(round(3 * math.sqrt(2 * ($res)**2)))")
    v.db.droprow river__lines__snapped where="length>$minl" output=shorts__
    v.net -c input=shorts__ output=shorts__net op=nodes
    # remove bridges/single lines
    v.net.bridge in=shorts__net out=shorts__bridges method=bridge
    v.db.addtable shorts__bridges
    tmpf=$(g.tempfile pid=$$)
    v.db.select shorts__bridges -c col=cat > $tmpf
    v.extract -r shorts__net file=$tmpf output=shorts__nobridges
    # create network and add node attributes (darea and components)
    v.net -c input=shorts__nobridges output=triangles__net op=nodes
    v.db.addtable triangles__net layer=2
    v.what.rast triangles__net layer=2 raster=$input_accumulation_sqkm column=darea
    v.net.components in=triangles__net out=triangles__net__comp method=weak
    v.db.addcolumn triangles__net layer=2 col="comp int"
    v.distance from=triangles__net from_layer=2 to=triangles__net__comp dmax=0 upload=to_attr col=comp to_column=comp
    # filter out lines in python
    river_network_utils.py prune_triangles triangles__net river__lines__snapped rivers


    # recompute half-width and width for rivgraph input, wrong outside of river mask
    g.region -p rast=$input_high_res_river_mask
    v.to.rast in=rivers out=rivers use=val value=1 memory=$(($MEMORY*$input_cpus))
    r.patch in=rivers,$input_high_res_river_mask out=river_mask_gaps
    r.mask rast=river_mask_gaps
    # river half-width
    r.grow.distance -n input=river_mask_gaps distance=rivers__halfwidth
    # width and correct for half width of single pixel wide rivers should be 15m
    r.mapcalc exp='river_width=rivers__halfwidth*2 - nsres()'
    r.mapcalc exp='river_halfwidth=river_width/2'
    r.mask -r
    
    clean_tmp
}


line_grwl_intersection_nodes(){
    input_network=$1
    input_grwl_raster=$2
    output_nodes=$3

    DIST_NEXT_NODE=90  # m

    r.to.vect in=$input_grwl_raster out=intersection__polygons type=area

    # add intersections as nodes to network
    v.overlay ain=$input_network atype=line bin=intersection__polygons out=line_polygon__overlay op=and
    v.net input=line_polygon__overlay output=line_polygon__overlay__net op=nodes -c

    v.db.addtable line_polygon__overlay__net col="dist double,cluster int" layer=2
    v.distance from=line_polygon__overlay__net from_layer=2 from_type=point to=intersection__polygons \
        to_type=boundary dmax=0 upload=dist column=dist
    # cluster nodes to be able to only pick out one
    v.extract line_polygon__overlay__net layer=2 out=line_polygon__overlay__nodes
    v.cluster line_polygon__overlay__nodes out=line_polygon__overlay__nodes__clustered layer=3 min=2 dist=$DIST_NEXT_NODE
    v.distance from=line_polygon__overlay__net from_layer=2 to=line_polygon__overlay__nodes__clustered to_layer=3 upload=cat col=cluster

    # pinch points like =><= will also create undesired overlaps with the polygon here, so we need to check the GRWL value on either side of point
    python -c '
from grass_attribute_table import GrassAttributeTable
from river_network_utils import network_nodes
nodes = network_nodes("line_polygon__overlay__net")
nodestbl = GrassAttributeTable("line_polygon__overlay__net", layer=2)
linestbl = GrassAttributeTable("line_polygon__overlay__net", layer=1)
# if no dist, then its a regular network node
nlcats = nodes.loc[nodestbl.dist.dropna().index, "line_cats"]
nvals = nlcats.apply(lambda c: linestbl.loc[c, "b_value"].unique())
nodestbl.loc[nvals.index, "n_grwl"] = nvals.apply(len)
nodestbl.loc[nvals.index, "grwl_vals"] = nvals.apply(lambda d: ",".join(map(str, d)))
# take first in cluster or if a regular node is part of cluster, dont take any
nodestbl.cluster=nodestbl.cluster.fillna(0).astype(int)
cln = nodestbl[(nodestbl.cluster>0) & (nodestbl.n_grwl>1)]
# first cats of clusters without real junction nodes
cats = cln.reset_index().groupby("cluster").first().cat
no_junction = cats.index.difference(set(nodestbl.cluster[(nodestbl.cluster > 0) & nodestbl.dist.isna()]))
nodestbl.loc[cats[no_junction].values, "cluster"] = 0
nodestbl.write()
'

    v.extract line_polygon__overlay__net layer=2 type=point where="n_grwl>1 AND cluster=0" output=line_polygon__nodes
    v.category line_polygon__nodes op=chlayer layer=2,1 output=$output_nodes
    v.db.connect $output_nodes layer=1 table=$output_nodes
    v.db.connect $output_nodes layer=2 -d
}


build_river_network(){
    # build segment network
    input_river_lines=$1
    input_river_width=$2
    input_grwl_mask=$3
    input_accumulation_sqkm=$4
    input_interbasin_nodes=$5
    input_coast=$6
    input_sinks=$7
    output_mapset=$8  # may exist already

    g.mapset -c $output_mapset
    g.region -p rast=$input_grwl_mask

    # add coast and sinks needed at nodes to avoid trimming links
    v.to.rast in=$input_sinks out=sinks__rast use=val memory=$MEMORY value=2
    r.patch in=sinks__rast,$input_coast out=outlet__ memory=$MEMORY
    # make vector again to later overlap with invalid network components
    r.to.vect in=outlet__ out=outlet_area type=area column=type

    # remove longest line of short loops
    build_network $input_river_lines segments__with__small_loops
    v.what.rast map=segments__with__small_loops layer=2 raster=$input_river_width col=river_width
    v.db.addcolumn segments__with__small_loops col='length double' layer=1
    v.to.db segments__with__small_loops option=length column="length"
    river_network_utils.py clean-small-loops segments__with__small_loops segments__no__loops \
        --factor=2
    toclean=segments__no__loops

    # clean more aggressively, multiple iterations to reconnect cleaned lines and check again
    toclean=segments__no__loops
    tmpf=$(g.tempfile pid=$$)
    echo start > $tmpf
    iteration=0
    while [ ! -z "$(cat $tmpf)" ] && (($iteration < $MAX_CLEANING_ITERATIONS)); do
        build_network $toclean segments__dirty
        v.what.rast map=segments__dirty layer=2 raster=$input_grwl_mask col=grwl_value
        v.what.rast map=segments__dirty layer=2 raster=$input_river_width col=river_width
        # avoid 3 pixel wide river bias by enforcing a minimum width
        v.db.update map=segments__dirty layer=2 col=river_width where="river_width < 90" value=90
        v.what.rast map=segments__dirty layer=2 raster=$input_accumulation_sqkm col=darea
        # preserve outlet nodes by giving them minimum stream darea
        v.what.rast map=segments__dirty layer=2 raster=outlet__ col=is_outlet
        v.db.update map=segments__dirty layer=2 column=is_outlet value=0 where="is_outlet IS NULL"
        #v.db.update segments__dirty layer=2 column=darea \
        #    where="darea<$MIN_STREAM_DAREA AND is_outlet=1" \
        #    value=$MIN_STREAM_DAREA
        v.db.addcolumn segments__dirty col='is_lake int DEFAULT 0' layer=2
        v.db.update segments__dirty layer=2 column=is_lake value=1 where="grwl_value=180 AND is_outlet=0"
        v.db.addcolumn segments__dirty col='length double' layer=1
        v.to.db segments__dirty option=length column="length"
        river_network_utils.py clean_network segments__dirty segments_cleaned__$iteration \
            --node_min_length_column=river_width --factor=2 --min_accum=$MIN_STREAM_DAREA \
            --node_lake_column=is_lake --lake_factor=10 --node_outlet_column=is_outlet \
            --cleaned_cats_file=$tmpf
        toclean=segments_cleaned__$iteration
        echo Finished cleaning iteration $iteration
        iteration=$(($iteration + 1))
    done

    # remove longest line of short loops again
    build_network $toclean segments__with__small_loops
    v.what.rast map=segments__with__small_loops layer=2 raster=$input_river_width col=river_width
    v.db.addcolumn segments__with__small_loops col='length double' layer=1
    v.to.db segments__with__small_loops option=length column="length"
    river_network_utils.py clean-small-loops segments__with__small_loops segments__no__loops \
        --factor=2
    toclean=segments__no__loops

    # remove short dangles
    build_network $toclean segments__with__short
    v.what.rast map=segments__with__short layer=2 raster=outlet__ col=is_outlet
    v.db.update segments__with__short col=is_outlet val=0 where="is_outlet IS NULL" la=2
    v.db.addcolumn segments__with__short col='length double' layer=1
    v.to.db segments__with__short option=length column="length"
    river_network_utils.py clean-short-dangles segments__with__short segments__no__short__dangles \
        --min-length=$(($MAX_REACH_LENGTH / 2))
    toclean=segments__no__short__dangles

    # remove sinks
    if [ ! -z "$(v.db.select -c $input_sinks col=cat)" ]; then
        v.buffer in=$input_sinks out=sinks__shrunk dist=$((-$(current_resolution)/2))
        v.overlay ain=$toclean bin=sinks__shrunk operator=not output=${toclean}__nosinks
        toclean=${toclean}__nosinks
    fi

    # build network again after cleaning
    build_network $toclean segments__poly

    # create dense segment points
    v.to.rast input=segments__poly use=val out=segments__vertices memory=$MEMORY
    r.to.vect in=segments__vertices out=segments__vertices type=point

    # potential outlet points to avoid confluence outlet problem
    r.mapcalc exp="outlet__segments=if(isnull(outlet__) + isnull(segments__vertices) == 0, segments__vertices, null())"
    r.to.vect in=outlet__segments out=potential_outlet_points type=point

    # add grwl change intersection nodes
    line_grwl_intersection_nodes segments__poly $input_grwl_mask grwl_boundary__nodes__unsnapped
    snap_points grwl_boundary__nodes__unsnapped segments__vertices grwl_boundary_nodes
    v.net -s in=segments__poly out=segments__with__grwl_breaks points=grwl_boundary_nodes operation=connect thresh=0
    # add interbasin nodes
    if [ -z "$(v.db.select -c $input_interbasin_nodes col=cat)" ]; then
        g.rename vect=segments__with__grwl_breaks,segments__badline_cats --o
    else
        snap_points $input_interbasin_nodes segments__vertices input_interbasin_nodes__snapped
        res=$(current_resolution)
        v.what.rast map=input_interbasin_nodes__snapped rast=$input_coast col=is_coast
        v.edit map=input_interbasin_nodes__snapped tool=delete where="snapped_distance>$res OR is_coast IS NOT NULL"
        v.net -s in=segments__with__grwl_breaks out=segments__badline_cats \
            points=input_interbasin_nodes__snapped operation=connect thresh=$(($res - 1))
    fi

    # clean categories after splits
    v.category in=segments__badline_cats option=del cat=-1 layer=1 output=segments__badline_cats__no_l1
    v.category in=segments__badline_cats__no_l1 option=del cat=-1 layer=2 output=segments__badline_cats__no_l12
    v.category input=segments__badline_cats__no_l12 option=add cat=1 step=1 layer=1 type=line output=segments__badcats__l1
    v.category input=segments__badcats__l1 option=add cat=1 step=1 layer=2 type=point output=segments
    v.db.droptable segments layer=1 -f
    v.db.addtable map=segments layer=1 table=segments_lines
    v.db.droptable segments layer=2 -f
    v.db.addtable map=segments layer=2 table=segments_nodes \
        column="type varchar(126),grwl_change varchar(126),interbasin_cat int"

    # add attributes
    v.distance from=segments from_layer=2 to=grwl_boundary_nodes dmax=0 upload=to_attr column=grwl_change to_column=grwl_vals
    v.db.update segments layer=2 where="grwl_change IS NOT NULL" col=type val=grwl_change
    if [ -z "$(v.db.select -c $input_interbasin_nodes col=cat)" ]; then
        echo No interbasin node added.
    else
        v.distance from=segments from_layer=2 to=$input_interbasin_nodes dmax=$((30-1)) upload=cat column=interbasin_cat
        v.db.update segments layer=2 where="interbasin_cat IS NOT NULL" col=type val=interbasin
        # drop from interbasin nodes outward facing lines to split network
        river_network_utils.py split_interbasin_network segments
        # make line cats continuous again
        v.category segments op=del out=seg__no__line_cats cat=-1 layer=1
        v.category seg__no__line_cats op=add layer=1 output=segments type=line --o
        v.db.droptable segments -f 
        v.db.addtable segments
    fi

    # add coast and sinks back to trimmed network
    v.what.rast map=segments layer=2 raster=$input_coast col=coast
    v.db.addcolumn segments layer=2 col="sink_id int"
    if [ ! -z "$(v.db.select -c $input_sinks col=cat)" ]; then
        v.what.vect segments layer=2 col=sink_id query_map=$input_sinks query_col=cat
    fi

    clean_tmp
}


select_network_components(){
    input_segments=$1
    input_river_width=$2
    input_grwl_mask=$3
    input_accumulation_sqkm=$4
    input_dem=$5
    input_domain=$6
    input_potential_outlet_points=$7  # to fix confluence at outlet problem
    output_mapset=$8  # may exist already

    g.mapset -c $output_mapset
    g.region -p rast=$input_grwl_mask

    g.copy vect=$input_segments,segments__all

    # add attributes
    river_network_utils.py upload_line_counts segments__all --column="n_segments"
    v.what.rast map=segments__all layer=2 raster=$input_accumulation_sqkm col=darea
    # decide what nodes are outlets
    # if dangle node with grwl = 126 is outlet:
    #  - leads to networks in costal zones producing many outlets and bifurcations
    #  - some real outlets that dont reach the ocean in GRWL mask are correctly classified
    #v.db.update segments__all layer=2 where="(coast!='' OR sink_id!='' OR grwl_value=126) AND n_segments=1" col=type value="outlet"
    # if only coastal and sink dangles are outlets:
    #  - requires costal GRWL gap be fixed as lots of networks wont have an outlet
    v.db.update segments__all layer=2 where="(coast IS NOT NULL OR sink_id IS NOT NULL) AND n_segments=1" col=type value="outlet"
    # add network components id (separate networks)
    v.net.components in=segments__all out=segments__components method=weak
    # v.net.components only retains layer 1/lines with cleaned categories (continuous 1...N)
    # i.e. the following join only works if segments__all lines also have clean cats
    cats=$(v.category segments__all op=report -g | grep "1 line")
    if (( $(echo $cats | cut -d" " -f4) != 1 || $(echo $cats | cut -d" " -f3) != $(echo $cats | cut -d" " -f5) )); then
        echo "segments__all doesnt seem to have continuous cats (1...N): $cats"
        exit 1
    fi
    v.db.join segments__all col=cat other_table=$(v.db.connect segments__components -g layer=1 | grep 1/ | cut -d"|" -f2) \
        other_col=cat subset=comp
    # transfer to nodes
    v.db.addcolumn segments__all layer=2 col="comp int"
    v.what.vect map=segments__all layer=2 column=comp query_map=segments__components query_col=comp
    # drop nodes that are not aligned with lines (likely from interbasin or grwl change nodes)
    echo Remove $(v.db.select segments__all where="comp IS NULL" la=2 -c | wc -l) nodes that dont have a component.
    v.edit map=segments__all tool='delete' layer=2 where="comp IS NULL"
    db.execute sql="DELETE FROM segments__all_segments_nodes WHERE comp IS NULL;"

    # drop components with not the majority of nodes in domain
    river_network_utils.py network_components_in_domain segments__all $input_domain --output=segments__domain

    # check if we need to include any components with confluence outlets
    g.copy vect=$input_potential_outlet_points,potential__outlet__points
    v.db.addcolumn potential__outlet__points col="comp int"
    v.what.vect map=potential__outlet__points query_map=segments__domain query_column=comp column=comp
    v.what.rast map=potential__outlet__points raster=$input_accumulation_sqkm col=darea
    # find potential outlets in components without outlets, take the one with the largest darea
    no_outlet_cats=$(python -c "import grass_attribute_table as gat; import warnings
comps = gat.GrassAttributeTable('segments__domain', layer=2).groupby('comp').type.agg(lambda x: 'outlet' not in x.unique())
outl = gat.GrassAttributeTable('potential__outlet__points').dropna(subset=['comp'])
outl.comp = outl.comp.astype(int)
outlno = outl[outl.comp.isin(comps.index[comps]) & (outl.darea >= $MIN_STREAM_DAREA)].groupby('comp').darea.idxmax()
print(','.join(outlno.astype(str)))
")
    nconoutlets=$(echo $no_outlet_cats | cut -d, --output-delimiter=" " -f1- | wc -w)
    echo Found $nconoutlets confluence outlets.
    if (($nconoutlets > 0)); then
        v.extract in=potential__outlet__points cats=$no_outlet_cats out=confluence_outlet
        catmax=$(v.category segments__domain op=report -g la=2 | grep "2 point" | cut -d" " -f5)
        v.category in=confluence_outlet out=confluence__outlet op=sum cat=$catmax
        # add points to network
        v.net -s input=segments__domain points=confluence__outlet output=segments__domain__with_conout operation=connect threshold=1
        # fix new points cats + attributes
        v.to.db --o map=segments__domain__with_conout layer=2 option=cat columns=cat query_layer=2
        v.db.update segments__domain__with_conout layer=2 column=type val=outlet where="n_segments IS NULL"
        v.what.vect segments__domain__with_conout layer=2 column=comp query_map=confluence_outlet query_col=comp
        v.what.vect segments__domain__with_conout layer=2 column=darea query_map=confluence_outlet query_col=darea
        v.db.update segments__domain__with_conout layer=2 column=n_segments val=2 where="n_segments IS NULL"
        # fix new lines cats + attributes
        duplines=$(v.category opt=print in=segments__domain__with_conout type=line | sort -n | uniq -d)
        tbl=$(v.db.connect segments__domain__with_conout -g | grep "^1/" | cut -d"|" -f2)
        tmpfile=$(g.tempfile $$)
        for c in $duplines; do
            newcat=$(($(v.category -g segments__domain__with_conout op=report | grep "1 line" | cut -d" " -f5) +1))
            oldcomp=$(v.db.select segments__domain__with_conout col=comp where="cat=$c" -c)
            # This seems to be buggy: v.edit segments__domain__with_conout tool=catadd id=$lid type=line cat=$newcat
            v.out.ascii segments__domain__with_conout form=standard cat=$c | head -n-1 > $tmpfile
            echo "1   $newcat" >> $tmpfile
            v.edit segments__domain__with_conout tool=delete type=line cat=$c
            v.edit segments__domain__with_conout tool=add in=$tmpfile
            db.execute "INSERT INTO $tbl (cat, comp) VALUES ($newcat, $oldcomp);"
        done
        g.rename vect=segments__domain,segments__domain__without_conout
        g.rename vect=segments__domain__with_conout,segments__domain
    else
        create_empty_vector confluence_outlet
    fi

    # add more attributes
    v.what.rast map=segments__domain layer=2 raster=$input_grwl_mask col=grwl_value
    v.what.rast map=segments__domain layer=2 raster=$input_river_width col=river_width
    v.what.rast map=segments__domain layer=2 raster=$input_dem col=elevation

    # drop invalid components
    river_network_utils.py valid_network_components segments__domain segments segments_invalid

    # make segments raster
    v.to.rast input=segments use=cat out=segments memory=$MEMORY

    # for ease of use/inspection
    v.generalize segments method=douglas out=segments_simple thres=$MAX_REACH_LENGTH

    clean_tmp
}


rivgraph_input(){
    input_segments_network=$1
    input_segments_raster=$2
    input_river_halfwidth=$3
    input_elevation=$4
    input_accumulation=$5
    output_mapset=$6  # may exist already
    output_dir=$7

    echo Output: $output_mapset, $output_dir

    g.mapset -c $output_mapset
    g.region -p rast=$input_segments_raster
    g.copy vect=$input_segments_network,segments__
    mkdir -p $output_dir
    # interger C-order cell index, python 0-based, double as it exceeds the integer range
    r.mapcalc exp="idx__ = double(col()-1) * nrows() + (row()-1)"

    # fill differences between GRWL and DEM
    r.mapcalc exp="DEM__nogaps=if(not(isnull($input_segments_raster)) & isnull($input_elevation), 0, $input_elevation)"
    r.mapcalc exp="accumulation__nogaps=if(not(isnull($input_segments_raster)) & isnull($input_accumulation), 1, $input_accumulation)"

    # distance from outlet raster
    v.to.rast segments__ layer=2 where="type='outlet'" use=val output=outlets__ memory=$MEMORY
    r.grow.distance outlets__ distance=outlet__distance
    v.what.rast map=segments__ layer=2 raster=outlet__distance column=outlet_distance

    r.mapcalc exp="mainchannel__accum=if($input_accumulation > $RIVGRAPH_MAINCHANNEL_DAREA, $input_accumulation, null())"
    r.grow.distance mainchannel__accum value=mainchannel__accum__gradient distance=mainchannel__distance
    v.what.rast map=segments__ layer=2 raster=mainchannel__accum__gradient column=mainchannel_darea

    river_network_utils.py rivgraph_nodes_input segments__ idx__ \
        - to_csv $output_dir/nodes.csv
    river_network_utils.py rivgraph_links_input segments__ $input_segments_raster idx__ \
        $input_river_halfwidth DEM__nogaps accumulation__nogaps $output_dir/nodes.csv - to_csv $output_dir/links.csv
    python -c "import grass.script as g; r=g.region(); print(r['cols'], ',', r['rows'])" > $output_dir/cols_rows
    
    clean_tmp true
}


rivgraph_output_to_grass(){
    input_segments=$1
    input_flip_links_csv=$2
    input_node_info=$3
    input_link_info=$4
    output_mapset=$5

    g.mapset -c $output_mapset

    river_network_utils.py flip_segments $input_segments $input_flip_links_csv segments_with_cycles \
        --node_attr=$input_node_info --link_attr=$input_link_info

    river_network_utils.py remove_directed_cycles segments_with_cycles segments 

    # not needed if Great Lakes not part of valid segments
    if [ "{wildcards.domain}" == "MISS__with__Great_Lakes" ]; then
        python -c "from river_network_utils import *; import networkx as nx; \
            stnode=58860; \
            G=river_networkx('segments', directed=False); G.remove_edge(stnode, 58832); nodes=nx.descendants(G, stnode) | {{stnode}}; \
            lines=nx.to_pandas_edgelist(G.subgraph(nodes))['id'].to_list(); \
            open('tmp_MISS_prune_lines', 'w').write(','.join(map(str, lines))); \
            open('tmp_MISS_prune_nodes', 'w').write(','.join(map(str, nodes - {{stnode}})))"
        v.edit segments tool=delete cats=$(cat tmp_MISS_prune_lines)
        v.edit segments la=2 tool=delete cats=$(cat tmp_MISS_prune_nodes)
        r="DELETE FROM segments_segments WHERE cat in ( $(cat tmp_MISS_prune_lines) );"
        db.execute "$r"
        r="DELETE FROM segments_segments_nodes WHERE cat in ( $(cat tmp_MISS_prune_nodes) );"
        db.execute "$r"
    fi

    river_network_utils.py upload_strahler_order segments
    v.extract segments_with_cycles where="cycle>0" layer=1 output=cycles
    v.extract segments layer=2 where="continuity_violated=1" output=discontinuity__la2
    if [ $(v.info discontinuity__la2 -t | grep points) == "points=0" ]; then
        g.rename vect=discontinuity__la2,discontinuities
    else
        # switch cats and table to layer 1
        v.category discontinuity__la2 layer=2,1 option=chlayer output=discontinuities
        v.db.connect discontinuities layer=2 -d
        v.db.connect discontinuities table=discontinuities
    fi
}


centerline_subbasins(){
    input_drainage=$1
    input_centerline_rast=$2
    input_memory=$3
    output_mapset=$4

    g.mapset -c $output_mapset
    g.region -p rast=$input_drainage

    r.to.vect input=$input_centerline_rast type=point output=centerline_points
    r.stream.basins direction=$input_drainage basins=centerline_subbasins points=centerline_points \
        memory=$input_memory
    # area of subbasins
    r.mapcalc exp="cell_area_sqkm=area()/10^(6)"
    r.stats.zonal cover=cell_area_sqkm base=centerline_subbasins out=subbasin_area method=sum
    r.mapcalc exp="one_cell=1"
    r.stats.zonal cover=one_cell base=centerline_subbasins out=subbasin_ncells method=sum
    
    clean_tmp
}

centerline_accumulation(){

    # all inputs should come from the r_watershed_sfd output
    input_drainage=$1
    input_centerlines=$2
    input_accumulation=$3

    g.region -p rast=$input_drainage

    # single route accumulation into centreline
    accum=" \
    1 \
    + if($input_drainage[-1, 1]==1 & isnull($input_centerlines[-1,1]), $input_accumulation[-1,1]) \
    + if($input_drainage[ 0, 1]==2 & isnull($input_centerlines[0, 1]), $input_accumulation[0,1]) \
    + if($input_drainage[ 1, 1]==3 & isnull($input_centerlines[1, 1]), $input_accumulation[1,1]) \
    + if($input_drainage[ 1, 0]==4 & isnull($input_centerlines[1, 0]), $input_accumulation[1,0]) \
    + if($input_drainage[ 1,-1]==5 & isnull($input_centerlines[1,-1]), $input_accumulation[1,-1]) \
    + if($input_drainage[ 0,-1]==6 & isnull($input_centerlines[0,-1]), $input_accumulation[0,-1]) \
    + if($input_drainage[-1,-1]==7 & isnull($input_centerlines[-1,-1]),$input_accumulation[-1,-1]) \
    + if($input_drainage[-1, 0]==8 & isnull($input_centerlines[-1,0]), $input_accumulation[-1,0]) \
    "
    r.mapcalc exp="sfd_centreline_accumulation=if(isnull(centerlines), 0, $accum)" --o
    
    clean_tmp
}

reaches(){
    input_segments_vect=$1
    input_segments_rast=$2
    input_grwl_mask=$3
    input_reach_id_rank=$4
    output_mapset=$5

    global_id_offset=$(($LOCAL_ID_SPACE*$input_reach_id_rank))

    g.mapset -c $output_mapset
    g.region rast=$input_grwl_mask -p

    v.split input=$input_segments_vect output=segments__split length=$MAX_REACH_LENGTH

    # get segment_id into the table
    v.category segments__split op=add layer=3 output=segments__split__reach_ids type=line
    v.db.addtable segments__split__reach_ids layer=3 column="segment_id int" table=reaches_lines
    v.to.db map=segments__split__reach_ids layer=3 option=query query_layer=1 \
        query_column="cat+$global_id_offset" columns=segment_id

    # change layer 3 to layer 1 again (reach lines)
    v.db.droptable segments__split__reach_ids layer=1 -f
    v.category in=segments__split__reach_ids op=del layer=1 out=reaches__no__l1 cat=-1
    v.category in=reaches__no__l1 op=chlayer layer=3,1 out=reaches__l1cats
    v.db.connect -d reaches__l1cats layer=3
    v.db.connect -o reaches__l1cats layer=1 table=reaches_lines
    # add reach nodes
    v.net -c in=reaches__l1cats out=reaches_zigzag operation=nodes

    v.db.droptable reaches_zigzag layer=2 -f
    v.db.addtable reaches_zigzag layer=2 column="segment_node_id int" table=reaches_nodes
    v.db.join reaches_zigzag layer=2 col=cat other_table=segments__split_segments_nodes other_col=cat subset_columns=n_segments
    v.db.update reaches_zigzag layer=2 column=segment_node_id where="n_segments IS NOT NULL" query_col="cat+$global_id_offset"

    # add global ids and add global segment ids
    v.db.addcolumn reaches_zigzag layer=1 col="global_id int"
    v.db.update reaches_zigzag layer=1 col=global_id query_col="cat+$global_id_offset"
    v.db.addcolumn reaches_zigzag layer=2 col="global_id int"
    v.db.update reaches_zigzag layer=2 col=global_id query_col="cat+$global_id_offset"

    # add component
    python -c "from grass_attribute_table import GrassAttributeTable; import pandas as pd
reaches=GrassAttributeTable('reaches_zigzag')
segm = GrassAttributeTable('$input_segments_vect')
reaches['catchment_id'] = segm.loc[reaches.segment_id - $global_id_offset, 'comp'].values + $global_id_offset
reaches.write()
"
    # add global topology ids
    river_network_utils.py add_global_topology_ids reaches_zigzag --global-id-offset=$global_id_offset

    # reach raster, guaranteed to be same as segments raster
    v.to.rast -d in=reaches_zigzag out=reaches__dense use=attr attr=global_id memory=$MEMORY
    r.mapcalc exp="reaches=if(isnull($input_segments_rast), null(), reaches__dense)"

    ## SMOOTH reaches
    # densify for better generalisation
    v.split reaches_zigzag out=reaches__dense -n len=30
    # smooth corners but preserve node locations
    v.generalize reaches__dense method=snake out=reaches__snakes thres=60  # thresh has no effect
    v.generalize reaches__snakes method=reumann out=reaches thresh=4.5
    v.generalize reaches_zigzag method=douglas out=reaches_simple thres=$MAX_REACH_LENGTH

    echo reaches line columns:
    v.db.connect -c reaches
    echo reaches node columns:
    v.db.connect -c reaches

    clean_tmp
}


reach_attributes(){
    input_reaches=$1
    input_grwl_stats=$2
    input_width_stats=$3
    output_mapset=$4

    g.mapset -c $output_mapset

    g.copy vect=$input_reaches,reaches

    # add component
    python -c "from grass_attribute_table import GrassAttributeTable; import pandas as pd
reaches = GrassAttributeTable('reaches')
catglid = reaches.reset_index().set_index('global_id').cat 
grwlinfo = pd.read_csv('$input_grwl_stats', index_col='zone')
grwlinfo.index = catglid[grwlinfo.index].values
width = pd.read_csv('$input_width_stats', index_col='zone')
width.index = catglid[width.index]
reaches['grwl_overlap'] = grwlinfo['non_null_cells']/grwlinfo[['non_null_cells', 'null_cells']].sum(axis=1)
reaches['grwl_value'] = grwlinfo['median'].fillna(0).round(0).astype(int)
reaches['grwl_width_median'] = width['median'].fillna(0)
reaches['grwl_width_std'] = width['stddev'].fillna(0)
reaches.write()
"
    # reach attributes
    v.to.db reaches option=length column=length
    v.to.db reaches option=sinuous column=sinuosity
    v.to.db reaches option=azimuth column=azimuth

    echo reaches line columns:
    v.db.connect -c reaches
    echo reaches node columns:
    v.db.connect -c reaches

    clean_tmp
}


segment_attributes(){
    input_segments=$1
    input_reach_id_rank=$2
    input_csv_list=$3
    output_mapset=$4

    g.mapset -c $output_mapset

    g.copy vect=$input_segments,segments

    for v in length sinuous azimuth ; do
        v.to.db segments option=$v column=$v
    done
    # make tidy columns
    river_network_utils.py tidy_segments_columns segments --global_id_offset=$(($LOCAL_ID_SPACE*$input_reach_id_rank))

    # add darea info
    python -c "from grass_attribute_table import GrassAttributeTable
import pandas as pd; import river_network_utils as ru

lines = GrassAttributeTable('segments', layer=1)
nodes = GrassAttributeTable('segments', layer=2)
csvs = pd.read_csv('$input_csv_list', index_col=0, header=None)[1]

# partitioned d-area
dainfo = pd.read_csv(csvs['darea_segment_info'], index_col='cat')
lines[[f'drainage_area_{io}' for io in dainfo.columns]] = dainfo
dainfo = pd.read_csv(csvs['darea_node_info'], index_col='cat')
nodes['drainage_area'] = dainfo['in']

# mainstem d-area
dainfo = pd.read_csv(csvs['darea_ms_segment_info'], index_col='cat')
lines[[f'drainage_area_mainstem_{io}' for io in dainfo.columns]] = dainfo
dainfo = pd.read_csv(csvs['darea_ms_node_info'], index_col='cat')
nodes['drainage_area_mainstem'] = dainfo['in']

# bifurcation balance
lines['bifurcation_balance_out'] = (lines.drainage_area_out - lines.drainage_area_mainstem_out) / \
        lines[['drainage_area_out', 'drainage_area_mainstem_out']].max(axis=1)

# width share at bifurcation nodes
lwid = lines.set_index('global_id').width_adjusted
dsws = {n: lwid[map(int, li.split(','))] for n, li in nodes.downstream_line_ids.items() if li}
dsws = {n: (w/w.sum()).round(6).astype(str) for n, w in dsws.items()}
nodes['downstream_width_share'] = pd.Series({n: ','.join(w) for n, w in dsws.items()})

# grwl values
grwlinfo = pd.read_csv(csvs['grwl_info'], index_col='zone')
lines['grwl_overlap'] = grwlinfo['non_null_cells']/grwlinfo[['non_null_cells', 'null_cells']].sum(axis=1)
lines['grwl_value'] = grwlinfo['median'].fillna(0).round(0).astype(int)
# OSM names
osmnames = pd.read_csv(csvs['osm_river_names'], index_col='cat')
lines['name'] = osmnames['en'].where(~osmnames['en'].isna(), osmnames['local'])
lines['name'] = ru.fillna_shortest_path('segments', 'name', data=lines, priority='length')
lines['name_local'] = osmnames['local']

nodes.write()
lines.write()
"
    river_network_utils.py upstream_bifurcations segments --output_column=n_bifurcations_upstream

    # make a simplified GRWL only network
    g.copy vect=segments,segments__grwl
    python -c "from grass_attribute_table import GrassAttributeTable
from river_network_utils import river_networkx; import networkx as nx

netw = river_networkx('segments__grwl')
lines = GrassAttributeTable('segments__grwl', layer=1)
nodes = GrassAttributeTable('segments__grwl', layer=2)

lines['purge'] = 1
# filter line criteria
lines.loc[(lines.grwl_value != 180) & (lines.grwl_overlap >= 0.95), 'purge'] = 0
# ignore short dangles
lines.loc[(lines.upstream_line_ids.apply(len)==0) & (lines.length<1000), 'purge'] = 1
# make sure all downstream links are also included
glnodes = lines.loc[lines.purge == 0, 'downstream_node_id'].unique()
ncats = nodes.reset_index().set_index('global_id').loc[glnodes, 'cat']
for n in ncats:
    edges = [netw[f][t][0]['id'] for f, t in nx.dfs_edges(netw, n)]
    lines.loc[edges, 'purge'] = 0
# filter nodes
valid = lines.loc[lines.purge == 0, ['upstream_node_id', 'downstream_node_id']].values.flatten()
nodes['purge']=1
nodes.loc[nodes.global_id.isin(valid), 'purge'] = 0
lines.write()
nodes.write()
"
    river_network_utils.py prune_network_components segments__grwl purge [1]
    v.db.dropcolumn map=segments__grwl layer=1 columns=purge
    v.db.dropcolumn map=segments__grwl layer=2 columns=purge
    v.generalize segments__grwl method=douglas out=segments_simple thres=250
}


subbasins2catchment(){
    input_reclass_rules=$1
    input_subbasins=$2
    input_added_subbasins=$3
    output_name=$4  # rast and vect

    r.reclass in=$input_subbasins out=$output_name rules=$input_reclass_rules

    if [ ! -z $input_added_subbasins ]; then
        r.patch in=$output_name,$input_added_subbasins out=patched__subbasins
        g.rename rast=patched__subbasins,$output_name
    fi
    r.to.vect -v -s -t in=$output_name out=${output_name} type=area
    #v.clean ${output_name}__with__sm__areas tool=rmarea thresh=$((32*30*30)) out=$output_name
    v.db.addtable map=$output_name col="global_id int"
    # copy cat to global_id for later indexing
    v.db.update $output_name column=global_id query_col=cat
    v.to.db $output_name option=area column=area units=kilometers
}

catchments(){
    # Create catchments for network components, segments and reaches by reclassifying the
    # centerline subbasins and converting them to vector
    input_centerline_tbl=$1  # csv tbl with centerline_id,reach_id
    input_centerline_subbasins=$2
    input_reaches_vect=$3  # must have segment_id and comp column
    input_land_mask=$4
    output_mapset=$5

    g.mapset -c $output_mapset
    g.region rast=$input_centerline_subbasins

    temp1=$(g.tempfile $$)
    temp2=$(g.tempfile $$)
    temp3=$(g.tempfile $$)

    # coastal catchments
    r.mapcalc exp="is__coastal__catch=if(not(isnull($input_land_mask)) & isnull($input_centerline_subbasins), 1, null())"
    r.clump in=is__coastal__catch out=coastal__catchments__local
    # check if not a component that is touching domain edges
    r.mapcalc exp="domain__edge=if(row()==1 | col()==1 | row()==nrows() | col()==ncols(), 1, null())"
    r.stats -n -c domain__edge,coastal__catchments__local sep="=" | { cut -d= -f2-3 ; echo "* = NULL"; } \
        | r.reclass in=coastal__catchments__local rules=- out=has__edge
    # offset by max component id, making it global and continues with catchment_ids
    max_comp=$(v.db.select $input_reaches_vect col=catchment_id -c | sort -n | tail -1)
    r.mapcalc exp="coastal_catchments=if(isnull(has__edge), coastal__catchments__local+$max_comp, null())"
    r.to.vect -v -s -t in=coastal_catchments out=coastal_catchments type=area
    v.db.addtable map=coastal_catchments

    # create reclass tables
    python -c "from grass_attribute_table import GrassAttributeTable
import pandas as pd
clp = pd.read_csv('$input_centerline_tbl', index_col='centerline_id')
reaches = GrassAttributeTable('$input_reaches_vect').set_index('global_id')
print(f'Centerline points without reach_id: {clp.reach_id.isna().sum()}')
clp.dropna(subset=['reach_id'], inplace=True)
clp['segment_id'] = reaches.loc[clp.reach_id, 'segment_id'].values
clp['catchment_id'] = reaches.loc[clp.reach_id, 'catchment_id'].values
clp['reach_id'].to_csv('$temp1', sep='=', header=False)
clp['segment_id'].to_csv('$temp2', sep='=', header=False)
clp['catchment_id'].to_csv('$temp3', sep='=', header=False)
"
    subbasins2catchment $temp1 $input_centerline_subbasins "" reach_catchments
    subbasins2catchment $temp2 $input_centerline_subbasins "" segment_catchments
    subbasins2catchment $temp3 $input_centerline_subbasins coastal_catchments component_catchments
    # save reclassed map to raster
    for i in reach segment component; do
        g.rename rast=${i}_catchments,${i}__reclass
        r.mapcalc exp="${i}_catchments=${i}__reclass"
    done

    # assign coastal column
    v.db.addcolumn component_catchments column="is_coastal int DEFAULT 0"
    v.db.update component_catchments value=1 col=is_coastal where="cat > $max_comp"

    clean_tmp
}


height_above_nearest_drainage(){
    input_segments_vect=$1
    input_segments_rast=$2
    input_network_cache=$3
    input_drainage=$4
    input_elevation=$5
    input_memory=$6
    output_mapset=$7

    g.mapset -c $output_mapset
    g.region rast=$input_segments_rast

    river_network_utils.py drainage_from_network \
        $input_segments_vect grit_drainage --network_cache=$input_network_cache

    r.mapcalc exp="drainage=if(isnull(grit_drainage), max($input_drainage, 0), grit_drainage)"

    # activate debug and use self compiled with additional debug to show current outlet for loops
    r.stream.distance stream=$input_segments_rast direction=drainage \
        elevation=$input_elevation dist=distance__all diff=hand__all mem=$input_memory

    r.mapcalc exp="hand=if(hand__all<=100, float(hand__all), null())"
    r.mapcalc exp="distance=if(isnull(distance__all), null(), int(distance__all))"

    clean_tmp
    # check loops: r.mapcalc exp="loops=if((d==1 & d[-1,1]==5) | (d==2 & d[-1,0]==6) | (d==3 & d[-1,-1]==7) | (d==4 & d[0,-1]==8) | (d==5 & d[1,-1]==1) | (d==6 & d[1,0]==2) | (d==7 & d[1,1]==3) | (d==8 & d[0,1]==4), d, null())" --o
    # r.stats -1ng loops
}


river_transects(){
    input_reaches=$1
    input_river_width_column=$2
    output_mapset=$3

    g.mapset -c $output_mapset

    bins="30 60 100 200 400 600 1000 2000 200000"
    floor="$input_river_width_column>0"
    for b in $bins ; do
        v.extract $input_reaches output=rivers__transects__width__$b where="$floor AND $input_river_width_column<=$b"
        bm=$(($b > 10000 ? 10000 : $b))
        v.transects input=rivers__transects__width__$b output=transects___width__$b \
            metric=straight transect_perpendicular=trend transect_spacing=100 dleft=$(($bm * 4)) dright=$(($bm * 4))
        # floor for next iteration
        floor="$input_river_width_column>$b"
        # add table
        v.db.addtable transects___width__$b columns="reach_id int"
        v.db.addtable transects___width__$b layer=2
        v.to.db transects___width__$b query_layer=2 query_column=cat option=query column=reach_id
    done
    v.patch -e input=$(g.list vect pat="transects___width__*" sep=comma) output=transects
}


snap_stations(){
    input_stations_vect=$1
    input_centerline_tbl=$2
    output_mapset=$3

    g.mapset -c $output_mapset

    # in case no stations found in domain, create empty output vectors
    if [ -z "$(v.db.select $input_stations_vect -c)" ]; then
        create_empty_vector stations_valid
        create_empty_vector stations_snapped
        create_empty_vector stations_invalid
        exit 0
    fi
    tmpf1=$(g.tempfile $$).gpkg
    tmpf2=$(g.tempfile $$).gpkg
    v.out.ogr in=$input_stations_vect out=$tmpf1 --o

    snap_stations.py $tmpf1 $input_centerline_tbl --max-distance=$STATION_MAX_SNAP_DISTANCE \
        --max_darea_deviation=$STATION_DAREA_MAX_DEVIATION --output=$tmpf2 - describe - to_string
    v.in.ogr in=$tmpf2 key="cat" out=stations_snapped

    filter="abs(deviation_snapped)>$STATION_DAREA_MAX_DEVIATION OR distance_snapped>$STATION_MAX_SNAP_DISTANCE"
    v.db.droprow in=stations_snapped out=stations_valid where="$filter"
    v.extract in=stations_snapped out=stations_invalid where="$filter"
}
snap_stations_id(){
    input_stations_vect=$1
    input_centerline_tbl=$2
    input_reaches_vect=$3
    output_mapset=$4

    g.mapset -c $output_mapset

    tmpf2=$(g.tempfile $$).gpkg
    v.out.ogr in=$input_reaches_vect out=$tmpf2 --o
    echo reach vector saved
    # in case no stations found in domain, create empty output vectors
    if [ -z "$(v.db.select $input_stations_vect -c)" ]; then
        create_empty_vector stations_valid
        create_empty_vector stations_snapped
        create_empty_vector stations_invalid
        exit 0
    fi
    tmpf1=$(g.tempfile $$).gpkg
    tmpf2=$(g.tempfile $$).gpkg
    tmpf3=$(g.tempfile $$).gpkg
    v.out.ogr in=$input_stations_vect out=$tmpf1 --o
    v.out.ogr in=$input_reaches_vect out=$tmpf2 --o

    snap_stations_p2p.py $tmpf1 $input_centerline_tbl $tmpf2 --max-distance=$STATION_MAX_SNAP_DISTANCE \
        --max_darea_deviation=$STATION_DAREA_MAX_DEVIATION --output=$tmpf3 - describe - to_string
    v.in.ogr in=$tmpf3 key="cat" out=stations_snapped

    filter="abs(deviation_snapped)>$STATION_DAREA_MAX_DEVIATION OR distance_snapped>$STATION_MAX_SNAP_DISTANCE"
    v.db.droprow in=stations_snapped out=stations_valid where="$filter"
    v.extract in=stations_snapped out=stations_invalid where="$filter"
}

match_stations_bkfw(){
    input_stations_vect=$1
    input_centerline_tbl=$2
    input_target_cols=$3
    input_csv_index_cols=$4
    output_mapset=$5

    g.mapset -c $output_mapset

    # in case no stations found in domain, create empty output vectors
    if [ -z "$(v.db.select $input_stations_vect -c)" ]; then
        create_empty_vector stations_valid
        create_empty_vector stations_snapped
        create_empty_vector stations_invalid
        exit 0
    fi
    tmpf1=$(g.tempfile $$).gpkg
    tmpf2=$(g.tempfile $$).gpkg
    v.out.ogr in=$input_stations_vect out=$tmpf1 --o

    snap_stations_p2p.py $tmpf1 $input_centerline_tbl --target_cols=$input_target_cols --csv_index_cols=$input_csv_index_cols \
        --max-distance=$STATION_MAX_SNAP_DISTANCE \
        --max_darea_deviation=$STATION_DAREA_MAX_DEVIATION --output=$tmpf2 - describe - to_string
    v.in.ogr in=$tmpf2 key="cat" out=stations_snapped

    filter="abs(deviation_snapped)>$STATION_DAREA_MAX_DEVIATION OR distance_snapped>$STATION_MAX_SNAP_DISTANCE"
    v.db.droprow in=stations_snapped out=stations_valid where="$filter"
    v.extract in=stations_snapped out=stations_invalid where="$filter"
}

snap_stations_points_binned(){
    input_stations_vect=$1
    input_river_lines=$2
    input_darea_rast=$3
    input_ref_darea_column=$4
    output_mapset=$5

    g.mapset -c $output_mapset
    g.region rast=$input_darea_rast

    # in case no stations found in domain, create empty output vectors
    if [ -z "$(v.db.select $input_stations_vect -c)" ]; then
        create_empty_vector stations_valid
        create_empty_vector stations_snapped
        create_empty_vector stations_invalid
        exit 0
    fi

    # loop over drainage area classes to snap stations only to rivers in same class
    patch=""
    sep=""
    lb=50
    for ub in 100 200 500 1000 2000 5000 10000 50000 100000 1000000 10000000; do
        v.extract in=$input_stations_vect out=stations__$ub \
            where="$input_ref_darea_column > $lb AND $input_ref_darea_column <= $ub"
        if [ ! -z "$(v.db.connect -g stations__$ub)" ]; then
        if [ ! -z "$(v.db.select stations__$ub -c)" ]; then
            v.extract in=$input_river_lines out=streams__$ub \
                where="drainage_area > $lb*0.8 AND drainage_area <= $ub*1.2"
            snap_points stations__$ub streams__$ub stations_snapped__$ub
            patch=$patch${sep}stations_snapped__${ub}
            sep=,
        fi
        fi
        lb=$ub
    done
    # either patch multiple classes or rename if only single class with stations
    if (( $(echo $patch | grep -o , | wc -l) > 0 )); then
        v.patch -e in=$patch out=stations_snapped
    else
        g.rename vect=$patch,stations_snapped
    fi

    # add darea and deviation
    v.what.rast stations_snapped raster=$input_darea_rast column=darea_derived
    v.db.addcolumn stations_snapped column='darea_deviation double'
    v.db.update stations_snapped column=darea_deviation query_col="(darea_derived/$input_ref_darea_column-1)*100"

    # filter stations
    filter="abs(darea_deviation)>$STATION_DAREA_MAX_DEVIATION OR snapped_distance>$STATION_MAX_SNAP_DISTANCE"
    v.db.droprow in=stations_snapped out=stations_valid where="$filter"
    v.extract in=stations_snapped out=stations_invalid where="$filter"

    clean_tmp
}


snap_stations_lines(){
    input_stations_vect=$1
    input_river_lines=$2
    input_darea_rast=$3
    input_ref_darea_column=$4
    output_mapset=$5

    g.mapset -c $output_mapset
    g.region rast=$input_darea_rast

    if [ -z "$(v.db.select $input_stations_vect -c)" ]; then
        create_empty_vector stations_valid
        create_empty_vector stations_snapped
        create_empty_vector stations_invalid
        exit 0
    fi

    snap_points $input_stations_vect $input_river_lines stations_snapped

    v.what.rast stations_snapped raster=$input_darea_rast column=darea_derived

    v.db.addcolumn stations_snapped column='darea_deviation double'
    v.db.update stations_snapped column=darea_deviation query_col="(darea_derived/$input_ref_darea_column-1)*100"

    # filter stations
    filter="abs(darea_deviation)>$STATION_DAREA_MAX_DEVIATION OR snapped_distance>$STATION_MAX_SNAP_DISTANCE"
    v.db.droprow in=stations_snapped out=stations_valid where="$filter"
    v.extract in=stations_snapped out=stations_invalid where="$filter"

    # try to snap invalid stations to second nearest segment
    # patch=""
    # for s in $(v.db.select -c stations_invalid col=cat,segment_id,area); do
    #     sid=$(echo $s | cut -d"|" -f1)
    #     v.extract stations_invalid cat=$sid out=station__invalid
    #     v.extract $input_river_lines cats=$(echo $s | cut -d"|" -f2) -r out=segments__bar1 --o
    #     snap_points station__invalid segments__bar1 station__snapped__again__$sid
    #     patch=$patch,station__snapped__again__$sid
    # done
    # v.patch -e in=$patch out=stations_invalid_resnapped
    # v.what.rast stations_invalid_resnapped raster=$input_darea_rast column=darea_derived_2nd

    clean_tmp
}

flow_accumulation(){
    input_dem_accumulation=$1
    input_network_darea=$2
    input_cpus=$3
    output_mapset=$4

    g.mapset -c $output_mapset
    g.region rast=$input_dem_accumulation

    # outside of network accumulation must not exceed minimal darea (patches inconsistencies between network and accumulation routing)
    r.mapcalc exp="drainage_area=if(isnull($input_network_darea), \
        if($input_dem_accumulation > $MIN_STREAM_DAREA, area()/1000000, $input_dem_accumulation), $input_network_darea)"
    r.colors -n -a map=drainage_area color=inferno

    # check where accumulation was reduced
    #r.mapcalc exp="spurious__accumulation=if(isnull($input_network_darea) & $input_dem_accumulation > $MIN_STREAM_DAREA, $input_dem_accumulation, null())"
    #r.to.vect in=spurious__accumulation type=point out=spurious_accumulation column=darea

    #r.patch in=$input_network_darea,$input_dem_accumulation nprocs=$input_cpus memory=$(($input_cpus*$MEMORY)) \
    #    output=drainage_area
    
    r.mapcalc exp="darea_diff=if(not(isnull($input_network_darea)), $input_network_darea - $input_dem_accumulation, null())"

    clean_tmp
}


stations_osm_riverbank(){
    input_stations=$1
    input_riverbank_polygons=$2
    input_reaches=$3
    output_mapset=$4
    # add GRASS extension
    # g.extension v.transects url="https://github.com/OSGeo/grass-addons/tree/grass8/src/vector/v.transects"
    echo printing available grass extension
    g.extension -a
    g.mapset -c $output_mapset

    TRANSECT_DISTANCE=50  # m
    TRANSECT_LENGTH=2000  # m from centerline
    if [ -z "$(v.db.select $input_stations -c)" ]; then
        create_empty_vector stations_osm_riverbank
        exit 0
    fi

    # select stations within water polygons
    v.select ain=$input_stations bin=$input_riverbank_polygons out=stations_valid_osm_riverbank
    # YL
    if [ -z "$(g.findfile vector file=stations_valid_osm_riverbank -n | head -n1 | cut -d"=" -f2)" ]; then
        create_empty_vector stations_osm_riverbank
        exit 0
    fi

    # get osm width into table
    tmpf=$(v.db.select stations_valid_osm_riverbank column=reach_id -c | paste -s -d, -)
    v.extract $input_reaches out=station_osm_reaches where="global_id in ($tmpf)"
    # create transects with reference to reach id
    v.transects in=station_osm_reaches out=reach_transects transect_spacing=$TRANSECT_DISTANCE \
        dleft=$TRANSECT_LENGTH dright=$TRANSECT_LENGTH
    # YL 
    # problem -s
    echo printing the info of output!
    v.info reach_transects
    echo add column reach_cat to table!
    v.db.addtable reach_transects --verbose columns="reach_cat int" key=cat layer=1
    v.db.addtable reach_transects --verbose layer=2 key=cat
    echo printing info of reach_transects layer2!
    v.info -c reach_transects layer=2
    v.to.db reach_transects query_layer=2 query_column=cat option=query column=reach_cat


    # reach_id value missing, the query does not copy the cat to reach_cat
    # problem -e
    v.db.join reach_transects col=reach_cat other_table=station_osm_reaches other_column=cat subset_columns=global_id
    v.db.renamecolumn reach_transects col=global_id,reach_id
 


    # overlay with polygons
    v.overlay ain=reach_transects bin=$input_riverbank_polygons op=and output=reach_transects_in_bank olayer=1,2,0 -t
    # YL 
    echo !
    v.db.connect -p reach_transects_in_bank
    v.db.addtable reach_transects_in_bank --verbose columns="transect_cat"
    v.db.addtable reach_transects_in_bank --verbose layer=2

    v.to.db reach_transects_in_bank query_layer=2 query_column=cat option=query column=transect_cat
    # YL
    echo printing reach_transects_in_bank after to.db!
    v.db.select reach_transects_in_bank
    v.select ain=reach_transects_in_bank bin=station_osm_reaches out=reach_transects_in_bank_reach
                # YL
    echo printing reach_transects_in_bank_reach after v.select!
    v.db.select reach_transects_in_bank_reach
    
    v.db.join reach_transects_in_bank_reach col=transect_cat other_table=reach_transects other_column=cat subset_columns=reach_id
    echo join done!
    # print-YL
    echo reach_transects_in_bank_reach,layer1!
    v.db.select reach_transects_in_bank_reach layer=1
    echo reach_transects_in_bank_reach,layer2!
    v.db.select reach_transects_in_bank_reach layer=2

    v.to.db reach_transects_in_bank_reach op=length col=width
    # print-YL
    echo reach_transects_in_bank_reach connections!
    v.db.connect -p reach_transects_in_bank_reach
    echo reach_transects_in_bank_reach,layer1 after todb!
    v.db.select reach_transects_in_bank_reach layer=1
    echo reach_transects_in_bank_reach,layer2 after todb!
    v.db.select reach_transects_in_bank_reach layer=2

    python -c "from grass_attribute_table import GrassAttributeTable
ts = GrassAttributeTable('reach_transects_in_bank_reach')
st = GrassAttributeTable('stations_valid_osm_riverbank')
ts['reach_id'] = ts.reach_id.astype(int)
gb = ts.groupby('reach_id').width.describe()
for i in list(set(st.reach_id) - set(gb.index)):
    gb.loc[i] = 0
st['osm_median_width'] = gb.loc[st.reach_id, '50%'].values
st['osm_n_transects'] = gb.loc[st.reach_id, 'count'].astype(int).values
st.write()
"

    # only keep stations that have more than 1 valid transect
    v.db.droprow stations_valid_osm_riverbank out=stations_osm_riverbank where="osm_n_transects<2"
}


raster_output_grid(){
    input_domain=$1
    output_mapset=$2
    output_tile_list=$3

    g.mapset $output_mapset

    python -c "import grass.script as grass; from grass_attribute_table import GrassAttributeTable
reg, res = grass.region(), $RASTER_OUTPUT_GRID_SIZE
llcorner = [int(reg['w']/res)*res, int(reg['s']/res)*res]
size = [int((reg['n'] - llcorner[1])/res) + 1, int((reg['e'] - llcorner[0])/res) + 1]
grass.run_command('v.mkgrid', map='grid__all', grid=size, box=[res, res], coor=llcorner, position='coor')
tbl = GrassAttributeTable('grid__all')
tbl['west'] = (tbl.col - 1)*res + llcorner[0]
tbl['west'] = tbl.west.apply(lambda t: 'W%09i' % (t * -1) if t < 0 else 'E%09i' % t)
tbl['south'] = (size[0] - tbl.row)*res + llcorner[1]
tbl['south'] = tbl.south.apply(lambda t: 'S%09i' % (t * -1) if t < 0 else 'N%09i' % t)
tbl.write()
"
    v.select ain=grid__all bin=$input_domain output=grid

    v.db.select grid col=cat,west,south sep=comma > $output_tile_list
}

regional_stream_density(){
    # This doesnt seem to work as correlation between width and drainage are is poor
    # output from r_watershed_sfd
    input_river_width=$1
    input_centerlines=$2
    input_accumulation_sqkm=$3
    input_subbasins=$4
    # from segments calculation, bar-filled river mask with GRWL values to filter rivers from other classes   
    input_grwl=$5

    g.region -p rast=$input_accumulation_sqkm
    r.mask $input_accumulation_sqkm

    # create subregions and cleaned subbasins
    # only use centerline cells of rivers narrower than MAX_WIDTH_THRESH and are rivers (i.e. not lakes or tidal or anything else)
    r.mapcalc exp="centreline__width__thresh=if(not(isnull($input_centerlines)) \
        & $input_river_width < $MAX_WIDTH_THRESHOLD & $input_grwl == 255, $input_river_width, null())"
    r.stats.zonal base=$input_subbasins cover=centreline__width__thresh method=count output=subbasin__centerline__count
    r.mapcalc exp="subbasins__large=if(subbasin__centerline__count>=$MIN_CENTRELINE_SUBBASIN_CELLS, $input_subbasins, null())"
    r.grow.distance in=subbasins__large value=subbasins__large__grown \
        maximum_dist=$(python -c "import grass.script as g; print(g.region()['nsres'] * 1.1)")
    r.mapcalc exp="subbasins__large__grown=int(subbasins__large__grown)" --o
    r.mapcalc exp="subbasins__small=if(isnull(subbasins__large), 1, null())"
    r.clump input=subbasins__small out=subbasins__small__clumped
    r.statistics base=subbasins__small__clumped cover=subbasins__large__grown method=mode output=subbasins__small__largeids  # largeids are in labels!
    r.mapcalc exp="subregions=if(isnull(subbasins__large), int(@subbasins__small__largeids), int(subbasins__large))"
    r.colors map=subregions color=random
    
    clean_tmp
}

# execute functions from commandline: script.sh function_name arguments ...
"$@"
