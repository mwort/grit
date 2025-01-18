#!/usr/bin/env bash
set -e -u
export GRASS_OVERWRITE=1


cut_vector(){
    input_region=$1
    input_vector=$2
    output=$3

    echo $input_region

    g.region $input_region -p

    v.in.region out=ROI__
    v.select ain=$input_vector bin=ROI__ alayer=1 atype=line out=cut__vect__ROI
    layers=$(v.db.connect -g $input_vector | cut -d"|" -f1 | sort | cut -d'/' -f2)
    echo $(v.db.connect -g $input_vector)
    echo Found layers: $layers
    # bug in v.select means attributes arent carried over/table is empty
    v.db.droptable cut__vect__ROI la=1 -f
    v.db.addtable cut__vect__ROI
    g.copy vect=$input_vector,tmp__
    tbl=$(v.db.connect -g tmp__ | cut -d"|" -f2 | sort | head -1)
    echo Found table: $tbl
    v.db.join cut__vect__ROI layer=1 col=cat other_table=$tbl other_col=cat

    v.out.ogr in=cut__vect__ROI out=$output layer=1 output_layer=$(echo $layers | cut -d" " -f1)

    if (( $(echo $layers | wc -w) == 2 )); then
            v.select ain=$input_vector bin=ROI__ alayer=2 out=cut__vect__ROI_2 --o
            v.out.ogr -u in=cut__vect__ROI_2 out=$output layer=2 output_layer=$(echo $layers | cut -d" " -f2)
    fi
    # clean up
    g.remove -f vect pat=*__*
    g.region -d
}

cut_raster(){
    input_region=$1  # needs to have the correct 
    input_raster=$2
    output=$3

    g.region $input_region align=$input_raster -p

    r.out.gdal in=$input_raster out=$output createopt="COMPRESS=DEFLATE" -m -c

}

# enable function execution as subcommands of script
$@