#!/usr/bin/env bash
set -e -u -x
export PS4='$(date --iso-8601=seconds): $(printf "%5d" "$LINENO") >> '
export GRASS_OVERWRITE=1

MAX_DISCRETISATION_LENGTH=400  # m


split_lines(){
    input_lines=$1
    input_length=$2
    output_lines=$3

    v.split in=$input_lines out=${input_lines}__split len=$input_length

    v.category ${input_lines}__split op=add layer=3 output=${input_lines}__split__ids type=line
    v.db.addtable ${input_lines}__split__ids layer=3 column="line_id int" table=${output_lines}_lines
    v.to.db map=${input_lines}__split__ids layer=3 option=query query_layer=1 \
        query_column="cat" columns=line_id

    # change layer 3 to layer 1 again (reach lines)
    v.db.droptable ${input_lines}__split__ids layer=1 -f
    v.category in=${input_lines}__split__ids op=del layer=1 out=${input_lines}__split__ids__no__l1 cat=-1
    v.category in=${input_lines}__split__ids__no__l1 op=chlayer layer=3,1 out=$output_lines
    v.db.connect -d $output_lines layer=3
    v.db.connect -o $output_lines layer=1 table=${output_lines}_lines
    tbl=$(v.db.connect -g $input_lines | grep "^1/" | cut -d"|" -f2)
    cols=$(v.db.connect -c $input_lines | cut -d"|" -f2 | tail -n +2 | paste -sd, -)
    echo $tbl
    v.db.join $output_lines col=line_id other_table=$tbl other_col=cat sub=$cols
}


select_features(){
    input_select_from=$1
    input_select_with=$2
    output_selected=$3

    v.select ain=$input_select_from bin=$input_select_with out=$output_selected op=within
    v.db.droptable -f $output_selected
    v.db.addtable $output_selected
    tbl=$(v.db.connect -g $input_select_from | grep "^1/" | cut -d"|" -f2)
    v.db.join $output_selected col=cat other_table=$tbl other_col=cat
}

didnt_grit_snap(){
    input_didnt=$1
    input_reaches=$2
    input_catchment=$3
    output_mapset=$4
    output_gpkg=$5

    g.mapset -c $output_mapset

    g.copy vect=$input_reaches,grit__reaches,$input_didnt,didnt_all

    # prepare didnt network
    v.extract $input_catchment where="is_coastal=0" out=non__coastal__catchments
    split_lines didnt_all $MAX_DISCRETISATION_LENGTH didnt_split
    select_features didnt_split non__coastal__catchments didnt_selected
    v.to.db didnt_selected option=azimuth col=azimuth
    v.to.db didnt_selected option=length col=length_split

    # subset reaches to avoid excessive vect processing
    g.region -a -p vect=didnt_selected grow=100
    v.in.region out=ROI

    # prepare GRIT
    select_features grit__reaches ROI grit_reaches_selected
    split_lines grit_reaches_selected $MAX_DISCRETISATION_LENGTH grit_reaches_split
    v.to.db grit_reaches_split option=azimuth col=azimuth

    # snapping
    v.db.addcolumn didnt_selected col="grit_cat int,snap_dist double,grit_azimuth double"
    v.distance from=didnt_selected from_type=line to=grit_reaches_split to_type=line \
        upload=cat,dist,to_attr col="grit_cat,snap_dist,grit_azimuth" to_column=azimuth
    # calculate minimum azimuth difference
    python -c "from grass_attribute_table import GrassAttributeTable; \
        didnt = GrassAttributeTable('didnt_selected'); \
        didnt['azimuth_diff'] = ((didnt.azimuth - didnt.grit_azimuth + 180) % 360 - 180).abs(); \
        didnt.write(); vald=didnt.snap_dist<didnt.avgwidth ; \
        print(f'Found {(~vald).sum()*100/len(didnt)}% invalid, {(vald & (didnt.azimuth_diff<90)).sum()*100/vald.sum()}% with right direction.')"

    v.out.ogr didnt_selected output=$output_gpkg

    # clean up
    g.remove -f vect pat="*__*"
}

nhd_grit_snap(){
    input_nhd=$1
    input_reaches=$2
    input_catchment=$3
    output_mapset=$4
    output_gpkg=$5

    g.mapset -c $output_mapset

    g.copy vect=$input_reaches,grit__reaches,$input_nhd,nhd__all

    # prepare NHD network
    v.extract $input_catchment where="is_coastal=0" out=non__coastal__catchments
    v.dissolve in=non__coastal__catchments col=is_coastal out=valid__catchment__area
    # remove all columns except those needed
    keep="cat|Divergence|TotDASqKM|DivDASqKM"
    todel=$(v.db.connect -c nhd__all | cut -d"|" -f2 | grep -v -E "$keep" | paste -sd, -)
    v.db.dropcolumn map=nhd__all col=$todel
    #select_features nhd__all valid__catchment__area nhd_selected
    v.overlay ain=nhd__all bin=valid__catchment__area out=nhd_selected op=and
    split_lines nhd_selected $MAX_DISCRETISATION_LENGTH nhd_split
    v.to.db nhd_split option=azimuth col=azimuth
    v.to.db nhd_split option=length col=length_split

    # subset reaches to avoid excessive vect processing
    g.region -a -p vect=nhd_split grow=100
    v.in.region out=ROI

    # prepare GRIT
    select_features grit__reaches ROI grit_reaches_selected
    split_lines grit_reaches_selected $MAX_DISCRETISATION_LENGTH grit_reaches_split
    v.to.db grit_reaches_split option=azimuth col=azimuth

    # snapping
    v.db.addcolumn nhd_split col="grit_cat int,snap_dist double,grit_azimuth double"
    v.distance from=nhd_split from_type=line to=grit_reaches_split to_type=line \
        upload=cat,dist,to_attr col="grit_cat,snap_dist,grit_azimuth" to_column=azimuth
    # calculate minimum azimuth difference
    python -c "from grass_attribute_table import GrassAttributeTable; \
        didnt = GrassAttributeTable('nhd_split'); \
        didnt['azimuth_diff'] = ((didnt.azimuth - didnt.grit_azimuth + 180) % 360 - 180).abs(); \
        didnt.write()"

    v.out.ogr nhd_split output=$output_gpkg

    # clean up
    g.remove -f vect pat="*__*"
}

"$@"