#!/usr/bin/env bash

input_segments=$1
input_cats=$2
output=$3

outvect=rg__fixlink__fixlinks
v.extract $input_segments cat=$input_cats out=$outvect
v.edit map=rg__fixlink__fixlinks tool=flip cats=$input_cats

if [ -f $output ]; then
    echo Appending to $output
    v.in.ogr $output out=rg__fixlink__output
    v.db.droptable rg__fixlink__output -f
    v.patch --o -a in=rg__fixlink__fixlinks out=rg__fixlink__output
    v.db.addtable rg__fixlink__output
    outvect=rg__fixlink__output
fi

v.db.select $outvect
v.out.ogr --o $outvect out=$output

g.remove vect pat=rg__fixlink__* -f