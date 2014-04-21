#!/bin/bash

if [ "$#" -ne 3 ] ; then 
    echo "usage: get_boundaries.sh file col1"
    exit 
fi

nrows=`liststruc $1 | grep rows | awk '{print $8}'`

funtable -s $2 $1[1] tmp_sorted.$1
filter='[col '$2'][#row==1||#row=='$nrows']'
dumpfits tmp_sorted.$1"$filter"

#funtable -s $3 $1[1] tmp_sorted.$1
#filter='[col '$3'][#row==1||#row=='$nrows']' 
#dumpfits tmp_sorted.$1"$filter"

rm tmp_sorted.$1

exit

