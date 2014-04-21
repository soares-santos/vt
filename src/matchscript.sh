#!/bin/bash

#set -x

if [ "$#" -ne 2 ] 
then  
    echo "Usage: $0 <file1> <file2>"
    echo " Files are fits tables with at least two columns named  id  and  host_id".
    echo " First file has the host_id repeated in all rows".
    echo " The program will return one line "
    echo "    ID_1 ID_2 fmatch"
    echo " where ID_1 is the host_id in the first file "
    echo "       ID_2 is the host_id of the matched cluster in the second file "
    echo "       fmatch is the fraction of rows of the first file that were matched " 
    exit
fi

ID_A=`dumpfits $1'[#row==1][col host_id]' | awk '{print $1}'`
NA=`liststruc $1 | grep rows | awk '{print $8}' `

bla=`funindex $1[1] id` 
bla=`funindex $2[1] id`
funjoin -j id -m 1 -a1 "id host_id" -a2 "host_id" $1[1] $2[1] tmp.fit 

NM=`liststruc tmp.fit | grep rows | awk '{print $8}' `

if [ $NM -eq 0 ] ; 
then
    echo "$ID_A -1 0.00"
    exit
fi

boundaries=`get_boundaries.sh $2 host_id id` 
binning=`echo $boundaries | awk '{print $1":"$2":1"}'` 
tmp=`funhist -w tmp.fit[1] host_id_2 $binning | grep -v " 0 " | tail +8 | awk '{print $2,int($4)}' | sort -nr | head -1`
ID_B=`echo $tmp | awk '{print $2}'`
NM=`echo $tmp | awk '{print $1}'`

output=`echo "$NA $NM $ID_A $ID_B" | awk '{print $3, $4, ($2 *1.0)/$1}'` 
echo $output
 
exit

