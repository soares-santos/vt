#!/bin/bash

#set -x 

if [ "$#" -ne 1 ] ; then echo "Usage: $0 <galaxy_cat>" ; exit ; fi   

echo "VT :: merge run1 and run2 catalogs"
echo `date`

ID="ID"
if [ -f params.info ] 
then source params.info
else echo "Warning: File params.info not found. Using defaults."
fi

clustercat=${1}.clustercat.fit
memberscat=${1}.memberscat.fit
rankedlist=${1}.clusters.fit

echo "ranking clusters by Ngals"

# rank clusters by nvt
funtable -s "nvt sig" ${clustercat}[1] tmpC.fit
fitscopy tmpC.fit'[col *;RANK=(NAXIS2+1-#row)]' \!${clustercat}
funtable -s "RANK" ${clustercat}[1] tmpC.fit && mv tmpC.fit ${clustercat}

# add rank column to members catalog
fitscopy ${clustercat}'[col RANK;id]' tmpC.fit
funindex tmpC.fit[1] id
funindex ${memberscat}[1] host_id
funjoin -m 1 -j1 host_id -j2 id ${memberscat}[1] tmpC.fit[1] tmpM.fit
fitscopy tmpM.fit'[col -id_2]' \!${memberscat} && rm tmpM.fit
rm -f tmpC_id.idx ${memberscat}_host_id.idx tmpC.fit

# create input file for matchscript.py
fitscopy ${clustercat}'[col RANK;ra_c;dec_c;redsh=z]' tmpC.fit
fitscopy ${memberscat}'[col gid=id;RANK]' tmpM.fit
fitscopy tmpC.fit'[col RANK;RA=ra_c;DEC=dec_c;Z=redsh]' tmp.fit && mv tmp.fit tmpC.fit
fitscopy tmpM.fit'[col ID=gid;RANK]' tmp.fit && mv tmp.fit tmpM.fit
funtable tmpC.fit[1] $rankedlist "RANK RA DEC Z" && rm tmpC.fit
funtable -a tmpM.fit[1] $rankedlist "ID RANK" && rm tmpM.fit
rm -f *.idx

echo "cluster catalog done."
echo `date`

exit 

