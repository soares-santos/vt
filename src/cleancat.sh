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

candidates1=${1}.candidates.fit 
membercand1=${1}.membercand.fit 
candidates2=${1}.candidates.run2.fit 
membercand2=${1}.membercand.run2.fit 
match_table=${1}.candidates_match_table.fit
clustercat=${1}.clustercat.fit
memberscat=${1}.memberscat.fit
rankedlist=${1}.clusters.fit

if [ ! -f $candidates1 ] ; then echo "File $candidates1 not found." ; exit ; fi  
if [ ! -f $membercand1 ] ; then echo "File $membercand1 not found." ; exit ; fi  
if [ ! -f $candidates2 ] ; then echo "File $candidates2 not found." ; exit ; fi  
if [ ! -f $membercand2 ] ; then echo "File $membercand2 not found." ; exit ; fi  
if [ ! -f $match_table ] ; then echo "File $match_table not found." ; exit ; fi  

echo "== select primary found in run1"

# select primary found in run1
rm -f tmp.fit filterC.txt filterM.txt
run1=`funtable ${match_table}[1]["primary==1&&run==1"] stdout | fundisp -n stdin[1] "id_cand"`
if [ ! "$run1" = "" ] 
then
filterC=`echo $run1 | awk '{ORS=".or.";i=0;while(i++<NF-1) print "id=="$i      ;ORS="\n"; print "id=="$i}'`
filterM=`echo $run1 | awk '{ORS=".or.";i=0;while(i++<NF-1) print "host_id=="$i ;ORS="\n"; print "host_id=="$i}'`
echo $filterC > filterC.txt
echo $filterM > filterM.txt
fitscopy $candidates1'[@filterC.txt]' tmp.fit && mv tmp.fit $clustercat && rm filterC.txt
fitscopy $membercand1'[@filterM.txt]' tmp.fit && mv tmp.fit $memberscat && rm filterM.txt
fitscopy $clustercat'[col *;dist2boxcenter(1D)=-999.;id_cand(K)=id]' tmp.fit
fitscopy tmp.fit'[col -id]' \!$clustercat
fitscopy $clustercat'[col *;id(K)=id_cand*10+1]' \!tmp.fit && mv tmp.fit $clustercat
fitscopy $memberscat'[col *;host_cand(K)=host_id]' \!tmp.fit
fitscopy tmp.fit'[col -host_id]' \!$memberscat
fitscopy $memberscat'[col *;host_id(K)=host_cand*10+1]' \!tmp.fit && mv tmp.fit $memberscat 
fi

echo "== select primary found in run2"

# select primary found in run2
rm -f tmpC.fit tmpM.fit tmp.fit filterC.txt filterM.txt 
run2=`funtable ${match_table}[1]["primary==1&&run==2"] stdout | fundisp -n stdin[1] "id_cand"`
if [ ! "$run2" = "" ] 
then
filterC=`echo $run2 | awk '{ORS=".or.";i=0;while(i++<NF-1) print "id=="$i      ;ORS="\n"; print "id=="$i}'`
filterM=`echo $run2 | awk '{ORS=".or.";i=0;while(i++<NF-1) print "host_id=="$i ;ORS="\n"; print "host_id=="$i}'`
echo $filterC > filterC.txt
echo $filterM > filterM.txt
fitscopy $candidates2'[@filterC.txt]' tmpC.fit && rm filterC.txt
fitscopy $membercand2'[@filterM.txt]' tmpM.fit && rm filterM.txt
fitscopy tmpC.fit'[col *;id_cand(K)=id]' tmp.fit
fitscopy tmp.fit'[col -id]' \!tmpC.fit
fitscopy tmpC.fit'[col *;id(K)=id_cand*10+2]' \!tmp.fit && mv tmp.fit tmpC.fit
fitscopy tmpM.fit'[col *;host_cand(K)=host_id]' \!tmp.fit
fitscopy tmp.fit'[col -host_id]' \!tmpM.fit
fitscopy tmpM.fit'[col *;host_id(K)=host_cand*10+2]' \!tmp.fit && mv tmp.fit tmpM.fit
fi

echo "== merge run1 and run2"

# merge run1 and run2
tabmerge tmpC.fit+1 ${clustercat}+1 && rm tmpC.fit
tabmerge tmpM.fit+1 ${memberscat}+1 && rm tmpM.fit
 
# join members list with the input catalog on ID
fitscopy ${memberscat}'[col id;local_dens;x_0;y_0;central;host_id]' tmpM.fit && rm ${memberscat}
funindex tmpM.fit[1] id
funindex ${1}.fit[1] $ID 
funjoin -m 1 -j1 $ID -j2 id ${1}.fit[1] tmpM.fit[1] ${memberscat} && rm tmpM.fit
fitscopy ${memberscat}'[col -id_2]' tmpM.fit && mv tmpM.fit ${memberscat}
rm -f *.idx

echo "== ranking clusters by Ngals"

# rank clusters by nvt
funtable -s "nvt sig" ${clustercat}[1] tmpC.fit
fitscopy tmpC.fit'[col *;RANK=(NAXIS2+1-#row)]' \!${clustercat}
funtable -s "RANK" ${clustercat}[1] tmpC.fit && mv tmpC.fit ${clustercat}

echo "== add rank column to members catalog"

# add rank column to members catalog
fitscopy ${clustercat}'[col RANK;id]' tmpC.fit
funindex tmpC.fit[1] id
funindex ${memberscat}[1] host_id
funjoin -m 1 -j1 host_id -j2 id ${memberscat}[1] tmpC.fit[1] tmpM.fit
fitscopy tmpM.fit'[col -id_2]' \!${memberscat} && rm tmpM.fit
rm -f tmpC_id.idx ${memberscat}_host_id.idx tmpC.fit

echo "== create input file for matchscript.py"

# create input file for matchscript.py
fitscopy ${clustercat}'[col RANK;ra_c;dec_c;redsh=z]' tmpC.fit
fitscopy ${memberscat}'[col gid(K)=id;RANK]' tmpM.fit
fitscopy tmpC.fit'[col RANK;RA(1D)=ra_c;DEC(1D)=dec_c;Z=redsh]' tmp.fit && mv tmp.fit tmpC.fit
fitscopy tmpM.fit'[col ID(K)=gid;RANK]' tmp.fit && mv tmp.fit tmpM.fit
funtable tmpC.fit[1] $rankedlist "RANK RA DEC Z" && rm tmpC.fit
funtable -a tmpM.fit[1] $rankedlist "ID RANK" && rm tmpM.fit
rm -f *.idx

echo "== cluster catalog done."
echo `date`

exit 

