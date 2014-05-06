#!/bin/bash

#set -x

if [ "$#" -lt 1 ] ; then echo "Usage: findclusters <galaxy_cat>" ; exit ; fi

echo `date`
echo "Starting the VT cluster finder."

glx_cat=$1  # input galaxy catalog filename
RA=RA       # RA colname
DEC=DEC     # DEC colname
Z=Z         # redshift colname  
ID=ID       # ID colname
rcl=0.0     # min confidence level 
boost=1.    # 1 - (boost * frac) is the threshold for detection
err_z=0.03  # err_z * (1+z) is the half width of the zphot shells in the second run
# frac is a function of redshift:
redshift=(0.0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0 1.1 1.2 1.3 1.4 1.5)
 frac[0]=0.0
 frac[1]=0.00002247
 frac[2]=0.00237922
 frac[3]=0.00556252
 frac[4]=0.00716849
 frac[5]=0.00930243
 frac[6]=0.01080613
 frac[7]=0.01289336
 frac[8]=0.01566860
 frac[9]=0.01646563
frac[10]=0.01788454
frac[11]=0.01613668
frac[12]=0.01201421
frac[13]=0.00872151
frac[14]=0.00491137
frac[15]=0.00150253
# user may want to modify the params above:
if [ -f params.info ] ; 
    then source params.info
    else echo "Warning: File params.info not found. Using defaults."
fi
echo "Parameters set:"
echo "RA=$RA"
echo "DEC=$DEC"
echo "Z=$Z"
echo "ID=$ID"
echo "rcl=$rcl"
echo "boost=$boost"
echo "err_z=$err_z"
echo "glx_cat=$glx_cat"
if [ ! -f boxes.info ] ; then echo "Error: File boxes.info not found.";exit ; fi
if [ ! -f $VT_DIR/bin/vt ] ; then echo "Error: Executable file vt not found.";exit; fi
if [ ! -f $glx_cat ] ; then echo "Error: File $glx_cat not found.";exit ; fi

cls_cat=${glx_cat}_VTclusters.fit
rm -f $cls_cat
run1=${glx_cat}.run1
run2=${glx_cat}.run2
rm -rf $run1 $run2
mkdir $run1 $run2

base=`pwd`
rm -f boxes
cat boxes.info | grep -v '#' > boxes
nboxes=`cat boxes | wc -l` 
echo "$nboxes boxes to process at Run 1."
echo "Processing..."
n=0
while [ $n -lt $nboxes ] 
do
    n=`expr $n + 1` 
    tmp=`head -$n boxes | tail -1`
    id=`echo $tmp | awk '{print $1}'`
    ra_min=`echo $tmp | awk '{print $2}'`
    ra_max=`echo $tmp | awk '{print $3}'`
    dec_min=`echo $tmp | awk '{print $4}'`
    dec_max=`echo $tmp | awk '{print $5}'`    
    z_min=`echo $tmp | awk '{print $6}'`
    z_max=`echo $tmp | awk '{print $7}'`
    A=`echo $tmp | awk '{print $8}'`
    g=`echo $tmp | awk '{print $9}'`
    tmp=$glx_cat".box_"$id
    rm -rf $tmp 
    mkdir $tmp
    filter="$RA>$ra_min-1.0/cos($dec_min*#deg)&&$RA<$ra_max+1.0/cos($dec_min*#deg)&&$DEC>$dec_min-1.0&&$DEC<$dec_max+1.0&&$Z>$z_min&&$Z<$z_max"    
    `fitscopy ''$1'['$filter']' $tmp/$tmp.fit`
    cd $tmp
    scl=0.9
    ired=0
    nred=15
    while [ $ired -lt $nred ]
    do
	jred=`expr $ired + 1`
	bool1=`echo "( $z_min + $z_max ) * 0.5 > ${redshift[$ired]}" | bc`
	bool2=`echo "( $z_min + $z_max ) * 0.5 < ${redshift[$jred]}" | bc`
 	if [ $bool1 -eq 1 ] ; then
	if [ $bool2 -eq 1 ] ; then
	    scl=`echo "1.0 - ( $boost * ( ${frac[$ired]} + ${frac[$jred]} ) * 0.5 )" | bc`
	    ired=$nred
	fi
	fi
	ired=`expr $ired + 1`
    done
    `vt -v -A $A -g $g -r $rcl -s $scl -i $id -F $ra_min $ra_max $dec_min $dec_max -C $ID $RA $DEC $Z $tmp.fit > $tmp.log`
    if [ -f ${tmp}_VTclusters.fit ] 
    then 
	cp ${tmp}_VTclusters.fit ${tmp}_VTgalaxies.fit $base/$run2
	cd $base/$run2 
	if [ -f ${glx_cat}.candidates.fit ]
	then 
	    tabmerge ${tmp}_VTclusters.fit+1 ${glx_cat}.candidates.fit+1
	    tabmerge ${tmp}_VTgalaxies.fit+1 ${glx_cat}.membercand.fit+1
 	    rm ${tmp}_VTclusters.fit ${tmp}_VTgalaxies.fit
	else 
	    mv ${tmp}_VTclusters.fit ${glx_cat}.candidates.fit
 	    mv ${tmp}_VTgalaxies.fit ${glx_cat}.membercand.fit
	fi
    fi
    cd $base
    mv $tmp $run1
done
rm boxes
echo "Run 1 complete."
echo `date`

cp $glx_cat $run2
cd $run2
cp ${glx_cat}.candidates.fit ${glx_cat}.membercand.fit $base
`fitscopy ${glx_cat}.'candidates.fit[col id;ra_min=ra_c-1.5/cos(dec_c*#deg);ra_max=ra_c+1.5/cos(dec_c*#deg);dec_min=dec_c-1.5;dec_max=dec_c+1.5;z_min=z_c-'$err_z'*(1+z_c);z_max=z_c+'$err_z'*(1+z_c);w_amp=A;w_pow=g]' run2.boxes.fit`
`dumpfits 'run2.boxes.fit[1]' > boxes`

nboxes=`cat boxes | wc -l` 
echo "$nboxes cluster candidates to process at Run 2."
if [ $nboxes -lt 1 ] ; then exit ; fi
echo "Processing..."
n=0
while [ $n -lt $nboxes ] 
do
    n=`expr $n + 1`
    tmp=`head -$n boxes | tail -1`
    id=`echo $tmp | awk '{print $1}'`
    ra_min=`echo $tmp | awk '{print $2}'`
    ra_max=`echo $tmp | awk '{print $3}'`
    dec_min=`echo $tmp | awk '{print $4}'`
    dec_max=`echo $tmp | awk '{print $5}'`    
    z_min=`echo $tmp | awk '{print $6}'`
    z_max=`echo $tmp | awk '{print $7}'`
    A=`echo $tmp | awk '{print $8}'`
    g=`echo $tmp | awk '{print $9}'`
    tmp=$glx_cat".box_"$id
    rm -rf $tmp 
    mkdir $tmp
    filter="$RA>$ra_min&&$RA<$ra_max&&$DEC>$dec_min&&$DEC<$dec_max&&$Z>$z_min&&$Z<$z_max"    
    `fitscopy ''$1'['$filter']' $tmp/$tmp.fit`
    cd $tmp
    scl=0.9
    ired=0
    nred=15
    while [ $ired -lt $nred ]
    do
	jred=`expr $ired + 1`
	bool1=`echo "( $z_min + $z_max ) * 0.5 > ${redshift[$ired]}" | bc`
	bool2=`echo "( $z_min + $z_max ) * 0.5 < ${redshift[$jred]}" | bc`
 	if [ $bool1 -eq 1 ] ; then
	if [ $bool2 -eq 1 ] ; then
	    scl=`echo "1.0 - ( $boost * ( ${frac[$ired]} + ${frac[$jred]} ) * 0.5 )" | bc`
	    ired=$nred
	fi
	fi
	ired=`expr $ired + 1`
    done
    `vt -v -A $A -g $g -r $rcl -s $scl -C $ID $RA $DEC $Z $tmp.fit > $tmp.log`
    if [ -f ${tmp}_VTclusters.fit ] 
    then 
	mv ${tmp}_VTclusters.fit VTclusters.tmp.1.fit
	filter="dist2boxcenter=angsep(ra_c,dec_c,0.5*($ra_min+$ra_max),0.5*($dec_min+$dec_max))"
	`fitscopy VTclusters.tmp.1.fit'[col *;'$filter']' VTclusters.tmp.2.fit`
	`funtable -s dist2boxcenter VTclusters.tmp.2.fit[1] VTclusters.tmp.1.fit`
	id2=`dumpfits 'VTclusters.tmp.1.fit[col id][#row==1]' | awk '{print $1'}`
	if [ "$id2" = "" ]
	then
	    echo "Warning: Run2 id variable is empty. No cluster found in in this box ?"
	else
	    filter="host_id==$id2"
	    `fitscopy ${tmp}_VTgalaxies.fit'['$filter']' VTgalaxies.tmp.1.fit`
	    filter="id_run1(K)=$id"
	    `fitscopy VTclusters.tmp.1.fit'[#row==1][col *;'$filter']' \!VTclusters.tmp.2.fit`
	    `fitscopy VTgalaxies.tmp.1.fit'[col *;'$filter']' VTgalaxies.tmp.2.fit`	
	    if [ -f $base/${glx_cat}.candidates.run2.fit ]
	    then 
		tabmerge VTclusters.tmp.2.fit+1 $base/${glx_cat}.candidates.run2.fit+1
		tabmerge VTgalaxies.tmp.2.fit+1 $base/${glx_cat}.membercand.run2.fit+1
	    else 
		mv VTclusters.tmp.2.fit $base/${glx_cat}.candidates.run2.fit
		mv VTgalaxies.tmp.2.fit $base/${glx_cat}.membercand.run2.fit
	    fi
	fi
	cd $base
    fi
    cd $base/$run2
done
cd $base
rm -f tmp.fit
fitscopy $base/${glx_cat}.candidates.run2.fit'[col -id]' tmp.fit && mv tmp.fit  $base/${glx_cat}.candidates.run2.fit
fitscopy $base/${glx_cat}.candidates.run2.fit'[col *;id(K)=id_run1]' tmp.fit && mv tmp.fit  $base/${glx_cat}.candidates.run2.fit
fitscopy $base/${glx_cat}.candidates.run2.fit'[col -id_run1]' tmp.fit && mv tmp.fit  $base/${glx_cat}.candidates.run2.fit
fitscopy $base/${glx_cat}.membercand.run2.fit'[col -host_id]' tmp.fit && mv tmp.fit  $base/${glx_cat}.membercand.run2.fit
fitscopy $base/${glx_cat}.membercand.run2.fit'[col *;host_id(K)=id_run1]' tmp.fit && mv tmp.fit  $base/${glx_cat}.membercand.run2.fit
fitscopy $base/${glx_cat}.membercand.run2.fit'[col -id_run1]' tmp.fit && mv tmp.fit  $base/${glx_cat}.membercand.run2.fit
echo "Run 2 complete."
echo `date`

match_clusters -v -C id host_id ${glx_cat}.candidates.fit ${glx_cat}.membercand.fit ${glx_cat}.candidates.run2.fit ${glx_cat}.membercand.run2.fit
echo `date`

exit
