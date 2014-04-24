#!/bin/bash

usage() {
    echo "
$0 is a install script for the vt cluster finder.

Usage: $0 [-h] [-q] [-o] [-r] [<target_directory>]

A vt directory will be created under the <target_directory>   
and the vt code will be installed there. If omitted, the 
<target_directory> is set to the current directory.

After installing, setup the environment for your vt runs
with the command:
source <target-directory>/vt/setup-vt.sh 

Options:
-h    Print this help.
-o    Install overwriting a previous installation.
-q    Quiet install.  
-r    Reinstall. Only vt code, not the libraries.
" 1>&2
    exit 1
}

OVERWRITE=0
VERBOSE=1 
REINSTALL=0

nargs="$#"

while getopts ":hoqr" opt; do
    case $opt in 
	h) 
	    usage
	    ;;
	o)
	    OVERWRITE=1
	    (( nargs -=1 ))
	    ;;
	q)
	    VERBOSE=0
	    (( nargs -=1 ))
	    ;;	    
	r)
	    REINSTALL=1
	    (( nargs -=1 ))
	    ;;	    
	*) 
	    echo "Invalid option."  
	    usage
	    ;;
    esac
done

if [ "$nargs" -eq 0 ]; then 
    TARGET_DIR=`pwd`
else
    TARGET_DIR="${@: -1}" 
    if [ ! -d $TARGET_DIR ]; then usage ; fi
fi
	    
if (( $VERBOSE )) ; then 
    echo "== set parameters == "
    echo " ...... This are your install.sh parameters: 
        TARGET_DIR = $TARGET_DIR
        VERBOSE = $VERBOSE
        OVERWRITE = $OVERWRITE
        REINSTALL = $REINSTALL"
fi

dir=$TARGET_DIR/vt
export VT_DIR=$dir
if (( $VERBOSE )) ; then echo "== create dir $dir  ==" ; fi
if [ -d "$dir" ]; then 
    if (( $REINSTALL )) ; then
	if (( $VERBOSE )) ; then echo " ...... $dir already exists. Reinstall vt here." ; fi
    else
	if (( $OVERWRITE )) ; then
	    if (( $VERBOSE )) ; then echo " ...... overwriting previous vt installation." ; fi
	    rm -rf $dir || exit 1
	    mkdir $dir 
	    mkdir $dir/bin
	else
	    echo " Error: Directory $dir already exists."
	    echo "        Use option -o to overwrite $dir and ALL its contents."
	    echo "        Use option -r to reinstall only the vt code."
	    exit 1
	fi
    fi
else
    mkdir $dir || exit 1
    mkdir $dir/bin
fi

INSTALL_LOG=$dir/install.log
rm -rf $INSTALL_LOG 
touch $INSTALL_LOG
if (( $VERBOSE )) ; then 
    echo "== perform installation  ==" 
    echo " ...... This is your install logfile:"
    echo "        INSTALL_LOG=$INSTALL_LOG"
fi

if (( $REINSTALL == 0 )) ; then
    if [ -d "$CFITSIO_DIR" ] ; then 
	if (( $VERBOSE )) ; then echo "== using cfitsio found in $CFITSIO_DIR ==" ; fi
    else
	if (( $VERBOSE )) ; then echo "== make cfitsio library ==" ; fi
	cd utils/cfitsio
	./configure --prefix=$dir >> $INSTALL_LOG 2>&1
	make >> $INSTALL_LOG 2>&1
	make install >> $INSTALL_LOG 2>&1
	make clean >> $INSTALL_LOG 2>&1
	if [ ! -f "$dir/lib/libcfitsio.a" ]; then  
	    echo "Error: libcifitsio.a not found. Did cfitsio install fail?"
	    exit 1
	fi
	cd - >> /dev/null 2>&1
    fi
    if [ -d "$FUNTOOLS_DIR" ] ; then
	if (( $VERBOSE )) ; then echo "== using funtools found in $FUNTOOLS_DIR ==" ; fi
    else	    
	if (( $VERBOSE )) ; then  echo "== make funtools library ==" ; fi
	cd utils/funtools
	./configure --prefix=$dir >> $INSTALL_LOG 2>&1
	make >> $INSTALL_LOG 2>&1
	make install >> $INSTALL_LOG 2>&1
	make clean >> $INSTALL_LOG 2>&1
	if [ ! -f "$dir/lib/libfuntools.a" ]; then  
	    echo "Error: libfuntools.a not found. Did funtools install fail?"
	    exit 1
	fi
	cd - >> /dev/null 2>&1
    fi
fi

if (( $VERBOSE )) ; then echo "== make vt ==" ; fi
cd src
make >> $INSTALL_LOG 2>&1
if [ ! -f "$dir/bin/fitscopy" ]; then  
    echo "Error: vtfind command not found. Did vt install fail?"
    exit 1
fi
cd - >> /dev/null 2>&1
cp setup-vt.sh $dir
sed -i "s|\`pwd\`|$VT_DIR|" $dir/setup-vt.sh
cp -r example $dir
cp README $dir

if (( $VERBOSE )) ; then echo "== install complete ==" ; fi

exit 0