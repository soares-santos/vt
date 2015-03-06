
## Usage: source setup-vt.sh
## This script will setup environment variables for vt cluster finder.

## say hello and get started

if (( $VERBOSE )) ; then echo " == set environment for vt runs ==" ; fi

## export path to your VT installation

export VT_DIR=`pwd`
if (( $VERBOSE )) ; then echo " ...... your vt working dir is: $VT_DIR" ; fi

## create unsetup script in case you want to remove vt from your system

CURRENT_PATH=$PATH
UNSETUP_SCRIPT="$VT_DIR/run/unsetup-vt.sh"
rm -rf $UNSETUP_SCRIPT
touch $UNSETUP_SCRIPT
echo "
# This is the unsetup-vt.sh script file. 
# Source this file to reset vt specific environment variables and 
# to remove vt specific bin directories from your PATH" >> $UNSETUP_SCRIPT
echo "export VT_DIR=\"\"" >> $UNSETUP_SCRIPT
echo "export PATH=$CURRENT_PATH" >> $UNSETUP_SCRIPT 
if (( $VERBOSE )); then echo " ...... file $UNSETUP_SCRIPT updated" ; fi

## add VT bin directory to your path

if [[ :$PATH: == *:"$VT_DIR/bin":* ]] ; then
    if (( $VERBOSE )) ; then echo " ...... dir $VT_DIR/bin found in PATH" ; fi
else
    export PATH=$VT_DIR/bin:$PATH
    if (( $VERBOSE )) ; then echo " ...... dir $VT_DIR/bin added to PATH" ; fi
fi

## add "." to your path

if [[ :$PATH: == *:".":* ]] ; then
    if (( $VERBOSE )) ; then echo " ...... dir ./ found in PATH" ; fi
else
    export PATH=.:$PATH
    if (( $VERBOSE )) ; then echo " ...... dir ./ added to PATH" ; fi
fi

## say goodbye and move on

if (( $VERBOSE )) ; then echo " == vt environment is set ==" ; fi
