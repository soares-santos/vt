#!/bin/bash

# Source this script to setup environment for vt cluster finder.

# Uncomment this line to run the script in verbose mode:
# VERBOSE=1  

if (( $VERBOSE )) ; then echo " == set environment for vt runs ==" ; fi

export VT_DIR=`pwd`
if (( $VERBOSE )) ; then echo " ...... your vt working dir is: $VT_DIR" ; fi

CURRENT_PATH=$PATH

UNSETUP_SCRIPT="$VT_DIR/unsetup-vt.sh"
rm -rf $UNSETUP_SCRIPT
touch $UNSETUP_SCRIPT
echo "
# This is the unsetup-vt.sh script file. 
# Source this file to reset vt specific environment variables and 
# to remove vt specific bin directories from your PATH" >> $UNSETUP_SCRIPT
echo "export VT_DIR=\"\"" >> $UNSETUP_SCRIPT
echo "export PATH=$CURRENT_PATH" >> $UNSETUP_SCRIPT 
if (( $VERBOSE )); then echo " ...... file $UNSETUP_SCRIPT updated" ; fi


if [[ :$PATH: == *:"$VT_DIR/bin":* ]] ; then
    if (( $VERBOSE )) ; then echo " ...... dir $VT_DIR/bin found in PATH" ; fi
else
    export PATH=$VT_DIR/bin:$PATH
    if (( $VERBOSE )) ; then echo " ...... dir $VT_DIR/bin added to PATH" ; fi
fi

if [[ :$PATH: == *:".":* ]] ; then
    if (( $VERBOSE )) ; then echo " ...... dir ./ found in PATH" ; fi
else
    export PATH=.:$PATH
    if (( $VERBOSE )) ; then echo " ...... dir ./ added to PATH" ; fi
fi

if (( $VERBOSE )) ; then echo " == vt environment is set ==" ; fi

