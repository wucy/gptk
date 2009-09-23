#!/bin/bash
# Check that source folder is set up - if not, set it up
# This should only be needed once.
SRC_FOLDER=`pwd`/psgp/src
CUR_FOLDER=`pwd`
if [[ ! -f $SRC_FOLDER/Makefile ]]
then
    echo Source folder not setup. Running setup script.
    cd $SRC_FOLDER
    ./setup.sh
    cd $CUR_FOLDER
    echo Done.
else
    echo Source folder setup.
fi

# Build the R package
rm psgp_*.tar.gz     # Remove previous versions
R CMD build psgp     # Build package

