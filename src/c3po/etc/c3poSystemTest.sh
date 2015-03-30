#!/bin/bash

#===================================================================#
# sytsem settings test routine for c3po project 
# Stefan Radl, February 2015
#===================================================================#

#- include functions
source $CFDEM_SRC_DIR/lagrangian/cfdemParticle/etc/functions.sh

#- show gcc settings
checkGPP="false"

#- system settings
clear
echo "**********************"
echo "CPPPO system settings:"
echo "**********************"

echo "CFDEM_VERSION=$CFDEM_VERSION"
echo "couple to OF_VERSION=$WM_PROJECT_VERSION"
echo "compile option=$WM_COMPILE_OPTION"

echo ""
echo "checking environment variables..."
checkEnvComment "$USEHDF5" '$USEHDF5' "yes"

echo ""
echo "checking if paths are set correctly..."
checkDirComment "$C3PO_SRC_DIR" '$C3PO_SRC_DIR' "yes"
checkDirComment "$C3PO_QT5_DIR" '$C3PO_QT5_DIR' "yes"
checkDirComment "$C3PO_HDF5_DIR" '$C3PO_HDF5_DIR' "yes"
checkDirComment "$C3PO_HDF5_LIB" '$C3PO_HDF5_LIB' "yes"
checkDirComment "$C3PO_CHEMKIN_SRC_DIR" '$C3PO_CHEMKIN_SRC_DIR' "yes"
checkDirComment "$C3PO_EXAMPLE_DIR" '$C3PO_EXAMPLE_DIR' "yes"
checkDirComment "$C3PO_QT5_LIB" '$C3PO_QT5_LIB' "yes"
checkDirComment "$C3PO_QT5_INC" '$C3PO_QT5_INC' "yes"
checkDirComment "$C3PO_HDF5_INC" '$C3PO_HDF5_INC' "yes"
echo ""

echo "relevant library names:"
echo "*******************"

echo ""
echo "relevant aliases for c3po can be found here:"
echo $C3PO_SRC_DIR/etc/bashrc
echo "*******************"
echo ""

if [ $checkGPP == "true" ]
  then
    echo "g++:"
    which g++
    g++ --version

    echo "gcc:"
    which gcc
    gcc --version

    echo "mpic++:"
    which mpic++
    mpic++ --version

    echo "mpirun:"
    which mpirun
    mpirun --version
fi

