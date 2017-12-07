#!/usr/bin/env bash

######################################################################
# 
# SCaMP setup script
#
# Sets up conda environment and installs prerequsitie packages
#
######################################################################

CONDA=$(which conda)
DIR="$( cd "$( dirname "$0" )" ; pwd  -P)"

if [ -z ${CONDA} ]; then
    echo
    echo "conda not found. Please install conda and bioconda before running this script"
    echo "See https://bioconda.github.io/#using-bioconda for installation instructions"
    echo
    exit 1
fi 

echo "Checking conda environments..."
ENV_EXISTS=$(conda env list|grep -c SCaMP)

if [ ${ENV_EXISTS} == 1 ]; then
    echo "Existing SCaMP conda environment found"
    echo "Do you wish to remove this environment and create a new one [y/N]?"
    while true; do
	read -rsn1 key
	if [ "${key,,}" = "y" ]; then
	    echo 
	    echo "Removing SCaMP environment..."
	    echo
	    #$(conda remove --name SCaMP --all)
	    #if [ $? != 0 ]; then
		#echo
		#echo "An error occured removing the SCaMP environment..."
		#echo
		#exit 1
	    #fi
	    break
	else
	    echo
	    echo "Leaving existing SCaMP environment and exiting..."
	    echo 
	    exit 1
	fi
    done 
fi


echo
echo "Creating SCaMP environment..."
echo
$(conda create -n SCaMP)
if [ $? != 0 ]; then
    echo
    echo "An error occured removing the SCaMP environment..."
    echo
    exit 1
fi

echo "Activating environment..."
$(source activate SCaMP)
echo
if [ $? != 0 ]; then
    echo
    echo "An error occured removing the SCaMP environment..."
    echo
    exit 1
fi

echo
echo "Installing prerequisite packages..."
echo
$(conda install -f ${DIR}../etc/conda_packages.txt)
if [ $? != 0 ]; then
    echo
    echo "An error occured removing the SCaMP environment..."
    echo
    exit 1
fi
