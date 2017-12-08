#!/usr/bin/env bash

######################################################################
# 
# SCaMP setup script
#
# Sets up conda environment and installs prerequisite packages
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
	    conda remove -y --name SCaMP --all
	    if [ $? != 0 ]; then
			echo
			echo "An error occured removing the SCaMP environment..."
			echo
			exit 1
	    fi
	    break
	else
	    echo
	    echo "Leaving existing SCaMP environment in place..."
	    echo 
		KEEPENV=1
	    break
	fi
    done 
fi

if [ -z "${KEEPENV}" ]; then
	echo
	echo "Creating SCaMP environment..."
	echo
	conda create -y -n SCaMP
	if [ $? != 0 ]; then
    	echo
    	echo "An error occured creating the SCaMP environment..."
    	echo
    	exit 1
	fi

	echo "Activating environment..."
	source activate SCaMP
	if [ $? != 0 ]; then
		echo
		echo "An error occured activating the SCaMP environment..."
		echo
		exit 1
	fi

	echo
	echo "Installing prerequisite packages..."
	echo
	conda install -y --file ${DIR}/../etc/conda_packages.txt
	if [ $? != 0 ]; then
		echo
    	echo "An error occured installing the SCaMP prerequisites"
    	echo
    	exit 1
	fi
fi

PATHSET=$(which SCaMP 2>/dev/null|grep -c SCaMP)
if [ "${PATHSET} eq 0" ]; then
	echo
	echo "Would you like the SCaMP bin directory appending to your default path [y/N]?"
	echo
	while true; do
		read -rsn1 key
		if [ "${key,,}" = "y" ]; then
			if [ -e ~/.bashrc ]; then
				echo "#Added by SCaMP setup.sh" >> ~/.bashrc
				echo "export PATH=\"${DIR}:\$PATH\"" >> ~/.bashrc
				echo "Path added to ~/.bashrc..."
				echo "Logout then login again to ensure PATH is correctly set"
				echo
				break
			fi
		else
			break
		fi
	done
else 
	echo
	echo "Path already correctly set..."
	echo
fi

echo "setup complete..."
echo
