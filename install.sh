#!/bin/bash

## This script will install the tools required for the STRling pipelines.
## It will fetched each tool from the web and placed into the tools/ subdirectory.
## Paths to all installed tools can be found in the file tools.groovy at the
## end of execution of this script. These paths can be changed if a different
## version of software is required.
##

installdir=$PWD
refdir=$PWD/reference-data
toolspec=$PWD/pipelines/pipeline_config.groovy
bpipeconfig=$PWD/pipelines/bpipe.config
bpipeconfig_template=$PWD/pipelines/config-examples/bpipe.config_template

mkdir -p tools/bin
cd tools

if [ ! -f $toolspec ] ; then
    echo "The configuration file pipelines/pipeline_config.groovy was not found, creating it."
    else
    echo -n "WARNING: pipelines/pipeline_config.groovy already exists so will be overwritten by this installation. "
    echo "Creating backup pipelines/pipeline_config.groovy.backup in case you wish to retreive the previous version of this file."
    cp $toolspec ${toolspec}.backup
fi

#a list of which programs need to be installed
commands="strling bpipe python"

function strling_install {
    wget --no-check-certificate https://github.com/quinlan-lab/STRling/releases/latest/download/strling
    chmod +x strling
    ln -s $PWD/strling $PWD/bin/
}

function bpipe_install {
    wget -O bpipe-0.9.9.5.tar.gz https://github.com/ssadedin/bpipe/releases/download/0.9.9.5/bpipe-0.9.9.5.tar.gz
    tar -zxvf bpipe-0.9.9.5.tar.gz ; rm bpipe-0.9.9.5.tar.gz
    ln -s $PWD/bpipe-0.9.9.5/bin/* $PWD/bin/
}

# Installs miniconda, Python 3 + required packages
# (and any other dependancies listed in environment.yml)
function python_install {
    wget -O miniconda.sh https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
    bash miniconda.sh -b -p $PWD/miniconda
    rm miniconda.sh
    $PWD/miniconda/bin/conda env create -f ../scripts/environment.yml
    ln -s $PWD/miniconda/envs/strling/bin/* $PWD/bin/
#    source activate STR
}

#populate toolspec
echo "// Bpipe pipeline config file" > $toolspec
echo "// Paths are relative to the directory the pipeline is running in, so absolute" >> $toolspec
echo "// paths are recommended." >> $toolspec
echo >> $toolspec
echo "// Adjust parameters" >> $toolspec
echo "REF=\"path/to/reference_genome.fasta\"" >> $toolspec
echo >> $toolspec

#set STRetch base directory
echo "// STRling installation location" >> $toolspec
echo "STRLING=\"$installdir\"" >> $toolspec
echo >> $toolspec

echo "// Paths to tools used by the pipeline" >> $toolspec

for c in $commands ; do
    c_path=`which $PWD/bin/$c 2>/dev/null`
    if [ -z $c_path ] ; then
	echo "$c not found, fetching it"
	${c}_install
	c_path=`which $PWD/bin/$c 2>/dev/null`
    fi
    echo "$c=\"$c_path\"" >> $toolspec
done

echo >> $toolspec
echo "// By default, uses other samples in the same batch as a control" >> $toolspec
echo "CONTROL=\"\"" >> $toolspec
echo "// Uncomment the line below to use a set of WGS samples as controls, or specify your own" >> $toolspec
echo "//CONTROL=\"$refdir/controls\"" >> $toolspec
echo >> $toolspec

if [ ! -f $bpipeconfig ] ; then
    echo "pipelines/bpipe.config not found, creating it"
    #copy bpipe.config template to pipeline directory
    cp $bpipeconfig_template $bpipeconfig
fi

#loop through commands to check they are all installed
echo "**********************************************************"
echo "Checking that all required tools were installed:"
Final_message="All commands installed successfully!"
for c in $commands ; do
    c_path=`which $PWD/bin/$c 2>/dev/null`
    if [ -z $c_path ] ; then
	echo -n "WARNING: $c could not be found!!!! "
	echo "You will need to download and install $c manually, then add its path to $toolspec"
	Final_message="WARNING: One or more command did not install successfully. See warning messages above. You will need to correct this before running STRetch."
    else
        echo "$c looks like it has been installed"
    fi
done

echo "**********************************************************"

#check for reference data
if [ ! -f $REF ] ; then
    echo -n "WARNING: reference files could not be found!!!! "
    echo "You will need to download them manually, then add the path to $toolspec"
else
    echo "Reference data has been found"
fi

echo "**********************************************************"

#check for config files
if [ ! -f $toolspec ] ; then
    echo -n "WARNING: pipelines/pipeline_config.groovy could not be found!!!! "
    echo "You will need to create this file manually."
else
    echo "It looks like pipelines/pipeline_config.groovy exists"
fi

if [ ! -f $bpipeconfig ] ; then
    echo -n "WARNING: pipelines/bpipe.config could not be found. "
    echo -n "This file is only required when using a job submission system. "
    echo "If required, you will need to create this file manually."
else
    echo "It looks like pipelines/bpipe.config exists"
fi

echo "**********************************************************"
echo $Final_message
