#!/bin/bash
# script for Quality Checking
# fastqc and multiqc

# activate environment
source ~/.bashrc
eval "$(micromamba shell hook --shell bash)"
micromamba activate bioinfo

# set paths
RAW=/Users/sls6550/work/hTBA_Alc_Singh/data

# create directory for quality if needed
QUALITY=${RAW}/quality 
mkdir -p ${QUALITY}

# run fastqc on all files 
fastqc ${RAW}/*.fastq.gz -o ${QUALITY}

# If you encounter an error when run multiQC, you should run the following commands:
#export LC_ALL=en_US.utf-8
#export LANG=en_US.utf-8

# run MultiQC
multiqc ${QUALITY}/*.zip -o ${QUALITY} 
