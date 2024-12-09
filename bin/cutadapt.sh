#!/bin/bash
## basic bash file for cutadapt on all samples

#activate environment using micromamba
source ~/.bashrc
eval "$(micromamba shell hook --shell bash)"
micromamba activate bioinfo
set -uex

#change to working directory
cd /Users/sls6550/"OneDrive - The Pennsylvania State University"/Ganda_Singh_Bioinfo/16s_hTBA

# working directory for sample names
WORKDIR=Sam_16S
mkdir -p $WORKDIR

#input directory for files
INDIR=$WORKDIR/Data

# output directory
OUTDIR=$WORKDIR/adapter-removed/
mkdir -p $OUTDIR

# fastqc directory
QUALDIR=$WORKDIR/quality_postadapt
mkdir -p $QUALDIR

# adapter sequences for cutadapt
AD1=AGATCGGAAGAGCACACGTCTGAACTCCAGTCA
AD2=AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT

# get list of sample names
ls $INDIR | awk -F '_L001_R' '{print $1}' | uniq > $WORKDIR/names.txt

## ---- cutadapt ----

# run in parallel

cat $WORKDIR/names.txt |\
parallel -j10 /Users/sls6550/micromamba/envs/bioinfo/bin/cutadapt \
-j 0 \
-a $AD1 -A $AD2 \
--max-n=0 \
-o $OUTDIR/{}_trimmed_1.fastq \
-p $OUTDIR/{}_trimmed_2.fastq \
$INDIR/{}_L001_R1_001.fastq.gz \
$INDIR/{}_L001_R2_001.fastq.gz

# print finished message
echo "cutadapt done"

## ---- fastqc on trimmed reads ----

#run fastqc on all reads
fastqc $OUTDIR/*.fastq -o $QUALDIR

# get multiqc report for all fastqc
multiqc $QUALDIR/* -o $QUALDIR -n adapter-removed.multiqc

# print finished message
echo "quality check done"
