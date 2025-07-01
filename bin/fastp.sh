#!/bin/bash
# Script for fastp.sh
# fastp

# Activate environment
source ~/.bashrc
eval "$(micromamba shell hook --shell bash)"
micromamba activate bioinfo

# Make working directory
WORKDIR="/Users/sls6550/work/hTBA_Alc_Singh"

# Input directory
INDIR=${WORKDIR}/data

# Write output directories
OUTDIR=${INDIR}/fastp_filt  # Passed reads
FAILDIR=${INDIR}/fastp_fail  # Failed reads
mkdir -p ${OUTDIR}
mkdir -p ${FAILDIR} 

# Get list of sample names
ls $INDIR/*.fastq.gz | xargs -n 1 basename | awk -F '_L001_R' '{print $1}' | uniq > $INDIR/names.txt

# Run in parallel
cat $INDIR/names.txt | \
parallel --verbose \
/Users/sls6550/micromamba/envs/bioinfo/bin/fastp \
  -i $INDIR/{}_L001_R1_001.fastq.gz \
  -I $INDIR/{}_L001_R2_001.fastq.gz \
  -o $OUTDIR/{}_fastp_1.fastq \
  -O $OUTDIR/{}_fastp_2.fastq \
  --failed_out ${FAILDIR}/{}_failed.fastq \
  -q 20 \
  --cut_tail \
  --cut_tail_window_size 4 \
  --cut_tail_mean_quality 20 \
  -x \
  -p \
  --thread 4 \
  -j ${OUTDIR}/fastp.json \
  -h ${OUTDIR}/fastp.html

echo "Done with trimming"
