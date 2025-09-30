#!/bin/bash
####
#  QC BAM alignment
# input; analysis folder path, sample name
#
####

echo "SAMPLE: $1"
echo "DIR: $2"

TXT_IN=${2}/raw_results/tables/${1}_identified_matched.txt
BAM_IN=${2}/aligned/${1}.st.bam
BAM_OUT=${2}/qc/${1}_identified_matched_sorted.bam
HIST=${2}/qc/${1}_identified_matched_coverage.txt
STATS=${2}/qc/${1}_aligned_stats.txt

cut -f1,2,3,4,5 $TXT_IN > ${2}/raw_results/tables/${1}_identified_matched.bed
samtools view -b -h -L ${2}/raw_results/tables/${1}_identified_matched.bed $BAM_IN | samtools sort > $BAM_OUT
samtools index $BAM_OUT
#samtools coverage -A --histogram $BAM_OUT > $HIST
samtools stats $BAM_IN | grep ^SN | cut -f 2- > $STATS

