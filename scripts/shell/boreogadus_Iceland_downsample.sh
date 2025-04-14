#!/bin/bash

#module unload bio/samtools/1.11
#module load bio/samtools/1.11

JOBS_FILE=/home/mchale/Documents/boreogadus/scripts/boreogadus_batch_downsample-bams_input.txt
IDS=$(cat ${JOBS_FILE})

for sample_line in ${IDS}
do
	sample_id=$(echo ${sample_line} | awk -F ":" '{print $2}')
	downsample_value=$(echo ${sample_line} | awk -F ":" '{print $3}')

	echo ${sample_id} ${downsample_value}
	
	samtools view -bo /home/mchale/Documents/boreogadus/downsample/boreogadus_${sample_id}_downsampled.bam -s ${downsample_value} /home/mchale/Documents/Arctic_cod/New_genome/Iceland_samples/${sample_id}.bam

	samtools depth /home/mchale/Documents/boreogadus/downsample/boreogadus_${sample_id}_downsampled.bam | awk -v name=${sample_id} '{sum+=$3} END {print name,"\t",sum/561633924}' > /home/mchale/Documents/boreogadus/depth/${sample_id}_downsampled_depth_20241214.csv
	
	samtools index /home/mchale/Documents/boreogadus/downsample/boreogadus_${sample_id}_downsampled.bam

done
