#!/bin/bash

# 2024-12-13
# depth with samtools

# set directories
BAMFILE=/home/mchale/Documents/boreogadus/scripts/boreogadus_bamlist_plus_downsampled.txt 	# path to textfile that lists all bamfiles with their path
OUTFILE=/home/mchale/Documents/boreogadus/bam/ # directory where depths per individual will be places
ARCTICBAM=/home/mchale/Documents/Arctic_cod/New_genome/Iceland_samples

# ARCTIC COD GENOME SIZE = 561633924

while read sample_id
do
	cd ${ARCTICBAM}
	echo ${sample_id}
	sample=$(basename -- $sample_id) # removed the filepath from the name
	filename="${sample%.*}" # removed the ".bam" from the name
	echo ${filename}
	samtools depth ${sample_id} | awk -v name=${filename} '{sum+=$3} END {print name,"\t",sum/561633924}' > ${OUTFILE}/${filename}_depth_20241213.csv
done < ${BAMFILE} # this is where it calls in the bamfile

# Create new CSV and then write all depths to it
touch ${OUTFILE}/boreogadus_depths_20241213.csv
cat ${OUTFILE}/*_depth_20241213.csv >> ${OUTFILE}/boreogadus_depths_20241213.csv

