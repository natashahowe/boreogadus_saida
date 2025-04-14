#!/bin/bash

# 2024-11-25
# depth with samtools

# you first need to create a textfile with 1 line per bamfile.
	# example line: /home/mchale/Documents/Arctic_cod/sample1528_sorted.bam
# to do this, you can type:
	# ls -d $PWD/*_sorted.bam > arctic_cod_bamfile_list.txt
# if bams are in multiple folders, you can do this within each folder and concat:
	# cat bamfile_list_folder1.txt bamfile_list_folder2.txt >> arctic_cod_bamfile_list.txt

# find depth
BAMFILE=/home/mchale/Documents/New_Arctic_cod/arctic_cod_bamfile_list.txt 	# path to textfile that lists all bamfiles with their path
OUTFILE=/home/mchale/Documents/Depth_calcs/ # directory where depths per individual will be places

# ARCTIC COD GENOME SIZE = 561633924

while read sample_id
do
	cd ${ARCTICBAM}
	echo ${sample_id}
	sample=$(basename -- $sample_id) # removed the filepath from the name
	filename="${sample%.*}" # removed the ".bam" from the name
	echo ${filename}
	samtools depth ${sample_id} | awk -v name=${filename} '{sum+=$3} END {print name,"\t",sum/561633924}' > ${OUTFILE}/${filename}_depth_2024-11-25.csv
done < ${BAMFILE} # this is where it calls in the bamfile

# Create new CSV and then write all depths to it
touch ${OUTFILE}/arctic_cod_depths_2024-11-25.csv
cat ${OUTFILE}/*_depth_2024-11-25.csv >> ${OUTFILE}/arctic_cod_depths_2024-11-25.csv

