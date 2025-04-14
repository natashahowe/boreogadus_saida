#!/bin/bash

JOBS_FILE=/home/mchale/Documents/boreogadus/scripts/boreogadus_alignARRAY_input.txt
IDS=$(cat ${JOBS_FILE})

for sample_line in ${IDS}
do
	job_index=$(echo ${sample_line} | awk -F ":" '{print $1}')
	fq_r1=$(echo ${sample_line} | awk -F ":" '{print $2}')
	fq_r2=$(echo ${sample_line} | awk -F ":" '{print $3}')

    sample_id=$(echo $fq_r1 | sed 's!^.*/!!')
    sample_id=${sample_id%%_*}

    bwa mem -M -t 10 /home/mchale/Documents/boreogadus/bwa/Arctic_cod_genome ${fq_r1} ${fq_r2} 2> /home/mchale/Documents/boreogadus/bwa/boreogadus_${sample_id}_bwa-mem.out > /home/mchale/Documents/boreogadus/bamtools/boreogadus_${sample_id}.sam

    samtools view -bS -F 4 /home/mchale/Documents/boreogadus/bamtools/boreogadus_${sample_id}.sam > /home/mchale/Documents/boreogadus/bamtools/boreogadus_${sample_id}.bam
    rm /home/mchale/Documents/boreogadus/bamtools/boreogadus_${sample_id}.sam

    samtools view -h /home/mchale/Documents/boreogadus/bamtools/boreogadus_${sample_id}.bam | samtools view -buS - | samtools sort -o /home/mchale/Documents/boreogadus/bamtools/boreogadus_${sample_id}_sorted.bam
    rm /home/mchale/Documents/boreogadus/bamtools/boreogadus_${sample_id}.bam

    samtools depth -aa /home/mchale/Documents/boreogadus/bamtools/boreogadus_${sample_id}_sorted.bam | cut -f 3 | gzip > /home/mchale/Documents/boreogadus/bamtools/boreogadus_${sample_id}_sorted.depth.gz

    samtools index /home/mchale/Documents/boreogadus/bamtools/boreogadus_${sample_id}_sorted.bam

done