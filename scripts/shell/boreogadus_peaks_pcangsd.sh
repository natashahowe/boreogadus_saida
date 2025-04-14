#!/bin/bash

# Natasha Howe

PCADIR=/home/mchale/Documents/boreogadus/pca
PCANGSD=/home/mchale/programsNH/pcangsd/pcangsd
GLSDIR=/home/mchale/Documents/boreogadus/gls

cd $PCANGSD

CHRNAME=OZ177918.1
CHR=15
STARTPOS=7500000
ENDPOS=20000000

cd ${PCANGSD} 

PEAKFILE=/home/mchale/Documents/boreogadus/scripts/boreogadus_peakfile1.txt 

#while IFS=$'\t' read CHRNAME CHR STARTPOS ENDPOS REST; 
#do

FILENAME=boreogadus_chr${CHR}_s${STARTPOS}_e${ENDPOS}
zcat ${GLSDIR}/boreogadus_${CHRNAME}_allpops_20241216.beagle.gz | awk -v s=$STARTPOS -v e=$ENDPOS -F'[\t_]' ' $2 >= s && $2 <= e' | gzip > ${GLSDIR}/peaks/${FILENAME}.beagle.gz

pcangsd --threads 4 --beagle ${GLSDIR}/peaks/${FILENAME}.beagle.gz --out ${PCADIR}/peaks/${FILENAME}

#done < ${PEAKFILE}