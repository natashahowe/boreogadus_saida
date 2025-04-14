#!/bin/bash

# Natasha Howe
# New boreogadus reference genome

# ANGSD Version 0.930
## Specify directories with data
REFGENOME=/home/mchale/Documents/Arctic_cod/New_genome/Arctic_cod_genome.fna
BAMDIR=/home/mchale/Documents/boreogadus
ANGSDDIR=/home/mchale/Desktop/angsd/

# tab-delimited text file with two columns, 1) pop name, 2) sample size
POP_FILE=${BAMDIR}/scripts/boreogadus_pop_samplesize.txt

#Filters:
FILTERS="-minQ 30 -minMapQ 30 -C 50 -remove_bads 1 -uniqueOnly 1 -skipTriallelic 1" # got rid of minMaf and snp-pval --> using filtered sites flag
#Options:
OPTIONS="-GL 1 -doSaf 1 -doMaf 2 -doMajorMinor 1 -doCounts 1 -doPost 1 -doGlf 2"

cd ${ANGSDDIR}
./angsd sites index ${BAMDIR}/gls/boreogadus_filtered_SNPs_filetwo.txt

while IFS=$'\t' read POP samplesize REST;
do echo $POP $samplesize;
    mindepth=$((samplesize*2))
    maxdepth=$((samplesize*20))

    ./angsd -nthreads 12 \
	-bam ${BAMDIR}/boreogadus_${POP}_bams.txt \
	-anc ${REFGENOME} -ref ${REFGENOME} \
	-sites ${BAMDIR}/gls/boreogadus_filtered_SNPs_filetwo.txt \
	-out ${BAMDIR}/gls/boreogadus_${POP}_20241210 \
	${FILTERS} ${OPTIONS} \
    	-setminDepth ${mindepth} -setmaxDepth ${maxdepth} -minInd 7 #chose minInd7 to do fst for Labrador, too

done < $POP_FILE
