#!/bin/bash

# Natasha Howe
# New boreogadus reference genome

# ANGSD Version 0.930
## Specify directories with data
REFGENOME=/home/mchale/Documents/Arctic_cod/New_genome/Arctic_cod_genome.fna
BAMDIR=/home/mchale/Documents/boreogadus
ANGSDDIR=/home/mchale/Desktop/angsd/

#Filters:
FILTERS="-minQ 30 -minMapQ 30 -minMaf 0.01 -C 50 -remove_bads 1 -uniqueOnly 1 -skipTriallelic 1 -only_proper_pairs 1 -SNP_pval 1e-6"
#Options:
OPTIONS="-GL 1 -doSaf 1 -doMaf 1 -doMajorMinor 1 -doCounts 1 -doPost 1 -doGlf 2"

samplesize=164

cd ${ANGSDDIR}

    #mindepth=$((samplesize*2))
    #maxdepth=$((samplesize*20))
    #minind=$((samplesize / 3))

    ./angsd -nthreads 12 \
	-bam ${BAMDIR}/boreogadus_filtered_bams.txt \
	-anc ${REFGENOME} -ref ${REFGENOME} \
	-out ${BAMDIR}/gls/boreogadus_allPops_downsampleIce_20241214 \
	-GL 1 -doSaf 1 -doMaf 1 -doMajorMinor 1 -doCounts 1 -doPost 1 -doGlf 2 \
	-minQ 30 -minMapQ 30 -minMaf 0.01 -remove_bads 1 -uniqueOnly 1 -skipTriallelic 1 -only_proper_pairs 1 -SNP_pval 1e-6 \
	-setMinDepth 328 -setMaxDepth 3280 -minInd 54

