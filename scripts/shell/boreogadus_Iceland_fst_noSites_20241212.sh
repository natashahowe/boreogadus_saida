#!/bin/bash

# Natasha Howe
# New boreogadus reference genome

# ANGSD Version 0.930
## Specify directories with data
REFGENOME=/home/mchale/Documents/Arctic_cod/New_genome/Arctic_cod_genome.fna
BAMDIR=/home/mchale/Documents/boreogadus
ANGSDDIR=/home/mchale/Desktop/angsd/
OUTDIR=/home/mchale/Documents/boreogadus/fst
GLSDIR=/home/mchale/Documents/boreogadus/gls

#Filters:
FILTERS="-minQ 30 -minMapQ 30 -C 50 -remove_bads 1 -uniqueOnly 1 -skipTriallelic 1" # got rid of minMaf and snp-pval --> using filtered sites flag
#Options:
OPTIONS="-GL 1 -doSaf 1 -doMaf 2 -doMajorMinor 1 -doCounts 1 -doPost 1 -doGlf 2"

cd ${ANGSDDIR}
mindepth=$((20))
maxdepth=$((200))

    ./angsd -nthreads 12 \
	-bam ${BAMDIR}/boreogadus_Iceland_bams.txt \
	-anc ${REFGENOME} -ref ${REFGENOME} \
	-out ${GLSDIR}/boreogadus_Iceland_nosites_20241212 \
	${FILTERS} ${OPTIONS} -snp_pval 1e-6 \
    	-setminDepth 20 -setmaxDepth 200 -minInd 7 #chose minInd7 to do fst for Labrador, too


    ./angsd -nthreads 12 \
	-bam ${BAMDIR}/boreogadus_Iceland_bams_nodownsample.txt \
        -anc ${REFGENOME} -ref ${REFGENOME} \
        -out ${GLSDIR}/boreogadus_Iceland_nosites_nodownsample_20241212 \
        ${FILTERS} ${OPTIONS} -snp_pval 1e-6 \
        -setminDepth 36 -setmaxDepth 360 -minInd 7 #chose minInd7 to do fst for Labrador, too

    ./angsd -nthreads 12 \
        -bam ${BAMDIR}/boreogadus_Iceland_bams_nodownsample.txt \
        -anc ${REFGENOME} -ref ${REFGENOME} \
	-sites ${GLSDIR}/boreogadus_filtered_SNPs_filetwo.txt \
        -out ${GLSDIR}/boreogadus_Iceland_nodownsample_20241212 \
        ${FILTERS} ${OPTIONS} \
        -setminDepth 36 -setmaxDepth 360 -minInd 7 #chose minInd7 to do fst for Labrador, too

#realSFS

./misc/realSFS -P 8 ${GLSDIR}/boreogadus_Iceland_nosites_20241212.saf.idx ${GLSDIR}/boreogadus_Chukchi_20241210.saf.idx -fold 1 > ${OUTDIR}/boreogadus_Iceland-Chukchi_nosites_20241212.sfs
./misc/realSFS fst index ${GLSDIR}/boreogadus_Iceland_nosites_20241212.saf.idx -fold 1 ${GLSDIR}/boreogadus_Chukchi_20241210.saf.idx -fold 1 -sfs ${OUTDIR}/boreogadus_Iceland-Chukchi_nosites_20241212.sfs -fstout ${OUTDIR}/boreogadus_Iceland-Chukchi_nosites_20241212.sfs
./misc/realSFS fst stats2 ${OUTDIR}/boreogadus_Iceland-Chukchi_nosites_20241212.sfs.fst.idx -win 1 -step 1 > ${OUTDIR}/boreogadus_Iceland-Chukchi_nosites_20241212.fst.SNP.txt

./misc/realSFS -P 8 ${GLSDIR}/boreogadus_Iceland_nodownsample_20241212.saf.idx ${GLSDIR}/boreogadus_Chukchi_20241210.saf.idx -fold 1 > ${OUTDIR}/boreogadus_Iceland-Chukchi_nodownsample_20241212.sfs
./misc/realSFS fst index ${GLSDIR}/boreogadus_Iceland_nodownsample_20241212.saf.idx -fold 1 ${GLSDIR}/boreogadus_Chukchi_20241210.saf.idx -fold 1 -sfs ${OUTDIR}/boreogadus_Iceland-Chukchi_nodownsample_20241212.sfs -fstout ${OUTDIR}/boreogadus_Iceland-Chukchi_nodownsample_20241212.sfs
./misc/realSFS fst stats2 ${OUTDIR}/boreogadus_Iceland-Chukchi_nodownsample_20241212.sfs.fst.idx -win 1 -step 1 > ${OUTDIR}/boreogadus_Iceland-Chukchi_nodownsample_20241212.fst.SNP.txt

./misc/realSFS -P 8 ${GLSDIR}/boreogadus_Iceland_nosites_nosites_nodownsample_20241212.saf.idx ${GLSDIR}/boreogadus_Chukchi_20241210.saf.idx -fold 1 > ${OUTDIR}/boreogadus_Iceland-Chukchi_nosites_nodownsample_20241212.sfs
./misc/realSFS fst index ${GLSDIR}/boreogadus_Iceland_nosites_nodownsample_20241212.saf.idx -fold 1 ${GLSDIR}/boreogadus_Chukchi_20241210.saf.idx -fold 1 -sfs ${OUTDIR}/boreogadus_Iceland-Chukchi_nosites_nodownsample_20241212.sfs -fstout ${OUTDIR}/boreogadus_Iceland-Chukchi_nosites_nodownsample_20241212.sfs
./misc/realSFS fst stats2 ${OUTDIR}/boreogadus_Iceland-Chukchi_nosites_nodownsample_20241212.sfs.fst.idx -win 1 -step 1 > ${OUTDIR}/boreogadus_Iceland-Chukchi_nosites_nodownsample_20241212.fst.SNP.txt
