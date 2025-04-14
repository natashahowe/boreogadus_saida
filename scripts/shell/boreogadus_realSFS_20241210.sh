#!/bin/bash

# Natasha Howe

# realSFS steps for Fst
## Specify directories with data
OUTDIR=/home/mchale/Documents/boreogadus/fst
GLSDIR=/home/mchale/Documents/boreogadus/gls
ANGSDDIR=/home/mchale/Desktop/angsd/
FST_FILE=/home/mchale/Documents/boreogadus/scripts/boreogadus_fst_comparisons.txt
#FST_FILE=/home/mchale/Documents/boreogadus/scripts/boreogadus_fst_subset_comparisons.txt

cd ${ANGSDDIR}

while IFS=$'\t' read POP1 POP2 REST;
do echo $POP1 $POP2;
./misc/realSFS -P 8 ${GLSDIR}/boreogadus_${POP1}_20241210.saf.idx -fold 1 ${GLSDIR}/boreogadus_${POP2}_20241210.saf.idx -fold 1 > ${OUTDIR}/boreogadus_${POP1}-${POP2}_20241210.sfs
./misc/realSFS fst index ${GLSDIR}/boreogadus_${POP1}_20241210.saf.idx -fold 1 ${GLSDIR}/boreogadus_${POP2}_20241210.saf.idx -fold 1 -sfs ${OUTDIR}/boreogadus_${POP1}-${POP2}_20241210.sfs -fstout ${OUTDIR}/boreogadus_${POP1}-${POP2}_20241210.sfs
./misc/realSFS fst stats ${OUTDIR}/boreogadus_${POP1}-${POP2}_20241210.sfs.fst.idx > ${OUTDIR}/boreogadus_${POP1}-${POP2}_20241210.sfs.global.fst
./misc/realSFS fst stats2 ${OUTDIR}/boreogadus_${POP1}-${POP2}_20241210.sfs.fst.idx -win 1 -step 1 > ${OUTDIR}/boreogadus_${POP1}-${POP2}_20241210.fst.SNP.txt
done < $FST_FILE
