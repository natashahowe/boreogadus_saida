#!/bin/bash

# Natasha Howe

# PCAngsd
# Using beagle file from something Matt already ran

## Specify directories with data
PCADIR=/home/mchale/Documents/boreogadus/pca
BEAGLE=/home/mchale/Desktop/angsd/Cod_filteredSNPs.beagle.gz # Matt ran this with a -sites flag
PCANGSD=/home/mchale/programsNH/pcangsd/pcangsd

cd ${PCANGSD}

#pcangsd -b ~/Desktop/angsd/Cod_filteredSNPs.beagle.gz -o ~/Documents/boreogadus/pca/boreogadus_wholegenome_20241210

#pcangsd -b ${BEAGLE} -o ${PCADIR}/boreogadus_wholegenome_20241210

# filter out Chukchi
#pcangsd --threads 4 --beagle ${BEAGLE} --out ${PCADIR}/boreogadus_wholegenome_noChukchi_20241210 --filter ${PCADIR}/noChukchi_pcafilter.txt

# filter out Chukchi & Coronation
pcangsd --threads 4 --beagle ${BEAGLE} --out ${PCADIR}/boreogadus_wholegenome_noChukchiCoronation_20241210 --filter ${PCADIR}/noChukchiCoronation_pcafilter.txt

# filter out Chukchi and Iceland
pcangsd --threads 4 --beagle ${BEAGLE} --out ${PCADIR}/boreogadus_wholegenome_noChukchiIceland_20241210 --filter ${PCADIR}/noChukchiIceland_pcafilter.txt

# filter out Iceland
pcangsd --threads 4 --beagle ${BEAGLE} --out ${PCADIR}/boreogadus_wholegenome_noIceland_20241210 --filter ${PCADIR}/noIceland_pcafilter.txt
