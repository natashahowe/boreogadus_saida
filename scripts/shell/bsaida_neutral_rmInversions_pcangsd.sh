#!/bin/bash

#SBATCH --cpus-per-task=8
#SBATCH --time=0-05:00:00
#SBATCH --job-name=boreogadus_neutral
#SBATCH --mail-type=FAIL
#SBATCH --partition=standard,medmem
#SBATCH --mail-user=natasha.howe@noaa.gov
#SBATCH --output=/home/nhowe/arctic/job_outfiles/boreogadus_neutral_pca_rmInv_%A.out

module purge
module load bio/pcangsd/1.36.4

PROJDIR="/home/nhowe/arctic"
SUFFIX="unfused_neutral_alpha0.05"
nEigen=4

NEWSUFFIX='unfused_neutral_alpha0.05_e4_rmPCAngsdOutliers_rmInversions'

NEWSNPS=${PROJDIR}/scripts/boreogadus_${NEWSUFFIX}.tsv

# remove beagle file if already exists
rm ${PROJDIR}/gls/boreogadus_${NEWSUFFIX}.beagle.gz

zcat ${PROJDIR}/gls/boreogadus_${SUFFIX}.beagle.gz | grep -Fwf ${NEWSNPS} | gzip > ${PROJDIR}/gls/boreogadus_${NEWSUFFIX}.beagle.gz


## ------ Run PCANGSD ------ ##
echo "start neutral dataset pcangsd"

pcangsd --threads 8 -b ${PROJDIR}/gls/boreogadus_${NEWSUFFIX}.beagle.gz -o ${PROJDIR}/pca/boreogadus_${NEWSUFFIX} --iter 500

#echo "end neutral dataset pcangsd"

## ------ Run PCANGSD: Loadings ------ ##

echo "Create pca, selection, & loading weights file with ${nEigen} pcs"
pcangsd --threads 8 -b ${PROJDIR}/gls/boreogadus_${NEWSUFFIX}.beagle.gz -o ${PROJDIR}/pca/boreogadus_${NEWSUFFIX}_e${nEigen} --sites-save --eig ${nEigen} --selection --snp-weights --iter 500

# pcangsd --threads 8 -b ${PROJDIR}/gls/boreogadus_${NEWSUFFIX}.beagle.gz -o ${PROJDIR}/pca/boreogadus_${NEWSUFFIX}_pcadapt --sites-save --pcadapt --iter 500
echo "end neutral dataset pcangsd loadings/selection"

echo "-- End of Script --"
