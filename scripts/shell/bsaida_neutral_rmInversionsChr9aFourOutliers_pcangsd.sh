#!/bin/bash

#SBATCH --cpus-per-task=8
#SBATCH --time=0-05:00:00
#SBATCH --job-name=boreogadus_neutral
#SBATCH --mail-type=FAIL
#SBATCH --partition=standard,medmem
#SBATCH --mail-user=natasha.howe@noaa.gov
#SBATCH --output=/home/nhowe/arctic/job_outfiles/boreogadus_neutral_pca_rmInv_%A.out

# Ran after bsaida_neutral_rmInversions_pluschr9a_pcangsd.sh

module purge
module load bio/pcangsd/1.36.4

PROJDIR="/home/nhowe/arctic"
nEigen=10

NEWSUFFIX='unfused_neutral_alpha0.05_e4_rmPCAngsdOutliers_rmInversions'

## ------ Run PCANGSD ------ ##
echo "start neutral dataset pcangsd"

pcangsd --threads 8 -b ${PROJDIR}/gls/boreogadus_${NEWSUFFIX}_rmchr9a.beagle.gz -o ${PROJDIR}/pca/boreogadus_${NEWSUFFIX}_rmchr9a_rmOutliers --iter 500 --filter boreogadus_neutral_pca_four_outliers.txt

echo "end neutral dataset pcangsd"

## ------ Run PCANGSD: Loadings ------ ##

echo "Create pca, selection, & loading weights file with ${nEigen} pcs"
# pcangsd --threads 8 -b ${PROJDIR}/gls/boreogadus_${NEWSUFFIX}_rmchr9a.beagle.gz -o ${PROJDIR}/pca/boreogadus_${NEWSUFFIX}_rmchr9a_e${nEigen} --sites-save --eig ${nEigen} --selection --snp-weights --iter 500 --filter boreogadus_neutral_pca_four_outliers.txt

echo "end neutral dataset pcangsd loadings/selection"

echo "-- End of Script --"
