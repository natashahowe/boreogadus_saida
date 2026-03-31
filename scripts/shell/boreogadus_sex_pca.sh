#!/bin/bash

#SBATCH --job-name=sexPCA
#SBATCH --time=0-1:00:00
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=natasha.howe@noaa.gov
#SBATCH --output=/home/nhowe/arctic/job_outfiles/sex_pca.%j.out 

module unload bio/pcangsd/0.99
module load bio/pcangsd/0.99
source /opt/bioinformatics/venv/pcangsd-0.99/bin/activate

path=/home/nhowe/arctic
prefix=boreogadus
chrom=OZ177908.1
firstpos=26000000
lastpos=27500000

# call in beagle file for pop and chromosome
BEAGLE=${path}/gls/${prefix}_${chrom}_20241222.beagle.gz

# cut 9a and 9b for the regions of interest
zcat ${BEAGLE} | awk -v s=$firstpos -v e=$lastpos -F'[\t_]' '$2 >= s && $2 <= e' | gzip > ${path}/gls/${prefix}_${chrom}_s${firstpos}_e${lastpos}.beagle.gz

pcangsd.py -threads 10 -beagle ${path}/gls/${prefix}_${chrom}_s${firstpos}_e${lastpos}.beagle.gz -o ${path}/pca/${prefix}_${chrom}_s${firstpos}_e${lastpos}
