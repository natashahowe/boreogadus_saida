#!/bin/bash

#SBATCH --cpus-per-task=10
#SBATCH --time=0-05:00:00
#SBATCH --job-name=bsaida-pca
#SBATCH --mail-type=FAIL,END
#SBATCH --partition=standard,medmem
#SBATCH --mail-user=natasha.howe@noaa.gov
#SBATCH --output=/home/nhowe/arctic/job_outfiles/boreogadus_fusion_adaptive-neutral_20241222.%A.out

module unload bio/pcangsd/0.99
module load bio/pcangsd/0.99
source /opt/bioinformatics/venv/pcangsd-0.99/bin/activate

zcat /home/nhowe/arctic/gls/boreogadus_unfused_chroms_20241222.beagle.gz | grep -Fwf /home/nhowe/arctic/fst/boreogadus_unfused_neutral_markerPos_alpha0.05.txt | gzip > /home/nhowe/arctic/gls/boreogadus_unfused_neutral_alpha0.05.beagle.gz

zcat /home/nhowe/arctic/gls/boreogadus_fused_chroms_20241222.beagle.gz | grep -Fwf /home/nhowe/arctic/fst/boreogadus_fused_adaptive_markerPos_alpha0.05.txt | gzip > /home/nhowe/arctic/gls/boreogadus_fused_adaptive_alpha0.05.beagle.gz

pcangsd.py -threads 10 -beagle /home/nhowe/arctic/gls/boreogadus_unfused_chroms_20241222.beagle.gz -o /home/nhowe/arctic/pca/boreogadus_unfused_chroms_20241222
pcangsd.py -threads 10 -beagle /home/nhowe/arctic/gls/boreogadus_unfused_neutral_alpha0.05.beagle.gz -o /home/nhowe/arctic/pca/boreogadus_unfused_neutral_alpha0.05

pcangsd.py -threads 10 -beagle /home/nhowe/arctic/gls/boreogadus_fused_chroms_20241222.beagle.gz -o /home/nhowe/arctic/pca/boreogadus_fused_chroms_20241222
pcangsd.py -threads 10 -beagle /home/nhowe/arctic/gls/boreogadus_fused_adaptive_alpha0.05.beagle.gz -o /home/nhowe/arctic/pca/boreogadus_fused_adaptive_alpha0.05
