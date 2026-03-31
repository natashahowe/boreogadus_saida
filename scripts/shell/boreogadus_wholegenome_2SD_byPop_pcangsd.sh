#!/bin/bash

#SBATCH --cpus-per-task=10
#SBATCH --time=0-05:00:00
#SBATCH --job-name=bsaida-pca
#SBATCH --mail-type=FAIL,END
#SBATCH --partition=medmem
#SBATCH --mail-user=natasha.howe@noaa.gov
#SBATCH --output=/home/nhowe/arctic/job_outfiles/boreogadus_wholegenome_adaptive-neutral_20241222.%A.out

module unload bio/pcangsd/0.99
module load bio/pcangsd/0.99
source /opt/bioinformatics/venv/pcangsd-0.99/bin/activate

POP='Pond'

zcat /home/nhowe/arctic/gls/boreogadus_wholegenome_20241222.beagle.gz | grep -Fwf /home/nhowe/arctic/fst/Pond_neutral_markerPos_alpha0.05.txt | gzip > /home/nhowe/arctic/gls/${POP}_boreogadus_wholegenome_20241222_neutral_alpha0.05.beagle.gz

zcat /home/nhowe/arctic/gls/boreogadus_wholegenome_20241222.beagle.gz | grep -Fwf /home/nhowe/arctic/fst/Pond_neutral_markerPos_alpha0.001.txt | gzip > /home/nhowe/arctic/gls/${POP}_boreogadus_wholegenome_20241222_neutral_alpha0.001.beagle.gz

zcat /home/nhowe/arctic/gls/boreogadus_wholegenome_20241222.beagle.gz | grep -Fwf /home/nhowe/arctic/fst/Pond_adaptive_markerPos_alpha0.05.txt | gzip > /home/nhowe/arctic/gls/${POP}_boreogadus_wholegenome_20241222_adaptive_alpha0.05.beagle.gz

zcat /home/nhowe/arctic/gls/boreogadus_wholegenome_20241222.beagle.gz | grep -Fwf /home/nhowe/arctic/fst/Pond_adaptive_markerPos_alpha0.001.txt | gzip > /home/nhowe/arctic/gls/${POP}_boreogadus_wholegenome_20241222_adaptive_alpha0.001.beagle.gz

pcangsd.py -threads 10  -beagle /home/nhowe/arctic/gls/${POP}_boreogadus_wholegenome_20241222_neutral_alpha0.05.beagle.gz -o /home/nhowe/arctic/pca/${POP}_boreogadus_wholegenome_20241222_neutral_alpha0.05

pcangsd.py -threads 10  -beagle /home/nhowe/arctic/gls/${POP}_boreogadus_wholegenome_20241222_neutral_alpha0.001.beagle.gz -o /home/nhowe/arctic/pca/${POP}_boreogadus_wholegenome_20241222_neutral_alpha0.001

pcangsd.py -threads 10  -beagle /home/nhowe/arctic/gls/${POP}_boreogadus_wholegenome_20241222_adaptive_alpha0.05.beagle.gz -o /home/nhowe/arctic/pca/${POP}_boreogadus_wholegenome_20241222_adaptive_alpha0.05

pcangsd.py -threads 10  -beagle /home/nhowe/arctic/gls/${POP}_boreogadus_wholegenome_20241222_adaptive_alpha0.001.beagle.gz -o /home/nhowe/arctic/pca/${POP}_boreogadus_wholegenome_20241222_adaptive_alpha0.001

