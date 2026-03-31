#!/bin/bash

#SBATCH --cpus-per-task=10
#SBATCH --time=0-05:00:00
#SBATCH --job-name=bsaida-pca
#SBATCH --mail-type=FAIL,END
#SBATCH --partition=standard
#SBATCH --mail-user=natasha.howe@noaa.gov
#SBATCH --output=/home/nhowe/arctic/job_outfiles/boreogadus_chroms_pca.%A-%a.out
#SBATCH --array=1-18%18

module unload bio/pcangsd/0.99
module load bio/pcangsd/0.99
source /opt/bioinformatics/venv/pcangsd-0.99/bin/activate

JOBS_FILE=/home/nhowe/arctic/scripts/boreogadus_chrom_array.txt
IDS=$(cat ${JOBS_FILE})

for sample_line in ${IDS}
do
      job_index=$(echo ${sample_line} | awk -F ":" '{print $1}')
      beagle_file=$(echo ${sample_line} | awk -F ":" '{print $2}')
      if [[ ${SLURM_ARRAY_TASK_ID} == ${job_index} ]]; then
      break
      fi
done

chrom=$(echo $beagle_file | sed 's!^.*/!!')
chrom=${chrom%.beagle.gz}

zcat /home/nhowe/arctic/gls/boreogadus_wholegenome_20241222.beagle.gz | awk -F'[_\t]' -v c=${chrom} '$1 == c' | gzip > /home/nhowe/arctic/gls/boreogadus_${chrom}_20241222.beagle.gz

# Run PCA
pcangsd.py -threads 10 -beagle /home/nhowe/arctic/gls/boreogadus_${chrom}_20241222.beagle.gz -o /home/nhowe/arctic/pca/boreogadus_${chrom}_20241222

