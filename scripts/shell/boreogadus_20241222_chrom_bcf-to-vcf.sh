#!/bin/bash

#SBATCH --cpus-per-task=10
#SBATCH --time=0-24:00:00
#SBATCH --job-name=bsaida-bcf
#SBATCH --mail-type=FAIL,END
#SBATCH --partition=standard,medmem
#SBATCH --mail-user=natasha.howe@noaa.gov
#SBATCH --output=/home/nhowe/arctic/job_outfiles/boreogadus_chrom_bcf-to-vcf.%A-%a.out
#SBATCH --array=1-18%9

module unload bio/bcftools/1.11
module load bio/bcftools/1.11

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

echo 'bcf to vcf'
bcftools view /home/nhowe/arctic/gls/boreogadus_${chrom}_allpops_20241222_bcf.bcf > /home/nhowe/arctic/gls/boreogadus_${chrom}_allpops_20241222_bcf.vcf

echo 'now bgzip'
bgzip /home/nhowe/arctic/gls/boreogadus_${chrom}_allpops_20241222_bcf.vcf

echo 'now index'
bcftools index /home/nhowe/arctic/gls/boreogadus_${chrom}_allpops_20241222_bcf.vcf.gz

