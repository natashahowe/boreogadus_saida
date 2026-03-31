#!/bin/bash

#SBATCH --cpus-per-task=5
#SBATCH --time=0-20:00:00
#SBATCH --job-name=bsaida-admix
#SBATCH --partition=standard
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=natasha.howe@noaa.gov
#SBATCH --output=/home/nhowe/arctic/job_outfiles/boreogadus_unfused_admix_%A-%a.out
#SBATCH --array=2-6%6

module unload bio/ngsadmix
module load bio/ngsadmix

for k_val in {2..6}
do
	if [[ ${SLURM_ARRAY_TASK_ID} == ${k_val} ]]; then
		break
	fi
done

NGSadmix -likes /home/nhowe/arctic/gls/boreogadus_unfused_chroms_20241222.beagle.gz -K ${k_val} -outfiles /home/nhowe/arctic/admixture/boreogadus_unfused_chroms_k${k_val}-1 -P 5 -minMaf 0
NGSadmix -likes /home/nhowe/arctic/gls/boreogadus_unfused_chroms_20241222.beagle.gz -K ${k_val} -outfiles /home/nhowe/arctic/admixture/boreogadus_unfused_chroms_k${k_val}-2 -P 5 -minMaf 0
NGSadmix -likes /home/nhowe/arctic/gls/boreogadus_unfused_chroms_20241222.beagle.gz -K ${k_val} -outfiles /home/nhowe/arctic/admixture/boreogadus_unfused_chroms_k${k_val}-0 -P 5 -minMaf 0

