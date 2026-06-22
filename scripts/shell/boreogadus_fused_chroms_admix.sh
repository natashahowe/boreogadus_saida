#!/bin/bash

#SBATCH --cpus-per-task=10
#SBATCH --time=0-20:00:00
#SBATCH --job-name=bsaida-admix
#SBATCH --partition=medmem
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=natasha.howe@noaa.gov
#SBATCH --output=/home/nhowe/arctic/job_outfiles/boreogadus_fused_admix_%A-%a.out
#SBATCH --array=2-6%6

module unload bio/ngsadmix
module load bio/ngsadmix

for k_val in {2..6}
do
	if [[ ${SLURM_ARRAY_TASK_ID} == ${k_val} ]]; then
		break
	fi
done

#NGSadmix -likes /home/nhowe/arctic/gls/boreogadus_fused_chroms_20241222.beagle.gz -K ${k_val} -outfiles /home/nhowe/arctic/admixture/boreogadus_fused_chroms_k${k_val}-1 -P 5 -minMaf 0
#NGSadmix -likes /home/nhowe/arctic/gls/boreogadus_fused_chroms_20241222.beagle.gz -K ${k_val} -outfiles /home/nhowe/arctic/admixture/boreogadus_fused_chroms_k${k_val}-2 -P 5 -minMaf 0
#NGSadmix -likes /home/nhowe/arctic/gls/boreogadus_fused_chroms_20241222.beagle.gz -K ${k_val} -outfiles /home/nhowe/arctic/admixture/boreogadus_fused_chroms_k${k_val}-0 -P 5 -minMaf 0

#NGSadmix -likes /home/nhowe/arctic/gls/boreogadus_fused_adaptive_alpha0.05.beagle.gz -K ${k_val} -outfiles /home/nhowe/arctic/admixture/boreogadus_fused_adaptive_alpha0.05_k${k_val}-1 -P 10 -minMaf 0
#NGSadmix -likes /home/nhowe/arctic/gls/boreogadus_fused_adaptive_alpha0.05.beagle.gz -K ${k_val} -outfiles /home/nhowe/arctic/admixture/boreogadus_fused_adaptive_alpha0.05_k${k_val}-2 -P 10 -minMaf 0
#NGSadmix -likes /home/nhowe/arctic/gls/boreogadus_fused_adaptive_alpha0.05.beagle.gz -K ${k_val} -outfiles /home/nhowe/arctic/admixture/boreogadus_fused_adaptive_alpha0.05_k${k_val}-0 -P 10 -minMaf 0
#NGSadmix -likes /home/nhowe/arctic/gls/boreogadus_fused_adaptive_alpha0.05.beagle.gz -K ${k_val} -outfiles /home/nhowe/arctic/admixture/boreogadus_fused_adaptive_alpha0.05_k${k_val}-3 -P 10 -minMaf 0
#NGSadmix -likes /home/nhowe/arctic/gls/boreogadus_fused_adaptive_alpha0.05.beagle.gz -K ${k_val} -outfiles /home/nhowe/arctic/admixture/boreogadus_fused_adaptive_alpha0.05_k${k_val}-4 -P 10 -minMaf 0
#NGSadmix -likes /home/nhowe/arctic/gls/boreogadus_fused_adaptive_alpha0.05.beagle.gz -K ${k_val} -outfiles /home/nhowe/arctic/admixture/boreogadus_fused_adaptive_alpha0.05_k${k_val}-5 -P 10 -minMaf 0

NGSadmix -likes /home/nhowe/arctic/gls/boreogadus_fused_adaptive_alpha0.05.beagle.gz -K ${k_val} -outfiles /home/nhowe/arctic/admixture/boreogadus_fused_adaptive_alpha0.05_k${k_val}-6 -P 10 -minMaf 0
NGSadmix -likes /home/nhowe/arctic/gls/boreogadus_fused_adaptive_alpha0.05.beagle.gz -K ${k_val} -outfiles /home/nhowe/arctic/admixture/boreogadus_fused_adaptive_alpha0.05_k${k_val}-7 -P 10 -minMaf 0
NGSadmix -likes /home/nhowe/arctic/gls/boreogadus_fused_adaptive_alpha0.05.beagle.gz -K ${k_val} -outfiles /home/nhowe/arctic/admixture/boreogadus_fused_adaptive_alpha0.05_k${k_val}-8 -P 10 -minMaf 0
NGSadmix -likes /home/nhowe/arctic/gls/boreogadus_fused_adaptive_alpha0.05.beagle.gz -K ${k_val} -outfiles /home/nhowe/arctic/admixture/boreogadus_fused_adaptive_alpha0.05_k${k_val}-9 -P 10 -minMaf 0
