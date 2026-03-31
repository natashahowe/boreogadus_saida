#!/bin/bash

#SBATCH --cpus-per-task=4
#SBATCH --time=0-20:00:00
#SBATCH --job-name=bsaida-admix
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=natasha.howe@noaa.gov
#SBATCH --output=/home/nhowe/arctic/job_outfiles/boreogadus_unfused_admix_%A-%a.out
#SBATCH --partition=standard,medmem
#SBATCH --array=2-6%6

module unload bio/ngsadmix
module load bio/ngsadmix

for k_val in {2..6}
do
	if [[ ${SLURM_ARRAY_TASK_ID} == ${k_val} ]]; then
		break
	fi
done

NGSadmix -likes /home/nhowe/arctic/gls/boreogadus_unfused_neutral_alpha0.05_e4_rmPCAngsdOutliers_rmInversions_rmchr9a.beagle.gz -K ${k_val} -outfiles /home/nhowe/arctic/admixture/boreogadus_unfused_neutral_alpha0.05_e4_rmPCAngsdOutliers_rmInversions_rmchr9a_k${k_val}-0 -P 4 -minMaf 0
NGSadmix -likes /home/nhowe/arctic/gls/boreogadus_unfused_neutral_alpha0.05_e4_rmPCAngsdOutliers_rmInversions_rmchr9a.beagle.gz -K ${k_val} -outfiles /home/nhowe/arctic/admixture/boreogadus_unfused_neutral_alpha0.05_e4_rmPCAngsdOutliers_rmInversions_rmchr9a_k${k_val}-1 -P 4 -minMaf 0
NGSadmix -likes /home/nhowe/arctic/gls/boreogadus_unfused_neutral_alpha0.05_e4_rmPCAngsdOutliers_rmInversions_rmchr9a.beagle.gz -K ${k_val} -outfiles /home/nhowe/arctic/admixture/boreogadus_unfused_neutral_alpha0.05_e4_rmPCAngsdOutliers_rmInversions_rmchr9a_k${k_val}-2 -P 4 -minMaf 0

#NGSadmix -likes /home/nhowe/arctic/gls/boreogadus_unfused_pruned_loci.beagle.gz -K ${k_val} -outfiles /home/nhowe/arctic/admixture/boreogadus_unfused_pruned_loci_k${k_val}-0 -P 4 -minMaf 0
#NGSadmix -likes /home/nhowe/arctic/gls/boreogadus_unfused_pruned_loci.beagle.gz -K ${k_val} -outfiles /home/nhowe/arctic/admixture/boreogadus_unfused_pruned_loci_k${k_val}-1 -P 4 -minMaf 0
#NGSadmix -likes /home/nhowe/arctic/gls/boreogadus_unfused_pruned_loci.beagle.gz -K ${k_val} -outfiles /home/nhowe/arctic/admixture/boreogadus_unfused_pruned_loci_k${k_val}-2 -P 4 -minMaf 0

