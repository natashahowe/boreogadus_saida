#!/bin/bash

#SBATCH --cpus-per-task=5
#SBATCH --time=0-20:00:00
#SBATCH --job-name=ld_bs
#SBATCH --partition=medmem,himem
#SBATCH --mem=75GB
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=natasha.howe@noaa.gov
#SBATCH --output=/home/nhowe/arctic/job_outfiles/boreogadus_plink_ld_and_plot_%A-%a.out
#SBATCH --array=1

module purge
module load bio/plink/1.90b6.23

#module load R/4.4.1
# I did something to the R environment when trying to install packages and it messed something up.
# --> have to use R through a mamba environment now...

# Naming
DATADIR='/home/nhowe/arctic'
TFILE="boreogadus_wholegenome_20241222"
ARRAYFILE="/home/nhowe/arctic/scripts/chr1-5_peak_input.txt"

CHROMS=$(< "$ARRAYFILE")

for peak in ${CHROMS}
do
        job_index=$(echo ${peak} | awk -F ":" '{print $1}')
        chrom=$(echo ${peak} | awk -F ":" '{print $2}')
        if [[ ${SLURM_ARRAY_TASK_ID} == ${job_index} ]]; then
                break
        fi
done

FILENAME='boreogadus_'${chrom}'_plink_2minr'

echo "-- Plink run ld"
plink --tfile ${DATADIR}'/gls/'${TFILE} --r2 inter-chr --allow-extra-chr \
	-chr ${chrom} \
	--out ${DATADIR}'/peaks/'${FILENAME} \
	--maf 0.05 --ld-window-r2 0.2 #\
	#--remove ${DATADIR}'/scripts/boreogadus_rmIceland_plink.txt'
echo "-- Plink done ld"

# remove column titles and replace spaces from plinks output with tabs, then only keep snp1, snp2, and r2
tail -n +2 ${DATADIR}'/peaks/'${FILENAME}'.ld' |  tr -s ' ' '\t' | cut -f2,5,7 > ${DATADIR}'/peaks/'${FILENAME}'_r2.ld'

# --- Mamba Environment Activation ---------------------------------------------

# The full path to mamba/conda installation.
CONDA_BASE=$(conda info --base)

# The path to your specific mamba environment
MAMBA_ENV_NAME="rockfish"

# Source the mamba initialization script to make 'conda' and 'mamba' commands available
source "${CONDA_BASE}/etc/profile.d/conda.sh"

# Activate the mamba environment
conda activate "${MAMBA_ENV_NAME}"

# ---- Reset the R Library Path --------------------------------------------------
# use R library set in the rockfish environment (with appropriate libraries installed)
export R_LIBS_USER="/home/nhowe/.conda/envs/rockfish/lib/R/library"

echo "Run ld plot"
Rscript --vanilla ${DATADIR}/scripts/ldplot_plink.R ${DATADIR}'/peaks' ${FILENAME} /figures/
echo " --- Done"
