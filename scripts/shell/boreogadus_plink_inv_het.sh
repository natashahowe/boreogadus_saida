#!/bin/bash

#SBATCH --cpus-per-task=1
#SBATCH --time=0-20:00:00
#SBATCH --job-name=bs_het
#SBATCH --partition=standard,medmem
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=natasha.howe@noaa.gov
#SBATCH --output=/home/nhowe/arctic/job_outfiles/boreogadus_plink_het_%A-%a.out
#SBATCH --array=1-23%12

module purge
module load bio/plink/1.90b6.23

#module load R/4.4.1
# I did something to the R environment when trying to install packages and it messed something up.
# --> have to use R through a mamba environment now...

# Naming
DATADIR='/home/nhowe/arctic'
TFILE="boreogadus_wholegenome_20241222"
ARRAYFILE="/home/nhowe/arctic/scripts/bsaida_inversions_input.txt"

INVERSIONS=$(< "$ARRAYFILE")

mkdir -p ${DATADIR}'/peaks/het'

for peak in ${INVERSIONS}
do
        job_index=$(echo ${peak} | awk -F ":" '{print $1}')
        chrom=$(echo ${peak} | awk -F ":" '{print $2}')
	startpos=$(echo ${peak} | awk -F ":" '{print $3}')
	endpos=$(echo ${peak} | awk -F ":" '{print $4}')
        if [[ ${SLURM_ARRAY_TASK_ID} == ${job_index} ]]; then
                break
        fi
done

FILENAME=boreogadus_${chrom}_s${startpos}_e${endpos}

echo $FILENAME

echo "-- Plink peak bed --> run"
plink --tfile ${DATADIR}/gls/${TFILE} --allow-extra-chr \
        -chr ${chrom} --from-bp ${startpos} --to-bp ${endpos} \
        --make-bed \
	--out ${DATADIR}/peaks/het/${FILENAME}

echo "-- Plink peak bed ----> done "


echo "-- Plink het --> run"
plink --bfile ${DATADIR}/peaks/het/${FILENAME} --allow-extra-chr \
	--het --out ${DATADIR}/peaks/het/${FILENAME}

echo "-- Plink het ----> done "

