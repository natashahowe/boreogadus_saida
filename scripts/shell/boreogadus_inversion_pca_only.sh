#!/bin/bash
#SBATCH --job-name=pca_inv
#SBATCH --time=1-12:00:00
#SBATCH -c 5
#SBATCH --partition=standard,medmem
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=natasha.howe@noaa.gov
#SBATCH --output=/home/nhowe/arctic/job_outfiles/boreogadus_localpca_202603.%A-%a.out
#SBATCH --array=1-23%12

# HOW TO RUN THE FOLLOWING
# NO SUBSETTING EXAMPLE: sbatch scripts/ld_pca_analysis_maturity_v5.sh -d /home/usr/chinook -k scripts/peak_analysis_test_repeat_peak.txt -b pollREF_wholegenome_polymorphic.beagle -f /figures/ -p pollREF -m pollREF_filtered_bamslist.txt -l no 

# SUBSETTING EXAMPLE: sbatch scripts/ld_pca_analysis_maturity_v5.sh -d /home/usr/chinook -k scripts/peak_analysis_test_repeat_peak.txt -b pollREF_wholegenome_polymorphic.beagle -f /figures/ -p pollREF -m pollREF_filtered_bamslist.txt -l no 

while getopts "d:k:b:f:p:m:l:e:" opt;
do
	case $opt in
		d) DATADIR=${OPTARG};;
		k) PEAKSFILE=${OPTARG};;
		b) BEAGLE=${OPTARG};;
		f) FIGOUT=${OPTARG};;
		p) PREFIX=${OPTARG};;
		m) BAMLIST=${OPTARG};;
		l) LDSUBSET=${OPTARG};;
		e) LDEVERY=${OPTARG};;
	esac
done

PEAKS=$(cat ${DATADIR}/${PEAKSFILE})

## peaksARRAY_input.txt will have XXX sets of data separated by colons: the job index, the chromosome of the peak of interest, the start postition, the end position, and the name of peak

for peak in ${PEAKS}
do
	job_index=$(echo ${peak} | awk -F ":" '{print $1}')
	chrom=$(echo ${peak} | awk -F ":" '{print $2}')
	s=$(echo ${peak} | awk -F ":" '{print $3}')
	e=$(echo ${peak} | awk -F ":" '{print $4}')
	if [[ ${SLURM_ARRAY_TASK_ID} == ${job_index} ]]; then
		break
	fi
done

##########################################
## Manipulate data for downstream analyses

# check to see if your beagle file is compressed and if not zip it
if [[ -f ${DATADIR}'/gls/'${BEAGLE}'.gz' ]]
then
	echo "Whole genome beagle already zipped."
elif [[ ${BEAGLE: -3} == '.gz' ]]
then
	BEAGLE=${BEAGLE%.*}
	echo "Whole genome beagle already zipped."
else
	gzip ${DATADIR}'/gls/'${BEAGLE}'.gz'
fi

# check to see if your output folders exist. if not, make them
if [[ -d ${DATADIR}/peaks ]]
then
	echo "Peaks directory already present."
else
	mkdir ${DATADIR}/peaks
fi

if [[ -d ${DATADIR}/peaks/figures ]]
then
	echo "Figures directory already present."
else
	mkdir ${DATADIR}/peaks/figures
fi

FILENAME=${PREFIX}_${chrom}_s${s}_e${e}

echo 'chrom:' $chrom
echo 'startpos:' $s
echo 'endpos:' $e
echo 'filename: ' $FILENAME

# check to see if your subset beagle file exists. If not, create it
if [[ -f ${DATADIR}'/peaks/'${FILENAME}'.beagle.gz' ]]
then
	echo "Beagle for peak already exists."
	echo "If rerunning with the same naming convention, remove files from" ${DATADIR}/peaks "with your prefix," ${FILENAME}
else
	# Determine if chromosome has underscore for proper subsetting
	if [[ ${chrom} =~ "_" ]]
	then
   		zcat ${DATADIR}'/gls/'${BEAGLE}'.gz' | grep ${chrom} | awk '{print $1}' | awk -F_ -v s=${s} -v e=${e} '$3 >= s && $3 <= e' > ${DATADIR}'/peaks/'${FILENAME}'.txt'
	else
		zcat ${DATADIR}'/gls/'${BEAGLE}'.gz' | grep ${chrom} | awk '{print $1}' | awk -F_ -v s=${s} -v e=${e} '$2 >= s && $2 <= e' > ${DATADIR}'/peaks/'${FILENAME}'.txt'
	fi
	zcat ${DATADIR}'/gls/'${BEAGLE}'.gz' | grep -Fwf ${DATADIR}'/peaks/'${FILENAME}'.txt' | gzip > ${DATADIR}'/peaks/'${FILENAME}'.beagle.gz'
fi

############### Run a local PCA

echo "Starting covariance matrix calculation."

mkdir -p ${DATADIR}/peaks/pca/

module purge
module load bio/pcangsd/1.36.4

# run pcangsd to get covariance matrix
#pcangsd --threads 5 --iter 500 \
#	-b ${DATADIR}'/peaks/'${FILENAME}'.beagle.gz' \
#	-o ${DATADIR}'/peaks/pca/'${FILENAME}

# Filter two consistently missing data indivs
pcangsd --threads 5 --iter 500 \
        -b ${DATADIR}'/peaks/'${FILENAME}'.beagle.gz' \
        -o ${DATADIR}'/peaks/pca/'${FILENAME}'_filter' \
	--filter ${DATADIR}'/scripts/boreogadus_pca_rm_missingness_inds.txt'

############### done with pcangsd calc

