#!/bin/bash
#SBATCH --job-name=satsuma
#SBATCH --mail-user=natasha.howe@noaa.gov
#SBATCH --mail-type=FAIL,END
#SBATCH --mem=120GB
#SBATCH --partition=himem
#SBATCH --time=6-12:00:00
#SBATCH --cpus-per-task=15
#SBATCH --output=/home/nhowe/arctic/job_outfiles/satsuma_synteny.%j.out

module load bio/satsuma/r41
module load bio/satsuma/r41

/opt/bioinformatics/bio/satsuma/satsuma-r41/SatsumaSynteny -n 15 -q /home/nhowe/reference_genomes/cod/V2/gadMor2.fasta -t /home/nhowe/reference_genomes/cod/V3/GCF_902167405.1_gadMor3.0_genomic.fna -o /home/nhowe/arctic/genome_synteny/