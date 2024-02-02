#!/bin/bash
#
#SBATCH --job-name=download
#SBATCH --cpus-per-task=1
#SBATCH --mem=8GB
#SBATCH --output=download.out
#SBATCH --time=120
#SBATCH --mail-type=ALL
#SBATCH --mail-user=daniel.malzl@imp.ac.at

url=https://biomedical-sequencing.at/projects/BSA_0347_TK_SRX_COHORT_f593ab2948624fe793db39a848dbbb2c/RAW
for sample in `cat accessions/srcx_samples.txt`;
do
	wget -r -np -nH --cut-dirs=3 -R html -e robots=off $url/${sample}_Transcriptome/ -P fastqs
	wget -r -np -nH --cut-dirs=3 -R html -e robots=off $url/${sample}_Antibody_Multiplex_CITE/ -P fastqs
done
