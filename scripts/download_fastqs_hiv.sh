#!/bin/bash
#
#SBATCH --job-name=download
#SBATCH --cpus-per-task=1
#SBATCH --mem=8GB
#SBATCH --output=download.out
#SBATCH --time=120
#SBATCH --mail-type=ALL
#SBATCH --mail-user=daniel.malzl@imp.ac.at

url=https://biomedical-sequencing.at/projects/BSA_0512_GSTARY_HIV_8c51b8f89482493c813a7d6a79615cf1/RAW
for sample in `cat accessions/hiv_samples.txt`;
do
	wget -r -np -nH --cut-dirs=3 -R html -e robots=off $url/${sample}_Transcriptome/ -P fastqs
	wget -r -np -nH --cut-dirs=3 -R html -e robots=off $url/${sample}_Antibody_Multiplex/ -P fastqs
done
