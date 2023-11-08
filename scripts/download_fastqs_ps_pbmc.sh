#!/bin/bash
#
#SBATCH --job-name=download
#SBATCH --cpus-per-task=1
#SBATCH --mem=8GB
#SBATCH --output=download.out
#SBATCH --time=120
#SBATCH --mail-type=ALL
#SBATCH --mail-user=daniel.malzl@imp.ac.at

url=https://biomedical-sequencing.at/projects/BSA_0581_CoR_HSCT1121_0ed14a948c0c4b9a99605e39891c6756/RAW
wget -r -np -nH --cut-dirs=3 -R html -e robots=off $url/Treg_blood_Transcriptome/ -P fastqs
wget -r -np -nH --cut-dirs=3 -R html -e robots=off $url/Treg_blood_Antibody_Multiplex/ -P fastqs
