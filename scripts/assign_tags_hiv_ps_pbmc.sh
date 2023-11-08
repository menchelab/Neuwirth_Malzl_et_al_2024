#!/bin/bash
#
#SBATCH --job-name=cellranger
#SBATCH --cpus-per-task=4
#SBATCH --mem=32GB
#SBATCH --output=assign_tags.out
#SBATCH --mail-type=ALL
#SBATCH --time=00-08:00:00
#SBATCH --mail-user=daniel.malzl@imp.ac.at

counhto --csv accessions/hto_samples.csv -p 4
