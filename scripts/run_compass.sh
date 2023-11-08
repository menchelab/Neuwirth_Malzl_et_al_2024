#!/bin/bash
#
#SBATCH --job-name=compass
#SBATCH --cpus-per-task=16
#SBATCH --mem=64GB
#SBATCH --array=1-45
#SBATCH --output=out/compass_%a.out
#SBATCH --mail-type=ALL
#SBATCH --qos=long
#SBATCH --time=14-00:00:00
#SBATCH --mail-user=daniel.malzl@imp.ac.at

sample=$(cat data/sample.txt | awk -v var=$SLURM_ARRAY_TASK_ID 'NR==var {print}')

mkdir compass/${sample}_tregs
compass --data data/skin_inflammatory_disease.${sample}.expression.tsv \
	--output-dir compass/${sample}_tregs \
	--num-processes 15 \
	--species homo_sapiens \
	--latent-space data/skin_inflammatory_disease.${sample}.latent.tsv
