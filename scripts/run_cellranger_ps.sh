#!/bin/bash
#
#SBATCH --job-name=cellranger
#SBATCH --cpus-per-task=8
#SBATCH --mem=40GB
#SBATCH --array=1-96
#SBATCH --output=out/cellranger_%a.out
#SBATCH --mail-type=ALL
#SBATCH --qos=medium
#SBATCH --time=01-00:00:00
#SBATCH --mail-user=daniel.malzl@imp.ac.at

module load cellranger/6.0.1

gsm=$(tail -n +2 accessions/PS_ENA_metadata.tsv | cut -f 1 -d$'\t' | sort | awk -v var=$SLURM_ARRAY_TASK_ID 'NR==var {print}')
echo $gsm
cellranger count --id=${gsm} \
		 --fastqs=/groups/pavri/bioinfo/daniel/scStary/fastqs \
		 --sample=${gsm} \
		 --transcriptome=/resources/references/10x/refdata-cellranger-GRCh38-3.0.0 \
		 --nosecondary \
		 --disable-ui \
		 --localcores=8 \
		 --localmem=40
