#!/bin/bash
#
#SBATCH --job-name=cellranger
#SBATCH --cpus-per-task=6
#SBATCH --mem=32GB
#SBATCH --array=1-36
#SBATCH --output=out/cellranger_%a.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=daniel.malzl@imp.ac.at

module load cellranger/6.0.1

gsm=$(tail -n +2 accessions/UC_SRA_metadata.txt | cut -f 21 -d "," | sort | uniq | awk -v var=$SLURM_ARRAY_TASK_ID 'NR==var {print}')
cd cellranger
cellranger count --id=${gsm} \
		 --fastqs=/scratch/daniel.malzl/scStary/fastqs \
		 --sample=${gsm} \
		 --transcriptome=/resources/references/10x/refdata-cellranger-GRCh38-3.0.0 \
		 --nosecondary \
		 --disable-ui \
		 --localcores=6 \
		 --localmem=32
