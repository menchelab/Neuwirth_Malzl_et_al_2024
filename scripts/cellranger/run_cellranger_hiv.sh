#!/bin/bash
#
#SBATCH --job-name=cellranger
#SBATCH --cpus-per-task=8
#SBATCH --mem=40GB
#SBATCH --array=1-7
#SBATCH --output=out/cellranger_%a.out
#SBATCH --mail-type=ALL
#SBATCH --time=00-08:00:00
#SBATCH --mail-user=daniel.malzl@imp.ac.at

module load cellranger/6.0.1

sample_id=$(cat /groups/pavri/bioinfo/daniel/scStary/cr_resources/hiv/sample_hto_demux.csv | cut -f 1 -d ',' | awk -v var=$SLURM_ARRAY_TASK_ID 'NR==var {print}')
feature_ref=$(cat /groups/pavri/bioinfo/daniel/scStary/cr_resources/hiv/sample_hto_demux.csv | cut -f 2 -d ',' | awk -v var=$SLURM_ARRAY_TASK_ID 'NR==var {print}')
libraries=$(cat /groups/pavri/bioinfo/daniel/scStary/cr_resources/hiv/sample_hto_demux.csv | cut -f 3 -d ',' | awk -v var=$SLURM_ARRAY_TASK_ID 'NR==var {print}')

echo $sample_id
echo $feature_ref
echo $libraries

cellranger count --id=${sample_id} \
		 --libraries=/groups/pavri/bioinfo/daniel/scStary/cr_resources/hiv/$libraries \
		 --feature-ref=/groups/pavri/bioinfo/daniel/scStary/cr_resources/hiv/$feature_ref \
		 --transcriptome=/resources/references/10x/refdata-cellranger-GRCh38-3.0.0 \
		 --nosecondary \
		 --disable-ui \
		 --localcores=8 \
		 --localmem=40
