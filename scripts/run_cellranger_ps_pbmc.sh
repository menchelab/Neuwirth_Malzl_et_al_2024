#!/bin/bash
#
#SBATCH --job-name=cellranger
#SBATCH --cpus-per-task=8
#SBATCH --mem=40GB
#SBATCH --output=out/cellranger_ps_pbmc.out
#SBATCH --mail-type=ALL
#SBATCH --time=00-08:00:00
#SBATCH --mail-user=daniel.malzl@imp.ac.at

module load cellranger/6.0.1

cellranger count --id=BSF_1113_Treg_blood_1 \
		 --libraries=/groups/pavri/bioinfo/daniel/scStary/cr_resources/ps_pbmc/hto_demux_config.csv \
		 --feature-ref=/groups/pavri/bioinfo/daniel/scStary/cr_resources/ps_pbmc/hto_cmo_ref_0.csv \
		 --transcriptome=/resources/references/10x/refdata-cellranger-GRCh38-3.0.0 \
		 --nosecondary \
		 --disable-ui \
		 --localcores=8 \
		 --localmem=40
