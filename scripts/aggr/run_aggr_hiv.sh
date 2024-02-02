#!/bin/bash
#
#SBATCH --job-name=aggr
#SBATCH --cpus-per-task=8
#SBATCH --mem=64GB
#SBATCH --output=aggr.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=daniel.malzl@imp.ac.at

module load cellranger/6.0.1
cellranger aggr --id=HIV \
	 	--csv=/groups/pavri/bioinfo/daniel/scStary/cr_resources/hiv/aggr.csv \
		--nosecondary \
		--disable-ui \
		--normalize=none

cp HIV/outs/count/filtered_feature_bc_matrix.h5 /groups/pavri/bioinfo/daniel/scStary/raw/hiv_cellranger.filtered.h5
