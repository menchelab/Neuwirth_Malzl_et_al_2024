#!/bin/bash
#
#SBATCH --job-name=aggr
#SBATCH --cpus-per-task=8
#SBATCH --mem=64GB
#SBATCH --output=aggr.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=daniel.malzl@imp.ac.at

module load cellranger/6.0.1
cd cellranger
cellranger aggr --id=UC \
	 	--csv=../accessions/aggr_uc.csv \
		--nosecondary \
		--disable-ui \
		--normalize=none

cp UC/outs/count/filtered_feature_bc_matrix.h5 ../raw/uc_cellranger.filtered.h5
python3 ../scripts/make_aggr_metadata.py -i ../raw/uc_cellranger.filtered.h5 \
				         -c ../accessions/aggr_uc.csv \
				         -o ../raw/uc_cellranger.metadata.tsv
