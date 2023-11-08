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
cellranger aggr --id=PS \
	 	--csv=accessions/aggr_ps.csv \
		--nosecondary \
		--disable-ui \
		--normalize=none

cp PS/outs/count/filtered_feature_bc_matrix.h5 raw/ps_cellranger.filtered.h5
python3 scripts/make_aggr_metadata.py -i raw/ps_cellranger.filtered.h5 \
				      	 -c accessions/aggr_ps.csv \
				         -o raw/ps_cellranger.metadata.tsv
