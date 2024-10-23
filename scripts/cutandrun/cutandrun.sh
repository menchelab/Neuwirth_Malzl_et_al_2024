module load build-env/f2022
module load nextflow/23.04.2

nextflow run nf-core/cutandrun \
	-c resource/cnr_config.conf \
	-profile cbe \
	-w /scratch-cbe/users/daniel.malzl/cpm \
	--input resource/cnr_samples.csv \
	--outdir results \
	--genome GRCh38 \
	--peakcaller seacr,macs2 \
	--blacklist resource/cnr_blacklist_hg38.bed \
	--normalisation_mode CPM \
	--seacr_stringent stringent \
	--seacr_norm norm \
	--igg_scale_factor 1 \
	--replicate_threshold 1 \
	-resume
