module load build-env/f2022
module load nextflow/23.10.1

nextflow run nf-core/atacseq \
	--input resource/atac_samples.csv \
	--outdir atacseq \
	--genome GRCh38 \
	--aligner bwa \
	--read_length 150 \
	--narrow_peak \
	-w /scratch-cbe/users/daniel.malzl/work/ \
	-profile cbe \
	-resume
