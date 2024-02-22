#!/bin/bash
#
#BATCH --job-name=count
#SBATCH --cpus-per-task=8
#SBATCH --mem=32GB
#SBATCH --output=count.out
#SBATCH --mail-type=ALL
#SBATCH --time=00-01:00:00
#SBATCH --mail-user=daniel.malzl@imp.ac.at

module load deeptools/3.3.1-foss-2018b-python-3.6.6

bamdir=results/02_alignment/bowtie2/target/markdup/

for name in h3k4me1 h3k27ac;
do
	for peaks in all reproducible;
	do
		multiBamSummary BED-file \
			-b $bamdir/${name}*bam \
			-o counts/${name}.${peaks}.npz \
			--BED merged_peaks/${name}.${peaks}.bed \
			--outRawCounts counts/${name}.${peaks}.tsv \
			--maxFragmentLength 1000 \
			--smartLabels \
			-p 8
	done
done

multiBamSummary BED-file \
	-b $bamdir/h*bam \
	-o counts/genes.extended.npz \
	--BED genes.extended.bed \
	--outRawCounts counts/genes.extended.tsv \
	--maxFragmentLength 1000 \
	--smartLabels \
	-p 8

multiBamSummary BED-file \
	-b $bamdir/h*bam \
	-o counts/genes_tss_promoters.npz \
	--BED genes_tss_promoters.bed \
	--outRawCounts counts/genes_tss_promoters.tsv \
	--maxFragmentLength 1000 \
	--smartLabels \
	-p 8
