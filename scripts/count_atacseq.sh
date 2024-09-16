#!/bin/bash
#
#BATCH --job-name=count
#SBATCH --cpus-per-task=16
#SBATCH --mem=64GB
#SBATCH --output=count.out
#SBATCH --mail-type=ALL
#SBATCH --time=00-02:00:00
#SBATCH --mail-user=daniel.malzl@imp.ac.at

module load deeptools/3.3.1-foss-2018b-python-3.6.6

multiBamSummary BED-file \
	-b atacseq/bwa/merged_library/*bam \
	-o atacseq/atacseq_counts.npz \
	--BED atacseq/bwa/merged_replicate/macs2/broad_peak/consensus/consensus_peaks.mRp.clN.bed \
	--outRawCounts atacseq/atacseq_counts.tsv \
	-p 16 \
	--ignoreDuplicates \
	--smartLabels
