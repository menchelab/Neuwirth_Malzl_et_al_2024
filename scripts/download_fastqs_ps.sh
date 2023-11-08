#!/bin/bash
#
#SBATCH --job-name=download
#SBATCH --cpus-per-task=2
#SBATCH --mem=16GB
#SBATCH --output=download.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=daniel.malzl@imp.ac.at

IFS=$'\n' # setting newline to sole separator to avoid generic splitting on whitespaces in for loops
echo -e 'sample_id\tpatient_id$Status\tTissue\tCell_fraction\tread1_url\tread2_url' > accessions/PS_ENA_metadata.tsv
for sample in `tail -n +2 accessions/PS_ENA_metadata.txt | cut -f 1,5,8,10,11,57,59 -d$'\t'`;
do
	IFS=$'\t' read sid state location fraction read1 read2 <<< $sample
	echo $sid	
	echo $sample >> accessions/PS_ENA_metadata.tsv
	for url in $read1 $read2;
	do
		file=$( echo $read1 | cut -f 11 -d / )
		if ! [ -f "fastqs/$file" ]; then
			wget $url -P fastqs &
		fi
	done
	wait
done
