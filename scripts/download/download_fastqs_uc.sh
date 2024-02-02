#!/bin/bash
#
#SBATCH --job-name=download
#SBATCH --cpus-per-task=10
#SBATCH --mem=16GB
#SBATCH --output=download.out
#SBATCH --qos=medium
#SBATCH --mail-type=ALL
#SBATCH --mail-user=daniel.malzl@imp.ac.at

module load sra-toolkit/2.9.6-1-centos_linux64
for sra in `tail -n +2 accessions/UC_SRA_metadata.txt | cut -f 1 -d ","`;
do
	echo $sra
	fasterq-dump $sra -S -p -P -L 3 -e 8 -O fastqs
	gzip fastqs/${sra}* &
done
wait # wait for gzip processes to finish

python3 scripts/rename_fastqs.py -i accessions/UC_SRA_metadata.txt -d fastqs
