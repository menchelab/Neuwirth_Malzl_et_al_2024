#!/bin/bash
#
#SBATCH --job-name=cellranger
#SBATCH --cpus-per-task=8
#SBATCH --mem=50GB
#SBATCH --array=1-38
#SBATCH --output=out/cellranger_%a.out
#SBATCH --mail-type=ALL
#SBATCH --qos=short
#SBATCH --time=00-04:00:00
#SBATCH --mail-user=daniel.malzl@imp.ac.at

module load cellranger/6.0.1

pvid=$(cat /groups/pavri/bioinfo/daniel/scStary/cr_resources/sirolimus/sample_hto_demux.csv | cut -f 1 -d ',' | awk -v var=$SLURM_ARRAY_TASK_ID 'NR==var {print}')
echo $pvid
cellranger multi \
	--id=${pvid} \
	--csv=/groups/pavri/bioinfo/daniel/scStary/cr_resources/sirolimus/hto_demux_config_${pvid}.conf \
	--disable-ui \
	--localcores=8 \
	--localmem=50
