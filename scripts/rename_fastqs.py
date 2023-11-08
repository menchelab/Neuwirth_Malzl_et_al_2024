import argparse as ap
import pandas as pd
import os

parser = ap.ArgumentParser()
parser.add_argument(
	'-i', 
	'--csv', 
	required = True, 
	help = 'csv file containing the metadata of each SRA file as downloaded from SRA'
)
parser.add_argument(
	'-d', 
	'--directory', 
	required = True, 
	help = 'directory where the downloaded fastqs reside'
)
args = parser.parse_args()
meta = pd.read_csv(args.csv)
for samplename, runs in meta.groupby('Sample Name'):
	for samplenum, runid in enumerate(runs.Run):
		samplenum += 1
		for i in [1, 2]:
			old = os.path.join(
				args.directory,
				f'{runid}_{i}.fastq.gz'
			)
			new = os.path.join(
				args.directory,
				f'{samplename}_S{samplenum}_L001_R{i}_001.fastq.gz'
			)
			if not os.path.exists(new):
				os.rename(old, new)
