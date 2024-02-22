module load bedtools/2.27.1-foss-2018b

peakdir=results/03_peak_calling/04_called_peaks/seacr

mkdir tmp

for name in h3k27ac h3k4me1;
do
	for gtype in nt sat1;
	do
		for i in 1 2;
		do
			cut -f 1,2,3 ${peakdir}/${name}_${gtype}_R${i}.seacr.peaks.stringent.bed > tmp/${name}_${gtype}_R${i}.bed
		done

		# reproducible peaks
		bedtools intersect -wa -a tmp/${name}_${gtype}_R1.bed -b tmp/${name}_${gtype}_R2.bed > tmp/${name}_${gtype}_12.bed
		bedtools intersect -wa -b tmp/${name}_${gtype}_R1.bed -a tmp/${name}_${gtype}_R2.bed > tmp/${name}_${gtype}_21.bed
		cat tmp/${name}_${gtype}_12.bed tmp/${name}_${gtype}_21.bed | sort -k1,1 -k2,2n > tmp/${name}_${gtype}.reproducible.sort.bed
		bedtools merge -i tmp/${name}_${gtype}.reproducible.sort.bed > merged_peaks/${name}_${gtype}.reproducible.bed

		# all peaks
		cat tmp/${name}_${gtype}_R1.bed tmp/${name}_${gtype}_R2.bed | sort -k1,1 -k2,2n > tmp/${name}_${gtype}.all.sort.bed
		bedtools merge -i tmp/${name}_${gtype}.all.sort.bed > merged_peaks/${name}_${gtype}.all.bed
	done

	for peaks in all reproducible;
	do
		cat merged_peaks/${name}_nt.${peaks}.bed merged_peaks/${name}_sat1.${peaks}.bed | sort -k1,1 -k2,2n > tmp/${name}.${peaks}.sort.bed
		bedtools merge -i tmp/${name}.${peaks}.sort.bed > merged_peaks/${name}.${peaks}.bed
	done
done

rm -r tmp
