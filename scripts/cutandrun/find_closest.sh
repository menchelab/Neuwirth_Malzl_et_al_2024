module load bedtools/2.27.1-foss-2018b

mkdir tmp

cat genes.bed | sort -k1,1 -k2,2n > tmp/genes.bed

for peaks in `ls merged_peaks`;
do
	cat merged_peaks/$peaks | sort -k1,1 -k2,2n > tmp/$peaks
	bedtools closest -a tmp/$peaks -b tmp/genes.bed -t all -D a > closest/${peaks%.*}.closestgenes.bed
done
