wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_45/gencode.v45.annotation.gtf.gz

python3 gtf_to_bed.py \
	-it gene \
	-a gene_id gene_name \
	-i gencode.v45.annotation.gtf.gz \
	-o genes.extended.bed \
	-u 1000

python3 gtf_to_bed.py \
	-it gene \
	-a gene_id gene_name \
	-i gencode.v45.annotation.gtf.gz \
	-o genes.bed \
	-tss \
	-promoter

cat genes.bed genes.promoter.bed genes.tss.bed | sort -k1,1 -k2,2n > genes_tss_promoters.bed
