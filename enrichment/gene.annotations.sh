#get promoter / gene body annotations.

mkdir -p regions

if [ ! -f "hg19.genome" ]; then
	 mysql --user=genome --host=genome-mysql.cse.ucsc.edu -A -e "select chrom, size from hg19.chromInfo" | grep -v _ | grep -v size | sort -k1,1 -k2,2n > hg19.genome
fi

sed '1d' NCBI_RefSeq_hg19.txt | awk 'BEGIN{FS="\t";OFS="\t"}{if($3=="+") {print $2,$4-1, $4} else {print $2,$5-1, $5}}' | grep -v "_" | sort -k1,1 -k2,2n | bedtools slop -i - -b 3000 -g hg19.genome  | bedtools merge -i - > regions/promoters.3kb.bed

sed '1d' NCBI_RefSeq_hg19.txt | awk 'BEGIN{FS="\t";OFS="\t"}{print $2,$4, $5}' | grep -v "_" | sort -k1,1 -k2,2n | bedtools merge -i - | bedtools subtract -a - -b regions/promoters.3kb.bed > regions/gene.body.no.promoter.bed

