
ref=$1 #differential peak calls in refernce tissue (eg, LuCaPs)
exclude=$2
outdir=$3

echo number of differential peaks with q less than 0.05
wc -l $ref

awk 'BEGIN{OFS="\t"} $6 > 3 && $7<0.000001 {print $1,$2,$3}' $ref > ${outdir}/reference_up.bed
echo peaks that are up in reference by the specified thresholds:
wc -l ${outdir}/reference_up.bed

bedtools intersect -v -a ${outdir}/reference_up.bed -b $exclude > ${outdir}/reference_up_filtered.bed
echo number of these peaks left after removing overlap with WBC MeDIP data:
wc -l ${outdir}/reference_up_filtered.bed

