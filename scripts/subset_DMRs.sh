
ref=$1 #differential peak calls in refernce tissue (eg, LuCaPs)
exclude=$2
outdir=$3

padj=0.000001
#padj=0.01
lfc=3
#lfc=0

echo number of differential peaks with q less than 0.05
wc -l $ref

#DMRs that are up in case compared to control reference tissues (ie, up in NEPC)
awk -v p=$padj -v c=$lfc 'BEGIN{OFS="\t"} $6 > c && $7<p {print $1,$2,$3}' $ref > ${outdir}/reference_up.bed
echo peaks that are up in reference by the specified thresholds:
wc -l ${outdir}/reference_up.bed

bedtools intersect -v -a ${outdir}/reference_up.bed -b $exclude > ${outdir}/reference_up_filtered.bed
echo number of these peaks left after removing overlap with WBC MeDIP data:
wc -l ${outdir}/reference_up_filtered.bed

#DMRs that are up in control compared to case reference tissues (ie, up in NEPC)
awk -v p=$padj -v c=$lfc 'BEGIN{OFS="\t"} $6 < -1*c && $7<p {print $1,$2,$3}' $ref > ${outdir}/reference_down.bed
echo peaks that are up in reference by the specified thresholds:
wc -l ${outdir}/reference_down.bed

bedtools intersect -v -a ${outdir}/reference_down.bed -b $exclude > ${outdir}/reference_down_filtered.bed
echo number of these peaks left after removing overlap with WBC MeDIP data:
wc -l ${outdir}/reference_down_filtered.bed


