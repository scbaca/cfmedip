#compare differentially methylated nucleotides from Misha Beltran's 2018 Nature Medicine paper to LuCaP DMRs
RBBS='NEPC_PRAD_RBBS.txt'
LuCaP='../analysis/dmrs/lucap_cf_nepc_vs_prad_train_test_final/reference/reference_DMRs.bed' #this has differential peaks between all LuCaP PRAD and NEPC included in this analysis.
OUT='out/overlap.bed'

mkdir -p out

printf 'chr\tstart\tend\tLuCaP\tRBBS\n' > $OUT
sed '1d' $RBBS | bedtools intersect -wo -a $LuCaP -b - | cut -f1-3,6,13 >> $OUT

Rscript plot.R
