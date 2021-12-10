#subset LuCaP DMRs based on whether they overlap promoters or gene bodies
mkdir -p out

#ensure that adjacent dmrs are merged:
#bedtools merge -i ../analysis/dmrs/lucap_cf_nepc_vs_prad_train_test_final/rms.padj0.001.lfc2/reference_up.bed -d 1 > tmp.up.bed
bedtools merge -i ../analysis/dmrs/lucap_cf_nepc_vs_prad_train_test_final/rms.padj0.001.lfc2/reference_up_filtered.bed -d 1 > tmp.up.bed
#bedtools merge -i ../analysis/dmrs/lucap_cf_nepc_vs_prad_train_test_final/rms.padj0.001.lfc2/reference_down.bed -d 1 > tmp.down.bed
bedtools merge -i ../analysis/dmrs/lucap_cf_nepc_vs_prad_train_test_final/rms.padj0.001.lfc2/reference_down_filtered.bed -d 1 > tmp.down.bed

#NEPC-up DMRs overlapping promoters
bedtools intersect -a tmp.up.bed -b regions/promoters.3kb.bed > out/up_promoter.bed

#NEPC-up DMRs overlapping gene bodies
bedtools intersect -a tmp.up.bed -b regions/gene.body.no.promoter.bed > out/up_genebody.bed

#PRAD-up DMRs overlapping promoters
bedtools intersect -a tmp.down.bed -b regions/promoters.3kb.bed > out/down_promoter.bed

#PRAD-up DMRs overlapping gene bodies
bedtools intersect -a tmp.down.bed -b regions/gene.body.no.promoter.bed > out/down_genebody.bed

#rm tmp.up.bed tmp.down.bed
