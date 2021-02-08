#summarize info about coefficients used in models for each iteartion (eg, LFC in cfDNA, LFC in LuCaPs, was it a signficant DMR, which iteration did it come frome
DIR=$1
OUTFILE=$2
reference=$3

printf '' > ${OUTFILE}.tmp
for coef in ${DIR}/iter.*/*coef.bed; do
	sig=`echo $coef | sed 's/coef/diff/'`
	#annotate the coefficient bed file wtih whether the DMR was signficant in cfDNA (1 = yes, 0 = no) and the iteration number
	iter=`echo $coef | sed 's/.*iter.//' | sed 's/\/.*//'`
	echo summarizing model coefficients for iteration $iter
	bedtools intersect -wao -a $coef -b $sig | \
		awk 'BEGIN{FS="\t";OFS="\t"} {if($NF==0){print $1, $2, $3, $4, $5, $6, $7, $8, 0} else {print $1, $2, $3, $4, $5, $6, $7, $8, 1}}' | \
		awk -v iter=$iter 'BEGIN{FS="\t";OFS="\t"} $8!=0 {print $0, iter}' >> ${OUTFILE}.tmp 
done
printf 'chr\tstart\tend\tmean1.plasma\tmean2.plasma\tlfc.plasma\tpadj.plasma\tcoef\tsig.plasma\titeration\tmean1.lucap\tmean2.lucap\tlfc.lucap\tpadj.lucap\n' > $OUTFILE
bedtools intersect -wo -a ${OUTFILE}.tmp -b $reference | cut --complement -f11-13,18 >> ${OUTFILE}


