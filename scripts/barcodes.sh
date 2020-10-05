for f in attic/*_1.fq.gz; do
	printf '%s\t%s\n' $f `zcat $f | head -n 1 | sed 's/.*://'`
done
