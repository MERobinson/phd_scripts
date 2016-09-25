#!/bin/bash

RESDIR=/media/mark/HARWOOD-WD3/Harwood/Resources
GTF=$RESDIR/GFFs/Dicty_12.09.14/Combined_12.09.14_v2-12_masked.gtf
INDEXDIR=/Volumes/HARWOOD-WD3/Tools/bowtie2/bowtie2-2.2.6/indexes
GTFINDEX=$INDEXDIR/tophat_gtf_index"
GENINDEX=$INDEXDIR/D_discoideum_Ax4_May_2009-masked_renamed_chr_ddb"

# tophat2 to align
FASTQ=$(find -type f -name "*.fastq")
for FILE in $FASTQ; do
	NAME=$(basename "$FILE")
	NAME=${NAME%.*}
	tophat -p 8 -o ./$NAME --transcriptome-index $GTFINDEX $GENINDEX $FILE
done

# HTseq-count for raw counts
BAM=$(find -type f -name "*accepted_hits.bam")
for FILE in $BAM; do
	NAME=${FILE%/*}
        NAME=$(basename "$NAME")
	echo "Counting reads for $NAME"
	htseq-count --stranded=no --format=bam -q $FILE $GTF > $NAME.counts
done

# stringtie on just WT samples for norm exprs
BAM=$(find -type f -name "*accepted_hits.bam" | grep "WT_0H")
for FILE in $BAM; do
	NAME=${FILE%/*}
	NAME=$(basename "$NAME")
	echo "Analysing expression for $NAME"
	mkdir -p ballgown
	stringtie -e -B -p 8 -G $GTF -o ballgown/$NAME/${NAME}_strng.gtf $FILE
done

exit 0
