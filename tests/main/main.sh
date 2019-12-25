#!/usr/bin/env bash
# vim: set noexpandtab tabstop=2:

xlab='Position in read (bp)'
ylab='Methylation level'
R --slave --no-save --no-restore --no-init-file \
	-e "width=6" \
	-e "height=3" \
	-e "infile='../result_mbias_strand.txt'" \
	-e "outfile='main.pdf'" \
	-e "xlab='$xlab'" \
	-e "ylab='$ylab'" \
	-e "source('bseqc2mbiasplot.R', echo=T)"

