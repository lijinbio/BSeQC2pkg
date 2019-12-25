#!/usr/bin/env bash
# vim: set noexpandtab tabstop=2:

xlab='Position in read (bp)'
ylab='Methylation level'
R --slave --no-save --no-restore --no-init-file \
	-e "width=6" \
	-e "height=3" \
	-e "pico=F" \
	-e "infile='../result_mbias_pe.txt'" \
	-e "outfile='main_pair.pdf'" \
	-e "xlab='$xlab'" \
	-e "ylab='$ylab'" \
	-e "source('bseqc2mbiasplot.R', echo=F)"

