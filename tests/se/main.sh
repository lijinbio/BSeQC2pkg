#!/usr/bin/env bash
# vim: set noexpandtab tabstop=2:

set -v
hg19=/mnt/folders/resource/genome/hg19/fasta/hg19_25chr/hg19.fa
f=T1L1.bam
../../src/bseqc2/bseqc2 -h
../../src/bseqc2/bseqc2 -i "$f" -l 50 -r "$hg19" -o outdir/result.txt
