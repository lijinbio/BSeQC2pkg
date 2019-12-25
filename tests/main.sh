#!/usr/bin/env bash
# vim: set noexpandtab tabstop=2:

set -v
f=LC1_chr_1k.bam
pefiltertag -i "$f" -s -t 4
../src/bseqc2/bseqc2 -h
tmpdir=$(mktemp -u -d)
../src/bseqc2/bseqc2 -i "$f" -l 160 -r /mnt/folders/resource/genome/hg38/fasta/hg38.fa -o "$tmpdir/result.txt"
tree "$tmpdir"
