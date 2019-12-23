bin_PROGRAMS = bseqc2

samtools_INCLUDE = $(top_srcdir)/lib/samtools-0.1.20
samtools_LIB = $(top_srcdir)/lib/samtools-0.1.20

CXXFLAGS = -g -O3 -std=c++11 -static
bseqc2_CPPFLAGS = -Wall -w -I$(samtools_INCLUDE)
bseqc2_LDFLAGS = -L$(samtools_LIB)
bseqc2_LDADD = -lbam -lz -lpthread -lboost_program_options
bseqc2_SOURCES = bseqc2.cpp
