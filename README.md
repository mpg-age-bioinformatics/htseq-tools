# README #

This repository hosts scripts and tools for RNAseq analysis

# tuxedo-slurm.sh #
this script runs an full RNAseq pipeline under a slurm jobs distribution system
\nfastQC - quality control
\nflexbar - adapters and quality trimming
\ntophat - aligner
\ncufflinks - transcripts assembly and quantification
\ncuffmerge - merging of assemblies
\ncuffcompare - corrrection of lost protein ids
\ncuffquant - transcript expression profiles
\ncuffdiff - differential expression analysis

# tuxedo-slurm.sh #
this script runs an full RNAseq pipeline under a slurm jobs distribution system \
using about 18 processes per file allowing full analysis of 20 libraries with \
50 - 150 M reads per libary to complete under 12h
fastQC - quality control
flexbar - adapters and quality trimming
hisat - aligner
stringtie - transcripts assembly and quantification
cuffmerge - merging of assemblies
cuffquant - transcript expression profiles
cuffdiff - differential expression analysis

# aDiff #
this python script takes cuffdiff results as inputs and generates excel report \
tables for each parallel comparison. It annotates each gene with the respective \
GO term. It also uses DAVID to perform GO enrichment analysis of biological process, \
cellular components,  and molecular function using significantly changed genes, \
transcripts, promoter usage, splicing, CDS  - biological process, cellular \
component, molecular function.