# README #

This repository hosts scripts and tools for RNAseq analysis

#### tuxedo-slurm.sh
this script runs a full RNAseq pipeline under a slurm jobs distribution system <br />
fastQC - quality control <br />
flexbar - adapters and quality trimming <br />
tophat - aligner <br />
cufflinks - transcripts assembly and quantification <br />
cuffmerge - merging of assemblies <br />
cuffcompare - corrrection of lost protein ids <br />
cuffquant - transcript expression profiles <br />
cuffdiff - differential expression analysis <br />

Please read the instructions inside the script for usage

#### tuxedo_v3-slurm.sh 
this script runs a full RNAseq pipeline under a slurm jobs distribution system <br />
using about 18 processes per file allowing full analysis of 20 libraries with <br />
50 - 150 M reads per library to complete under 12h <br />
fastQC - quality control <br />
flexbar - adapters and quality trimming <br />
hisat - aligner <br />
stringtie - transcripts assembly and quantification <br />
cuffmerge - merging of assemblies <br />
cuffquant - transcript expression profiles <br />
cuffdiff - differential expression analysis <br />

Please read the instructions inside the script for usage <br />

####  aDiff
this python script takes cuffdiff results as inputs and generates excel report <br />
tables for each parallel comparison. It annotates each gene with the respective <br />
GO term. It also uses DAVID to perform GO enrichment analysis of biological process, <br />
cellular components,  and molecular function using significantly changed genes, <br />
transcripts, promoter usage, splicing, CDS.

Please read the instructions inside the script for usage <br />