# rnaseq-tools

This repository hosts scripts and tools for RNAseq analysis.

## Tools

#### ENSrefs

*Made for SLURM*

This is the script we use for downloading genome assemblies from ENSEMBL and 
build the respective indexes and folders structure. 
 
Usage: 
```
srun ENSrefs destination_folder
```



### tuxedo-slurm.sh

*Made for SLURM*

This script runs a full RNAseq pipeline under a slurm jobs distribution system

* fastQC - quality control 
* flexbar - adapters and quality trimming 
* tophat - aligner 
* cufflinks - transcripts assembly and quantification 
* cuffmerge - merging of assemblies 
* cuffcompare - corrrection of lost protein ids 
* cuffquant - transcript expression profiles 
* cuffdiff - differential expression analysis 

Please read the instructions inside the script for usage.

### tuxedo_v3-slurm.sh

*Made for SLURM*

This script runs a full RNAseq pipeline under a slurm jobs distribution system 
using about 18 processes per file allowing full analysis of 20 libraries with 
50 - 150 M reads per library to complete under 12h 

* fastQC - quality control 
* flexbar - adapters and quality trimming 
* hisat - aligner 
* stringtie - transcripts assembly and quantification 
* cuffmerge - merging of assemblies 
* cuffquant - transcript expression profiles 
* cuffdiff - differential expression analysis 

Please read the instructions inside the script for usage. 

####  aDiff

This python script takes cuffdiff results as inputs and generates excel report 
tables for each parallel comparison. It annotates each gene with the respective 
GO and KEGG term. It also uses DAVID to perform GO enrichment analysis of biological process, 
cellular components, molecular function and other DAVID annotations using significantly changed genes, 
transcripts, promoter usage, splicing, CDS.

Usage: 

```
aDiff -h
```

#### QC.R

This R script performs a minimal quality control analysis on cuffdiff outputs. 
The output folder needs to be created before using this script. 
Not all plot allways perform with every dataset - you might therefore need to 
comment out some plots.

This R scripts requires the `cummeRbund` package. To install it, run the following
in an `R` console:

```
source('http://www.bioconductor.org/biocLite.R')
biocLite('cummeRbund', ask=FALSE)
```

Usage:

```
mkdir cummeRbund_output_folder
Rscript QC.R cuffdiff_output_folder cummeRbund_output_folder
```

## License

> Copyright 2015
> * Jorge Boucas <JBoucas@age.mpg.de>
> * Sven Templer <templer@age.mpg.de>

See file `LICENSE`.
