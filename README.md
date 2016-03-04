# HTSeq-tools

This repository hosts scripts and tools for analysis of high throughput sequencing data.

## Tools

#### ensref

*Made for SLURM & Environment Modules Project*

This is the script we use for downloading genome assemblies from ENSEMBL and 
build the respective indexes and folders structure. 
 
For usage check the help output: 
```
./ensref --help
```

The former script ENSrefs is deprecated. See the history with

```
git log -- ENSrefs build-ref.sh
git log -p -- ENSrefs build-ref.sh # with diff
```

#### CCG

this script takes a wget command from CCG donwload instructions, 
downloads the respective files to a destination of choice checking file integrity afterwards.

usage: 
```
CCG -l "wget ..." -o destination_folder
```

#### getChromosome

This script allows you to generate a fasta file for a chromosome contained in a multifasta file.

Usage:
```
getChromosome -h
# For >2 dna:chromosome chromosome:GRCm38:2:1:182113224:1 REF use:
getChromosome -i GRCm38.dna.primary_assembly.fa -o output_folder -c 2
```

#### make.genome

Create genome file for bedtools. The output is a tab-delimited file of chromosome name followed by its length.

Usage:
```
make.genome -h
```

#### GATKlift

This scripts lifts GTF and fasta files from GATK reconstructions (ie. FastaAlternateReferenceMaker).

Usage:
```
GATKlift -h
```

#### tuxedo-slurm.sh

*Made for SLURM & Environment Modules Project*

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

#### tuxedo_v3-slurm.sh

*Made for SLURM & Environment Modules Project*

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

#### mart

A command line tool to annotate ids in a column of tables (.tsv, .csv, .xls[x] format)
with the R/biomaRt tool and optionally merge the result with the original table.

Dependencies:
* python >= 2.7
* python/argparse
* R (tested on version 3.2)
* R/WriteXLS
* R/argparse
* R/biomaRt
* R/plyr
* R/readxl

Usage:
```
mart -h
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

#### SNPsFilter

This executable interactive python scripts lets you merge `*.tabular files 
outputed by [CloudMap](http://usegalaxy.org/cloudmap) and filter SNPs of 
interest based on the maximum number of samples containing the SNP and minimum
number of SNPs per gene.

Usage:

```
SNPsFilter
```

#### meran-slurm.sh

A HPC/slurm pipeline to perform RNA methylation analysis.
Modify the sample sheet and path variables at the beginning of the script.
Run script on trimmed read data. Submit in slurm environment.

#### GATK_lofreq_slurm.sh

*Made for SLURM & Environment Modules Project*

This is a standard GATK and lofreq variant call pipeline with snpEff annotations.

It is made to make calls on two lines eg. wt and mutant.

It is set for single end sequencing.

Please read the instructions inside the script for usage.


## License

> Copyright 2015
> * Jorge Boucas <JBoucas@age.mpg.de>
> * Sven Templer <templer@age.mpg.de>

See file `LICENSE`.
