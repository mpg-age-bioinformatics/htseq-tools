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

This scripts lifts fasta files from GATK reconstructions (ie. FastaAlternateReferenceMaker).

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

#### chains_slurm.sh

*Made for SLURM & Environment Modules Project*

This template lets you create chain files for 2 fasta files of the same species
after genome reconstruction or de-novo assembly.

This script should be run from project/scripts/chains_slurm.sh

More info can be found here: http://hgwdev.cse.ucsc.edu/~kent/src/unzipped/hg/doc/liftOver.txt and here: http://genomewiki.ucsc.edu/index.php/Minimal_Steps_For_LiftOver

Usage:
```
./chains_slurm.sh /path/to/old/fasta.fa /path/to/new/fasta.fa
```

The respective chains file can be found in chains_output/oldFasta_To_newFasta.chain.gz

#### CellPlot

Working Python implementation of CellPlot and SymPlot from the CellPlot package for R.

```
usage: CellPlot [-h] [-i INPUT] [-n NTERMS] [-o OUTPUT]

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        Tab separated tables with the following collumns from
                        an enrichment analysis: 'Enrichment', 'Significant',
                        'Annotated', 'Term', and 'log2fc'. For log2fc each
                        cell must contain a comma separated string with the
                        log2fc for the genes enriched in the respective term.
                        eg. '-Inf,-1,2,3.4,3.66,Inf' (default: None)
  -n NTERMS, --nterms NTERMS
                        Number of terms to display (default: 10)
  -o OUTPUT, --output OUTPUT
                        /path/to/output/prefix This will create the files
                        prefix.CellPlot.svg, prefix.CellPlot.png,
                        prefix.SymPlot.svg, prefix.SymPlot.png (default: None)
``` 

#### topgo

Run gene ontology enrichment analysis using R/Bioconductor's topGO package.

Usage:

```
usage: topgo [-h] [-i NAME] [-e NAME] [-n NAME [NAME ...]] [-f FORMAT]
             [-o OUTPUT_PREFIX] [-x] [-v] [-w]
             organism table [table ...]

topgo - perform gene ontology enrichment analysis using Bioconductor/topGO

positional arguments:
  organism              Select organism [second field as identifier from
                        Bioconductor AnnotationDbi 'org' databases]. E.g. 'Dm'
                        for 'org.Dm.eg.db'. Supported are Ag = Anopheles, Bt =
                        Bovine, Ce = Worm, Cf = Canine, Dm = Fly, Dr =
                        Zebrafish, Gg = Chicken, Hs = Human, Mm = Mouse, Mmu =
                        Rhesus, Pt = Chimp, Rn = Rat. See also: http://biocond
                        uctor.org/packages/release/BiocViews.html#___OrgDb
                        Currently unsupported are At, EcK12, EcSakai, Pf, Ss,
                        Xl.
  table                 Path to tabular input file(s).

optional arguments:
  -h, --help            show this help message and exit
  -i NAME, --id-column NAME
                        Column name of gene identifiers (Ensembl gene ids) to
                        query. (default: gene)
  -e NAME, --expression-column NAME
                        Column name of gene expression values. Merged into
                        output column 'GenesSignificantExpression'. (default:
                        None)
  -n NAME [NAME ...], --name-columns NAME [NAME ...]
                        Column name(s) for other values to be coerced as
                        lists. Merged into output column
                        'GenesSignificantNAME'. (default: None)
  -f FORMAT, --input-format FORMAT
                        Select the input file format. (default: tsv)
  -o OUTPUT_PREFIX, --output-prefix OUTPUT_PREFIX
                        Output file prefix. Can contain slashes for (sub)
                        folders. Prepended to input file base names without
                        suffix and sheet names. (default: topgo.)
  -x                    Also output a .xlsx version. (default: False)
  -v, --verbose         Be more verbose on what is done. (default: False)
  -w, --show-warnings   Show 'warnings()' after reading tables. (default:
                        False)
```

#### goplots

Create plots from gene ontology or other enrichment analysis.

Usage:

```
usage: goplots [-h] [-s SHEETS] [-d DELIMITER] [--col-term COL_TERM]
               [--col-enrich COL_ENRICH] [--col-nsign COL_NSIGN]
               [--col-nannotated COL_NANNOTATED] [--col-deg COL_DEG] [-f]
               [-n NTERMS] [-o OUTPUT_FOLDER]
               table [table ...]

positional arguments:
  table                 Input table

optional arguments:
  -h, --help            show this help message and exit
  -s SHEETS, --sheets SHEETS
                        Sheet for xl input (default: None)
  -d DELIMITER, --delimiter DELIMITER
                        For text input files, the column separator (default: )
  --col-term COL_TERM   Column name for term names/ids (default: termName)
  --col-enrich COL_ENRICH
                        Column name for (log) fold enrichment values (default:
                        foldEnrichment)
  --col-nsign COL_NSIGN
                        Column name for number of significant genes (default:
                        listHits)
  --col-nannotated COL_NANNOTATED
                        Column name for number of annotated genes (default:
                        Annotated)
  --col-deg COL_DEG     Column name for differential gene expression (e.g.
                        logFC-) values (default: geneExpression)
  -f, --format-deg      Replace missing gene expression values by mean of
                        values per term (default: False)
  -n NTERMS, --nterms NTERMS
                        Number of terms to plot (default: 20)
  -o OUTPUT_FOLDER, --output-folder OUTPUT_FOLDER
                        Output directory (default: goplots_output)
```

## License

> Copyright 2015
> * Jorge Boucas <JBoucas@age.mpg.de>
> * Sven Templer <templer@age.mpg.de>

See file `LICENSE`.
