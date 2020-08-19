#!/bin/bash

# TODO: Set General pipeline settings
# MPI-Age Settings
HOMESOURCE="source ~/.bashrc"
SLURMPARTITION="blade,himem,hugemem"
SING="singularity exec /beegfs/common/singularity/bioinformatics_software.v3.0.2.sif /bin/bash"

# reference genome
genome_folder=/beegfs/common/genomes/caenorhabditis_elegans/89/
gtf=${genome_folder}original.gtf
chromosomes=${genome_folder}original.toplevel.genome
fasta=${genome_folder}original.toplevel.fa


# project folder
top=/beegfs/group_bit/data/projects/departments/Bioinformatics/bit_pipelines/kallisto_deseq2/
raw=/beegfs/group_bit/data/projects/departments/Bioinformatics/bit_tmp_circRNA_pipeline/raw_data/linear_reads/
series='RNP6_circRNA'


# circRNA folder if circRNA counts are supposed to be included, 'None' otherwise
circRNA="/beegfs/group_bit/data/projects/departments/Bioinformatics/bit_pipelines/circRNA/dcc_out/RNP6/"

## If ERCC spikes are not used, comment the following lines until file sufixes
ref=${top}/references/
mkdir -p ${ref}

ercc_folder=/beegfs/common/genomes/ercc_spikes/                                                                                                                                                                                          
ercc_gtf=${ercc_folder}ERCC92.gtf
ercc_fa=${ercc_folder}ERCC92.fa 
ercc_sizes=${ercc_folder}ERCC92.genome

# Concatenating Ref genome with ERCC spikes
cat ${gtf} ${ercc_gtf} > ${ref}ERCC_REF.gtf
cat ${chromosomes} ${ercc_sizes} > ${ref}ERCC_REF.genome
cat ${fasta} ${ercc_fa} > ${ref}ERCC_REF.fa

gtf=${ref}ERCC_REF.gtf
chromosomes=${ref}ERCC_REF.genome
fasta=${ref}ERCC_REF.fa

# files sufixes
read1_sufix="-READ_1.fastq.gz"
read2_sufix="-READ_2.fastq.gz" # if single, write "none"

# "caenorhabditis elegans", "drosophila melanogaster", "mus musculus", "homo sapiens", "saccharomyces cerevisiae"
species="caenorhabditis elegans"

# http://www.ensembl.org/info/website/archives/index.html
biomart_host="http://jan2019.archive.ensembl.org/biomart/"

# datasets "celegans_gene_ensembl" "dmelanogaster_gene_ensembl","mmusculus_gene_ensembl","hsapiens_gene_ensembl", "scerevisiae_gene_ensembl"
biomart_dataset="celegans_gene_ensembl"

# for single end reads in kalisto
fragment_size=200
fragment_size_standard_deviation=30

# DAVID
daviddatabase="ENSEMBL_GENE_ID"
DAVIDUSER="franziska.metge@age.mpg.de"

# cytoscape
cytoscape_host="10.20.16.43"

# samples, if paired only READ is required
samplestable=$(readlink -f .)/sample_sheet.xlsx
# and excel file with one sheet of the form
#Files   daf2    treat
#S_001-F_And2-L____wt-___-cont-REP_1-READ_1.fastq.gz wt  control
#S_002-F_And2-L____wt-___-druA-REP_1-READ_1.fastq.gz wt  drug_A
#S_003-F_And2-L____wt-___-druB-REP_1-READ_1.fastq.gz wt  drug_B
#S_004-F_And2-L__daf2-___-cont-REP_1-READ_1.fastq.gz mut control
#S_005-F_And2-L__daf2-___-druA-REP_1-READ_1.fastq.gz mut drug_A
#S_006-F_And2-L__daf2-___-druB-REP_1-READ_1.fastq.gz mut drug_B
#S_007-F_And2-L____wt-___-cont-REP_2-READ_1.fastq.gz wt  control
#S_008-F_And2-L____wt-___-druA-REP_2-READ_1.fastq.gz wt  drug_A
#S_009-F_And2-L____wt-___-druB-REP_2-READ_1.fastq.gz wt  drug_B
#S_010-F_And2-L__daf2-___-cont-REP_2-READ_1.fastq.gz mut control
#S_011-F_And2-L__daf2-___-druA-REP_2-READ_1.fastq.gz mut drug_A
#S_012-F_And2-L__daf2-___-druB-REP_2-READ_1.fastq.gz mut drug_B
#S_013-F_And2-L____wt-___-cont-REP_3-READ_1.fastq.gz wt  control
#S_014-F_And2-L____wt-___-druA-REP_3-READ_1.fastq.gz wt  drug_A
#S_015-F_And2-L____wt-___-druB-REP_3-READ_1.fastq.gz wt  drug_B
#S_016-F_And2-L__daf2-___-cont-REP_3-READ_1.fastq.gz mut control
#S_017-F_And2-L__daf2-___-druA-REP_3-READ_1.fastq.gz mut drug_A
#S_018-F_And2-L__daf2-___-druB-REP_3-READ_1.fastq.gz mut drug_B

samples_treatment=$(readlink -f .)/samples_MasterTable.txt
# A tsv file with the following format
#	Treat
#S_001-F_And2-L____wt-___-cont-REP_1	control
#S_002-F_And2-L____wt-___-druA-REP_1	drug_A
#S_003-F_And2-L____wt-___-druB-REP_1	drug_B
#S_004-F_And2-L__daf2-___-cont-REP_1	control
#S_005-F_And2-L__daf2-___-druA-REP_1	drug_A
#S_006-F_And2-L__daf2-___-druB-REP_1	drug_B
#S_007-F_And2-L____wt-___-cont-REP_2	control
#S_008-F_And2-L____wt-___-druA-REP_2	drug_A
#S_009-F_And2-L____wt-___-druB-REP_2	drug_B
#S_010-F_And2-L__daf2-___-cont-REP_2	control
#S_011-F_And2-L__daf2-___-druA-REP_2	drug_A
#S_012-F_And2-L__daf2-___-druB-REP_2	drug_B
#S_013-F_And2-L____wt-___-cont-REP_3	control
#S_014-F_And2-L____wt-___-druA-REP_3	drug_A
#S_015-F_And2-L____wt-___-druB-REP_3	drug_B
#S_016-F_And2-L__daf2-___-cont-REP_3	control
#S_017-F_And2-L__daf2-___-druA-REP_3	drug_A
#S_018-F_And2-L__daf2-___-druB-REP_3	drug_B

#############################################################################

echo "Creating required folders"

mkdir -p ${top}

tmp=${top}/tmp/
logs=${top}/slurm_logs/
fqc=${top}/fastqc_output/
kalout=${top}/kallisto_output/
deg=${top}/${series}/deseq2_output/
mqc=${top}/${series}/multiqc_out/
qcp=${top}/${series}/qc_plots/
ind=${top}/kallisto_index/
scripts=$(readlink -f .)/
cdna_fasta=${ind}transcripts.norRNA.fa
kallisto_index=${ind}transcripts.norRNA.idx
feature_out=${top}${series}/featureCounts_output_kal/

mkdir -p ${tmp} ${logs} ${fqc} ${kalout} ${deg} ${mqc} ${qcp} ${ind} ${feature_out}

#############################################################################

#echo "Prepare gene id to gene name reference table"
echo "Extracting non rRNA transcripts"

############################################################################# 


if [ ! -e ${cdna_fasta} ] ; then 

grep -v -i rrna ${gtf} > ${tmp}no.rRNA.gtf
module load cufflinks
# echo "Done with grep"
# module load gff
# needs samtools faidx ${fasta}
# https://github.com/gpertea/gffread
gffread -w ${cdna_fasta} -g ${fasta} ${tmp}no.rRNA.gtf

fi

############################################################################# 

echo "Determine paired vs single end sequencing"
cd ${raw}
test_read_1=$(ls *${read1_sufix}* | head -n 1)
test_read_2=${test_read_1%${read1_sufix}}${read2_sufix}
echo ${test_read_1} 
if [ -e ${test_read_2} ] ; then seq="paired" ; echo ${test_read_2} ; else seq="single" ; test_read_2="" ; fi
echo "This is ${seq} end sequencing."

module load slurm

#############################################################################

echo "Starting FASTQC"

#############################################################################

cd ${raw} 
for file in $(ls *.fastq.gz | grep -v tmp ); do 

if [ ! -e ${fqc}${file%.fastq.gz}_fastqc.html ] ; then

rm -rf ${logs}fastqc.${file%.fastq.gz}.*.*

sbatch --partition $SLURMPARTITION << EOF
#!/bin/bash
#SBATCH --output ${logs}fastqc.${file%.fastq.gz}.%j.out
#SBATCH --error ${logs}fastqc.${file%.fastq.gz}.%j.err
#SBATCH -c 4
#SBATCH --job-name='fastqc'

${SING} << SIN
#!/bin/bash
${HOMESOURCE}

module load jdk fastqc
cd ${raw}

# FASTQC call
fastqc -t 4 -o ${fqc} ${file}

SIN
EOF

fi

done

#############################################################################

echo "Building required index."

#############################################################################

if [ ! -e ${kallisto_index} ] ; then  

rm -rf ${logs}kallisto.index.*.out

id=$(sbatch --partition $SLURMPARTITION --parsable << EOF
#!/bin/bash
#SBATCH --output ${logs}kallisto.index.%j.out
#SBATCH -c 4
#SBATCH --job-name="index"

${SING} << SIN
#!/bin/bash
${HOMESOURCE}

module load kallisto 
module load python/3.8.0
cd ${ind}

kallisto index -i ${kallisto_index} ${cdna_fasta}

SIN
EOF
)

srun -d afterok:${id} echo "Done building Kallisto's index."

fi

#############################################################################

echo "Determining strandness"

#############################################################################

if [ ! -e ${raw}${series}.strandness.txt ] ; then

rm -rf ${logs}strand.${series}.${test_read_1%${read1_sufix}}.*.out

python << PYOF
inGTF="${gtf}"
outbed="${ind}gene.model.bed"
bed=[]
with open(inGTF, "r") as fin:
    for line in fin:
        if line[0] != "#":
            l=line.split("\t")
            if l[2] == "gene":
                gene_id=l[8].split('gene_id "')[1].split('";')[0]
                bedline=[l[0],l[3],l[4],gene_id, ".",l[6]]
                bedline="\t".join(bedline)
                bed.append(bedline)
with open(outbed, "w") as fout:
    fout.write("\n".join(bed))
PYOF

id=$(sbatch --partition $SLURMPARTITION --parsable << EOF
#!/bin/bash
#SBATCH --output ${logs}strand.${series}.${test_read_1%${read1_sufix}}.%j.out
#SBATCH -c 18
#SBATCH --job-name="strand"

${SING} << SIN
#!/bin/bash
${HOMESOURCE}

module load kallisto
module load python/3.8.0
cd ${raw}

zcat ${test_read_1} | head -n 16000000 > tmp.${test_read_1}

if [ "${seq}" == "paired" ]
then 

echo ${seq}
zcat ${test_read_2} | head -n 16000000 > tmp.${test_read_2}
kallisto quant -t 18 -i ${kallisto_index} --genomebam -g ${gtf} -c ${chromosomes} -o $(pwd)/tmp.${test_read_1%${read1_sufix}} -b 100 tmp.${test_read_1} tmp.${test_read_2}
rm -rf tmp.${test_read_1} tmp.${test_read_2}
else

echo ${seq}
kallisto quant -t 18 -i ${kallisto_index} --genomebam -g ${gtf} -c ${chromosomes} -o $(pwd)/tmp.${test_read_1%${read1_sufix}} -b 100 --single -l ${fragment_size} -s ${fragment_size_standard_deviation} tmp.${test_read_1}
rm -rf tmp.${test_read_1}

fi

module load python/3.8.0
cd ${top}
unset PYTHONHOME
pip3 install virtualenv --user
rm -rf RSeQC  
virtualenv RSeQC
unset PYTHONUSERBASE
source RSeQC/bin/activate
pip install numpy
pip install RSeQC
infer_experiment.py -i ${raw}tmp.${test_read_1%${read1_sufix}}/pseudoalignments.bam -r ${ind}gene.model.bed > ${raw}${series}.infer_experiment.txt

python > ${raw}${series}.strandness.txt << PYOF
import sys
filein="${raw}${series}.infer_experiment.txt"
text=open(filein, "r").readlines()
text=[ s.split("\n")[0] for s in text ]
text=[ s for s in text if len(s) > 0 ]
strand={"fr-secondstrand":text[2].split(": ")[-1], \
       "fr-firststrand":text[3].split(": ")[-1]}
if float(strand["fr-secondstrand"]) > 0.75:
    print("fr-secondstrand")
    print("\n".join(text))
    sys.stdout.flush()
    sys.exit(0)
elif float(strand["fr-firststrand"]) > 0.75:
    print("fr-firststrand")
    print("\n".join(text))
    sys.stdout.flush()
    sys.exit(0)
elif ( float(strand["fr-firststrand"]) < 0.6 ) & ( float(strand["fr-secondstrand"]) < 0.6 ):
    print("unstranded")
    print("\n".join(text))
    sys.stdout.flush()
    sys.exit(0)
else:
    print("unable to determine strand")
    print("\n".join(text))
    sys.stdout.flush()
    sys.exit(1)
PYOF
SIN
EOF
)

echo "Waiting for strandness determination job ${id} to complete."
srun --partition $SLURMPARTITION -d afterok:${id} echo "Strandness determination complete."

rm -rf ${raw}tmp.${test_read_1%${read1_sufix}}  

fi


#############################################################################

strand=$(head -n 1 ${raw}${series}.strandness.txt) 
echo "This is ${strand} data. Starting Kallisto."

#############################################################################

# kallisto

#############################################################################

cd ${raw}

if [ "${strand}" == "fr-firststrand" ] ; then
    kallisto_strand="--rf-stranded"
elif [ "${strand}" == "fr-secondstrand" ] ; then
    kallisto_strand="--fr-stranded"
elif [ "${strand}" == "unstranded" ] ; then
    kallisto_strand=""
else
    exit
fi

ids=""

for f in $(ls *${read1_sufix}* | grep -v tmp ) ; do

if [ ! -e ${kalout}${f%${read1_sufix}}/pseudoalignments.bam ] ; then

rm -rf ${logs}kallisto.${f%${read1_sufix}}.*.out

ids=${ids}:$(sbatch --partition $SLURMPARTITION --parsable << EOF
#!/bin/bash
#SBATCH --output ${logs}kallisto.${f%${read1_sufix}}.%j.out
#SBATCH -c 18
#SBATCH --job-name="kallisto"

${SING} << SIN
#!/bin/bash
${HOMESOURCE}

module load kallisto

cd ${raw}

if [ "${seq}" == "paired" ]
then 

kallisto quant -t 18 -i ${kallisto_index} ${kallisto_strand} --genomebam -g ${gtf} -c ${chromosomes} -o ${kalout}${f%${read1_sufix}} -b 100 ${f} ${f%${read1_sufix}}${read2_sufix}

else

kallisto quant -t 18 -i ${kallisto_index} ${kallisto_strand} --genomebam -g ${gtf} -c ${chromosomes} -o ${kalout}${f%${read1_sufix}} -b 100 --single -l ${fragment_size} -s ${fragment_size_standard_deviation} ${f}

module load samtools
cd ${kalout}${f%${read1_sufix}}
samtools flagstat pseudoalignments.bam > flagstats.txt

fi
SIN
EOF
)

fi

done

ids=${ids}:$(sbatch --partition $SLURMPARTITION --parsable << EOF
#!/bin/bash
#SBATCH --output ${logs}empty.%j.out
EOF
)

echo "Waiting for kallisto jobs ${ids} to complete."
srun --partition $SLURMPARTITION -d afterok${ids} echo "Finished Kallisto. Starting featureCounts and multiqc."

#############################################################################

# featureCounts

#############################################################################

if [ "${strand}" == "fr-firststrand" ] ; then
    featureCounts_direction=2
elif [ "${strand}" == "fr-secondstrand" ] ; then
    featureCounts_direction=1
elif [ "${strand}" == "unstranded" ] ; then
    featureCounts_direction=0
else
    exit
fi

ids=""
cd ${kalout}
echo $(pwd)
for file in $(ls) ; do 

if [ ! -e ${feature_out}/${file}_biotype_counts_mqc.txt ] ; then

rm -rf ${logs}featureCount_kal_${file}.*.out
ids=${ids}:$(sbatch --partition $SLURMPARTITION --parsable << EOF
#!/bin/bash
#SBATCH --cpus-per-task=18
#SBATCH -o ${logs}featureCount_kal_${file}.%j.out
#SBATCH --job-name="featureCounts"

${SING} << SIN
#!/bin/bash
${HOMESOURCE}

module load subread

cd ${kalout}/${file}

echo "featureCounts -a $gtf -T 18 -g gene_id -o ${feature_out}/${file}_gene.featureCounts.txt -p -s $featureCounts_direction pseudoalignments.bam"
featureCounts -a $gtf -T 18 -g gene_id -o ${feature_out}/${file}_gene.featureCounts.txt -p -s $featureCounts_direction pseudoalignments.bam

echo "featureCounts -a $gtf -T 18 -g gene_biotype -o ${feature_out}/${file}_biotype.featureCounts.txt -p -s $featureCounts_direction pseudoalignments.bam"
featureCounts -a $gtf -T 18 -g gene_biotype -o ${feature_out}/${file}_biotype.featureCounts.txt -p -s $featureCounts_direction pseudoalignments.bam

cut -f 1,7 ${feature_out}/${file}_biotype.featureCounts.txt | tail -n +3 | (echo "$biotypes_header" && cat) >> ${feature_out}/${file}_biotype_counts_mqc.tmp.txt
grep -v '^\s' ${feature_out}/${file}_biotype_counts_mqc.tmp.txt > ${feature_out}/${file}_biotype_counts_mqc.txt

SIN
EOF
)

fi

done

ids=${ids}:$(sbatch --partition $SLURMPARTITION --parsable << EOF
#!/bin/bash
#SBATCH --output ${logs}empty.%j.out
EOF
)

echo "Waiting for FeatureCounts jobs ${ids} to complete."
srun --partition $SLURMPARTITION -d afterok${ids} echo "Starting Renaming FeatureCounts headers."

cat > ${scripts}featureCounts_Headers.py << PYOF
import os
import csv
import pandas as pd
import numpy as np

os.chdir("${feature_out}")
print(os.getcwd())

files_path="${feature_out}"
print(files_path)
files=os.listdir(files_path)
files=[s for s in files if "summary" not in s and "mqc" not in s]

for f in files:
    
    print(f)
    print(files_path+f)
    file=pd.read_csv(files_path+f, sep="\t")
    header=file.columns[0]
    file_=pd.read_csv(files_path+f, sep="\t", comment="#")
    sample_name=f.split(".",1)[0]
    file_=file_.rename(columns={'pseudoalignments.bam':sample_name})
    
    text_file = open(files_path+f, 'w')
    text_file.write(header+"\n")
    text_file.close()
     
    file_.to_csv(files_path+f, mode='a', header=True, sep="\t", index=None)
    
    file_sum=pd.read_csv(files_path+f+".summary", sep="\t")
    file_sum=file_sum.rename(columns={'pseudoalignments.bam':f})
    file_sum.to_csv(files_path+f+".summary", index=None, sep="\t")

PYOF

ids=""

rm -rf ${logs}fc_headers.*.out

ids=${ids}:$(sbatch --partition $SLURMPARTITION --parsable << EOF
#!/bin/bash
#SBATCH --output ${logs}fc_headers.%j.out
#SBATCH -c 4
#SBATCH --job-name="fc_headers"

${SING} << SIN
#!/bin/bash
${HOMESOURCE}

module purge
module load python/3.8.0
python3 ${scripts}featureCounts_Headers.py

SIN

EOF
)


ids=${ids}:$(sbatch --partition $SLURMPARTITION --parsable << EOF
#!/bin/bash
#SBATCH --output ${logs}empty.%j.out
EOF
)

echo "Waiting for Ranaming headers jobs ${ids} to complete."
srun --partition $SLURMPARTITION -d afterok${ids} echo "Finished Renaming FeatureCounts headers."


#############################################################################

# MultiQC

#############################################################################

if [ ! -e ${mqc}/multiqc_report.html ] ; then

rm -rf ${logs}multiqc.*.out
sbatch --partition $SLURMPARTITION -d afterok${ids} --parsable << EOF
#!/bin/bash
#SBATCH --cpus-per-task=1
#SBATCH -o ${logs}multiqc.%j.out
#SBATCH --job-name="multiqc"

${SING} << SIN
#!/bin/bash
${HOMESOURCE}

# Install multiqc
module load python/3.8.0
#pip install multiqc --user --ignore-installed

cd ${top}
pip3 install virtualenv --user
unset PYTHONHOME
virtualenv multiqc_kallisto
unset PYTHONUSERBASE
source multiqc_kallisto/bin/activate
pip install multiqc --ignore-installed  
multiqc ${fqc} ${feature_out} ${kalout} -f -o ${mqc}

SIN
EOF

fi

#############################################################################

# DESeq2 part 1

#############################################################################


if [ ! -e ${deg}deseq2.part1.Rdata ] ; then  

rm -rf ${logs}deseq2.*.out

cat > ${deg}deseq2.R << ROF
print(Sys.getenv("R_LIBS_USER"))
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager", repos="http://cloud.r-project.org/" )
if (!requireNamespace("tximportData", quietly = TRUE))
  BiocManager::install("tximportData", version = "3.8")
if (!requireNamespace("tximport", quietly = TRUE))
  BiocManager::install("tximport")
if (!requireNamespace("DESeq2", quietly = TRUE))
  BiocManager::install("DESeq2")
if (!requireNamespace("apeglm", quietly = TRUE))
  BiocManager::install("apeglm")
if (!requireNamespace("readr", quietly = TRUE))
  BiocManager::install("readr")
if (!requireNamespace("rhdf5", quietly = TRUE))
  BiocManager::install("rhdf5")

library(tximportData)
library(tximport)
library(DESeq2)
library(rhdf5)
library(readr)
library(apeglm)

# generate tx2gene for gene based counts (to compare to Berenice data)
tx2gene <- read_csv("${deg}tx2gene.csv")
tx2gene <- tx2gene[rowSums(is.na(tx2gene)) == 0,]
tx2gene <- droplevels(tx2gene)
save.image("${deg}deseq2.part1.Rdata")
sessionInfo()
ROF

cat > ${deg}read.template.py << PYOF
import pandas as pd
from biomart import BiomartServer
import itertools
import AGEpy as age
sfile="${samplestable}"
sdf=pd.read_excel(sfile)
fs=sdf.columns.tolist()[0]
sdf[fs]=sdf[fs].apply(lambda x: x.split("${read1_sufix}")[0] )
sdf.index=sdf[fs].tolist()
sdf=sdf.drop([fs],axis=1)
cols=sdf.columns.tolist()
mods=[ [x,y] for x in cols for y in cols ]
interactions=[]
for c in mods:
    if c not in interactions:
        if [c[1], c[0]] not in interactions:
            if c[0] != c[1]:
                interactions.append(c)

single_models=cols

textout=[]
for m in single_models:
    variants=[ s for s in cols if s != m ] 
    tmp=sdf.copy()
    tmp["_group_"]=m
    for c in variants:
        tmp["_group_"]=tmp["_group_"]+"."+c+"_"+tmp[c]
    for g in list(set(tmp["_group_"].tolist())):
        outdf=tmp[tmp["_group_"]==g]
        model_data=list(set(outdf[m].tolist()))
        model_pairs=[ list(set([x,y])) for x in model_data for y in model_data ]
        model_pairs_=[ pair for pair in model_pairs if len(pair) > 1 ]
        model_pairs=[]
        for pair in model_pairs_:
            if pair not in model_pairs:
                model_pairs.append(pair)
        for pair in model_pairs:
            outdf_=outdf[outdf[m].isin(pair)]
            outdf_=outdf_.drop(["_group_"],axis=1)
            coef=outdf_[m].tolist()
            ref=coef[0]
            target=[ t for t in coef if t != ref ][0]
            coef=m+"_"+target+"_vs_"+ref
            filename=coef+g.split(m)[-1]
            outdf_.to_csv("${deg}"+filename+".input.tsv", sep="\t")
            text=filename+"\t"+m+"\t"+g+"\t"+ref+"\t"+coef
            textout.append(text)
        
with open("${deg}models.txt", "w") as mout:
     mout.write("\n".join(textout) + "\n")

GTF=age.readGTF("${gtf}")
GTF["gene_id"]=age.retrieve_GTF_field(field="gene_id",gtf=GTF)
GTF["transcript_id"]=age.retrieve_GTF_field(field="transcript_id",gtf=GTF)
GTF["gene_biotype"]=age.retrieve_GTF_field(field="gene_biotype",gtf=GTF)
tx2gene=GTF[["transcript_id","gene_id"]].drop_duplicates()
tx2gene.columns=["TXNAME","GENEID"]
tx2gene[["TXNAME","GENEID"]].to_csv("${deg}"+"tx2gene.csv", quoting=1, index=None)

attributes=["ensembl_gene_id","go_id","name_1006"]
server = BiomartServer( "${biomart_host}" )
organism=server.datasets["${biomart_dataset}"]
response=organism.search({"attributes":attributes})
response=response.content.decode().split("\n")
response=[s.split("\t") for s in response ]
bio_go=pd.DataFrame(response,columns=attributes)
bio_go.to_csv("${deg}"+"annotated/biotypes_go_raw.txt", index=None, sep="\t")
bio_go.columns = ["ensembl_gene_id","GO_id","GO_term"]

def CombineAnn(df):
    return pd.Series(dict(ensembl_gene_id = "; ".join([ str(s) for s in list(set(df["ensembl_gene_id"]))  if str(s) != "nan" ] ) ,\
                       GO_id = "; ".join([ str(s) for s in list(set(df["GO_id"])) if str(s) != "nan" ] ) ,\
                       GO_term = "; ".join([ str(s) for s in list(set(df["GO_term"])) if str(s) != "nan" ] ) ,\
                      ) )

bio_go=bio_go.groupby(by="ensembl_gene_id", as_index=False).apply(CombineAnn)

bio_go.reset_index(inplace=True, drop=True)
#bio_go.to_csv("${deg}"+"annotated/biotypes_go.txt", sep= "\t", index=None)

GTF=GTF[["gene_id","gene_biotype"]].drop_duplicates()
GTF.columns=["ensembl_gene_id","gene_biotype"]
bio_go=pd.merge(GTF,bio_go,on=["ensembl_gene_id"],how="outer")
bio_go.to_csv("${deg}"+"annotated/biotypes_go.txt", sep= "\t", index=None)

PYOF

# see tx2gene.sh

mkdir -p ${deg}annotated

id=$(sbatch --partition $SLURMPARTITION --parsable << EOF
#!/bin/bash
#SBATCH --output ${logs}deseq2.%j.out
#SBATCH -c 4
#SBATCH --job-name="deseq2"

${SING} << SIN
#!/bin/bash
${HOMESOURCE}
module purge

module load python/3.8.0
pip3 install pandas --user
pip3 install xlrd --user
pip3 install biomart --user
pip3 install git+https://github.com/mpg-age-bioinformatics/AGEpy.git --user

python3 ${deg}read.template.py

module load rlang
Rscript ${deg}deseq2.R

SIN
EOF
)


echo "Waiting for DESeq2 part 1 job ${id} to complete."
srun --partition $SLURMPARTITION -d afterok:${id} echo "Finished DESeq2 part 1. Starting DESeq2 pair wise comparisons."

fi

#############################################################################

# DESeq2 part 2 (pair-wise comparisons)

#############################################################################

ids=""

while read line ; do
#filename+"\t"+m+"\t"+g+"\t"+ref+"\t"+coef
infile=$(echo ${line} | awk '{ print $1 }').input.tsv
model=$(echo ${line} | awk '{ print $2 }')
ref=$(echo ${line} | awk '{ print $4 }')
coef=$(echo ${line} | awk '{ print $5 }')
out=${deg}$(echo ${line} | awk '{ print $1 }')
outs=$(echo ${line} | awk '{ print $1 }')

cat > ${out}.deseq2.R << EOF 
load("${deg}deseq2.part1.Rdata")
library(tximportData)
library(tximport)
library(DESeq2)
library(rhdf5)
library(readr)
library(apeglm)

# filein
sampleTable<-read.delim2("${deg}${infile}",sep = "\t", row.names = 1)
samples<-row.names(sampleTable)

# dir
dir<-"${kalout}"
files<-file.path(dir, samples, "abundance.h5")
names(files)<-samples
txi <- tximport(files, type = "kallisto", txOut = TRUE)
txi <- summarizeToGene(txi, tx2gene = tx2gene)

# model
dds <- DESeqDataSetFromTximport(txi, sampleTable, ~${model} )

# if circRNA file is present, add circRNA counts to dds, else proceed
# here, get count table, 
# add circRNA counts
circRNA_folder="${circRNA}"
if(circRNA_folder != "None"){
  # gene count table
  gene_count = as.data.frame(counts(dds, normalized=FALSE))

  # circRNA table
  circRNA = read.delim('${circRNA}/CircRNACount', check.names = FALSE, as.is = TRUE )
  coord = read.delim('${circRNA}/CircCoordinates', check.names = FALSE, as.is = TRUE)
  
  # reformat circRNA header
  names(circRNA) <- gsub('.Chimeric.out.junction', '', names(circRNA))
  row.names(circRNA) <- paste0('circ_', coord[, 'Chr'], ':', coord[,'Start'], '-', coord[,'End'], '|', coord[,'Strand'], '|', coord[,'Gene'])
  
  circRNA = circRNA[,names(gene_count)]
  
  # add circRNAs to gene count table
  gene_count = rbind(gene_count, circRNA)
  
  # create DESeqDataSetFromMatrix
  dds <- DESeqDataSetFromMatrix(gene_count, sampleTable, ~${model})
} 

# model
dds\$${model} <- relevel(dds\$${model}, ref = "${ref}")

keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
dds <- DESeq(dds)
res.counts<-counts(dds, normalized=TRUE)

# coef
resLFC <- lfcShrink(dds,  coef="${coef}", type="apeglm")
counts.dge<-merge(res.counts,resLFC, by=0, all=FALSE)
row.names(counts.dge)<-counts.dge\$Row.names
counts.dge<-counts.dge[, !(colnames(counts.dge) %in% c("Row.names"))]
counts.dge <- counts.dge[order(counts.dge\$pvalue),]

# file out
write.table(counts.dge, "${out}.results.tsv", sep="\t")

EOF

if [ ! -e ${out}.results.tsv ] ; then

rm -rf ${logs}deseq2.${outs}.*.log

ids=${ids}:$(sbatch --partition $SLURMPARTITION --parsable << EOF
#!/bin/bash
#SBATCH --output ${logs}deseq2.${outs}.%j.out
#SBATCH -c 4
#SBATCH --job-name="deseq2"

${SING} << SIN
#!/bin/bash
${HOMESOURCE}

module purge
module load rlang
which R

Rscript ${out}.deseq2.R

SIN
EOF
)

fi

done < ${deg}models.txt

ids=${ids}:$(sbatch --partition $SLURMPARTITION --parsable << EOF
#!/bin/bash
#SBATCH --output ${logs}empty.%j.out
EOF
)

echo "Waiting for DESeq2 part 2 jobs ${ids} to complete."
srun --partition $SLURMPARTITION -d afterok${ids} echo "Finished DESeq2 part 2. Starting Master Table generation."

#############################################################################

# DESeq2 MASTER TABLE

#############################################################################

if [ ! -e ${deg}all_results_stats.tsv ] ; then  

rm -rf ${logs}res_counts.*.out

cat > ${scripts}res_counts.R << EOF 
setwd("${top}")

library("readxl")
library(dplyr)
#library(factoextra)
library(ggplot2)
library(gridExtra)
library(tximportData)
library(tximport)
library(DESeq2)
library(rhdf5)
library(readr)
library(apeglm)
library(WriteXLS)
library(gplots)

load("${deg}deseq2.part1.Rdata")

sampleTable<-read.delim2("${samples_treatment}",sep = "\t", row.names = 1)
samples<-row.names(sampleTable)

# dir
dir <- "${kalout}"
files<-file.path(dir, samples, "abundance.h5")
names(files)<-samples

txi <- tximport(files, type = "kallisto", txOut = TRUE)
txi <- summarizeToGene(txi, tx2gene = tx2gene)

# model
dds <- DESeqDataSetFromTximport(txi, sampleTable, ~ Treat )


# if circRNA file is present, add circRNA counts to dds, else proceed
# here, get count table, 
# add circRNA counts
circRNA_folder="${circRNA}"
if(circRNA_folder != "None"){
  # gene count table
  gene_count = as.data.frame(counts(dds, normalized=FALSE))

  # circRNA table
  circRNA = read.delim('${circRNA}/CircRNACount', check.names = FALSE, as.is = TRUE )
  coord = read.delim('${circRNA}/CircCoordinates', check.names = FALSE, as.is = TRUE)
  
  # reformat circRNA header
  names(circRNA) <- gsub('.Chimeric.out.junction', '', names(circRNA))
  row.names(circRNA) <- paste0('circ_', coord[, 'Chr'], ':', coord[,'Start'], '-', coord[,'End'], '|', coord[,'Strand'], '|', coord[,'Gene'])
  
  circRNA = circRNA[,names(gene_count)]
  
  # add circRNAs to gene count table
  gene_count = rbind(gene_count, circRNA)
  
  # create DESeqDataSetFromMatrix
  dds <- DESeqDataSetFromMatrix(gene_count, sampleTable, ~${model})
} 

dds <- estimateSizeFactors(dds)
res.counts <- counts(dds, normalized=TRUE)

res_counts <- as.data.frame(res.counts)

WriteXLS(res_counts, ExcelFileName = "${deg}all_res_counts.xlsx", row.names = TRUE, col.names = TRUE)

result_tables <- list.files('${deg}', pattern = '.results.tsv')

for(f in result_tables){
  tmp <- read.delim(paste0('${deg}', f))
  tmp <- tmp[,c('log2FoldChange', 'pvalue', 'padj')]
  names(tmp) <- paste(names(tmp), gsub('.results.tsv', '', gsub('Group_', '', f)), sep = '.')
  res_counts <- merge(res_counts, tmp, by.x = 'row.names', by.y = 'row.names')
  print(nrow(res_counts))
  names(res_counts)[names(res_counts) == "Row.names"] <- "ensembl_gene_id"
  print(length(res_counts[,'ensembl_gene_id']))
  row.names(res_counts) <- res_counts[,'ensembl_gene_id']
  res_counts <- res_counts[, !duplicated(colnames(res_counts))]
}

write.table(res_counts, '${deg}all_results_stats.tsv', sep = '\t', quote = F, row.names = F)
WriteXLS(res_counts, ExcelFileName = "${deg}all_results_stats.xlsx", row.names = FALSE, col.names = TRUE)

EOF

id=$(sbatch --partition $SLURMPARTITION --parsable << EOF
#!/bin/bash
#SBATCH --output ${logs}res_counts.%j.out
#SBATCH -c 4
#SBATCH --job-name="res_counts"

#${SING} << SIN
##!/bin/bash
#${HOMESOURCE}

module load rlang

Rscript ${scripts}res_counts.R

#SIN

EOF
)


echo "Waiting for DESeq2 Master Table  job ${id} to complete."
srun --partition $SLURMPARTITION -d afterok:${id} echo "Finished generating Master Table. Starting Annotation."

fi

#############################################################################

# DESeq2 part 3 (annotations)

#############################################################################

if [ ! -e ${deg}annotated/significant.xlsx ] ; then

rm -rf ${logs}annotate.deseq2.*.out

cat > ${deg}annotate.deseq2.py << PYOF
import pandas as pd
import os
import AGEpy as age

GTF=age.readGTF("${gtf}")
GTF["gene_id"]=age.retrieve_GTF_field(field="gene_id",gtf=GTF)
GTF["gene_name"]=age.retrieve_GTF_field(field="gene_name",gtf=GTF)
id_name=GTF[["gene_id","gene_name"]].drop_duplicates()
id_name.reset_index(inplace=True, drop=True)
id_name.columns=["ensembl_gene_id","gene_name"]

#id_name=pd.read_table("${ind}cdna.norRNA.tsv")
#id_name=id_name[["gene_id","gene_symbol","description"]]
#id_name.columns=["ensembl_gene_id","gene_name","description"]
#id_name=id_name.drop_duplicates()
#id_name.reset_index(inplace=True,drop=True)

bio_go=pd.read_csv("${deg}"+"annotated/biotypes_go.txt", sep="\t")

deg_files=os.listdir("${deg}")
deg_files=[ s for s in deg_files if "results.tsv" in s ]
i=1
s=[]
dfs={}
for f in deg_files:
    df=pd.read_table("${deg}"+f)
    df=pd.merge(id_name,df,left_on=["ensembl_gene_id"],right_index=True, how="right") # change to gene_id
    df=pd.merge(df,bio_go,on=["ensembl_gene_id"],how="left")
    df=df.sort_values(by=["padj"],ascending=True)
    df.to_csv("${deg}"+"annotated/"+f, sep="\t",index=None)
    df.to_excel("${deg}"+"annotated/"+f.replace('.tsv', '.xlsx'), index=None)
    n=f.split(".results.tsv")[0]
    s.append([i,n])
    df=df[df["padj"]<0.05]
    df.reset_index(inplace=True, drop=True)
    dfs[i]=df
    i=i+1
sdf=pd.DataFrame(s,columns=["sheet","comparison"])
EXC=pd.ExcelWriter("${deg}"+"annotated/significant.xlsx")
sdf.to_excel(EXC,"summary",index=None)
for k in list(dfs.keys()):
    dfs[k].to_excel(EXC, str(k),index=None)
EXC.close()

mt=pd.read_csv("${deg}"+"all_results_stats.tsv", sep="\t")

mt_ann=pd.merge(id_name,mt,on=["ensembl_gene_id"], how="right")
mt_ann=pd.merge(mt_ann,bio_go,on=["ensembl_gene_id"],how="left")

mt_ann.to_csv("${deg}"+"annotated/masterTable_annotated.tsv", sep="\t",index=None)
mt_ann.to_excel("${deg}"+"annotated/masterTable_annotated.xlsx", index=None)

PYOF

id=$(sbatch --partition $SLURMPARTITION --parsable << EOF
#!/bin/bash
#SBATCH --output ${logs}annotate.deseq2.%j.out
#SBATCH -c 4
#SBATCH --job-name="adeseq2"

${SING} << SIN
#!/bin/bash
${HOMESOURCE}
module purge

module load python/3.8.0

python3 ${deg}annotate.deseq2.py

SIN
EOF
)

echo "Waiting for DESeq2 annotation part 1 job ${id} to complete."
srun --partition $SLURMPARTITION -d afterok:${id} echo "Finished DESeq2 annotation. Starting  enrichment analysis."

fi

#############################################################################

# DAVID

#############################################################################

cd ${deg}annotated

ids=""

for f in $(ls *results.tsv); do

    if [ ! -e ${deg}/annotated/${f%results.tsv}DAVID.tsv ] ; then   

    PY=${f%results.tsv}deseq2.p2.py
    cat > ${deg}${PY} << PYOF
import pandas as pd
import AGEpy as age
import os 
import sys

deseq2="${deg}/annotated/"
f="${f}"

if os.path.isfile(deseq2+f.replace("results.tsv","DAVID.xlsx")):
    sys.exit()

df=pd.read_csv(deseq2+f,sep="\t")
df=df[df["padj"]<0.05]

dics=df[["ensembl_gene_id","gene_name","log2FoldChange"]]
dics["ensembl_gene_id"]=dics["ensembl_gene_id"].apply(lambda x: x.upper())
dics.index=dics["ensembl_gene_id"].tolist()
names_dic=dics[["gene_name"]].to_dict()["gene_name"]
exp_dic=dics[["log2FoldChange"]].to_dict()["log2FoldChange"]

genes=df["ensembl_gene_id"].tolist()

DAVID=age.DAVIDenrich(database="${daviddatabase}",\
                categories="GOTERM_BP_FAT,GOTERM_CC_FAT,GOTERM_MF_FAT,KEGG_PATHWAY,PFAM,PROSITE,GENETIC_ASSOCIATION_DB_DISEASE,OMIM_DISEASE",\
               user="${DAVIDUSER}",\
               ids=genes, verbose=True)

if type(DAVID) == type(pd.DataFrame()):
    #for c in DAVID.columns.tolist():
    #    DAVID[c]=DAVID[c].apply(lambda x: x.decode())
    #print(DAVID.head(),DAVID["geneIds"].tolist(),names_dic )
    DAVID["genes name"]=DAVID["geneIds"].apply(lambda x: ", ".join([ str(names_dic[s.upper()]) for s in x.split(", ") ] ) )
    DAVID["log2fc"]=DAVID["geneIds"].apply(lambda x: ", ".join([ str(exp_dic[s.upper()]) for s in x.split(", ") ] ) )

    DAVID.to_csv(deseq2+f.replace("results.tsv","DAVID.tsv"), sep="\t", index=None)
    EXC=pd.ExcelWriter(deseq2+f.replace("results.tsv","DAVID.xlsx"))
    df.to_excel(EXC,"genes",index=None)
    for cat in DAVID["categoryName"].tolist():
        tmp=DAVID[DAVID["categoryName"]==cat]
        tmp.to_excel(EXC,cat,index=None)
    EXC.close()
PYOF

rm -rf ${logs}${PY%py}*.out

ids=${ids}:$(sbatch --partition $SLURMPARTITION --parsable << EOF
#!/bin/bash
#SBATCH --output ${logs}${PY%py}%j.out
#SBATCH -c 4
#SBATCH --job-name="deseq2.p2"
module load shifter

${SING} << SIN
#!/bin/bash
${HOMESOURCE}
module purge

module load python/3.8.0

python3 ${deg}${PY}

SIN
EOF
)

fi

done

ids=${ids}:$(sbatch --partition $SLURMPARTITION --parsable << EOF
#!/bin/bash
#SBATCH --output ${logs}empty.%j.out
EOF
)

echo "Waiting for DAVID annotation  job ${ids} to complete."
srun --partition $SLURMPARTITION -d afterok${ids} echo "Finished DAVID annotation. Starting cytoscape and drawing DAIVD plots."


#############################################################################

# Draw DAVID plots

#############################################################################

cd ${scripts}

wget https://raw.githubusercontent.com/mpg-age-bioinformatics/htseq-tools/master/david_to_cellplot.R
chmod 775 david_to_cellplot.R

ids=""

ids=${ids}:$(sbatch --partition $SLURMPARTITION --parsable << EOF
#!/bin/bash
#SBATCH --output ${logs}david_plot%j.out
#SBATCH -c 2
#SBATCH --job-name="david_plot"

for f in ${deg}/annotated/*DAVID.tsv ; do

${SING} << SIN
#!/bin/bash
${HOMESOURCE}

module load rlang

Rscript ${scripts}david_to_cellplot.R \${f} tsv GOTERM_BP_FAT 15
Rscript ${scripts}david_to_cellplot.R \${f} tsv KEGG_PATHWAY 15


SIN
done
EOF
)

ids=${ids}:$(sbatch --partition $SLURMPARTITION --parsable << EOF
#!/bin/bash
#SBATCH --output ${logs}empty.%j.out
EOF
)

#echo "Waiting for DAVID plots job ${ids} to complete."
#srun --partition $SLURMPARTITION -d afterok:${ids} echo "Finished pipeline. Done"


#############################################################################

# PPIs with String and Cytoscape

#############################################################################

rm -rf ${logs}PPIs.*.out

cat > ${deg}PPIs.py << PYOF
import pandas as pd
import numpy as np
import AGEpy as age
import sys
import os
from py2cytoscape import cyrest
from py2cytoscape.cyrest.base import *
import paramiko
from time import sleep
import matplotlib
import matplotlib.pyplot as plt
import tempfile

################# in values ################################

fin=sys.argv[1]
species="${species}"
host="${cytoscape_host}"
biomarthost="${biomart_host}"

###########################################################

python_output="/".join(fin.split("/")[:-1])
target=fin.replace("results.tsv","cytoscape")

if os.path.isfile(target+".cys"):
    sys.exit()


taxons={"caenorhabditis elegans":"6239","drosophila melanogaster":"7227",\
       "mus musculus":"10090","homo sapiens":"9606", "saccharomyces cerevisiae": "4932"}
tags={"caenorhabditis elegans":"CEL","drosophila melanogaster":"DMEL",\
       "mus musculus":"MUS","homo sapiens":"HSA"}

taxon_id=taxons[species]

### ATTENTION ### if you are using yeast, you will need to uncomment the follwing lines 
organismtag=tags[species]

if not os.path.isfile(python_output+"/homdf.txt"):
    print("Could not find ageing evidence table. Using biomart to create one.")
    sys.stdout.flush()
    homdf,HSA,MUS,CEL,DMEL=age.FilterGOstring(host=biomarthost)
    homdf.to_csv(python_output+"/homdf.txt", index=None,sep="\t")
else:
    print("Found existing ageing evidence table.")
    sys.stdout.flush()
homdf=pd.read_csv(python_output+"/homdf.txt", sep="\t")
aging_genes=homdf[[organismtag+"_ensembl_gene_id","evidence"]].dropna()
aging_genes=aging_genes[aging_genes[organismtag+"_ensembl_gene_id"]!="None"]
aging_genes=aging_genes[organismtag+"_ensembl_gene_id"].tolist()
### till here

dfin=pd.read_csv(fin, sep="\\t")

cytoscape=cyrest.cyclient(host=host)
cytoscape.version()
cytoscape.session.new()
cytoscape.vizmap.apply(styles="default")

# Annotate aging evindence
def CheckEvidence(x,aging_genes=aging_genes):
    if x in aging_genes:
        res="aging_gene"
    else:
        res="no"
    return res

### also comment this line
dfin["evidence"]=dfin["ensembl_gene_id"].apply(lambda x:CheckEvidence(x) )

dfin["baseMean"]=dfin["baseMean"].apply(lambda x: np.log10(x))

qdf=dfin[dfin["padj"]<0.05]
qdf=qdf.sort_values(by=["padj"],ascending=True)
query_genes=qdf["ensembl_gene_id"].tolist()[:1000]
limit=int(len(query_genes)*.25)
response=api("string", "protein query",\
                      {"query":",".join(query_genes),\
                       "cutoff":str(0.4),\
                       "species":species,\
                       "limit":str(limit),\
                       "taxonID":taxon_id},\
                       host=host, port="1234")

cytoscape.layout.force_directed(defaultSpringCoefficient=".000004", defaultSpringLength="5")
defaults_dic={"NODE_SHAPE":"ellipse",\
               "NODE_SIZE":"60",\
               "NODE_FILL_COLOR":"#AAAAAA",\
               "EDGE_TRANSPARENCY":"120"}
defaults_list=cytoscape.vizmap.simple_defaults(defaults_dic)

NODE_LABEL=cytoscape.vizmap.mapVisualProperty(visualProperty="NODE_LABEL",\
                                              mappingType="passthrough",\
                                              mappingColumn="display name")

cytoscape.vizmap.create_style(title="dataStyle",\
                              defaults=defaults_list,\
                              mappings=[NODE_LABEL])
sleep(2)
cytoscape.vizmap.apply(styles="dataStyle")


uploadtable=dfin[dfin["padj"]<0.05][["ensembl_gene_id","baseMean","log2FoldChange","evidence"]].dropna()
# uploadtable=dfin[dfin["padj"]<0.05][["ensembl_gene_id","baseMean","log2FoldChange"]].dropna() ### use this line if you are using yeast

cytoscape.table.loadTableData(uploadtable,df_key="ensembl_gene_id",table_key_column="query term")
sleep(10)


cmap = matplotlib.cm.get_cmap("bwr")
norm = matplotlib.colors.Normalize(vmin=-4, vmax=4)
min_color=matplotlib.colors.rgb2hex(cmap(norm(-4)))
center_color=matplotlib.colors.rgb2hex(cmap(norm(0)))
max_color=matplotlib.colors.rgb2hex(cmap(norm(4)))  

NODE_FILL_COLOR=cytoscape.vizmap.mapVisualProperty(visualProperty="NODE_FILL_COLOR",mappingType="continuous",\
                                                   mappingColumn="log2FoldChange",\
                                                 lower=[-4,min_color],\
                                                   center=[0.0,center_color],\
                                                   upper=[4,max_color])

### do not do this if you are using yeast ...
# apply diamond shape and increase node size to nodes with aging evidence
NODE_SHAPE=cytoscape.vizmap.mapVisualProperty(visualProperty="NODE_SHAPE",mappingType="discrete",mappingColumn="evidence",\
                                              discrete=[ ["aging_gene","no"], ["DIAMOND", "ellipse"] ])
NODE_SIZE=cytoscape.vizmap.mapVisualProperty(visualProperty="NODE_SIZE",mappingType="discrete",mappingColumn="evidence",\
                                             discrete=[ ["aging_gene","no"], ["100.0","60.0"] ])
###
cytoscape.vizmap.update_style("dataStyle",mappings=[NODE_SIZE,NODE_SHAPE,NODE_FILL_COLOR])
# cytoscape.vizmap.update_style("dataStyle",mappings=[NODE_FILL_COLOR]) # if using yeast

cytoscape.vizmap.apply(styles="dataStyle")

network = "current"
namespace='default'
PARAMS=set_param(["columnList","namespace","network"],["SUID",namespace,network])
network=api(namespace="network", command="get attribute",PARAMS=PARAMS, host=host,port='1234',version='v1')
network=int(network[0]["SUID"])


basemean = cytoscape.table.getTable(table="node",columns=["baseMean"], network = network)

min_NormInt = min(basemean.dropna()["baseMean"].tolist())
max_NormInt = max(basemean.dropna()["baseMean"].tolist())
cent_NormInt = np.mean([min_NormInt,max_NormInt])

cmap = matplotlib.cm.get_cmap("Reds")
norm = matplotlib.colors.Normalize(vmin=min_NormInt, vmax=max_NormInt)
min_color=matplotlib.colors.rgb2hex(cmap(norm(np.mean([min_NormInt,max_NormInt]))))
center_color=matplotlib.colors.rgb2hex(cmap(norm(cent_NormInt)))
max_color=matplotlib.colors.rgb2hex(cmap(norm(max_NormInt)))  

NODE_BORDER_PAINT=cytoscape.vizmap.mapVisualProperty(visualProperty="NODE_BORDER_PAINT",\
                                                     mappingType="continuous",\
                                                     mappingColumn="baseMean",\
                                                     lower=[min_NormInt,min_color],\
                                                     center=[np.mean([min_NormInt,max_NormInt]),center_color],\
                                                     upper=[max_NormInt,max_color])
cytoscape.vizmap.update_style("dataStyle",mappings=[NODE_BORDER_PAINT])

NODE_BORDER_WIDTH=cytoscape.vizmap.mapVisualProperty(visualProperty="NODE_BORDER_WIDTH",\
                                                     mappingType="continuous",\
                                                     mappingColumn="baseMean",\
                                                     lower=[min_NormInt,2],\
                                                     center=[np.mean([min_NormInt,max_NormInt]),4],\
                                                     upper=[max_NormInt,8])
cytoscape.vizmap.update_style("dataStyle",mappings=[NODE_BORDER_WIDTH])

cytoscape.vizmap.apply(styles="dataStyle")
cytoscape.network.rename(name="main String network")

cytoscape.network.select(edgeList="all", extendEdges="true")
cytoscape.network.create(source="current",nodeList="selected")
cytoscape.network.rename(name="main String network (edges only)")

cytoscape.network.set_current(network="main String network (edges only)")

log2FoldChange = cytoscape.table.getTable(table="node",columns=["log2FoldChange"])
if int(len(log2FoldChange)*.10) > 0:
    log2FoldChange["log2FoldChange"]=log2FoldChange["log2FoldChange"].apply(lambda x: abs(x))
    log2FoldChange=log2FoldChange.sort_values(by=["log2FoldChange"],ascending=False)
    top_nodes=log2FoldChange.index.tolist()[:int(len(log2FoldChange)*.10)]

    cytoscape.network.set_current(network="main String network (edges only)")
    cytoscape.network.select(nodeList="name:"+",".join(top_nodes))
    cytoscape.network.select(firstNeighbors="any",network="current")
    sleep(5)
    cytoscape.network.create(source="current",nodeList="selected")
    cytoscape.network.rename(name="top "+str(int(len(log2FoldChange)*.10))+" changed firstNeighbors")

def MAKETMP():
    (fd, f) = tempfile.mkstemp()
    f="/tmp/"+f.split("/")[-1]
    return f

cys=MAKETMP()
cyjs=MAKETMP()
main_png=MAKETMP()
main_pdf=MAKETMP()
edg_png=MAKETMP()
edg_pdf=MAKETMP()
neig_png=MAKETMP()
neig_pdf=MAKETMP()

cytoscape.session.save_as(session_file=cys)
cytoscape.network.export(options="CYJS",OutputFile=cyjs)

cytoscape.network.set_current(network="main String network")
cytoscape.network.deselect(edgeList="all",nodeList="all")
cytoscape.view.export(options="PNG",outputFile=main_png)
cytoscape.view.export(options="PDF",outputFile=main_pdf)

cytoscape.network.set_current(network="main String network (edges only)")
cytoscape.network.deselect(edgeList="all",nodeList="all")
cytoscape.view.export(options="PNG",outputFile=edg_png)
cytoscape.view.export(options="PDF",outputFile=edg_pdf)

if int(len(log2FoldChange)*.10) > 0:
    cytoscape.network.set_current(network="top "+str(int(len(log2FoldChange)*.10))+" changed firstNeighbors")
    cytoscape.network.deselect(edgeList="all",nodeList="all")
    sleep(5)
    cytoscape.view.export(options="PNG",outputFile=neig_png)
    cytoscape.view.export(options="PDF",outputFile=neig_pdf)

ssh = paramiko.SSHClient()
ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())
ssh.connect(host)

ftp_client=ssh.open_sftp()

for f, extension, local in zip([cys,cyjs,main_png,main_pdf,edg_png,edg_pdf,neig_png,neig_pdf],\
                                [".cys",".cyjs",".png",".pdf",".png",".pdf",".png",".pdf" ],\
                                [target+".cys",target+".cyjs",target+".main.png",target+".main.pdf",\
                                target+".main.edges.png",target+".main.edges.pdf",\
                                target+".topFirstNeighbors.png",target+".topFirstNeighbors.pdf"]):

    try:
        ftp_client.get(f+extension,local)
        ssh_stdin, ssh_stdout, ssh_stderr = ssh.exec_command("rm "+f+extension )
    except:
        print("No "+local)
        sys.stdout.flush()
print("Done with cytoscape.")

PYOF

cd ${deg}annotated
list_of_results="$(ls  *.results.tsv)"

sbatch --partition $SLURMPARTITION --parsable << EOF
#!/bin/bash
#SBATCH --output ${logs}PPIs.%j.out
#SBATCH -c 4
#SBATCH --job-name="ppis"

cd ${deg}annotated
for f in *.results.tsv ; do
    if [ ! -e ${deg}/annotated/${f%results.tsv}cytoscape.cys ] ; then

${SING} << SIN
#!/bin/bash
${HOMESOURCE}
module purge

module load python/3.8.0

python3 ${deg}PPIs.py ${deg}annotated/\${f}

SIN
fi
done
EOF


#############################################################################

# Draw QC plots

#############################################################################

if [ ! -e ${qcp}pca_all_samples.pdf ] ; then

rm -rf ${logs}qc_plots.*.out

cat > ${scripts}QC_plots.py << PYOF
import os
import AGEpy as age
import pandas as pd
import numpy as np
import matplotlib
import seaborn as sns

import matplotlib.pyplot as plt
import scipy
from  matplotlib.colors import LinearSegmentedColormap
from scipy import stats
from sklearn.cluster import KMeans
from sklearn import metrics
from scipy.spatial.distance import cdist
import multiprocessing as mp
from collections import defaultdict
from sklearn.metrics.pairwise import euclidean_distances
from matplotlib.backends.backend_pdf import PdfPages
from scipy.cluster.hierarchy import fcluster
import scipy.stats as stats
from matplotlib.pyplot import rc
import matplotlib
from scipy.stats import hypergeom
import itertools
import sys
import statsmodels.stats.multitest as multi
from sklearn.decomposition import PCA
from sklearn import preprocessing
from itertools import cycle
import scipy.spatial.distance

from matplotlib import pyplot as plt
from scipy.cluster.hierarchy import dendrogram, linkage
import numpy as np

#%matplotlib inline

os.chdir("${top}")
os.getcwd()

df_original=pd.read_excel("${deg}"+"all_res_counts.xlsx", index_col=0)
print(len(df_original))
df_original.head()

samp_df=pd.read_csv("${samples_treatment}", sep="\t")
samp_df=samp_df.rename(columns={'Unnamed: 0':'Samples'})

samples_dict={}
conditions=list(set([s for s in samp_df['Treat'].tolist()]))
for c in conditions:
    samples=samp_df.loc[samp_df['Treat'] == c , 'Samples'].tolist()
    samples_dict[c]=samples

#conditions_var="${cond_var}"
#conditions=conditions_var.split(" ")
#conditions

df_norm=np.log10(df_original+1)
df_norm.head()

font = {'family' : 'Serif',
            'weight' : 'semibold',
            'size'   : 12}
font_ = {'family' : 'Serif',
            'weight' : 'semibold',
            'size'   : 14}


#### GROUPED KDE ####

fig = plt.figure( frameon=False,figsize=(7,7))

for key in samples_dict.keys():
    dict_value=samples_dict[key]
    values=df_norm[dict_value].values.tolist()
    values=[item for sublist in values for item in sublist]
    values=[s for s in values if s != 0]
    sns.set()
    sns.kdeplot(values, shade=True, label=key)

plt.xlabel("log10(counts+1)",font_)
plt.ylabel("Density",font_)
plt.legend(loc='best', prop=font)
plt.title('Genes', font_)
plt.savefig("${qcp}grouped.KDE.pdf", dpi=600,bbox_inches='tight', pad_inches=0.1,format='pdf')   
plt.show()

#### SAMPLE KDE ####

matplotlib.rc('font', **font)

fig = plt.figure( frameon=False,figsize=(7,7))

for f in df_norm.columns.tolist():
    values=[s for s in df_norm[f].tolist()]
    values=[s for s in values if s != 0]    
    sns.set()
    sns.kdeplot(values, shade=True, label=f)

plt.xlabel("log10(counts+1)",font_)
plt.ylabel("Density",font_)
plt.title('Genes', font_)
plt.legend(bbox_to_anchor=(1, 1), prop=font_)
plt.savefig("${qcp}sample.KDE.pdf", dpi=600,bbox_inches='tight', pad_inches=0.1,format='pdf')   
plt.show()

#### GROUPED BAR PLOTS ####

matplotlib.rc('font', **font)

fig = plt.figure( frameon=False,figsize=(7,7))

medianprops = dict(color='black', linewidth=1)
#flierprops = dict(marker='D', markerfacecolor='black', markersize=6,linestyle='none')

cond_values=[]
for key in samples_dict.keys():
    dict_value=samples_dict[key]
    values=df_norm[dict_value].values.tolist()
    values=[item for sublist in values for item in sublist]
    values=[s for s in values if s != 0]
    cond_values.append(values)
  
sns.set()
bplot=plt.boxplot(cond_values, labels=conditions, patch_artist=True, showfliers=True, medianprops=medianprops)
plt.xticks(fontfamily='Serif')
plt.xlabel("condition",font_)
plt.ylabel("log10(counts+1)",font_)
palette = sns.color_palette("husl", len(conditions))
for patch, color in zip(bplot['boxes'], palette):
        patch.set_facecolor(color)

plt.savefig("${qcp}grouped.barPlots.pdf", dpi=600,bbox_inches='tight', pad_inches=0.1,format='pdf')   
plt.show()

#### SAMPLE BAR PLOTS ####

samples=[col for col in df_norm.columns.tolist()]
values=[]

for col in df_norm.columns.tolist():
    tmp=df_norm[col].values.tolist()
    tmp=[s for s in tmp if s != 0]
    values.append(tmp)

matplotlib.rc('font', **font)

fig = plt.figure(frameon=False,figsize=(10,7))

sns.set()
bplot=plt.boxplot(values, labels=samples, patch_artist=True, showfliers=True, medianprops=medianprops)
plt.xticks(rotation=90, fontfamily='Serif')
plt.xlabel("condition",font_)
plt.ylabel("log10(counts+1)",font_)

palette = sns.color_palette("husl", len(samples))
for patch, color in zip(bplot['boxes'], palette):
        patch.set_facecolor(color)

plt.savefig("${qcp}sample.barPlots.pdf", dpi=600,bbox_inches='tight', pad_inches=0.1,format='pdf')
plt.show()

#### SCATTER PLOT MATRIX ####

file_name="${qcp}count.matrix.scatter.plot.pdf"

pdf = PdfPages(file_name)

w=len(samples_dict.keys())
l=len(samples_dict.keys())
i=1

fig = plt.figure(figsize=(w*12,l*12))

for comp1 in samples_dict.keys():
    value_1=samples_dict[comp1]
    for comp2 in samples_dict.keys():
        value_2=samples_dict[comp2]
        if comp1 == comp2:            
            sns.set()
            ax = fig.add_subplot(l,w,i)
            values=df_norm[value_1].values.tolist()
            values=[item for sublist in values for item in sublist]
            values=[s for s in values if s != 0]
            sns.kdeplot(values)
            #ax.set_xlabel("log10(counts+1)",font_)
            ax.set_ylabel(comp1,font_)
            ax.set_title(comp2, font)
            i=i+1
        else:
            cols_y=value_1
            cols_x=value_2
            tmp = df_norm[value_1 + value_2]
            tmp = tmp[ tmp != 0]
            ax = fig.add_subplot(l,w,i)
            x=tmp[cols_x].mean(axis = 1)
            y=tmp[cols_y].mean(axis = 1)
            ax.scatter( x, y , color="tab:blue", alpha=0.5, s=1,)
            ax.set_xlabel("log10(counts+1)",font_)
            ax.set_ylabel(comp1+"\n\nlog10(counts+1)",font_)
            #ax.legend(loc='best', prop=font_)
            ax.set_title(comp2, font)
            i=i+1
        if i == w*l+1:
            plt.savefig(pdf, dpi=300, bbox_inches='tight', pad_inches=0.1,format='pdf')
	    #plt.show()
            #plt.close()
            i=1
            fig = plt.figure(figsize=(w*12,l*12))
if i != 1:
    plt.savefig(pdf, dpi=300, bbox_inches='tight', pad_inches=0.1,format='pdf')
    #plt.show()
    #plt.close()
    i=1
    fig = plt.figure(figsize=(w*12,l*12))
pdf.close()

#### DENDROGRAM ####

df_norm.T.index

# plt.figure(figsize=(10, 10))
# plt.ylabel('distance')
# dendrogram(linkage(df_norm.T, 'ward'), truncate_mode='lastp', p=len(conditions))
# plt.show()

plt.figure(figsize=(10, 10))
plt.ylabel('distance')
dendrogram(linkage(df_norm.T, 'ward'), labels=df_norm.T.index, leaf_rotation=90.0)
plt.savefig("${qcp}sample.dendrogram.pdf", dpi=600,bbox_inches='tight', pad_inches=0.1,format='pdf')   
plt.show()

path="${deg}"+"annotated/"
data_files=os.listdir(path)
data_files=[s for s in data_files if ".results.tsv" in s]

conditions_new=[s.split("_",1)[1].split(".",1)[0].split("_vs_") for s in data_files]
conditions_new=[item for sublist in conditions_new for item in sublist]
conditions_new=list(set(conditions_new))
conditions_new

#### MA PLOTS ####

file_name="${qcp}MA.plots.pdf"

pdf = PdfPages(file_name)

w=len(conditions_new)
l=len(conditions_new)
i=1

fig = plt.figure(figsize=(w*12,l*12))

for comp1 in conditions_new:
    for comp2 in conditions_new:
        if comp1 == comp2:
            
            sns.set_style("whitegrid")
            ax = fig.add_subplot(l,w,i)
            
            ax.set_xlabel("log10(baseMean)",font_)
            ax.set_ylabel(comp1+"\n\nlog2FoldChange",font_)
            ax.set_title(comp2, font)
            
            i=i+1
            
        else:
            file=[s for s in data_files if comp1 in s and comp2 in s]
            df_tmp=pd.read_csv(path+file[0], sep="\t")
            
            ax = fig.add_subplot(l,w,i)
            
            tmp=df_tmp[df_tmp["padj"]>=0.05]
            x=tmp["baseMean"].tolist()
            y=tmp["log2FoldChange"].tolist()
            x=[ np.log10(s) for s in x ]
            ax.scatter( x, y , color="k", alpha=0.5, s=1, label='NonSignificant')

            tmp=df_tmp[df_tmp["padj"]<0.05]
            x=tmp["baseMean"].tolist()
            y=tmp["log2FoldChange"].tolist()
            x=[ np.log10(s) for s in x ]
            ax.scatter( x, y , color="r", alpha=0.5, s=1, label='Significant')

            ax.set_xlabel("log10(baseMean)",font_)
            ax.set_ylabel(comp1+"\n\nlog2FoldChange",font_)
            ax.legend(loc='best', prop=font_)
            ax.set_title(comp2, font)


            i=i+1
        
        if i == w*l+1:
            plt.savefig(pdf, dpi=300, bbox_inches='tight', pad_inches=0.1,format='pdf')
	    #plt.show()
            #plt.close()
            i=1
            fig = plt.figure(figsize=(w*12,l*12))

if i != 1:
    plt.savefig(pdf, dpi=300, bbox_inches='tight', pad_inches=0.1,format='pdf')
    #plt.show()
    #plt.close()
    i=1
    fig = plt.figure(figsize=(w*12,l*12))

pdf.close()

#### VOLCANO PLOTS ####

file_name="${qcp}volcano.plots.pdf"

pdf = PdfPages(file_name)

w=len(conditions_new)
l=len(conditions_new)
i=1

fig = plt.figure(figsize=(w*12,l*12))

for comp1 in conditions_new:
    for comp2 in conditions_new:
        if comp1 == comp2:
            
            sns.set_style("whitegrid")
            ax = fig.add_subplot(l,w,i)
            
            ax.set_xlabel("log2FoldChange",font_)
            ax.set_ylabel(comp1+"\n\n-log10(Adj.P.value)",font_)
            ax.set_title(comp2, font)
            
            i=i+1
            
        else:
            file=[s for s in data_files if comp1 in s and comp2 in s]
            df_tmp=pd.read_csv(path+file[0], sep="\t")
            
            ax = fig.add_subplot(l,w,i)
            
            tmp=df_tmp[df_tmp["padj"]>=0.05]
            x=tmp["log2FoldChange"].tolist()
            y=tmp["padj"].tolist()
            y=[ -np.log10(s) for s in y ]
            ax.scatter( x, y , color="k", alpha=0.5, s=1, label='NonSignificant')

            tmp=df_tmp[df_tmp["padj"]<0.05]
            x=tmp["log2FoldChange"].tolist()
            y=tmp["padj"].tolist()
            y=[ -np.log10(s) for s in y ]
            ax.scatter( x, y , color="r", alpha=0.5, s=1, label='Significant')

            ax.set_xlabel("log2FoldChange",font_)
            ax.set_ylabel(comp1+"\n\n-log10(Adj.P.value)",font_)
            ax.legend(loc='best', prop=font_)
            ax.set_title(comp2, font)

            i=i+1
        
        if i == w*l+1:
            plt.savefig(pdf, dpi=300, bbox_inches='tight', pad_inches=0.1,format='pdf')
	    #plt.show()
            #plt.close()
            i=1
            fig = plt.figure(figsize=(w*12,l*12))

if i != 1:
    plt.savefig(pdf, dpi=300, bbox_inches='tight', pad_inches=0.1,format='pdf')
    #plt.show()
    #plt.close()
    i=1
    fig = plt.figure(figsize=(w*12,l*12))

pdf.close()

#### P.VALUE DISTRIBUTION ####

file_name="${qcp}p.value.dist.pdf"

pdf = PdfPages(file_name)

w=len(conditions_new)
l=len(conditions_new)
i=1

fig = plt.figure(figsize=(w*12,l*12))

for comp1 in conditions_new:
    for comp2 in conditions_new:
        if comp1 == comp2:
            
            sns.set_style("whitegrid")
            ax = fig.add_subplot(l,w,i)
            
            ax.set_xlabel("P.Value",font_)
            ax.set_ylabel(comp1+"\n\nFrequency",font_)
            ax.set_title(comp2, font)
            
            i=i+1
            
        else:
            file=[s for s in data_files if comp1 in s and comp2 in s]
            df_tmp=pd.read_csv(path+file[0], sep="\t")
            
            ax = fig.add_subplot(l,w,i)
            
            df_tmp['pvalue'].hist(bins=30)
            
            ax.set_xlabel("P.Value",font_)
            ax.set_ylabel(comp1+"\n\nFrequency",font_)
            ax.set_title(comp2, font)

            i=i+1
        
        if i == w*l+1:
            plt.savefig(pdf, dpi=300, bbox_inches='tight', pad_inches=0.1,format='pdf')
	    #plt.show()
            #plt.close()
            i=1
            fig = plt.figure(figsize=(w*12,l*12))

if i != 1:
    plt.savefig(pdf, dpi=300, bbox_inches='tight', pad_inches=0.1,format='pdf')
    #plt.show()
    #plt.close()
    i=1
    fig = plt.figure(figsize=(w*12,l*12))

pdf.close()

#### Q.VALUE DISTRIBUTION ####

file_name="${qcp}q.value.dist.pdf"

pdf = PdfPages(file_name)

w=len(conditions_new)
l=len(conditions_new)
i=1

fig = plt.figure(figsize=(w*12,l*12))

for comp1 in conditions_new:
    for comp2 in conditions_new:
        if comp1 == comp2:
            
            sns.set_style("whitegrid")
            ax = fig.add_subplot(l,w,i)
            
            ax.set_xlabel("Q.Value",font_)
            ax.set_ylabel(comp1+"\n\nFrequency",font_)
            ax.set_title(comp2, font)
            
            i=i+1
            
        else:
            file=[s for s in data_files if comp1 in s and comp2 in s]
            df_tmp=pd.read_csv(path+file[0], sep="\t")
            
            ax = fig.add_subplot(l,w,i)
            
            df_tmp['padj'].hist(bins=30)
            
            ax.set_xlabel("Q.Value",font_)
            ax.set_ylabel(comp1+"\n\nFrequency",font_)
            ax.set_title(comp2, font)

            i=i+1
        
        if i == w*l+1:
            plt.savefig(pdf, dpi=300, bbox_inches='tight', pad_inches=0.1,format='pdf')
	    #plt.show()
            #plt.close()
            i=1
            fig = plt.figure(figsize=(w*12,l*12))

if i != 1:
    plt.savefig(pdf, dpi=300, bbox_inches='tight', pad_inches=0.1,format='pdf')
    #plt.show()
    #plt.close()
    i=1
    fig = plt.figure(figsize=(w*12,l*12))

pdf.close()

#### GENESET LEVEL PLOTS ####

genes=[]
for f in data_files:
    tmp=pd.read_csv(path+f, sep="\t")
    tmp_genes=tmp.loc[tmp['padj'] < 0.05 , 'ensembl_gene_id'].tolist()
    for g in tmp_genes:
        if g not in genes:
            genes.append(g)
            
len(genes)

df_heat=df_norm.loc[df_norm.index.isin(genes),]
df_heat.head()

#print(conditions)

for key in samples_dict.keys():
    value=samples_dict[key]
    df_heat['mean('+key+')']=df_heat[value].mean(axis=1)
df_heat.head()

#### GROUPED HEATMAP ####

matplotlib.rc('font', **font)

#fig = plt.figure( frameon=False,figsize=(14,14))
cm=sns.clustermap(df_heat[[s for s in df_heat.columns.tolist() if 'mean' in s]],figsize=(14,14))
cm.ax_row_dendrogram.set_visible(False)
cm.ax_col_dendrogram.set_visible(False)
plt.savefig("${qcp}grouped.heatMap.pdf", dpi=600,bbox_inches='tight', pad_inches=0.1,format='pdf')   
plt.show()

#### SAMPLE HEATMAP ####

matplotlib.rc('font', **font)

#fig = plt.figure( frameon=False,figsize=(14,14))
cm=sns.clustermap(df_heat[[s for s in df_heat.columns.tolist() if 'mean' not in s]],figsize=(14,16))
cm.ax_row_dendrogram.set_visible(False)
cm.ax_col_dendrogram.set_visible(False)
plt.savefig("${qcp}sample.heatMap.pdf", dpi=600,bbox_inches='tight', pad_inches=0.1,format='pdf')   
plt.show()

#### SIGNIFICANT FEATURES MATRIX ####

df_sig=pd.DataFrame(columns=conditions_new, index=conditions_new)

for comp1 in conditions_new:
    for comp2 in conditions_new:
        if comp1 == comp2:
            sig_fearures=0
            df_sig.loc[comp1,comp2]=sig_fearures
        else:
            file=[s for s in data_files if comp1 in s and comp2 in s]
            tmp=pd.read_csv(path+file[0], sep="\t")
            sig_fearures=len(tmp.loc[tmp['padj'] < 0.05,])
            df_sig.loc[comp1,comp2]=sig_fearures
        
df_sig

mask = np.triu(np.ones_like(df_sig.astype(float), dtype=np.bool))

matplotlib.rc('font', **font)

fig = plt.figure( frameon=False,figsize=(9,7))

sns.heatmap(df_sig.astype(float),annot=True,fmt='.10g', mask=mask)

plt.title('No of significant features', font_)
plt.yticks(rotation=0) 

plt.savefig("${qcp}sigFeatures.matrix.pdf", dpi=600,bbox_inches='tight', pad_inches=0.1,format='pdf')   
plt.show()

df_corr=df_norm.copy()
for key in samples_dict.keys():
    value=samples_dict[key]
    df_corr['mean('+key+')']=df_corr[value].mean(axis=1)

df_corr.head()

#### GROUP DISTANCE MATRIX ####

dist=metrics.pairwise_distances(np.array(df_corr[[s for s in df_corr.columns.tolist() if 'mean' in s]].T))

df_dist=pd.DataFrame(dist, columns=[s.split("(")[1].strip(")") for s in df_corr.columns.tolist() if 'mean' in s], index=[s.split("(")[1].strip(")") for s in df_corr.columns.tolist() if 'mean' in s])
matplotlib.rc('font', **font)

annot = {'family' : 'Serif',
         'weight' : 'semibold',
         'size'   : 12}

fig = plt.figure( frameon=False,figsize=(12,10))

sns.heatmap(df_dist.T,annot=True,fmt='.5g',annot_kws=annot)

plt.title('Group Distance Matrix', font_)

plt.savefig("${qcp}group.distance.matrix.pdf", dpi=600,bbox_inches='tight', pad_inches=0.1,format='pdf')   
plt.show()

#### SAMPLE DISTANCE MATRIX ####

dist_all=metrics.pairwise_distances(np.array(df_corr[[s for s in df_corr.columns.tolist() if 'mean' not in s]].T))

df_dist_all=pd.DataFrame(dist_all, columns=[s for s in df_corr.columns.tolist() if 'mean' not in s], index=[s for s in df_corr.columns.tolist() if 'mean' not in s])
matplotlib.rc('font', **font)

annot = {'family' : 'Serif',
         'weight' : 'semibold',
         'size'   : 10}

fig = plt.figure( frameon=False,figsize=(12,10))

sns.heatmap(df_dist_all.T,annot=True,fmt='.3g',annot_kws=annot)

plt.title('Sample Distance Matrix', font_)

plt.savefig("${qcp}sample.distance.matrix.pdf", dpi=600,bbox_inches='tight', pad_inches=0.1,format='pdf')
plt.show()

#### PCA ####

file_name="${qcp}pca.pdf"

pdf = PdfPages(file_name)

w=len(conditions_new)
l=len(conditions_new)
i=1

fig = plt.figure(figsize=(w*12,l*12))

for comp1 in conditions_new:
    for comp2 in conditions_new:
        if comp1 == comp2:
            
            sns.set_style("darkgrid")
            ax = fig.add_subplot(l,w,i)
            
            ax.set_xlabel("Component1",font_)
            ax.set_ylabel(comp1+"\n\nComponent2",font_)
            ax.set_title(comp2, font)
            
            i=i+1
            
        else:
            file=[s for s in data_files if comp1 in s and comp2 in s]
            df_tmp=pd.read_csv(path+file[0], sep="\t")
            
            ax = fig.add_subplot(l,w,i)
            
            data=df_tmp.copy()
            
            cols=[s for s in data.columns if "_R" in s]
            data_pca=data[cols]
            data_pca=np.log10(data_pca + 1)
            data_pca=data_pca.set_index(data['ensembl_gene_id'])
            #data_pca.head()

            df_pca=data_pca.T.reset_index()
            #df_pca

            #print(df_pca.shape)

            df_pca.set_index('index', inplace=True)
            df_pca.index.names = ['Sample']

            pca = PCA(copy=True, iterated_power='auto', n_components=2, random_state=None,svd_solver='auto', tol=0.0, whiten=False)

            #scaling the values
            df_pca_scaled = preprocessing.scale(df_pca, axis = 1)

            projected=pca.fit_transform(df_pca_scaled)

            #print(pca.explained_variance_ratio_)

            tmp=pd.DataFrame(projected)
            tmp.rename(columns={0: 'Component 1',1: 'Component 2'}, inplace=True)
            #tmp.head()

            final_pca=pd.merge(df_pca.reset_index(),tmp,left_index=True,right_index=True)
            #final_pca.head()

            for s in final_pca['Sample'].tolist():
                for key in samples_dict.keys():
                    value=samples_dict[key]
                    if s in value:
                        final_pca.loc[final_pca['Sample'] == s,'expt'] = key

            font = {'family' : 'serif',
                    'weight' : 'bold',
                    'size'   : 12}


            color_gen = cycle(('blue', 'red', 'green', 'pink','yellow'))
            
            for lab in set(final_pca['expt']):
                plt.scatter(final_pca.loc[final_pca['expt'] == lab, 'Component 1'], 
                            final_pca.loc[final_pca['expt'] == lab, 'Component 2'], 
                            c=next(color_gen),
                            label=lab)

            ax.set_xlabel('Component 1  - ' + str(pca.explained_variance_ratio_[0]*100)+ " % ", fontdict=font)
            ax.set_ylabel(comp1+"\n\nComponent 2  - "+ str(pca.explained_variance_ratio_[1]*100)+ "  % ", fontdict=font)
            #ax.set_xticks(fontsize=14)
            #ax.set_yticks(fontsize=14)
            ax.legend(loc='best', fontsize=14, prop=font_)
            ax.set_title(comp2, fontdict=font)

            i=i+1
        
        if i == w*l+1:
            plt.savefig(pdf, dpi=300, bbox_inches='tight', pad_inches=0.1,format='pdf')
	    #plt.show()
            #plt.close()
            i=1
            fig = plt.figure(figsize=(w*12,l*12))

if i != 1:
    plt.savefig(pdf, dpi=300, bbox_inches='tight', pad_inches=0.1,format='pdf')
    #plt.show()
    #plt.close()
    i=1
    fig = plt.figure(figsize=(w*12,l*12))

pdf.close()



### PCA all ####

file_name="${qcp}pca_all_samples.pdf"

pdf = PdfPages(file_name)

plt.figure(figsize=(14,14))

#sns.set_style("darkgrid")

pca_data = df_heat[[s for s in df_heat.columns.tolist() if 'mean' not in s]]
df_pca=pca_data.T.reset_index()

df_pca.set_index('index', inplace=True)
df_pca.index.names = ['Sample']

pca = PCA(copy=True, iterated_power='auto', n_components=2, random_state=None,svd_solver='auto', tol=0.0, whiten=False)

#scaling the values
df_pca_scaled = preprocessing.scale(df_pca, axis = 1)

projected=pca.fit_transform(df_pca_scaled)
#print(pca.explained_variance_ratio_)

tmp=pd.DataFrame(projected)
tmp.rename(columns={0: 'Component 1',1: 'Component 2'}, inplace=True)
#tmp.head()

final_pca=pd.merge(df_pca.reset_index(),tmp,left_index=True,right_index=True)
#final_pca.head()

for s in final_pca['Sample'].tolist():
    for key in samples_dict.keys():
        value=samples_dict[key]
        if s in value:
            final_pca.loc[final_pca['Sample'] == s,'expt'] = key

font = {'family' : 'serif',
        'weight' : 'bold',
        'size'   : 12}

color_gen = cycle(('blue', 'red', 'green', 'pink','yellow'))
color_gen = sns.color_palette(None, len(conditions))

plot_data = final_pca[["Sample", "Component 1", "Component 2", "expt"]]

# map group color
plot_data['color'] = plot_data['expt'].map(dict(zip(set(plot_data["expt"]), color_gen)))

for lab in set(plot_data['expt']):
    scatter = plt.scatter(plot_data.loc[plot_data['expt'] == lab, 'Component 1'],
    plot_data.loc[plot_data['expt'] == lab, 'Component 2'],
    c=plot_data.loc[plot_data['expt'] == lab, 'color'],
    label=lab)

plt.legend(handles=scatter.legend_elements()[0], labels=set(plot_data["expt"]), title = 'Group')

plt.xlabel('Component 1  - ' + str(pca.explained_variance_ratio_[0]*100)+ " % ", fontdict=font)
plt.ylabel('Component 2  - ' + str(pca.explained_variance_ratio_[1]*100)+ "  % ", fontdict=font)
plt.title("PCA of all Samples", font)

plt.savefig(pdf, dpi=300, bbox_inches='tight', pad_inches=0.1,format='pdf')

pdf.close()
 

PYOF

id=$(sbatch --partition $SLURMPARTITION --parsable << EOF
#!/bin/bash
#SBATCH --output ${logs}qc_plots.%j.out
#SBATCH -c 4
#SBATCH --job-name="QCplots"

${SING} << SIN
#!/bin/bash
${HOMESOURCE}

module purge
module load python/3.8.0

python3 ${scripts}QC_plots.py

SIN
EOF
)

echo "Waiting for QC plots job ${id} to complete."
srun --partition $SLURMPARTITION -d afterok:${id} echo "Finished pipeline. Done."

fi


exit

