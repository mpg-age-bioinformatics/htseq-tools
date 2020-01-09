#!/bin/bash

# TODO: Set General pipeline settings
# MPI-Age Settings
HOMESOURCE="source ~/.bashrc"
SLURMPARTITION="blade,himem,hugemem"
SHIFTER="shifter --image=hub.age.mpg.de/bioinformatics/software:v2.0.6 /bin/bash"
SHIFTER_FLEXBAR="shifter --image=hub.age.mpg.de/flexbar:3.5.0--hb7ba0dd_3 /bin/bash"

# reference genome
genome_folder=/beegfs/common/genomes/mus_musculus/91/
gtf=${genome_folder}original.gtf
chromosomes=${genome_folder}original.toplevel.genome
fasta=${genome_folder}original.primary_assembly.fa

# barcode fasta files
barcodes=$(readlink -f .)/barcodes.fa

# raw files
umi_file="/beegfs/group_bit/data/raw_data/departments/Martin_Denzel/MD_Kira_a6CD34/mpg_L18583-1_1-16seq_S126_R1_001.fastq.gz"
rna_file="/beegfs/group_bit/data/raw_data/departments/Martin_Denzel/MD_Kira_a6CD34/mpg_L18583-1_1-16seq_S126_R2_001.fastq.gz"

# demultiplexed files sufixes
read1_sufix="-READ_2.fastq.gz"
read2_sufix="none" # if single, write "none"

# "caenorhabditis elegans", "drosophila melanogaster", "mus musculus", "homo sapiens", "saccharomyces cerevisiae"
species="mus musculus"

# http://www.ensembl.org/info/website/archives/index.html
biomart_host="http://jan2019.archive.ensembl.org/biomart/"

# datasets "celegans_gene_ensembl" "dmelanogaster_gene_ensembl","mmusculus_gene_ensembl","hsapiens_gene_ensembl", "scerevisiae_gene_ensembl"
biomart_dataset="mmusculus_gene_ensembl"

# for single end reads in kalisto
fragment_size=75
fragment_size_standard_deviation=10

# DAVID
daviddatabase="ENSEMBL_GENE_ID"
DAVIDUSER="franziska.metge@age.mpg.de"

# cytoscape
cytoscape_host="192.168.50.166"

# samples, if paired only READ is required
samplestable="/beegfs/group_bit/data/projects/departments/Martin_Denzel/MD_Kira_a6CD34/scripts.FMetge/RNAseq_sample_sheet.xlsx"
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

#############################################################################

echo "Creating required folders"
mkdir -p ../slurm_logs
mkdir -p ../fastqc_output
mkdir -p ../tmp
mkdir -p ../kallisto_output
mkdir -p ../deseq2_output
mkdir -p ../featureCounts_output_kal
mkdir -p ../multiqc_output
mkdir -p ../kallisto_index

top=$(readlink -f ../)/
tmp=$(readlink -f ../tmp)/
raw=$(readlink -f ../raw_data)/
logs=$(readlink -f ../slurm_logs)/
fqc=$(readlink -f ../fastqc_output)/
kalout=$(readlink -f ../kallisto_output)/
deg=$(readlink -f ../deseq2_output)/
mqc=$(readlink -f ../multiqc_output)/
ind=$(readlink -f ../kallisto_index)/
scripts=$(readlink -f .)/
cdna_fasta=${ind}transcripts.norRNA.fa
kallisto_index=${ind}transcripts.norRNA.idx



#############################################################################

echo "Demultiplexing reads."

#############################################################################

if [ ! -e ${raw}/demultiplexed.txt ] ; then 

rm -rf ${logs}demultiplex.*.out

id=$(sbatch --partition $SLURMPARTITION --parsable << EOF
#!/bin/bash
#SBATCH --output ${logs}demultiplex.%j.out
#SBATCH -c 4
#SBATCH --job-name="demult"
module load shifter

${SHIFTER_FLEXBAR} << SHI
#!/bin/bash
${HOMESOURCE}

which flexbar

cd ${raw}

flexbar -r ${umi_file} -p ${rna_file} -b $barcodes --umi-tags -n 4 -t ${raw} -u 5 -q TAIL -qt 25 -qf i1.8 -z GZ


SHI
EOF
)

srun -p $SLURMPARTITION -d afterok:${id} echo "Done demultiplexing."


#############################################################################

echo "rename files."

#############################################################################

cd $raw

for f in $(ls *.fastq.gz);

do

        mv $f ${f#_barcode_}

done


#############################################################################

echo "reformat read name."

#############################################################################

cd ${raw}

for f in $(ls *.fastq.gz);

do
echo $f

zcat $f | sed 's/ /:/' > ${f%.gz}
gzip -f ${f%.gz}

done
fi

touch demultiplexed.txt


#############################################################################

#echo "Prepare gene id to gene name reference table"
echo "Extracting non rRNA transcripts"

############################################################################# 


if [ ! -e ${cdna_fasta} ] ; then 

grep -v -i rrna ${gtf} > ${tmp}no.rRNA.gtf
#echo "Done with grep"
#module load gff
# needs samtools faidx ${fasta}
# https://github.com/gpertea/gffread
~/gffread/gffread -w ${cdna_fasta} -g ${fasta} ${tmp}no.rRNA.gtf

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

if [ ! -e ../fastqc_output/${file%.fastq.gz}_fastqc.html ] ; then

rm -rf ${logs}fastqc.${file%.fastq.gz}.*.*

sbatch --partition $SLURMPARTITION << EOF
#!/bin/bash
#SBATCH --output ${logs}fastqc.${file%.fastq.gz}.%j.out
#SBATCH --error ${logs}fastqc.${file%.fastq.gz}.%j.err
#SBATCH -c 4
#SBATCH --job-name='fastqc'
module load shifter

${SHIFTER} << SHI
#!/bin/bash
${HOMESOURCE}
module load jdk fastqc
cd ${raw}
# FASTQC call
fastqc -t 4 -o ../fastqc_output ${file}
SHI
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
module load shifter

${SHIFTER} << SHI
#!/bin/bash
${HOMESOURCE}

module load kallisto 

cd ${ind}

kallisto index -i ${kallisto_index} ${cdna_fasta}

SHI
EOF
)

srun -d afterok:${id} echo "Done building Kallisto's index."

fi

#############################################################################

echo "Determining strandness"

#############################################################################

if [ ! -e ${raw}strandness.txt ] ; then

rm -rf ${logs}strand.${test_read_1%${read1_sufix}}.*.out

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
#SBATCH --output ${logs}strand.${test_read_1%${read1_sufix}}.%j.out
#SBATCH -c 18
#SBATCH --job-name="strand"
module load shifter

${SHIFTER} << SHI
#!/bin/bash
${HOMESOURCE}

module load kallisto

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

module load python
cd ../
unset PYTHONHOME
pip install virtualenv --user
rm -rf RSeQC  
virtualenv RSeQC
unset PYTHONUSERBASE
source RSeQC/bin/activate
pip install numpy
pip install RSeQC
infer_experiment.py -i raw_data/tmp.${test_read_1%${read1_sufix}}/pseudoalignments.bam -r ${ind}gene.model.bed > raw_data/infer_experiment.txt

python > raw_data/strandness.txt << PYOF
import sys
filein="raw_data/infer_experiment.txt"
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
SHI
EOF
)

echo "Waiting for strandness determination job ${id} to complete."
srun --partition $SLURMPARTITION -d afterok:${id} echo "Strandness determination complete."

rm -rf ${raw}tmp.${test_read_1%${read1_sufix}}  

fi


#############################################################################

strand=$(head -n 1 ${raw}strandness.txt) 
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
module load shifter

${SHIFTER} << SHI
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
SHI
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
srun --partition $SLURMPARTITION -d afterok${ids} echo "Finished Kallisto. Starting to extract reads using bedtools."


#############################################################################

# extracting transcriptomic reads
echo "extracting transcriptomic reads"

############################################################################# 


ids=""
cd ${top}/kallisto_output
echo $(pwd)

umi_tag='$5'
gene_name='$17'
for file in $(ls) ; do

rm -rf ${logs}bedtools_kal_${file}.*.out
ids=${ids}:$(sbatch --partition $SLURMPARTITION --parsable << EOF
#!/bin/bash
#SBATCH --cpus-per-task=18
#SBATCH -o ${logs}bedtools_kal_${file}.%j.out
#SBATCH --job-name="bedtools"
module load shifter

module load bedtools
cd ${top}kallisto_output/${file}

bedtools intersect -a pseudoalignments.bam -b ${ind}gene.model.bed -bed -wo > tmp.bed
sed 's/_/\t/' tmp.bed |awk '{print $umi_tag "\t" $gene_name}' |sort -u > gene_reads.bed

EOF
)

done

ids=${ids}:$(sbatch --partition $SLURMPARTITION --parsable << EOF
#!/bin/bash
#SBATCH --output ${logs}empty.%j.out
EOF
)


echo "Waiting for bedtools jobs ${ids} to complete."
srun --partition $SLURMPARTITION -d afterok${ids} echo "Finished Bedtools. Starting featureCounts and multiqc."


#############################################################################

# featureCounts

#############################################################################

if [ "${strand}" == "fr-firststrand" ] ; then
    featureCounts_direction=2
elif [ "${strand}" == "fr-secondstrand" ] ; then
    featureCounts_direction=1
elif [ "${strand}" == "unstranded" ] ; then
    featureCounts_direction=0
else:
    exit
fi

ids=""
cd ${top}kallisto_output
echo $(pwd)
for file in $(ls) ; do 

if [ ! -e ${top}featureCounts_output_kal/${file}_biotype_counts_mqc.txt ] ; then

rm -rf ${logs}featureCount_kal_${file}.*.out
ids=${ids}:$(sbatch --partition $SLURMPARTITION --parsable << EOF
#!/bin/bash
#SBATCH --cpus-per-task=18
#SBATCH -o ${logs}featureCount__kal_${file}.%j.out
#SBATCH --job-name="featureCounts"
module load shifter

${SHIFTER} << SHI
#!/bin/bash
${HOMESOURCE}
module load subread
cd ${top}kallisto_output/${file}
echo "featureCounts -a $gtf -T 18 -g gene_id -o ${top}featureCounts_output_kal/${file}_gene.featureCounts.txt -p -s $featureCounts_direction pseudoalignments.bam"
featureCounts -a $gtf -T 18 -g gene_id -o ${top}featureCounts_output_kal/${file}_gene.featureCounts.txt -p -s $featureCounts_direction pseudoalignments.bam
echo "featureCounts -a $gtf -T 18 -g gene_biotype -o ${top}featureCounts_output_kal/${file}_biotype.featureCounts.txt -p -s $featureCounts_direction pseudoalignments.bam"
featureCounts -a $gtf -T 18 -g gene_biotype -o ${top}featureCounts_output_kal/${file}_biotype.featureCounts.txt -p -s $featureCounts_direction pseudoalignments.bam
cut -f 1,7 ${top}featureCounts_output_kal/${file}_biotype.featureCounts.txt | tail -n +3 | (echo "$biotypes_header" && cat) >> ${top}featureCounts_output_kal/${file}_biotype_counts_mqc.txt

SHI
EOF
)

fi

done

ids=${ids}:$(sbatch --partition $SLURMPARTITION --parsable << EOF
#!/bin/bash
#SBATCH --output ${logs}empty.%j.out
EOF
)

#############################################################################

# MultiQC

#############################################################################

if [ ! -e ${top}multiqc_kallisto_output/multiqc_report.html ] ; then

rm -rf ${logs}multiqc.*.out
sbatch --partition $SLURMPARTITION -d afterok${ids} --parsable << EOF
#!/bin/bash
#SBATCH --cpus-per-task=1
#SBATCH -o ${logs}multiqc.%j.out
#SBATCH --job-name="multiqc"

${SHIFTER} << SHI
#!/bin/bash
${HOMESOURCE}

# Install multiqc
module load python
#pip install multiqc --user --ignore-installed

cd ${top}
pip install virtualenv --user
unset PYTHONHOME
virtualenv multiqc_kallisto
unset PYTHONUSERBASE
source multiqc_kallisto/bin/activate
pip install multiqc --ignore-installed  
multiqc fastqc_output/ featureCounts_output_kal/ kallisto_output/ -f -o multiqc_kallisto_output

SHI
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

GTF=age.readGTF("${gtf}")
GTF["gene_id"]=age.retrieve_GTF_field(field="gene_id",gtf=GTF)
GTF["transcript_id"]=age.retrieve_GTF_field(field="transcript_id",gtf=GTF)
GTF["gene_biotype"]=age.retrieve_GTF_field(field="gene_biotype",gtf=GTF)
tx2gene=GTF[["transcript_id","gene_id"]].drop_duplicates()
tx2gene.columns=["TXNAME","GENEID"]
tx2gene[["TXNAME","GENEID"]].to_csv("${deg}"+"tx2gene.csv", quoting=1, index=None)

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
module load shifter

${SHIFTER} << SHI
#!/bin/bash
${HOMESOURCE}
module purge

module load python/3.6.5
pip3 install pandas --user
pip3 install xlrd --user
pip3 install biomart --user
pip3 install git+https://github.com/mpg-age-bioinformatics/AGEpy.git --user

python3 ${deg}read.template.py

module load rlang
Rscript ${deg}deseq2.R

SHI
EOF
)


echo "Waiting for DESeq2 part 1 job ${id} to complete."
srun --partition $SLURMPARTITION -d afterok:${id} echo "Finished DESeq2 part 1. Starting DESeq2 pair wise comparisons."

fi

#############################################################################

# DESeq2 part 2 (pair-wise comparisons)
# based on UMI counts per gene

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
files<-file.path(dir, samples, "gene_reads.bed")
names(files)<-samples

# while reading in, count umi per gene, merge to bigger table
gene_counts = read.delim(files[1], header = F, as.is = T)
gene_counts = as.data.frame(table(gene_counts[,2]))
names(gene_counts)[2] = names(files)[1]

for(f in 2:length(files)){
  tmp = read.delim(files[f], header = F, as.is = T)
  tmp = as.data.frame(table(tmp[,2]))
  gene_counts = merge(gene_counts, tmp, by = 'Var1', all = TRUE)
  names(gene_counts)[ncol(gene_counts)] = names(files)[f]
  write.table(gene_counts, "${out}.umi_gene_counts.tsv", sep = "\t", quote = F, row.names = F)
}

gene_counts[is.na(gene_counts)] <- 0
row.names(gene_counts) = gene_counts[,1]
gene_counts = gene_counts[,-1]

# model
dds <- DESeqDataSetFromMatrix(gene_counts, sampleTable, ~${model} )
# model
dds\$${model} <- relevel(dds\$${model}, ref = "${ref}")

keep <- rowSums(counts(dds) >= 10) == nrow(sampleTable)
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
module load shifter

${SHIFTER} << SHI
#!/bin/bash
${HOMESOURCE}
module purge
module load rlang
which R

Rscript ${out}.deseq2.R

SHI
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
srun --partition $SLURMPARTITION -d afterok${ids} echo "Finished DESeq2 part 2. Starting annotation."

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
PYOF

id=$(sbatch --partition $SLURMPARTITION --parsable << EOF
#!/bin/bash
#SBATCH --output ${logs}annotate.deseq2.%j.out
#SBATCH -c 4
#SBATCH --job-name="adeseq2"
module load shifter

${SHIFTER} << SHI
#!/bin/bash
${HOMESOURCE}
module purge

module load python/3.6.5

python3 ${deg}annotate.deseq2.py

SHI
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

${SHIFTER} << SHI
#!/bin/bash
${HOMESOURCE}
module purge

module load python/3.6.5

python3 ${deg}${PY}

SHI
EOF
)

fi

done

echo "Waiting for DAVID annotation  job ${ids} to complete."
srun --partition $SLURMPARTITION -d afterok:${ids} echo "Finished DAVID annotation. Starting cytoscape and drawing DAIVD plots."


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
module load shifter

cd ${deg}annotated
for f in *.results.tsv ; do

${SHIFTER} << SHI
#!/bin/bash
${HOMESOURCE}
module purge

module load python/3.6.5

python3 ${deg}PPIs.py ${deg}annotated/\${f}

SHI
done
EOF

#############################################################################

# Draw DAVID plots

#############################################################################

cd ${scripts}

wget https://raw.githubusercontent.com/mpg-age-bioinformatics/htseq-tools/master/david_to_cellplot.R
chmod 775 david_to_cellplot.R


sbatch --partition $SLURMPARTITION --parsable << EOF
#!/bin/bash
#SBATCH --output ${logs}david_plot%j.out
#SBATCH -c 2
#SBATCH --job-name="david_plot"
module load shifter


for f in ${deg}/annotated/*DAVID.tsv ; do

${SHIFTER} << SHI
#!/bin/bash
${HOMESOURCE}

module load rlang

Rscript david_to_cellplot.R ${f} tsv GOTERM_BP_FAT 15
Rscript david_to_cellplot.R ${f} tsv KEGG_PATHWAY 15


SHI
done
EOF

exit




