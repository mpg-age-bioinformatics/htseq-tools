#!/bin/bash
#source /beegfs/common/software/2017/age-bioinformatics.2017.only.rc

# this script starts the variant calling pipeline.
# usage: 
# ./complete.pipe.S_xxx.sh S_012 blade,himem mm9 150
# in order to start the pipeline for sample S_012, make sure sample S_012 exists as 1.fastq.gz and 2.fastq.gz in the raw_data folder
# the second parameter specifies the partition you want to submit your job to
# the thrid parameter specifies which organims the sample is coming from. mm9 and hg38 are supported right now
# note that there is neither a wt nor a synthetic wt for hg38, outputfiles are not filtered for any mutations, needs to be adjusted once there is enough data
# when uploading data using bit always upload the scripts using --scripts in order to track with which paramters the programs were run

sample=$1
partition=$2
genome=$3 # either mm9 or hg38
mem=$4 # specify the number of GB memory
top=$(readlink -f ../)/
raw=$(readlink -f ../raw_data)/
fastqc_out=$(readlink -f ../fastqc_output)/
bwa=$(readlink -f ../bwa_output)/
ss=$(readlink -f ../sorted_sam)/
md=$(readlink -f ../mark_duplicates)/
rg=$(readlink -f ../replace_groups)/
rs=$(readlink -f ../reorder_sam)/
others=$(readlink -f ../others)/
rec=$(readlink -f ../recalibrated_output)/
pu=$(readlink -f ../pileup_output)/
tmp=$(readlink -f ../tmp)/
logs=$(readlink -f ../slurm_logs)/
temp=${tmp}/
real=$(readlink -f ../realign_output)/
ssb=$(readlink -f ../sorted_sam_b)/
#echo $synthetic_background
mkdir -p ${fastqc_out}  ${bwa} ${ss} ${md} ${rg} ${rs} ${others} ${rec} ${refs} ${tmp} ${logs} ${pu} ${real} ${ssb} ${logs}

echo ${top}
echo "starting pipeline on sample ${sample}" > ${logs}/${sample}.pipeline.log
date >> ${logs}/${sample}.pipeline.log

# set reference genome
if [ "${genome}" == 'hg38' ] ;
    then
        kegg_genome="hsa" # test
        refs=$(readlink -f ../references)/hg38
        cd ${refs}
        bwa_index=$(readlink -f  bwa_index/hg38.fa)
        vcf=$(readlink -f  common_all_20170710.fixed.sorted.vcf) 
        faref=$(readlink -f hg38.fa)
        exomebed=$(readlink -f S04380110_Covered.hg38.bed) # get correct one from moritz
        genoem_flag="H1"
        synthetic_background=${top}/synthetic_wildtype/S_043_human_synthetic_wt.vcf
        #wt_file=${top}/filter/S_044-F_HORN-L____H1-___-____-REP_1.SNPs.vcf
        wt_file=${top}/filter/S_054-F_HORN-L____H1-___-____-REP_1.SNPs.vcf
        echo "yes"
elif [ "${genome}" == 'mm9' ] ;
    then
        kegg_genome="mmu"
        refs=$(readlink -f ../references)/mm9
        cd ${refs}
        bwa_index=$(readlink -f bwa_index/mm9.fa)
        vcf=$(readlink -f C57BL_6NJ.mgp.v5.snps.dbSNP142.ucsc.mm9.s.vcf)
        faref=$(readlink -f mm9.fa)
        exomebed=$(readlink -f S0276129_Covered_format.bed)
        genoem_flag="M1"
        synthetic_background=${top}/synthetic_wildtype/variantes_synthetic_S_013-S_035.txt
        wt_file=${top}/filter/S_016-F_HORN-L____M1-___-____-REP_1.SNPs.vcf
        lucky_file=${top}/filter/S_002-F_HORN-L____M1-___-____-REP_1.SNPs.vcf
        echo "yes"
else
    echo "ERROR: You did not specify an available genome. Available genomes are hg38 and mm9"
    exit 1
fi


cd ${raw}

IDS=

for f in $( ls ${sample}*1.fastq.gz);
	do echo "#!/bin/bash
	module load bwa/0.7.15
	module load samtools/1.3.1
	module load picard/2.8.1
	module load gatk/3.4.46
	
	fastqc -t 10 -o ${fastqc_out} ${sample}*.fastq.gz
	echo 'BWA'
	bwa mem -t 18 -M ${bwa_index} ${f} ${f%1.fastq.gz}2.fastq.gz > ${bwa}${f%-READ_1.fastq.gz}.sam 
	echo 'sam to bam'
	samtools view -bS ${bwa}${f%-READ_1.fastq.gz}.sam > ${bwa}${f%-READ_1.fastq.gz}.bam
	
	samtools sort -@ 10 -o ${bwa}${f%-READ_1.fastq.gz}.sorted.bam ${bwa}${f%-READ_1.fastq.gz}.bam
	
	samtools index ${bwa}${f%-READ_1.fastq.gz}.sorted.bam

	samtools flagstat ${bwa}${f%-READ_1.fastq.gz}.bam > ${bwa}${f%-READ_1.fastq.gz}.bam.stat 
	echo 'sort SAM picard'
	java -jar \${PICARD} SortSam SORT_ORDER=coordinate INPUT=${bwa}${f%-READ_1.fastq.gz}.bam \
OUTPUT=${ss}${f%-READ_1.fastq.gz}.bam TMP_DIR=${temp}
	echo 'MarkDuplicates'
	java -Xmx16g -jar \${PICARD} MarkDuplicates I=${ss}${f%-READ_1.fastq.gz}.bam \
O=${md}${f%-READ_1.fastq.gz}.bam CREATE_INDEX=true VALIDATION_STRINGENCY=LENIENT \
TMP_DIR=${temp} M=${md}metrics.txt

	rm -rf ${rg}${f%-READ_1.fastq.gz}.bam; java -jar \${PICARD} AddOrReplaceReadGroups I=${ss}${f%-READ_1.fastq.gz}.bam O=${rg}${f%-READ_1.fastq.gz}.bam SO=coordinate \
RGID=${f%-READ_1.fastq.gz} RGLB=${f%-READ_1.fastq.gz} RGPL=illumina RGPU=${f%-READ_1.fastq.gz} \
RGSM=${f%-READ_1.fastq.gz} CREATE_INDEX=true  

	java -jar \${PICARD} ReorderSam INPUT=${rg}${f%-READ_1.fastq.gz}.bam \
OUTPUT=${rs}${f%-READ_1.fastq.gz}.bam REFERENCE=${faref}
	
	java -jar \${PICARD} BuildBamIndex INPUT=${rs}${f%-READ_1.fastq.gz}.bam
	echo 'IndelRealigner'
	java -jar \${GATK} -nt 18 -T RealignerTargetCreator -R ${faref} -I ${rs}${f%-READ_1.fastq.gz}.bam -o ${others}${f%-READ_1.fastq.gz}.list 

	java -jar \${GATK} -T IndelRealigner -R ${faref} -I ${rs}${f%-READ_1.fastq.gz}.bam \
-targetIntervals ${others}${f%-READ_1.fastq.gz}.list -o ${real}${f%-READ_1.fastq.gz}.bam

	java -jar \${PICARD} BuildBamIndex INPUT=${real}${f%-READ_1.fastq.gz}.bam
	echo 'BaseRecalibrator'
	java -Xmx36g -Djava.io.tmpdir=${temp}${f%-READ_1.fastq.gz} -jar \${GATK} -nct 18 \
-T BaseRecalibrator -knownSites ${vcf} -lowMemoryMode -R ${faref} \
-I ${real}${f%-READ_1.fastq.gz}.bam -o ${others}${f%-READ_1.fastq.gz}.table -L ${exomebed}

	java -jar \${GATK} -T PrintReads -R ${faref} -I ${real}${f%-READ_1.fastq.gz}.bam \
-BQSR ${others}${f%-READ_1.fastq.gz}.table -o ${rec}${f%-READ_1.fastq.gz}.bam \
--validation_strictness LENIENT
	echo 'samtools mpileup'
	rm -rf ${pu}${f%-READ_1.fastq.gz}.vcf; samtools mpileup -l ${exomebed} -f ${faref} -v -t AD,ADR,INFO/AD -u ${rec}${f%-READ_1.fastq.gz}.bam \
-o ${pu}${f%-READ_1.fastq.gz}.vcf

	" > ${tmp}${f%-READ_1.fastq.gz}.sh 
	chmod 755 ${tmp}${f%-READ_1.fastq.gz}.sh
	rm -rf ${logs}${f%-READ_1.fastq.gz}.*.out
	ID=$(sbatch --parsable -o ${logs}${f%-READ_1.fastq.gz}.%j.out -c 18 --mem=${mem}GB -p ${partition} ${tmp}${f%-READ_1.fastq.gz}.sh)
	IDS=:${ID}
done

echo "Done with submissions"
srun --mem=100 -p ${partition} -d afterok${IDS} echo "Waiting for pileup to finish"

filter=$(readlink -f ../filter)/
pileup=$(readlink -f ../pileup_output)/

mkdir -p ${filter}

cd ${pileup}
for f in $(ls *.vcf | grep ${sample}); do egrep -v '(#|INDEL)' ${f} | awk '{if ($5 !="<*>") {print $0}}' > ${filter}${f%vcf}SNPs.vcf; done


if [ "${genome}" == 'mm9' ] || [ "${genome}" == 'hg38' ] 
    then
cd ${filter}
echo "selecting.."

module load python
python - << EOF

import pandas as pd
import numpy as np


def Calcs(x):
    l=x.split("AD=")[1].split(";")[0].split(",")
    l=[ int(s) for s in l ]
    #l=l[:-1]
    sums=np.sum(l)
    smaller=l[1]
    res=str(sums)+","+str(smaller)
    return res

def FILTER(D, min_freq):
    D["tmp"]=D[7].apply(lambda x: Calcs(x) )
    D["sums"]=D["tmp"].apply(lambda x: x.split(",")[0])
    D["smaller alternative al."]=D["tmp"].apply(lambda x: x.split(",")[1] )
    D["freq. (smaller alternative al.) (%)"]=D["smaller alternative al."].astype(float)/D["sums"].astype(float)*100.0
    D=D[[0,1, "sums", "freq. (smaller alternative al.) (%)"]]
    # filter by allele frequency
    D=D[D["freq. (smaller alternative al.) (%)"] > min_freq]
    return D

filterfolder="${filter}"
backgroundfile="${synthetic_background}"
wtfile="${wt_file}"
lucky="${lucky_file}"
selected="${sample}-F_ACUS-L____${genoem_flag}-___-____-REP_1.SNPs.vcf"


outNoBack="${sample}-F_ACUS-L____${genoem_flag}-___-____-REP_1.SNPs.noback.vcf"
outWT="${sample}-F_ACUS-L____${genoem_flag}-___-____-REP_1.SNPs.nowt.vcf"
outLucky="${sample}-F_ACUS-L____${genoem_flag}-___-____-REP_1.SNPs.lucky.vcf"

background=pd.read_table(backgroundfile,header=0, usecols=[0, 1])
wt=pd.read_table(wtfile,header=None, usecols=[0, 1, 7])
lucky_wt = pd.read_table(lucky, header=None, usecols=[0,1])
df_selected=pd.read_table(filterfolder+selected, header=None)

wt = FILTER(wt, 3)

background=list(background["Chr"].astype(str)+"_"+background["Pos"].astype(str))
wt=list(wt[0].astype(str)+"_"+wt[1].astype(str))
lucky_wt = list(lucky_wt[0].astype(str)+"_"+lucky_wt[1].astype(str))

df_selected["ref"]=df_selected[0].astype(str)+"_"+df_selected[1].astype(str)

df_selected_noback=df_selected[~df_selected["ref"].isin(background)]
df_selected_wt=df_selected[~df_selected["ref"].isin(wt)]
df_selected_lucky = df_selected[~df_selected["ref"].isin(lucky_wt)]


df_selected_noback=df_selected_noback.drop(["ref"],axis=1)
df_selected_wt=df_selected_wt.drop(['ref'], axis=1)
df_selected_lucky=df_selected_lucky.drop(['ref'], axis=1)


df_selected_noback.to_csv(filterfolder+outNoBack, sep='\t', index=None, header=None)
df_selected_wt.to_csv(filterfolder+outWT, sep='\t', index=None, header=None)
df_selected_lucky.to_csv(filterfolder+outLucky, sep='\t', index=None, header=None)


EOF

else
    echo "skipping selection"
    cd ${filter}
    rm ${sample}-F_ACUS-L____${genoem_flag}-___-____-REP_1.SNPs.noback.vcf
    rm ${sample}-F_ACUS-L____${genoem_flag}-___-____-REP_1.SNPs.nowt.vcf
    rm ${sample}-F_ACUS-L____${genoem_flag}-___-____-REP_1.SNPs.lucky.vcf
    ln -s ${sample}-F_ACUS-L____${genoem_flag}-___-____-REP_1.SNPs.vcf ${sample}-F_ACUS-L____${genoem_flag}-___-____-REP_1.SNPs.lucky.vcf
    ln -s ${sample}-F_ACUS-L____${genoem_flag}-___-____-REP_1.SNPs.vcf ${sample}-F_ACUS-L____${genoem_flag}-___-____-REP_1.SNPs.noback.vcf
    ln -s ${sample}-F_ACUS-L____${genoem_flag}-___-____-REP_1.SNPs.vcf ${sample}-F_ACUS-L____${genoem_flag}-___-____-REP_1.SNPs.nowt.vcf
fi

echo "annotation"

cd ${filter}
module load snpeff/4.3.i
java -Xmx4g -jar ${SNPEFF} -ud 0 $genome ${sample}-F_ACUS-L____${genoem_flag}-___-____-REP_1.SNPs.noback.vcf | grep protein_coding | egrep "MODERAT|HIGH"  > ${sample}-F_ACUS-L____${genoem_flag}-___-____-REP_1.SNPs.noback.protein.vcf
java -Xmx4g -jar ${SNPEFF} -ud 0 $genome ${sample}-F_ACUS-L____${genoem_flag}-___-____-REP_1.SNPs.nowt.vcf | grep protein_coding | egrep "MODERAT|HIGH"  > ${sample}-F_ACUS-L____${genoem_flag}-___-____-REP_1.SNPs.nowt.protein.vcf
java -Xmx4g -jar ${SNPEFF} -ud 0 $genome ${sample}-F_ACUS-L____${genoem_flag}-___-____-REP_1.SNPs.lucky.vcf | grep protein_coding | egrep "MODERAT|HIGH"  > ${sample}-F_ACUS-L____${genoem_flag}-___-____-REP_1.SNPs.lucky.protein.vcf

echo 'filtering'
python - << EOF

import pandas as pd
import numpy as np
import scipy.stats as stats


inFolder="${filter}"
nobackProtein=pd.read_table(inFolder+"${sample}-F_ACUS-L____${genoem_flag}-___-____-REP_1.SNPs.noback.protein.vcf",header=None)
nowtProtein=pd.read_table(inFolder+"${sample}-F_ACUS-L____${genoem_flag}-___-____-REP_1.SNPs.nowt.protein.vcf",header=None)
luckyProtein=pd.read_table(inFolder+"${sample}-F_ACUS-L____${genoem_flag}-___-____-REP_1.SNPs.lucky.protein.vcf",header=None)


def Calcs(x):
    l=x.split("AD=")[1].split(";")[0].split(",")
    l=[ int(s) for s in l ]
    #l=l[:-1]
    sums=np.sum(l)
    smaller=l[1]
    res=str(sums)+","+str(smaller)
    return res

def ANNOTATE(D, sampleID, outname, min_reads, min_freq, hits_per_gene):
    D["ANN"]=D[7].apply(lambda x: x.split("ANN=")[1] )
    D["tmp"]=D[7].apply(lambda x: Calcs(x) )
    D["sums"]=D["tmp"].apply(lambda x: x.split(",")[0] )
    D["smaller alternative al."]=D["tmp"].apply(lambda x: x.split(",")[1] )
    D["freq. (smaller alternative al.) (%)"]=D["smaller alternative al."].astype(float)/D["sums"].astype(float)*100.0
    D['genes'] = D["ANN"].apply(lambda x: x.split("|")[3])
    D['sampleID'] = [sampleID]*D.shape[0]
    # keep only desired columns
    D=D[['sampleID', 0,1,3,4,"sums","smaller alternative al.", "freq. (smaller alternative al.) (%)", "genes","ANN"]]
    D.columns=[['sampleID', "chrom","pos","REF","ALT","sums","smaller alternative al.","freq. (smaller alternative al.) (%)", "genes", "ANN"]]
    # filter by allele frequency
    D=D[D["freq. (smaller alternative al.) (%)"]>min_freq]
    # filter by read depths
    D=D[D["sums"].astype(float)>min_reads]
    # filter by number of hits per gene
    sel_genes = list(D['genes'].value_counts()[D['genes'].value_counts() >= hits_per_gene].index)
    D = D[D['genes'].isin(sel_genes)]
    # write table to file
    D.to_csv(inFolder+outname,index=None,sep="\t")
    return D


nobackProtein_lenient=ANNOTATE(nobackProtein, "${sample}", "${sample}-F_ACUS-L____${genoem_flag}-___-____-REP_1.SNPs.noback.protein.freq.lenient.vcf", 30, 3, 1)
nowtProtein_lenient=ANNOTATE(nowtProtein, "${sample}", "${sample}-F_ACUS-L____${genoem_flag}-___-____-REP_1.SNPs.nowt.protein.freq.lenient.vcf", 30, 3, 1)
luckyProtein_lenient=ANNOTATE(luckyProtein, "${sample}", "${sample}-F_ACUS-L____${genoem_flag}-___-____-REP_1.SNPs.lucky.protein.freq.lenient.vcf", 30, 3, 1)


nobackProtein_stringent=ANNOTATE(nobackProtein, "${sample}", "${sample}-F_ACUS-L____${genoem_flag}-___-____-REP_1.SNPs.noback.protein.freq.stringent.vcf", 50, 5, 2)
nowtProtein_stringent=ANNOTATE(nowtProtein, "${sample}", "${sample}-F_ACUS-L____${genoem_flag}-___-____-REP_1.SNPs.nowt.protein.freq.stringent.vcf", 50, 5, 2)
luckyProtein_stringent=ANNOTATE(luckyProtein, "${sample}", "${sample}-F_ACUS-L____${genoem_flag}-___-____-REP_1.SNPs.lucky.protein.freq.stringent.vcf", 50, 5, 2)


EXC=pd.ExcelWriter(inFolder+"results.${sample}.xlsx")
nowtProtein_lenient.to_excel(EXC,sheet_name="minus_wildtype_lenient",index=None) # change according to human once relevant
nowtProtein_stringent.to_excel(EXC,sheet_name="minus_wildtype_stringent",index=None)
nobackProtein_lenient.to_excel(EXC,sheet_name="minus_syn_wt_lenient",index=None)
nobackProtein_stringent.to_excel(EXC,sheet_name="minus_syn_wt_stringent",index=None)
luckyProtein_lenient.to_excel(EXC,sheet_name="lucky_lenient",index=None) 
luckyProtein_stringent.to_excel(EXC,sheet_name="lucky_stringent",index=None)

#kegg.to_excel(EXC,sheet_name="KEGG Pathway analysis",index=None)
EXC.save()

EOF

echo "finished pipeline on sample ${sample}" >> ${logs}/${sample}.pipeline.log
date >> ${logs}/${sample}.pipeline.log

exit
