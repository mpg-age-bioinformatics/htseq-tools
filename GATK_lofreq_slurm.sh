#!/bin/bash

"
This is a standard GATK and lofreq variant call pipeline with snpEff annotations.

It is made to make calls on two lines eg. wt and mutant.

It is set for single end sequencing.

This script needs to run form inside the folder scripts in a working project with the following structure:

project/scripts
project/raw_data

Raw data needs to be labeled in the following fashion:

Sample_serial-Folder-Line-Time_point/day_of_life(day_post_Treament)-treament-REPlicate-READ

and with the exact number of characters as in the example bellow:

S_001-F_HaTS-L____N2-__0-____-REP_1-READ_1.fastq.gz

S_XXX-F_XXXX-L_XXXXX-XXX-XXXX-REP_X-READ_x.fastq.gz

Please notice that for paired samples, the S_XXX is the same.
" > /dev/null 2>&1

LineA="-L____wt-" # this needs to be the same string that you used in -L_XXXXX- for the 1st line
LineB="-L____mt-" # this needs to be the same string that you used in -L_XXXXX- for the 2nd line

LabelA="wt" # Pick a label for the 1st line
LabelB="mt" # Pick a label for the second line

adapters_file=/path/to/adapters.txt # TruSeq Adapters
oriFasta=/path/to/organism.primary_assembly.fa
faref=pick_a_name_for_the_fasta_reference.fa # name to be attributed to the reference fasta on your project folder
fadic=pick_a_name_for_the_fasta_reference.dict # name to be attributed to the reference dictionary on your project folder
pic=/beegfs/common/software/picard-tools/1.137/picard.jar
vcf=/path/to/known.vcf # VCF file download from ENSEMBL
snpEffdb=snpEff_database

# Create output folders
mkdir ../fastqc_output > /dev/null 2>&1
mkdir ../flexbar_output > /dev/null 2>&1
mkdir ../slurm_logs > /dev/null 2>&1
mkdir ../tmp > /dev/null 2>&1
mkdir ../star_output > /dev/null 2>&1
mkdir ../star_2pass > /dev/null 2>&1
mkdir ../replace_groups > /dev/null 2>&1
mkdir ../sortsam_output > /dev/null 2>&1
mkdir ../markduplicates_output > /dev/null 2>&1
mkdir ../others > /dev/null 2>&1
mkdir ../realign_output > /dev/null 2>&1
mkdir ../recalibrate_output > /dev/null 2>&1
mkdir ../variantscall_output > /dev/null 2>&1
mkdir ../star_2pass_index > /dev/null 2>&1
mkdir ../variant_filtration > /dev/null 2>&1
mkdir ../split_cigars > /dev/null 2>&1
mkdir ../lofreq_output > /dev/null 2>&1


# Paths
top=$(readlink -f ../)/
raw=$(readlink -f ../raw_data)/
filt=$(readlink -f ../variant_filtration)/
qc=$(readlink -f ../fastqc_output)/
logs=$(readlink -f ../slurm_logs)/
tmp=$(readlink -f ../tmp)/
fb=$(readlink -f ../flexbar_output)/
bm=$(readlink -f ../star_output)/
sti2=$(readlink -f ../star_2pass_index)/
st2=$(readlink -f ../star_2pass)/
rg=$(readlink -f ../replace_groups)/
ss=$(readlink -f ../sortsam_output)/
sp=$(readlink -f ../split_cigars)/
md=$(readlink -f ../markduplicates_output)/
others=$(readlink -f ../others)/
real=$(readlink -f ../realign_output)/
rec=$(readlink -f ../recalibrate_output)/
var=$(readlink -f ../variantscall_output)/
lofreq=$(readlink -f ../lofreq_output)/


# Load required software
module load vcftools
module load IGV
module load Java
module load SAMtools
module load STAR
module load FastQC
module load Flexbar
module load R
module load lofreq

snpEff=/path/to/snpEff/4.2/snpEff.jar
GA=/path/to/GATK/3.4-46/GenomeAnalysisTK.jar

# Generate indexes
printf '\nGenerating indexes\n\n'  
mkdir ../others > /dev/null 2&>1
cd ../others

srun vcftools --vcf ${vcf} --out ${vcf::(-4)}_indels --recode --recode-INFO-all --keep-only-indels
deletions=$(readlink -f ${vcf::(-4)}_indels.recode.vcf)

srun igvtools index ${vcf::(-4)}_indels.recode.vcf

srun cp ${oriFasta} ${faref}
srun java -jar ${pic} CreateSequenceDictionary R=${faref} O=${fadic}
srun samtools faidx ${faref}

rm -rf STAR
mkdir STAR
ref=$(readlink -f STAR) 
cd STAR
srun cp ${oriFasta} .
fas=$(ls *.*)
srun --cpus-per-task=18 STARshort --runMode genomeGenerate --genomeDir ${ref} --genomeFastaFiles ${fas} --runThreadN 18
cd ..

faref=$(readlink -f ${faref})
IDS=
cd ${raw}
for f in $(ls *1.fastq.gz);

    do echo "#!/bin/bash
    
    printf '\nFastQC\n\n'
    fastqc -t 4 -o ${qc} ${f}
    
    printf '\nFlexbar\n\n'
    flexbar -r ${f} -t ${fb}${f::(-9)} -n 18 -a ${adapters_file} -ao 10 -u 5 -q 20 -f i1.8 -ae ANY
    
    printf '\nFirst STAR mapping\n\n'
    mkdir ${bm}${f::(-16)}
    cd ${bm}${f::(-16)}
    STARshort --genomeDir ${ref} --readFilesIn ${fb}${f::(-3)} --runThreadN 18
    
    printf '\nSecond STAR genome index\n\n'  
    mkdir ${sti2}${f::(-16)}
    cp ${ref}/${fas} ${sti2}${f::(-16)}
    cd ${sti2}${f::(-16)}
    STARshort --runMode genomeGenerate --genomeDir ${sti2}${f::(-16)} --genomeFastaFiles ${fas} --sjdbFileChrStartEnd ${bm}${f::(-16)}/SJ.out.tab --sjdbOverhang 75 --runThreadN 18 
    mkdir ${st2}${f::(-16)}
    
    printf '\nSecond STAR mapping\n\n' 
    cd  ${st2}${f::(-16)}
    STARshort --genomeDir ${sti2}${f::(-16)} --readFilesIn ${fb}${f::(-3)} --runThreadN 18 
    
    printf '\nSort sam file\n\n' 
    java -jar ${pic} SortSam INPUT=${st2}${f::(-16)}/Aligned.out.sam OUTPUT=${ss}${f::(-16)}.bam SORT_ORDER=coordinate
    
    printf '\nRead groups\n\n' 
    java -jar ${pic} AddOrReplaceReadGroups I=${ss}${f::(-16)}.bam O=${rg}${f::(-16)}.bam SO=coordinate RGID=${f::5} RGLB=lib1 RGPL=illumina RGPU=unit1 RGSM=${f::5}    
    
    printf '\nMark duplicates\n\n' 
    java -jar ${pic} MarkDuplicates I=${rg}${f::(-16)}.bam O=${md}${f::(-16)}.bam CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT M=${md}metrics.txt
    cd ${md}
    
    printf '\nIndexing bam\n\n' 
    java -jar ${pic} BuildBamIndex INPUT=${f::(-16)}.bam
    
    printf '\nSplit Cigar Reads\n\n' 
    java -jar ${GA} -T SplitNCigarReads -R ${faref} -I ${f::(-16)}.bam -o ${sp}${f::(-16)}.bam -rf ReassignOneMappingQuality -RMQF 255 -RMQT 60 -U ALLOW_N_CIGAR_READS
    
    printf '\nRealigner Target Creator\n\n' 
    java -jar ${GA} -nt 18 -T RealignerTargetCreator -known ${deletions} -R ${faref} -I ${sp}${f::(-16)}.bam -o ${others}${f::(-16)}.list  
    
    rm ${tmp}${f}_p1.sh
    " > ${tmp}${f}_p1.sh
    chmod 755 ${tmp}${f}_p1.sh
    rm ${logs}${f}_p1.*.out > /dev/null 2>&1
    ID=$(sbatch -o ${logs}${f}_p1.%j.out --cpus-per-task=18 ${tmp}${f}_p1.sh)
    echo ${ID}
 
    echo "#!/bin/bash
    
    # IndelRealigner has a marginal effect on the number of indels which are recovered,
    # it is memory greedy and therefore it tends to fail and it can be skipped.  
    printf '\nIndel Realigner\n\n' 
    java -Xmx60g -jar ${GA} -T IndelRealigner -known ${deletions} -R ${faref} -I ${sp}${f::(-16)}.bam -targetIntervals ${others}${f::(-16)}.list -o ${real}${f::(-16)}.bam 
    ####rm ${real}${f::(-16)}*
    ####cp ${sp}${f::(-16)}.bam ${real}${f::(-16)}.bam
    cd ${real}
    printf '\nBuild Bam index\n\n'
    java -jar ${pic} BuildBamIndex INPUT=${f::(-16)}.bam

    printf '\nFirst Base Recalibrator\n\n'
    java -Xmx60g -jar ${GA} -T BaseRecalibrator -lowMemoryMode -R ${faref} -I ${real}${f::(-16)}.bam -knownSites ${vcf} -o ${others}${f::(-16)}.table
    
    printf '\nSecond Base Recalibrator\n\n'     
    java -jar ${GA} -T BaseRecalibrator -R ${faref} -I ${real}${f::(-16)}.bam -knownSites ${vcf} -BQSR ${others}${f::(-16)}.table -o ${others}${f::(-16)}_post.table
    
    printf '\nAnalyze Covariates\n\n'    
    java -jar ${GA} -T AnalyzeCovariates -R ${faref} -l DEBUG -before ${others}${f::(-16)}.table -after ${others}${f::(-16)}_post.table -plots ${others}${f::(-16)}.recal_plots.pdf
    
    printf '\nPrint Reads\n\n'  
    java -jar ${GA} -T PrintReads -R ${faref} -I ${real}${f::(-16)}.bam -BQSR ${others}${f::(-16)}.table -o ${rec}${f::(-16)}.bam # check the table input
    
    printf '\nHaplotype Caller\n\n'  
    java -jar ${GA} -T HaplotypeCaller -R ${faref} -I ${rec}${f::(-16)}.bam -dontUseSoftClippedBases -stand_call_conf 20.0 -ERC GVCF --variant_index_type LINEAR --variant_index_parameter 128000 -stand_emit_conf 20.0 -o ${var}${f::(-16)}.vcf  #--genotyping_mode DISCOVERY -stand_emit_conf 10 -stand_call_conf 30 -o ${var}${f::(-16)}.vcf -bamout ${var}${f::(-16)}.bam -ERC GVCF --variant_index_type LINEAR --variant_index_parameter 128000
    
    #printf '\nVariant Filtration\n\n'   
    #java -jar ${GA} -T VariantFiltration -R ${faref} -V ${var}${f::(-16)}.vcf -window 35 -cluster 3 -filterName FS -filter 'FS > 30.0' -filterName QD -filter 'QD < 2.0' -o ${filt}${f::(-16)}.vcf 
    
    printf '\nLoFreq\n'
    lofreq call --call-indels --no-default-filter --verbose -f ${faref} -o ${lofreq}${f::(-16)}.vcf ${rec}${f::(-16)}.bam

    printf '\nsnpEff on LoFreq output\n'
    java -Xmx4g -jar ${snpEff} -ud 0 -v ${snpEffdb} ${lofreq}${f::(-16)}.vcf > ${lofreq}${f::(-16)}.ann.vcf

    rm ${tmp}${f}_p2.sh
    " > ${tmp}${f}_p2.sh
    chmod 755 ${tmp}${f}_p2.sh
    rm ${logs}${f}_p2.*.out > /dev/null 2>&1
    id=$(sbatch -d afterok:${ID:20} --mem=60gb -o ${logs}${f}_p2.%j.out ${tmp}${f}_p2.sh)
    IDS=${IDS}:${id:20}
done

cd ${var} 
GVCFs1=$(for i in $(ls *${LineA}*.vcf); do echo -n "-V ${var}${i} "; done)
GVCFs2=$(for i in $(ls *${LineB}*.vcf); do echo -n "-V ${var}${i} "; done)


# current issues with --filterExpression. As this is not required for reconstruction we skipped it for now. 

echo "#!/bin/bash
cd ${others}
java -jar ${GA} -T GenotypeGVCFs -R ${faref} ${GVCFs1::(-1)} -o ${var}${LabelA}.vcf
java -jar ${GA} -T SelectVariants -R ${faref} -V ${var}${LabelA}.vcf -selectType SNP -o ${var}${LabelA}_raw_snps.vcf
java -jar ${GA} -T SelectVariants -R ${faref} -V ${var}${LabelA}.vcf -selectType INDEL -o ${var}${LabelA}_raw_INDELS.vcf
#java -jar ${GA} -T VariantFiltration -R ${faref} -V ${var}${LabelA}_raw_snps.vcf --filterExpression 'QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5' --filterName '${LabelA}_snps_a' --filterExpression 'ReadPosRankSum < -8.0' --filterName '${LabelA}_snps_b' -o ${var}${LabelA}_filt_snps.vcf
#java -jar ${GA} -T VariantFiltration -R ${faref} -V ${var}${LabelA}_raw_INDELS.vcf --filterExpression 'QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0' --filterName '${LabelA}_indels' -o ${var}${LabelA}_filt_INDELS.vcf
java -jar ${GA} -T SelectVariants -R ${faref} --variant ${var}${LabelA}_raw_snps.vcf -select 'QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0' -o ${var}${LabelA}_PASS_snps.vcf
java -jar ${GA} -T SelectVariants -R ${faref} --variant ${var}${LabelA}_raw_INDELS.vcf -select 'QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0' -o ${var}${LabelA}_PASS_INDELS.vcf
java -cp ${GA} org.broadinstitute.gatk.tools.CatVariants -R ${faref} -V ${var}${LabelA}_PASS_snps.vcf -V ${var}${LabelA}_PASS_INDELS.vcf -out ${LabelA}.pass.vcf -assumeSorted
java -Xmx4g -jar ${snpEff} -ud 0 -v ${snpEffdb} ${LabelA}.pass.vcf > ${LabelA}.pass.ann.vcf
java -jar ${GA} -T FastaAlternateReferenceMaker -R ${faref} -o ${LabelA}.g.fa -V ${LabelA}.pass.vcf
rm ${tmp}${LabelA}.joint.sh
" > ${tmp}${LabelA}.joint.sh
chmod 755 ${tmp}${LabelA}.joint.sh
rm ${logs}${LabelA}.*.out > /dev/null 2>&1
sbatch -d afterok${IDS} -o ${logs}${LabelA}.%j.out ${tmp}${LabelA}.joint.sh

# Line 2
echo "#!/bin/bash
cd ${others}
java -jar ${GA} -T GenotypeGVCFs -R ${faref} ${GVCFs1::(-1)} -o ${var}${LabelB}.vcf
java -jar ${GA} -T SelectVariants -R ${faref} -V ${var}${LabelB}.vcf -selectType SNP -o ${var}${LabelB}_raw_snps.vcf
java -jar ${GA} -T SelectVariants -R ${faref} -V ${var}${LabelB}.vcf -selectType INDEL -o ${var}${LabelB}_raw_INDELS.vcf
#java -jar ${GA} -T VariantFiltration -R ${faref} -V ${var}${LabelB}_raw_snps.vcf --filterExpression 'QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5' --filterName '${LabelB}_snps_a' --filterExpression 'ReadPosRankSum < -8.0' --filterName '${LabelB}_snps_b'  -o ${var}${LabelB}_filt_snps.vcf
#java -jar ${GA} -T VariantFiltration -R ${faref} -V ${var}${LabelB}_raw_INDELS.vcf --filterExpression 'QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0' --filterName '${LabelB}_indels' -o ${var}${LabelB}_filt_INDELS.vcf
java -jar ${GA} -T SelectVariants -R ${faref} --variant ${var}${LabelB}_raw_snps.vcf -select 'QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0' -o ${var}${LabelB}_PASS_snps.vcf
java -jar ${GA} -T SelectVariants -R ${faref} --variant ${var}${LabelB}_raw_INDELS.vcf -select 'QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0' -o ${var}${LabelB}_PASS_INDELS.vcf
java -cp ${GA} org.broadinstitute.gatk.tools.CatVariants -R ${faref} -V ${var}${LabelB}_PASS_snps.vcf -V ${var}${LabelB}_PASS_INDELS.vcf -out ${LabelB}.pass.vcf -assumeSorted
java -Xmx4g -jar ${snpEff} -ud 0 -v ${snpEffdb} ${LabelB}.pass.vcf > ${LabelB}.pass.ann.vcf
java -jar ${GA} -T FastaAlternateReferenceMaker -R ${faref} -o ${LabelB}.g.fa -V ${LabelB}.pass.vcf
rm ${tmp}${LabelB}.joint.sh
" > ${tmp}${LabelB}.joint.sh
chmod 755 ${tmp}${LabelB}.joint.sh
rm ${logs}${LabelB}.*.out
sbatch -d afterok${IDS} -o ${logs}${LabelB}.%j.out ${tmp}${LabelB}.joint.sh

exit

