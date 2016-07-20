#!/bin/bash

"This script needs to run form inside the folder scripts in a working project with the following structure:
project/scripts
project/raw_data


Raw data needs to be labeled in the following fashion:

Sample_serial-Folder-Line-Time_point/day_of_life(day_post_Treament)-treament-REPlicate-READ

and with the exact number of characters as in the example bellow:

S_001-F_HaTS-L____N2-__0-____-REP_1-READ_1.fastq.gz

S_XXX-F_XXXX-L_XXXXX-XXX-XXXX-REP_X-READ_x.fastq.gz

Please notice that for paired samples, the S_XXX is the same.


Make sure you have edited the last section of this script - cuffdiff - before you execute this script." > /dev/null 2>&1


#############################################################################

# Define series as SE or PE and stranded or unstranded

SE_unstr=("YiTS" "YiDR" "YiIS" "ShTe")
SE_str=("Yid1" "OeDA" "AgMi")
PE_str=("RoSt" "HaTS" "HaIS")
PE_uns=("XHFC")
mix=("Yid3")

unstr=("YiTS" "YiDR" "YiIS" "ShTe" "XHFC")
str=("Yid1" "OeDA" "RoSt" "HaTS" "HaIS" "AgMi" )
#mix=("Yid3")


# Which series do you which to work on:

series="XHFC"

# Reference genome

ann=/data/genomes/mus_musculus/GRCm38_79
ori_GTF=$(readlink -f ${ann}/GRCm38.79.gtf)
hisat_index=${ann}/hisat/GRCm38.dna.toplevel.fa
adapters_file=/beegfs/group_bit/home/JBoucas/documents/TruSeqAdapters.txt
genome=${ann}/hisat/GRCm38.dna.toplevel.fa


#############################################################################


echo "Creating required folders"
mkdir ../slurm_logs
mkdir ../fastqc_output
mkdir ../tmp
mkdir ../flexbar_output
mkdir ../hisat_output
mkdir ../stringtie_output
mkdir ../cuffmerge_output
mkdir ../cuffdiff_output
mkdir ../cuffquant_output


top=$(readlink -f ../)/
tmp=$(readlink -f ../tmp)/
raw=$(readlink -f ../raw_data)/
rawt=$(readlink -f ../flexbar_output)/
merg=$(readlink -f ../cuffmerge_output)/ 
qua=$(readlink -f ../cuffquant_output)/ 



# Required function

function contains() {
    local n=$#
    local value=${!n}
    for ((i=1;i < $#;i++)) {
        if [ "${!i}" == "${value}" ]; then
            echo "y"
            return 0
        fi
    }
    echo "n"
    return 1
}

#############################################################################

echo "Starting FASTQC"

ids=

cd ${raw} 
for serie in $series; do
    cd ${raw}
    for file in $(ls *${serie}*.fastq.gz); do 
        echo "#!/bin/bash
        module load pigz
        module load FastQC
        cp ${raw}${file} ${tmp}
        cd ${tmp}
        unpigz -p 4 ${file}
        # FASTQC call
        fastqc -t 4 -o ../fastqc_output ${file::(-3)}
        rm ${tmp}fastqc_${file::(-9)}.sh
        " > ${tmp}fastqc_${file::(-9)}.sh

        cd ${tmp} 
        chmod 755 ${tmp}fastqc_${file::(-9)}.sh 
        rm ../slurm_logs/fastqc_${file::(-9)}.*.out > /dev/null 2>&1  
        id=$(sbatch -p blade,himem,hugemem --cpus-per-task=4 -o ../slurm_logs/fastqc_${file::(-9)}.%j.out ${tmp}fastqc_${file::(-9)}.sh)
        sleep 2
        ids=${ids}:${id:20}
    done
done

echo "Waiting for FASTQC jobs${ids} to complete"
srun -p blade,himem,hugemem -d afterok${ids} echo "FASTQC done. Starting Flexbar."

#############################################################################

ids=

cd ${raw}
for serie in $series; do
    cd ${raw}
    for file in $(ls *${serie}*1.fastq.gz); do
        cd ${raw}
        if [[ -e ${file::(-10)}2.fastq.gz ]]; then
            echo "#!/bin/bash 
            module load pigz
            module load Flexbar

            # Flexbar call for paired end reads

            flexbar -r ${tmp}${file::(-10)}1.fastq \
            -p ${tmp}${file::(-10)}2.fastq -t ${top}flexbar_output/${file::(-11)} \
            -n 18 -a ${adapters_file} \
            -ao 10 -u 5 -q 20 -m 20 -f i1.8 -ae ANY
            cd ${top}flexbar_output
            pigz -p 18 ${file::(-10)}1.fastq
            pigz -p 18 ${file::(-10)}2.fastq
            rm ${tmp}flexbar_${file::(-8)}sh
            " > ${tmp}flexbar_${file::(-8)}sh
        else
            echo "#!/bin/bash
            module load pigz
            module load Flexbar

            # Flexbar call for single end reads

            flexbar -r ${tmp}${file::(-10)}1.fastq \
            -t ${top}flexbar_output/${file::(-11)}_1 \
            -n 18 -a ${adapters_file} \
            -ao 10 -u 5 -q 20 -m 20 -f i1.8 -ae ANY
            cd ${top}flexbar_output
            pigz -p 18 ${file::(-3)}
            rm ${tmp}flexbar_${file::(-8)}sh
            " > ${tmp}flexbar_${file::(-8)}sh

        fi

    cd ${tmp}
    chmod 755 ${tmp}flexbar_${file::(-8)}sh
    rm ../slurm_logs/flexbar_${file::(-8)}*.out > /dev/null 2>&1  
    id=$(sbatch -p blade,himem,hugemem --cpus-per-task=18 -o ../slurm_logs/flexbar_${file::(-8)}%j.out ${tmp}flexbar_${file::(-8)}sh)
    sleep 2
    ids=${ids}:${id:20}    
    done
done

echo "Waiting for Flexbar jobs${ids} to complete"
srun -p blade,himem,hugemem -d afterok${ids} echo "FLEXBAR done. Starting HiSat and StringTie"


#############################################################################

ids=

cd ${rawt}
for serie in $series; do
    cd ${rawt}
    for file in $(ls *${serie}*1.fastq.gz); do
        
        # Libraries and paired end vs. single end settings for HISAT    
 
        if [[ $(contains "${SE_unstr[@]}" "$serie") == "y" ]]; then
            lib=
            files="-U ${file}"
        elif [[ $(contains "${PE_uns[@]}" "$serie") == "y" ]]; then
            lib=
            files="-1 ${file} -2 ${file::(-10)}2.fastq.gz"
        elif [[ $(contains "${SE_str[@]}" "$serie") == "y" ]]; then
            lib="--rna-strandness R"
            files="-U ${file}"
        elif [[ $(contains "${PE_str[@]}" "$serie") == "y" ]]; then
            lib="--rna-strandness RF"
            files="-1 ${file} -2 ${file::(-10)}2.fastq.gz"
        elif [[ $(contains "${mix[@]}" "$serie") == "y" ]]; then
            files=-U ${file}
            REP=${file:30:5}
            if [[ ${REP} == REP_3 ]]; then
                lib="--rna-strandness R"
            else
                lib=
            fi
        fi

        echo "#!/bin/bash
        cd ${rawt}
        module load Bowtie2
        module load HISAT

        # HISAT call 

        hisat2 -p 18 ${lib} --dta --met-file ${top}hisat_output/${file::(-16)}.stats \
        -x ${hisat_index} -S ${top}hisat_output/${file::(-16)}.sam \
        ${files}

        cd ${top}hisat_output
        module load SAMtools
        
        # Use samtools to select mapped reads and sort them

        samtools view -@ 18 -bhS -F 4 ${file::(-16)}.sam | samtools sort -@ 18 ${file::(-16)}.bam -
        mkdir ${top}stringtie_output/${file::(-16)}

        module load StringTie

        # StringTie call

        stringtie ${file::(-16)}.bam -o ${top}stringtie_output/${file::(-16)}.gtf \
        -p 18 -G ${ori_GTF} -f 0.99 \
        -C ${top}stringtie_output/${file::(-16)}_full_cov.gtf \
        -b ${top}stringtie_output/${file::(-16)} 
        rm ${tmp}HS_ST_${file::(-16)}.sh
        " > ${tmp}HS_ST_${file::(-16)}.sh

        cd ${tmp}
        chmod 755 ${tmp}HS_ST_${file::(-16)}.sh 
        rm ../slurm_logs/HS_ST_${file::(-16)}.*.out > /dev/null 2>&1
        id=$(sbatch -p blade,himem,hugemem --cpus-per-task=18 -o ../slurm_logs/HS_ST_${file::(-16)}.%j.out ${tmp}HS_ST_${file::(-16)}.sh)
        sleep 2
        ids=${ids}:${id:20}
    done
done

echo "Waiting for HISAT and StringTie jobs${ids} to complete"
srun -p blade,himem,hugemem -d afterok${ids} echo "HiSat and StringTie done. Starting cuffmerge"
 
#############################################################################

for serie in $series; do
    rm ${tmp}assemblies_${serie}.txt

    cd ${top}stringtie_output
    mkdir full_coverage
    mv *_full_cov.gtf full_coverage
    
    # Select only transcripts which have full coverage

    cd full_coverage
    for gtf in $(ls *${serie}*.gtf); do
        readlink -f ${gtf} >> ${tmp}assemblies_${serie}.txt
    done

    cd ${top}
    mkdir cuffmerge_output/${serie}
    cmout=$(readlink -f cuffmerge_output/${serie})/
    echo ${serie}

    cd ${top}
    module load Cufflinks
    
    # Cuffmerge call

    srun -p blade,himem,hugemem --cpus-per-task=2 cuffmerge -p 2 \
    -o ${cmout} --min-isoform-fraction 1.0 \
    -g ${ori_GTF} -s ${genome} ${tmp}assemblies_${serie}.txt
done

cd ${tmp}
cd ../scripts

#############################################################################

echo "Starting cuffquant"

ids=

for serie in $series; do

    # Library settings for cuffquant

    if [[ $(contains "${unstr[@]}" "$serie") == "y" ]]; then
        lib="fr-unstranded"
    elif [[ $(contains "${str[@]}" "$serie") == "y" ]]; then
        lib="fr-firststrand"
    elif [[ $(contains "${mix[@]}" "$serie") == "y" ]]; then
        lib="fr-unstranded"
    fi

    cd ${top}hisat_output
    for file in $(ls *${serie}*.bam); do 
        echo "#!/bin/bash
        cd ${top}cuffquant_output
        mkdir ${serie}
        cd ${serie}
        module load Cufflinks

        # Cuffquant call

        cuffquant -p 18 --library-type ${lib} \
        -o ${file::(-4)} \
        ${top}cuffmerge_output/${serie}/merged.gtf \
        ${top}hisat_output/${file}
        rm ${tmp}quant_${file::(-4)}.sh
        " > ${tmp}quant_${file::(-4)}.sh
        cd ${tmp}
        chmod 755 ${tmp}quant_${file::(-4)}.sh
        rm ../slurm_logs/quant_${file::(-4)}.*.out > /dev/null 2>&1  
        id=$(sbatch -p blade,himem,hugemem --cpus-per-task=18 -o ../slurm_logs/quant_${file::(-4)}.%j.out ${tmp}quant_${file::(-4)}.sh)
        sleep 2
        ids=${ids}:${id:20}
    done
done

echo "Waiting for cuffquant jobs${ids} to complete"
srun -p blade,himem,hugemem -d afterok${ids} echo "Cuffquant done. Starting cuffdiff."


#############################################################################


#### cuff diff >>>> one section per serie ######

serie=XHFC
mkdir ${top}cuffdiff_output/${serie}
dout=$(readlink -f ${top}cuffdiff_output/${serie})
lib="fr-unstranded"

echo "#!/bin/bash
cd ${qua}${serie}

module load Cufflinks

# Cuffdiff call

cuffdiff -p 18 --library-type ${lib} \
-L Fema_Y,Male_Y,Male_O \
-o ${dout} --dispersion-method per-condition \
${top}cuffmerge_output/${series}/merged.gtf \
S_002-F_XHFC-L_____F-__Y-____-REP_2/abundances.cxb,S_003-F_XHFC-L_____F-__Y-____-REP_3/abundances.cxb,S_004-F_XHFC-L_____F-__Y-____-REP_4/abundances.cxb,S_006-F_XHFC-L_____F-__Y-____-REP_6/abundances.cxb,S_007-F_XHFC-L_____F-__Y-____-REP_7/abundances.cxb,S_008-F_XHFC-L_____F-__Y-____-REP_8/abundances.cxb \
S_009-F_XHFC-L_____M-__Y-____-REP_1/abundances.cxb,S_010-F_XHFC-L_____M-__Y-____-REP_2/abundances.cxb,S_011-F_XHFC-L_____M-__Y-____-REP_3/abundances.cxb,S_012-F_XHFC-L_____M-__Y-____-REP_4/abundances.cxb,S_013-F_XHFC-L_____M-__Y-____-REP_5/abundances.cxb,S_014-F_XHFC-L_____M-__Y-____-REP_6/abundances.cxb,S_015-F_XHFC-L_____M-__Y-____-REP_7/abundances.cxb \
S_017-F_XHFC-L_____M-__O-____-REP_1/abundances.cxb,S_018-F_XHFC-L_____M-__O-____-REP_2/abundances.cxb,S_019-F_XHFC-L_____M-__O-____-REP_3/abundances.cxb,S_020-F_XHFC-L_____M-__O-____-REP_4/abundances.cxb,S_021-F_XHFC-L_____M-__O-____-REP_5/abundances.cxb

rm ${tmp}diff_${serie}.sh" > ${tmp}diff_${serie}.sh

#### END section

for serie in ${series}; do
    cd ${tmp}
    chmod 755 ${tmp}diff_${serie}.sh
    rm ../slurm_logs/diff_${serie}.*.out > /dev/null 2>&1  
    sbatch -p blade,himem,hugemem --mem=724gb --cpus-per-task=18 -o ../slurm_logs/diff_${serie}.%j.out ${tmp}diff_${serie}.sh
    sleep 2
done

exit
