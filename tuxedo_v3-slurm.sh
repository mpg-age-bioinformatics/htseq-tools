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


Make sure you have edited the last section of this script - cuffdiff - before you execute this script."


#############################################################################

# Define series as SE or PE and stranded or unstranded

SE_unstr=("YiTS" "YiDR" "YiIS" "ShTe")
SE_str=("Yid1" "OeDA" "AgMi")
PE_str=("RoSt" "HaTS" "HaIS")
PE_uns=("XuFC")
mix=("Yid3")

unstr=("YiTS" "YiDR" "YiIS" "ShTe" "XuFC")
str=("Yid1" "OeDA" "RoSt" "HaTS" "HaIS" "AgMi" )
#mix=("Yid3")


# Which series do you which to work on:

series="XuFC"

# Reference genome

ann=/data/genomes/mus_musculus/GRCm38_79
ori_GTF=$(readlink -f  ${ann}/GRCm38.79.gtf)
hisat_index=${ann}/hisat/GRCm38.dna.toplevel.fa
adapters_file=/beegfs/group_bit/home/JBoucas/documents/TruSeqAdapters.txt
genome=${ann}/hisat/GRCm38.dna.toplevel.fa


#############################################################################


echo "Creating required folders"
mkdir ../slurm_logs
mkdir ../fastqc
mkdir ../tmp
mkdir ../raw_trimmed
mkdir ../V3_hisat_output
mkdir ../V3_stringtie_output
mkdir ../V3_cuffmerge_output
mkdir ../V3_cuffdiff_output
mkdir ../V3_cuffquant_output


top=$(readlink -f ../)/
tmp=$(readlink -f ../tmp)/
raw=$(readlink -f ../raw_data)/
rawt=$(readlink -f ../raw_trimmed)/
merg=$(readlink -f ../V3_cuffmerge_output)/ 
qua=$(readlink -f ../V3_cuffquant_output)/ 



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

cd ${raw} 

if [[ -e ${tmp}fastqc.ids ]]; then
rm ${tmp}fastqc.ids
fi

for serie in $series; do
for file in $(ls *${serie}*.fastq.gz); do echo "#!/bin/bash
module load pigz
module load FastQC
cp ${raw}${file} ${tmp}
cd ${tmp}
unpigz -p 4 ${file}
fastqc -t 4 -o ../fastqc ${file::(-3)}
rm ${tmp}fastqc_${file::(-9)}.sh" > ${tmp}fastqc_${file::(-9)}.sh

cd ${tmp} 
chmod 755 ${tmp}fastqc_${file::(-9)}.sh 
rm ../slurm_logs/fastqc_${file::(-9)}.*.out
sbatch -p blade,himem,hugemem --cpus-per-task=4 -o ../slurm_logs/fastqc_${file::(-9)}.%j.out ${tmp}fastqc_${file::(-9)}.sh 2>&1 | tee ${tmp}fastqc_${file::(-2)}id

id=$(cat ${tmp}fastqc_${file::(-2)}id | grep 'Submitted batch job')

echo -n :${id:20} >> ${tmp}fastqc.ids
rm ${tmp}fastqc_${file::(-2)}id
done
done

fastqc_ids=$(cat ${tmp}fastqc.ids)
srun -p blade,himem,hugemem -d afterok${fastqc_ids} echo "FASTQC done"

#############################################################################

cd ${raw}

if [[ -e ${tmp}flexbar.ids ]]; then
rm ${tmp}flexbar.ids
fi


for serie in $series; do
for file in $(ls *${serie}*1.fastq.gz); do

if [[ -e ${file::(-10)}2.fastq.gz ]]; then

echo "#!/bin/bash 
module load pigz
module load Flexbar
flexbar -r ${tmp}${file::(-10)}1.fastq \
-p ${tmp}${file::(-10)}2.fastq -t ${top}raw_trimmed/${file::(-11)} \
-n 18 -a ${adapters_file} \
-ao 10 -u 5 -q 20 -m 20 -f i1.5 -ae ANY
cd ${top}raw_trimmed
pigz -p 18 ${file::(-10)}1.fastq
pigz -p 18 ${file::(-10)}2.fastq
rm ${tmp}flexbar_${file::(-8)}sh" > ${tmp}flexbar_${file::(-8)}sh

else

echo "#!/bin/bash
module load pigz
module load Flexbar
flexbar -r ${tmp}${file::(-10)}1.fastq \
-t ${top}raw_trimmed/${file::(-11)}_1 \
-n 18 -a ${adapters_file} \
-ao 10 -u 5 -q 20 -m 20 -f i1.5 -ae ANY
cd ${top}raw_trimmed
pigz -p 18 ${file::(-3)}
rm ${tmp}flexbar_${file::(-8)}sh" > ${tmp}flexbar_${file::(-8)}sh

fi

cd ${tmp}
chmod 755 ${tmp}flexbar_${file::(-8)}sh
rm ../slurm_logs/flexbar_${file::(-8)}*.out
sbatch -p blade,himem,hugemem --cpus-per-task=18 -o ../slurm_logs/flexbar_${file::(-8)}%j.out ${tmp}flexbar_${file::(-8)}sh 2>&1 | tee ${tmp}flexbar_${file::(-8)}id

id=$(cat ${tmp}flexbar_${file::(-8)}id | grep 'Submitted batch job')

echo -n :${id:20} >> ${tmp}flexbar.ids
rm ${tmp}flexbar_${file::(-8)}id
done
done

flexbar_ids=$(cat ${tmp}flexbar.ids)
srun -p blade,himem,hugemem -d afterok${flexbar_ids} echo "FLEXBAR done"


#############################################################################


if [[ -e ${tmp}V3_HS_ST.ids ]]; then
rm ${tmp}V3_HS_ST.ids
fi


cd ${rawt}


for serie in $series; do
for file in $(ls *${serie}*1.fastq.gz); do 


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
hisat -p 18 ${lib} --met-file ${top}V3_hisat_output/${file::(-16)}.stats \
-x ${hisat_index} -S ${top}V3_hisat_output/${file::(-16)}.sam \
${files}

cd ${top}V3_hisat_output
module load SAMtools

printf 'Filtering mapped reads and sorting bam file'

samtools view -@ 18 -bhS -F 4 ${file::(-16)}.sam | samtools sort -@ 18 - ${file::(-16)}
mkdir ${top}V3_stringtie_output/${file::(-16)}

printf 'Starting StringTie'

module load StringTie
stringtie ${file::(-16)}.bam -o ${top}V3_stringtie_output/${file::(-16)}.gtf \
-p 18 -G ${ori_GTF} \
-C ${top}V3_stringtie_output/${file::(-16)}_full_cov.gtf \
-b ${top}V3_stringtie_output/${file::(-16)} 
rm ${tmp}V3_HS_ST_${file::(-16)}.sh" > ${tmp}V3_HS_ST_${file::(-16)}.sh

cd ${tmp}
chmod 755 ${tmp}V3_HS_ST_${file::(-16)}.sh 
rm ../slurm_logs/V3_HS_ST_${file::(-16)}.*.out
sbatch -p blade,himem,hugemem --cpus-per-task=18 -o ../slurm_logs/V3_HS_ST_${file::(-16)}.%j.out ${tmp}V3_HS_ST_${file::(-16)}.sh 2>&1 | tee ${tmp}V3_HS_ST_${file::(-16)}.id

id=$(cat ${tmp}V3_HS_ST_${file::(-16)}.id | grep 'Submitted batch job')

echo -n :${id:20} >> ${tmp}V3_HS_ST.ids
rm ${tmp}V3_HS_ST_${file::(-16)}.id

done
done

HS_ST_ids=$(cat ${tmp}V3_HS_ST.ids)
srun -p blade,himem,hugemem -d afterok${HS_ST_ids} echo "HiSat and StringTie done"
 
#############################################################################

echo "Starting cuffmerge"


for serie in $series; do

rm ${tmp}V3_assemblies_${serie}.txt

cd ${top}V3_stringtie_output
mkdir full_coverage
mv *_full_cov.gtf full_coverage

for gtf in $(ls *${serie}*.gtf); do
readlink -f ${gtf} >> ${tmp}V3_assemblies_${serie}.txt
done

cd ${top}
mkdir V3_cuffmerge_output/${serie}
cmout=$(readlink -f V3_cuffmerge_output/${serie})/
echo ${serie}

cd ${top}
module load Cufflinks
srun -p blade,himem,hugemem --cpus-per-task=2 cuffmerge -p 2 \
-o ${cmout} --min-isoform-fraction 1.0 \
-g ${ori_GTF} -s ${genome} ${tmp}V3_assemblies_${serie}.txt

done

cd ${tmp}
cd ../scripts

#############################################################################

echo "Starting cuffquant"

if [[ -e ${tmp}V3_quant.ids ]]; then
rm ${tmp}V3_quant.ids
fi


for serie in $series; do

if [[ $(contains "${unstr[@]}" "$serie") == "y" ]]; then
lib="fr-unstranded"

elif [[ $(contains "${str[@]}" "$serie") == "y" ]]; then
lib="fr-firststrand"


elif [[ $(contains "${mix[@]}" "$serie") == "y" ]]; then
lib="fr-unstranded"

fi

cd ${top}V3_hisat_output

for file in $(ls *${serie}*.bam); do echo "#!/bin/bash

cd ${top}V3_cuffquant_output
mkdir ${serie}
cd ${serie}
module load Cufflinks
cuffquant -p 18 --library-type ${lib} \
-o ${file::(-4)} \
${top}V3_cuffmerge_output/${serie}/merged.gtf \
${top}V3_hisat_output/${file}
rm ${tmp}V3_quant_${file::(-4)}.sh" > ${tmp}V3_quant_${file::(-4)}.sh

cd ${tmp}
chmod 755 ${tmp}V3_quant_${file::(-4)}.sh
rm ../slurm_logs/V3_quant_${file::(-4)}.*.out
sbatch -p blade,himem,hugemem --cpus-per-task=18 -o ../slurm_logs/V3_quant_${file::(-4)}.%j.out ${tmp}V3_quant_${file::(-4)}.sh 2>&1 | tee ${tmp}V3_quant_${file::(-4)}.id
id=$(cat ${tmp}V3_quant_${file::(-4)}.id | grep 'Submitted batch job')
echo -n :${id:20} >> ${tmp}V3_quant.ids
rm ${tmp}V3_quant_${file::(-4)}.id
done
done

quant_ids=$(cat ${tmp}V3_quant.ids)

srun -p blade,himem,hugemem -d afterok${quant_ids} echo "Starting cuffdiff"


#############################################################################


#### cuff diff >>>> one section per serie ######

serie=XuFC
mkdir ${top}V3_cuffdiff_output/${serie}
dout=$(readlink -f ${top}V3_cuffdiff_output/${serie})
lib="fr-unstranded"

echo "#!/bin/bash
cd ${qua}${serie}

module load Cufflinks
cuffdiff -p 18 --library-type ${lib} \
-L Fema_Y,Male_Y,Male_O \
-o ${dout} --dispersion-method per-condition \
${top}V3_cuffmerge_output/${series}/merged.gtf \
S_002-F_XuFC-L_____F-__Y-____-REP_2/abundances.cxb,S_003-F_XuFC-L_____F-__Y-____-REP_3/abundances.cxb,S_004-F_XuFC-L_____F-__Y-____-REP_4/abundances.cxb,S_006-F_XuFC-L_____F-__Y-____-REP_6/abundances.cxb,S_007-F_XuFC-L_____F-__Y-____-REP_7/abundances.cxb,S_008-F_XuFC-L_____F-__Y-____-REP_8/abundances.cxb \
S_009-F_XuFC-L_____M-__Y-____-REP_1/abundances.cxb,S_010-F_XuFC-L_____M-__Y-____-REP_2/abundances.cxb,S_011-F_XuFC-L_____M-__Y-____-REP_3/abundances.cxb,S_012-F_XuFC-L_____M-__Y-____-REP_4/abundances.cxb,S_013-F_XuFC-L_____M-__Y-____-REP_5/abundances.cxb,S_014-F_XuFC-L_____M-__Y-____-REP_6/abundances.cxb,S_015-F_XuFC-L_____M-__Y-____-REP_7/abundances.cxb \
S_017-F_XuFC-L_____M-__O-____-REP_1/abundances.cxb,S_018-F_XuFC-L_____M-__O-____-REP_2/abundances.cxb,S_019-F_XuFC-L_____M-__O-____-REP_3/abundances.cxb,S_020-F_XuFC-L_____M-__O-____-REP_4/abundances.cxb,S_021-F_XuFC-L_____M-__O-____-REP_5/abundances.cxb

rm ${tmp}V3_diff_${serie}.sh" > ${tmp}V3_diff_${serie}.sh

#### END section

for serie in ${series}; do
cd ${tmp}
chmod 755 ${tmp}V3_diff_${serie}.sh
rm ../slurm_logs/V3_diff_${serie}.*.out
sbatch -p blade,himem,hugemem --mem=724gb --cpus-per-task=18 -o ../slurm_logs/V3_diff_${serie}.%j.out ${tmp}V3_diff_${serie}.sh
done

exit
