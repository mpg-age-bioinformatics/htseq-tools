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
mix=("Yid3")

unstr=("YiTS" "YiDR" "YiIS" "ShTe")
str=("Yid1" "OeDA" "RoSt" "HaTS" "HaIS" "AgMi" )
#mix=("Yid3")


# Which series do you which to work on:

series="AgMi"

# Reference genome

ann=/data/genomes/mus_musculus/GRCm38_79
ori_GTF=$(readlink -f  ${ann}/GRCm38.*.gtf)
GTF_file=$(readlink -f ${ann}/cuffcmp_GTF.*.gtf)
GTF_index=$(readlink -f ${ann}/cuffcmp_GTF_index/*.gff)
GTF_index=${GTF_index::(-4)}
genome=$(ls ${ann}/bowtie2/*dna.toplevel.fa)
genomeN=${genome::(-3)}

adapters_file=full_path_to_adapters_file_for_flexbar_to_use


#############################################################################



echo "Creating required folders"
mkdir ../slurm_logs
mkdir ../fastqc
mkdir ../tmp
mkdir ../tophat_output
mkdir ../cufflinks_output
mkdir ../raw_trimmed
mkdir ../cuffmerge_output
mkdir ../cuffquant_output
mkdir ../cuffdiff_output

top=$(readlink -f ../)/
tmp=$(readlink -f ../tmp)/
raw=$(readlink -f ../raw_data)/
rawt=$(readlink -f ../raw_trimmed)/
cli=$(readlink -f ../cufflinks_output)/
qua=$(readlink -f ../cuffquant_output)/
merg=$(readlink -f ../cuffmerge_output)/

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
unpigz -p 10 ${file}
fastqc -t 10 -o ../fastqc ${file::(-3)}
rm ${tmp}fastqc_${file::(-9)}.sh" > ${tmp}fastqc_${file::(-9)}.sh

cd ${tmp} 
chmod 755 ${tmp}fastqc_${file::(-9)}.sh 
rm ../slurm_logs/fastqc_${file::(-9)}.*.out
sbatch -p blade,himem,hugemem --cpus-per-task=10 -o ../slurm_logs/fastqc_${file::(-9)}.%j.out ${tmp}fastqc_${file::(-9)}.sh 2>&1 | tee ${tmp}fastqc_${file::(-2)}id

id=$(cat ${tmp}fastqc_${file::(-2)}id | grep 'Submitted batch job')
echo -n :${id:20} >> ${tmp}fastqc.ids
rm ${tmp}fastqc_${file::(-2)}id
done
done

fastqc_ids=$(cat ${tmp}fastqc.ids)
srun -p blade,himem,hugemem -d afterok${fastqc_ids} echo "FASTQC done"

#############################################################################

cd ${tmp}

if [[ -e ${tmp}flexbar.ids ]]; then
rm ${tmp}flexbar.ids
fi


for serie in $series; do
for file in $(ls *${serie}*1.fastq); do

if [[ -e ${file::(-7)}2.fastq ]]; then

echo "#!/bin/bash 
module load pigz
module load Flexbar

flexbar -r ${tmp}${file::(-7)}1.fastq \
-p ${tmp}${file::(-7)}2.fastq -t ${top}raw_trimmed/${file::(-8)} \
-n 20 -a ${adapters_file} \
-ao 10 -u 5 -q 20 -m 20 -f i1.8 -ae ANY
cd ${top}raw_trimmed
pigz -p 10 ${file::(-7)}1.fastq
pigz -p 10 ${file::(-7)}2.fastq

rm ${tmp}flexbar_${file::(-5)}sh" > ${tmp}flexbar_${file::(-5)}sh

else

echo "#!/bin/bash
module load pigz
module load Flexbar

flexbar -r ${tmp}${file::(-7)}1.fastq \
-t ${top}raw_trimmed/${file::(-6)} \
-n 20 -a -a ${adapters_file} \
-ao 10 -u 5 -q 20 -m 20 -f i1.8 -ae ANY
cd ${top}raw_trimmed
pigz -p 10 ${file}
rm ${tmp}flexbar_${file::(-5)}sh" > ${tmp}flexbar_${file::(-5)}sh

fi

cd ${tmp}
chmod 755 ${tmp}flexbar_${file::(-5)}sh
rm ../slurm_logs/flexbar_${file::(-5)}*.out
sbatch -p blade,himem,hugemem --cpus-per-task=20 -o ../slurm_logs/flexbar_${file::(-5)}%j.out ${tmp}flexbar_${file::(-5)}sh 2>&1 | tee ${tmp}flexbar_${file::(-5)}id
id=$(cat ${tmp}flexbar_${file::(-5)}id | grep 'Submitted batch job')
echo -n :${id:20} >> ${tmp}flexbar.ids
rm ${tmp}flexbar_${file::(-5)}id
done
done

flexbar_ids=$(cat ${tmp}flexbar.ids)
srun -p blade,himem,hugemem -d afterok${flexbar_ids} echo "FLEXBAR done"


#############################################################################


if [[ -e ${tmp}assemblies.txt ]]; then
rm ${tmp}assemblies.txt
fi
if [[ -e ${tmp}th_cl.ids ]]; then
rm ${tmp}${tmp}th_cl.ids
fi


cd ${rawt}


for serie in $series; do
for file in $(ls *${serie}*1.fastq.gz); do 


if [[ $(contains "${SE_unstr[@]}" "$serie") == "y" ]]; then
lib="fr-unstranded"
files=${file}

elif [[ $(contains "${SE_str[@]}" "$serie") == "y" ]]; then
lib="fr-firststrand"
files=${file}

elif [[ $(contains "${PE_str[@]}" "$serie") == "y" ]]; then
lib="fr-firststrand"
files="${file} ${file::(-10)}2.fastq.gz"

elif [[ $(contains "${mix[@]}" "$serie") == "y" ]]; then
files=${file}
REP=${file:30:5}

if [[ ${REP} == REP_3 ]]; then
lib="fr-firststrand"
else
lib="fr-unstranded"
fi
fi

echo "#!/bin/bash
cd ${rawt}
module load Bowtie2
module load TopHat
tophat -p 20 --library-type ${lib} \
--transcriptome-index ${GTF_index} \
-o ${top}tophat_output/${file::(-16)} \
${genomeN} \
${files}

module load Cufflinks
cufflinks -p 20 --library-type ${lib} \
-g ${GTF_file} --no-faux-reads \
-o ${top}cufflinks_output/${file::(-16)} \
${top}tophat_output/${file::(-16)}/accepted_hits.bam

cat ${top}cufflinks_output/${file::(-16)}/transcripts.gtf | grep yes > ${top}cufflinks_output/${file::(-16)}/transcripts_full_read.gtf

rm ${tmp}th_cl_${file::(-16)}.sh" > ${tmp}th_cl_${file::(-16)}.sh

cd ${tmp}
chmod 755 ${tmp}th_cl_${file::(-16)}.sh 
rm ../slurm_logs/th_cl_${file::(-16)}.*.out
sbatch -p blade,himem,hugemem --cpus-per-task=20 -o ../slurm_logs/th_cl_${file::(-16)}.%j.out ${tmp}th_cl_${file::(-16)}.sh 2>&1 | tee ${tmp}th_cl_${file::(-16)}.id
id=$(cat ${tmp}th_cl_${file::(-16)}.id | grep 'Submitted batch job')
echo -n :${id:20} >> ${tmp}th_cl.ids
rm ${tmp}th_cl_${file::(-16)}.id
echo "cufflinks_output/${file::(-16)}/transcripts_full_read.gtf" >> ${tmp}assemblies.txt

done
done

th_cl_ids=$(cat ${tmp}th_cl.ids)
srun -p blade,himem,hugemem -d afterok${th_cl_ids} echo "TopHat and Cufflinks done"

#############################################################################


echo "Starting cuffmerge"

for serie in $series; do

cat ${tmp}assemblies.txt | grep ${serie} > ${tmp}assemblies_${serie}.txt

cd ${top}
mkdir cuffmerge_output/${serie}
cmout=$(readlink -f cuffmerge_output/${serie})/
echo ${serie}

cd ${top}
module load Cufflinks
srun -p blade,himem,hugemem --cpus-per-task=2 cuffmerge -p 2 -o ${cmout} \
-g ${GTF_file} -s ${genome} ${tmp}assemblies_${serie}.txt

cd ${cmout}
mkdir cuffcompare_output
cd cuffcompare_output
echo "STARTING CUFFCOMPARE"
module load Cufflinks
srun -p himem,hugemem,blade --cpus-per-task=2 cuffcompare -C -r ${ori_GTF} -s ${ann}/chromosomes ${cmout}merged.gtf
done

cd ${tmp}
cd ../scripts

#############################################################################

echo "Starting cuffquant"

if [[ -e ${tmp}quant.ids ]]; then
rm ${tmp}quant.ids
fi


for serie in $series; do

if [[ $(contains "${unstr[@]}" "$serie") == "y" ]]; then
lib="fr-unstranded"

elif [[ $(contains "${str[@]}" "$serie") == "y" ]]; then
lib="fr-firststrand"


elif [[ $(contains "${mix[@]}" "$serie") == "y" ]]; then
lib="fr-unstranded"

fi

cd ${top}tophat_output

for file in $(ls -d *${serie}*); do echo "#!/bin/bash
cd ${top}cuffquant_output
mkdir ${serie}
cd ${serie}
module load Cufflinks
cuffquant -p 20 --library-type ${lib} \
-o ${file} \
${top}cuffmerge_output/${serie}/cuffcompare_output/cuffcmp.combined.gtf \
${top}tophat_output/${file}/accepted_hits.bam
rm ${tmp}quant_${file}.sh" > ${tmp}quant_${file}.sh

cd ${tmp}
chmod 755 ${tmp}quant_${file}.sh
rm ../slurm_logs/quant_${file}.*.out
sbatch -p blade,himem,hugemem --cpus-per-task=20 -o ../slurm_logs/quant_${file}.%j.out ${tmp}quant_${file}.sh 2>&1 | tee ${tmp}quant_${file}.id
id=$(cat ${tmp}quant_${file}.id | grep 'Submitted batch job')
echo -n :${id:20} >> ${tmp}quant.ids
rm ${tmp}quant_${file}.id
done
done
quant_ids=$(cat ${tmp}quant.ids)

srun -p blade,himem,hugemem -d afterok${quant_ids} echo "Starting cuffdiff"


#############################################################################

echo "Starting cuffdiff"

#### cuff diff >>>> one section per serie ######

serie=AgMi
mkdir ${top}cuffdiff_output/${serie}
dout=$(readlink -f ${top}cuffdiff_output/${serie})
lib="fr-firststrand"

echo "#!/bin/bash
cd ${qua}${serie}

module load Cufflinks

cuffdiff -p 20 --library-type ${lib} \
-L liver_10M,heart_10M,cereb_10M,hippo_10M,liver_27M,heart_27M,cereb_27M,hippo_27M \
-o ${dout} --dispersion-method per-condition \
${top}cuffmerge_output/${serie}/cuffcompare_output/cuffcmp.combined.gtf \
S_001-F_AgMi-L_live-_10-____-REP_1/abundances.cxb,S_005-F_AgMi-L_live-_10-____-REP_2/abundances.cxb,S_009-F_AgMi-L_live-_10-____-REP_3/abundances.cxb \
S_002-F_AgMi-L_hear-_10-____-REP_1/abundances.cxb,S_006-F_AgMi-L_hear-_10-____-REP_2/abundances.cxb,S_010-F_AgMi-L_hear-_10-____-REP_3/abundances.cxb \
S_003-F_AgMi-L_cere-_10-____-REP_1/abundances.cxb,S_007-F_AgMi-L_cere-_10-____-REP_2/abundances.cxb,S_011-F_AgMi-L_cere-_10-____-REP_3/abundances.cxb \
S_004-F_AgMi-L_hipp-_10-____-REP_1/abundances.cxb,S_008-F_AgMi-L_hipp-_10-____-REP_2/abundances.cxb,S_012-F_AgMi-L_hipp-_10-____-REP_3/abundances.cxb \
S_013-F_AgMi-L_live-_27-____-REP_1/abundances.cxb,S_017-F_AgMi-L_live-_27-____-REP_2/abundances.cxb,S_021-F_AgMi-L_live-_27-____-REP_3/abundances.cxb \
S_014-F_AgMi-L_hear-_27-____-REP_1/abundances.cxb,S_018-F_AgMi-L_hear-_27-____-REP_2/abundances.cxb,S_022-F_AgMi-L_hear-_27-____-REP_3/abundances.cxb \
S_015-F_AgMi-L_cere-_27-____-REP_1/abundances.cxb,S_019-F_AgMi-L_cere-_27-____-REP_2/abundances.cxb,S_023-F_AgMi-L_cere-_27-____-REP_3/abundances.cxb \
S_016-F_AgMi-L_hipp-_27-____-REP_1/abundances.cxb,S_020-F_AgMi-L_hipp-_27-____-REP_2/abundances.cxb,S_024-F_AgMi-L_hipp-_27-____-REP_3/abundances.cxb

rm ${tmp}diff_${serie}.sh" > ${tmp}diff_${serie}.sh

#### END section



for serie in ${series}; do
cd ${tmp}
chmod 755 ${tmp}diff_${serie}.sh
rm ../slurm_logs/diff_${serie}.*.out
sbatch -p blade,himem,hugemem --mem=512gb --cpus-per-task=20 -o ../slurm_logs/diff_${serie}.%j.out ${tmp}diff_${serie}.sh
done

exit
