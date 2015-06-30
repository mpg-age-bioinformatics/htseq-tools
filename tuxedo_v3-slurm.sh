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

series="HaIS"

# Reference genome

ann=/data/genomes/caenorhabditis_elegans/WBcel235_79
ori_GTF=$(readlink -f  ${ann}/WBcel235.*.gtf)
hisat_index=${ann}/hisat/WBcel235.dna.toplevel
adapters_file=full_path_to_adapters_file_for_flexbar_to_use
genome=${ann}/hisat/WBcel235.dna.toplevel.fa

#############################################################################



echo "Creating required folders"
mkdir ../slurm_logs
mkdir ../fastqc
mkdir ../tmp
mkdir ../raw_trimmed
mkdir ../hisat_output
mkdir ../stringtie_output
mkdir ../cuffmerge_output
mkdir ../cuffdiff_output
mkdir ../cuffquant_output


top=$(readlink -f ../)/
tmp=$(readlink -f ../tmp)/
raw=$(readlink -f ../raw_data)/
rawt=$(readlink -f ../raw_trimmed)/
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


if [[ -e ${tmp}HS_ST.ids ]]; then
rm ${tmp}${tmp}HS_ST.ids
fi


cd ${rawt}


for serie in $series; do
for file in $(ls *${serie}*1.fastq.gz); do 


if [[ $(contains "${SE_unstr[@]}" "$serie") == "y" ]]; then
lib=
files="-U ${file}"

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
hisat -p 20 ${lib} --met-file ${top}hisat_output/${file::(-16)}.stats \
-x ${hisat_index} -S ${top}hisat_output/${file::(-16)}.sam \
${files}

cd ${top}hisat_output
module load SAMtools

printf 'Filtering mapped reads and sorting bam file'

samtools view -@ 20 -bhS -F 4 ${file::(-16)}.sam | samtools sort -@ 20 - ${file::(-16)}

mkdir ${top}stringtie_output/${file::(-16)}

printf 'Starting StringTie'

module load StringTie
stringtie ${file::(-16)}.bam -o ${top}stringtie_output/${file::(-16)}.gtf \
-p 20 -G ${ori_GTF} \
-C ${top}stringtie_output/${file::(-16)}_full_cov.gtf \
-b ${top}stringtie_output/${file::(-16)} 

rm ${tmp}HS_ST_${file::(-16)}.sh" > ${tmp}HS_ST_${file::(-16)}.sh

cd ${tmp}
chmod 755 ${tmp}HS_ST_${file::(-16)}.sh 
rm ../slurm_logs/HS_ST_${file::(-16)}.*.out
sbatch -p blade,himem,hugemem --cpus-per-task=20 -o ../slurm_logs/HS_ST_${file::(-16)}.%j.out ${tmp}HS_ST_${file::(-16)}.sh 2>&1 | tee ${tmp}HS_ST_${file::(-16)}.id
id=$(cat ${tmp}HS_ST_${file::(-16)}.id | grep 'Submitted batch job')
echo -n :${id:20} >> ${tmp}HS_ST.ids
rm ${tmp}HS_ST_${file::(-16)}.id

done
done

HS_ST_ids=$(cat ${tmp}HS_ST.ids)
srun -p blade,himem,hugemem -d afterok${HS_ST_ids} echo "HiSat and StringTie done"

#############################################################################

echo "Starting cuffmerge"


for serie in $series; do

cd ${top}stringtie_output

for full_cov in $(ls *${serie}*_full_cov.gtf); do
readlink -f ${full_cov} >> ${tmp}assemblies_${serie}.txt
done

cd ${top}
mkdir cuffmerge_output/${serie}
cmout=$(readlink -f cuffmerge_output/${serie})/
echo ${serie}

cd ${top}
module load Cufflinks
srun -p blade,himem,hugemem --cpus-per-task=2 cuffmerge -p 2 --min-isoform-fraction 1.0 -o ${cmout} \
-g ${ori_GTF} -s ${genome} ${tmp}assemblies_${serie}.txt

cd ${cmout}
mkdir cuffcompare_output
cd cuffcompare_output
echo "STARTING CUFFCOMPARE"
module load Cufflinks
#srun -p himem,hugemem,blade --cpus-per-task=2 cuffcompare -C -r ${ori_GTF} -s ${ann}/chromosomes ${cmout}merged.gtf
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

cd ${top}hisat_output

for file in $(ls *${serie}*.bam); do echo "#!/bin/bash
cd ${top}cuffquant_output
mkdir ${serie}
cd ${serie}
module load Cufflinks
cuffquant -p 20 --library-type ${lib} \
-o ${file::(-4)} \
${top}cuffmerge_output/${serie}/merged.gtf \
${top}hisat_output/${file}
rm ${tmp}quant_${file::(-4)}.sh" > ${tmp}quant_${file::(-4)}.sh

cd ${tmp}
chmod 755 ${tmp}quant_${file::(-4)}.sh
rm ../slurm_logs/quant_${file::(-4)}.*.out
sbatch -p blade,himem,hugemem --cpus-per-task=20 -o ../slurm_logs/quant_${file::(-4)}.%j.out ${tmp}quant_${file::(-4)}.sh 2>&1 | tee ${tmp}quant_${file::(-4)}.id
id=$(cat ${tmp}quant_${file::(-4)}.id | grep 'Submitted batch job')
echo -n :${id:20} >> ${tmp}quant.ids
rm ${tmp}quant_${file::(-4)}.id
done
done
quant_ids=$(cat ${tmp}quant.ids)

srun -p blade,himem,hugemem -d afterok${quant_ids} echo "Starting cuffdiff"


#############################################################################


#### cuff diff >>>> one section per serie ######

serie=HaIS
mkdir ${top}cuffdiff_output/${serie}
dout=$(readlink -f ${top}cuffdiff_output/${serie})
lib="fr-firststrand"

echo "#!/bin/bash
cd ${qua}${serie}

module load Cufflinks

cuffdiff -p 20 --library-type ${lib} \
-L N2,daf2 \
-o ${dout} --dispersion-method per-condition \
${top}cuffmerge_output/${series}/merged.gtf \
S_160-F_HaIS-L____N2-___-____-REP_1/abundances.cxb,S_161-F_HaIS-L____N2-___-____-REP_2/abundances.cxb,S_162-F_HaIS-L____N2-___-____-REP_3/abundances.cxb \
S_163-F_HaIS-L__daf2-___-____-REP_1/abundances.cxb,S_164-F_HaIS-L__daf2-___-____-REP_2/abundances.cxb,S_165-F_HaIS-L__daf2-___-____-REP_3/abundances.cxb

rm ${tmp}diff_${serie}.sh" > ${tmp}diff_${serie}.sh

#### END section

for serie in ${series}; do
cd ${tmp}
chmod 755 ${tmp}diff_${serie}.sh
rm ../slurm_logs/diff_${serie}.*.out
sbatch -p blade,himem,hugemem --mem=512gb --cpus-per-task=20 -o ../slurm_logs/diff_${serie}.%j.out ${tmp}diff_${serie}.sh
done

exit
