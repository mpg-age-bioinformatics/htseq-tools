#!/bin/bash

# meran-slurm.sh - a meRanTK pipeline for running RNA methylation analysis on a slurm-HPC system
#
# copyright (c) 2015 - Sven E. Templer <sven.templer@gmail.com>

############
# USER INPUT
############

# root project directory to put output folder into
home=/path/to/base/directory

# notification mail address
email=mail@foo.bar

# path to meRanTK software v 1.1.1a
path_meran=software/meRanTK
# path to samtools binary, script is compatible with version 0.1.18
bin_sam=$path_meran/extutil/x86_64-linux-thread-multi/tophat2/samtools_0.1.18

# reference genome fasta and gff files from tophat2 transcriptome assembly
index_fa=/beegfs/common/genomes/homo_sapiens/83/toplevel_tophat2/index.fa
index_gtf=/beegfs/common/genomes/homo_sapiens/83/toplevel_tophat2/index.fa.gff

# sample 'sheet' - to create use e.g. vim with ':r! ls reads_trim/*' for example
samples="
# comment lines start with '#' and are ignored
# empty lines are ignored
# fields are white-space separeated

# sample information as
# 1) group = group identifier, must contain 2 different in total
# 2) read_fwd = path (relative to \$home or absolute) to forward reads
# 3) read_rev = optional path to reverse reads

young reads_trim/s49_1.fastq reads_trim/s49_2.fastq
young reads_trim/s50_1.fastq reads_trim/s50_2.fastq
young reads_trim/s51_1.fastq reads_trim/s51_2.fastq
young reads_trim/s52_1.fastq reads_trim/s52_2.fastq
young reads_trim/s53_1.fastq reads_trim/s53_2.fastq
young reads_trim/s54_1.fastq reads_trim/s54_2.fastq

old reads_trim/s55_1.fastq reads_trim/s55_2.fastq
#old reads_trim/s56_1.fastq reads_trim/s56_2.fastq
old reads_trim/s57_1.fastq reads_trim/s57_2.fastq
old reads_trim/s58_1.fastq reads_trim/s58_2.fastq
old reads_trim/s59_1.fastq reads_trim/s59_2.fastq
old reads_trim/s60_1.fastq reads_trim/s60_2.fastq
"

############
# environment
############

file_check () {
  [[ ! -f $1 ]] && echo "error: file '$1' does not exist" && exit
}

cat_array () { 
  local IFS=','
  echo "$*"
}

submit_sbatch () {
  local pid
  if [[ -f ${1}.pid ]]
  then
    echo -n "* skipped (already submitted) " >&2
  else
    echo "$2" > ${1}.sbatch
    sbatch --parsable ${1}.sbatch > ${1}.pid
    echo -n "* submitted " >&2
  fi
  pid=$(cat ${1}.pid)
  echo "${1}.sbatch ($pid)" >&2
  echo $pid
  return
}

unset dep dep2
unset samples_group samples_name samples_fq1 samples_fq2
unset samples_name_uniq samples_name_a samples_name_b

cd $home
folders=(meran_index meran_align meran_call meran_compare)
mkdir -p ${folders[@]} ${folders[@]/#/slurm/}

while read -r line
do
  [[ -z $line ]] && continue
  [[ -n $(echo "$line" | egrep -o "^\ *#|^\ *$") ]] && continue
  line=($line)
  group=${line[0]}
  fq1=${line[1]}
  fq2=${line[2]}
  name=${fq1%.gz}
  name=${name%.fastq}
  name=${name%.fq}
  name=$(basename $name)
  file_check $fq1
  [[ -n $fq2 ]] && file_check $fq2
  samples_group+=("$group")
  samples_fq1+=("$fq1")
  samples_fq2+=("$fq2")
  samples_name+=("$name")
done <<< "$samples"

groups=($(printf '%s\n' "${samples_group[@]}" | sort | uniq))
[[ ${#groups[@]} -ne 2 ]] && \
  echo "error: expecting 2 samples groups, but found ${#groups[@]}" && \
  exit
comparison=${groups[0]}_vs_${groups[1]}
mkdir -p meran_compare/$comparison
samples=("name group file_fwd file_rev")
for i in ${!samples_name[@]}
do
  samples+=("${samples_name[$i]} ${samples_group[$i]} ${samples_fq1[$i]} ${samples_fq2[$i]}")
done
echo "* samples to analyse:"
printf '%s\n' "${samples[@]}" | column -t | tee meran_compare/$comparison/samples

############
# index
############

# !! note: for refseq use:
# wget -q -r ftp://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/mRNA_Prot/human.*.rna.fna.g
# gunzip gunzip ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/mRNA_Prot/* 
# cat ftp.../*fna > meran_index/index.fa
# $path_meran/util/mkRefSeq2GeneMap.pl -f meran_index/index.fa -m meran_index/index.fa.map

sbatch="#!/bin/bash
#SBATCH -n 18
#SBATCH -o slurm/meran_index/index.out
#SBATCH -J meran_index
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=$email
ln -s $index_fa meran_index/index.fa
ln -s $index_gtf meran_index/index.gtf
awk '\$3==\"transcript\"{print \$14, \$18, \$4}' meran_index/index.gtf | tr -d '\";' > meran_index/index.map
$path_meran/meRanT mkbsidx -fa meran_index/index.fa -id meran_index -t 18
"

dep=$(submit_sbatch slurm/meran_index/index "$sbatch")

############
# align, call
############

for i in ${!samples_name[@]}
do
  
  name=${samples_name[$i]}

  ###########
  # align
  ##########
  
  fq1="-f ${samples_fq1[$i]}"
  fq2=${samples_fq2[$i]}
  [[ -n $fq2 ]] && fq2="-r $fq2"
  mkdir -p meran_align/$name

  sbatch="#!/bin/bash
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 18
#SBATCH -o slurm/meran_align/${name}.out
#SBATCH -J meran_align_$name
#SBATCH -d afterok:$dep
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=$email

$path_meran/meRanT align -t 18 \
  -o meran_align/$name -S aligned.sam \
  -un -ud meran_align/$name/unaligned \
  -MM -mbp -k 10 -ra \
  -i2g meran_index/index.map \
  -x meran_index/index_C2T \
  $fq1 $fq2

if [[ ! -f meran_align/$name/aligned_cut.sam ]]; then cut -f1-11 meran_align/$name/aligned.sam > meran_align/$name/aligned_cut.sam; fi
if [[ ! -f meran_align/$name/aligned.bam ]]; then $bin_sam view -bS -o meran_align/$name/aligned.bam meran_align/$name/aligned_cut.sam; fi
if [[ ! -f meran_align/$name/aligned_sorted.bam ]]; then $bin_sam sort meran_align/$name/aligned.bam meran_align/$name/aligned_sorted; fi
$bin_sam index meran_align/$name/aligned_sorted.bam meran_align/$name/aligned_sorted.bai
"

  pid=$(submit_sbatch slurm/meran_align/$name "$sbatch")

  ###########
  # call
  ##########

  mkdir -p meran_call/$name

  sbatch="#!/bin/bash
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 38
#SBATCH -o slurm/meran_call/${name}.out
#SBATCH -J meran_call_$name
#SBATCH -d afterok:$pid
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=$email

$path_meran/meRanCall -p 38 \
  -o meran_call/$name \
  -bam meran_align/$name/aligned_sorted.bam \
  -f meran_index/index.fa \
  -fs5 5 -rs5 5 -rl 50 -sc 10 \
  -ei 0.05 -cr 0.99 -fdr 0.05 -tref
"

  pid=$(submit_sbatch slurm/meran_call/$name "$sbatch")
  dep2="${dep2}:$pid"

done

###########
# compare
##########

group_a=${groups[0]}
group_b=${groups[1]}

files_a=($(awk -v group=$group_a 'NR>1 && $2==group {print $1}' meran_compare/$comparison/samples))
files_a=($(printf 'meran_call/%s_FDR_0.05.txt ' ${files_a[@]}))
files_a=$(cat_array ${files_a[@]})

files_b=($(awk -v group=$group_b 'NR>1 && $2==group {print $1}' meran_compare/$comparison/samples))
files_b=($(printf 'meran_call/%s_FDR_0.05.txt ' ${files_b[@]}))
files_b=$(cat_array ${files_b[@]})

sbatch="#!/bin/bash
#SBATCH -o slurm/meran_compare/${comparison}.out
#SBATCH -J meran_compare_$comparison
#SBATCH -d afterok$dep2
#SBATCH --mail-type=END
#SBATCH --mail-user=$email

cd meran_compare/$comparison
$path_meran/meRanCompare \
  -mr 3 -s 0.01 -fdr 0.01 \
  -na $group_a -nb $group_b \
  -fa $files_a \
  -fb $files_b
"

pid=$(submit_sbatch slurm/meran_compare/$comparison "$sbatch")

echo "* done - notification email will be send to '$email'"
echo "  (after last job ended and if any job failed)"

